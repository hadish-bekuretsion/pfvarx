"""Statistical models for allele-count evidence.

All likelihoods are evaluated in log space using ``scipy.special.gammaln`` /
``betaln`` so that zero counts, deep sites and extreme allele fractions remain
finite and deterministic.

Models
------
A. Binomial allele-count model.
B. Beta-binomial overdispersion model, parameterised by mean ``mu`` and
   intra-class correlation ``rho`` (``rho -> 0`` recovers the binomial).
C. Likelihood-ratio comparison of reference-only (H0) and alternate-allele
   (H1) hypotheses, with parametric-bootstrap or boundary chi-square null
   calibration.
D. Optional Bayesian posterior, conditionally enabled.

Scientific caveat: none of these numbers is a validated probability of a true
biological variant. They are evidence measures that must be calibrated against
independent control and truth data before any report-level interpretation.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import optimize, stats
from scipy.special import betaln, gammaln, logsumexp

from pfvarx.errors import ModelError

_EPS = 1e-12


# ---------------------------------------------------------------------------
# Basic log-likelihoods
# ---------------------------------------------------------------------------
def log_binomial_pmf(k: int, n: int, p: float) -> float:
    """log P(K=k | n, p), stable at p in {0, 1} and k in {0, n}."""
    if n < 0 or k < 0 or k > n:
        return -math.inf
    p = min(max(float(p), 0.0), 1.0)
    if p == 0.0:
        return 0.0 if k == 0 else -math.inf
    if p == 1.0:
        return 0.0 if k == n else -math.inf
    return (
        gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
        + k * math.log(p) + (n - k) * math.log1p(-p)
    )


def log_beta_binomial_pmf(k: int, n: int, mu: float, rho: float) -> float:
    """log P(K=k | n, mu, rho) for the beta-binomial.

    ``mu`` is the mean allele probability and ``rho`` the intra-class
    correlation in (0, 1). ``rho == 0`` returns the binomial limit.
    """
    if n < 0 or k < 0 or k > n:
        return -math.inf
    if not (0.0 <= mu <= 1.0):
        raise ModelError(f"beta-binomial mu must be in [0,1], got {mu}")
    if not (0.0 <= rho < 1.0):
        raise ModelError(f"beta-binomial rho must be in [0,1), got {rho}")
    if rho < _EPS:
        return log_binomial_pmf(k, n, mu)
    if mu <= 0.0:
        return 0.0 if k == 0 else -math.inf
    if mu >= 1.0:
        return 0.0 if k == n else -math.inf
    alpha = mu * (1.0 - rho) / rho
    beta = (1.0 - mu) * (1.0 - rho) / rho
    return (
        gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
        + betaln(k + alpha, n - k + beta) - betaln(alpha, beta)
    )


def log_dirichlet_multinomial_pmf(counts, alpha) -> float:
    """log P(counts | alpha) for the Dirichlet-multinomial (multi-allelic sites)."""
    counts = np.asarray(counts, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    if counts.shape != alpha.shape:
        raise ModelError("counts and alpha must have the same shape")
    if np.any(counts < 0) or np.any(alpha <= 0):
        raise ModelError("counts must be non-negative and alpha strictly positive")
    n = counts.sum()
    a0 = alpha.sum()
    return float(
        gammaln(n + 1) - gammaln(counts + 1).sum()
        + gammaln(a0) - gammaln(n + a0)
        + (gammaln(counts + alpha) - gammaln(alpha)).sum()
    )


# ---------------------------------------------------------------------------
# Hypothesis construction
# ---------------------------------------------------------------------------
def error_mixed_alt_probability(true_fraction: float, error_rate: float) -> float:
    """Observed alternate-read probability under a true allele fraction.

    A read supports ALT either because it came from an ALT-bearing genome
    (probability ``f``) or because of a technical error at a REF template
    (probability ``(1-f) * e``). This keeps the biological allele fraction
    distinct from the read-support probability.
    """
    f = min(max(float(true_fraction), 0.0), 1.0)
    e = min(max(float(error_rate), 0.0), 1.0)
    return f + (1.0 - f) * e


@dataclass
class HypothesisLikelihoods:
    n: int
    k: int
    error_rate: float
    alt_fraction: float
    p0: float
    p1: float
    log_l0: float
    log_l1: float
    model: str
    rho: float = 0.0
    alt_fraction_mle: float | None = None
    p1_mle: float | None = None

    @property
    def lrt_statistic(self) -> float:
        return 2.0 * (self.log_l1 - self.log_l0)

    @property
    def log_likelihood_ratio(self) -> float:
        return self.log_l1 - self.log_l0

    def as_dict(self) -> dict:
        d = asdict(self)
        d["lrt_statistic"] = self.lrt_statistic
        d["log_likelihood_ratio"] = self.log_likelihood_ratio
        return d


def evaluate_hypotheses(
    k: int,
    n: int,
    *,
    error_rate: float,
    alt_fraction: float,
    model: str = "binomial",
    rho: float = 0.0,
) -> HypothesisLikelihoods:
    """Evaluate H0 against the composite H1: true AF >= ``alt_fraction``.

    The alternative likelihood is maximised over the allowed AF interval. In
    the binomial model the constrained MLE is analytic; for beta-binomial data
    a bounded scalar search maximises the concave likelihood in its mean
    parameter and compares both interval endpoints. ``alt_fraction`` remains
    the configured lower boundary; ``alt_fraction_mle`` records the fitted value.
    """
    if n <= 0 or k < 0 or k > n:
        raise ModelError(f"Invalid counts k={k}, n={n} (depth must be > 0)")
    p0 = min(max(error_rate, 0.0), 1.0)
    e = min(max(float(error_rate), 0.0), 1.0)
    f_min = min(max(float(alt_fraction), 0.0), 1.0)
    p_min = error_mixed_alt_probability(f_min, e)
    if model == "binomial":
        p1 = max(p_min, k / n)
        f_mle = f_min if e >= 1.0 else min(1.0, max(f_min, (p1 - e) / (1.0 - e)))
        ll0 = log_binomial_pmf(k, n, p0)
        ll1 = log_binomial_pmf(k, n, p1)
    elif model == "beta_binomial":
        def alt_log_likelihood(f: float) -> float:
            return log_beta_binomial_pmf(
                k, n, error_mixed_alt_probability(f, e), rho
            )

        candidates = [f_min, 1.0]
        if f_min < 1.0:
            result = optimize.minimize_scalar(
                lambda f: -alt_log_likelihood(f),
                bounds=(f_min, 1.0),
                method="bounded",
                options={"xatol": 1e-10},
            )
            if result.success:
                candidates.append(float(result.x))
        f_mle = max(candidates, key=alt_log_likelihood)
        p1 = error_mixed_alt_probability(f_mle, e)
        ll0 = log_beta_binomial_pmf(k, n, p0, rho)
        ll1 = log_beta_binomial_pmf(k, n, p1, rho)
    else:
        raise ModelError(f"Unknown count model {model!r}")
    return HypothesisLikelihoods(
        n=n, k=k, error_rate=error_rate, alt_fraction=alt_fraction,
        p0=p0, p1=p1, log_l0=ll0, log_l1=ll1, model=model, rho=rho,
        alt_fraction_mle=f_mle, p1_mle=p1,
    )


def one_sided_count_pvalue(
    k: int, n: int, *, error_rate: float, model: str = "binomial", rho: float = 0.0
) -> float:
    """Exact upper-tail p-value for excess ALT counts under the configured null.

    Validity is conditional on the supplied null parameters. The caller applies
    a Bonferroni bound over the complete bounded hypothesis universe, so
    candidate filtering does not redefine the testing family.
    """
    if n <= 0 or k < 0 or k > n:
        raise ModelError(f"Invalid counts k={k}, n={n} (depth must be > 0)")
    p0 = min(max(float(error_rate), 0.0), 1.0)
    if model == "binomial" or (model == "beta_binomial" and rho < _EPS):
        p = float(stats.binom.sf(k - 1, n, p0))
    elif model == "beta_binomial":
        if not 0.0 <= rho < 1.0:
            raise ModelError("rho must be in [0,1)")
        if p0 in (0.0, 1.0):
            p = 1.0 if (p0 == 0.0 and k == 0) or p0 == 1.0 else 0.0
        else:
            concentration = (1.0 - rho) / rho
            p = float(stats.betabinom.sf(
                k - 1, n, p0 * concentration, (1.0 - p0) * concentration
            ))
    else:
        raise ModelError(f"Unknown count model {model!r}")
    # Prevent numerical underflow of a positive tail from fabricating zero.
    if p == 0.0 and 0.0 < p0 < 1.0:
        p = float(np.nextafter(0.0, 1.0))
    return min(1.0, max(0.0, p))


# ---------------------------------------------------------------------------
# Null calibration
# ---------------------------------------------------------------------------
def lrt_pvalue_bootstrap(
    k: int,
    n: int,
    *,
    error_rate: float,
    alt_fraction: float,
    model: str = "binomial",
    rho: float = 0.0,
    replicates: int = 200,
    seed: int = 20240101,
) -> float:
    """Parametric-bootstrap p-value for the LRT.

    The null distribution is generated under H0 with the *same* n and error
    process, and the LRT is recomputed for each simulated site using the same
    constrained composite H1. This avoids relying on chi-square asymptotics when the null lies
    on the boundary, parameters are estimated, or sites were selected.
    """
    if replicates <= 0:
        raise ModelError("bootstrap replicates must be positive")
    obs = evaluate_hypotheses(k, n, error_rate=error_rate, alt_fraction=alt_fraction,
                              model=model, rho=rho).lrt_statistic
    rng = np.random.default_rng(seed)
    if model == "binomial":
        sims = rng.binomial(n, min(max(error_rate, 0.0), 1.0), size=replicates)
    else:
        mu = min(max(error_rate, 0.0), 1.0)
        if rho < _EPS:
            sims = rng.binomial(n, mu, size=replicates)
        else:
            alpha = mu * (1.0 - rho) / rho
            beta = (1.0 - mu) * (1.0 - rho) / rho
            probs = rng.beta(alpha, beta, size=replicates)
            sims = rng.binomial(n, probs)
    exceed = 0
    for kk in sims:
        stat = evaluate_hypotheses(
            int(kk), n, error_rate=error_rate, alt_fraction=alt_fraction,
            model=model, rho=rho,
        ).lrt_statistic
        if stat >= obs - 1e-12:
            exceed += 1
    if exceed == 0:
        # The observed statistic exceeds every null replicate: the Monte-Carlo
        # estimate has hit its resolution floor of 1/(B+1).  Rather than switch
        # to a *different* null calibration (which would mix two incompatible
        # reference distributions), compute the exact conditional null tail by
        # enumerating K ~ null over 0..n.  This is exact for the stated null and
        # avoids anti-conservative chi-square extrapolation.
        return lrt_pvalue_exact(
            k, n, error_rate=error_rate, alt_fraction=alt_fraction,
            model=model, rho=rho,
        )
    return (1.0 + exceed) / (1.0 + replicates)


def lrt_pvalue_exact(
    k: int,
    n: int,
    *,
    error_rate: float,
    alt_fraction: float,
    model: str = "binomial",
    rho: float = 0.0,
) -> float:
    """Exact conditional null tail of the LRT statistic.

    Computes ``P_null(LRT(K) >= LRT(k))`` by summing the null probability mass
    over every possible count ``K`` in ``0..n``.  The null is the same
    error-only model used for H0 (binomial or beta-binomial), so this is exact
    conditional on ``n``, ``error_rate`` and ``rho``.  It is used to refine the
    Monte-Carlo bootstrap when the observed statistic exceeds all replicates.
    """
    if n <= 0 or k < 0 or k > n:
        raise ModelError(f"Invalid counts k={k}, n={n} (depth must be > 0)")
    obs = evaluate_hypotheses(
        k, n, error_rate=error_rate, alt_fraction=alt_fraction, model=model, rho=rho
    ).lrt_statistic
    total = 0.0
    for kk in range(n + 1):
        stat = evaluate_hypotheses(
            kk, n, error_rate=error_rate, alt_fraction=alt_fraction,
            model=model, rho=rho,
        ).lrt_statistic
        if stat >= obs - 1e-12:
            total += math.exp(_logpmf(kk, n, error_rate, model, rho))
    return min(1.0, max(0.0, total))


def _logpmf(k: int, n: int, p: float, model: str, rho: float) -> float:
    if model == "binomial":
        return log_binomial_pmf(k, n, p)
    return log_beta_binomial_pmf(k, n, p, rho)


def lrt_pvalue_chi2_mixture(statistic: float) -> float:
    """Boundary-corrected p-value: 0.5 * chi2_1 tail (Self & Liang 1987).

    Offered as an analysis option only; not the default because selection and
    estimated nuisance parameters violate the regularity conditions.
    """
    if statistic <= 0:
        return 1.0
    return float(0.5 * stats.chi2.sf(statistic, df=1))


# ---------------------------------------------------------------------------
# Bayesian posterior (optional)
# ---------------------------------------------------------------------------
def bayesian_posterior(
    k: int,
    n: int,
    *,
    error_rate: float,
    alt_fraction: float,
    prior_h1: float,
    model: str = "binomial",
    rho: float = 0.0,
) -> dict:
    """Posterior P(H1 | D) with an explicit, user-supplied prior.

    The prior is kept separate from the likelihood in the returned record so
    that reports never conflate prior belief with observed read evidence.
    """
    if not (0.0 <= prior_h1 <= 1.0):
        raise ModelError("prior_h1 must be in [0,1]")
    hyp = evaluate_hypotheses(k, n, error_rate=error_rate, alt_fraction=alt_fraction,
                              model=model, rho=rho)
    log_prior = np.array([math.log1p(-prior_h1) if prior_h1 < 1 else -math.inf,
                          math.log(prior_h1) if prior_h1 > 0 else -math.inf])
    log_post = np.array([hyp.log_l0, hyp.log_l1]) + log_prior
    norm = logsumexp(log_post)
    post = np.exp(log_post - norm)
    return {
        "posterior_h0": float(post[0]),
        "posterior_h1": float(post[1]),
        "prior_h1": prior_h1,
        "log_marginal": float(norm),
        "bayes_factor_h1_h0": float(math.exp(min(hyp.log_l1 - hyp.log_l0, 700.0))),
        "model": model,
    }


# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------
@dataclass
class CalibrationProfile:
    mode: str
    error_rate: float
    rho: float
    model: str
    n_sites: int
    n_observations: int
    error_rate_ci: tuple[float, float]
    rho_ci: tuple[float, float] | None
    overdispersion_tested: bool
    overdispersion_present: bool
    goodness_of_fit_pvalue: float | None
    dataset_id: str
    status: str
    warnings: list[str] = field(default_factory=list)
    version: str = "1.0.0"

    def as_dict(self) -> dict:
        return asdict(self)


def fit_error_rate(alt_counts, depths, *, confidence: float = 0.95) -> dict:
    """MLE of a single error probability from negative-control sites.

    Uses the pooled binomial MLE ``sum(k)/sum(n)`` with a Jeffreys
    (Beta(0.5, 0.5)) equal-tailed credible interval, which is well behaved at
    zero counts.
    """
    k = np.asarray(alt_counts, dtype=float)
    n = np.asarray(depths, dtype=float)
    if k.shape != n.shape:
        raise ModelError("alt_counts and depths must have equal length")
    if np.any(k < 0) or np.any(n < 0) or np.any(k > n):
        raise ModelError("require 0 <= alt_counts <= depths")
    total_k = float(k.sum())
    total_n = float(n.sum())
    if total_n <= 0:
        raise ModelError("Cannot estimate an error rate from zero total depth")
    p = total_k / total_n
    alpha = 0.5 + total_k
    beta = 0.5 + (total_n - total_k)
    lo = float(stats.beta.ppf((1 - confidence) / 2, alpha, beta))
    hi = float(stats.beta.ppf(1 - (1 - confidence) / 2, alpha, beta))
    return {
        "error_rate": p,
        "ci_low": lo,
        "ci_high": hi,
        "total_alt": total_k,
        "total_depth": total_n,
        "n_sites": int(len(k)),
        "method": "pooled_binomial_mle_jeffreys_ci",
    }


def fit_beta_binomial_rho(alt_counts, depths, mu: float | None = None) -> dict:
    """Fit overdispersion ``rho`` by maximising the beta-binomial likelihood.

    A profile-likelihood interval is returned. If the fitted rho is at the
    boundary (0) the binomial model is preferred.
    """
    k = np.asarray(alt_counts, dtype=float)
    n = np.asarray(depths, dtype=float)
    if k.shape != n.shape or len(k) < 3:
        raise ModelError("beta-binomial calibration requires >=3 paired observations")
    if mu is None:
        mu = float(k.sum() / max(n.sum(), 1.0))
    mu = min(max(mu, 1e-6), 1 - 1e-6)

    def negll(rho: float) -> float:
        if rho <= 0:
            return -sum(log_binomial_pmf(int(kk), int(nn), mu) for kk, nn in zip(k, n, strict=True))
        return -sum(
            log_beta_binomial_pmf(int(kk), int(nn), mu, rho) for kk, nn in zip(k, n, strict=True)
        )

    grid = np.concatenate([[0.0], np.linspace(1e-4, 0.95, 240)])
    vals = np.array([negll(r) for r in grid])
    best = float(grid[int(np.argmin(vals))])
    min_ll = float(vals.min())
    # Profile-likelihood 95% interval: drop <= 1.92 log-likelihood units.
    keep = grid[vals <= min_ll + 1.92]
    rho_ci = (float(keep.min()), float(keep.max())) if len(keep) else (best, best)

    # Likelihood-ratio test of rho=0 (boundary, so 0.5*chi2_1).
    # ll0 is the log-likelihood under H0 (rho=0); -min_ll is the maximised
    # log-likelihood.  The statistic is 2*(maxLL - ll0) and must be >= 0.
    ll0 = -negll(0.0)
    ll_max = -min_ll
    lr = max(0.0, 2.0 * (ll_max - ll0))
    p_overdisp = lrt_pvalue_chi2_mixture(lr)

    return {
        "rho": best,
        "mu": mu,
        "rho_ci_low": rho_ci[0],
        "rho_ci_high": rho_ci[1],
        "overdispersion_lr": lr,
        "overdispersion_pvalue": p_overdisp,
        "overdispersion_present": bool(p_overdisp < 0.05 and best > 0.0),
        "n_observations": int(len(k)),
        "method": "beta_binomial_mle_profile_likelihood",
    }


def goodness_of_fit(alt_counts, depths, *, mu: float, rho: float) -> dict:
    """Pearson-style dispersion statistic against the fitted model.

    For the binomial, the Pearson chi-square/df ratio indicates extra-binomial
    variation. For the beta-binomial the same statistic is reported as a
    diagnostic (not a formal p-value).
    """
    k = np.asarray(alt_counts, dtype=float)
    n = np.asarray(depths, dtype=float)
    if len(k) == 0:
        raise ModelError("goodness_of_fit requires observations")
    if rho <= 0:
        var = n * mu * (1 - mu)
        resid = (k - n * mu) ** 2 / np.maximum(var, _EPS)
    else:
        var = n * mu * (1 - mu) * (n + (1 - rho) / rho) / (1 + (1 - rho) / rho)
        resid = (k - n * mu) ** 2 / np.maximum(var, _EPS)
    df = max(len(k) - 1, 1)
    stat = float(resid.sum())
    return {
        "pearson_chi2": stat,
        "df": df,
        "dispersion_ratio": stat / df,
        "pvalue": float(stats.chi2.sf(stat, df)) if rho <= 0 else None,
    }


def calibrate_from_controls(
    alt_counts,
    depths,
    *,
    mode: str = "unspecified",
    dataset_id: str = "unspecified",
    model: str = "auto",
    confidence: float = 0.95,
) -> CalibrationProfile:
    """Build a mode-specific calibration profile from negative-control sites.

    ``model='auto'`` selects the beta-binomial only when a likelihood-ratio
    test demonstrates overdispersion; otherwise the binomial is retained.
    """
    err = fit_error_rate(alt_counts, depths, confidence=confidence)
    warnings: list[str] = []
    rho = 0.0
    rho_ci = None
    over_tested = False
    over_present = False
    gof_p = None
    selected = "binomial"
    if model in {"auto", "beta_binomial"} and len(np.asarray(alt_counts)) >= 3:
        over_tested = True
        try:
            fit = fit_beta_binomial_rho(alt_counts, depths, mu=err["error_rate"])
            rho = fit["rho"]
            rho_ci = (fit["rho_ci_low"], fit["rho_ci_high"])
            over_present = fit["overdispersion_present"]
            if model == "beta_binomial" or over_present:
                selected = "beta_binomial"
            gof = goodness_of_fit(alt_counts, depths, mu=err["error_rate"], rho=rho)
            gof_p = gof["pvalue"]
            if not over_present:
                warnings.append("No significant overdispersion detected; binomial retained.")
        except ModelError as exc:
            warnings.append(f"Overdispersion not calibrated: {exc}")
            selected = "binomial"
    elif model == "beta_binomial":
        warnings.append("Insufficient control sites for beta-binomial calibration; falling back to binomial.")

    status = "CALIBRATED_DEVELOPMENT" if err["n_sites"] >= 30 else "INSUFFICIENT_CONTROL_DATA"
    if err["n_sites"] < 30:
        warnings.append(
            f"Only {err['n_sites']} control sites; error-rate uncertainty is wide. "
            "Do not treat this profile as a validated LoD."
        )
    return CalibrationProfile(
        mode=mode,
        error_rate=err["error_rate"],
        rho=rho,
        model=selected,
        n_sites=err["n_sites"],
        n_observations=int(err["total_depth"]),
        error_rate_ci=(err["ci_low"], err["ci_high"]),
        rho_ci=rho_ci,
        overdispersion_tested=over_tested,
        overdispersion_present=over_present,
        goodness_of_fit_pvalue=gof_p,
        dataset_id=dataset_id,
        status=status,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Multiple testing and intervals
# ---------------------------------------------------------------------------
def bh_fdr(pvalues) -> np.ndarray:
    """Benjamini-Hochberg adjusted q-values (step-up, monotone)."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1:
        raise ModelError("pvalues must be one-dimensional")
    n = len(p)
    if n == 0:
        return np.array([])
    if np.any(~np.isfinite(p)):
        raise ModelError("pvalues must be finite")
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty(n)
    out[order] = q
    return out


def wilson_interval(successes: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion (Wilson 1927)."""
    if trials <= 0:
        return (0.0, 1.0)
    if successes < 0 or successes > trials:
        raise ModelError("successes must satisfy 0 <= successes <= trials")
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    phat = successes / trials
    denom = 1 + z**2 / trials
    centre = (phat + z**2 / (2 * trials)) / denom
    half = z * math.sqrt(phat * (1 - phat) / trials + z**2 / (4 * trials**2)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def allele_fraction_interval(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Jeffreys equal-tailed interval for an allele fraction.

    Preferred over Wald intervals at low counts and extreme fractions.
    """
    if n <= 0:
        return (0.0, 1.0)
    lo = float(stats.beta.ppf((1 - confidence) / 2, k + 0.5, n - k + 0.5))
    hi = float(stats.beta.ppf(1 - (1 - confidence) / 2, k + 0.5, n - k + 0.5))
    return (max(0.0, lo), min(1.0, hi))


def expected_lod(depth: int, error_rate: float, *, sensitivity: float = 0.95, alpha: float = 0.01) -> float:
    """Analytic planning estimate of the detectable allele fraction.

    Solves for the smallest ``f`` such that the probability of observing at
    least the critical count under H1 is >= ``sensitivity``, where the critical
    count is the smallest ``k`` whose binomial tail under H0 is <= ``alpha``.
    This is a *planning* calculation, not a validated LoD.
    """
    if depth <= 0:
        return float("nan")
    crit = int(stats.binom.isf(alpha, depth, error_rate)) + 1
    crit = min(max(crit, 1), depth)
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        p1 = error_mixed_alt_probability(mid, error_rate)
        power = float(stats.binom.sf(crit - 1, depth, p1))
        if power >= sensitivity:
            hi = mid
        else:
            lo = mid
    return hi


def required_depth_for_fraction(
    fraction: float, error_rate: float, *, sensitivity: float = 0.95, alpha: float = 0.01, max_depth: int = 100000
) -> int:
    """Planning estimate of depth needed to detect ``fraction`` at target power."""
    for depth in (10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000, 2000, 5000, 10000):
        lod = expected_lod(depth, error_rate, sensitivity=sensitivity, alpha=alpha)
        if lod <= fraction:
            return depth
    return max_depth

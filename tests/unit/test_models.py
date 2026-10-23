"""Statistical model tests: binomial, beta-binomial, LRT, FDR, intervals."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pfvarx.errors import ModelError
from pfvarx.models.count_models import (
    allele_fraction_interval,
    bayesian_posterior,
    bh_fdr,
    error_mixed_alt_probability,
    evaluate_hypotheses,
    expected_lod,
    fit_beta_binomial_rho,
    fit_error_rate,
    log_beta_binomial_pmf,
    log_binomial_pmf,
    lrt_pvalue_bootstrap,
    lrt_pvalue_chi2_mixture,
    one_sided_count_pvalue,
    wilson_interval,
)


def test_binomial_pmf_sums_to_one():
    total = sum(math.exp(log_binomial_pmf(k, 50, 0.3)) for k in range(51))
    assert total == pytest.approx(1.0, abs=1e-9)


def test_beta_binomial_reduces_to_binomial_at_rho_zero():
    for k in range(0, 21):
        bb = log_beta_binomial_pmf(k, 20, 0.25, 0.0)
        bi = log_binomial_pmf(k, 20, 0.25)
        assert bb == pytest.approx(bi, abs=1e-9)


def test_error_mixed_probability_bounds():
    assert error_mixed_alt_probability(0.0, 0.01) == pytest.approx(0.01)
    assert error_mixed_alt_probability(1.0, 0.01) == pytest.approx(1.0)
    assert 0.01 < error_mixed_alt_probability(0.5, 0.01) < 0.51


def test_zero_alt_observations_has_negative_lrt():
    h = evaluate_hypotheses(0, 100, error_rate=0.01, alt_fraction=0.2)
    assert h.lrt_statistic <= 1e-9


def test_single_alt_observation_is_weak():
    h = evaluate_hypotheses(1, 100, error_rate=0.01, alt_fraction=0.2)
    assert h.lrt_statistic < 5


def test_strong_alt_has_large_lrt():
    h = evaluate_hypotheses(40, 100, error_rate=0.01, alt_fraction=0.4)
    assert h.lrt_statistic > 50


def test_composite_h1_fits_high_af_instead_of_fixed_detection_boundary():
    k, n, error, boundary = 80, 100, 0.005, 0.05
    composite = evaluate_hypotheses(
        k, n, error_rate=error, alt_fraction=boundary
    )
    fixed_p1 = error_mixed_alt_probability(boundary, error)
    fixed_ll1 = log_binomial_pmf(k, n, fixed_p1)
    fixed_lrt = 2 * (fixed_ll1 - log_binomial_pmf(k, n, error))
    assert composite.alt_fraction_mle == pytest.approx(
        (k / n - error) / (1 - error), abs=1e-10
    )
    assert composite.lrt_statistic > fixed_lrt
    assert composite.p1_mle == pytest.approx(k / n)


def test_beta_binomial_composite_mle_respects_detection_boundary():
    h = evaluate_hypotheses(
        70, 100, error_rate=0.01, alt_fraction=0.10,
        model="beta_binomial", rho=0.02,
    )
    assert h.alt_fraction_mle >= 0.10
    assert h.alt_fraction_mle == pytest.approx((0.70 - 0.01) / 0.99, abs=1e-2)


def test_exact_upper_tail_pvalue_is_monotone_in_alt_count():
    p = [one_sided_count_pvalue(k, 100, error_rate=0.01) for k in (0, 1, 5, 20)]
    assert p == sorted(p, reverse=True)
    assert p[0] == pytest.approx(1.0)


def test_null_simulation_exposes_discovery_selection_and_unconditional_control():
    # A threshold-selected candidate set has a highly enriched raw-p tail.
    # The same exact null p-values still control the unconditional event of a
    # selected, nominally significant site in this one-test family.
    rng = np.random.default_rng(82461)
    n, error, discovery_min, reps, alpha = 60, 0.01, 3, 20_000, 0.05
    counts = rng.binomial(n, error, size=reps)
    selected = counts >= discovery_min
    pvals = np.array([
        one_sided_count_pvalue(int(k), n, error_rate=error) for k in counts
    ])
    nominal_rejections = selected & (pvals <= alpha)
    assert selected.mean() > 0.01
    assert nominal_rejections.sum() / selected.sum() > 0.5
    # Family-wise probability includes all possible outcomes, not only the
    # discovered subset; a Bonferroni bound for this one-hypothesis family is
    # therefore at most alpha. Simulated rate should be near the exact ~2%.
    assert nominal_rejections.mean() < alpha


def test_zero_depth_raises():
    with pytest.raises(ModelError):
        evaluate_hypotheses(0, 0, error_rate=0.01, alt_fraction=0.2)


def test_counts_exceeding_depth_raise():
    with pytest.raises(ModelError):
        evaluate_hypotheses(5, 3, error_rate=0.01, alt_fraction=0.2)


def test_bootstrap_pvalue_in_unit_interval():
    p = lrt_pvalue_bootstrap(40, 100, error_rate=0.01, alt_fraction=0.4, replicates=100)
    assert 0.0 < p <= 1.0


def test_bootstrap_strong_variant_gets_small_p():
    p = lrt_pvalue_bootstrap(40, 100, error_rate=0.01, alt_fraction=0.4, replicates=100)
    assert p < 1e-3


def test_chi2_mixture_pvalue_monotone():
    assert lrt_pvalue_chi2_mixture(0.0) == 1.0
    assert lrt_pvalue_chi2_mixture(10.0) < lrt_pvalue_chi2_mixture(1.0)


def test_bayesian_posterior_increases_with_support():
    low = bayesian_posterior(2, 100, error_rate=0.01, alt_fraction=0.2, prior_h1=0.01)
    high = bayesian_posterior(40, 100, error_rate=0.01, alt_fraction=0.4, prior_h1=0.01)
    assert high["posterior_h1"] > low["posterior_h1"]
    assert 0.0 <= low["posterior_h1"] <= 1.0
    assert 0.0 <= high["posterior_h1"] <= 1.0


def test_bh_fdr_is_monotone_and_bounded():
    p = [0.001, 0.01, 0.02, 0.5, 0.9]
    q = bh_fdr(p)
    assert all(0 <= x <= 1 for x in q)
    assert q[0] <= q[-1]


def test_bh_fdr_matches_manual_for_simple_case():
    q = bh_fdr([0.01, 0.04, 0.03, 0.20])
    assert q[0] == pytest.approx(0.04)


def test_wilson_interval_contains_estimate():
    lo, hi = wilson_interval(5, 10)
    assert lo < 0.5 < hi


def test_wilson_interval_edge_cases():
    lo, hi = wilson_interval(0, 10)
    assert lo == 0.0 and hi > 0
    lo, hi = wilson_interval(10, 10)
    assert hi == pytest.approx(1.0) and lo < 1.0


def test_allele_fraction_interval():
    lo, hi = allele_fraction_interval(10, 100)
    assert lo < 0.1 < hi


def test_fit_error_rate_recovers_low_rate():
    alt = [0, 1, 0, 0, 1, 0, 0, 0, 1, 0]
    depths = [100] * 10
    fit = fit_error_rate(alt, depths)
    assert 0.0 <= fit["error_rate"] < 0.05


def test_fit_rho_in_unit_interval():
    alt = [5, 8, 2, 10, 7]
    depths = [100] * 5
    fit = fit_beta_binomial_rho(alt, depths, mu=0.06)
    assert 0.0 <= fit["rho"] <= 1.0


def test_expected_lod_positive():
    lod = expected_lod(100, 0.01)
    assert lod > 0

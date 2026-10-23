"""Versioned variant decision policy.

This module is the *single* place where statistical evidence is converted into
a decision.  It implements the staged process required by the PfVarX
specification:

1.  Hard gates (depth, callability, mapping, reference concordance, UMI model
    validity) produce an explicit ``NO_CALL``/``ALT_UNCERTAIN`` with named
    filters -- never a silent drop.
2.  Candidate support gates (minimum alternate observations / fraction).
3.  Model evaluation (binomial or beta-binomial) against the null of sequencing
    error only.
4.  Conservative family-wise Bonferroni adjustment over the caller's bounded
    hypothesis universe (candidate filtering does not define the family).
5.  Optional, explicitly-enabled Bayesian posterior (disabled by default).
6.  Decision assignment into the five decision categories plus minority
    fraction bins.

Design rules that are enforced here and asserted in tests:

* A raw likelihood-ratio statistic is **never** reported as a posterior.
* An uncalibrated model can never yield ``ALT_PASS_VALIDATED``.
* ``ALT_PASS_VALIDATED`` requires a ``VALIDATED`` policy status *and* a
  posterior/LRT result that clears the configured threshold.
* Every decision carries a machine-readable list of filter reasons.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pfvarx.calling.evidence import SiteEvidence
from pfvarx.config.schema import PfVarXConfig
from pfvarx.models.count_models import (
    bayesian_posterior,
    evaluate_hypotheses,
    one_sided_count_pvalue,
)

# ---------------------------------------------------------------------------
# Controlled vocabularies (spec section 6.3 / 7.2).  These strings are part of
# the public contract and are asserted in the test-suite.
# ---------------------------------------------------------------------------

POLICY_VERSION = "pfvarx-policy-1.0.0"

FILTER_LOW_EFFECTIVE_DEPTH = "LowEffectiveDepth"
FILTER_LOW_ALT_SUPPORT = "LowAltSupport"
FILTER_UNCALIBRATED_MODEL = "UncalibratedModel"
FILTER_STRAND_BIAS_REVIEW = "StrandBiasReview"
FILTER_HOMOPOLYMER_INDEL_REVIEW = "HomopolymerIndelReview"
FILTER_OUTSIDE_CALLABLE = "OutsideCallableRegion"
FILTER_AMBIGUOUS_MAPPING = "AmbiguousMapping"
FILTER_UMI_MODEL_INVALID = "UMIModelInvalid"
FILTER_REFERENCE_MISMATCH = "ReferenceMismatch"

ALL_FILTERS = (
    FILTER_LOW_EFFECTIVE_DEPTH,
    FILTER_LOW_ALT_SUPPORT,
    FILTER_UNCALIBRATED_MODEL,
    FILTER_STRAND_BIAS_REVIEW,
    FILTER_HOMOPOLYMER_INDEL_REVIEW,
    FILTER_OUTSIDE_CALLABLE,
    FILTER_AMBIGUOUS_MAPPING,
    FILTER_UMI_MODEL_INVALID,
    FILTER_REFERENCE_MISMATCH,
)

DECISION_NO_CALL = "NO_CALL"
DECISION_REFERENCE_SUPPORTED = "REFERENCE_SUPPORTED"
DECISION_ALT_EXPLORATORY = "ALT_DETECTED_EXPLORATORY"
DECISION_ALT_PASS = "ALT_PASS_VALIDATED"
DECISION_ALT_UNCERTAIN = "ALT_UNCERTAIN"

ALL_DECISIONS = (
    DECISION_NO_CALL,
    DECISION_REFERENCE_SUPPORTED,
    DECISION_ALT_EXPLORATORY,
    DECISION_ALT_PASS,
    DECISION_ALT_UNCERTAIN,
)

BIN_BELOW_THRESHOLD = "BELOW_THRESHOLD"
BIN_LOW_FRACTION = "LOW_FRACTION_EXPLORATORY"
BIN_MINORITY = "MINORITY_CANDIDATE"
BIN_INTERMEDIATE = "INTERMEDIATE_FRACTION"
BIN_PREDOMINANT = "PREDOMINANT_ALT"

ALL_BINS = (
    BIN_BELOW_THRESHOLD,
    BIN_LOW_FRACTION,
    BIN_MINORITY,
    BIN_INTERMEDIATE,
    BIN_PREDOMINANT,
)


def minority_bin(alt_fraction: float | None, cfg: PfVarXConfig) -> str:
    """Assign a minority fraction bin using configured boundaries.

    Imported callsets may lack a computed allele fraction; such sites are
    reported as ``BELOW_THRESHOLD`` rather than raising, since the bin is a
    descriptive annotation and not a decision input.
    """
    if alt_fraction is None:
        return BIN_BELOW_THRESHOLD
    rep = cfg.reporting
    if alt_fraction < cfg.candidate.exploratory_min_alt_fraction:
        return BIN_BELOW_THRESHOLD
    if alt_fraction < rep.minority_candidate_fraction:
        return BIN_LOW_FRACTION
    if alt_fraction < rep.intermediate_fraction:
        return BIN_MINORITY
    if alt_fraction < rep.predominant_fraction:
        return BIN_INTERMEDIATE
    return BIN_PREDOMINANT


@dataclass
class SiteCall:
    """The complete, auditable outcome for one candidate site."""

    evidence: SiteEvidence
    decision: str
    filters: list[str] = field(default_factory=list)
    bin: str = BIN_BELOW_THRESHOLD

    model: str = "binomial"
    model_version: str = POLICY_VERSION
    log_likelihood_h0: float | None = None
    log_likelihood_h1: float | None = None
    lrt_statistic: float | None = None
    lrt_pvalue: float | None = None
    lrt_pvalue_adjusted: float | None = None
    posterior_h1: float | None = None
    calibration_profile: str = "UNCALIBRATED"

    policy_status: str = "PROVISIONAL_DEVELOPMENT"
    policy_id: str | None = None
    notes: list[str] = field(default_factory=list)

    @classmethod
    def from_evidence_row(cls, row: dict) -> SiteCall:
        """Reconstruct a SiteCall from a serialised evidence-table row.

        Numeric fields are coerced defensively; missing values become ``None``
        rather than raising, so that partially written tables still render.
        """

        def _f(key: str) -> float | None:
            raw = row.get(key)
            if raw in (None, "", "None", "nan", "NA"):
                return None
            try:
                return float(raw)
            except (TypeError, ValueError):
                return None

        def _i(key: str) -> int | None:
            val = _f(key)
            return None if val is None else int(val)

        def _b(key: str) -> bool | None:
            raw = row.get(key)
            if raw in (None, "", "None"):
                return None
            return str(raw).strip().lower() in {"true", "1", "yes", "y"}

        evidence = SiteEvidence(
            sample_id=str(row.get("sample_id", "sample")),
            chrom=str(row.get("chrom", "")),
            pos=int(_i("pos") or 0),
            ref=str(row.get("ref", "")),
            alt=str(row.get("alt", "")),
            var_type=str(row.get("var_type", "snp")),
            ref_count=int(_i("ref_count") or 0),
            alt_count=int(_i("alt_count") or 0),
            other_count=int(_i("other_count") or 0),
            depth=int(_i("depth") or 0),
            informative_depth=int(_i("informative_depth") or 0),
            alt_fraction=_f("alt_fraction"),
            af_ci_low=_f("af_ci_low"),
            af_ci_high=_f("af_ci_high"),
            mean_base_quality_alt=_f("mean_base_quality_alt"),
            mean_base_quality_ref=_f("mean_base_quality_ref"),
            mean_mapping_quality=_f("mean_mapping_quality"),
            alt_forward=_i("alt_forward"),
            alt_reverse=_i("alt_reverse"),
            ref_forward=_i("ref_forward"),
            ref_reverse=_i("ref_reverse"),
            strand_bias_pvalue=_f("strand_bias_pvalue"),
            homopolymer_length=_i("homopolymer_length"),
            left_aligned=_b("left_aligned"),
            duplicate_fragments_excluded=_i("duplicate_fragments_excluded"),
            low_mapq_excluded=_i("low_mapq_excluded"),
            low_bq_excluded=_i("low_bq_excluded"),
            secondary_excluded=_i("secondary_excluded"),
            sequencing_mode=str(row.get("sequencing_mode", "")),
            umi_status=str(row.get("umi_status", "")),
        )
        filters_raw = row.get("filter_reasons") or row.get("filters") or ""
        if isinstance(filters_raw, str):
            filters = [f for f in filters_raw.replace(";", ",").split(",") if f.strip()]
        else:
            filters = list(filters_raw)
        notes_raw = row.get("policy_notes") or ""
        notes = [n for n in str(notes_raw).split(";") if n.strip()]
        return cls(
            evidence=evidence,
            decision=str(row.get("decision", "NO_CALL")),
            filters=filters,
            bin=str(row.get("minority_bin", BIN_BELOW_THRESHOLD)),
            model=str(row.get("model", "binomial")),
            model_version=str(row.get("model_version", POLICY_VERSION)),
            log_likelihood_h0=_f("log_likelihood_h0"),
            log_likelihood_h1=_f("log_likelihood_h1"),
            lrt_statistic=_f("lrt_statistic"),
            lrt_pvalue=_f("lrt_pvalue"),
            lrt_pvalue_adjusted=_f("lrt_pvalue_adjusted"),
            posterior_h1=_f("posterior_h1"),
            calibration_profile=str(row.get("calibration_profile", "UNCALIBRATED")),
            policy_status=str(row.get("policy_status", "PROVISIONAL_DEVELOPMENT")),
            policy_id=row.get("policy_id"),
            notes=notes,
        )

    @property
    def is_pass(self) -> bool:
        return self.decision == DECISION_ALT_PASS

    @property
    def is_alt(self) -> bool:
        return self.decision in (
            DECISION_ALT_EXPLORATORY,
            DECISION_ALT_PASS,
            DECISION_ALT_UNCERTAIN,
        )

    def as_dict(self) -> dict:
        d = self.evidence.as_dict()
        d.update(
            {
                "decision": self.decision,
                "filter_status": "PASS" if not self.filters else "REVIEW",
                "filter_reasons": ";".join(self.filters) if self.filters else ".",
                "minority_bin": self.bin,
                "model": self.model,
                "model_version": self.model_version,
                "log_likelihood_h0": self.log_likelihood_h0,
                "log_likelihood_h1": self.log_likelihood_h1,
                "lrt_statistic": self.lrt_statistic,
                "lrt_pvalue": self.lrt_pvalue,
                "lrt_pvalue_adjusted": self.lrt_pvalue_adjusted,
                "posterior_h1": self.posterior_h1,
                "calibration_profile": self.calibration_profile,
                "policy_status": self.policy_status,
                "policy_id": self.policy_id,
                "policy_notes": ";".join(self.notes) if self.notes else ".",
            }
        )
        return d


def _callability_threshold(cfg: PfVarXConfig) -> int:
    return cfg.callability.min_depth_for_mode(cfg.run.mode)


def evaluate_site(
    evidence: SiteEvidence,
    cfg: PfVarXConfig,
    *,
    in_callable_region: bool = True,
    reference_matches: bool = True,
    ambiguous_mapping: bool = False,
    umi_model_valid: bool = True,
    compute_bootstrap: bool = False,
    seed: int = 20240101,
) -> SiteCall:
    """Evaluate a single site and return a fully annotated :class:`SiteCall`."""
    call = SiteCall(
        evidence=evidence,
        decision=DECISION_NO_CALL,
        policy_status=cfg.run.policy_status,
        policy_id=cfg.policy_id(),
    )
    call.bin = minority_bin(evidence.alt_fraction, cfg)

    filters: list[str] = []

    # --- Stage 1: hard gates -------------------------------------------------
    if not reference_matches:
        filters.append(FILTER_REFERENCE_MISMATCH)
        call.decision = DECISION_NO_CALL
        call.filters = filters
        call.notes.append("Candidate reference allele disagrees with the pinned FASTA.")
        return call

    if not in_callable_region:
        filters.append(FILTER_OUTSIDE_CALLABLE)

    if ambiguous_mapping:
        filters.append(FILTER_AMBIGUOUS_MAPPING)

    if cfg.run.mode == "tbs_umi" and not umi_model_valid:
        filters.append(FILTER_UMI_MODEL_INVALID)

    effective_depth = evidence.informative_depth
    min_depth = _callability_threshold(cfg)
    if effective_depth < min_depth:
        filters.append(FILTER_LOW_EFFECTIVE_DEPTH)

    hard_gate_failed = any(
        f in filters
        for f in (
            FILTER_OUTSIDE_CALLABLE,
            FILTER_AMBIGUOUS_MAPPING,
            FILTER_UMI_MODEL_INVALID,
            FILTER_LOW_EFFECTIVE_DEPTH,
        )
    )

    # --- Stage 2: candidate support -----------------------------------------
    is_indel = evidence.var_type.upper() in {"INDEL", "INS", "DEL"}
    min_alt = (
        cfg.candidate.indel_min_alt_observations
        if is_indel
        else cfg.candidate.snp_min_alt_observations
    )
    low_support = (
        evidence.alt_count < min_alt
        or evidence.alt_fraction < cfg.statistics.min_alt_fraction_h1
    )
    if low_support:
        filters.append(FILTER_LOW_ALT_SUPPORT)

    # --- Stage 3: model evaluation ------------------------------------------
    error_rate = cfg.statistics.error_rate
    rho = cfg.statistics.beta_binomial_rho or 0.0
    model = cfg.statistics.count_model
    if model == "auto":
        model = "beta_binomial" if rho > 0 else "binomial"

    call.model = model
    call.calibration_profile = "UNCALIBRATED_UNVERIFIED"
    if cfg.run.policy_status == "VALIDATED":
        call.notes.append(
            "Validated certification is unavailable: no verifiable assay-validation/LoD "
            "provenance is consumed by this caller."
        )

    if error_rate is None:
        filters.append(FILTER_UNCALIBRATED_MODEL)
        call.notes.append(
            "No empirically calibrated sequencing error rate; model output is "
            "exploratory and cannot support a validated call."
        )
        error_rate = 0.005  # provisional placeholder for exploratory reporting only

    n = evidence.ref_count + evidence.alt_count
    k = evidence.alt_count
    if n > 0:
        lik = evaluate_hypotheses(
            k=k,
            n=n,
            error_rate=error_rate,
            alt_fraction=cfg.statistics.min_alt_fraction_h1,
            model=model,
            rho=rho,
        )
        call.log_likelihood_h0 = lik.log_l0
        call.log_likelihood_h1 = lik.log_l1
        call.lrt_statistic = lik.lrt_statistic

        # Exact one-sided count tail avoids asymptotic LRT calibration and is
        # independent of the fitted composite-H1 nuisance value. The complete
        # scan universe is supplied to apply_fdr for conservative correction.
        call.lrt_pvalue = one_sided_count_pvalue(
            k=k, n=n, error_rate=error_rate, model=model, rho=rho
        )
        call.notes.append(
            "Exact one-sided count-tail p-value conditional on informative depth and "
            "the configured null; retained in the legacy lrt_pvalue output field."
        )
        if compute_bootstrap:
            call.notes.append(
                "Exact one-sided count-tail calibration supersedes the requested "
                "bootstrap LRT calibration."
            )

        if (
            cfg.statistics.bayesian_enabled
            and cfg.statistics.posterior_threshold_if_validated is not None
        ):
            call.notes.append(
                "Posterior is exploratory and uses the configured point AF boundary; "
                "no composite-AF prior or assay-validation/LoD provenance is available."
            )
            post = bayesian_posterior(
                k=k,
                n=n,
                error_rate=error_rate,
                alt_fraction=cfg.statistics.min_alt_fraction_h1,
                prior_h1=cfg.statistics.bayesian_prior_h1,
                model=model,
                rho=rho,
            )
            call.posterior_h1 = post["posterior_h1"]

    # --- Strand bias / homopolymer review -----------------------------------
    if (
        evidence.strand_bias_pvalue is not None
        and evidence.strand_bias_pvalue < cfg.candidate.strand_bias_alpha
    ):
        filters.append(FILTER_STRAND_BIAS_REVIEW)
    if (
        is_indel
        and evidence.homopolymer_length >= cfg.candidate.homopolymer_review_length
    ):
        filters.append(FILTER_HOMOPOLYMER_INDEL_REVIEW)

    # --- Stage 6: decision assignment ---------------------------------------
    call.filters = filters
    call.decision = _assign_decision(
        evidence=evidence,
        cfg=cfg,
        filters=filters,
        hard_gate_failed=hard_gate_failed,
        low_support=low_support,
        min_depth=min_depth,
        posterior_h1=call.posterior_h1,
    )
    return call


def _assign_decision(
    *,
    evidence: SiteEvidence,
    cfg: PfVarXConfig,
    filters: list[str],
    hard_gate_failed: bool,
    low_support: bool,
    min_depth: int,
    posterior_h1: float | None = None,
) -> str:
    if hard_gate_failed:
        return DECISION_NO_CALL

    if evidence.alt_count == 0:
        if evidence.informative_depth >= min_depth:
            return DECISION_REFERENCE_SUPPORTED
        return DECISION_NO_CALL

    if low_support or evidence.informative_depth < min_depth:
        return DECISION_ALT_UNCERTAIN

    review_only = {
        FILTER_STRAND_BIAS_REVIEW,
        FILTER_HOMOPOLYMER_INDEL_REVIEW,
        FILTER_UNCALIBRATED_MODEL,
    }
    blocking = set(filters) - review_only

    if blocking:
        return DECISION_ALT_UNCERTAIN

    # No assay-validation/LoD provenance object is currently represented or
    # verified. Configuration status, an error-rate value or a posterior must
    # therefore never certify a PASS.
    return DECISION_ALT_EXPLORATORY


def apply_fdr(
    calls: list[SiteCall], cfg: PfVarXConfig, *, family_size: int | None = None
) -> list[SiteCall]:
    """Apply Bonferroni over the complete pre-filter hypothesis universe.

    ``family_size`` must bound every possible tested allele in this sample's
    scan, not just emitted candidates. The correction is per sample, not across
    a multi-sample cohort. When the bound is unavailable, inferential calls fail
    closed instead of applying BH to a discovery-selected subset.
    """
    for c in calls:
        if c.decision == DECISION_ALT_PASS:
            c.decision = DECISION_ALT_EXPLORATORY
            c.notes.append(
                "Validated PASS downgraded: no verifiable assay-validation/LoD "
                "provenance is available to this policy."
            )
    indexed = [(i, c.lrt_pvalue) for i, c in enumerate(calls) if c.lrt_pvalue is not None]
    if family_size is None or family_size < 1:
        for c in calls:
            if c.is_alt:
                c.decision = DECISION_ALT_UNCERTAIN
                c.notes.append(
                    "No complete pre-selection test-family bound supplied; "
                    "candidate-selected p-values are not used for inference."
                )
        return calls
    if not indexed:
        return calls
    for idx, p in indexed:
        calls[idx].lrt_pvalue_adjusted = min(1.0, float(p) * family_size)

    for c in calls:
        if c.lrt_pvalue_adjusted is None or not c.is_alt:
            continue
        threshold = (
            cfg.statistics.sensitivity_fdr_q
            if c.decision == DECISION_ALT_PASS
            else cfg.statistics.exploratory_fdr_q
        )
        if c.lrt_pvalue_adjusted > threshold:
            if c.decision == DECISION_ALT_PASS:
                c.decision = DECISION_ALT_EXPLORATORY
            elif c.decision == DECISION_ALT_EXPLORATORY:
                # An exploratory ALT still requires statistical support at the
                # exploratory FDR.  A site that fails it is uncertain, not a call.
                c.decision = DECISION_ALT_UNCERTAIN
            c.notes.append(
                f"Bonferroni family-adjusted p={c.lrt_pvalue_adjusted:.3g} exceeds "
                f"{threshold:.3g}; call downgraded."
            )
    return calls


def summarise_calls(calls: list[SiteCall]) -> dict:
    """Aggregate decision and filter counts for QC and reporting."""
    decisions: dict[str, int] = dict.fromkeys(ALL_DECISIONS, 0)
    bins: dict[str, int] = dict.fromkeys(ALL_BINS, 0)
    filters: dict[str, int] = dict.fromkeys(ALL_FILTERS, 0)
    for c in calls:
        decisions[c.decision] = decisions.get(c.decision, 0) + 1
        bins[c.bin] = bins.get(c.bin, 0) + 1
        for f in c.filters:
            filters[f] = filters.get(f, 0) + 1
    return {
        "n_sites": len(calls),
        "decisions": decisions,
        "minority_bins": bins,
        "filter_counts": filters,
        "policy_version": POLICY_VERSION,
    }

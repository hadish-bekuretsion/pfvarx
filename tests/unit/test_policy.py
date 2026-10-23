"""Decision-policy tests: gates, decisions, filters, FDR, minority bins."""

from __future__ import annotations

from pfvarx.calling.evidence import SiteEvidence
from pfvarx.calling.policy import (
    ALL_DECISIONS,
    ALL_FILTERS,
    DECISION_ALT_EXPLORATORY,
    DECISION_ALT_PASS,
    DECISION_ALT_UNCERTAIN,
    DECISION_NO_CALL,
    DECISION_REFERENCE_SUPPORTED,
    apply_fdr,
    evaluate_site,
    minority_bin,
    summarise_calls,
)
from tests.conftest import make_config


def ev(alt_count, depth, *, alt_fraction=None, **kw):
    ref_count = depth - alt_count
    af = alt_count / depth if depth else 0.0
    base = {
        "sample_id": "S",
        "chrom": "Pf3D7_01_v3",
        "pos": 100,
        "ref": "A",
        "alt": "T",
        "var_type": "snp",
        "ref_count": ref_count,
        "alt_count": alt_count,
        "depth": depth,
        "informative_depth": depth,
        "alt_fraction": af if alt_fraction is None else alt_fraction,
    }
    base.update(kw)
    return SiteEvidence(**base)


def test_filter_and_decision_vocabularies_are_complete():
    assert set(ALL_FILTERS) >= {
        "LowEffectiveDepth", "LowAltSupport", "UncalibratedModel",
        "StrandBiasReview", "HomopolymerIndelReview", "OutsideCallableRegion",
        "AmbiguousMapping", "UMIModelInvalid", "ReferenceMismatch",
    }
    assert set(ALL_DECISIONS) == {
        "NO_CALL", "REFERENCE_SUPPORTED", "ALT_DETECTED_EXPLORATORY",
        "ALT_PASS_VALIDATED", "ALT_UNCERTAIN",
    }


def test_low_depth_is_no_call():
    cfg = make_config()
    call = evaluate_site(ev(1, 3), cfg)
    assert call.decision == DECISION_NO_CALL
    assert "LowEffectiveDepth" in call.filters


def test_reference_mismatch_is_no_call():
    cfg = make_config()
    call = evaluate_site(ev(50, 100), cfg, reference_matches=False)
    assert call.decision == DECISION_NO_CALL
    assert "ReferenceMismatch" in call.filters


def test_outside_callable_is_no_call():
    cfg = make_config()
    call = evaluate_site(ev(50, 100), cfg, in_callable_region=False)
    assert call.decision == DECISION_NO_CALL
    assert "OutsideCallableRegion" in call.filters


def test_ambiguous_mapping_is_uncertain():
    cfg = make_config()
    call = evaluate_site(ev(50, 100), cfg, ambiguous_mapping=True)
    assert call.decision in (DECISION_NO_CALL, DECISION_ALT_UNCERTAIN)
    assert "AmbiguousMapping" in call.filters


def test_umi_invalid_is_uncertain():
    cfg = make_config(run={"mode": "tbs_umi"}, preprocessing={"umi_design": "D", "umi_min_family_size": 2})
    call = evaluate_site(ev(50, 100, sequencing_mode="tbs_umi"), cfg, umi_model_valid=False)
    assert "UMIModelInvalid" in call.filters


def test_insertions_and_deletions_use_indel_policy_thresholds():
    cfg = make_config(
        candidate={
            "snp_min_alt_observations": 1,
            "indel_min_alt_observations": 3,
            "homopolymer_review_length": 6,
        }
    )
    for variant_type in ("INS", "DEL"):
        call = evaluate_site(
            ev(2, 100, var_type=variant_type, homopolymer_length=6), cfg
        )
        assert "LowAltSupport" in call.filters
        assert "HomopolymerIndelReview" in call.filters


def test_reference_supported_when_no_alt():
    cfg = make_config()
    call = evaluate_site(ev(0, 100), cfg)
    assert call.decision == DECISION_REFERENCE_SUPPORTED


def test_strong_alt_is_exploratory_not_validated_in_provisional_mode():
    cfg = make_config()
    call = evaluate_site(ev(50, 100), cfg)
    assert call.decision == DECISION_ALT_EXPLORATORY
    assert "UncalibratedModel" in call.filters
    # A raw LRT must never be reported as a posterior.
    assert call.posterior_h1 is None


def test_alt_pass_requires_validated_policy():
    cfg = make_config(
        run={"mode": "wgs", "policy_status": "VALIDATED"},
        reference={"assembly_id": "Pf3D7", "fasta": "f", "sha256": "a" * 64},
        statistics={
            "error_rate": 0.005,
            "bayesian_enabled": True,
            "posterior_threshold_if_validated": 0.99,
            "beta_binomial_rho": 0.01,
        },
    )
    call = evaluate_site(ev(50, 100), cfg)
    # User-supplied status/error/Bayes settings are not verified assay/LoD
    # provenance, so the certification state must fail closed.
    assert call.decision != DECISION_ALT_PASS
    assert call.posterior_h1 is not None
    assert call.calibration_profile == "UNCALIBRATED_UNVERIFIED"
    assert any("no verifiable assay-validation/LoD provenance" in n for n in call.notes)


def test_fdr_downgrades_non_significant_exploratory():
    cfg = make_config()
    strong = evaluate_site(ev(50, 100), cfg)
    weak = evaluate_site(ev(2, 100), cfg)
    # Force both to exploratory to test the FDR step in isolation.
    strong.decision = DECISION_ALT_EXPLORATORY
    weak.decision = DECISION_ALT_EXPLORATORY
    strong.lrt_pvalue = 1e-12
    weak.lrt_pvalue = 0.4
    apply_fdr([strong, weak], cfg, family_size=2)
    assert strong.decision == DECISION_ALT_EXPLORATORY
    assert weak.decision == DECISION_ALT_UNCERTAIN


def test_selection_threshold_is_not_used_as_the_test_family():
    cfg = make_config()
    selected = evaluate_site(ev(10, 100), cfg)
    selected.decision = DECISION_ALT_EXPLORATORY
    selected.lrt_pvalue = 0.001
    # Candidate-only BH would leave this apparently strong, but the supplied
    # complete pre-filter universe has 1,000 hypotheses.
    apply_fdr([selected], cfg, family_size=1000)
    assert selected.lrt_pvalue_adjusted == 1.0
    assert selected.decision == DECISION_ALT_UNCERTAIN
    assert any("Bonferroni family-adjusted" in n for n in selected.notes)


def test_missing_family_bound_fails_closed_without_dropping_evidence():
    cfg = make_config()
    selected = evaluate_site(ev(50, 100), cfg)
    assert selected.lrt_pvalue is not None
    apply_fdr([selected], cfg)
    assert selected.decision == DECISION_ALT_UNCERTAIN
    assert selected.lrt_pvalue is not None
    assert any("No complete pre-selection test-family bound" in n for n in selected.notes)


def test_existing_pass_object_is_downgraded_without_validation_provenance():
    cfg = make_config()
    call = evaluate_site(ev(50, 100), cfg)
    call.decision = DECISION_ALT_PASS
    call.posterior_h1 = 1.0
    call.lrt_pvalue = 1e-12
    apply_fdr([call], cfg, family_size=1)
    assert call.decision == DECISION_ALT_EXPLORATORY
    assert any("no verifiable assay-validation/LoD provenance" in n for n in call.notes)


def test_minority_bins_cover_expected_boundaries():
    cfg = make_config()
    assert minority_bin(0.005, cfg) == "BELOW_THRESHOLD"
    assert minority_bin(0.02, cfg) == "LOW_FRACTION_EXPLORATORY"
    assert minority_bin(0.10, cfg) == "MINORITY_CANDIDATE"
    assert minority_bin(0.50, cfg) == "INTERMEDIATE_FRACTION"
    assert minority_bin(0.95, cfg) == "PREDOMINANT_ALT"


def test_summarise_calls_counts_all_decisions():
    cfg = make_config()
    calls = [evaluate_site(ev(50, 100), cfg), evaluate_site(ev(0, 100), cfg)]
    summary = summarise_calls(calls)
    assert summary["n_sites"] == 2
    assert set(summary["decisions"]) == set(ALL_DECISIONS)

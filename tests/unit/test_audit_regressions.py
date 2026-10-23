"""Regression tests for defects found during the full pipeline audit.

Each test pins a specific bug that was identified and repaired so that it
cannot silently reappear.
"""

from __future__ import annotations

import math

import pysam

from pfvarx.calling.evidence import SiteEvidence
from pfvarx.calling.policy import (
    DECISION_ALT_PASS,
    DECISION_NO_CALL,
    SiteCall,
    evaluate_site,
    minority_bin,
)
from pfvarx.calling.writer import write_vcf
from pfvarx.downstream.diversity import diversity_summary
from pfvarx.downstream.mixed_infection import mixed_infection_summary
from pfvarx.models.count_models import (
    evaluate_hypotheses,
    fit_beta_binomial_rho,
    lrt_pvalue_bootstrap,
    lrt_pvalue_exact,
)
from pfvarx.release import RELEASE_STATUS, overall_status
from tests.conftest import make_config


def _evidence(**kw) -> SiteEvidence:
    base = {
        "sample_id": "S1", "chrom": "Pf3D7", "pos": 100, "ref": "A", "alt": "T",
        "var_type": "SNP", "ref_count": 10, "alt_count": 5, "other_count": 0,
        "depth": 15, "informative_depth": 15, "alt_fraction": 1 / 3,
        "af_ci_low": 0.15, "af_ci_high": 0.55, "mean_base_quality_alt": 35,
        "mean_base_quality_ref": 36, "mean_mapping_quality": 55, "alt_forward": 3,
        "alt_reverse": 2, "ref_forward": 5, "ref_reverse": 5,
        "strand_bias_pvalue": 0.5, "homopolymer_length": 0, "left_aligned": True,
    }
    base.update(kw)
    return SiteEvidence(**base)


def _row(**kw) -> dict:
    ev = _evidence()
    row = {k: getattr(ev, k) for k in ev.__dataclass_fields__}
    row.update(
        decision="ALT_UNCERTAIN", filter_reasons="LowAltSupport;StrandBiasReview",
        minority_bin="INTERMEDIATE_FRACTION", model="binomial", model_version="v1",
        log_likelihood_h0=-1.0, log_likelihood_h1=-2.0, lrt_statistic=2.0,
        lrt_pvalue=0.1, lrt_pvalue_adjusted=0.1, posterior_h1=None,
        calibration_profile="UNCALIBRATED", policy_status="PROVISIONAL_DEVELOPMENT",
        policy_notes="",
    )
    row.update(kw)
    return row


def test_overall_release_status_is_the_weakest_scientific_status():
    order = [
        "NOT_IMPLEMENTED",
        "IMPLEMENTED_NOT_VALIDATED",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "VALIDATED_WITH_RESTRICTIONS",
        "RESEARCH_USE_READY",
    ]
    weakest = min(
        (component["scientific_status"] for component in RELEASE_STATUS.values()),
        key=order.index,
    )
    assert overall_status() == weakest


# --- VCF writer -----------------------------------------------------------

def test_vcf_writes_multiple_filters(tmp_path):
    cfg = make_config()
    call = SiteCall.from_evidence_row(_row(policy_id=cfg.policy_id()))
    out = tmp_path / "x.vcf"
    write_vcf([call], cfg, str(out), sample_id="S1")
    with pysam.VariantFile(str(out)) as vf:
        rec = next(iter(vf))
        assert set(rec.filter.keys()) == {"LowAltSupport", "StrandBiasReview"}


def test_vcf_uses_haploid_genotype(tmp_path):
    cfg = make_config()
    call = SiteCall.from_evidence_row(_row(policy_id=cfg.policy_id()))
    out = tmp_path / "x.vcf"
    write_vcf([call], cfg, str(out), sample_id="S1")
    with pysam.VariantFile(str(out)) as vf:
        rec = next(iter(vf))
        assert len(rec.samples["S1"]["GT"]) == 1
        assert rec.samples["S1"]["GT"] == (1,)


# --- policy ---------------------------------------------------------------

def test_minority_bin_handles_none():
    cfg = make_config()
    assert minority_bin(None, cfg) == "BELOW_THRESHOLD"


def test_posterior_gate_blocks_low_posterior():
    cfg = make_config(
        run={"policy_status": "VALIDATED"},
        statistics={
            "bayesian_enabled": True,
            "posterior_threshold_if_validated": 0.99,
            "error_rate": 0.005,
            "min_alt_fraction_h1": 0.05,
        },
    )
    call = SiteCall.from_evidence_row(_row(policy_id=cfg.policy_id()))
    call.posterior_h1 = 0.5
    from pfvarx.calling.policy import _assign_decision

    decision = _assign_decision(
        evidence=call.evidence, cfg=cfg, filters=[], hard_gate_failed=False,
        low_support=False, min_depth=1, posterior_h1=call.posterior_h1,
    )
    assert decision != DECISION_ALT_PASS


def test_reference_mismatch_fails_closed():
    cfg = make_config()
    call = evaluate_site(_evidence(), cfg, reference_matches=False)
    assert call.decision == DECISION_NO_CALL
    assert "ReferenceMismatch" in call.filters


def test_umi_mode_fails_closed_without_consensus():
    cfg = make_config(
        run={"mode": "tbs_umi"},
        preprocessing={"umi_design": "M12S+T"},
    )
    call = evaluate_site(_evidence(sequencing_mode="tbs_umi"), cfg, umi_model_valid=False)
    assert "UMIModelInvalid" in call.filters


# --- statistical models ---------------------------------------------------

def test_overdispersion_lrt_detects_strong_overdispersion():
    import random

    rng = random.Random(7)
    alts, depths = [], []
    for _ in range(200):
        n = 80
        # Beta(2, 20) overdispersion around mu ~ 0.09.
        p = rng.betavariate(2, 20)
        k = sum(1 for _ in range(n) if rng.random() < p)
        alts.append(k)
        depths.append(n)
    fit = fit_beta_binomial_rho(alts, depths, mu=0.09)
    assert fit["rho"] > 0
    assert fit["overdispersion_lr"] > 0
    assert fit["overdispersion_pvalue"] < 0.05


def test_bootstrap_tail_uses_exact_null_not_chi_square():
    # When no simulated null statistic exceeds the observed value, the
    # bootstrap must fall back to the exact conditional null tail rather than
    # silently returning the chi-square value.
    from pfvarx.models.count_models import lrt_pvalue_chi2_mixture

    k, n = 50, 60
    exact = lrt_pvalue_exact(k, n, error_rate=0.005, alt_fraction=0.05)
    boot = lrt_pvalue_bootstrap(k, n, error_rate=0.005, alt_fraction=0.05, replicates=200, seed=1)
    chi = lrt_pvalue_chi2_mixture(
        evaluate_hypotheses(k, n, error_rate=0.005, alt_fraction=0.05).lrt_statistic
    )
    assert math.isclose(exact, boot, rel_tol=1e-12)
    assert not math.isclose(exact, chi, rel_tol=1e-12)


# --- downstream honesty ---------------------------------------------------

def test_diversity_does_not_report_watterson_theta():
    calls = [SiteCall.from_evidence_row(_row())]
    s = diversity_summary(calls, sample_id="S1", callable_bases=1_000_000)
    assert s.watterson_theta is None


def test_mixed_infection_reports_no_numeric_coi():
    calls = [SiteCall.from_evidence_row(_row(minority_bin="INTERMEDIATE_FRACTION")) for _ in range(4)]
    m = mixed_infection_summary(calls, sample_id="S1")
    assert m.coi_point_estimate is None


# --- variant import -------------------------------------------------------

def test_import_handles_missing_dp_header(tmp_path):
    import pysam

    from pfvarx.io.variant_import import import_variants

    cfg = make_config()
    vcf = tmp_path / "in.vcf"
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.contigs.add("Pf3D7")
    header.add_sample("S1")
    header.add_meta("FORMAT", items=[("ID", "GT"), ("Number", "1"), ("Type", "String"), ("Description", "Genotype")])
    with pysam.VariantFile(str(vcf), "w", header=header) as vf:
        rec = vf.new_record(contig="Pf3D7", start=99, stop=100, alleles=("A", "T"))
        rec.samples["S1"]["GT"] = (1,)
        vf.write(rec)
    res = import_variants(str(vcf), cfg, sample_id="S1")
    assert res.n_records == 1
    assert res.calls[0].evidence.var_type == "SNP"


def test_import_marks_indels_uppercase(tmp_path):
    import pysam

    from pfvarx.io.variant_import import import_variants

    cfg = make_config()
    vcf = tmp_path / "in.vcf"
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.contigs.add("Pf3D7")
    header.add_sample("S1")
    header.add_meta("FORMAT", items=[("ID", "GT"), ("Number", "1"), ("Type", "String"), ("Description", "Genotype")])
    with pysam.VariantFile(str(vcf), "w", header=header) as vf:
        rec = vf.new_record(contig="Pf3D7", start=99, stop=101, alleles=("AT", "A"))
        rec.samples["S1"]["GT"] = (1,)
        vf.write(rec)
    res = import_variants(str(vcf), cfg, sample_id="S1")
    assert res.calls[0].evidence.var_type == "INDEL"

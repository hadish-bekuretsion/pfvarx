"""Evidence-extraction tests: pileup counting, scope guards, discovery."""

from __future__ import annotations

import pytest

from pfvarx.calling.evidence import discover_snps, extract_site_evidence
from pfvarx.errors import InputValidationError


def test_extract_snp_evidence_counts_alt(synthetic, cfg):
    v = synthetic.variants[0]  # 50% SNP
    ev = extract_site_evidence(
        synthetic.bam_path,
        synthetic.contig,
        v.pos,
        v.ref,
        v.alt,
        cfg,
        sample_id="S1",
        fasta_path=synthetic.reference_path,
    )
    assert ev.var_type == "SNP"
    assert ev.depth > 0
    assert ev.alt_count > 0
    assert ev.ref_count > 0
    assert 0.0 <= ev.alt_fraction <= 1.0
    assert ev.sequencing_mode == "wgs"
    assert ev.umi_status == "NOT_APPLICABLE"


def test_extract_snp_alt_fraction_near_truth(synthetic, cfg):
    v = synthetic.variants[0]
    ev = extract_site_evidence(
        synthetic.bam_path, synthetic.contig, v.pos, v.ref, v.alt, cfg,
        sample_id="S1", fasta_path=synthetic.reference_path,
    )
    # 50% truth variant should be recovered in a broad band
    assert 0.30 <= ev.alt_fraction <= 0.70


def test_extract_rejects_out_of_scope_allele(synthetic, cfg):
    with pytest.raises(InputValidationError):
        extract_site_evidence(
            synthetic.bam_path,
            synthetic.contig,
            300,
            "A" * 60,
            "T",
            cfg,
            sample_id="S1",
        )


def test_extract_rejects_missing_bam(cfg):
    with pytest.raises(InputValidationError):
        extract_site_evidence("/nonexistent/x.bam", "c", 1, "A", "T", cfg)


def test_discover_snps_finds_truth_variant(synthetic, cfg):
    cands = discover_snps(
        synthetic.bam_path, cfg, fasta_path=synthetic.reference_path
    )
    positions = {c.pos for c in cands}
    for v in synthetic.variants:
        assert v.pos in positions, f"truth SNP at {v.pos} not discovered"


def test_site_evidence_zero_depth_defaults(cfg):
    from pfvarx.calling.evidence import SiteEvidence

    ev = SiteEvidence(sample_id="x", chrom="c", pos=1, ref="A", alt="T", var_type="SNP")
    assert ev.depth == 0
    assert ev.alt_count == 0
    assert ev.alt_fraction == 0.0

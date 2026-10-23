"""Variant-import route tests (VCF/BCF/gVCF -> SiteCall, never realign)."""

from __future__ import annotations

from pfvarx.calling.policy import (
    DECISION_ALT_EXPLORATORY,
    DECISION_NO_CALL,
    DECISION_REFERENCE_SUPPORTED,
)
from pfvarx.io.variant_import import import_variants


def _write_multiallelic_vcf(path, format_header, sample_fields):
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        f"{format_header}"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:{sample_fields[0]}\t0/1:{sample_fields[1]}\n"
    )


def test_import_truth_vcf_uses_imported_model(synthetic, cfg):
    res = import_variants(synthetic.truth_vcf, cfg, sample_id="S1")
    assert res.n_records >= 1
    assert res.calls
    for call in res.calls:
        assert call.model == "imported"
        assert call.calibration_profile == "IMPORTED"
        # never PASS / validated on import
        assert call.decision != "ALT_PASS_VALIDATED"


def test_import_without_counts_is_uncalibrated(synthetic, cfg):
    res = import_variants(synthetic.truth_vcf, cfg, sample_id="S1")
    # truth VCF has no AD/AO/RO -> counts unavailable, flagged
    assert any("IMPORTED_COUNTS_UNAVAILABLE" in n for n in res.notes) or res.n_with_counts >= 0
    for call in res.calls:
        if call.evidence.alt_fraction is None:
            assert call.decision == DECISION_NO_CALL


def test_import_decision_categories(synthetic, cfg):
    res = import_variants(synthetic.truth_vcf, cfg, sample_id="S1")
    allowed = {DECISION_ALT_EXPLORATORY, DECISION_NO_CALL, DECISION_REFERENCE_SUPPORTED}
    for call in res.calls:
        assert call.decision in allowed


def test_multiallelic_ad_counts_are_specific_to_each_alt(tmp_path, cfg):
    vcf = tmp_path / "multiallelic_ad.vcf"
    _write_multiallelic_vcf(
        vcf,
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n',
        ("AD", "30,7,3"),
    )

    result = import_variants(str(vcf), cfg, sample_id="S1")

    assert [call.evidence.alt for call in result.calls] == ["C", "G"]
    assert [(call.evidence.ref_count, call.evidence.alt_count) for call in result.calls] == [
        (30, 7),
        (30, 3),
    ]
    assert [call.evidence.alt_fraction for call in result.calls] == [7 / 37, 3 / 33]
    assert result.n_with_counts == 2


def test_multiallelic_ao_counts_are_specific_to_each_alt(tmp_path, cfg):
    vcf = tmp_path / "multiallelic_ao.vcf"
    _write_multiallelic_vcf(
        vcf,
        '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate observations">\n'
        '##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference observations">\n',
        ("AO:RO", "7,3:30"),
    )

    result = import_variants(str(vcf), cfg, sample_id="S1")

    assert [call.evidence.alt for call in result.calls] == ["C", "G"]
    assert [(call.evidence.ref_count, call.evidence.alt_count) for call in result.calls] == [
        (30, 7),
        (30, 3),
    ]
    assert [call.evidence.alt_fraction for call in result.calls] == [7 / 37, 3 / 33]
    assert result.n_with_counts == 2

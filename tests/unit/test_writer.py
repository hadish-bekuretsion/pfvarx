"""Output-contract tests: evidence schema, FILTER vocabulary, QUAL policy."""

from __future__ import annotations

import csv

import pysam
import pytest

from pfvarx.calling.caller import call_sample
from pfvarx.calling.policy import (
    ALL_FILTERS,
    DECISION_ALT_PASS,
    SiteCall,
    SiteEvidence,
)
from pfvarx.calling.writer import _FILTERS, EVIDENCE_COLUMNS, _qual
from tests.conftest import make_config


def test_evidence_schema_is_stable():
    assert len(EVIDENCE_COLUMNS) == 46
    assert EVIDENCE_COLUMNS[0] == "sample_id"
    assert "posterior_h1" in EVIDENCE_COLUMNS
    assert "policy_id" in EVIDENCE_COLUMNS


def test_filter_vocabulary_matches_policy():
    assert set(_FILTERS) == set(ALL_FILTERS)
    assert len(_FILTERS) == 9


def test_evidence_table_and_vcf_written(synthetic, tmp_path):
    cfg = make_config()
    res = call_sample(
        synthetic.bam_path,
        cfg,
        sample_id="S1",
        fasta_path=synthetic.reference_path,
        compute_bootstrap=False,
        write_outputs=True,
        output_dir=str(tmp_path / "out"),
    )
    assert res.evidence_path and res.vcf_path

    with open(res.evidence_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
    assert header == EVIDENCE_COLUMNS

    with pysam.VariantFile(res.vcf_path) as vf:
        assert "PVX_DECISION" in vf.header.info
        assert "PVX_POLICY_ID" in vf.header.info
        for rec in vf:
            # provisional run => never a QUAL value
            assert rec.qual is None
            for f in rec.filter:
                assert f in ALL_FILTERS or f == "PASS"


def test_qual_only_with_validated_posterior():
    ev = SiteEvidence(sample_id="s", chrom="c", pos=1, ref="A", alt="T", var_type="SNP")
    call = SiteCall(evidence=ev, decision=DECISION_ALT_PASS, posterior_h1=0.99)
    provisional = make_config()
    assert _qual(call, provisional) is None

    validated = make_config(
        run={"policy_status": "VALIDATED"},
        reference={"assembly_id": "Pf3D7_01_v3"},
        statistics={"error_rate": 0.005},
    )
    q = _qual(call, validated)
    assert q is not None and q > 0

    call.posterior_h1 = None
    assert _qual(call, validated) is None


def test_vcf_never_encodes_lrt_as_qual(synthetic, tmp_path):
    cfg = make_config()
    res = call_sample(
        synthetic.bam_path, cfg, sample_id="S1",
        fasta_path=synthetic.reference_path, compute_bootstrap=False,
        write_outputs=True, output_dir=str(tmp_path / "o"),
    )
    with pysam.VariantFile(res.vcf_path) as vf:
        for rec in vf:
            if rec.qual is not None:
                pytest.fail("provisional record must not carry QUAL")

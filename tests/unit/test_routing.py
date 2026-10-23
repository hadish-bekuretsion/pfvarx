"""Input routing tests: FASTQ, BAM, VCF and error paths."""

from __future__ import annotations

import pytest

from pfvarx.errors import RoutingError
from pfvarx.io.routing import (
    ROUTE_BAM,
    ROUTE_FASTQ,
    ROUTE_VARIANT,
    plan_run,
    route_sample,
    summarise_plan,
)
from pfvarx.io.validators import SampleRecord
from tests.conftest import make_config


def _touch(path, text="x"):
    path.write_text(text)
    return str(path)


def test_bam_route_reuses_alignment(tmp_path):
    bam = _touch(tmp_path / "s.bam", "not-a-real-bam-but-nonempty")
    s = SampleRecord(sample_id="S1", input_type="bam", bam=bam)
    d = route_sample(s, make_config())
    assert d.route == ROUTE_BAM
    assert "realignment" in " ".join(d.reasons).lower()


def test_fastq_route_requires_alignment(tmp_path):
    r1 = _touch(tmp_path / "s_R1.fastq", "@r1\nACGT\n+\nIIII\n")
    r2 = _touch(tmp_path / "s_R2.fastq", "@r2\nACGT\n+\nIIII\n")
    s = SampleRecord(sample_id="S1", input_type="fastq", fastq_1=r1, fastq_2=r2)
    d = route_sample(s, make_config())
    assert d.route == ROUTE_FASTQ
    assert "preprocess" in d.stages


def test_vcf_route_imports(tmp_path):
    vcf = _touch(tmp_path / "s.vcf", "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    s = SampleRecord(sample_id="S1", input_type="vcf", vcf=vcf)
    d = route_sample(s, make_config())
    assert d.route == ROUTE_VARIANT
    assert "normalize" in d.stages


def test_declared_type_mismatch_raises(tmp_path):
    bam = _touch(tmp_path / "s.bam", "x")
    s = SampleRecord(sample_id="S1", input_type="vcf", bam=bam)
    with pytest.raises(RoutingError):
        route_sample(s, make_config())


def test_missing_evidence_raises():
    s = SampleRecord(sample_id="S1", input_type="bam")
    with pytest.raises(RoutingError):
        route_sample(s, make_config())


def test_missing_file_raises(tmp_path):
    s = SampleRecord(sample_id="S1", input_type="bam", bam=str(tmp_path / "nope.bam"))
    with pytest.raises(RoutingError):
        route_sample(s, make_config())


def test_bam_rejected_when_not_accepted(tmp_path):
    bam = _touch(tmp_path / "s.bam", "x")
    cfg = make_config(inputs={"samplesheet": "s.tsv", "accept_existing_bam": False})
    s = SampleRecord(sample_id="S1", input_type="bam", bam=bam)
    with pytest.raises(RoutingError):
        route_sample(s, cfg)


def test_plan_run_detects_duplicate_ids(tmp_path):
    bam = _touch(tmp_path / "s.bam", "x")
    s1 = SampleRecord(sample_id="S1", input_type="bam", bam=bam)
    s2 = SampleRecord(sample_id="S1", input_type="bam", bam=bam)
    with pytest.raises(RoutingError):
        plan_run([s1, s2], make_config())


def test_summary_counts_routes(tmp_path):
    bam = _touch(tmp_path / "s.bam", "x")
    vcf = _touch(tmp_path / "v.vcf", "##fileformat=VCFv4.3\n")
    samples = [
        SampleRecord(sample_id="A", input_type="bam", bam=bam),
        SampleRecord(sample_id="B", input_type="vcf", vcf=vcf),
    ]
    summary = summarise_plan(plan_run(samples, make_config()))
    assert summary["n_samples"] == 2
    assert summary["reuses_alignment"] == 1
    assert summary["imports_variants"] == 1

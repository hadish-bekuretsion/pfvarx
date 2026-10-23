"""End-to-end and workflow integration tests.

These exercise the real code paths (no mocks): synthetic data generation,
calling, evidence/VCF writing, import, QC and reporting, plus the atomic-write
and idempotent-rerun guarantees that a resumable workflow depends on.
"""

from __future__ import annotations

import csv
import json
import os

from pfvarx.calling.caller import call_sample
from pfvarx.cli import main
from pfvarx.io.variant_import import import_variants
from pfvarx.qc.integrator import run_qc, sample_qc
from pfvarx.reporting.exporters import export_pdf, export_workbook
from tests.conftest import make_config


def _write_config(
    tmp_path, reference_path: str = "reference.fa", checksum: str = "REQUIRED_SHA256"
) -> str:
    cfg = tmp_path / "run.yaml"
    cfg.write_text(
        "run:\n"
        "  mode: wgs\n"
        "  intent: research\n"
        "  policy_status: PROVISIONAL_DEVELOPMENT\n"
        "reference:\n"
        "  assembly_id: Pf3D7\n"
        f"  fasta: {reference_path}\n"
        f"  sha256: {checksum}\n"
        "inputs:\n"
        "  samplesheet: samples.tsv\n"
    )
    return str(cfg)


def test_cli_demo_end_to_end(tmp_path, capsys):
    cfg = _write_config(tmp_path)
    out = tmp_path / "demo"
    rc = main(["demo", "--config", cfg, "--out", str(out), "--coverage", "60"])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    assert "RESEARCH USE ONLY" in payload["banner"]
    assert os.path.exists(payload["vcf"])
    evidence = out / "results" / "DEMO" / "DEMO.evidence.tsv"
    assert evidence.exists()
    with open(evidence, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    decisions = {r["decision"] for r in rows}
    assert "ALT_DETECTED_EXPLORATORY" in decisions


def test_atomic_write_leaves_no_tmp(synthetic, tmp_path):
    cfg = make_config()
    res = call_sample(
        synthetic.bam_path, cfg, sample_id="S1",
        fasta_path=synthetic.reference_path, compute_bootstrap=False,
        write_outputs=True, output_dir=str(tmp_path / "o"),
    )
    assert res.vcf_path and os.path.exists(res.vcf_path)
    leftovers = list((tmp_path / "o").rglob("*.tmp"))
    assert leftovers == []


def test_rerun_is_idempotent(synthetic, tmp_path):
    cfg = make_config()
    out = str(tmp_path / "o")
    r1 = call_sample(
        synthetic.bam_path, cfg, sample_id="S1",
        fasta_path=synthetic.reference_path, compute_bootstrap=False,
        write_outputs=True, output_dir=out,
    )
    with open(r1.evidence_path) as fh:
        first = fh.read()
    r2 = call_sample(
        synthetic.bam_path, cfg, sample_id="S1",
        fasta_path=synthetic.reference_path, compute_bootstrap=False,
        write_outputs=True, output_dir=out,
    )
    with open(r2.evidence_path) as fh:
        second = fh.read()
    assert first == second


def test_full_chain_calling_qc_report(synthetic, tmp_path):
    cfg = make_config()
    res = call_sample(
        synthetic.bam_path, cfg, sample_id="S1",
        fasta_path=synthetic.reference_path, compute_bootstrap=False,
        write_outputs=True, output_dir=str(tmp_path / "o"),
    )
    metrics = sample_qc("S1", res.calls, cfg)
    qc = run_qc(metrics, cfg, n_samples=1)
    assert qc["summary"]["overall_status"] in {"PASS", "WARN", "FAIL"}
    pdf = export_pdf(
        res.calls, qc, str(tmp_path / "report.pdf"), cfg=cfg, sample_id="S1"
    )
    xlsx = export_workbook(res.calls, qc, str(tmp_path / "report.xlsx"))
    assert os.path.exists(pdf) and os.path.getsize(pdf) > 0
    assert os.path.exists(xlsx)


def test_import_route_writes_outputs(synthetic, tmp_path):
    from pfvarx.calling.writer import write_evidence_table, write_vcf

    cfg = make_config()
    res = import_variants(synthetic.truth_vcf, cfg, sample_id="S1")
    vcf = write_vcf(res.calls, cfg, str(tmp_path / "imp.vcf"), sample_id="S1")
    ev = write_evidence_table(res.calls, str(tmp_path / "imp.tsv"))
    assert os.path.exists(vcf)
    assert os.path.exists(ev)
    assert all(c.model == "imported" for c in res.calls)


def test_cli_validate_routing_plan(tmp_path, capsys, synthetic):
    from pfvarx.io.checksums import sha256_file

    cfg = _write_config(
        tmp_path,
        reference_path=synthetic.reference_path,
        checksum=sha256_file(synthetic.reference_path),
    )
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample_id\tinput_type\tassay\tfastq_1\tfastq_2\tbam\tvcf\n"
        f"S1\tBAM\tWGS\t\t\t{synthetic.bam_path}\t\n"
    )
    out = tmp_path / "routing.json"
    rc = main([
        "validate", "--config", cfg,
        "--samplesheet", str(sheet),
        "--out", str(out),
    ])
    assert rc == 0
    plan = json.loads(out.read_text())
    assert plan["samples"][0]["sample_id"] == "S1"
    assert plan["samples"][0]["route"] == "alignment_reuse"

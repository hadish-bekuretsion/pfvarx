"""CLI smoke tests: parser wiring and JSON-emitting subcommands."""

from __future__ import annotations

import json

from pfvarx.cli import build_parser, main


def test_parser_has_all_subcommands():
    parser = build_parser()
    actions = [a for a in parser._actions if a.dest == "command"]
    assert actions
    names = set(actions[0].choices)
    for expected in {
        "init", "validate", "run", "call", "calibrate", "benchmark", "report",
        "parameters", "status", "demo", "import", "annotate",
        "downstream", "relatedness", "qc",
    }:
        assert expected in names


def test_status_emits_json(capsys):
    rc = main(["status"])
    assert rc == 0
    out = capsys.readouterr().out
    payload = json.loads(out)
    assert "banner" in payload
    assert "RESEARCH USE ONLY" in payload["banner"]
    assert "overall" in payload
    assert "components" in payload


def test_annotate_writes_header_when_no_annotations(tmp_path, monkeypatch, capsys):
    import csv
    from types import SimpleNamespace

    import pfvarx.annotation.annotator as annotator
    import pfvarx.cli as cli

    evidence = tmp_path / "empty.evidence.tsv"
    evidence.write_text("sample_id\tchrom\tpos\tref\talt\n")
    output = tmp_path / "annotated.tsv"
    from pfvarx.io.checksums import sha256_file

    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n")
    cfg = SimpleNamespace(
        reference=SimpleNamespace(fasta=str(reference), sha256=sha256_file(reference))
    )
    monkeypatch.setattr(cli, "_load_cfg", lambda args: cfg)
    monkeypatch.setattr(annotator, "annotate_calls", lambda *args, **kwargs: [])

    rc = main([
        "annotate", "--config", "unused.yaml", "--evidence", str(evidence), "--out", str(output)
    ])

    assert rc == 0
    capsys.readouterr()
    with output.open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert "annotation_engine" in rows[0]
    assert "resistance_evidence" in rows[0]


def test_run_output_guard_rejects_existing_summary(tmp_path):
    from types import SimpleNamespace

    import pytest

    from pfvarx.cli import _check_run_output_collisions
    from pfvarx.errors import ConfigError

    (tmp_path / "run_summary.json").write_text("prior run")
    with pytest.raises(ConfigError, match="Refusing to overwrite"):
        _check_run_output_collisions(
            tmp_path, [SimpleNamespace(sample_id="S1")], ["diversity"],
            evidence_reuse=False, overwrite=False,
        )


def test_run_output_guard_allows_snakefile_staged_evidence(tmp_path):
    from types import SimpleNamespace

    from pfvarx.cli import _check_run_output_collisions

    (tmp_path / "routing_plan.json").write_text("upstream plan")
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    (sample_dir / "S1.evidence.tsv").write_text("upstream evidence")
    (tmp_path / "logs").mkdir()
    (tmp_path / "logs" / "call.log").write_text("upstream log")

    _check_run_output_collisions(
        tmp_path, [SimpleNamespace(sample_id="S1")], ["diversity"],
        evidence_reuse=True, overwrite=False,
    )


def test_reference_checksum_must_match_exact_fasta(tmp_path):
    from types import SimpleNamespace

    import pytest

    from pfvarx.cli import _verify_reference_checksum
    from pfvarx.errors import ConfigError
    from pfvarx.io.checksums import sha256_file

    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTACGT\n")
    cfg = SimpleNamespace(
        reference=SimpleNamespace(fasta=str(reference), sha256=sha256_file(reference))
    )
    _verify_reference_checksum(cfg)

    cfg.reference.sha256 = "0" * 64
    with pytest.raises(ConfigError, match="checksum does not match"):
        _verify_reference_checksum(cfg)

    cfg.reference.sha256 = "REQUIRED_SHA256"
    with pytest.raises(ConfigError, match="setup placeholder"):
        _verify_reference_checksum(cfg)


def test_import_cli_rejects_reference_digest_mismatch_before_writing(tmp_path):
    import yaml

    from pfvarx.cli import main

    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTACGT\n")
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({
        "run": {"mode": "wgs"},
        "reference": {
            "assembly_id": "test",
            "fasta": str(reference),
            "sha256": "0" * 64,
        },
        "inputs": {"samplesheet": str(tmp_path / "samples.tsv")},
    }, sort_keys=False))
    vcf = tmp_path / "input.vcf"
    vcf.write_text("not read because preflight must fail first\n")
    out = tmp_path / "out"

    rc = main([
        "import", "--config", str(config), "--vcf", str(vcf),
        "--sample-id", "S1", "--out", str(out),
    ])

    assert rc == 3
    assert not out.exists()


def test_catalogue_override_must_exist(tmp_path):
    from pfvarx.cli import _catalogue_available
    from tests.conftest import make_config

    cfg = make_config(downstream={"resistance_catalogue": str(tmp_path / "missing.yaml")})
    assert not _catalogue_available(cfg)


def test_drug_resistance_gff_is_checked_before_run(tmp_path):
    from types import SimpleNamespace

    import pytest

    from pfvarx.cli import _resolve_objectives
    from pfvarx.errors import ConfigError
    from tests.conftest import make_config

    with pytest.raises(ConfigError, match="annotation GFF3"):
        _resolve_objectives(
            SimpleNamespace(objectives=["drug_resistance"], metadata=None),
            make_config(),
            [],
        )


def test_parameters_emits_json(capsys):
    rc = main(["parameters"])
    assert rc == 0
    records = json.loads(capsys.readouterr().out)
    assert isinstance(records, list)
    assert len(records) > 0
    assert "unit" in records[0]


def test_parameters_mode_filter(capsys):
    rc = main(["parameters", "--mode", "wgs"])
    assert rc == 0
    records = json.loads(capsys.readouterr().out)
    assert all("wgs" in r.get("modes", []) for r in records)


def test_missing_config_returns_error(tmp_path):
    rc = main([
        "validate",
        "--config", str(tmp_path / "nope.yaml"),
        "--samplesheet", str(tmp_path / "samples.tsv"),
    ])
    assert rc in (2, 3)


def test_init_creates_only_configured_input_directories(tmp_path, capsys):
    import yaml

    data = tmp_path / "inputs"
    results = tmp_path / "final-results"
    logs = tmp_path / "run-logs"
    work = tmp_path / "run-work"
    config = tmp_path / "run.yaml"
    config.write_text(yaml.safe_dump({
        "run": {"mode": "wgs"},
        "reference": {
            "assembly_id": "Pf3D7",
            "fasta": str(data / "reference" / "Pf3D7.fa"),
            "sha256": "REQUIRED_SHA256",
        },
        "inputs": {
            "root": str(data),
            "raw_dir": str(data / "raw"),
            "prepared_dir": str(data / "prepared"),
            "metadata_dir": str(data / "metadata"),
            "samplesheet": str(data / "prepared" / "samples.tsv"),
        },
        "outputs": {
            "root": str(results),
            "logs_dir": str(logs),
            "work_dir": str(work),
        },
    }, sort_keys=False))

    rc = main(["init", "--config", str(config), "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert rc == 0
    assert payload["stage"] == "input_preparation"
    assert payload["output_directories_created"] is False
    assert (data / "raw").is_dir()
    assert (data / "prepared").is_dir()
    assert (data / "metadata").is_dir()
    assert (data / "reference").is_dir()
    assert not results.exists()
    assert not logs.exists()
    assert not work.exists()


def test_run_evidence_dir_executes_only_selected_objectives(tmp_path, synthetic):
    """`run --evidence-dir` must reuse evidence and run only the selected objectives."""
    import yaml

    from pfvarx.calling.caller import call_sample
    from pfvarx.config.schema import PfVarXConfig
    from pfvarx.io.checksums import sha256_file

    cfg_data = {
        "run": {"mode": "wgs", "intent": "research", "policy_status": "PROVISIONAL_DEVELOPMENT"},
        "reference": {
            "assembly_id": "Pf3D7",
            "fasta": synthetic.reference_path,
            "sha256": sha256_file(synthetic.reference_path),
        },
        "inputs": {"samplesheet": str(tmp_path / "samples.tsv")},
        "outputs": {
            "root": str(tmp_path / "out"),
            "logs_dir": str(tmp_path / "run-logs"),
            "work_dir": str(tmp_path / "run-work"),
        },
    }
    cfg = PfVarXConfig.model_validate(cfg_data)
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg_data, sort_keys=False))

    (tmp_path / "samples.tsv").write_text(
        "sample_id\tinput_type\tbam\nS1\tbam\t" + synthetic.bam_path + "\n"
    )

    # Produce evidence once (as the Snakemake calling stage would).
    evidence_root = tmp_path / "evidence"
    call_sample(
        synthetic.bam_path,
        cfg,
        sample_id="S1",
        fasta_path=synthetic.reference_path,
        output_dir=str(evidence_root / "variants"),
    )
    assert (evidence_root / "variants" / "S1" / "S1.evidence.tsv").exists()

    rc = main([
        "run",
        "--config", str(cfg_path),
        "--samplesheet", str(tmp_path / "samples.tsv"),
        "--objectives", "diversity",
        "--evidence-dir", str(evidence_root),
        "--out", str(tmp_path / "out"),
    ])
    assert rc == 0

    manifest = json.loads((tmp_path / "out" / "run_manifest.json").read_text())
    assert {o["objective_id"] for o in manifest["outcomes"]} == {"diversity"}
    assert (tmp_path / "out" / "objectives" / "diversity").is_dir()
    assert not (tmp_path / "out" / "objectives" / "relatedness").exists()
    assert (tmp_path / "run-logs" / "run.log").is_file()
    assert (tmp_path / "run-work").is_dir()


def test_run_fastq_without_workflow_exits_nonzero(tmp_path, synthetic):
    import yaml

    from pfvarx.io.checksums import sha256_file

    cfg_data = {
        "run": {"mode": "wgs", "intent": "research", "policy_status": "PROVISIONAL_DEVELOPMENT"},
        "reference": {
            "assembly_id": "Pf3D7",
            "fasta": synthetic.reference_path,
            "sha256": sha256_file(synthetic.reference_path),
        },
        "inputs": {"samplesheet": str(tmp_path / "samples.tsv")},
        "objectives": ["diversity"],
        "outputs": {"root": str(tmp_path / "out")},
    }
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg_data, sort_keys=False))
    samplesheet = tmp_path / "samples.tsv"
    samplesheet.write_text(
        "sample_id\tinput_type\tfastq_1\tfastq_2\n"
        f"S1\tfastq\t{synthetic.fastq_r1}\t{synthetic.fastq_r2}\n"
    )

    rc = main([
        "run", "--config", str(cfg_path), "--samplesheet", str(samplesheet),
        "--objectives", "diversity", "--out", str(tmp_path / "out"),
    ])

    assert rc != 0


def test_demo_writes_run_manifest(tmp_path):
    """The 60-second demo must leave a manifest so the dashboard is populated."""
    import yaml

    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        yaml.safe_dump(
            {
                "run": {
                    "mode": "wgs",
                    "intent": "research",
                    "policy_status": "PROVISIONAL_DEVELOPMENT",
                },
                "reference": {
                    "assembly_id": "Pf3D7",
                    "fasta": "reference.fa",
                    "sha256": "REQUIRED_SHA256",
                },
                "inputs": {"samplesheet": "samples.tsv"},
                "outputs": {"root": str(tmp_path / "demo" / "results")},
            },
            sort_keys=False,
        )
    )
    rc = main(
        ["demo", "--config", str(cfg_path), "--out", str(tmp_path / "demo"), "--coverage", "30"]
    )
    assert rc == 0
    manifest_path = tmp_path / "demo" / "results" / "run_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text())
    assert set(manifest["objectives"]) == {"diversity", "mixed_infection", "coi"}
    assert manifest["status"] == "completed"

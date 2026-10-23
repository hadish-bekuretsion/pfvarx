"""End-to-end smoke test for objective-driven execution.

Builds two synthetic BAM samples, then runs ``pfvarx run`` three times with
different objective selections and asserts that only the requested artefacts
appear on disk.

Run from the repository root::

    python scripts/smoke_objectives.py /tmp/pfvarx-smoke
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

from pfvarx.synthetic import Variant, simulate_dataset


def build_dataset(root: Path) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    import pysam

    ref = root / "reference.fa"
    datasets = []
    for i, (seed, cov) in enumerate([(3, 100), (5, 90)]):
        ds = simulate_dataset(
            str(root / f"sample{i + 1}"),
            length=1500,
            coverage=cov,
            variants=None,
            seed=seed,
            write_fastq=False,
        )
        datasets.append(ds)
    # one shared reference for both samples
    shutil.copy(datasets[0].reference_path, ref)
    return {"reference": str(ref), "datasets": datasets}


def write_samplesheet(root: Path, datasets: list) -> Path:
    path = root / "samples.tsv"
    lines = ["sample_id\tinput_type\tbam\tmetadata_json"]
    for i, ds in enumerate(datasets):
        meta = {"country": "Kenya" if i == 0 else "Uganda", "collection_date": f"2023-0{i + 1}-15"}
        lines.append(
            f"sample{i + 1}\tbam\t{ds.bam_path}\t{json.dumps(meta)}"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def write_config(root: Path, reference: str, samplesheet: Path) -> Path:
    cfg = {
        "schema_version": 1,
        "project": "PfVarX-smoke",
        "run": {"mode": "wgs", "intent": "research", "policy_status": "PROVISIONAL_DEVELOPMENT"},
        "reference": {
            "assembly_id": "Pf3D7",
            "fasta": reference,
            "sha256": "REQUIRED_SHA256",
        },
        "inputs": {"samplesheet": str(samplesheet)},
        "outputs": {"root": str(root / "results")},
    }
    path = root / "config.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def run(cfg: Path, objectives: str, out: Path) -> subprocess.CompletedProcess:
    cmd = [
        sys.executable, "-m", "pfvarx.cli", "run",
        "--config", str(cfg),
        "--samplesheet", str(cfg.parent / "samples.tsv"),
        "--objectives", objectives,
        "--out", str(out),
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/pfvarx-smoke")
    if root.exists():
        shutil.rmtree(root)
    built = build_dataset(root)
    sheet = write_samplesheet(root, built["datasets"])
    cfg = write_config(root, built["reference"], sheet)

    failures: list[str] = []

    # --- Case 1: single-sample objectives only -----------------------------
    out1 = root / "results_single"
    proc = run(cfg, "diversity,mixed_infection,coi", out1)
    print("CASE 1 stdout:", proc.stdout[-1500:])
    if proc.returncode != 0:
        failures.append(f"case1 exit={proc.returncode} stderr={proc.stderr[-800:]}")
    manifest = json.loads((out1 / "run_manifest.json").read_text())
    got = {o["objective_id"] for o in manifest["outcomes"]}
    if got != {"diversity", "mixed_infection", "coi"}:
        failures.append(f"case1 objectives {got}")
    for oid in ("diversity", "mixed_infection", "coi"):
        for sample in ("sample1", "sample2"):
            if not (out1 / "objectives" / oid / sample).exists():
                failures.append(f"case1 missing dir {oid}/{sample}")
    for absent in ("relatedness", "phylogenetics", "population_structure", "spatial_temporal", "drug_resistance"):
        if (out1 / "objectives" / absent).exists():
            failures.append(f"case1 unexpected dir {absent}")
    if (out1 / "alignments").exists():
        failures.append("case1 alignments/ created for a BAM route")

    # --- Case 2: multi-sample objectives -----------------------------------
    out2 = root / "results_multi"
    proc = run(cfg, "relatedness,phylogenetics", out2)
    print("CASE 2 stdout:", proc.stdout[-1500:])
    if proc.returncode != 0:
        failures.append(f"case2 exit={proc.returncode} stderr={proc.stderr[-800:]}")
    manifest2 = json.loads((out2 / "run_manifest.json").read_text())
    got2 = {o["objective_id"] for o in manifest2["outcomes"]}
    if got2 != {"relatedness", "phylogenetics"}:
        failures.append(f"case2 objectives {got2}")
    if not (out2 / "objectives" / "phylogenetics" / "tree.nwk").exists():
        failures.append("case2 no tree.nwk")
    if (out2 / "objectives" / "diversity").exists():
        failures.append("case2 unexpected diversity dir")

    # --- Case 3: metadata objective without metadata must fail fast ---------
    cfg_nometa = root / "config_nometa.yaml"
    data = yaml.safe_load(cfg.read_text())
    cfg_nometa.write_text(yaml.safe_dump(data, sort_keys=False))
    sheet_nometa = root / "samples_nometa.tsv"
    sheet_nometa.write_text(
        "\n".join(
            ["sample_id\tinput_type\tbam"]
            + [f"sample{i + 1}\tbam\t{ds.bam_path}" for i, ds in enumerate(built["datasets"])]
        )
        + "\n"
    )
    out3 = root / "results_nometa"
    cmd = [
        sys.executable, "-m", "pfvarx.cli", "run",
        "--config", str(cfg_nometa),
        "--samplesheet", str(sheet_nometa),
        "--objectives", "population_structure",
        "--out", str(out3),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        failures.append("case3 should have failed without metadata")
    if "metadata" not in (proc.stderr + proc.stdout).lower():
        failures.append("case3 error message does not mention metadata")
    if out3.exists() and (out3 / "objectives").exists():
        failures.append("case3 created output before validating prerequisites")

    # --- Case 4: unknown objective -----------------------------------------
    proc = run(cfg, "not_a_real_objective", root / "results_bad")
    if proc.returncode == 0:
        failures.append("case4 accepted an unknown objective")

    print()
    if failures:
        print("FAILURES:")
        for f in failures:
            print("  -", f)
        return 1
    print("ALL SMOKE CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

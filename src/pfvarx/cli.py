"""PfVarX command-line interface.

Subcommands
-----------
    pfvarx validate     Validate inputs, route samples, write a routing plan.
    pfvarx run          Execute the full pipeline for a configured run.
    pfvarx call         Call variants for a single BAM/CRAM.
    pfvarx calibrate    Fit error rate / over-dispersion from control sites.
    pfvarx benchmark    Run the benchmark / ablation framework.
    pfvarx report       Export CSV/XLSX/PNG/PDF or launch the dashboard.
    pfvarx parameters   Print or export the parameter dictionary.
    pfvarx status       Print the release-status matrix.
    pfvarx demo         Generate a synthetic dataset and run an end-to-end call.

Every command prints the mandatory research-use banner on start-up.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from pfvarx import __version__
from pfvarx.errors import PfVarXError
from pfvarx.release import overall_status, status_banner

BANNER = status_banner()


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _emit(payload: Any, *, as_json: bool) -> None:
    if as_json:
        print(json.dumps(payload, indent=2, default=str))
    else:
        print(payload)


def _load_cfg(args: argparse.Namespace):
    from pfvarx.config.loader import load_config

    return load_config(
        args.config,
        profile=getattr(args, "profile", None),
        validate_paths=getattr(args, "check_paths", False),
    )


def _verify_reference_checksum(cfg, reference_path: str | Path | None = None) -> None:
    """Require the configured digest to match the FASTA used for analysis."""
    from pfvarx.errors import ConfigError
    from pfvarx.io.checksums import sha256_file

    expected = str(cfg.reference.sha256)
    path = Path(reference_path or cfg.reference.fasta)
    if not path.is_absolute():
        path = Path.cwd() / path
    if expected.upper().startswith("REQUIRED_"):
        raise ConfigError(
            "Reference checksum is still a setup placeholder.",
            cause=f"reference.fasta={path}",
            action=(
                "Run sha256sum on the exact reference FASTA, set reference.sha256 to the "
                "64-hex digest, then validate/run again."
            ),
        )
    if not path.is_file():
        raise ConfigError(
            f"Reference FASTA not found: {path}",
            action="Correct reference.fasta or the direct-call --fasta path.",
        )
    try:
        actual = sha256_file(path)
    except OSError as exc:
        raise ConfigError(
            f"Could not read reference FASTA for checksum verification: {path}",
            cause=str(exc),
            action="Check file permissions and storage availability.",
        ) from exc
    if actual.lower() != expected.lower():
        raise ConfigError(
            "Reference FASTA checksum does not match reference.sha256.",
            cause=f"configured={expected}; observed={actual}; file={path}",
            action="Verify the assembly and update reference.sha256 from the exact FASTA.",
        )


def _write_json(path: str | Path, payload: Any) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, default=str))
    return out


# --------------------------------------------------------------------------- #
# Commands
# --------------------------------------------------------------------------- #
def _selected_objectives(args: argparse.Namespace, cfg) -> list[str]:
    """CLI ``--objectives`` wins over the configuration file."""
    raw = getattr(args, "objectives", None)
    if raw:
        if isinstance(raw, str):
            return [part.strip() for part in raw.split(",") if part.strip()]
        return list(raw)
    return list(getattr(cfg, "objectives", []) or [])


def _catalogue_available(cfg) -> bool:
    from pfvarx.annotation.catalogue import DEFAULT_CATALOGUE

    path = getattr(cfg.downstream, "resistance_catalogue", None)
    return Path(path).is_file() if path else DEFAULT_CATALOGUE.is_file()


def _resolve_objectives(args, cfg, samples):
    """Resolve objectives with every prerequisite checked *before* execution."""
    from pfvarx.execution import load_metadata
    from pfvarx.objectives import resolve

    selected = _selected_objectives(args, cfg)
    metadata = None
    if selected:
        metadata = load_metadata(samples, cfg, explicit_path=getattr(args, "metadata", None))
    resolved = resolve(
        selected,
        n_samples=len(samples),
        have_metadata=metadata is not None,
        have_catalogue=_catalogue_available(cfg),
    )
    if "drug_resistance" in resolved.selected:
        gff = getattr(cfg.annotation, "gff", None)
        if not gff or not Path(gff).is_file():
            from pfvarx.errors import ConfigError

            raise ConfigError(
                "The drug-resistance objective requires an existing annotation GFF3 file.",
                cause=f"annotation.gff is unset or unavailable: {gff!r}",
                action="Set annotation.gff to the matching reference's GFF3 gene model.",
            )
    return resolved, metadata


def _check_run_output_collisions(out_root: Path, decisions, objectives, *,
                                 evidence_reuse: bool, overwrite: bool,
                                 logs_dir: Path | None = None) -> None:
    """Reject known output-file collisions before creating or writing anything.

    Existing directories alone are not collisions: the Snakemake workflow
    stages evidence and logs under the same results root before invoking this
    command.  Existing files this command will replace are collisions unless
    the user explicitly enables ``outputs.overwrite``.
    """
    if overwrite:
        return
    candidates = [
        out_root / "config" / "effective_config.json",
        out_root / "objective_plan.json",
        out_root / "run_manifest.json",
        out_root / "run_summary.json",
        (logs_dir or out_root / "logs") / "run.log",
    ]
    if not evidence_reuse:
        candidates.append(out_root / "routing_plan.json")
        for decision in decisions:
            sid = decision.sample_id
            candidates.extend(
                (
                    out_root / sid / f"{sid}.evidence.tsv",
                    out_root / sid / f"{sid}.pfvarx.vcf",
                    out_root / "variants" / sid / f"{sid}.evidence.tsv",
                    out_root / "variants" / sid / f"{sid}.pfvarx.vcf",
                )
            )
    collisions = [str(path) for path in candidates if path.is_file()]
    for objective in objectives:
        objective_root = out_root / "objectives" / objective
        if objective_root.is_dir():
            collisions.extend(str(path) for path in objective_root.rglob("*") if path.is_file())
    if collisions:
        from pfvarx.errors import ConfigError

        raise ConfigError(
            "Refusing to overwrite existing run outputs: " + ", ".join(sorted(set(collisions))),
            action="Choose a new outputs.root or set outputs.overwrite=true.",
        )


def cmd_objectives(args: argparse.Namespace) -> int:
    from pfvarx.objectives import describe, format_table

    if args.json:
        _emit({"banner": BANNER, "objectives": describe()}, as_json=True)
    else:
        print(format_table())
    return 0


def cmd_init(args: argparse.Namespace) -> int:
    """Create configured input directories only; never create run outputs."""
    from pfvarx.workspace import create_input_directories, plan_input_directories

    cfg = _load_cfg(args)
    planned = plan_input_directories(cfg)
    created = create_input_directories(cfg)
    created_set = set(created)
    existing = [path for path in planned if path not in created_set]
    payload = {
        "stage": "input_preparation",
        "created": [str(path) for path in created],
        "already_present": [str(path) for path in existing],
        "output_directories_created": False,
    }
    if args.json:
        _emit(payload, as_json=True)
    else:
        print("Input directories ready:")
        for path in planned:
            marker = "created" if path in created_set else "already present"
            print(f"  {marker}: {path}")
        print("No results, logs, or work directories were created.")
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    """Validate inputs, route samples and preview the objective plan.

    Nothing heavy runs here.  The command exists so that every problem a run
    could hit (missing reference, unknown objective, absent metadata, absent
    catalogue, too few samples for a multi-sample objective) is reported before
    a single byte of output is written.
    """
    from pfvarx.io.routing import route_all, summarise_plan
    from pfvarx.io.validators import read_samplesheet
    from pfvarx.objectives import OBJECTIVES
    from pfvarx.workspace import plan_layout, resolve_output_directories

    cfg = _load_cfg(args)
    samples, sheet_result = read_samplesheet(args.samplesheet, base_dir=args.base_dir)
    if not sheet_result.ok:
        _emit(sheet_result.as_dict(), as_json=args.json)
        return 2

    decisions, routing_errors = route_all(samples, cfg)
    if routing_errors:
        payload = {
            "stage": "routing",
            "n_errors": len(routing_errors),
            "errors": routing_errors,
            "corrective_action": (
                "Fix every item above, then re-run 'pfvarx validate'. "
                "No computation has been performed and no output was written."
            ),
        }
        if args.json:
            _emit(payload, as_json=True)
        else:
            print(f"\n{len(routing_errors)} sample(s) could not be routed:\n")
            for err in routing_errors:
                print(f"  - {err.get('message', err)}")
                if err.get("corrective_action"):
                    print(f"    fix: {err['corrective_action']}")
            print("\nNo computation has been performed and no output was written.")
        return 2
    _verify_reference_checksum(cfg)
    summary = summarise_plan(decisions)
    summary["samplesheet"] = sheet_result.as_dict()

    resolved, metadata = _resolve_objectives(args, cfg, samples)
    summary["objectives"] = resolved.as_dict()
    if metadata is not None:
        summary["metadata"] = {
            "n_rows": int(len(metadata)),
            "columns": [str(c) for c in metadata.columns],
        }
    out_root = Path(args.out or cfg.outputs.root)
    logs_dir, work_dir = resolve_output_directories(cfg, out_root)
    layout = plan_layout(
        out_root,
        objectives=resolved.selected,
        routes=[d.route for d in decisions],
        needs_annotation=any(o in resolved.selected for o in ("drug_resistance",)),
        logs_dir=logs_dir,
        work_dir=work_dir,
    )
    summary["workspace"] = layout.as_dict()

    if not args.json:
        _print_preflight(summary, resolved, layout, OBJECTIVES)

    if args.out:
        _write_json(args.out, summary)
    if args.json:
        _emit(summary, as_json=True)
    return 0


def _print_preflight(summary, resolved, layout, objectives_registry) -> None:
    print("")
    print("=" * 78)
    print("PfVarX pre-execution summary")
    print("=" * 78)
    print(f"samples        : {summary['n_samples']}")
    for route, count in summary["routes"].items():
        print(f"  route {route:<22}: {count}")
    print(f"objectives     : {', '.join(resolved.selected) or '(none selected)'}")
    print(f"required stages: {', '.join(resolved.stages) or '(none)'}")
    if resolved.selected:
        print("objective detail:")
        for oid in resolved.selected:
            obj = objectives_registry[oid]
            print(f"  - {oid:<22} {obj.summary}")
            print(f"    artefacts: {', '.join(obj.artefacts)}")
    print("directories that will be created:")
    for d in layout.directories:
        print(f"  + {d:<32} {layout.rationale.get(d, '')}")
    if resolved.warnings:
        print("warnings:")
        for w in resolved.warnings:
            print(f"  ! {w}")
    print("=" * 78)
    print("")


def cmd_run(args: argparse.Namespace) -> int:
    """Execute the pipeline for exactly the objectives the user selected.

    Ten-step sequence:
      1 read config  2 validate  3 determine start stage  4 prepare inputs
      5 build execution plan  6 create directories  7 execute
      8 record logs/status  9 generate relevant results  10 final summary
    """
    from pfvarx.config.loader import write_effective_config
    from pfvarx.execution import (
        RunLogger,
        execute_objectives,
        format_summary,
        load_calls,
        write_manifest,
    )
    from pfvarx.io.routing import route_all, summarise_plan
    from pfvarx.io.validators import read_samplesheet
    from pfvarx.workspace import prepare, resolve_output_directories

    # 1. read config
    cfg = _load_cfg(args)

    # 2. validate inputs
    samples, sheet_result = read_samplesheet(args.samplesheet, base_dir=args.base_dir)
    if not sheet_result.ok:
        _emit(sheet_result.as_dict(), as_json=args.json)
        return 2
    decisions, routing_errors = route_all(samples, cfg)
    if routing_errors:
        for err in routing_errors:
            print(f"[routing] {err.get('message', err)}", file=sys.stderr)
            if err.get("corrective_action"):
                print(f"          fix: {err['corrective_action']}", file=sys.stderr)
        print(
            "No computation has been performed and no output was written.",
            file=sys.stderr,
        )
        return 2
    _verify_reference_checksum(cfg)

    # 3-4. objectives + prerequisite checks before any heavy work
    resolved, metadata = _resolve_objectives(args, cfg, samples)

    out_root = Path(args.out or cfg.outputs.root)
    logs_dir, work_dir = resolve_output_directories(cfg, out_root)
    evidence_reuse = bool(getattr(args, "evidence_dir", None))
    evidence_files: dict[str, Path] = {}
    if evidence_reuse:
        evidence_root = Path(args.evidence_dir)
        missing = []
        for decision in decisions:
            sid = decision.sample_id
            candidates = (
                evidence_root / "variants" / sid / f"{sid}.evidence.tsv",
                evidence_root / sid / f"{sid}.evidence.tsv",
            )
            found = next((path for path in candidates if path.is_file()), None)
            if found is None:
                missing.append(sid)
            else:
                evidence_files[sid] = found
        if missing:
            print(
                "Evidence reuse failed for: " + ", ".join(missing)
                + ". Run the calling stage first or drop --evidence-dir.",
                file=sys.stderr,
            )
            return 2
    _check_run_output_collisions(
        out_root,
        decisions,
        resolved.selected,
        evidence_reuse=evidence_reuse,
        overwrite=cfg.outputs.overwrite,
        logs_dir=logs_dir,
    )

    # 6. create only the directories this run needs
    layout, created, ws_warnings = prepare(
        out_root,
        objectives=resolved.selected,
        routes=[d.route for d in decisions],
        overwrite=cfg.outputs.overwrite,
        logs_dir=logs_dir,
        work_dir=work_dir,
    )
    logger = RunLogger(logs_dir / "run.log")
    for warning in ws_warnings:
        logger.warning(warning)

    # 5. execution plan is recorded before execution for auditability
    write_effective_config(cfg, out_root / "config" / "effective_config.json")
    routing_plan_path = out_root / "routing_plan.json"
    if not (evidence_reuse and routing_plan_path.is_file()):
        _write_json(routing_plan_path, summarise_plan(decisions))
    _write_json(out_root / "objective_plan.json", resolved.as_dict())

    # 7. execute per-sample evidence production
    calls_by_sample: dict[str, list] = {}
    sample_results: list[dict] = []

    if evidence_reuse:
        # Reuse mode: evidence tables were produced by an upstream orchestrator
        # (the Snakemake workflow).  Only the selected objectives are executed.
        for decision in decisions:
            sid = decision.sample_id
            found = evidence_files[sid]
            calls_by_sample[sid] = load_calls(found)
            sample_results.append({
                "sample_id": sid, "route": decision.route,
                "status": "evidence_reused", "evidence": str(found),
            })
            logger.info(f"sample {sid}: reused evidence {found}")
    else:
        for decision in decisions:
            sample = next(s for s in samples if s.sample_id == decision.sample_id)
            logger.info(f"sample {sample.sample_id}: route={decision.route}")
            evidence_path: str | None = None
            if decision.route == "variant_import":
                from pfvarx.calling.writer import write_evidence_table, write_vcf
                from pfvarx.io.variant_import import import_variants

                imp = import_variants(sample.vcf, cfg, sample_id=sample.sample_id)
                sample_dir = out_root / "variants" / sample.sample_id
                sample_dir.mkdir(parents=True, exist_ok=True)
                vcf_path = write_vcf(
                    imp.calls, cfg,
                    str(sample_dir / f"{sample.sample_id}.pfvarx.vcf"),
                    sample_id=sample.sample_id,
                )
                evidence_path = write_evidence_table(
                    imp.calls, str(sample_dir / f"{sample.sample_id}.evidence.tsv")
                )
                sample_results.append({
                    "sample_id": sample.sample_id, "route": decision.route,
                    "status": "imported", "vcf": vcf_path, "evidence": evidence_path,
                    "n_records": imp.n_records, "n_with_counts": imp.n_with_counts,
                })
            elif decision.route == "fastq_alignment":
                logger.warning(
                    f"sample {sample.sample_id}: FASTQ alignment is performed by the Snakemake "
                    "workflow; the in-process driver cannot align reads"
                )
                sample_results.append({
                    "sample_id": sample.sample_id, "route": decision.route,
                    "status": "not_executed_requires_alignment",
                    "note": "FASTQ alignment is performed only by the Snakemake workflow.",
                })
            else:
                from pfvarx.calling.caller import call_sample

                res = call_sample(
                    sample.bam,
                    cfg,
                    sample_id=sample.sample_id,
                    fasta_path=cfg.reference.fasta,
                    callable_bed=cfg.reference.callable_mask_bed,
                    output_dir=str(out_root / "variants"),
                )
                evidence_path = res.evidence_path
                sample_results.append({
                    "sample_id": sample.sample_id, "route": decision.route,
                    "status": "called", "vcf": res.vcf_path,
                    "evidence": res.evidence_path, "summary": res.summary,
                })
            if evidence_path and Path(evidence_path).exists():
                calls_by_sample[sample.sample_id] = load_calls(evidence_path)

    # 7-9. execute the selected objectives only
    manifest = execute_objectives(
        root=out_root,
        resolved=resolved,
        calls_by_sample=calls_by_sample,
        samples=samples,
        cfg=cfg,
        metadata=metadata,
        logger=logger,
        run_id=getattr(cfg.run, "run_id", None),
    )
    manifest.directories = list(layout.directories)
    manifest.routes = summarise_plan(decisions)["routes"]

    # 10. record everything and summarise
    write_manifest(manifest, out_root)
    _write_json(
        out_root / "run_summary.json",
        {"banner": BANNER, "samples": sample_results, "manifest": manifest.as_dict()},
    )

    if args.json:
        _emit({"banner": BANNER, "samples": sample_results, "manifest": manifest.as_dict()}, as_json=True)
    else:
        print(format_summary(manifest))

    return 0 if manifest.status == "completed" else 1



def cmd_call(args: argparse.Namespace) -> int:
    from pfvarx.calling.caller import call_sample

    cfg = _load_cfg(args)
    _verify_reference_checksum(cfg, args.fasta or cfg.reference.fasta)
    res = call_sample(
        args.bam,
        cfg,
        sample_id=args.sample_id,
        fasta_path=args.fasta or cfg.reference.fasta,
        callable_bed=args.callable_bed or cfg.reference.callable_mask_bed,
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        max_candidates=args.max_candidates,
        output_dir=args.out,
    )
    payload = {
        "banner": BANNER,
        "sample_id": res.sample_id,
        "vcf": res.vcf_path,
        "evidence": res.evidence_path,
        "summary": res.summary,
    }
    _emit(payload, as_json=args.json)
    return 0


def cmd_calibrate(args: argparse.Namespace) -> int:
    from pfvarx.models.count_models import (
        calibrate_from_controls,
        fit_beta_binomial_rho,
        fit_error_rate,
    )

    alt_counts, depths = [], []
    cfg = _load_cfg(args)
    path = Path(args.controls)
    if not path.exists():
        print(f"Controls table not found: {path}", file=sys.stderr)
        return 2
    header = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if header is None:
            header = parts
            continue
        row = dict(zip(header, parts, strict=False))
        depths.append(int(row.get("depth", 0)))
        alt_counts.append(int(row.get("alt_count", 0)))
    if not depths:
        print("No control rows found", file=sys.stderr)
        return 2
    error_fit = fit_error_rate(alt_counts, depths)
    rho_fit = fit_beta_binomial_rho(alt_counts, depths, mu=max(error_fit["error_rate"], 1e-6))
    profile = calibrate_from_controls(
        alt_counts, depths, mode=cfg.run.mode, dataset_id=str(path)
    )
    payload = {
        "banner": BANNER,
        "n_control_sites": len(depths),
        "error_rate": error_fit,
        "overdispersion": rho_fit,
        "calibration_profile": profile.as_dict() if hasattr(profile, "as_dict") else profile,
    }
    if args.out:
        _write_json(args.out, payload)
    _emit(payload, as_json=args.json)
    return 0


def cmd_benchmark(args: argparse.Namespace) -> int:
    from pfvarx.benchmark.runner import run_benchmark

    cfg = _load_cfg(args)
    result = run_benchmark(
        cfg,
        out_dir=args.out,
        replicates=args.replicates,
        coverage=args.coverage,
        ablation=not args.no_ablation,
        seed=args.seed,
    )
    _emit(result, as_json=args.json)
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    from pfvarx.calling.policy import SiteCall
    from pfvarx.reporting import export_figures, export_pdf, export_workbook

    if args.what == "dashboard":
        import subprocess

        cmd = [
            sys.executable, "-m", "streamlit", "run",
            str(Path(__file__).with_name("reporting") / "dashboard.py"),
            "--", "--results", args.results,
        ]
        return subprocess.call(cmd)

    # Objective-aware export: when the results tree carries a run manifest, the
    # exporters emit only the objectives that were selected for that run.
    if args.what in {"workbook", "pdf", "all"}:
        results_root = Path(args.results)
        if (results_root / "run_manifest.json").exists():
            from pfvarx.reporting import report as objective_report

            produced: dict[str, object] = {}
            if args.what in {"workbook", "all"}:
                produced["workbook"] = objective_report.export_run_workbook(results_root)
                produced["csv"] = objective_report.export_run_csv(results_root)
            if args.what in {"pdf", "all"}:
                produced["pdf"] = objective_report.export_run_pdf(results_root)
            produced["figures"] = objective_report.export_run_figures(results_root)
            _emit(
                {"banner": BANNER, "mode": "objective_aware", "produced": produced},
                as_json=args.json,
            )
            return 0

    # Rebuild calls from one or more evidence tables for export.
    import csv

    evidence_paths: list[Path] = []
    if args.evidence:
        evidence_paths = [Path(args.evidence)]
    elif args.results:
        evidence_paths = sorted(Path(args.results).rglob("*.evidence.tsv"))
    if not evidence_paths:
        print("No evidence table found; pass --evidence or --results", file=sys.stderr)
        return 2
    missing = [str(p) for p in evidence_paths if not p.exists()]
    if missing:
        print(f"Evidence table not found: {', '.join(missing)}", file=sys.stderr)
        return 2
    rows: list[dict] = []
    for path in evidence_paths:
        with path.open() as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    calls = [SiteCall.from_evidence_row(r) for r in rows]
    cfg = _load_cfg(args)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    figs = export_figures(calls, str(out_dir / "figures"))
    if args.what in {"workbook", "all"}:
        export_workbook(calls, {}, str(out_dir / "report.xlsx"), metadata={"sample_id": args.sample_id})
    if args.what in {"pdf", "all"}:
        export_pdf(calls, {}, str(out_dir / "report.pdf"), cfg=cfg,
                   sample_id=args.sample_id, figures=figs)
    _emit({"banner": BANNER, "figures": figs, "out_dir": str(out_dir)}, as_json=args.json)
    return 0


def cmd_parameters(args: argparse.Namespace) -> int:
    from pfvarx.config.loader import export_parameter_dictionary
    from pfvarx.config.parameters import registry_records

    if args.out:
        export_parameter_dictionary(args.out)
        print(f"Parameter dictionary written to {args.out}")
        return 0
    records = registry_records()
    if args.mode:
        records = [r for r in records if args.mode in r.get("modes", [])]
    _emit(records, as_json=True)
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    from pfvarx.release import RELEASE_STATUS

    components = {
        name: {
            "software_status": rec["software_status"],
            "scientific_status": rec["scientific_status"],
            "note": rec["note"],
        }
        for name, rec in RELEASE_STATUS.items()
    }
    _emit(
        {"banner": BANNER, "overall": overall_status(), "components": components},
        as_json=True,
    )
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    from pfvarx.calling.caller import call_sample
    from pfvarx.synthetic import Variant, simulate_dataset

    out = Path(args.out)
    data_dir = out / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    # Build the reference first so REF alleles are guaranteed to match.
    from pfvarx.synthetic import write_reference

    ref_path = str(data_dir / "reference.fa")
    write_reference(ref_path, "Pf3D7_01_v3", 2000, seed=args.seed)
    import pysam

    with pysam.FastaFile(ref_path) as fa:
        seq = fa.fetch("Pf3D7_01_v3").upper()

    def _alt(base: str) -> str:
        return {"A": "T", "T": "A", "C": "G", "G": "C"}[base]

    specs = [(300, 0.50), (700, 0.20), (1200, 0.08)]
    variants = [Variant(pos, seq[pos - 1], _alt(seq[pos - 1]), frac) for pos, frac in specs]
    ds = simulate_dataset(
        str(data_dir),
        length=2000,
        coverage=args.coverage,
        variants=variants,
        seed=args.seed,
    )
    cfg = _load_cfg(args)
    res = call_sample(
        ds.bam_path, cfg, sample_id="DEMO",
        fasta_path=ds.reference_path, output_dir=str(out / "results"),
    )

    # Run a small, fixed objective set so the demo produces a real manifest and
    # a dashboard that is populated out of the box.
    from pfvarx.execution import RunLogger, execute_objectives, load_calls, write_manifest
    from pfvarx.objectives import resolve

    demo_objectives = ["diversity", "mixed_infection", "coi"]
    resolved = resolve(
        demo_objectives, n_samples=1, have_metadata=False, have_catalogue=False
    )
    manifest = execute_objectives(
        root=out / "results",
        resolved=resolved,
        calls_by_sample={"DEMO": load_calls(res.evidence_path)},
        samples=[],
        cfg=cfg,
        logger=RunLogger(out / "results" / "logs" / "run.log"),
    )
    write_manifest(manifest, out / "results")

    payload = {"banner": BANNER, "dataset": str(out / "data"),
               "vcf": res.vcf_path, "summary": res.summary,
               "objectives": demo_objectives,
               "manifest": str(out / "results" / "run_manifest.json")}
    _emit(payload, as_json=True)
    return 0


def _load_calls(evidence_path: str) -> list:
    import csv

    from pfvarx.calling.policy import SiteCall

    with open(evidence_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [SiteCall.from_evidence_row(r) for r in rows]


def cmd_import(args: argparse.Namespace) -> int:
    from pfvarx.calling.writer import write_evidence_table, write_vcf
    from pfvarx.io.variant_import import import_variants

    cfg = _load_cfg(args)
    _verify_reference_checksum(cfg)
    res = import_variants(args.vcf, cfg, sample_id=args.sample_id, max_records=args.max_records)
    out = Path(args.out) / args.sample_id
    out.mkdir(parents=True, exist_ok=True)
    vcf_path = write_vcf(res.calls, cfg, str(out / f"{args.sample_id}.pfvarx.vcf"), sample_id=args.sample_id)
    ev_path = write_evidence_table(res.calls, str(out / f"{args.sample_id}.evidence.tsv"))
    payload = {
        "banner": BANNER,
        "sample_id": args.sample_id,
        "n_records": res.n_records,
        "n_with_counts": res.n_with_counts,
        "vcf": vcf_path,
        "evidence": ev_path,
        "notes": res.notes,
    }
    _emit(payload, as_json=args.json)
    return 0


def cmd_annotate(args: argparse.Namespace) -> int:
    from pfvarx.annotation.annotator import annotate_calls
    from pfvarx.annotation.catalogue import load_catalogue

    cfg = _load_cfg(args)
    _verify_reference_checksum(cfg, args.fasta or cfg.reference.fasta)
    calls = _load_calls(args.evidence)
    catalogue = None
    if args.catalogue:
        catalogue = load_catalogue(args.catalogue)
    annotated = annotate_calls(
        calls, cfg, gff_path=args.gff, fasta_path=args.fasta or cfg.reference.fasta,
        catalogue=catalogue,
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    import csv

    rows = [a.as_dict() for a in annotated]
    if rows:
        fieldnames = list(rows[0].keys())
    else:
        from pfvarx.calling.evidence import SiteEvidence

        fieldnames = list(SiteEvidence.__dataclass_fields__) + [
            "decision", "filter_status", "filter_reasons", "minority_bin", "model",
            "model_version", "log_likelihood_h0", "log_likelihood_h1", "lrt_statistic",
            "lrt_pvalue", "lrt_pvalue_adjusted", "posterior_h1", "calibration_profile",
            "policy_status", "policy_id", "policy_notes", "annotation_engine", "consequence",
            "gene_id", "gene_name", "amino_acid_change", "resistance_drugs",
            "resistance_markers", "resistance_evidence",
        ]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    _emit({"banner": BANNER, "n_annotated": len(annotated), "out": str(out)}, as_json=args.json)
    return 0


def cmd_downstream(args: argparse.Namespace) -> int:
    from pfvarx.downstream.diversity import diversity_summary
    from pfvarx.downstream.mixed_infection import mixed_infection_summary

    calls = _load_calls(args.evidence)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    div = diversity_summary(calls, sample_id=args.sample_id, callable_bases=args.callable_bases)
    mix = mixed_infection_summary(calls, sample_id=args.sample_id)
    _write_json(out / f"{args.sample_id}.diversity.json", div.as_dict())
    _write_json(out / f"{args.sample_id}.mixed_infection.json", mix.as_dict())
    _emit({"banner": BANNER, "diversity": div.as_dict(), "mixed_infection": mix.as_dict()}, as_json=args.json)
    return 0


def cmd_relatedness(args: argparse.Namespace) -> int:
    from pfvarx.downstream.phylogenetics import newick
    from pfvarx.downstream.relatedness import relatedness_matrix

    root = Path(args.results)
    samples: dict[str, list] = {}
    for ev in sorted(root.rglob("*.evidence.tsv")):
        sid = ev.name.replace(".evidence.tsv", "")
        samples[sid] = _load_calls(str(ev))
    if len(samples) < 2:
        _emit({"banner": BANNER, "note": "Need at least two samples for relatedness"}, as_json=args.json)
        return 0
    results = relatedness_matrix(samples, min_depth=args.min_depth)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    _write_json(out / "relatedness.json", [r.as_dict() for r in results])
    tree = newick(sorted(samples), results)
    (out / "tree.nwk").write_text(tree)
    _emit({"banner": BANNER, "n_pairs": len(results), "tree": tree}, as_json=args.json)
    return 0


def cmd_qc(args: argparse.Namespace) -> int:
    from pfvarx.qc.integrator import run_qc, sample_qc

    cfg = _load_cfg(args)
    root = Path(args.results)
    all_metrics = []
    for ev in sorted(root.rglob("*.evidence.tsv")):
        sid = ev.name.replace(".evidence.tsv", "")
        calls = _load_calls(str(ev))
        all_metrics.extend(sample_qc(sid, calls, cfg))
    summary = run_qc(all_metrics, cfg, n_samples=len(list(root.rglob("*.evidence.tsv"))))
    _write_json(args.out, summary)
    _emit(summary, as_json=args.json)
    return 0


# --------------------------------------------------------------------------- #
# Parser
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="pfvarx", description="PfVarX — P. falciparum SNP/INDEL platform")
    p.add_argument("--version", action="version", version=f"pfvarx {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    def add_common(sp):
        sp.add_argument("--config", required=True, help="Path to run configuration YAML")
        sp.add_argument("--profile", default=None, help="Mode profile to merge under the config")
        sp.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
        sp.add_argument("--check-paths", action="store_true", help="Verify referenced paths exist")

    init = sub.add_parser("init", help="Create configured input directories only (Stage 1)")
    init.add_argument("--config", required=True, help="Path to run configuration YAML")
    init.add_argument("--profile", default=None, help="Mode profile to merge under the config")
    init.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    init.set_defaults(func=cmd_init)

    v = sub.add_parser("validate", help="Validate inputs, route samples and preview the objective plan")
    add_common(v)
    v.add_argument("--samplesheet", required=True)
    v.add_argument("--base-dir", default=None)
    v.add_argument("--out", default=None, help="Write the pre-execution plan JSON here")
    v.add_argument(
        "--objectives",
        default=None,
        help="Comma-separated analysis objectives (e.g. drug_resistance,diversity). "
        "Overrides the config file. Run 'pfvarx objectives' to list them.",
    )
    v.add_argument("--metadata", default=None, help="TSV/CSV sample metadata table")
    v.set_defaults(func=cmd_validate)

    r = sub.add_parser("run", help="Execute the pipeline for the selected objectives")
    add_common(r)
    r.add_argument("--samplesheet", required=True)
    r.add_argument("--base-dir", default=None)
    r.add_argument("--out", default=None)
    r.add_argument(
        "--objectives",
        default=None,
        help="Comma-separated analysis objectives (e.g. drug_resistance,diversity). "
        "Overrides the config file. Run 'pfvarx objectives' to list them.",
    )
    r.add_argument("--metadata", default=None, help="TSV/CSV sample metadata table")
    r.add_argument(
        "--evidence-dir",
        default=None,
        help="Reuse evidence tables already present under this results tree and run "
        "only the selected objectives. Used by the Snakemake workflow after it has "
        "produced evidence for FASTQ/BAM/VCF inputs.",
    )
    r.set_defaults(func=cmd_run)

    obj = sub.add_parser("objectives", help="List the available analysis objectives")
    obj.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    obj.set_defaults(func=cmd_objectives)

    c = sub.add_parser("call", help="Call variants for one BAM/CRAM")
    add_common(c)
    c.add_argument("--bam", required=True)
    c.add_argument("--sample-id", required=True)
    c.add_argument("--fasta", default=None)
    c.add_argument("--callable-bed", default=None)
    c.add_argument("--chrom", default=None)
    c.add_argument("--start", type=int, default=None)
    c.add_argument("--end", type=int, default=None)
    c.add_argument("--max-candidates", type=int, default=None)
    c.add_argument("--out", default="results")
    c.set_defaults(func=cmd_call)

    cal = sub.add_parser("calibrate", help="Fit error rate and over-dispersion")
    add_common(cal)
    cal.add_argument("--controls", required=True, help="TSV with depth and alt_count columns")
    cal.add_argument("--out", default=None)
    cal.set_defaults(func=cmd_calibrate)

    b = sub.add_parser("benchmark", help="Run benchmark and ablation study")
    add_common(b)
    b.add_argument("--out", default="benchmarks/results")
    b.add_argument("--replicates", type=int, default=5)
    b.add_argument("--coverage", type=int, default=60)
    b.add_argument("--seed", type=int, default=11)
    b.add_argument("--no-ablation", action="store_true")
    b.set_defaults(func=cmd_benchmark)

    rep = sub.add_parser("report", help="Export reports or launch the dashboard")
    add_common(rep)
    rep.add_argument("what", choices=["workbook", "pdf", "all", "dashboard"])
    rep.add_argument("--evidence", default=None)
    rep.add_argument("--results", default="results")
    rep.add_argument("--sample-id", default="sample")
    rep.add_argument("--out", default="reports")
    rep.set_defaults(func=cmd_report)

    par = sub.add_parser("parameters", help="Print or export the parameter dictionary")
    par.add_argument("--out", default=None)
    par.add_argument("--mode", default=None)
    par.set_defaults(func=cmd_parameters)

    s = sub.add_parser("status", help="Print the release-status matrix")
    s.set_defaults(func=cmd_status)

    d = sub.add_parser("demo", help="Generate synthetic data and run an end-to-end call")
    add_common(d)
    d.add_argument("--out", default="demo")
    d.add_argument("--coverage", type=int, default=120)
    d.add_argument("--seed", type=int, default=7)
    d.set_defaults(func=cmd_demo)

    imp = sub.add_parser("import", help="Import a pre-called VCF/BCF/gVCF (no realignment)")
    add_common(imp)
    imp.add_argument("--vcf", required=True)
    imp.add_argument("--sample-id", required=True)
    imp.add_argument("--max-records", type=int, default=None)
    imp.add_argument("--out", default="results")
    imp.set_defaults(func=cmd_import)

    ann = sub.add_parser("annotate", help="Annotate an evidence table")
    add_common(ann)
    ann.add_argument("--evidence", required=True)
    ann.add_argument("--sample-id", default="sample")
    ann.add_argument("--gff", default=None)
    ann.add_argument("--fasta", default=None)
    ann.add_argument("--catalogue", default=None)
    ann.add_argument("--out", required=True)
    ann.set_defaults(func=cmd_annotate)

    dn = sub.add_parser("downstream", help="Diversity and mixed-infection summaries")
    add_common(dn)
    dn.add_argument("--evidence", required=True)
    dn.add_argument("--sample-id", default="sample")
    dn.add_argument("--callable-bases", type=int, default=None)
    dn.add_argument("--out", default="results")
    dn.set_defaults(func=cmd_downstream)

    rel = sub.add_parser("relatedness", help="Pairwise relatedness and NJ tree")
    add_common(rel)
    rel.add_argument("--results", default="results")
    rel.add_argument("--min-depth", type=int, default=10)
    rel.add_argument("--out", default="results")
    rel.set_defaults(func=cmd_relatedness)

    qc = sub.add_parser("qc", help="Unified QC across a results tree")
    add_common(qc)
    qc.add_argument("--results", default="results")
    qc.add_argument("--out", default="results/run_qc.json")
    qc.set_defaults(func=cmd_qc)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    print(f"[PfVarX] {BANNER}", file=sys.stderr)
    try:
        return args.func(args)
    except PfVarXError as exc:
        print(json.dumps(exc.as_dict(), indent=2), file=sys.stderr)
        return 3
    except FileNotFoundError as exc:
        print(f"File not found: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

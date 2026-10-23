"""Objective-driven execution engine.

The engine executes *only* the analyses the user selected.  It never runs a
downstream analysis because it happens to exist; it runs the objectives in
``ResolvedObjectives.selected`` and records an explicit outcome for each of
them, including the ones that could not run.

Design contract
---------------
* Every objective produces a :class:`ObjectiveOutcome` with one of four
  statuses: ``completed``, ``failed``, ``skipped`` or ``unavailable``.
  ``unavailable`` means a prerequisite (reference GFF, metadata table,
  catalogue, a second sample) was absent; it is *not* the same as "no signal".
* Per-sample objectives write into ``<root>/objectives/<objective_id>/<sample>/``.
  Run-scope objectives write into ``<root>/objectives/<objective_id>/``.
* The manifest is the single source of truth consumed by the dashboard and the
  exporters.  Nothing in the report layer infers what "should" be present.
"""

from __future__ import annotations

import csv
import json
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from pfvarx.errors import PfVarXError
from pfvarx.objectives import OBJECTIVES, ResolvedObjectives
from pfvarx.workspace import DIR_LOGS, DIR_OBJECTIVES

STATUS_COMPLETED = "completed"
STATUS_FAILED = "failed"
STATUS_SKIPPED = "skipped"
STATUS_UNAVAILABLE = "unavailable"

ALL_STATUSES = (STATUS_COMPLETED, STATUS_FAILED, STATUS_SKIPPED, STATUS_UNAVAILABLE)


# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #
class RunLogger:
    """Append-only run log with in-memory warning/error capture."""

    def __init__(self, path: Path | None = None) -> None:
        self.path = path
        self.warnings: list[str] = []
        self.errors: list[str] = []
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)

    def _write(self, level: str, message: str) -> None:
        line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} [{level}] {message}"
        if level == "WARNING":
            self.warnings.append(message)
        elif level == "ERROR":
            self.errors.append(message)
        if self.path is not None:
            with self.path.open("a", encoding="utf-8") as fh:
                fh.write(line + "\n")

    def info(self, message: str) -> None:
        self._write("INFO", message)

    def warning(self, message: str) -> None:
        self._write("WARNING", message)

    def error(self, message: str) -> None:
        self._write("ERROR", message)


# --------------------------------------------------------------------------- #
# Outcomes and manifest
# --------------------------------------------------------------------------- #
@dataclass
class ObjectiveOutcome:
    objective_id: str
    label: str
    status: str
    scope: str
    artefacts: dict[str, str] = field(default_factory=dict)
    samples: list[str] = field(default_factory=list)
    duration_s: float = 0.0
    message: str = ""
    error: dict[str, Any] | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "objective_id": self.objective_id,
            "label": self.label,
            "status": self.status,
            "scope": self.scope,
            "artefacts": dict(self.artefacts),
            "samples": list(self.samples),
            "duration_s": round(self.duration_s, 3),
            "message": self.message,
            "error": self.error,
        }


@dataclass
class RunManifest:
    run_id: str
    root: str
    started_at: str
    finished_at: str = ""
    duration_s: float = 0.0
    objectives: list[str] = field(default_factory=list)
    stages: list[str] = field(default_factory=list)
    routes: dict[str, int] = field(default_factory=dict)
    samples: list[str] = field(default_factory=list)
    outcomes: list[ObjectiveOutcome] = field(default_factory=list)
    directories: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    status: str = "running"

    def counts(self) -> dict[str, int]:
        out = dict.fromkeys(ALL_STATUSES, 0)
        for o in self.outcomes:
            out[o.status] = out.get(o.status, 0) + 1
        return out

    def as_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "root": self.root,
            "status": self.status,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "duration_s": round(self.duration_s, 3),
            "objectives": list(self.objectives),
            "stages": list(self.stages),
            "routes": dict(self.routes),
            "samples": list(self.samples),
            "directories": list(self.directories),
            "counts": self.counts(),
            "outcomes": [o.as_dict() for o in self.outcomes],
            "warnings": list(self.warnings),
            "errors": list(self.errors),
        }


# --------------------------------------------------------------------------- #
# Small IO helpers
# --------------------------------------------------------------------------- #
def _write_json(path: Path, payload: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return str(path)


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return str(path)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _scalar(v) for k, v in row.items()})
    return str(path)


def _scalar(value: Any) -> Any:
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(v) for v in value)
    if isinstance(value, dict):
        return json.dumps(value, default=str)
    return value


def load_calls(evidence_path: str | Path) -> list:
    """Reconstruct :class:`SiteCall` objects from an evidence table."""
    from pfvarx.calling.policy import SiteCall

    with open(evidence_path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [SiteCall.from_evidence_row(r) for r in rows]


def load_metadata(
    samples: list,
    cfg: Any,
    *,
    explicit_path: str | None = None,
) -> pd.DataFrame | None:
    """Build a sample metadata frame from the samplesheet or a metadata table.

    Precedence: an explicitly supplied metadata table wins; otherwise the
    per-sample ``metadata_json`` column of the samplesheet is used.  Returns
    ``None`` when no metadata is available at all.
    """
    path = explicit_path or getattr(cfg.downstream, "metadata_table", None)
    if path:
        p = Path(path)
        if not p.exists():
            raise PfVarXError(
                f"Metadata table not found: {p}",
                cause="downstream.metadata_table points at a missing file",
                action="Correct the path or remove the setting",
            )
        sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(p, sep=sep)
        if "sample_id" not in frame.columns:
            raise PfVarXError(
                f"Metadata table {p} has no sample_id column",
                action="Add a sample_id column that matches the samplesheet",
            )
        return frame

    rows = [{"sample_id": s.sample_id, **dict(getattr(s, "metadata", {}) or {})} for s in samples]
    if not any(len(getattr(s, "metadata", {}) or {}) for s in samples):
        return None
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Per-sample objective implementations
# --------------------------------------------------------------------------- #
def _obj_drug_resistance(calls: list, cfg: Any, sample_id: str, out: Path, logger: RunLogger) -> tuple[str, dict, str]:
    from pfvarx.annotation.annotator import annotate_calls
    from pfvarx.annotation.catalogue import load_catalogue

    gff = getattr(cfg.annotation, "gff", None)
    if not gff:
        return (
            STATUS_UNAVAILABLE,
            {},
            "No gene model supplied: set annotation.gff (GFF3) so amino-acid changes can be "
            "derived. Without it a resistance table would be empty for the wrong reason.",
        )
    catalogue = load_catalogue(getattr(cfg.downstream, "resistance_catalogue", None))
    annotated = annotate_calls(
        calls, cfg, gff_path=gff, fasta_path=cfg.reference.fasta, catalogue=catalogue
    )
    hits: list[dict] = []
    for site in annotated:
        for marker in site.resistance:
            row = marker.copy()
            row.update(
                {
                    "sample_id": sample_id,
                    "chrom": site.call.evidence.chrom,
                    "pos": site.call.evidence.pos,
                    "ref": site.call.evidence.ref,
                    "alt": site.call.evidence.alt,
                    "alt_fraction": site.call.evidence.alt_fraction,
                    "depth": site.call.evidence.informative_depth,
                    "decision": site.call.decision,
                }
            )
            hits.append(row)

    hits_path = _write_tsv(out / "resistance_hits.tsv", hits)
    by_drug: dict[str, int] = {}
    by_gene: dict[str, int] = {}
    for row in hits:
        by_drug[row["drug"]] = by_drug.get(row["drug"], 0) + 1
        by_gene[row["gene_name"]] = by_gene.get(row["gene_name"], 0) + 1
    summary = {
        "sample_id": sample_id,
        "catalogue_version": catalogue.version,
        "catalogue_assembly": catalogue.assembly,
        "clinical_use": catalogue.clinical_use,
        "n_calls_screened": len(calls),
        "n_marker_hits": len(hits),
        "hits_by_drug": by_drug,
        "hits_by_gene": by_gene,
        "note": (
            "Marker presence is not proof of treatment failure and this catalogue is "
            "not for clinical decision-making."
        ),
    }
    summary_path = _write_json(out / "resistance_summary.json", summary)
    msg = f"{len(hits)} resistance-marker hit(s) across {len(by_gene)} gene(s)"
    return STATUS_COMPLETED, {"resistance_hits": hits_path, "resistance_summary": summary_path}, msg


def _obj_mixed_infection(calls: list, cfg: Any, sample_id: str, out: Path, logger: RunLogger) -> tuple[str, dict, str]:
    from pfvarx.downstream.mixed_infection import mixed_infection_summary

    summary = mixed_infection_summary(calls, sample_id=sample_id)
    path = _write_json(out / "mixed_infection.json", summary.as_dict())
    return STATUS_COMPLETED, {"mixed_infection": path}, f"flag={summary.mixed_infection_flag}"


def _obj_coi(calls: list, cfg: Any, sample_id: str, out: Path, logger: RunLogger) -> tuple[str, dict, str]:
    from pfvarx.downstream.mixed_infection import mixed_infection_summary

    mix = mixed_infection_summary(calls, sample_id=sample_id)
    payload = {
        "sample_id": sample_id,
        "coi_point_estimate": mix.coi_point_estimate,
        "coi_method": "unavailable",
        "ordinal_indicator": mix.mixed_infection_flag,
        "n_intermediate_sites": mix.n_intermediate,
        "n_minority_sites": mix.n_minority,
        "max_minor_fraction": mix.max_minor_fraction,
        "note": (
            "No validated COI model is implemented. A numeric estimate derived from "
            "intermediate-fraction counts would be fabricated; use a validated method "
            "such as DEploid for a numeric COI."
        ),
    }
    path = _write_json(out / "coi.json", payload)
    return STATUS_COMPLETED, {"coi": path}, "ordinal indicator only; numeric COI not implemented"


def _obj_diversity(calls: list, cfg: Any, sample_id: str, out: Path, logger: RunLogger) -> tuple[str, dict, str]:
    from pfvarx.downstream.diversity import diversity_summary

    summary = diversity_summary(calls, sample_id=sample_id)
    path = _write_json(out / "diversity.json", summary.as_dict())
    return STATUS_COMPLETED, {"diversity": path}, f"pi={summary.nucleotide_diversity:.4g}"


# --------------------------------------------------------------------------- #
# Run-scope objective implementations
# --------------------------------------------------------------------------- #
def _obj_population_structure(
    calls_by_sample: dict[str, list], metadata: pd.DataFrame | None, cfg: Any, out: Path, logger: RunLogger
) -> tuple[str, dict, str]:
    from pfvarx.downstream.diversity import diversity_summary
    from pfvarx.downstream.mixed_infection import mixed_infection_summary
    from pfvarx.downstream.spatial_temporal import aggregate_by_metadata

    if metadata is None:
        return STATUS_UNAVAILABLE, {}, "No sample metadata found (samplesheet metadata_json or downstream.metadata_table)."

    per_sample = []
    for sid, calls in calls_by_sample.items():
        div = diversity_summary(calls, sample_id=sid)
        mix = mixed_infection_summary(calls, sample_id=sid)
        per_sample.append(
            {
                "sample_id": sid,
                "n_sites": div.n_sites,
                "proportion_polymorphic": div.proportion_polymorphic,
                "nucleotide_diversity": div.nucleotide_diversity,
                "expected_heterozygosity": div.expected_heterozygosity,
                "n_intermediate_sites": mix.n_intermediate,
            }
        )
    group_cols = [c for c in metadata.columns if c != "sample_id"][:2]
    if not group_cols:
        return STATUS_UNAVAILABLE, {}, "Metadata has no grouping columns beyond sample_id."
    table = aggregate_by_metadata(per_sample, metadata, group_by=group_cols)
    rows = table.to_dict(orient="records")
    tsv = _write_tsv(out / "population_structure.tsv", rows)
    payload = {
        "group_by": group_cols,
        "n_samples": len(per_sample),
        "n_groups": len(rows),
        "table": rows,
        "note": "Descriptive aggregation only; no population-genetic model is fitted.",
    }
    js = _write_json(out / "population_structure.json", payload)
    return STATUS_COMPLETED, {"population_structure_tsv": tsv, "population_structure_json": js}, f"{len(rows)} group(s) by {group_cols}"


def _obj_relatedness(
    calls_by_sample: dict[str, list], metadata: pd.DataFrame | None, cfg: Any, out: Path, logger: RunLogger
) -> tuple[str, dict, str]:
    from pfvarx.downstream.relatedness import relatedness_matrix

    if len(calls_by_sample) < 2:
        return STATUS_UNAVAILABLE, {}, "Relatedness requires at least two samples with called variants."
    results = relatedness_matrix(calls_by_sample)
    rows = [r.as_dict() for r in results]
    js = _write_json(out / "relatedness.json", rows)
    tsv = _write_tsv(out / "relatedness.tsv", rows)
    return STATUS_COMPLETED, {"relatedness_json": js, "relatedness_tsv": tsv}, f"{len(rows)} pair(s)"


def _obj_phylogenetics(
    calls_by_sample: dict[str, list], metadata: pd.DataFrame | None, cfg: Any, out: Path, logger: RunLogger
) -> tuple[str, dict, str]:
    from pfvarx.downstream.phylogenetics import distance_matrix, newick
    from pfvarx.downstream.relatedness import relatedness_matrix

    if len(calls_by_sample) < 2:
        return STATUS_UNAVAILABLE, {}, "Phylogenetics requires at least two samples with called variants."
    samples = sorted(calls_by_sample)
    results = relatedness_matrix(calls_by_sample)
    tree = newick(samples, results)
    tree_path = out / "tree.nwk"
    tree_path.parent.mkdir(parents=True, exist_ok=True)
    tree_path.write_text(tree, encoding="utf-8")
    dist = distance_matrix(results, samples)
    rows = [{"sample_a": a, "sample_b": b, "distance": d} for (a, b), d in dist.items()]
    tsv = _write_tsv(out / "distance_matrix.tsv", rows)
    return STATUS_COMPLETED, {"tree": str(tree_path), "distance_matrix": tsv}, f"{len(samples)} taxa"


def _obj_spatial_temporal(
    calls_by_sample: dict[str, list], metadata: pd.DataFrame | None, cfg: Any, out: Path, logger: RunLogger
) -> tuple[str, dict, str]:
    from pfvarx.downstream.diversity import diversity_summary
    from pfvarx.downstream.spatial_temporal import aggregate_by_metadata, temporal_series

    if metadata is None:
        return STATUS_UNAVAILABLE, {}, "No sample metadata found (samplesheet metadata_json or downstream.metadata_table)."
    per_sample = []
    for sid, calls in calls_by_sample.items():
        div = diversity_summary(calls, sample_id=sid)
        per_sample.append(
            {
                "sample_id": sid,
                "n_sites": div.n_sites,
                "nucleotide_diversity": div.nucleotide_diversity,
                "proportion_polymorphic": div.proportion_polymorphic,
            }
        )
    artefacts: dict[str, str] = {}
    messages: list[str] = []
    group_cols = [c for c in metadata.columns if c != "sample_id"][:2]
    if group_cols:
        agg = aggregate_by_metadata(per_sample, metadata, group_by=group_cols)
        artefacts["spatial_temporal_tsv"] = _write_tsv(out / "spatial_temporal.tsv", agg.to_dict(orient="records"))
        messages.append(f"grouped by {group_cols}")
    else:
        messages.append("no grouping columns in metadata")
    if "collection_date" in metadata.columns:
        series = temporal_series(per_sample, metadata)
        artefacts["temporal_series_tsv"] = _write_tsv(out / "temporal_series.tsv", series.to_dict(orient="records"))
        messages.append("temporal series written")
    payload = {
        "n_samples": len(per_sample),
        "group_by": group_cols,
        "has_collection_date": "collection_date" in metadata.columns,
        "artefacts": artefacts,
        "note": "Descriptive aggregation only; no spatio-temporal model is fitted.",
    }
    artefacts["spatial_temporal_json"] = _write_json(out / "spatial_temporal.json", payload)
    return STATUS_COMPLETED, artefacts, "; ".join(messages)


# --------------------------------------------------------------------------- #
# Dispatch tables
# --------------------------------------------------------------------------- #
PER_SAMPLE_RUNNERS: dict[str, Callable[..., tuple[str, dict, str]]] = {
    "drug_resistance": _obj_drug_resistance,
    "mixed_infection": _obj_mixed_infection,
    "coi": _obj_coi,
    "diversity": _obj_diversity,
}

RUN_SCOPE_RUNNERS: dict[str, Callable[..., tuple[str, dict, str]]] = {
    "population_structure": _obj_population_structure,
    "relatedness": _obj_relatedness,
    "phylogenetics": _obj_phylogenetics,
    "spatial_temporal": _obj_spatial_temporal,
}


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def execute_objectives(
    *,
    root: str | Path,
    resolved: ResolvedObjectives,
    calls_by_sample: dict[str, list],
    samples: list,
    cfg: Any,
    metadata: pd.DataFrame | None = None,
    logger: RunLogger | None = None,
    run_id: str | None = None,
) -> RunManifest:
    """Execute every selected objective and return the run manifest.

    Objectives whose prerequisites are missing are recorded as ``unavailable``
    with an actionable message rather than silently producing empty output.
    """
    root = Path(root)
    logger = logger or RunLogger(root / DIR_LOGS / "run.log")
    started = time.time()
    manifest = RunManifest(
        run_id=run_id or getattr(cfg.run, "run_id", None) or "run",
        root=str(root),
        started_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        objectives=list(resolved.selected),
        stages=list(resolved.stages),
        samples=sorted(calls_by_sample),
    )

    for objective_id in resolved.selected:
        spec = OBJECTIVES[objective_id]
        t0 = time.time()
        if not calls_by_sample:
            logger.warning(f"objective {objective_id}: no sample produced called variants")
            manifest.outcomes.append(
                ObjectiveOutcome(
                    objective_id=objective_id,
                    label=spec.label,
                    status=STATUS_UNAVAILABLE,
                    scope=spec.scope,
                    duration_s=0.0,
                    message=(
                        "No sample produced called variants, so this objective had nothing to "
                        "analyse. Check the per-sample statuses in run_summary.json."
                    ),
                )
            )
            continue
        logger.info(f"objective start: {objective_id} (scope={spec.scope})")
        try:
            if spec.scope == "per_sample":
                manifest.outcomes.extend(
                    _run_per_sample(objective_id, spec, calls_by_sample, cfg, root, logger)
                )
            else:
                status, artefacts, message = RUN_SCOPE_RUNNERS[objective_id](
                    calls_by_sample, metadata, cfg, root / DIR_OBJECTIVES / objective_id, logger
                )
                manifest.outcomes.append(
                    ObjectiveOutcome(
                        objective_id=objective_id,
                        label=spec.label,
                        status=status,
                        scope=spec.scope,
                        artefacts=artefacts,
                        samples=sorted(calls_by_sample),
                        duration_s=time.time() - t0,
                        message=message,
                    )
                )
        except PfVarXError as exc:
            logger.error(f"objective {objective_id} failed: {exc}")
            manifest.outcomes.append(
                ObjectiveOutcome(
                    objective_id=objective_id,
                    label=spec.label,
                    status=STATUS_FAILED,
                    scope=spec.scope,
                    samples=sorted(calls_by_sample),
                    duration_s=time.time() - t0,
                    message=str(exc),
                    error=exc.as_dict(),
                )
            )
        except Exception as exc:  # noqa: BLE001 - recorded, never swallowed silently
            logger.error(f"objective {objective_id} failed unexpectedly: {exc!r}")
            manifest.outcomes.append(
                ObjectiveOutcome(
                    objective_id=objective_id,
                    label=spec.label,
                    status=STATUS_FAILED,
                    scope=spec.scope,
                    samples=sorted(calls_by_sample),
                    duration_s=time.time() - t0,
                    message=f"{type(exc).__name__}: {exc}",
                )
            )
        logger.info(f"objective end: {objective_id} -> {manifest.outcomes[-1].status}")

    manifest.duration_s = time.time() - started
    manifest.finished_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    manifest.warnings = list(logger.warnings)
    manifest.errors = list(logger.errors)
    manifest.status = _overall_status(manifest)
    return manifest


def _run_per_sample(
    objective_id: str, spec: Any, calls_by_sample: dict[str, list], cfg: Any, root: Path, logger: RunLogger
) -> list[ObjectiveOutcome]:
    runner = PER_SAMPLE_RUNNERS[objective_id]
    outcomes: list[ObjectiveOutcome] = []
    for sample_id in sorted(calls_by_sample):
        t0 = time.time()
        out = root / DIR_OBJECTIVES / objective_id / sample_id
        try:
            status, artefacts, message = runner(calls_by_sample[sample_id], cfg, sample_id, out, logger)
        except PfVarXError as exc:
            logger.error(f"{objective_id}/{sample_id} failed: {exc}")
            outcomes.append(
                ObjectiveOutcome(
                    objective_id=objective_id,
                    label=spec.label,
                    status=STATUS_FAILED,
                    scope="per_sample",
                    samples=[sample_id],
                    duration_s=time.time() - t0,
                    message=str(exc),
                    error=exc.as_dict(),
                )
            )
            continue
        except Exception as exc:  # noqa: BLE001
            logger.error(f"{objective_id}/{sample_id} failed unexpectedly: {exc!r}")
            outcomes.append(
                ObjectiveOutcome(
                    objective_id=objective_id,
                    label=spec.label,
                    status=STATUS_FAILED,
                    scope="per_sample",
                    samples=[sample_id],
                    duration_s=time.time() - t0,
                    message=f"{type(exc).__name__}: {exc}",
                )
            )
            continue
        outcomes.append(
            ObjectiveOutcome(
                objective_id=objective_id,
                label=spec.label,
                status=status,
                scope="per_sample",
                artefacts=artefacts,
                samples=[sample_id],
                duration_s=time.time() - t0,
                message=message,
            )
        )
    return outcomes


def _overall_status(manifest: RunManifest) -> str:
    statuses = {o.status for o in manifest.outcomes}
    if not statuses:
        return STATUS_SKIPPED
    if statuses == {STATUS_COMPLETED}:
        return STATUS_COMPLETED
    if STATUS_FAILED in statuses:
        return "partial_failure"
    if STATUS_COMPLETED in statuses:
        return "partial"
    return STATUS_UNAVAILABLE


def write_manifest(manifest: RunManifest, root: str | Path) -> str:
    return _write_json(Path(root) / "run_manifest.json", manifest.as_dict())


def format_summary(manifest: RunManifest) -> str:
    """Human-readable execution summary for the terminal."""
    counts = manifest.counts()
    lines = [
        "",
        "=" * 78,
        f"PfVarX execution summary  (run_id={manifest.run_id})",
        "=" * 78,
        f"root            : {manifest.root}",
        f"status          : {manifest.status}",
        f"duration        : {manifest.duration_s:.2f} s",
        f"samples         : {len(manifest.samples)}",
        f"objectives      : {', '.join(manifest.objectives) or '(none)'}",
        f"stages          : {', '.join(manifest.stages) or '(none)'}",
        f"outcomes        : completed={counts['completed']} failed={counts['failed']} "
        f"skipped={counts['skipped']} unavailable={counts['unavailable']}",
        "-" * 78,
        f"{'OBJECTIVE':<24}{'STATUS':<14}{'SCOPE':<12}DETAIL",
    ]
    for outcome in manifest.outcomes:
        detail = outcome.message or ""
        if len(detail) > 60:
            detail = detail[:57] + "..."
        lines.append(f"{outcome.objective_id:<24}{outcome.status:<14}{outcome.scope:<12}{detail}")
    if manifest.warnings:
        lines.append("-" * 78)
        lines.append("warnings:")
        lines.extend(f"  - {w}" for w in manifest.warnings[:20])
    lines.append("=" * 78)
    lines.append(f"Manifest: {Path(manifest.root) / 'run_manifest.json'}")
    lines.append("")
    return "\n".join(lines)

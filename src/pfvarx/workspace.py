"""Filesystem planning for input preparation and run outputs.

Stage 1 uses :func:`create_input_directories` to create only configured input
locations. Stage 2 uses :func:`plan_layout`, :func:`validate`, and :func:`create`
to prepare outputs, logs, and work directories after run preflight succeeds.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pfvarx.errors import ConfigError

DIR_LOGS = "logs"
DIR_CONFIG = "config"
DIR_QC = "qc"
DIR_INTERMEDIATE = "intermediate"
DIR_ALIGNMENTS = "alignments"
DIR_VARIANTS = "variants"
DIR_OBJECTIVES = "objectives"
DIR_REPORTS = "reports"
DIR_DASHBOARD = "dashboard"


def _absolute(path: str | Path, base_dir: str | Path | None = None) -> Path:
    path = Path(path).expanduser()
    if path.is_absolute():
        return path.resolve(strict=False)
    return (Path(base_dir or Path.cwd()) / path).resolve(strict=False)


def resolve_output_directories(
    cfg: Any,
    root: str | Path | None = None,
    *,
    base_dir: str | Path | None = None,
) -> tuple[Path, Path]:
    """Return the effective logs and work directories for a run.

    Explicit paths are interpreted relative to ``base_dir`` (the current
    working directory by default). Otherwise they are derived from the selected
    final-results root.
    """
    out_root = Path(root if root is not None else cfg.outputs.root).expanduser()
    if not out_root.is_absolute():
        out_root = Path(base_dir or Path.cwd()) / out_root
    out_root = out_root.resolve(strict=False)
    logs_value = getattr(cfg.outputs, "logs_dir", None)
    work_value = getattr(cfg.outputs, "work_dir", None)
    logs = _absolute(logs_value, base_dir) if logs_value else out_root / DIR_LOGS
    work = _absolute(work_value, base_dir) if work_value else out_root / DIR_INTERMEDIATE
    return logs, work


def plan_input_directories(cfg: Any, *, base_dir: str | Path | None = None) -> list[Path]:
    """Plan directories needed to stage inputs, without creating anything.

    The plan includes the configured raw/prepared/metadata folders and parent
    folders for declared input files. It deliberately excludes every output,
    log, and work path and rejects input locations placed inside those paths.
    """
    base = Path(base_dir or Path.cwd()).resolve()
    inputs = cfg.inputs
    input_root = _absolute(getattr(inputs, "root", "data"), base)

    configured = (
        (getattr(inputs, "raw_dir", None), input_root / "raw", True),
        (getattr(inputs, "prepared_dir", None), input_root / "prepared", True),
        (getattr(inputs, "metadata_dir", None), input_root / "metadata", False),
    )
    planned: list[Path] = []
    for value, default, create_default in configured:
        if value is None and not create_default:
            continue
        planned.append(_absolute(value, base) if value else default)

    file_values = [
        getattr(inputs, "samplesheet", None),
        getattr(cfg.reference, "fasta", None),
        getattr(cfg.reference, "callable_mask_bed", None),
        getattr(cfg.preprocessing, "primer_bed", None),
        getattr(cfg.annotation, "gff", None),
        getattr(cfg.downstream, "metadata_table", None),
        getattr(cfg.downstream, "resistance_catalogue", None),
    ]
    for value in file_values:
        if value:
            planned.append(_absolute(value, base).parent)

    output_root = _absolute(cfg.outputs.root, base)
    logs_dir, work_dir = resolve_output_directories(cfg, output_root, base_dir=base)
    output_locations = (output_root, logs_dir, work_dir)
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in planned:
        path = path.resolve(strict=False)
        if path in seen:
            continue
        seen.add(path)
        for output in output_locations:
            if path == output or output in path.parents:
                raise ConfigError(
                    f"Input directory overlaps a configured output directory: {path}",
                    cause=f"Output location: {output}",
                    action="Move input files outside outputs.root, outputs.logs_dir, and outputs.work_dir.",
                    consequence="Input initialization was stopped; no directories were created.",
                )
        unique.append(path)
    return unique


def create_input_directories(cfg: Any, *, base_dir: str | Path | None = None) -> list[Path]:
    """Create only the Stage 1 input directories and return newly created paths."""
    paths = plan_input_directories(cfg, base_dir=base_dir)
    created: list[Path] = []
    for path in paths:
        if path.exists() and not path.is_dir():
            raise ConfigError(
                f"Input directory path exists but is not a directory: {path}",
                action="Choose a directory path or move the conflicting file.",
            )
        if not path.exists():
            try:
                path.mkdir(parents=True, exist_ok=True)
            except OSError as exc:
                raise ConfigError(
                    f"Cannot create input directory: {path}",
                    cause=str(exc),
                    action="Choose a writable input location and rerun `pfvarx init`.",
                ) from exc
            created.append(path)
    return created


@dataclass
class WorkspaceLayout:
    """The directories required for one Stage 2 run."""

    root: Path
    #: Paths relative to ``root`` or absolute external log/work paths.
    directories: list[str] = field(default_factory=list)
    objective_dirs: dict[str, str] = field(default_factory=dict)
    rationale: dict[str, str] = field(default_factory=dict)

    def path(self, directory: str) -> Path:
        value = Path(directory)
        return value if value.is_absolute() else self.root / value

    def as_dict(self) -> dict:
        return {
            "root": str(self.root),
            "directories": list(self.directories),
            "objective_dirs": dict(self.objective_dirs),
            "rationale": dict(self.rationale),
        }


def _layout_key(root: Path, path: str | Path) -> str:
    target = Path(path).resolve(strict=False)
    root_abs = root.resolve(strict=False)
    if target == root_abs:
        raise ConfigError("A log or work directory cannot equal outputs.root.")
    try:
        return target.relative_to(root_abs).as_posix()
    except ValueError:
        return str(target)


def plan_layout(
    root: str | Path,
    *,
    objectives: list[str] | None = None,
    routes: list[str] | None = None,
    needs_annotation: bool = False,
    logs_dir: str | Path | None = None,
    work_dir: str | Path | None = None,
) -> WorkspaceLayout:
    """Compute Stage 2 output paths without touching the filesystem."""
    root = Path(root)
    objectives = list(objectives or [])
    routes = set(routes or [])
    layout = WorkspaceLayout(root=root)

    def add(name: str, why: str) -> None:
        if name not in layout.directories:
            layout.directories.append(name)
            layout.rationale[name] = why

    if logs_dir is None:
        add(DIR_LOGS, "Run logs are written here (or to outputs.logs_dir if configured).")
    else:
        key = _layout_key(root, logs_dir)
        add(key, "Run logs are written to configured outputs.logs_dir.")
    add(DIR_CONFIG, "The resolved configuration and routing plan are snapshotted for reproducibility.")
    add(DIR_QC, "Read-level and variant-level QC is produced for every input type.")

    if "fastq_alignment" in routes:
        add(DIR_ALIGNMENTS, "FASTQ input: alignment products are written under the work directory.")
    if routes.intersection({"fastq_alignment", "alignment_reuse"}):
        if work_dir is None:
            add(DIR_INTERMEDIATE, "Intermediate files use the default work directory.")
        else:
            key = _layout_key(root, work_dir)
            add(key, "Intermediate files use configured outputs.work_dir.")

    if routes:
        add(DIR_VARIANTS, "VCF and evidence tables for called or imported variants.")
    if objectives:
        add(DIR_OBJECTIVES, "One subdirectory per selected analysis objective.")
        for objective in objectives:
            rel = f"{DIR_OBJECTIVES}/{objective}"
            layout.objective_dirs[objective] = rel
            add(rel, f"Selected objective: {objective}.")
    add(DIR_REPORTS, "Run-level report exports.")
    add(DIR_DASHBOARD, "Dashboard configuration describing what this run produced.")
    return layout


def validate(layout: WorkspaceLayout, *, overwrite: bool = False) -> list[str]:
    """Check all planned Stage 2 paths without creating them."""
    warnings: list[str] = []
    if layout.root.exists() and not layout.root.is_dir():
        raise ConfigError(
            f"Output root exists but is not a directory: {layout.root}",
            action="Choose a different outputs.root or remove the conflicting file.",
        )
    for directory in layout.directories:
        target = layout.path(directory)
        if target.exists() and not target.is_dir():
            raise ConfigError(
                f"Output directory exists but is not a directory: {target}",
                action="Choose another outputs.root, outputs.logs_dir, or outputs.work_dir.",
            )
        probe = target
        while not probe.exists() and probe != probe.parent:
            probe = probe.parent
        if probe.exists() and not os.access(probe, os.W_OK):
            raise ConfigError(
                f"Cannot write to planned output directory: {target}",
                cause=f"Nearest existing ancestor {probe} is not writable",
                action="Choose a writable output, log, and work location.",
            )
        if target.exists() and not overwrite:
            warnings.append(
                f"Output directory already exists and outputs.overwrite is false: {target}. "
                "Existing files are not removed; file-level collisions are checked by the CLI."
            )
    return warnings


def create(layout: WorkspaceLayout) -> list[str]:
    """Create exactly the planned Stage 2 directories."""
    created: list[str] = []
    for directory in layout.directories:
        target = layout.path(directory)
        if not target.exists():
            target.mkdir(parents=True, exist_ok=True)
            created.append(directory)
    return created


def prepare(
    root: str | Path,
    *,
    objectives: list[str] | None = None,
    routes: list[str] | None = None,
    overwrite: bool = False,
    logs_dir: str | Path | None = None,
    work_dir: str | Path | None = None,
) -> tuple[WorkspaceLayout, list[str], list[str]]:
    """Plan, validate, and create a Stage 2 workspace."""
    layout = plan_layout(
        root,
        objectives=objectives,
        routes=routes,
        logs_dir=logs_dir,
        work_dir=work_dir,
    )
    warnings = validate(layout, overwrite=overwrite)
    created = create(layout)
    return layout, created, warnings

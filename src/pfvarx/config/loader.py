"""Configuration loading, profile resolution and provenance capture."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import yaml

from pfvarx.config import parameters as P
from pfvarx.config.schema import PfVarXConfig
from pfvarx.errors import ConfigError

PROFILE_DIR = Path(__file__).resolve().parents[3] / "config" / "profiles"


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def read_yaml(path: str | Path) -> dict:
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"Configuration file not found: {p}", action="Check the path")
    try:
        data = yaml.safe_load(p.read_text())
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in {p}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"Configuration root must be a mapping, got {type(data).__name__}")
    return data


def available_profiles() -> list[str]:
    if not PROFILE_DIR.exists():
        return []
    return sorted(p.stem for p in PROFILE_DIR.glob("*.yaml"))


def load_profile(name_or_path: str) -> dict:
    p = Path(name_or_path)
    if not p.exists():
        p = PROFILE_DIR / f"{name_or_path}.yaml"
    if not p.exists():
        raise ConfigError(
            f"Unknown profile {name_or_path!r}; available: {available_profiles()}",
            action="Use --profile with a known profile name or an explicit YAML path",
        )
    return read_yaml(p)


def load_config(
    path: str | Path,
    *,
    profile: str | None = None,
    overrides: dict[str, Any] | None = None,
    validate_paths: bool = False,
    base_dir: str | Path | None = None,
) -> PfVarXConfig:
    """Load, merge and validate a configuration.

    Precedence: explicit config file > profile defaults > schema/registry defaults.
    """
    data = read_yaml(path)
    # A ``profile:`` key inside the run configuration selects the base profile.
    # It is a loader directive, not a schema field, so it must be removed before
    # validation (``extra="forbid"`` would otherwise reject the shipped example).
    declared_profile = data.pop("profile", None)
    profile = profile or (str(declared_profile) if declared_profile else None)
    if profile:
        data = _deep_merge(load_profile(profile), data)
    if overrides:
        data = _deep_merge(data, overrides)
    try:
        cfg = PfVarXConfig.model_validate(data)
    except ConfigError:
        raise
    except Exception as exc:  # pydantic ValidationError
        raise ConfigError(
            f"Configuration validation failed: {exc}",
            cause=str(exc),
            action="Correct the offending field(s); see docs/PARAMETERS.md",
            consequence="No stage will run",
        ) from exc
    if validate_paths:
        cfg.validate_paths(base_dir=base_dir)
    return cfg


def write_effective_config(cfg: PfVarXConfig, out_path: str | Path) -> Path:
    """Write a fully resolved config plus provenance sidecar."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "config": cfg.model_dump(mode="json"),
        "provenance": {
            "config_hash": cfg.config_hash(),
            "policy_id": cfg.policy_id(),
            "parameter_registry_hash": P.registry_hash(),
            "effective_parameters": cfg.effective_parameters(),
            "unvalidated_parameters": P.unvalidated_parameters(),
            "scientific_warnings": cfg.scientific_warnings(),
        },
    }
    tmp = out.with_suffix(out.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=False))
    tmp.replace(out)
    return out


def export_parameter_dictionary(out_path: str | Path) -> Path:
    """Write the machine-readable parameter dictionary (CSV + JSON sidecar)."""
    import csv

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    records = P.registry_records()
    fields = [
        "name", "category", "unit", "datatype", "default", "minimum", "maximum",
        "allowed", "modes", "rationale", "source", "policy_status",
        "calibration_dataset_id", "version",
    ]
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for rec in records:
            rec = dict(rec)
            rec["allowed"] = ";".join(map(str, rec["allowed"])) if rec["allowed"] else ""
            rec["modes"] = ";".join(rec["modes"])
            writer.writerow(rec)
    out.with_suffix(".json").write_text(json.dumps(records, indent=2, default=str))
    return out

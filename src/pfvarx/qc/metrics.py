"""Unified QC metric model and evaluation.

Every QC metric is emitted with the full provenance required by the spec:
value, unit, source, method, stage, sample, provisional-or-validated flag,
threshold, direction, status, reason and policy id.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from pfvarx.errors import ConfigError

ONTOLOGY_PATH = Path(__file__).resolve().parents[3] / "config" / "schemas" / "qc_ontology.json"

STATUS_PASS = "PASS"
STATUS_WARN = "WARN"
STATUS_FAIL = "FAIL"
STATUS_INFO = "INFO"


@dataclass
class QCMetric:
    metric: str
    value: float | str | None
    unit: str
    source: str
    method: str
    stage: str
    sample: str
    provisional_or_validated: str
    threshold: float | str | None
    direction: str
    status: str
    reason: str
    policy_id: str | None = None
    extra: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        d = {
            "metric": self.metric,
            "value": self.value,
            "unit": self.unit,
            "source": self.source,
            "method": self.method,
            "stage": self.stage,
            "sample": self.sample,
            "provisional_or_validated": self.provisional_or_validated,
            "threshold": self.threshold,
            "direction": self.direction,
            "status": self.status,
            "reason": self.reason,
            "policy_id": self.policy_id,
        }
        d.update(self.extra)
        return d


def load_ontology(path: str | None = None) -> dict:
    p = Path(path) if path else ONTOLOGY_PATH
    try:
        with open(p) as fh:
            return json.load(fh)
    except FileNotFoundError as exc:
        raise ConfigError(
            f"QC ontology not found: {p}",
            cause=str(exc),
            action="Reinstall package data or set an explicit ontology path.",
        ) from exc


def evaluate_metric(
    metric: str,
    value: float | str | None,
    *,
    sample: str = "run",
    stage: str | None = None,
    unit: str | None = None,
    source: str = "pfvarx",
    method: str = "computed",
    threshold: float | str | None = None,
    direction: str | None = None,
    policy_id: str | None = None,
    provisional: bool = True,
    warn_fraction: float = 0.1,
) -> QCMetric:
    """Evaluate a metric against a threshold and return a fully annotated record."""
    ontology = load_ontology()
    spec = ontology["metrics"].get(metric, {})
    unit = unit or spec.get("unit", "unknown")
    stage = stage or spec.get("stage", "unknown")
    direction = direction or spec.get("direction", "none")

    status = STATUS_INFO
    reason = "No threshold configured; informational only."
    if threshold is not None and isinstance(value, (int, float)) and direction in ("min", "max"):
        if direction == "min":
            if value >= threshold:
                status = STATUS_PASS
                reason = f"value {value:.4g} >= minimum {threshold:.4g}"
            elif value >= threshold * (1 - warn_fraction):
                status = STATUS_WARN
                reason = f"value {value:.4g} within {warn_fraction:.0%} below minimum {threshold:.4g}"
            else:
                status = STATUS_FAIL
                reason = f"value {value:.4g} < minimum {threshold:.4g}"
        else:
            if value <= threshold:
                status = STATUS_PASS
                reason = f"value {value:.4g} <= maximum {threshold:.4g}"
            elif value <= threshold * (1 + warn_fraction):
                status = STATUS_WARN
                reason = f"value {value:.4g} within {warn_fraction:.0%} above maximum {threshold:.4g}"
            else:
                status = STATUS_FAIL
                reason = f"value {value:.4g} > maximum {threshold:.4g}"

    return QCMetric(
        metric=metric,
        value=value,
        unit=unit,
        source=source,
        method=method,
        stage=stage,
        sample=sample,
        provisional_or_validated="PROVISIONAL" if provisional else "VALIDATED",
        threshold=threshold,
        direction=direction,
        status=status,
        reason=reason,
        policy_id=policy_id,
    )


def summarise_metrics(metrics: list[QCMetric]) -> dict:
    counts = {STATUS_PASS: 0, STATUS_WARN: 0, STATUS_FAIL: 0, STATUS_INFO: 0}
    for m in metrics:
        counts[m.status] = counts.get(m.status, 0) + 1
    overall = STATUS_PASS
    if counts[STATUS_FAIL]:
        overall = STATUS_FAIL
    elif counts[STATUS_WARN]:
        overall = STATUS_WARN
    return {"overall_status": overall, "counts": counts, "n_metrics": len(metrics)}

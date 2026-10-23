"""Run-level QC integration.

Combines per-sample calling summaries, downstream summaries and run metadata
into a single machine-readable QC table.
"""

from __future__ import annotations

from pfvarx.calling.policy import (
    DECISION_ALT_EXPLORATORY,
    DECISION_ALT_PASS,
    DECISION_ALT_UNCERTAIN,
    summarise_calls,
)
from pfvarx.config.schema import PfVarXConfig
from pfvarx.qc.metrics import QCMetric, evaluate_metric, summarise_metrics


def sample_qc(
    sample_id: str,
    calls,
    cfg: PfVarXConfig,
    *,
    mixed_infection_flag: str | None = None,
    n_polymorphic: int | None = None,
) -> list[QCMetric]:
    summary = summarise_calls(calls)
    decisions = summary["decisions"]
    n = summary["n_sites"]
    # ``UncalibratedModel`` is a policy-level marker applied uniformly to every
    # site in a PROVISIONAL run; it is not a site-specific QC failure, so it is
    # excluded from the per-site filter rate.
    filtered = sum(1 for c in calls if [f for f in c.filters if f != "UncalibratedModel"])
    strand = sum(1 for c in calls if "StrandBiasReview" in c.filters)
    n_alt = decisions.get(DECISION_ALT_EXPLORATORY, 0) + decisions.get(DECISION_ALT_PASS, 0)
    # In a provisional (uncalibrated) run every exploratory call is expected to
    # fall back to ALT_UNCERTAIN, so an uncertainty threshold would always fail.
    provisional = cfg.run.policy_status != "VALIDATED"
    uncertain_threshold = None if provisional else 0.5
    filter_threshold = None if provisional else 0.9

    metrics = [
        evaluate_metric(
            "n_candidates", n, sample=sample_id, policy_id=cfg.policy_id(),
            provisional=provisional,
        ),
        evaluate_metric(
            "n_alt_exploratory", decisions.get(DECISION_ALT_EXPLORATORY, 0),
            sample=sample_id, policy_id=cfg.policy_id(), provisional=provisional,
        ),
        evaluate_metric(
            "n_alt_pass", decisions.get(DECISION_ALT_PASS, 0),
            sample=sample_id, policy_id=cfg.policy_id(), provisional=provisional,
        ),
        evaluate_metric(
            "n_uncertain", decisions.get(DECISION_ALT_UNCERTAIN, 0),
            sample=sample_id, threshold=uncertain_threshold, policy_id=cfg.policy_id(),
            provisional=provisional,
        ),
        evaluate_metric(
            "filter_rate", (filtered / n) if n else 0.0,
            sample=sample_id, threshold=filter_threshold, policy_id=cfg.policy_id(),
            provisional=provisional,
        ),
        evaluate_metric(
            "strand_bias_rate", (strand / n_alt) if n_alt else 0.0,
            sample=sample_id, threshold=0.2, policy_id=cfg.policy_id(),
            provisional=provisional,
        ),
    ]
    if mixed_infection_flag is not None:
        metrics.append(
            evaluate_metric(
                "mixed_infection_flag", mixed_infection_flag, sample=sample_id,
                policy_id=cfg.policy_id(), provisional=provisional,
            )
        )
    if n_polymorphic is not None:
        metrics.append(
            evaluate_metric(
                "n_polymorphic", n_polymorphic, sample=sample_id,
                policy_id=cfg.policy_id(), provisional=provisional,
            )
        )
    return metrics


def run_qc(sample_metrics: list[QCMetric], cfg: PfVarXConfig, *, n_samples: int) -> dict:
    all_metrics = list(sample_metrics)
    all_metrics.append(
        evaluate_metric(
            "n_samples", n_samples, policy_id=cfg.policy_id(),
            provisional=cfg.run.policy_status != "VALIDATED",
        )
    )
    summary = summarise_metrics(all_metrics)
    return {
        "metrics": [m.as_dict() for m in all_metrics],
        "summary": summary,
        "policy_id": cfg.policy_id(),
        "policy_status": cfg.run.policy_status,
    }

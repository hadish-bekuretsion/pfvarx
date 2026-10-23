"""QC metric evaluation and run aggregation tests."""

from __future__ import annotations

from pfvarx.qc.metrics import (
    STATUS_FAIL,
    STATUS_INFO,
    STATUS_PASS,
    STATUS_WARN,
    evaluate_metric,
    load_ontology,
    summarise_metrics,
)


def test_ontology_loads():
    onto = load_ontology()
    assert "metrics" in onto
    assert len(onto["metrics"]) > 0


def test_min_direction_pass_warn_fail():
    assert evaluate_metric("x", 10, threshold=10, direction="min").status == STATUS_PASS
    assert evaluate_metric("x", 9.5, threshold=10, direction="min").status == STATUS_WARN
    assert evaluate_metric("x", 5, threshold=10, direction="min").status == STATUS_FAIL


def test_max_direction_pass_warn_fail():
    assert evaluate_metric("x", 5, threshold=10, direction="max").status == STATUS_PASS
    assert evaluate_metric("x", 10.5, threshold=10, direction="max").status == STATUS_WARN
    assert evaluate_metric("x", 20, threshold=10, direction="max").status == STATUS_FAIL


def test_no_threshold_is_info():
    m = evaluate_metric("x", 5)
    assert m.status == STATUS_INFO


def test_summarise_counts():
    metrics = [
        evaluate_metric("x", 10, threshold=10, direction="min"),
        evaluate_metric("x", 5, threshold=10, direction="min"),
        evaluate_metric("x", 1),
    ]
    s = summarise_metrics(metrics)
    assert s["counts"][STATUS_PASS] == 1
    assert s["counts"][STATUS_FAIL] == 1
    assert s["counts"][STATUS_INFO] == 1


def test_provisional_flag_rendered():
    m = evaluate_metric("x", 1)
    assert m.provisional_or_validated == "PROVISIONAL"

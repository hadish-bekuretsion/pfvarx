"""Tests for the objective-aware dashboard and exporters."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from pfvarx.objectives import OBJECTIVES
from pfvarx.reporting.dashboard import (
    available_pages,
    discover_run,
    make_objective_page,
    page_func_for,
)
from pfvarx.reporting.report import (
    export_run_csv,
    export_run_figures,
    export_run_pdf,
    export_run_workbook,
)


def _write_manifest(root: Path, objectives: list[str]) -> dict:
    """Create a minimal but realistic results tree for the given objectives."""
    root.mkdir(parents=True, exist_ok=True)
    outcomes = []
    for obj_id in objectives:
        obj = OBJECTIVES[obj_id]
        out_dir = root / "objectives" / obj_id
        out_dir.mkdir(parents=True, exist_ok=True)
        artefacts = {}
        if obj.scope == "per_sample":
            sample_dir = out_dir / "S1"
            sample_dir.mkdir()
            table = sample_dir / "result.tsv"
            pd.DataFrame({"site": [1, 2, 3], "value": [0.1, 0.2, 0.3]}).to_csv(
                table, sep="\t", index=False
            )
            artefacts["result"] = str(table)
        else:
            tree = out_dir / "tree.nwk"
            tree.write_text("(S1:0.1,S2:0.2);\n")
            artefacts["tree"] = str(tree)
        outcomes.append(
            {
                "objective_id": obj_id,
                "label": obj.label,
                "status": "completed",
                "scope": obj.scope,
                "artefacts": artefacts,
                "samples": ["S1"],
                "duration_s": 0.5,
                "message": "ok",
                "error": None,
            }
        )
    manifest = {
        "run_id": "TEST",
        "root": str(root),
        "status": "completed",
        "started_at": "2024-01-01T00:00:00Z",
        "finished_at": "2024-01-01T00:00:01Z",
        "duration_s": 1.0,
        "objectives": objectives,
        "stages": ["call"],
        "routes": {"alignment_reuse": 2},
        "samples": ["S1", "S2"],
        "directories": [],
        "counts": {"completed": len(objectives), "failed": 0, "skipped": 0, "unavailable": 0},
        "outcomes": outcomes,
        "warnings": [],
        "errors": [],
    }
    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    (root / "variants" / "S1").mkdir(parents=True, exist_ok=True)
    (root / "variants" / "S1" / "S1.evidence.tsv").write_text("site\tdecision\n1\tALT\n")
    return manifest


def test_dashboard_shows_only_selected_objectives(tmp_path):
    _write_manifest(tmp_path, ["diversity"])
    inv = discover_run(tmp_path)
    pages = available_pages(inv)
    assert OBJECTIVES["diversity"].label in pages
    assert OBJECTIVES["relatedness"].label not in pages
    assert OBJECTIVES["drug_resistance"].label not in pages
    assert "Objectives Summary" in pages


def test_dashboard_adds_pages_for_multi_sample_objectives(tmp_path):
    _write_manifest(tmp_path, ["relatedness", "phylogenetics"])
    pages = available_pages(discover_run(tmp_path))
    assert OBJECTIVES["relatedness"].label in pages
    assert OBJECTIVES["phylogenetics"].label in pages
    assert OBJECTIVES["diversity"].label not in pages


def test_page_func_for_returns_objective_renderer(tmp_path):
    _write_manifest(tmp_path, ["diversity"])
    inv = discover_run(tmp_path)
    func = page_func_for(OBJECTIVES["diversity"].label, inv)
    assert callable(func)
    assert func.__name__ == "page_objective_diversity"


def test_make_objective_page_is_named_after_the_objective():
    assert make_objective_page("coi").__name__ == "page_objective_coi"


def test_discover_run_reads_manifest_objectives(tmp_path):
    _write_manifest(tmp_path, ["coi"])
    inv = discover_run(tmp_path)
    assert inv["manifest"]["objectives"] == ["coi"]


def test_dashboard_falls_back_to_plan_without_manifest(tmp_path):
    (tmp_path / "objective_plan.json").write_text(
        json.dumps({"selected": ["mixed_infection"], "stages": ["call"]})
    )
    pages = available_pages(discover_run(tmp_path))
    assert OBJECTIVES["mixed_infection"].label in pages


def test_export_csv_contains_only_selected_objectives(tmp_path):
    _write_manifest(tmp_path, ["diversity"])
    path = export_run_csv(tmp_path)
    df = pd.read_csv(path)
    assert set(df["objective"]) == {"diversity"}


def test_export_workbook_has_summary_sheet(tmp_path):
    _write_manifest(tmp_path, ["diversity", "relatedness"])
    path = export_run_workbook(tmp_path)
    sheets = pd.read_excel(path, sheet_name=None)
    assert "Summary" in sheets
    assert set(sheets["Summary"]["objective"]) == {"diversity", "relatedness"}


def test_export_pdf_is_a_valid_pdf(tmp_path):
    _write_manifest(tmp_path, ["diversity"])
    path = export_run_pdf(tmp_path)
    data = Path(path).read_bytes()
    assert data.startswith(b"%PDF")
    assert len(data) > 1000


def test_export_figures_writes_pngs(tmp_path):
    _write_manifest(tmp_path, ["diversity"])
    figures = export_run_figures(tmp_path)
    assert figures
    assert all(f.endswith(".png") for f in figures)


def test_export_without_manifest_raises_clear_error(tmp_path):
    with pytest.raises(FileNotFoundError) as exc:
        export_run_csv(tmp_path)
    assert "run_manifest.json" in str(exc.value)

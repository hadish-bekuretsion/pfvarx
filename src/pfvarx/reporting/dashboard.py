"""PfVarX Streamlit dashboard — objective-aware.

The dashboard is generated from the *actual* run, not from a fixed template:

* ``run_manifest.json`` is the single source of truth.  It records which
  objectives the user selected, which stages were planned, and an explicit
  outcome (``completed`` / ``failed`` / ``skipped`` / ``unavailable``) for each
  objective together with the artefacts it produced.
* Pages for objectives the user did **not** select are never created, so an
  unselected analysis cannot appear as an empty or misleading section.
* If no manifest exists yet (a run has not been executed) the dashboard falls
  back to ``objective_plan.json`` and then to the effective configuration, so
  that a pre-execution preview is still possible.

The dashboard is read-only: it never mutates results.  Every page renders the
release-status banner so that an exploratory call can never be mistaken for a
clinically validated one.

Run with::

    streamlit run -m pfvarx.reporting.dashboard -- --results results

or via the CLI::

    pfvarx report dashboard --results results
"""

from __future__ import annotations

import json
import os
from collections.abc import Callable
from pathlib import Path
from typing import Any

try:  # pragma: no cover - streamlit is an optional extra
    import streamlit as st
except Exception:  # pragma: no cover
    st = None  # type: ignore[assignment]

import pandas as pd

from pfvarx.objectives import OBJECTIVES
from pfvarx.release import overall_status, status_banner

BANNER = status_banner()

#: Human-readable marker for each outcome status.
STATUS_ICON = {
    "completed": "✅ completed",
    "failed": "❌ failed",
    "skipped": "⏭️ skipped",
    "unavailable": "⚠️ unavailable",
    "running": "⏳ running",
}

#: Pages that are meaningful for every run, in display order.
BASE_PAGES = [
    "Run Overview",
    "Objectives Summary",
    "Input & Validation",
    "Stage QC",
    "Variant Calling",
    "Variant Annotation",
    "Software & Environment",
    "Logs & Troubleshooting",
]


# --------------------------------------------------------------------------- #
# Data loading helpers
# --------------------------------------------------------------------------- #
def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _read_tsv(path: Path) -> pd.DataFrame | None:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def _selected_objectives(inv: dict[str, Any]) -> list[str]:
    """Objectives the user asked for, from the manifest or the plan."""
    manifest = inv.get("manifest") or {}
    if manifest.get("objectives"):
        return list(manifest["objectives"])
    plan = inv.get("objective_plan") or {}
    if plan.get("selected"):
        return list(plan["selected"])
    cfg = inv.get("run_config") or {}
    return list(cfg.get("objectives") or (cfg.get("downstream") or {}).get("enabled") or [])


def _outcomes(inv: dict[str, Any]) -> dict[str, dict[str, Any]]:
    manifest = inv.get("manifest") or {}
    return {o["objective_id"]: o for o in manifest.get("outcomes", [])}


def discover_run(results_dir: str | Path) -> dict[str, Any]:
    """Scan a results directory and return a structured, objective-aware inventory."""
    root = Path(results_dir)
    inventory: dict[str, Any] = {
        "root": str(root),
        "exists": root.exists(),
        "samples": [],
        "run_config": None,
        "run_qc": None,
        "manifest": None,
        "objective_plan": None,
        "routing_plan": None,
        "logs": [],
        "reports": [],
    }
    if not root.exists():
        return inventory

    inventory["run_config"] = _read_json(root / "config" / "effective_config.json") or _read_json(
        root / "effective_config.json"
    )
    inventory["manifest"] = _read_json(root / "run_manifest.json")
    inventory["objective_plan"] = _read_json(root / "objective_plan.json")
    inventory["routing_plan"] = _read_json(root / "routing_plan.json")
    inventory["run_qc"] = _read_json(root / "run_qc.json")

    # Evidence lives under variants/<sample>/ in the objective-driven layout, but
    # older runs wrote it directly under <sample>/.  Support both.
    sample_dirs: dict[str, Path] = {}
    variants_root = root / "variants"
    if variants_root.is_dir():
        for d in sorted(p for p in variants_root.iterdir() if p.is_dir()):
            sample_dirs[d.name] = d
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        if d.name in {"config", "logs", "qc", "objectives", "reports", "dashboard",
                      "alignments", "intermediate", "variants"}:
            continue
        sample_dirs.setdefault(d.name, d)

    for sample_id, sample_dir in sample_dirs.items():
        sample: dict[str, Any] = {"sample_id": sample_id, "dir": str(sample_dir)}
        for pattern, key in (
            ("*.pfvarx.vcf", "vcf"),
            ("*.pfvarx.vcf.gz", "vcf"),
            ("*.evidence.tsv", "evidence"),
            ("*.annotated.tsv", "annotation"),
            ("sample_qc.json", "sample_qc"),
            ("*.diversity.json", "diversity"),
            ("*.mixed_infection.json", "mixed_infection"),
        ):
            for hit in sorted(sample_dir.glob(pattern)):
                sample[key] = str(hit)
        inventory["samples"].append(sample)

    inventory["logs"] = [str(p) for p in sorted(root.rglob("*.log"))]
    inventory["reports"] = [
        str(p)
        for p in sorted(root.rglob("*"))
        if p.suffix.lower() in {".pdf", ".xlsx", ".csv", ".png", ".svg"}
    ]
    return inventory


def _load_evidence(sample: dict[str, Any]) -> pd.DataFrame | None:
    if "evidence" in sample:
        return _read_tsv(Path(sample["evidence"]))
    return None


# --------------------------------------------------------------------------- #
# Generic artefact renderer
# --------------------------------------------------------------------------- #
def _render_artefact(path_str: str) -> None:
    """Render one artefact according to its file type."""
    path = Path(path_str)
    if not path.exists():
        st.warning(f"Artefact recorded but missing on disk: {path}")
        return
    suffix = path.suffix.lower()
    if suffix == ".json":
        st.json(_read_json(path) or {}, expanded=False)
    elif suffix in {".tsv", ".csv"}:
        sep = "\t" if suffix == ".tsv" else ","
        try:
            df = pd.read_csv(path, sep=sep)
        except Exception as exc:  # pragma: no cover - defensive
            st.warning(f"Could not read {path.name}: {exc}")
            return
        st.dataframe(df, use_container_width=True, height=360)
        st.download_button(
            f"Download {path.name}",
            data=path.read_bytes(),
            file_name=path.name,
            mime="text/csv",
            key=f"dl-{path}",
        )
    elif suffix == ".nwk":
        st.code(path.read_text(), language="text")
    elif suffix in {".png", ".svg", ".jpg", ".jpeg"}:
        st.image(str(path), use_container_width=True)
    else:
        st.code("\n".join(path.read_text(errors="replace").splitlines()[:200]))


# --------------------------------------------------------------------------- #
# Base page renderers
# --------------------------------------------------------------------------- #
def page_overview(inv: dict[str, Any]) -> None:
    st.header("Run Overview")
    if not inv["exists"]:
        st.error(f"Results directory not found: {inv['root']}")
        return
    manifest = inv.get("manifest") or {}
    cfg = inv.get("run_config") or {}
    counts = manifest.get("counts") or {}

    cols = st.columns(5)
    cols[0].metric("Samples", len(inv["samples"]))
    cols[1].metric("Objectives", len(_selected_objectives(inv)))
    cols[2].metric("Completed", counts.get("completed", 0))
    cols[3].metric("Failed", counts.get("failed", 0))
    cols[4].metric("Duration (s)", manifest.get("duration_s", "n/a"))

    st.subheader("Run status")
    status = manifest.get("status", "not executed")
    st.info(f"Overall status: **{status}**")

    st.subheader("Release status")
    st.warning(BANNER)

    if manifest:
        st.subheader("Execution summary")
        st.write(
            {
                "run_id": manifest.get("run_id"),
                "started_at": manifest.get("started_at"),
                "finished_at": manifest.get("finished_at"),
                "routes": manifest.get("routes"),
                "stages": manifest.get("stages"),
            }
        )
    else:
        st.info(
            "No run_manifest.json found — this run has not been executed yet. "
            "The dashboard is showing the planned configuration only."
        )

    if cfg:
        with st.expander("Effective configuration"):
            st.json(cfg, expanded=False)

    if inv["samples"]:
        st.subheader("Samples")
        st.dataframe(
            pd.DataFrame(
                [
                    {
                        "sample_id": s["sample_id"],
                        "artefacts": ", ".join(k for k in s if k not in {"sample_id", "dir"}),
                    }
                    for s in inv["samples"]
                ]
            ),
            use_container_width=True,
        )


def page_objectives_summary(inv: dict[str, Any]) -> None:
    st.header("Objectives Summary")
    selected = _selected_objectives(inv)
    if not selected:
        st.info(
            "No objectives were selected for this run. Select objectives with "
            "`--objectives` or the `objectives:` key in the configuration."
        )
        return

    outcomes = _outcomes(inv)
    rows = []
    for obj_id in selected:
        meta = OBJECTIVES.get(obj_id)
        outcome = outcomes.get(obj_id, {})
        rows.append(
            {
                "objective": obj_id,
                "label": meta.label if meta else outcome.get("label", obj_id),
                "status": STATUS_ICON.get(outcome.get("status", "not executed"), "— not executed"),
                "scope": outcome.get("scope", meta.scope if meta else ""),
                "duration_s": outcome.get("duration_s", ""),
                "detail": outcome.get("message") or outcome.get("error") or "",
            }
        )
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.subheader("Export")
    st.caption(
        "Export the run-level results. Only the selected objectives are included; "
        "unselected analyses are never represented."
    )
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("Export CSV"):
            try:
                from pfvarx.reporting.exporters import export_calls_csv  # noqa: F401
                from pfvarx.reporting.report import export_run_csv

                path = export_run_csv(inv["root"])
                st.success(f"Written: {path}")
            except Exception as exc:
                st.error(f"CSV export failed: {exc}")
    with col2:
        if st.button("Export Excel"):
            try:
                from pfvarx.reporting.report import export_run_workbook

                path = export_run_workbook(inv["root"])
                st.success(f"Written: {path}")
            except Exception as exc:
                st.error(f"Excel export failed: {exc}")
    with col3:
        if st.button("Export PDF"):
            try:
                from pfvarx.reporting.report import export_run_pdf

                path = export_run_pdf(inv["root"])
                st.success(f"Written: {path}")
            except Exception as exc:
                st.error(f"PDF export failed: {exc}")


def page_input(inv: dict[str, Any]) -> None:
    st.header("Input & Validation")
    cfg = inv.get("run_config") or {}
    ref = cfg.get("reference") or {}
    st.subheader("Reference")
    st.write(
        {
            "assembly_id": ref.get("assembly_id"),
            "fasta": ref.get("fasta"),
            "sha256": ref.get("sha256"),
            "callable_mask_bed": ref.get("callable_mask_bed"),
        }
    )
    st.subheader("Routing plan")
    plan = inv.get("routing_plan")
    if plan:
        st.dataframe(pd.DataFrame(plan.get("samples", plan)), use_container_width=True)
    else:
        st.info("No routing_plan.json found. Run `pfvarx validate` to generate one.")


def page_stage_qc(inv: dict[str, Any]) -> None:
    st.header("Stage QC")
    run_qc = inv.get("run_qc")
    if run_qc:
        st.subheader("Run-level QC")
        st.json(run_qc, expanded=False)
    rows = []
    for sample in inv["samples"]:
        if "sample_qc" in sample:
            data = _read_json(Path(sample["sample_qc"])) or {}
            for metric in data.get("metrics", []):
                rows.append(metric)
    if rows:
        df = pd.DataFrame(rows)
        st.subheader("Per-sample metrics")
        if "status" in df.columns:
            st.dataframe(df.sort_values(["status", "metric"]), use_container_width=True)
        else:
            st.dataframe(df, use_container_width=True)
    elif not run_qc:
        st.info("No QC artefacts found.")


def page_calling(inv: dict[str, Any]) -> None:
    st.header("Variant Calling")
    if not inv["samples"]:
        st.info("No samples available.")
        return
    selected = st.selectbox("Sample", [s["sample_id"] for s in inv["samples"]])
    sample = next((s for s in inv["samples"] if s["sample_id"] == selected), None)
    if not sample:
        st.info("No samples available.")
        return
    df = _load_evidence(sample)
    if df is None or df.empty:
        st.info("No evidence table found for this sample.")
        return
    st.caption(f"{len(df)} candidate sites")
    if "decision" in df.columns:
        st.subheader("Decision distribution")
        st.bar_chart(df["decision"].value_counts())
    if "minority_bin" in df.columns:
        st.subheader("Minority bins")
        st.bar_chart(df["minority_bin"].value_counts())
    st.subheader("Evidence table")
    st.dataframe(df, use_container_width=True, height=420)


def page_annotation(inv: dict[str, Any]) -> None:
    st.header("Variant Annotation")
    found = False
    for sample in inv["samples"]:
        if "annotation" in sample:
            found = True
            st.subheader(sample["sample_id"])
            df = _read_tsv(Path(sample["annotation"]))
            if df is not None:
                st.dataframe(df, use_container_width=True, height=320)
    if not found:
        st.info("No annotated tables found for the selected objectives.")


def page_software(inv: dict[str, Any]) -> None:
    st.header("Software & Environment")
    from pfvarx.preprocess.tools import tool_matrix

    st.subheader("Tool matrix")
    st.dataframe(pd.DataFrame(tool_matrix()), use_container_width=True)
    st.subheader("PfVarX release status")
    st.json(overall_status(), expanded=False)
    env_path = Path(inv["root"]) / "environment.json"
    if env_path.exists():
        st.subheader("Recorded environment")
        st.json(_read_json(env_path) or {}, expanded=False)


def page_logs(inv: dict[str, Any]) -> None:
    st.header("Logs & Troubleshooting")
    manifest = inv.get("manifest") or {}
    warnings = manifest.get("warnings") or []
    errors = manifest.get("errors") or []
    if warnings:
        st.subheader("Warnings")
        for w in warnings:
            st.warning(w)
    if errors:
        st.subheader("Errors")
        for e in errors:
            st.error(e)
    if not inv["logs"]:
        st.info("No log files found.")
        return
    choice = st.selectbox("Log file", inv["logs"])
    text = Path(choice).read_text(errors="replace")
    st.code("\n".join(text.splitlines()[-400:]))
    st.subheader("Common remedies")
    st.markdown(
        "- **ReferenceMismatch** — confirm the FASTA md5/contig names match the BAM header.\n"
        "- **UncalibratedModel** — supply an empirical error profile or run `pfvarx calibrate`.\n"
        "- **UMIModelInvalid** — verify the UMI design and minimum family size.\n"
        "- **OutsideCallableRegion** — extend the callable mask or lower the callable fraction.\n"
        "- **LowEffectiveDepth** — increase coverage or relax the depth threshold for the mode.\n"
        "- **Objective unavailable** — read the objective's detail message; it names the "
        "missing prerequisite (GFF, metadata table, catalogue or a second sample)."
    )


# --------------------------------------------------------------------------- #
# Objective pages — one per selected objective, generated from the registry
# --------------------------------------------------------------------------- #
def make_objective_page(objective_id: str) -> Callable[[dict[str, Any]], None]:
    """Return a renderer for one objective, driven by the run manifest."""

    def _page(inv: dict[str, Any]) -> None:
        meta = OBJECTIVES.get(objective_id)
        label = meta.label if meta else objective_id
        st.header(label)
        if meta:
            st.caption(meta.summary)

        outcome = _outcomes(inv).get(objective_id)
        if outcome is None:
            st.warning(
                "This objective was selected but has no recorded outcome. The run may "
                "not have reached it; check the log for the reason."
            )
            return

        status = outcome.get("status", "unknown")
        st.markdown(f"**Status:** {STATUS_ICON.get(status, status)}  ·  "
                    f"**Scope:** {outcome.get('scope', '')}  ·  "
                    f"**Duration:** {outcome.get('duration_s', 'n/a')} s")
        if outcome.get("message"):
            st.info(outcome["message"])
        if outcome.get("error"):
            st.error(outcome["error"])

        artefacts = outcome.get("artefacts") or {}
        if not artefacts:
            st.warning("No artefacts were produced for this objective.")
            return

        for name, path in artefacts.items():
            st.subheader(name.replace("_", " ").title())
            _render_artefact(path)

    _page.__name__ = f"page_objective_{objective_id}"
    return _page


# --------------------------------------------------------------------------- #
# Page registry — computed from the actual run
# --------------------------------------------------------------------------- #
def available_pages(inv: dict[str, Any]) -> list[str]:
    """Pages relevant to this run, in display order.

    Objectives the user did not select are absent, so the navigation cannot
    show an unselected analysis.
    """
    selected = _selected_objectives(inv)
    manifest = inv.get("manifest") or {}
    stages = set(manifest.get("stages") or (inv.get("objective_plan") or {}).get("stages") or [])

    pages = ["Run Overview", "Objectives Summary"]
    pages += [OBJECTIVES[o].label for o in selected if o in OBJECTIVES]
    pages += ["Input & Validation", "Stage QC"]
    if "call" in stages or any("evidence" in s for s in inv["samples"]):
        pages.append("Variant Calling")
    if "annotate" in stages or any("annotation" in s for s in inv["samples"]):
        pages.append("Variant Annotation")
    pages += ["Software & Environment", "Logs & Troubleshooting"]
    # De-duplicate while preserving order.
    return list(dict.fromkeys(pages))


def page_func_for(title: str, inv: dict[str, Any]) -> Callable[[dict[str, Any]], None]:
    """Resolve a page title to its renderer for this particular run."""
    base = {
        "Run Overview": page_overview,
        "Objectives Summary": page_objectives_summary,
        "Input & Validation": page_input,
        "Stage QC": page_stage_qc,
        "Variant Calling": page_calling,
        "Variant Annotation": page_annotation,
        "Software & Environment": page_software,
        "Logs & Troubleshooting": page_logs,
    }
    if title in base:
        return base[title]
    for obj_id, meta in OBJECTIVES.items():
        if meta.label == title:
            return make_objective_page(obj_id)
    return page_overview


def main(results_dir: str | None = None) -> None:  # pragma: no cover - UI entry
    if st is None:
        raise RuntimeError("streamlit is not installed; install the 'report' extra")
    st.set_page_config(page_title="PfVarX", layout="wide", page_icon="🧬")
    st.sidebar.title("PfVarX")
    st.sidebar.caption(BANNER)
    results_dir = results_dir or os.environ.get("PFVARX_RESULTS", "results")
    results_dir = st.sidebar.text_input("Results directory", results_dir)
    inv = discover_run(results_dir)
    pages = available_pages(inv)
    page = st.sidebar.radio("Navigate", pages)
    st.title(f"PfVarX — {page}")
    st.caption(BANNER)
    page_func_for(page, inv)(inv)


if __name__ == "__main__":  # pragma: no cover
    import sys

    args = sys.argv[1:]
    target = None
    if "--results" in args:
        target = args[args.index("--results") + 1]
    main(target)

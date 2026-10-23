"""Objective-aware result export.

Every exporter here reads ``run_manifest.json`` and exports *only* the
objectives the user selected.  An unselected analysis is never represented,
not even as an empty sheet or section.

Three formats are supported:

* :func:`export_run_csv`      — a summary CSV plus one CSV per tabular artefact.
* :func:`export_run_workbook` — a single ``.xlsx`` with a summary sheet and one
  sheet per tabular artefact.
* :func:`export_run_pdf`      — a paginated PDF report with the objective
  status table and each objective's textual results.
* :func:`export_run_figures`  — PNG bar charts for numeric artefact columns.

All files are written under ``<root>/reports/`` and the path is returned.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from pfvarx.objectives import OBJECTIVES
from pfvarx.release import status_banner

STATUS_LABEL = {
    "completed": "completed",
    "failed": "failed",
    "skipped": "skipped",
    "unavailable": "unavailable",
    "running": "running",
}


def _load_manifest(root: str | Path) -> dict[str, Any]:
    path = Path(root) / "run_manifest.json"
    if not path.exists():
        raise FileNotFoundError(
            f"No run_manifest.json in {root}. Execute the run before exporting results."
        )
    return json.loads(path.read_text())


def _summary_frame(manifest: dict[str, Any]) -> pd.DataFrame:
    rows = []
    for outcome in manifest.get("outcomes", []):
        obj_id = outcome.get("objective_id", "")
        meta = OBJECTIVES.get(obj_id)
        rows.append(
            {
                "objective": obj_id,
                "label": meta.label if meta else outcome.get("label", obj_id),
                "status": outcome.get("status", ""),
                "scope": outcome.get("scope", ""),
                "duration_s": outcome.get("duration_s", ""),
                "n_artefacts": len(outcome.get("artefacts") or {}),
                "message": outcome.get("message", ""),
                "error": outcome.get("error", ""),
            }
        )
    return pd.DataFrame(rows)


def _tabular_artefacts(manifest: dict[str, Any]) -> dict[str, pd.DataFrame]:
    """Return ``{sheet_name: dataframe}`` for every tabular artefact."""
    frames: dict[str, pd.DataFrame] = {}
    for outcome in manifest.get("outcomes", []):
        obj_id = outcome.get("objective_id", "objective")
        for name, path_str in (outcome.get("artefacts") or {}).items():
            path = Path(path_str)
            if not path.exists():
                continue
            suffix = path.suffix.lower()
            try:
                if suffix == ".tsv":
                    df = pd.read_csv(path, sep="\t")
                elif suffix == ".csv":
                    df = pd.read_csv(path)
                elif suffix == ".json":
                    payload = json.loads(path.read_text())
                    if isinstance(payload, list) and payload and isinstance(payload[0], dict):
                        df = pd.DataFrame(payload)
                    elif isinstance(payload, dict):
                        df = pd.DataFrame(
                            [{"key": k, "value": json.dumps(v) if isinstance(v, (dict, list)) else v}
                             for k, v in payload.items()]
                        )
                    else:
                        continue
                else:
                    continue
            except Exception:
                continue
            if df.empty:
                continue
            sheet = f"{obj_id}.{name}"[:31]
            frames[sheet] = df
    return frames


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)


def export_run_csv(root: str | Path) -> str:
    """Write a summary CSV and one CSV per tabular artefact."""
    root = Path(root)
    manifest = _load_manifest(root)
    out_dir = root / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / "objective_summary.csv"
    _summary_frame(manifest).to_csv(summary_path, index=False)

    for sheet, df in _tabular_artefacts(manifest).items():
        df.to_csv(out_dir / f"{_safe(sheet)}.csv", index=False)
    return str(summary_path)


def export_run_workbook(root: str | Path) -> str:
    """Write a single Excel workbook with a summary sheet and artefact sheets."""
    root = Path(root)
    manifest = _load_manifest(root)
    out_dir = root / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "pfvarx_results.xlsx"

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        _summary_frame(manifest).to_excel(writer, sheet_name="Summary", index=False)
        for sheet, df in _tabular_artefacts(manifest).items():
            df.to_excel(writer, sheet_name=sheet, index=False)
    return str(path)


def export_run_pdf(root: str | Path) -> str:
    """Write a paginated PDF report covering only the selected objectives."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    root = Path(root)
    manifest = _load_manifest(root)
    out_dir = root / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "pfvarx_results.pdf"

    styles = getSampleStyleSheet()
    story: list[Any] = []
    story.append(Paragraph("PfVarX — Run Results", styles["Title"]))
    story.append(Paragraph(status_banner(), styles["Italic"]))
    story.append(Spacer(1, 6 * mm))
    story.append(
        Paragraph(
            f"Run ID: {manifest.get('run_id', 'n/a')} · "
            f"Status: {manifest.get('status', 'n/a')} · "
            f"Samples: {len(manifest.get('samples', []))} · "
            f"Duration: {manifest.get('duration_s', 0)} s",
            styles["Normal"],
        )
    )
    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph("Objective outcomes", styles["Heading2"]))

    table_rows = [["Objective", "Status", "Scope", "Duration (s)"]]
    for outcome in manifest.get("outcomes", []):
        obj_id = outcome.get("objective_id", "")
        meta = OBJECTIVES.get(obj_id)
        table_rows.append(
            [
                meta.label if meta else obj_id,
                STATUS_LABEL.get(outcome.get("status", ""), outcome.get("status", "")),
                outcome.get("scope", ""),
                str(outcome.get("duration_s", "")),
            ]
        )
    table = Table(table_rows, hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f4e79")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef3f8")]),
            ]
        )
    )
    story.append(table)

    for outcome in manifest.get("outcomes", []):
        story.append(PageBreak())
        obj_id = outcome.get("objective_id", "")
        meta = OBJECTIVES.get(obj_id)
        story.append(Paragraph(meta.label if meta else obj_id, styles["Heading1"]))
        if meta:
            story.append(Paragraph(meta.summary, styles["Italic"]))
        story.append(Spacer(1, 3 * mm))
        story.append(
            Paragraph(
                f"Status: {STATUS_LABEL.get(outcome.get('status', ''), outcome.get('status', ''))}",
                styles["Normal"],
            )
        )
        if outcome.get("message"):
            story.append(Paragraph(f"Detail: {outcome['message']}", styles["Normal"]))
        if outcome.get("error"):
            story.append(Paragraph(f"Error: {outcome['error']}", styles["Normal"]))
        for name, path_str in (outcome.get("artefacts") or {}).items():
            story.append(Spacer(1, 3 * mm))
            story.append(Paragraph(name.replace("_", " ").title(), styles["Heading3"]))
            p = Path(path_str)
            if not p.exists():
                story.append(Paragraph("Artefact missing on disk.", styles["Normal"]))
                continue
            if p.suffix.lower() == ".json":
                text = p.read_text()
                story.append(Paragraph(f"<font size=7><pre>{text[:4000]}</pre></font>", styles["Normal"]))
            elif p.suffix.lower() in {".tsv", ".csv"}:
                df = pd.read_csv(p, sep="\t" if p.suffix.lower() == ".tsv" else ",").head(25)
                data = [list(df.columns)] + df.astype(str).values.tolist()
                t = Table(data, hAlign="LEFT")
                t.setStyle(
                    TableStyle(
                        [
                            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#d9e2f3")),
                            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                            ("FONTSIZE", (0, 0), (-1, -1), 7),
                            ("GRID", (0, 0), (-1, -1), 0.2, colors.grey),
                        ]
                    )
                )
                story.append(t)
            elif p.suffix.lower() == ".nwk":
                story.append(Paragraph(f"<font size=7>{p.read_text()[:2000]}</font>", styles["Normal"]))

    doc = SimpleDocTemplate(str(path), pagesize=A4)
    doc.build(story)
    return str(path)


def export_run_figures(root: str | Path) -> list[str]:
    """Write one PNG bar chart per numeric artefact column."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = Path(root)
    manifest = _load_manifest(root)
    out_dir = root / "reports" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []

    for sheet, df in _tabular_artefacts(manifest).items():
        numeric = df.select_dtypes("number")
        if numeric.empty:
            continue
        for column in numeric.columns[:4]:
            fig, ax = plt.subplots(figsize=(6, 3))
            df[column].plot(kind="bar", ax=ax)
            ax.set_title(f"{sheet} — {column}")
            ax.set_ylabel(column)
            fig.tight_layout()
            path = out_dir / f"{_safe(sheet)}_{_safe(str(column))}.png"
            fig.savefig(path, dpi=120)
            plt.close(fig)
            written.append(str(path))
    return written

"""Report exporters: CSV, XLSX, PNG figures and PDF summary.

All exporters are optional-dependency tolerant: if a library is missing the
exporter raises a clear error rather than silently producing a partial report.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from pfvarx.calling.policy import SiteCall, summarise_calls
from pfvarx.config.schema import PfVarXConfig


def _ensure_parent(path: str) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def calls_to_dataframe(calls: list[SiteCall]) -> pd.DataFrame:
    return pd.DataFrame([c.as_dict() for c in calls])


def export_calls_csv(calls: list[SiteCall], path: str) -> str:
    p = _ensure_parent(path)
    calls_to_dataframe(calls).to_csv(p, index=False)
    return str(p)


def export_workbook(
    calls: list[SiteCall],
    qc: dict,
    path: str,
    *,
    metadata: dict | None = None,
) -> str:
    """Write a multi-sheet Excel workbook (calls, decisions, QC, provenance)."""
    p = _ensure_parent(path)
    df = calls_to_dataframe(calls)
    summary = summarise_calls(calls)
    decisions = pd.DataFrame(
        [{"decision": k, "count": v} for k, v in summary["decisions"].items()]
    )
    bins = pd.DataFrame(
        [{"minority_bin": k, "count": v} for k, v in summary["minority_bins"].items()]
    )
    filters = pd.DataFrame(
        [{"filter": k, "count": v} for k, v in summary["filter_counts"].items()]
    )
    qc_df = pd.DataFrame(qc.get("metrics", []))
    prov = pd.DataFrame([metadata or {}])
    with pd.ExcelWriter(p, engine="openpyxl") as xw:
        df.to_excel(xw, sheet_name="calls", index=False)
        decisions.to_excel(xw, sheet_name="decisions", index=False)
        bins.to_excel(xw, sheet_name="minority_bins", index=False)
        filters.to_excel(xw, sheet_name="filters", index=False)
        qc_df.to_excel(xw, sheet_name="qc", index=False)
        prov.to_excel(xw, sheet_name="provenance", index=False)
    return str(p)


def export_figures(calls: list[SiteCall], out_dir: str) -> list[str]:
    """Render QC figures as PNG.  Returns the list of written paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    summary = summarise_calls(calls)

    # Decision counts
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = list(summary["decisions"].keys())
    values = list(summary["decisions"].values())
    ax.bar(labels, values, color="#3b6ea5")
    ax.set_ylabel("Sites")
    ax.set_title("Decision categories")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    path = out / "decisions.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    written.append(str(path))

    # Allele-fraction distribution
    afs = [c.evidence.alt_fraction for c in calls if c.evidence.alt_count > 0]
    if afs:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(afs, bins=25, color="#c0504d", edgecolor="white")
        ax.set_xlabel("Alternate allele fraction")
        ax.set_ylabel("Sites")
        ax.set_title("Allele-fraction distribution")
        fig.tight_layout()
        path = out / "allele_fraction.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(str(path))

    # Depth vs AF
    if calls:
        fig, ax = plt.subplots(figsize=(7, 4))
        xs = [c.evidence.informative_depth for c in calls]
        ys = [c.evidence.alt_fraction for c in calls]
        colors = ["#4f81bd" if c.is_alt else "#bfbfbf" for c in calls]
        ax.scatter(xs, ys, c=colors, s=12, alpha=0.7)
        ax.set_xlabel("Informative depth")
        ax.set_ylabel("Alternate allele fraction")
        ax.set_title("Depth vs allele fraction")
        fig.tight_layout()
        path = out / "depth_vs_af.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(str(path))
    return written


def export_pdf(
    calls: list[SiteCall],
    qc: dict,
    path: str,
    *,
    cfg: PfVarXConfig,
    sample_id: str,
    figures: list[str] | None = None,
) -> str:
    """Write a concise PDF report with the mandatory research-use banner."""
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import Image, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

    p = _ensure_parent(path)
    styles = getSampleStyleSheet()
    doc = SimpleDocTemplate(str(p), pagesize=A4, leftMargin=2 * cm, rightMargin=2 * cm)
    story = [
        Paragraph("PfVarX Variant Analysis Report", styles["Title"]),
        Paragraph(
            "<b>RESEARCH USE ONLY.</b> This report is not a clinical device and "
            f"must not be used for diagnosis or treatment decisions. Validation "
            f"status: {cfg.run.policy_status}.",
            styles["BodyText"],
        ),
        Spacer(1, 12),
        Paragraph(f"Sample: {sample_id}", styles["Heading2"]),
        Paragraph(f"Policy ID: {cfg.policy_id()}", styles["BodyText"]),
        Paragraph(f"Sequencing mode: {cfg.run.mode}", styles["BodyText"]),
        Spacer(1, 12),
    ]

    summary = summarise_calls(calls)
    data = [["Decision", "Sites"]] + [[k, str(v)] for k, v in summary["decisions"].items()]
    table = Table(data, colWidths=[8 * cm, 4 * cm])
    table.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.5, "#999999"), ("BACKGROUND", (0, 0), (-1, 0), "#dce6f1")]))
    story += [Paragraph("Decision summary", styles["Heading2"]), table, Spacer(1, 12)]

    qc_summary = qc.get("summary", {})
    story.append(Paragraph(f"Overall QC status: {qc_summary.get('overall_status', 'INFO')}", styles["BodyText"]))
    story.append(Spacer(1, 12))

    for fig in figures or []:
        try:
            story.append(Image(fig, width=15 * cm, height=8.5 * cm))
            story.append(Spacer(1, 8))
        except Exception:
            continue

    doc.build(story)
    return str(p)

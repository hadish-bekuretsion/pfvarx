"""Reporting exporter tests: CSV, XLSX, figures, PDF."""

from __future__ import annotations

import os

from pfvarx.calling.policy import SiteCall, SiteEvidence
from pfvarx.reporting.exporters import (
    calls_to_dataframe,
    export_calls_csv,
    export_figures,
    export_pdf,
    export_workbook,
)


def _calls():
    out = []
    for pos, (ref, alt) in enumerate([(50, 50), (90, 10), (100, 0)], start=1):
        depth = ref + alt
        ev = SiteEvidence(
            sample_id="S1", chrom="c", pos=pos, ref="A", alt="T", var_type="SNP",
            ref_count=ref, alt_count=alt, depth=depth, informative_depth=depth,
            alt_fraction=(alt / depth) if depth else 0.0,
        )
        out.append(SiteCall(evidence=ev, decision="ALT_DETECTED_EXPLORATORY"))
    return out


def test_calls_to_dataframe():
    df = calls_to_dataframe(_calls())
    assert len(df) == 3
    assert "decision" in df.columns


def test_export_csv(tmp_path):
    p = export_calls_csv(_calls(), str(tmp_path / "calls.csv"))
    assert os.path.exists(p)


def test_export_workbook(tmp_path):
    p = export_workbook(_calls(), {"metrics": []}, str(tmp_path / "r.xlsx"))
    assert os.path.exists(p)


def test_export_figures(tmp_path):
    figs = export_figures(_calls(), str(tmp_path / "figs"))
    assert len(figs) >= 1
    for f in figs:
        assert os.path.exists(f)


def test_export_pdf(tmp_path):
    from tests.conftest import make_config

    p = export_pdf(
        _calls(), {"metrics": []}, str(tmp_path / "r.pdf"),
        cfg=make_config(), sample_id="S1",
    )
    assert os.path.exists(p)
    assert os.path.getsize(p) > 0

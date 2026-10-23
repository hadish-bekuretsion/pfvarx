"""Reporting layer: tabular, workbook, figure, PDF and dashboard outputs.

All reporting artefacts carry the mandatory research-use banner.  No report may
present an uncalibrated likelihood ratio as a posterior probability or encode an
uncalibrated metric in the VCF QUAL column.
"""

from pfvarx.reporting.exporters import (
    calls_to_dataframe,
    export_calls_csv,
    export_figures,
    export_pdf,
    export_workbook,
)

__all__ = [
    "calls_to_dataframe",
    "export_calls_csv",
    "export_figures",
    "export_pdf",
    "export_workbook",
]

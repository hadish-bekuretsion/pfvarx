"""Spatial and temporal aggregation of sample-level summaries.

Metadata is supplied in the samplesheet; this module groups sample summaries by
arbitrary metadata columns and produces tidy tables suitable for plotting.
"""

from __future__ import annotations

import pandas as pd


def aggregate_by_metadata(
    sample_summaries: list[dict],
    metadata: pd.DataFrame,
    *,
    group_by: list[str],
    value_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Merge sample summaries with metadata and aggregate by group."""
    if not group_by:
        raise ValueError("group_by must contain at least one column")
    df = pd.DataFrame(sample_summaries)
    if "sample_id" not in df.columns:
        raise ValueError("sample_summaries must include a sample_id column")
    merged = df.merge(metadata, on="sample_id", how="left", suffixes=("", "_meta"))
    missing = [g for g in group_by if g not in merged.columns]
    if missing:
        raise ValueError(f"grouping columns not present in metadata: {missing}")
    if value_columns is None:
        value_columns = [c for c in df.columns if c != "sample_id" and pd.api.types.is_numeric_dtype(df[c])]
    agg = merged.groupby(group_by, dropna=False)[value_columns].agg(["mean", "median", "count"])
    agg.columns = ["_".join(map(str, c)).strip("_") for c in agg.columns]
    return agg.reset_index()


def temporal_series(sample_summaries: list[dict], metadata: pd.DataFrame, *, date_column: str = "collection_date") -> pd.DataFrame:
    """Return sample summaries ordered by a collection date column."""
    df = pd.DataFrame(sample_summaries).merge(metadata, on="sample_id", how="left", suffixes=("", "_meta"))
    if date_column not in df.columns:
        raise ValueError(f"date column {date_column!r} not found in metadata")
    df[date_column] = pd.to_datetime(df[date_column], errors="coerce")
    return df.sort_values(date_column)

"""PyRanges1-backed interval normalization and overlap helpers."""

from __future__ import annotations

import pandas as pd
import pyranges1 as pr


SEQ_REGION_COL = "seq_region_name"
START_COL = "seq_region_start"
END_COL = "seq_region_end"
STRAND_COL = "seq_region_strand"


def _normalise_strand_symbol(value: object) -> str:
    if value in ("+", "+1", "1", 1):
        return "+"
    if value in ("-", "-1", -1):
        return "-"
    return "."


def to_pyranges_frame(  # pylint: disable=too-many-arguments
    df: pd.DataFrame,
    *,
    feature_id_col: str,
    feature_prefix: str,
    source_name: str = "",
    source_role: str = "",
    source_class: str = "",
) -> pd.DataFrame:
    """Return a PyRanges-ready dataframe using 0-based half-open coordinates."""
    columns = [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        f"{feature_prefix}_id",
        f"{feature_prefix}_row_index",
        "source_name",
        "source_role",
        "source_class",
    ]
    if df.empty:
        return pd.DataFrame(columns=columns)

    out = df.copy()
    for col in (SEQ_REGION_COL, START_COL, END_COL, STRAND_COL, feature_id_col):
        if col not in out.columns:
            out[col] = pd.NA

    out[f"{feature_prefix}_row_index"] = out.index
    out[f"{feature_prefix}_id"] = out[feature_id_col].astype("string")
    out["Chromosome"] = out[SEQ_REGION_COL].astype("string").fillna("")
    starts = pd.to_numeric(out[START_COL], errors="coerce").fillna(0).astype(int)
    ends = pd.to_numeric(out[END_COL], errors="coerce").fillna(0).astype(int)
    out["Start"] = (starts - 1).clip(lower=0)
    out["End"] = ends.clip(lower=0)
    out["Strand"] = out[STRAND_COL].map(_normalise_strand_symbol)
    out["source_name"] = source_name
    out["source_role"] = source_role
    out["source_class"] = source_class
    out = out[(out["Chromosome"] != "") & (out["End"] > out["Start"])]
    return out[columns]


def to_pyranges(df: pd.DataFrame) -> pr.PyRanges:
    """Wrap a PyRanges-ready dataframe as a pyranges1 object."""
    if df.empty:
        return pr.PyRanges(
            pd.DataFrame(columns=["Chromosome", "Start", "End", "Strand"])
        )
    return pr.PyRanges(df)


def overlap_matches(
    query: pd.DataFrame,
    target: pd.DataFrame,
    *,
    same_strand: bool = True,
    suffix: str = "_target",
) -> pd.DataFrame:
    """Join two PyRanges-ready frames and report overlap bp per local match."""
    if query.empty or target.empty:
        return pd.DataFrame()
    strand_behavior = "same" if same_strand else "ignore"
    joined = to_pyranges(query).join_overlaps(
        to_pyranges(target),
        strand_behavior=strand_behavior,
        report_overlap_column="Overlap",
        suffix=suffix,
        preserve_input_order=False,
    )
    return pd.DataFrame(joined)


def best_overlap_by_query(
    query: pd.DataFrame,
    target: pd.DataFrame,
    *,
    query_id_col: str,
    target_id_col: str,
    same_strand: bool = True,
) -> pd.DataFrame:
    """Return best target and overlap count for each query interval."""
    columns = [
        query_id_col,
        target_id_col,
        "target_row_index",
        "overlap_bp",
        "overlap_count",
    ]
    matches = overlap_matches(query, target, same_strand=same_strand)
    if matches.empty:
        return pd.DataFrame(columns=columns)

    target_row_col = "target_row_index"
    if target_row_col not in matches.columns and "target_row_index_target" in matches:
        target_row_col = "target_row_index_target"

    counts = (
        matches.groupby(query_id_col, dropna=False)
        .size()
        .rename("overlap_count")
        .reset_index()
    )
    ranked = matches.sort_values(
        [query_id_col, "Overlap", target_id_col], ascending=[True, False, True]
    )
    best = ranked.drop_duplicates(query_id_col, keep="first")[
        [query_id_col, target_id_col, target_row_col, "Overlap"]
    ].rename(
        columns={
            target_row_col: "target_row_index",
            "Overlap": "overlap_bp",
        }
    )
    return best.merge(counts, on=query_id_col, how="left")[columns]

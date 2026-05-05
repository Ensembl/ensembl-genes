from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from dataclasses import dataclass
from typing import Iterable, Optional
from urllib.parse import unquote

import pandas as pd

from .db import (
    DBParams,
    connect,
    extract_all_genes,
    extract_all_transcripts,
    extract_all_translations,
    list_seq_regions,
)
from .io import write_manifest
from .loci import build_loci
from .utils import ensure_dir, make_run_id


PRIORITY_RANK = {"P0": 0, "P1": 1, "P2": 2, "P3": 3, "P4": 4}
CONFIDENCE_RANK = {"high": 0, "medium": 1, "low": 2, "unknown": 3}


@dataclass(frozen=True)
class AuditPaths:
    output_dir: str
    run_id: str
    fmt: str = "tsv"

    @property
    def root(self) -> str:
        return os.path.join(self.output_dir, self.run_id)

    @property
    def extract_dir(self) -> str:
        return os.path.join(self.root, "extract")

    @property
    def loci_dir(self) -> str:
        return os.path.join(self.root, "loci")

    @property
    def observability_dir(self) -> str:
        return os.path.join(self.root, "observability")


def _paths_from_args(args: argparse.Namespace) -> AuditPaths:
    run_dir = getattr(args, "run_dir", None)
    if run_dir:
        run_dir = os.path.abspath(run_dir)
        return AuditPaths(os.path.dirname(run_dir), os.path.basename(run_dir), args.format)
    output_dir = getattr(args, "output_dir", None)
    run_id = getattr(args, "run_id", None)
    if not output_dir or not run_id:
        raise ValueError("Provide either --run_dir or both --output_dir and --run_id")
    return AuditPaths(output_dir, run_id, args.format)


def _empty(columns: Iterable[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(columns))


def _read_table(path: str, fmt: Optional[str] = None) -> pd.DataFrame:
    if fmt is None:
        _, ext = os.path.splitext(path)
        fmt = ext.lstrip(".").lower()
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    if fmt == "parquet":
        return pd.read_parquet(path)
    sep = "\t" if fmt in ("tsv", "txt") else ","
    try:
        return pd.read_csv(path, sep=sep, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _load_named_table(base_dir: str, base: str, prefer_fmt: str, required: bool = True) -> pd.DataFrame:
    for fmt in (prefer_fmt, "tsv", "csv", "parquet"):
        path = os.path.join(base_dir, f"{base}.{fmt}")
        if os.path.exists(path):
            return _read_table(path, fmt)
    if required:
        raise FileNotFoundError(f"Could not find {base} in {base_dir}")
    return pd.DataFrame()


def _write_table(df: pd.DataFrame, path: str, fmt: str) -> None:
    ensure_dir(os.path.dirname(path))
    if fmt == "parquet":
        df.to_parquet(path, index=False)
    elif fmt == "csv":
        df.to_csv(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False)


def _span_bp(start: object, end: object) -> int:
    try:
        return max(0, int(end) - int(start) + 1)
    except (TypeError, ValueError):
        return 0


def _overlap_bp(a_start: object, a_end: object, b_start: object, b_end: object) -> int:
    try:
        start = max(int(a_start), int(b_start))
        end = min(int(a_end), int(b_end))
    except (TypeError, ValueError):
        return 0
    return max(0, end - start + 1)


def _normalise_strand(value: object) -> int:
    if pd.isna(value):
        return 0
    if value in ("+", "+1", "1", 1):
        return 1
    if value in ("-", "-1", -1):
        return -1
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _string_value(value: object, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    text = str(value)
    if text.lower() == "nan":
        return default
    return text


def _numeric_value(value: object, default: float = 0.0) -> float:
    if value is None or pd.isna(value):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _int_value(value: object, default: int = 0) -> int:
    return int(_numeric_value(value, float(default)))


def _normalise_gene_table(df: pd.DataFrame, source_kind: str) -> pd.DataFrame:
    cols = [
        "gene_id",
        "stable_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "biotype",
        "canonical_transcript_id",
        "logic_name",
    ]
    out = df.copy()
    for col in cols:
        if col not in out.columns:
            out[col] = pd.NA
    out["source_kind"] = source_kind
    if not out.empty:
        out["seq_region_name"] = out["seq_region_name"].astype(str)
        out["seq_region_strand"] = out["seq_region_strand"].map(_normalise_strand)
        out["seq_region_start"] = out["seq_region_start"].map(lambda x: _int_value(x, 0))
        out["seq_region_end"] = out["seq_region_end"].map(lambda x: _int_value(x, 0))
    return out[cols + ["source_kind"]]


def _normalise_transcript_table(df: pd.DataFrame, source_kind: str) -> pd.DataFrame:
    cols = [
        "transcript_id",
        "gene_id",
        "stable_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "biotype",
        "logic_name",
    ]
    out = df.copy()
    for col in cols:
        if col not in out.columns:
            out[col] = pd.NA
    out["source_kind"] = source_kind
    if not out.empty:
        out["seq_region_name"] = out["seq_region_name"].astype(str)
        out["seq_region_strand"] = out["seq_region_strand"].map(_normalise_strand)
        out["seq_region_start"] = out["seq_region_start"].map(lambda x: _int_value(x, 0))
        out["seq_region_end"] = out["seq_region_end"].map(lambda x: _int_value(x, 0))
    return out[cols + ["source_kind"]]


def _translation_ids(df: pd.DataFrame) -> set[int]:
    if df.empty or "transcript_id" not in df.columns:
        return set()
    ids: set[int] = set()
    for value in df["transcript_id"].dropna().tolist():
        try:
            ids.add(int(value))
        except (TypeError, ValueError):
            continue
    return ids


def _gene_support_by_id(tx_df: pd.DataFrame, tr_df: pd.DataFrame) -> pd.DataFrame:
    columns = ["gene_id", "transcript_count", "coding_transcript_count"]
    if tx_df.empty or "gene_id" not in tx_df.columns:
        return _empty(columns)
    tx = tx_df.copy()
    tr_ids = _translation_ids(tr_df)
    tx["has_translation"] = tx["transcript_id"].map(lambda x: int(_int_value(x, -1) in tr_ids))
    grouped = (
        tx.groupby("gene_id", dropna=False)
        .agg(transcript_count=("transcript_id", "count"), coding_transcript_count=("has_translation", "sum"))
        .reset_index()
    )
    return grouped[columns]


def _best_overlap(query: pd.Series, targets: pd.DataFrame, same_strand: bool = True) -> tuple[Optional[pd.Series], int]:
    if targets.empty:
        return None, 0
    seq = _string_value(query.get("seq_region_name"))
    strand = _normalise_strand(query.get("seq_region_strand"))
    subset = targets[targets["seq_region_name"].astype(str) == seq]
    if same_strand:
        subset = subset[subset["seq_region_strand"].map(_normalise_strand) == strand]
    if subset.empty:
        return None, 0
    best_row: Optional[pd.Series] = None
    best_bp = 0
    for _, target in subset.iterrows():
        bp = _overlap_bp(
            query.get("seq_region_start"),
            query.get("seq_region_end"),
            target.get("seq_region_start"),
            target.get("seq_region_end"),
        )
        if bp > best_bp:
            best_bp = bp
            best_row = target
    return best_row, best_bp


def _overlap_count(query: pd.Series, targets: pd.DataFrame, same_strand: Optional[bool]) -> int:
    if targets.empty:
        return 0
    seq = _string_value(query.get("seq_region_name"))
    strand = _normalise_strand(query.get("seq_region_strand"))
    subset = targets[targets["seq_region_name"].astype(str) == seq]
    if same_strand is True:
        subset = subset[subset["seq_region_strand"].map(_normalise_strand) == strand]
    elif same_strand is False:
        subset = subset[subset["seq_region_strand"].map(_normalise_strand) == -strand]
    count = 0
    for _, target in subset.iterrows():
        if _overlap_bp(
            query.get("seq_region_start"),
            query.get("seq_region_end"),
            target.get("seq_region_start"),
            target.get("seq_region_end"),
        ):
            count += 1
    return count


def _coverage_class(coverage: float) -> str:
    if coverage >= 0.95:
        return "span_represented_high"
    if coverage >= 0.80:
        return "span_represented_near"
    if coverage > 0:
        return "span_represented_partial"
    return "no_span_match"


def _evidence_fate_class(
    layer_row: pd.Series,
    best_core: Optional[pd.Series],
    coverage: float,
    layer_has_coding: bool,
    same_strand_core_count: int,
    opposite_strand_core_count: int,
    collapse_group_size: int,
) -> tuple[str, str, str, str]:
    if best_core is None or same_strand_core_count == 0:
        if opposite_strand_core_count > 0:
            return (
                "orphan_evidence",
                "opposite_strand_core_only",
                "P2" if layer_has_coding else "P3",
                "inspect_browser",
            )
        return (
            "orphan_evidence",
            "no_core_gene_built",
            "P1" if layer_has_coding else "P2",
            "try_candidate_rescue",
        )

    core_biotype = _string_value(best_core.get("biotype"), "unknown")
    layer_biotype = _string_value(layer_row.get("biotype"), "unknown")

    if layer_has_coding and core_biotype != "protein_coding":
        return (
            "represented_partial" if coverage < 0.95 else "represented_near",
            "coding_evidence_not_protein_coding",
            "P1",
            "check_biotype_assignment",
        )
    if coverage < 0.50:
        return ("represented_partial", "core_span_underrepresents_evidence", "P2", "inspect_browser")
    if coverage < 0.80:
        return ("represented_partial", "partial_core_representation", "P2", "inspect_browser")
    if collapse_group_size > 1 and layer_biotype == core_biotype:
        return ("collapsed_into_core", "possible_transcript_or_model_collapse", "P3", "inspect_collapse_group")
    if coverage >= 0.95:
        return ("represented_exact_or_near_span", "represented_span_only", "P4", "no_action")
    return ("represented_near", "represented_span_only", "P4", "no_action")


def audit_evidence_fate(
    core_genes: pd.DataFrame,
    layer_genes: pd.DataFrame,
    core_tx: pd.DataFrame,
    layer_tx: pd.DataFrame,
    core_tr: pd.DataFrame,
    layer_tr: pd.DataFrame,
) -> pd.DataFrame:
    core = _normalise_gene_table(core_genes, "core")
    layer = _normalise_gene_table(layer_genes, "layer")
    core_tx_norm = _normalise_transcript_table(core_tx, "core")
    layer_tx_norm = _normalise_transcript_table(layer_tx, "layer")

    layer_support = _gene_support_by_id(layer_tx_norm, layer_tr)
    core_support = _gene_support_by_id(core_tx_norm, core_tr)
    if not layer.empty:
        layer = layer.merge(layer_support, on="gene_id", how="left")
    if not core.empty:
        core = core.merge(core_support, on="gene_id", how="left")
    for df in (layer, core):
        if "transcript_count" not in df.columns:
            df["transcript_count"] = 0
        if "coding_transcript_count" not in df.columns:
            df["coding_transcript_count"] = 0
        df["transcript_count"] = pd.to_numeric(df["transcript_count"], errors="coerce").fillna(0).astype(int)
        df["coding_transcript_count"] = (
            pd.to_numeric(df["coding_transcript_count"], errors="coerce").fillna(0).astype(int)
        )

    preliminary: list[dict] = []
    for idx, row in layer.iterrows():
        best_core, overlap = _best_overlap(row, core, same_strand=True)
        layer_span = _span_bp(row.get("seq_region_start"), row.get("seq_region_end"))
        coverage = overlap / layer_span if layer_span else 0.0
        same_count = _overlap_count(row, core, same_strand=True)
        opp_count = _overlap_count(row, core, same_strand=False)
        preliminary.append(
            {
                "_row_index": idx,
                "layer_gene_id": row.get("gene_id"),
                "layer_stable_id": _string_value(row.get("stable_id")),
                "seq_region_name": _string_value(row.get("seq_region_name")),
                "seq_region_start": _int_value(row.get("seq_region_start")),
                "seq_region_end": _int_value(row.get("seq_region_end")),
                "seq_region_strand": _normalise_strand(row.get("seq_region_strand")),
                "layer_biotype": _string_value(row.get("biotype"), "unknown"),
                "layer_logic_name": _string_value(row.get("logic_name"), "unknown"),
                "layer_transcript_count": _int_value(row.get("transcript_count")),
                "layer_coding_transcript_count": _int_value(row.get("coding_transcript_count")),
                "best_core_gene_id": best_core.get("gene_id") if best_core is not None else pd.NA,
                "best_core_stable_id": _string_value(best_core.get("stable_id")) if best_core is not None else "",
                "best_core_biotype": (
                    _string_value(best_core.get("biotype"), "unknown") if best_core is not None else ""
                ),
                "best_core_transcript_count": (
                    _int_value(best_core.get("transcript_count")) if best_core is not None else 0
                ),
                "best_core_coding_transcript_count": (
                    _int_value(best_core.get("coding_transcript_count")) if best_core is not None else 0
                ),
                "same_strand_core_count": same_count,
                "opposite_strand_core_count": opp_count,
                "overlap_bp": overlap,
                "layer_span_bp": layer_span,
                "layer_span_coverage_by_core": coverage,
                "representation_class": _coverage_class(coverage),
                "_layer_has_coding": bool(_int_value(row.get("coding_transcript_count")) > 0),
            }
        )

    if not preliminary:
        return _empty(
            [
                "layer_gene_id",
                "layer_stable_id",
                "seq_region_name",
                "seq_region_start",
                "seq_region_end",
                "seq_region_strand",
                "layer_biotype",
                "layer_logic_name",
                "fate_class",
                "failure_class",
                "review_priority",
                "suggested_action",
            ]
        )

    out = pd.DataFrame(preliminary)
    best_counts = out["best_core_gene_id"].dropna().astype(str).value_counts().to_dict()
    fate_rows = []
    core_by_id = {str(r.get("gene_id")): r for _, r in core.iterrows()}
    layer_by_index = {idx: row for idx, row in layer.iterrows()}
    for _, row in out.iterrows():
        best_id = _string_value(row.get("best_core_gene_id"))
        best_core = core_by_id.get(best_id)
        layer_row = layer_by_index.get(row["_row_index"], pd.Series(dtype=object))
        collapse_size = best_counts.get(best_id, 0) if best_id else 0
        fate, failure, priority, action = _evidence_fate_class(
            layer_row,
            best_core,
            float(row["layer_span_coverage_by_core"]),
            bool(row["_layer_has_coding"]),
            int(row["same_strand_core_count"]),
            int(row["opposite_strand_core_count"]),
            collapse_size,
        )
        fate_rows.append((fate, failure, priority, action, collapse_size))

    out["fate_class"] = [x[0] for x in fate_rows]
    out["failure_class"] = [x[1] for x in fate_rows]
    out["review_priority"] = [x[2] for x in fate_rows]
    out["suggested_action"] = [x[3] for x in fate_rows]
    out["collapse_group_size"] = [x[4] for x in fate_rows]
    out["review_reason"] = out.apply(
        lambda r: f"{r['failure_class']}; coverage={float(r['layer_span_coverage_by_core']):.2f}; "
        f"logic_name={r['layer_logic_name']}",
        axis=1,
    )
    return out.drop(columns=["_row_index", "_layer_has_coding"]).sort_values(
        ["review_priority", "failure_class", "seq_region_name", "seq_region_start"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
    )


def _normalise_expected_genes(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "expected_gene_id",
        "expected_source",
        "reference_stable_id",
        "symbol",
        "biotype",
        "orthogroup_id",
        "busco_id",
        "expected_copy_number",
        "confidence",
    ]
    out = df.copy()
    for col in cols:
        if col not in out.columns:
            out[col] = pd.NA
    if "expected_gene_id" not in df.columns and "reference_stable_id" in df.columns:
        out["expected_gene_id"] = out["reference_stable_id"]
    out["confidence"] = out["confidence"].map(lambda x: _string_value(x, "unknown").lower())
    out["expected_copy_number"] = out["expected_copy_number"].map(lambda x: _int_value(x, 1) or 1)
    return out[cols]


def _normalise_expected_projections(df: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "target_seq_region_name": "seq_region_name",
        "target_start": "seq_region_start",
        "target_end": "seq_region_end",
        "target_strand": "seq_region_strand",
        "start": "seq_region_start",
        "end": "seq_region_end",
        "strand": "seq_region_strand",
    }
    out = df.rename(columns={k: v for k, v in rename.items() if k in df.columns}).copy()
    cols = [
        "expected_gene_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "projection_status",
        "projection_identity",
        "projection_coverage",
        "assembly_gap_overlap_bp",
        "repeat_overlap_bp",
    ]
    for col in cols:
        if col not in out.columns:
            out[col] = pd.NA
    if not out.empty:
        out["seq_region_name"] = out["seq_region_name"].map(lambda x: _string_value(x))
        out["seq_region_start"] = out["seq_region_start"].map(lambda x: _int_value(x))
        out["seq_region_end"] = out["seq_region_end"].map(lambda x: _int_value(x))
        out["seq_region_strand"] = out["seq_region_strand"].map(_normalise_strand)
        out["projection_status"] = out["projection_status"].map(lambda x: _string_value(x, "mapped").lower())
        out["projection_coverage"] = out["projection_coverage"].map(lambda x: _numeric_value(x, 0.0))
        out["projection_identity"] = out["projection_identity"].map(lambda x: _numeric_value(x, 0.0))
        out["assembly_gap_overlap_bp"] = out["assembly_gap_overlap_bp"].map(lambda x: _int_value(x, 0))
        out["repeat_overlap_bp"] = out["repeat_overlap_bp"].map(lambda x: _int_value(x, 0))
    return out[cols]


def _expected_presence_class(
    expected: pd.Series,
    projection: pd.Series,
    best_core: Optional[pd.Series],
    best_layer: Optional[pd.Series],
    core_overlap: int,
    layer_overlap: int,
    core_overlap_count: int,
    core_expected_count_for_best_core: int,
) -> tuple[str, str, str, str, str]:
    confidence = _string_value(expected.get("confidence"), "unknown").lower()
    expected_biotype = _string_value(expected.get("biotype"), "unknown")
    projection_status = _string_value(projection.get("projection_status"), "mapped").lower()
    projection_coverage = _numeric_value(projection.get("projection_coverage"), 0.0)
    gap_bp = _int_value(projection.get("assembly_gap_overlap_bp"), 0)
    expected_span = _span_bp(projection.get("seq_region_start"), projection.get("seq_region_end"))
    core_cov = core_overlap / expected_span if expected_span else 0.0

    if projection_status in ("unmapped", "missing") or expected_span == 0:
        if confidence == "high":
            return ("unresolved", "projection_unmapped", "P2", "add_or_check_projection", "unassessable")
        return ("unresolved", "projection_unmapped", "P3", "add_or_check_projection", "unassessable")

    if gap_bp > 0 and best_core is None:
        return ("assembly_limited", "assembly_gap_limited", "P2", "check_assembly_gap", "gap_limited")

    if best_core is None:
        if best_layer is not None or layer_overlap > 0:
            priority = "P1" if confidence == "high" else "P2"
            return ("missing_with_evidence", "no_core_gene_built", priority, "try_candidate_rescue", "intact")
        if projection_coverage >= 0.50 or projection_status in ("mapped", "partial"):
            priority = "P1" if confidence == "high" else "P2"
            return ("projection_only", "projected_expected_gene_not_built", priority, "inspect_browser", "intact")
        return ("missing_no_evidence", "evidence_absent_on_assembly", "P3", "mark_possible_true_loss", "uncertain")

    core_biotype = _string_value(best_core.get("biotype"), "unknown")
    if core_overlap_count > 1:
        return ("split", "split_gene_candidate", "P2", "inspect_split_candidate", "intact")
    if core_expected_count_for_best_core > 1:
        return ("fused", "fused_gene_candidate", "P2", "inspect_fusion_candidate", "intact")
    if expected_biotype == "protein_coding" and core_biotype not in ("protein_coding", "unknown"):
        return ("present_wrong_biotype", "wrong_biotype", "P1", "check_biotype_assignment", "intact")
    if core_cov >= 0.80:
        return ("present_clean", "represented_span_only", "P4", "no_action", "intact")
    if core_cov > 0:
        priority = "P1" if confidence == "high" else "P2"
        return ("present_degraded", "core_span_underrepresents_expected_gene", priority, "inspect_browser", "intact")
    return ("unresolved", "no_span_match_after_overlap_selection", "P3", "inspect_browser", "uncertain")


def audit_expected_presence(
    expected_genes: pd.DataFrame,
    expected_projections: pd.DataFrame,
    core_genes: pd.DataFrame,
    layer_genes: pd.DataFrame,
) -> pd.DataFrame:
    expected = _normalise_expected_genes(expected_genes)
    projections = _normalise_expected_projections(expected_projections)
    core = _normalise_gene_table(core_genes, "core")
    layer = _normalise_gene_table(layer_genes, "layer")

    if expected.empty:
        return _empty(
            [
                "expected_gene_id",
                "presence_class",
                "failure_class",
                "review_priority",
                "suggested_action",
            ]
        )

    if projections.empty:
        projections = expected[["expected_gene_id"]].copy()
        projections["seq_region_name"] = ""
        projections["seq_region_start"] = 0
        projections["seq_region_end"] = 0
        projections["seq_region_strand"] = 0
        projections["projection_status"] = "unmapped"
        projections["projection_identity"] = 0.0
        projections["projection_coverage"] = 0.0
        projections["assembly_gap_overlap_bp"] = 0
        projections["repeat_overlap_bp"] = 0
    else:
        projections = expected[["expected_gene_id"]].merge(projections, on="expected_gene_id", how="left")
        missing_projection = projections["projection_status"].isna()
        projections.loc[missing_projection, "seq_region_name"] = ""
        projections.loc[missing_projection, "seq_region_start"] = 0
        projections.loc[missing_projection, "seq_region_end"] = 0
        projections.loc[missing_projection, "seq_region_strand"] = 0
        projections.loc[missing_projection, "projection_status"] = "unmapped"
        projections.loc[missing_projection, "projection_identity"] = 0.0
        projections.loc[missing_projection, "projection_coverage"] = 0.0
        projections.loc[missing_projection, "assembly_gap_overlap_bp"] = 0
        projections.loc[missing_projection, "repeat_overlap_bp"] = 0

    merged = projections.merge(expected, on="expected_gene_id", how="left")
    preliminary: list[dict] = []
    best_core_by_expected: dict[str, str] = {}
    for _, row in merged.iterrows():
        query = pd.Series(
            {
                "seq_region_name": row.get("seq_region_name"),
                "seq_region_start": row.get("seq_region_start"),
                "seq_region_end": row.get("seq_region_end"),
                "seq_region_strand": row.get("seq_region_strand"),
            }
        )
        best_core, core_overlap = _best_overlap(query, core, same_strand=True)
        best_layer, layer_overlap = _best_overlap(query, layer, same_strand=True)
        core_count = _overlap_count(query, core, same_strand=True)
        expected_id = _string_value(row.get("expected_gene_id"))
        best_core_id = _string_value(best_core.get("gene_id")) if best_core is not None else ""
        if expected_id and best_core_id:
            best_core_by_expected[expected_id] = best_core_id
        preliminary.append(
            {
                "expected_gene_id": expected_id,
                "expected_source": _string_value(row.get("expected_source"), "unknown"),
                "reference_stable_id": _string_value(row.get("reference_stable_id")),
                "symbol": _string_value(row.get("symbol")),
                "expected_biotype": _string_value(row.get("biotype"), "unknown"),
                "orthogroup_id": _string_value(row.get("orthogroup_id")),
                "busco_id": _string_value(row.get("busco_id")),
                "expected_copy_number": _int_value(row.get("expected_copy_number"), 1),
                "confidence": _string_value(row.get("confidence"), "unknown").lower(),
                "seq_region_name": _string_value(row.get("seq_region_name")),
                "seq_region_start": _int_value(row.get("seq_region_start")),
                "seq_region_end": _int_value(row.get("seq_region_end")),
                "seq_region_strand": _normalise_strand(row.get("seq_region_strand")),
                "projection_status": _string_value(row.get("projection_status"), "mapped"),
                "projection_identity": _numeric_value(row.get("projection_identity")),
                "projection_coverage": _numeric_value(row.get("projection_coverage")),
                "assembly_gap_overlap_bp": _int_value(row.get("assembly_gap_overlap_bp")),
                "repeat_overlap_bp": _int_value(row.get("repeat_overlap_bp")),
                "best_core_gene_id": best_core.get("gene_id") if best_core is not None else pd.NA,
                "best_core_stable_id": _string_value(best_core.get("stable_id")) if best_core is not None else "",
                "best_core_biotype": (
                    _string_value(best_core.get("biotype"), "unknown") if best_core is not None else ""
                ),
                "best_layer_gene_id": best_layer.get("gene_id") if best_layer is not None else pd.NA,
                "best_layer_stable_id": _string_value(best_layer.get("stable_id")) if best_layer is not None else "",
                "best_layer_logic_name": (
                    _string_value(best_layer.get("logic_name"), "unknown") if best_layer is not None else ""
                ),
                "core_overlap_bp": core_overlap,
                "layer_overlap_bp": layer_overlap,
                "core_overlap_count": core_count,
                "_best_core_obj": best_core,
                "_best_layer_obj": best_layer,
                "_source_row": row,
            }
        )

    core_expected_counts = pd.Series(list(best_core_by_expected.values())).value_counts().to_dict()
    classified: list[dict] = []
    for row in preliminary:
        best_core_id = _string_value(row.get("best_core_gene_id"))
        presence, failure, priority, action, assembly_state = _expected_presence_class(
            row["_source_row"],
            row["_source_row"],
            row["_best_core_obj"],
            row["_best_layer_obj"],
            int(row["core_overlap_bp"]),
            int(row["layer_overlap_bp"]),
            int(row["core_overlap_count"]),
            int(core_expected_counts.get(best_core_id, 0)),
        )
        row = {k: v for k, v in row.items() if not k.startswith("_")}
        expected_span = _span_bp(row["seq_region_start"], row["seq_region_end"])
        row["expected_span_bp"] = expected_span
        row["expected_span_coverage_by_core"] = (
            float(row["core_overlap_bp"]) / float(expected_span) if expected_span else 0.0
        )
        row["presence_class"] = presence
        row["failure_class"] = failure
        row["assembly_state"] = assembly_state
        row["review_priority"] = priority
        row["suggested_action"] = action
        row["review_reason"] = (
            f"{failure}; confidence={row['confidence']}; "
            f"projection={row['projection_status']}:{row['projection_coverage']:.2f}; "
            f"core_coverage={row['expected_span_coverage_by_core']:.2f}"
        )
        classified.append(row)

    out = pd.DataFrame(classified)
    return out.sort_values(
        ["review_priority", "confidence", "seq_region_name", "seq_region_start"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9)
        if s.name == "review_priority"
        else (s.map(CONFIDENCE_RANK).fillna(9) if s.name == "confidence" else s),
    )


def build_or_load_loci(
    paths: AuditPaths,
    core_genes: pd.DataFrame,
    layer_genes: pd.DataFrame,
    gap_bp: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    loci = _load_named_table(paths.loci_dir, "loci.strict", paths.fmt, required=False)
    mapping = _load_named_table(paths.loci_dir, "gene_to_locus", paths.fmt, required=False)
    if not loci.empty and not mapping.empty:
        return loci, mapping
    loci, _, mapping = build_loci(core_genes, layer_genes, gap_bp=gap_bp)
    return loci, mapping


def build_audit_loci(loci: pd.DataFrame, expected_presence: pd.DataFrame, review_loci: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "locus_id",
        "seq_region_name",
        "seq_region_strand",
        "locus_start",
        "locus_end",
        "locus_length",
        "core_gene_count",
        "layer_gene_count",
        "expected_gene_count",
        "review_locus_count",
        "p1_review_locus_count",
    ]
    if loci.empty:
        return _empty(columns)

    out = loci.copy()
    if "locus_id" not in out.columns:
        out["locus_id"] = (
            out["seq_region_name"].astype(str)
            + ":"
            + out["seq_region_strand"].astype(str)
            + ":"
            + out["locus_start"].astype(int).astype(str)
            + ":"
            + out["locus_end"].astype(int).astype(str)
        )

    expected_counts = []
    review_counts = []
    p1_counts = []
    for _, locus in out.iterrows():
        expected_counts.append(_interval_row_count(locus, expected_presence, "expected_gene_id"))
        review_counts.append(_interval_row_count(locus, review_loci, "subject_id"))
        if review_loci.empty:
            p1_counts.append(0)
        else:
            p1_counts.append(
                _interval_row_count(locus, review_loci[review_loci["review_priority"] == "P1"], "subject_id")
            )
    out["expected_gene_count"] = expected_counts
    out["review_locus_count"] = review_counts
    out["p1_review_locus_count"] = p1_counts
    for col in columns:
        if col not in out.columns:
            out[col] = 0
    return out[columns]


def _interval_row_count(locus: pd.Series, rows: pd.DataFrame, id_col: str) -> int:
    if rows.empty or id_col not in rows.columns:
        return 0
    seq = _string_value(locus.get("seq_region_name"))
    strand = _normalise_strand(locus.get("seq_region_strand"))
    start = _int_value(locus.get("locus_start"))
    end = _int_value(locus.get("locus_end"))
    subset = rows[rows["seq_region_name"].astype(str) == seq]
    if "seq_region_strand" in subset.columns:
        subset = subset[subset["seq_region_strand"].map(_normalise_strand) == strand]
    ids = set()
    for _, row in subset.iterrows():
        if _overlap_bp(start, end, row.get("seq_region_start"), row.get("seq_region_end")):
            ids.add(_string_value(row.get(id_col), default=str(len(ids))))
    return len(ids)


def _failure_summary(evidence_fate: pd.DataFrame, expected_presence: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    if not evidence_fate.empty:
        grouped = (
            evidence_fate.groupby(["failure_class", "review_priority"], dropna=False)
            .size()
            .reset_index(name="count")
        )
        for _, row in grouped.iterrows():
            rows.append(
                {
                    "audit_track": "evidence_fate",
                    "class": row["failure_class"],
                    "review_priority": row["review_priority"],
                    "count": int(row["count"]),
                }
            )
    if not expected_presence.empty:
        grouped = (
            expected_presence.groupby(["presence_class", "failure_class", "review_priority"], dropna=False)
            .size()
            .reset_index(name="count")
        )
        for _, row in grouped.iterrows():
            rows.append(
                {
                    "audit_track": "expected_presence",
                    "class": f"{row['presence_class']}|{row['failure_class']}",
                    "review_priority": row["review_priority"],
                    "count": int(row["count"]),
                }
            )
    return pd.DataFrame(rows, columns=["audit_track", "class", "review_priority", "count"])


def build_source_profile(evidence_fate: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "layer_logic_name",
        "layer_biotype",
        "n_layer_models",
        "n_coding_models",
        "n_orphan_models",
        "n_p1_models",
        "median_core_coverage",
        "top_failure_class",
    ]
    if evidence_fate.empty:
        return _empty(columns)
    rows = []
    grouped = evidence_fate.groupby(["layer_logic_name", "layer_biotype"], dropna=False)
    for (logic_name, biotype), group in grouped:
        top_failure = group["failure_class"].value_counts().index[0] if not group.empty else ""
        rows.append(
            {
                "layer_logic_name": logic_name,
                "layer_biotype": biotype,
                "n_layer_models": int(len(group)),
                "n_coding_models": int((group["layer_coding_transcript_count"] > 0).sum()),
                "n_orphan_models": int((group["fate_class"] == "orphan_evidence").sum()),
                "n_p1_models": int((group["review_priority"] == "P1").sum()),
                "median_core_coverage": float(group["layer_span_coverage_by_core"].median()),
                "top_failure_class": top_failure,
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["n_p1_models", "n_orphan_models", "n_layer_models"], ascending=[False, False, False]
    )


def build_expected_source_profile(expected_presence: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "expected_source",
        "confidence",
        "n_expected_genes",
        "n_present_clean",
        "n_missing_with_evidence",
        "n_projection_only",
        "n_assembly_limited",
        "n_p1_genes",
        "top_presence_class",
    ]
    if expected_presence.empty:
        return _empty(columns)
    rows = []
    grouped = expected_presence.groupby(["expected_source", "confidence"], dropna=False)
    for (source, confidence), group in grouped:
        top_presence = group["presence_class"].value_counts().index[0] if not group.empty else ""
        rows.append(
            {
                "expected_source": source,
                "confidence": confidence,
                "n_expected_genes": int(len(group)),
                "n_present_clean": int((group["presence_class"] == "present_clean").sum()),
                "n_missing_with_evidence": int((group["presence_class"] == "missing_with_evidence").sum()),
                "n_projection_only": int((group["presence_class"] == "projection_only").sum()),
                "n_assembly_limited": int((group["presence_class"] == "assembly_limited").sum()),
                "n_p1_genes": int((group["review_priority"] == "P1").sum()),
                "top_presence_class": top_presence,
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["n_p1_genes", "n_missing_with_evidence", "n_expected_genes"], ascending=[False, False, False]
    )


def build_busco_crosswalk(expected_presence: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "busco_id",
        "expected_gene_id",
        "confidence",
        "projection_status",
        "presence_class",
        "failure_class",
        "review_priority",
        "best_core_stable_id",
        "expected_span_coverage_by_core",
        "suggested_action",
    ]
    if expected_presence.empty or "busco_id" not in expected_presence.columns:
        return _empty(columns)
    busco = expected_presence[expected_presence["busco_id"].map(lambda x: _string_value(x) != "")].copy()
    if busco.empty:
        return _empty(columns)
    return busco[columns].sort_values(
        ["review_priority", "busco_id"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
    )


def _has_busco(row: pd.Series) -> bool:
    return bool(_string_value(row.get("busco_id")))


def _is_high_confidence(row: pd.Series) -> bool:
    return _string_value(row.get("confidence"), "unknown").lower() == "high"


def _is_present_any(row: pd.Series) -> bool:
    return _string_value(row.get("presence_class")) in {
        "present_clean",
        "present_degraded",
        "present_wrong_canonical",
        "present_wrong_biotype",
        "split",
        "fused",
    }


def _is_missing_actionable(row: pd.Series) -> bool:
    return _string_value(row.get("presence_class")) in {
        "missing_with_evidence",
        "projection_only",
        "present_degraded",
        "present_wrong_biotype",
        "split",
        "fused",
    }


def _panel_rows(expected_presence: pd.DataFrame) -> list[tuple[str, str, pd.DataFrame]]:
    if expected_presence.empty:
        return []
    panels: list[tuple[str, str, pd.DataFrame]] = [
        ("all_expected", "All expected genes in the catalogue", expected_presence),
    ]
    high = expected_presence[expected_presence.apply(_is_high_confidence, axis=1)]
    if not high.empty:
        panels.append(("high_confidence", "High-confidence expected genes", high))
        high_non_busco = high[~high.apply(_has_busco, axis=1)]
        if not high_non_busco.empty:
            panels.append(
                (
                    "high_confidence_non_busco",
                    "High-confidence expected genes outside BUSCO",
                    high_non_busco,
                )
            )
    busco = expected_presence[expected_presence.apply(_has_busco, axis=1)]
    if not busco.empty:
        panels.append(("busco_linked", "Expected genes linked to BUSCO markers", busco))
    same_species = expected_presence[
        expected_presence["expected_source"].astype(str).isin(["prior_ensembl", "same_species", "previous_ensembl"])
    ]
    if not same_species.empty:
        panels.append(("same_species_reference", "Previous or same-species expected genes", same_species))
    for source, group in expected_presence.groupby("expected_source", dropna=False):
        source_name = _string_value(source, "unknown")
        panels.append((f"source:{source_name}", f"Expected genes from source {source_name}", group))
    return panels


def build_completeness_profile(expected_presence: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "panel_id",
        "panel_description",
        "n_expected",
        "n_present_any",
        "present_any_fraction",
        "n_present_clean",
        "present_clean_fraction",
        "n_missing_with_evidence",
        "missing_with_evidence_fraction",
        "n_projection_only",
        "projection_only_fraction",
        "n_assembly_limited",
        "assembly_limited_fraction",
        "n_actionable_loss",
        "actionable_loss_fraction",
        "n_p1",
        "p1_fraction",
    ]
    rows = []
    for panel_id, description, panel in _panel_rows(expected_presence):
        n_expected = len(panel)
        if n_expected == 0:
            continue
        n_present_any = int(panel.apply(_is_present_any, axis=1).sum())
        n_present_clean = int((panel["presence_class"] == "present_clean").sum())
        n_missing_with_evidence = int((panel["presence_class"] == "missing_with_evidence").sum())
        n_projection_only = int((panel["presence_class"] == "projection_only").sum())
        n_assembly_limited = int((panel["presence_class"] == "assembly_limited").sum())
        n_actionable_loss = int(panel.apply(_is_missing_actionable, axis=1).sum())
        n_p1 = int((panel["review_priority"] == "P1").sum())
        rows.append(
            {
                "panel_id": panel_id,
                "panel_description": description,
                "n_expected": n_expected,
                "n_present_any": n_present_any,
                "present_any_fraction": n_present_any / n_expected,
                "n_present_clean": n_present_clean,
                "present_clean_fraction": n_present_clean / n_expected,
                "n_missing_with_evidence": n_missing_with_evidence,
                "missing_with_evidence_fraction": n_missing_with_evidence / n_expected,
                "n_projection_only": n_projection_only,
                "projection_only_fraction": n_projection_only / n_expected,
                "n_assembly_limited": n_assembly_limited,
                "assembly_limited_fraction": n_assembly_limited / n_expected,
                "n_actionable_loss": n_actionable_loss,
                "actionable_loss_fraction": n_actionable_loss / n_expected,
                "n_p1": n_p1,
                "p1_fraction": n_p1 / n_expected,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def build_non_busco_high_confidence_losses(expected_presence: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "expected_gene_id",
        "expected_source",
        "reference_stable_id",
        "symbol",
        "expected_biotype",
        "orthogroup_id",
        "confidence",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "projection_status",
        "presence_class",
        "failure_class",
        "review_priority",
        "best_core_stable_id",
        "best_layer_stable_id",
        "best_layer_logic_name",
        "expected_span_coverage_by_core",
        "suggested_action",
        "review_reason",
    ]
    if expected_presence.empty:
        return _empty(columns)
    mask = (
        expected_presence.apply(_is_high_confidence, axis=1)
        & ~expected_presence.apply(_has_busco, axis=1)
        & expected_presence.apply(_is_missing_actionable, axis=1)
    )
    out = expected_presence[mask].copy()
    if out.empty:
        return _empty(columns)
    for col in columns:
        if col not in out.columns:
            out[col] = ""
    return out[columns].sort_values(
        ["review_priority", "expected_source", "seq_region_name", "seq_region_start"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
    )


def build_busco_proxy_calibration(expected_presence: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "comparison",
        "panel_id",
        "n_expected",
        "present_clean_fraction",
        "present_any_fraction",
        "actionable_loss_fraction",
        "assembly_limited_fraction",
        "interpretation",
    ]
    completeness = build_completeness_profile(expected_presence)
    if completeness.empty:
        return _empty(columns)

    wanted = completeness[
        completeness["panel_id"].isin(["busco_linked", "high_confidence_non_busco", "high_confidence"])
    ].copy()
    rows = []
    for _, row in wanted.iterrows():
        panel_id = row["panel_id"]
        if panel_id == "busco_linked":
            interpretation = "BUSCO-linked expected genes; sentinel conserved coding panel"
        elif panel_id == "high_confidence_non_busco":
            interpretation = "High-confidence expected genes that BUSCO does not cover"
        else:
            interpretation = "All high-confidence expected genes"
        rows.append(
            {
                "comparison": "busco_proxy_context",
                "panel_id": panel_id,
                "n_expected": int(row["n_expected"]),
                "present_clean_fraction": float(row["present_clean_fraction"]),
                "present_any_fraction": float(row["present_any_fraction"]),
                "actionable_loss_fraction": float(row["actionable_loss_fraction"]),
                "assembly_limited_fraction": float(row["assembly_limited_fraction"]),
                "interpretation": interpretation,
            }
        )

    by_panel = {row["panel_id"]: row for _, row in completeness.iterrows()}
    if "busco_linked" in by_panel and "high_confidence_non_busco" in by_panel:
        busco = by_panel["busco_linked"]
        non_busco = by_panel["high_confidence_non_busco"]
        delta = float(non_busco["actionable_loss_fraction"]) - float(busco["actionable_loss_fraction"])
        rows.append(
            {
                "comparison": "non_busco_minus_busco_actionable_loss",
                "panel_id": "high_confidence_non_busco_vs_busco_linked",
                "n_expected": int(non_busco["n_expected"]),
                "present_clean_fraction": pd.NA,
                "present_any_fraction": pd.NA,
                "actionable_loss_fraction": delta,
                "assembly_limited_fraction": pd.NA,
                "interpretation": (
                    "Positive values mean high-confidence non-BUSCO genes "
                    "are failing more often than BUSCO-linked genes"
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def build_copy_number_audit(expected_presence: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "copy_group_id",
        "group_source",
        "expected_copy_number",
        "observed_core_copy_number",
        "copy_number_class",
        "n_expected_rows",
        "review_priority",
        "member_expected_gene_ids",
    ]
    if expected_presence.empty:
        return _empty(columns)

    rows = []
    work = expected_presence.copy()
    work["copy_group_id"] = work.apply(_copy_group_id, axis=1)
    for group_id, group in work.groupby("copy_group_id", dropna=False):
        if not _string_value(group_id):
            continue
        expected_copy_number = max(1, int(group["expected_copy_number"].max()))
        observed_core_ids = {
            _string_value(x) for x in group["best_core_stable_id"].tolist() if _string_value(x)
        }
        observed = len(observed_core_ids)
        if observed == 0:
            klass = "missing_all_copies"
            priority = "P1" if (group["confidence"] == "high").any() else "P2"
        elif observed < expected_copy_number:
            klass = "collapsed_copy_number"
            priority = "P1" if (group["confidence"] == "high").any() else "P2"
        elif observed > expected_copy_number:
            klass = "expanded_copy_number"
            priority = "P2"
        else:
            klass = "copy_number_as_expected"
            priority = "P4"
        rows.append(
            {
                "copy_group_id": group_id,
                "group_source": _copy_group_source(group.iloc[0]),
                "expected_copy_number": expected_copy_number,
                "observed_core_copy_number": observed,
                "copy_number_class": klass,
                "n_expected_rows": int(len(group)),
                "review_priority": priority,
                "member_expected_gene_ids": ",".join(sorted(group["expected_gene_id"].astype(str).tolist())),
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["review_priority", "copy_number_class", "copy_group_id"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
    )


def _copy_group_id(row: pd.Series) -> str:
    for col in ("orthogroup_id", "busco_id"):
        value = _string_value(row.get(col))
        if value:
            return value
    return _string_value(row.get("expected_gene_id"))


def _copy_group_source(row: pd.Series) -> str:
    if _string_value(row.get("orthogroup_id")):
        return "orthogroup_id"
    if _string_value(row.get("busco_id")):
        return "busco_id"
    return "expected_gene_id"


def build_biotype_transition(evidence_fate: pd.DataFrame, expected_presence: pd.DataFrame) -> pd.DataFrame:
    columns = ["transition_track", "source_biotype", "core_biotype", "count", "p1_count"]
    rows = []
    if not evidence_fate.empty:
        grouped = evidence_fate.groupby(["layer_biotype", "best_core_biotype"], dropna=False)
        for (source_bt, core_bt), group in grouped:
            rows.append(
                {
                    "transition_track": "layer_to_core",
                    "source_biotype": source_bt,
                    "core_biotype": core_bt or "no_core",
                    "count": int(len(group)),
                    "p1_count": int((group["review_priority"] == "P1").sum()),
                }
            )
    if not expected_presence.empty:
        grouped = expected_presence.groupby(["expected_biotype", "best_core_biotype"], dropna=False)
        for (source_bt, core_bt), group in grouped:
            rows.append(
                {
                    "transition_track": "expected_to_core",
                    "source_biotype": source_bt,
                    "core_biotype": core_bt or "no_core",
                    "count": int(len(group)),
                    "p1_count": int((group["review_priority"] == "P1").sum()),
                }
            )
    if not rows:
        return _empty(columns)
    return pd.DataFrame(rows, columns=columns).sort_values(["p1_count", "count"], ascending=[False, False])


GFFCOMPARE_TMAP_COLUMNS = [
    "reference_gene_id",
    "reference_transcript_id",
    "class_code",
    "query_gene_id",
    "query_transcript_id",
    "num_exons",
    "fpkm",
    "tpm",
    "coverage",
    "length",
    "major_isoform_id",
    "reference_match_length",
]

GFFCOMPARE_CLASS_RANK = {
    "=": 0,
    "c": 1,
    "k": 2,
    "j": 2,
    "m": 3,
    "n": 3,
    "e": 4,
    "o": 4,
    "s": 4,
    "x": 4,
    "i": 4,
    "y": 5,
    "p": 5,
    "r": 5,
    "u": 6,
    "-": 7,
    "": 7,
}


def read_gffcompare_tmap(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
    if df.empty:
        return _empty(GFFCOMPARE_TMAP_COLUMNS)
    df = df.iloc[:, : len(GFFCOMPARE_TMAP_COLUMNS)]
    df.columns = GFFCOMPARE_TMAP_COLUMNS[: len(df.columns)]
    for col in GFFCOMPARE_TMAP_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA
    first = _string_value(df.iloc[0].get("reference_gene_id")).lower()
    if first in {"reference_gene_id", "ref_gene_id", "ref_gene_name"}:
        df = df.iloc[1:].reset_index(drop=True)
    return df[GFFCOMPARE_TMAP_COLUMNS]


def _gffcompare_structure_class(class_code: object) -> tuple[str, str, str]:
    code = _string_value(class_code, "-")
    if code == "=":
        return ("exact_intron_chain", "P4", "no_action")
    if code == "c":
        return ("contained_or_compatible", "P3", "inspect_terminal_exons")
    if code in {"j", "k"}:
        return ("splice_compatible_partial", "P2", "inspect_isoform_structure")
    if code in {"m", "n", "e", "o"}:
        return ("exonic_or_locus_overlap_only", "P2", "inspect_gene_structure")
    if code in {"s", "x", "i", "y", "p", "r"}:
        return ("ambiguous_or_problematic_overlap", "P2", "inspect_overlap_context")
    if code == "u":
        return ("query_intergenic_to_reference", "P3", "check_novel_or_false_positive")
    return ("no_reference_match", "P3", "check_reference_or_query_id")


def build_same_assembly_structure_audit(tmap: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "reference_gene_id",
        "reference_transcript_id",
        "query_gene_id",
        "query_transcript_id",
        "class_code",
        "structure_match_class",
        "review_priority",
        "suggested_action",
        "num_exons",
        "length",
        "reference_match_length",
        "review_reason",
    ]
    if tmap.empty:
        return _empty(columns)

    rows = []
    for _, row in tmap.iterrows():
        structure_class, priority, action = _gffcompare_structure_class(row.get("class_code"))
        reference_gene_id = _string_value(row.get("reference_gene_id"))
        reference_transcript_id = _string_value(row.get("reference_transcript_id"))
        query_gene_id = _string_value(row.get("query_gene_id"))
        query_transcript_id = _string_value(row.get("query_transcript_id"))
        rows.append(
            {
                "reference_gene_id": "" if reference_gene_id == "-" else reference_gene_id,
                "reference_transcript_id": "" if reference_transcript_id == "-" else reference_transcript_id,
                "query_gene_id": "" if query_gene_id == "-" else query_gene_id,
                "query_transcript_id": "" if query_transcript_id == "-" else query_transcript_id,
                "class_code": _string_value(row.get("class_code"), "-"),
                "structure_match_class": structure_class,
                "review_priority": priority,
                "suggested_action": action,
                "num_exons": _int_value(row.get("num_exons"), 0),
                "length": _int_value(row.get("length"), 0),
                "reference_match_length": _int_value(row.get("reference_match_length"), 0),
                "review_reason": (
                    f"gffcompare_class={_string_value(row.get('class_code'), '-')}; "
                    f"reference={reference_gene_id}|{reference_transcript_id}; "
                    f"query={query_gene_id}|{query_transcript_id}"
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["review_priority", "class_code", "reference_gene_id", "query_gene_id"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9)
        if s.name == "review_priority"
        else (s.map(GFFCOMPARE_CLASS_RANK).fillna(9) if s.name == "class_code" else s),
    )


def build_same_assembly_structure_summary(same_assembly_structure: pd.DataFrame) -> pd.DataFrame:
    columns = ["class_code", "structure_match_class", "review_priority", "count"]
    if same_assembly_structure.empty:
        return _empty(columns)
    return (
        same_assembly_structure.groupby(["class_code", "structure_match_class", "review_priority"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(
            ["review_priority", "count"],
            ascending=[True, False],
            key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
        )
    )


EXPECTED_PROTEIN_COLUMNS = [
    "expected_gene_id",
    "query_protein_id",
    "expected_source",
    "reference_stable_id",
    "symbol",
    "biotype",
    "orthogroup_id",
    "busco_id",
    "confidence",
]

REFERENCE_PROTEIN_HIT_COLUMNS = [
    "expected_gene_id",
    "query_protein_id",
    "target_gene_id",
    "target_stable_id",
    "target_transcript_id",
    "seq_region_name",
    "seq_region_start",
    "seq_region_end",
    "seq_region_strand",
    "aligner",
    "hit_rank",
    "percent_identity",
    "query_coverage",
    "target_coverage",
    "alignment_score",
    "frameshift_count",
    "stop_codon_count",
]


def _normalise_expected_proteins(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in EXPECTED_PROTEIN_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    if "expected_gene_id" not in df.columns and "query_protein_id" in df.columns:
        out["expected_gene_id"] = out["query_protein_id"]
    out["expected_gene_id"] = out.apply(
        lambda r: _string_value(r.get("expected_gene_id")) or _string_value(r.get("query_protein_id")),
        axis=1,
    )
    out["confidence"] = out["confidence"].map(lambda x: _string_value(x, "unknown").lower())
    return out[EXPECTED_PROTEIN_COLUMNS]


def _normalise_reference_protein_hits(df: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "target_seq_region_name": "seq_region_name",
        "target_start": "seq_region_start",
        "target_end": "seq_region_end",
        "target_strand": "seq_region_strand",
        "identity": "percent_identity",
        "qcov": "query_coverage",
        "tcov": "target_coverage",
    }
    out = df.rename(columns={k: v for k, v in rename.items() if k in df.columns}).copy()
    for col in REFERENCE_PROTEIN_HIT_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    for col in ("percent_identity", "query_coverage", "target_coverage", "alignment_score"):
        out[col] = out[col].map(lambda x: _numeric_value(x, 0.0))
    for col in ("hit_rank", "seq_region_start", "seq_region_end", "frameshift_count", "stop_codon_count"):
        out[col] = out[col].map(lambda x: _int_value(x, 0))
    out["seq_region_strand"] = out["seq_region_strand"].map(_normalise_strand)
    out["expected_gene_id"] = out.apply(
        lambda r: _string_value(r.get("expected_gene_id")) or _string_value(r.get("query_protein_id")),
        axis=1,
    )
    return out[REFERENCE_PROTEIN_HIT_COLUMNS]


def _best_protein_hit(group: pd.DataFrame) -> Optional[pd.Series]:
    if group.empty:
        return None
    work = group.copy()
    work["_rank_sort"] = work["hit_rank"].map(lambda x: _int_value(x, 999999) or 999999)
    work = work.sort_values(
        ["_rank_sort", "query_coverage", "percent_identity", "alignment_score"],
        ascending=[True, False, False, False],
    )
    return work.iloc[0]


def _protein_hit_query(best_hit: pd.Series) -> pd.Series:
    return pd.Series(
        {
            "seq_region_name": best_hit.get("seq_region_name"),
            "seq_region_start": best_hit.get("seq_region_start"),
            "seq_region_end": best_hit.get("seq_region_end"),
            "seq_region_strand": best_hit.get("seq_region_strand"),
        }
    )


def _protein_hit_class(
    best_hit: Optional[pd.Series],
    best_core: Optional[pd.Series],
    min_identity: float,
    min_query_coverage: float,
) -> tuple[str, str, str]:
    if best_hit is None:
        return ("no_protein_hit", "P2", "check_protein_mapping_or_true_absence")
    identity = _numeric_value(best_hit.get("percent_identity"), 0.0)
    query_coverage = _numeric_value(best_hit.get("query_coverage"), 0.0)
    frameshift_count = _int_value(best_hit.get("frameshift_count"), 0)
    stop_codon_count = _int_value(best_hit.get("stop_codon_count"), 0)
    if query_coverage < min_query_coverage or identity < min_identity:
        return ("weak_protein_hit", "P2", "inspect_low_confidence_homology")
    if frameshift_count > 0 or stop_codon_count > 0:
        return ("protein_hit_degraded", "P2", "inspect_frameshift_or_stop")
    if best_core is None:
        return ("protein_supported_no_core_gene", "P1", "try_candidate_rescue_from_protein")
    return ("protein_supported_built", "P4", "no_action")


def audit_reference_protein_set(
    expected_proteins: pd.DataFrame,
    protein_hits: pd.DataFrame,
    core_genes: pd.DataFrame,
    min_identity: float,
    min_query_coverage: float,
) -> pd.DataFrame:
    columns = [
        "expected_gene_id",
        "query_protein_id",
        "expected_source",
        "reference_stable_id",
        "symbol",
        "biotype",
        "orthogroup_id",
        "busco_id",
        "confidence",
        "best_target_gene_id",
        "best_target_stable_id",
        "best_target_transcript_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "aligner",
        "percent_identity",
        "query_coverage",
        "target_coverage",
        "frameshift_count",
        "stop_codon_count",
        "best_core_stable_id",
        "core_overlap_bp",
        "protein_hit_class",
        "review_priority",
        "suggested_action",
        "review_reason",
    ]
    hits = _normalise_reference_protein_hits(protein_hits)
    expected = _normalise_expected_proteins(expected_proteins)
    if expected.empty and not hits.empty:
        expected = hits[["expected_gene_id", "query_protein_id"]].drop_duplicates().copy()
        expected["expected_source"] = "reference_protein_set"
        expected["reference_stable_id"] = expected["expected_gene_id"]
        expected["symbol"] = ""
        expected["biotype"] = "protein_coding"
        expected["orthogroup_id"] = ""
        expected["busco_id"] = ""
        expected["confidence"] = "unknown"
        expected = expected[EXPECTED_PROTEIN_COLUMNS]
    if expected.empty:
        return _empty(columns)

    core = _normalise_gene_table(core_genes, "core")
    hits_by_expected = {str(k): v for k, v in hits.groupby("expected_gene_id", dropna=False)}
    rows = []
    for _, expected_row in expected.iterrows():
        expected_id = _string_value(expected_row.get("expected_gene_id"))
        query_protein_id = _string_value(expected_row.get("query_protein_id")) or expected_id
        group = hits_by_expected.get(expected_id, pd.DataFrame())
        if group.empty and query_protein_id:
            group = hits[hits["query_protein_id"].astype(str) == query_protein_id]
        best_hit = _best_protein_hit(group)
        best_core = None
        core_overlap = 0
        if best_hit is not None:
            target_stable_id = _string_value(best_hit.get("target_stable_id"))
            target_gene_id = _string_value(best_hit.get("target_gene_id"))
            if target_stable_id and "stable_id" in core.columns:
                matches = core[core["stable_id"].astype(str) == target_stable_id]
                if not matches.empty:
                    best_core = matches.iloc[0]
            if best_core is None and target_gene_id and "gene_id" in core.columns:
                matches = core[core["gene_id"].astype(str) == target_gene_id]
                if not matches.empty:
                    best_core = matches.iloc[0]
            if best_core is None:
                best_core, core_overlap = _best_overlap(_protein_hit_query(best_hit), core, same_strand=True)
            elif _string_value(best_hit.get("seq_region_name")):
                core_overlap = _overlap_bp(
                    best_hit.get("seq_region_start"),
                    best_hit.get("seq_region_end"),
                    best_core.get("seq_region_start"),
                    best_core.get("seq_region_end"),
                )
        protein_class, priority, action = _protein_hit_class(
            best_hit,
            best_core,
            min_identity=min_identity,
            min_query_coverage=min_query_coverage,
        )
        rows.append(
            {
                "expected_gene_id": expected_id,
                "query_protein_id": query_protein_id,
                "expected_source": _string_value(expected_row.get("expected_source"), "reference_protein_set"),
                "reference_stable_id": _string_value(expected_row.get("reference_stable_id")),
                "symbol": _string_value(expected_row.get("symbol")),
                "biotype": _string_value(expected_row.get("biotype"), "protein_coding"),
                "orthogroup_id": _string_value(expected_row.get("orthogroup_id")),
                "busco_id": _string_value(expected_row.get("busco_id")),
                "confidence": _string_value(expected_row.get("confidence"), "unknown").lower(),
                "best_target_gene_id": _string_value(best_hit.get("target_gene_id")) if best_hit is not None else "",
                "best_target_stable_id": (
                    _string_value(best_hit.get("target_stable_id")) if best_hit is not None else ""
                ),
                "best_target_transcript_id": (
                    _string_value(best_hit.get("target_transcript_id")) if best_hit is not None else ""
                ),
                "seq_region_name": _string_value(best_hit.get("seq_region_name")) if best_hit is not None else "",
                "seq_region_start": _int_value(best_hit.get("seq_region_start")) if best_hit is not None else 0,
                "seq_region_end": _int_value(best_hit.get("seq_region_end")) if best_hit is not None else 0,
                "seq_region_strand": _normalise_strand(best_hit.get("seq_region_strand"))
                if best_hit is not None
                else 0,
                "aligner": _string_value(best_hit.get("aligner")) if best_hit is not None else "",
                "percent_identity": _numeric_value(best_hit.get("percent_identity")) if best_hit is not None else 0.0,
                "query_coverage": _numeric_value(best_hit.get("query_coverage")) if best_hit is not None else 0.0,
                "target_coverage": _numeric_value(best_hit.get("target_coverage")) if best_hit is not None else 0.0,
                "frameshift_count": _int_value(best_hit.get("frameshift_count")) if best_hit is not None else 0,
                "stop_codon_count": _int_value(best_hit.get("stop_codon_count")) if best_hit is not None else 0,
                "best_core_stable_id": _string_value(best_core.get("stable_id")) if best_core is not None else "",
                "core_overlap_bp": core_overlap,
                "protein_hit_class": protein_class,
                "review_priority": priority,
                "suggested_action": action,
                "review_reason": (
                    f"{protein_class}; qcov="
                    f"{(_numeric_value(best_hit.get('query_coverage')) if best_hit is not None else 0.0):.2f}; "
                    f"pid={(_numeric_value(best_hit.get('percent_identity')) if best_hit is not None else 0.0):.2f}"
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["review_priority", "protein_hit_class", "expected_gene_id"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
    )


def build_reference_protein_summary(reference_protein_audit: pd.DataFrame) -> pd.DataFrame:
    columns = ["protein_hit_class", "review_priority", "count"]
    if reference_protein_audit.empty:
        return _empty(columns)
    return (
        reference_protein_audit.groupby(["protein_hit_class", "review_priority"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(
            ["review_priority", "count"],
            ascending=[True, False],
            key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
        )
    )


def build_feature_profile(
    evidence_fate: pd.DataFrame,
    expected_presence: pd.DataFrame,
    review_loci: pd.DataFrame,
    copy_number_audit: pd.DataFrame,
    busco_crosswalk: pd.DataFrame,
    completeness_profile: Optional[pd.DataFrame] = None,
    non_busco_high_confidence_losses: Optional[pd.DataFrame] = None,
    reference_protein_audit: Optional[pd.DataFrame] = None,
    same_assembly_structure: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    rows: list[dict] = []

    def add(group: str, name: str, value: float, denominator: Optional[float], description: str) -> None:
        fraction = value / denominator if denominator else pd.NA
        rows.append(
            {
                "metric_group": group,
                "metric_name": name,
                "value": value,
                "denominator": denominator if denominator is not None else pd.NA,
                "fraction": fraction,
                "description": description,
            }
        )

    n_layer = float(len(evidence_fate))
    add("evidence_fate", "layer_model_count", n_layer, None, "Layer models audited")
    if n_layer:
        add(
            "evidence_fate",
            "orphan_layer_model_count",
            float((evidence_fate["fate_class"] == "orphan_evidence").sum()),
            n_layer,
            "Layer models with no same-strand core representation",
        )
        add(
            "evidence_fate",
            "coding_orphan_layer_model_count",
            float(
                (
                    (evidence_fate["fate_class"] == "orphan_evidence")
                    & (evidence_fate["layer_coding_transcript_count"] > 0)
                ).sum()
            ),
            n_layer,
            "Coding layer models with no same-strand core representation",
        )
        add(
            "evidence_fate",
            "p1_evidence_issue_count",
            float((evidence_fate["review_priority"] == "P1").sum()),
            n_layer,
            "Layer evidence rows requiring high-priority review",
        )

    n_expected = float(len(expected_presence))
    add("expected_presence", "expected_gene_count", n_expected, None, "Expected genes audited")
    if n_expected:
        for presence_class in (
            "present_clean",
            "missing_with_evidence",
            "projection_only",
            "assembly_limited",
            "present_degraded",
        ):
            add(
                "expected_presence",
                f"{presence_class}_count",
                float((expected_presence["presence_class"] == presence_class).sum()),
                n_expected,
                f"Expected genes classified as {presence_class}",
            )
        add(
            "expected_presence",
            "p1_expected_issue_count",
            float((expected_presence["review_priority"] == "P1").sum()),
            n_expected,
            "Expected-gene rows requiring high-priority review",
        )

    n_review = float(len(review_loci))
    add("review", "review_locus_count", n_review, None, "P1/P2 review loci emitted")
    if n_review:
        add(
            "review",
            "p1_review_locus_count",
            float((review_loci["review_priority"] == "P1").sum()),
            n_review,
            "Review loci with P1 priority",
        )

    n_copy = float(len(copy_number_audit))
    add("copy_number", "copy_group_count", n_copy, None, "Copy-number groups audited")
    if n_copy:
        add(
            "copy_number",
            "copy_number_issue_count",
            float((copy_number_audit["copy_number_class"] != "copy_number_as_expected").sum()),
            n_copy,
            "Copy-number groups not matching expectation",
        )

    n_busco = float(len(busco_crosswalk))
    add("busco", "busco_expected_gene_count", n_busco, None, "Expected genes with BUSCO IDs")
    if n_busco:
        add(
            "busco",
            "busco_p1_or_p2_issue_count",
            float(busco_crosswalk["review_priority"].isin(["P1", "P2"]).sum()),
            n_busco,
            "BUSCO-linked expected genes requiring P1/P2 review",
        )

    if completeness_profile is not None and not completeness_profile.empty:
        for panel_id in ("high_confidence", "high_confidence_non_busco", "busco_linked"):
            panel = completeness_profile[completeness_profile["panel_id"] == panel_id]
            if panel.empty:
                continue
            panel_row = panel.iloc[0]
            n_panel = float(panel_row["n_expected"])
            add(
                "completeness",
                f"{panel_id}_present_clean_count",
                float(panel_row["n_present_clean"]),
                n_panel,
                f"Cleanly represented genes in completeness panel {panel_id}",
            )
            add(
                "completeness",
                f"{panel_id}_actionable_loss_count",
                float(panel_row["n_actionable_loss"]),
                n_panel,
                f"Actionable losses in completeness panel {panel_id}",
            )
    if non_busco_high_confidence_losses is not None:
        add(
            "completeness",
            "non_busco_high_confidence_loss_count",
            float(len(non_busco_high_confidence_losses)),
            None,
            "High-confidence non-BUSCO expected genes with actionable loss/degradation",
        )

    if reference_protein_audit is not None:
        n_reference_proteins = float(len(reference_protein_audit))
        add(
            "reference_protein",
            "expected_protein_count",
            n_reference_proteins,
            None,
            "Reference proteins audited with BUSCO-like protein-set logic",
        )
        if n_reference_proteins:
            for klass in ("protein_supported_built", "protein_supported_no_core_gene", "protein_hit_degraded"):
                add(
                    "reference_protein",
                    f"{klass}_count",
                    float((reference_protein_audit["protein_hit_class"] == klass).sum()),
                    n_reference_proteins,
                    f"Reference proteins classified as {klass}",
                )

    if same_assembly_structure is not None:
        n_same_assembly = float(len(same_assembly_structure))
        add(
            "same_assembly_structure",
            "gffcompare_transcript_count",
            n_same_assembly,
            None,
            "GffCompare query transcript rows audited",
        )
        if n_same_assembly:
            add(
                "same_assembly_structure",
                "exact_intron_chain_count",
                float((same_assembly_structure["structure_match_class"] == "exact_intron_chain").sum()),
                n_same_assembly,
                "Same-assembly query transcripts with exact intron-chain matches",
            )
            add(
                "same_assembly_structure",
                "p1_or_p2_structure_issue_count",
                float(same_assembly_structure["review_priority"].isin(["P1", "P2"]).sum()),
                n_same_assembly,
                "Same-assembly transcript rows requiring structure review",
            )

    return pd.DataFrame(
        rows,
        columns=["metric_group", "metric_name", "value", "denominator", "fraction", "description"],
    )


def build_release_readiness(
    expected_presence: pd.DataFrame,
    copy_number_audit: pd.DataFrame,
    reference_protein_audit: pd.DataFrame,
    same_assembly_structure: pd.DataFrame,
    completeness_profile: pd.DataFrame,
    busco_complete_percent: Optional[float],
    busco_floor_percent: float,
    max_high_confidence_actionable_loss_fraction: float,
    max_copy_number_issue_fraction: float,
    require_expected_genes: bool,
) -> pd.DataFrame:
    columns = [
        "gate_id",
        "status",
        "priority",
        "observed_value",
        "threshold",
        "rationale",
        "required_action",
        "target_table",
        "target_filter",
    ]
    rows: list[dict] = []

    def add(
        gate_id: str,
        status: str,
        priority: str,
        observed_value: object,
        threshold: object,
        rationale: str,
        required_action: str,
        target_table: str,
        target_filter: str,
    ) -> None:
        rows.append(
            {
                "gate_id": gate_id,
                "status": status,
                "priority": priority,
                "observed_value": observed_value,
                "threshold": threshold,
                "rationale": rationale,
                "required_action": required_action,
                "target_table": target_table,
                "target_filter": target_filter,
            }
        )

    if busco_complete_percent is not None:
        status = "FAIL" if busco_complete_percent < busco_floor_percent else "PASS"
        add(
            "busco_floor",
            status,
            "P1" if status == "FAIL" else "P4",
            f"{busco_complete_percent:.2f}",
            f">={busco_floor_percent:.2f}",
            "BUSCO completeness is below the configured floor.",
            "Do not release on BUSCO alone; run expected-gene and protein-set recovery before release.",
            "external_busco",
            "complete_percent",
        )

    if require_expected_genes and expected_presence.empty:
        add(
            "expected_catalogue_required",
            "FAIL",
            "P1",
            0,
            ">0",
            "No expected-gene catalogue was audited, so reference-informed completeness cannot be claimed.",
            "Provide prior/same-species/reference expected genes and projections before release sign-off.",
            "expected_gene_presence.tsv",
            "missing table",
        )

    if not completeness_profile.empty:
        high = completeness_profile[completeness_profile["panel_id"] == "high_confidence"]
        if not high.empty:
            fraction = _numeric_value(high.iloc[0].get("actionable_loss_fraction"), 0.0)
            status = "FAIL" if fraction > max_high_confidence_actionable_loss_fraction else "PASS"
            add(
                "high_confidence_actionable_loss",
                status,
                "P1" if status == "FAIL" else "P4",
                f"{fraction:.4f}",
                f"<={max_high_confidence_actionable_loss_fraction:.4f}",
                "High-confidence expected genes are missing, degraded, split, fused, or wrong-biotype.",
                "Resolve or explicitly classify high-confidence losses as assembly-limited before release.",
                "expected_gene_presence.tsv",
                "confidence=high;presence_class in actionable classes",
            )

    if not expected_presence.empty:
        missing_with_evidence = int((expected_presence["presence_class"] == "missing_with_evidence").sum())
        if missing_with_evidence:
            add(
                "expected_genes_missing_with_evidence",
                "FAIL",
                "P1",
                missing_with_evidence,
                0,
                "Expected genes have local layer support but no core gene.",
                "Run candidate rescue or adjust layer/source selection before release.",
                "expected_gene_presence.tsv",
                "presence_class=missing_with_evidence",
            )

    if not copy_number_audit.empty:
        copy_issues = copy_number_audit[copy_number_audit["copy_number_class"] != "copy_number_as_expected"]
        issue_fraction = len(copy_issues) / len(copy_number_audit) if len(copy_number_audit) else 0.0
        status = "FAIL" if issue_fraction > max_copy_number_issue_fraction else "PASS"
        add(
            "copy_number_regression",
            status,
            "P1" if status == "FAIL" else "P4",
            f"{issue_fraction:.4f}",
            f"<={max_copy_number_issue_fraction:.4f}",
            "Expected copy-number groups do not match observed core copy number.",
            "Review collapsed, missing, and expanded families; do not release unresolved key copy losses.",
            "copy_number_audit.tsv",
            "copy_number_class!=copy_number_as_expected",
        )

    if not reference_protein_audit.empty:
        protein_missing = int(
            (reference_protein_audit["protein_hit_class"] == "protein_supported_no_core_gene").sum()
        )
        if protein_missing:
            add(
                "protein_supported_missing_core",
                "FAIL",
                "P1",
                protein_missing,
                0,
                "Reference proteins align confidently to the assembly but lack overlapping core genes.",
                "Use these protein-supported loci as rescue candidates or evidence-ingest failures.",
                "reference_protein_audit.tsv",
                "protein_hit_class=protein_supported_no_core_gene",
            )

    if not same_assembly_structure.empty:
        p2_structure = int(same_assembly_structure["review_priority"].isin(["P1", "P2"]).sum())
        if p2_structure:
            add(
                "same_assembly_structure_regression",
                "WARN",
                "P2",
                p2_structure,
                0,
                "Same-assembly transcript structures differ from the reference annotation.",
                "Review structural regressions before claiming annotation consistency.",
                "same_assembly_structure.tsv",
                "review_priority in P1,P2",
            )

    if not rows:
        add(
            "release_readiness",
            "PASS",
            "P4",
            "no_blocking_signals",
            "configured gates",
            "No configured release gate failed.",
            "Proceed with normal review, preserving audit outputs for release notes.",
            "feature_profile.tsv",
            "all rows",
        )
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["status", "priority", "gate_id"],
        key=lambda s: s.map({"FAIL": 0, "WARN": 1, "PASS": 2}).fillna(3)
        if s.name == "status"
        else (s.map(PRIORITY_RANK).fillna(9) if s.name == "priority" else s),
    )


def build_recommendations(
    evidence_fate: pd.DataFrame,
    expected_presence: pd.DataFrame,
    copy_number_audit: pd.DataFrame,
    busco_proxy_calibration: pd.DataFrame,
    non_busco_high_confidence_losses: pd.DataFrame,
    reference_protein_audit: Optional[pd.DataFrame] = None,
    same_assembly_structure: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    columns = [
        "recommendation_id",
        "priority",
        "scope",
        "recommendation_type",
        "trigger_count",
        "rationale",
        "next_action",
        "target_table",
        "target_filter",
        "confidence",
    ]
    rows: list[dict] = []

    def add(
        recommendation_id: str,
        priority: str,
        scope: str,
        recommendation_type: str,
        trigger_count: int,
        rationale: str,
        next_action: str,
        target_table: str,
        target_filter: str,
        confidence: str = "medium",
    ) -> None:
        if trigger_count <= 0:
            return
        rows.append(
            {
                "recommendation_id": recommendation_id,
                "priority": priority,
                "scope": scope,
                "recommendation_type": recommendation_type,
                "trigger_count": int(trigger_count),
                "rationale": rationale,
                "next_action": next_action,
                "target_table": target_table,
                "target_filter": target_filter,
                "confidence": confidence,
            }
        )

    if not evidence_fate.empty:
        coding_orphans = evidence_fate[
            (evidence_fate["failure_class"] == "no_core_gene_built")
            & (evidence_fate["layer_coding_transcript_count"] > 0)
        ]
        add(
            "candidate_rescue_from_layer",
            "P1",
            "locus",
            "candidate_rescue",
            len(coding_orphans),
            "Coding layer models have no same-strand core representation.",
            "Inspect P1 loci; test whether layer candidates should be rescued or source/layer ordering changed.",
            "evidence_fate.tsv",
            "failure_class=no_core_gene_built;layer_coding_transcript_count>0",
            "high",
        )
        wrong_biotype = evidence_fate[evidence_fate["failure_class"] == "coding_evidence_not_protein_coding"]
        add(
            "coding_evidence_biotype_review",
            "P1",
            "locus",
            "biotype_or_cds_review",
            len(wrong_biotype),
            "Coding layer evidence is represented by a non-protein-coding core model.",
            "Review CDS, biotype assignment, and pseudogene/non-coding classification at these loci.",
            "evidence_fate.tsv",
            "failure_class=coding_evidence_not_protein_coding",
            "high",
        )
        partial = evidence_fate[
            evidence_fate["failure_class"].isin(
                ["core_span_underrepresents_evidence", "partial_core_representation"]
            )
        ]
        add(
            "underrepresented_layer_models",
            "P2",
            "locus",
            "structure_review",
            len(partial),
            "Core genes only partially span layer evidence.",
            "Inspect for fragmented models, split genes, truncation, or overly aggressive competition/collapse.",
            "evidence_fate.tsv",
            "failure_class in core_span_underrepresents_evidence,partial_core_representation",
        )

    if not expected_presence.empty:
        missing_with_evidence = expected_presence[
            expected_presence["presence_class"] == "missing_with_evidence"
        ]
        add(
            "expected_gene_rescue",
            "P1",
            "expected_gene",
            "expected_gene_rescue",
            len(missing_with_evidence),
            "Expected genes have local evidence but no core model.",
            "Prioritize these as build-limited losses; inspect source evidence and candidate generation.",
            "expected_gene_presence.tsv",
            "presence_class=missing_with_evidence",
            "high",
        )
        projection_only = expected_presence[expected_presence["presence_class"] == "projection_only"]
        add(
            "projection_without_local_support",
            "P2",
            "expected_gene",
            "evidence_ingest_or_projection_review",
            len(projection_only),
            "Expected genes project to the assembly but lack core/layer support.",
            "Check projection quality, evidence ingest, and whether source data should have entered the layer DB.",
            "expected_gene_presence.tsv",
            "presence_class=projection_only",
        )
        assembly_limited = expected_presence[expected_presence["presence_class"] == "assembly_limited"]
        add(
            "assembly_limited_expected_genes",
            "P2",
            "assembly",
            "assembly_quality_review",
            len(assembly_limited),
            "Expected genes overlap assembly gap signal and lack core representation.",
            "Separate these from genebuild failures; report assembly-limited completeness separately.",
            "expected_gene_presence.tsv",
            "presence_class=assembly_limited",
            "high",
        )

    if not copy_number_audit.empty:
        copy_issues = copy_number_audit[
            copy_number_audit["copy_number_class"] != "copy_number_as_expected"
        ]
        add(
            "copy_number_review",
            "P2",
            "family",
            "copy_number_review",
            len(copy_issues),
            "Observed core copy number does not match expected copy number.",
            "Review paralogue families for collapsed copies, expansions, and projection ambiguity.",
            "copy_number_audit.tsv",
            "copy_number_class!=copy_number_as_expected",
        )

    add(
        "non_busco_completeness_review",
        "P1",
        "completeness_panel",
        "beyond_busco_review",
        len(non_busco_high_confidence_losses),
        "High-confidence non-BUSCO expected genes show actionable loss or degradation.",
        "Do not rely on BUSCO alone; inspect non-BUSCO losses and source-specific patterns.",
        "non_busco_high_confidence_losses.tsv",
        "all rows",
        "high" if len(non_busco_high_confidence_losses) else "medium",
    )

    if not busco_proxy_calibration.empty:
        delta_rows = busco_proxy_calibration[
            busco_proxy_calibration["comparison"] == "non_busco_minus_busco_actionable_loss"
        ]
        if not delta_rows.empty:
            delta = _numeric_value(delta_rows.iloc[0].get("actionable_loss_fraction"), 0.0)
            if delta > 0:
                add(
                    "busco_false_reassurance_risk",
                    "P1",
                    "completeness_panel",
                    "busco_proxy_warning",
                    1,
                    "High-confidence non-BUSCO genes fail more often than BUSCO-linked genes.",
                    "Use completeness_profile.tsv as the primary completeness claim; treat BUSCO as a sentinel panel.",
                    "busco_proxy_calibration.tsv",
                    "comparison=non_busco_minus_busco_actionable_loss",
                    "medium",
                )

    if reference_protein_audit is not None and not reference_protein_audit.empty:
        protein_rescue = reference_protein_audit[
            reference_protein_audit["protein_hit_class"] == "protein_supported_no_core_gene"
        ]
        add(
            "reference_protein_rescue",
            "P1",
            "reference_protein",
            "protein_supported_candidate_rescue",
            len(protein_rescue),
            "Reference proteins align confidently to the assembly but lack an overlapping core gene.",
            "Inspect protein-supported loci; use them as candidate-rescue evidence or evidence-ingest checks.",
            "reference_protein_audit.tsv",
            "protein_hit_class=protein_supported_no_core_gene",
            "high",
        )
        degraded = reference_protein_audit[
            reference_protein_audit["protein_hit_class"] == "protein_hit_degraded"
        ]
        add(
            "reference_protein_degradation_review",
            "P2",
            "reference_protein",
            "protein_integrity_review",
            len(degraded),
            "Reference protein alignments contain frameshift or stop-codon signals.",
            "Inspect whether the target assembly, projection, or final model truncates the coding sequence.",
            "reference_protein_audit.tsv",
            "protein_hit_class=protein_hit_degraded",
        )

    if same_assembly_structure is not None and not same_assembly_structure.empty:
        structure_issues = same_assembly_structure[same_assembly_structure["review_priority"].isin(["P1", "P2"])]
        add(
            "same_assembly_structure_review",
            "P2",
            "same_assembly_annotation",
            "transcript_structure_review",
            len(structure_issues),
            "Same-assembly query transcripts differ structurally from the reference annotation.",
            "Use GffCompare class codes to separate exact intron-chain matches from partial, overlap-only, or novel rows.",
            "same_assembly_structure.tsv",
            "review_priority in P1,P2",
        )

    if not rows:
        return _empty(columns)
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["priority", "trigger_count"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "priority" else -s,
    )


def _review_rows_from_evidence(df: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for _, row in df.iterrows():
        priority = _string_value(row.get("review_priority"), "P4")
        if PRIORITY_RANK.get(priority, 9) > 2:
            continue
        rows.append(
            {
                "audit_track": "evidence_fate",
                "review_priority": priority,
                "seq_region_name": row["seq_region_name"],
                "seq_region_start": int(row["seq_region_start"]),
                "seq_region_end": int(row["seq_region_end"]),
                "seq_region_strand": int(row["seq_region_strand"]),
                "subject_id": row["layer_stable_id"] or row["layer_gene_id"],
                "class": row["failure_class"],
                "suggested_action": row["suggested_action"],
                "review_reason": row["review_reason"],
            }
        )
    return rows


def _review_rows_from_expected(df: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for _, row in df.iterrows():
        priority = _string_value(row.get("review_priority"), "P4")
        if PRIORITY_RANK.get(priority, 9) > 2:
            continue
        rows.append(
            {
                "audit_track": "expected_presence",
                "review_priority": priority,
                "seq_region_name": row["seq_region_name"],
                "seq_region_start": int(row["seq_region_start"]),
                "seq_region_end": int(row["seq_region_end"]),
                "seq_region_strand": int(row["seq_region_strand"]),
                "subject_id": row["expected_gene_id"],
                "class": row["presence_class"],
                "suggested_action": row["suggested_action"],
                "review_reason": row["review_reason"],
            }
        )
    return rows


def build_review_loci(evidence_fate: pd.DataFrame, expected_presence: pd.DataFrame, top_n: int) -> pd.DataFrame:
    rows = _review_rows_from_evidence(evidence_fate) + _review_rows_from_expected(expected_presence)
    if not rows:
        return _empty(
            [
                "audit_track",
                "review_priority",
                "seq_region_name",
                "seq_region_start",
                "seq_region_end",
                "seq_region_strand",
                "subject_id",
                "class",
                "suggested_action",
                "review_reason",
            ]
        )
    out = pd.DataFrame(rows)
    out = out[out["seq_region_name"].astype(str) != ""]
    out = out[out["seq_region_end"].astype(int) >= out["seq_region_start"].astype(int)]
    out = out.sort_values(
        ["review_priority", "audit_track", "seq_region_name", "seq_region_start"],
        key=lambda s: s.map(PRIORITY_RANK).fillna(9) if s.name == "review_priority" else s,
    )
    return out.head(top_n).reset_index(drop=True)


def write_bed(review_loci: pd.DataFrame, path: str, pad_bp: int = 0) -> None:
    ensure_dir(os.path.dirname(path))
    if review_loci.empty:
        open(path, "w", encoding="utf-8").close()
        return
    rows = []
    for _, row in review_loci.iterrows():
        start0 = max(0, int(row["seq_region_start"]) - 1 - pad_bp)
        end0 = int(row["seq_region_end"]) + pad_bp
        name = "|".join(
            [
                _string_value(row.get("review_priority"), "P4"),
                _string_value(row.get("audit_track"), "audit"),
                _string_value(row.get("class"), "class"),
                _string_value(row.get("subject_id"), "subject"),
            ]
        )
        score = max(0, 1000 - PRIORITY_RANK.get(_string_value(row.get("review_priority"), "P4"), 4) * 200)
        strand = {1: "+", -1: "-"}.get(_normalise_strand(row.get("seq_region_strand")), ".")
        rows.append([row["seq_region_name"], start0, end0, name, score, strand])
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)


def write_review_bed_sets(review_loci: pd.DataFrame, output_dir: str, pad_bp: int = 0) -> list[str]:
    bed_dir = os.path.join(output_dir, "review_beds")
    ensure_dir(bed_dir)
    written = []
    class_names = [
        "no_core_gene_built",
        "missing_with_evidence",
        "projection_only",
        "assembly_limited",
        "present_degraded",
        "present_wrong_biotype",
        "split",
        "fused",
    ]
    for class_name in class_names:
        subset = review_loci[review_loci["class"] == class_name] if not review_loci.empty else review_loci
        path = os.path.join(bed_dir, f"{class_name}.bed")
        write_bed(subset, path, pad_bp=pad_bp)
        written.append(path)
    return written


def write_expected_presence_bed(expected_rows: pd.DataFrame, path: str, pad_bp: int = 0) -> None:
    ensure_dir(os.path.dirname(path))
    if expected_rows.empty:
        open(path, "w", encoding="utf-8").close()
        return
    rows = []
    for _, row in expected_rows.iterrows():
        seq = _string_value(row.get("seq_region_name"))
        start = _int_value(row.get("seq_region_start"))
        end = _int_value(row.get("seq_region_end"))
        if not seq or end < start:
            continue
        start0 = max(0, start - 1 - pad_bp)
        end0 = end + pad_bp
        name = "|".join(
            [
                _string_value(row.get("review_priority"), "P4"),
                _string_value(row.get("presence_class"), "expected"),
                _string_value(row.get("expected_gene_id"), "expected_gene"),
            ]
        )
        score = max(0, 1000 - PRIORITY_RANK.get(_string_value(row.get("review_priority"), "P4"), 4) * 200)
        strand = {1: "+", -1: "-"}.get(_normalise_strand(row.get("seq_region_strand")), ".")
        rows.append([seq, start0, end0, name, score, strand])
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)


def _top_counts(df: pd.DataFrame, col: str, limit: int = 8) -> list[tuple[str, int]]:
    if df.empty or col not in df.columns:
        return []
    return [(str(k), int(v)) for k, v in df[col].value_counts().head(limit).items()]


def write_actionable_summary(
    path: str,
    evidence_fate: pd.DataFrame,
    expected_presence: pd.DataFrame,
    review_loci: pd.DataFrame,
    feature_profile: Optional[pd.DataFrame] = None,
    recommendations: Optional[pd.DataFrame] = None,
    release_readiness: Optional[pd.DataFrame] = None,
) -> None:
    ensure_dir(os.path.dirname(path))
    lines = [
        "# Genebuild observability summary",
        "",
        "This report is generated from the implemented draft audit. It is intended",
        "to point reviewers at actionable loci, not to replace manual review.",
        "",
        "## High-priority review counts",
        "",
    ]
    if review_loci.empty:
        lines.append("No P0/P1/P2 review loci were produced.")
    else:
        counts = review_loci["review_priority"].value_counts().sort_index()
        for priority, count in counts.items():
            lines.append(f"- {priority}: {int(count)} loci")
    lines.extend(["", "## Evidence fate signals", ""])
    if evidence_fate.empty:
        lines.append("No layer evidence rows were available.")
    else:
        for key, count in _top_counts(evidence_fate, "failure_class"):
            lines.append(f"- {key}: {count}")
    lines.extend(["", "## Expected gene presence signals", ""])
    if expected_presence.empty:
        lines.append("No expected-gene catalogue was provided.")
    else:
        for key, count in _top_counts(expected_presence, "presence_class"):
            lines.append(f"- {key}: {count}")
    if feature_profile is not None and not feature_profile.empty:
        lines.extend(["", "## Feature profile highlights", ""])
        highlight_names = {
            "orphan_layer_model_count",
            "coding_orphan_layer_model_count",
            "missing_with_evidence_count",
            "projection_only_count",
            "assembly_limited_count",
            "copy_number_issue_count",
            "busco_p1_or_p2_issue_count",
        }
        highlights = feature_profile[feature_profile["metric_name"].isin(highlight_names)]
        for _, row in highlights.iterrows():
            fraction = row["fraction"]
            fraction_text = "" if pd.isna(fraction) else f" ({float(fraction):.1%})"
            lines.append(f"- {row['metric_name']}: {int(row['value'])}{fraction_text}")
    if recommendations is not None and not recommendations.empty:
        lines.extend(["", "## Recommended actions", ""])
        for _, row in recommendations.head(12).iterrows():
            lines.append(
                f"- {row['priority']} {row['recommendation_type']}: "
                f"{row['next_action']} ({int(row['trigger_count'])} triggers)"
            )
    if release_readiness is not None and not release_readiness.empty:
        lines.extend(["", "## Release readiness gates", ""])
        for _, row in release_readiness.head(12).iterrows():
            lines.append(
                f"- {row['status']} {row['gate_id']}: {row['required_action']} "
                f"(observed={row['observed_value']}, threshold={row['threshold']})"
            )
    lines.extend(["", "## Top review loci", ""])
    if review_loci.empty:
        lines.append("No review loci to list.")
    else:
        for _, row in review_loci.head(20).iterrows():
            lines.append(
                f"- {row['review_priority']} "
                f"{row['seq_region_name']}:{row['seq_region_start']}-{row['seq_region_end']} "
                f"{row['audit_track']} {row['class']} {row['subject_id']} "
                f"-> {row['suggested_action']}"
            )
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def run_audit(args: argparse.Namespace) -> int:
    paths = _paths_from_args(args)
    ensure_dir(paths.observability_dir)

    core_genes = _load_named_table(paths.extract_dir, "core_genes", paths.fmt)
    layer_genes = _load_named_table(paths.extract_dir, "layer_genes", paths.fmt)
    core_tx = _load_named_table(paths.extract_dir, "core_transcripts", paths.fmt, required=False)
    layer_tx = _load_named_table(paths.extract_dir, "layer_transcripts", paths.fmt, required=False)
    core_tr = _load_named_table(paths.extract_dir, "core_translations", paths.fmt, required=False)
    layer_tr = _load_named_table(paths.extract_dir, "layer_translations", paths.fmt, required=False)

    loci, gene_to_locus = build_or_load_loci(paths, core_genes, layer_genes, args.locus_gap_bp)

    evidence_fate = audit_evidence_fate(core_genes, layer_genes, core_tx, layer_tx, core_tr, layer_tr)

    target_seq_regions = sorted(
        set(core_genes.get("seq_region_name", pd.Series(dtype=str)).dropna().astype(str))
        | set(layer_genes.get("seq_region_name", pd.Series(dtype=str)).dropna().astype(str))
    )
    generated_expected = populate_expected_inputs(args, paths, target_seq_regions)

    expected_presence = pd.DataFrame()
    if args.expected_genes:
        expected_genes = _read_table(args.expected_genes)
        expected_projections = _read_table(args.expected_projections) if args.expected_projections else pd.DataFrame()
        expected_presence = audit_expected_presence(expected_genes, expected_projections, core_genes, layer_genes)

    same_assembly_structure = pd.DataFrame()
    same_assembly_structure_summary = pd.DataFrame()
    if getattr(args, "gffcompare_tmap", None):
        same_assembly_structure = build_same_assembly_structure_audit(read_gffcompare_tmap(args.gffcompare_tmap))
        same_assembly_structure_summary = build_same_assembly_structure_summary(same_assembly_structure)

    reference_protein_audit = pd.DataFrame()
    reference_protein_summary = pd.DataFrame()
    if getattr(args, "expected_proteins", None) or getattr(args, "reference_protein_hits", None):
        expected_proteins = (
            _read_table(args.expected_proteins) if getattr(args, "expected_proteins", None) else pd.DataFrame()
        )
        reference_protein_hits = (
            _read_table(args.reference_protein_hits)
            if getattr(args, "reference_protein_hits", None)
            else pd.DataFrame()
        )
        reference_protein_audit = audit_reference_protein_set(
            expected_proteins,
            reference_protein_hits,
            core_genes,
            min_identity=args.protein_min_identity,
            min_query_coverage=args.protein_min_query_coverage,
        )
        reference_protein_summary = build_reference_protein_summary(reference_protein_audit)

    failure_summary = _failure_summary(evidence_fate, expected_presence)
    review_loci = build_review_loci(evidence_fate, expected_presence, args.top_n)
    audit_loci = build_audit_loci(loci, expected_presence, review_loci)
    source_profile = build_source_profile(evidence_fate)
    expected_source_profile = build_expected_source_profile(expected_presence)
    busco_crosswalk = build_busco_crosswalk(expected_presence)
    completeness_profile = build_completeness_profile(expected_presence)
    non_busco_high_confidence_losses = build_non_busco_high_confidence_losses(expected_presence)
    busco_proxy_calibration = build_busco_proxy_calibration(expected_presence)
    copy_number_audit = build_copy_number_audit(expected_presence)
    biotype_transition = build_biotype_transition(evidence_fate, expected_presence)
    recommendations = build_recommendations(
        evidence_fate,
        expected_presence,
        copy_number_audit,
        busco_proxy_calibration,
        non_busco_high_confidence_losses,
        reference_protein_audit,
        same_assembly_structure,
    )
    feature_profile = build_feature_profile(
        evidence_fate,
        expected_presence,
        review_loci,
        copy_number_audit,
        busco_crosswalk,
        completeness_profile,
        non_busco_high_confidence_losses,
        reference_protein_audit,
        same_assembly_structure,
    )
    release_readiness = build_release_readiness(
        expected_presence,
        copy_number_audit,
        reference_protein_audit,
        same_assembly_structure,
        completeness_profile,
        getattr(args, "busco_complete_percent", None),
        args.busco_floor_percent,
        args.max_high_confidence_actionable_loss_fraction,
        args.max_copy_number_issue_fraction,
        args.require_expected_genes,
    )

    suffix = args.format
    _write_table(audit_loci, os.path.join(paths.observability_dir, f"audit_loci.{suffix}"), args.format)
    _write_table(
        gene_to_locus,
        os.path.join(paths.observability_dir, f"gene_to_locus.{suffix}"),
        args.format,
    )
    _write_table(evidence_fate, os.path.join(paths.observability_dir, f"evidence_fate.{suffix}"), args.format)
    if not expected_presence.empty or args.expected_genes:
        _write_table(
            expected_presence,
            os.path.join(paths.observability_dir, f"expected_gene_presence.{suffix}"),
            args.format,
        )
    _write_table(source_profile, os.path.join(paths.observability_dir, f"source_profile.{suffix}"), args.format)
    _write_table(
        expected_source_profile,
        os.path.join(paths.observability_dir, f"expected_source_profile.{suffix}"),
        args.format,
    )
    _write_table(
        busco_crosswalk,
        os.path.join(paths.observability_dir, f"busco_expected_crosswalk.{suffix}"),
        args.format,
    )
    _write_table(
        completeness_profile,
        os.path.join(paths.observability_dir, f"completeness_profile.{suffix}"),
        args.format,
    )
    _write_table(
        non_busco_high_confidence_losses,
        os.path.join(paths.observability_dir, f"non_busco_high_confidence_losses.{suffix}"),
        args.format,
    )
    _write_table(
        busco_proxy_calibration,
        os.path.join(paths.observability_dir, f"busco_proxy_calibration.{suffix}"),
        args.format,
    )
    _write_table(
        copy_number_audit,
        os.path.join(paths.observability_dir, f"copy_number_audit.{suffix}"),
        args.format,
    )
    _write_table(
        biotype_transition,
        os.path.join(paths.observability_dir, f"biotype_transition.{suffix}"),
        args.format,
    )
    _write_table(
        recommendations,
        os.path.join(paths.observability_dir, f"recommendations.{suffix}"),
        args.format,
    )
    if not same_assembly_structure.empty or getattr(args, "gffcompare_tmap", None):
        _write_table(
            same_assembly_structure,
            os.path.join(paths.observability_dir, f"same_assembly_structure.{suffix}"),
            args.format,
        )
        _write_table(
            same_assembly_structure_summary,
            os.path.join(paths.observability_dir, f"same_assembly_structure_summary.{suffix}"),
            args.format,
        )
    if not reference_protein_audit.empty or getattr(args, "expected_proteins", None) or getattr(
        args, "reference_protein_hits", None
    ):
        _write_table(
            reference_protein_audit,
            os.path.join(paths.observability_dir, f"reference_protein_audit.{suffix}"),
            args.format,
        )
        _write_table(
            reference_protein_summary,
            os.path.join(paths.observability_dir, f"reference_protein_summary.{suffix}"),
            args.format,
        )
    _write_table(
        release_readiness,
        os.path.join(paths.observability_dir, f"release_readiness.{suffix}"),
        args.format,
    )
    _write_table(feature_profile, os.path.join(paths.observability_dir, f"feature_profile.{suffix}"), args.format)
    _write_table(failure_summary, os.path.join(paths.observability_dir, f"failure_mode_summary.{suffix}"), args.format)
    _write_table(review_loci, os.path.join(paths.observability_dir, f"review_loci.{suffix}"), args.format)
    write_bed(review_loci, os.path.join(paths.observability_dir, "review_loci.bed"), pad_bp=args.bed_pad_bp)
    review_bed_paths = write_review_bed_sets(review_loci, paths.observability_dir, pad_bp=args.bed_pad_bp)
    completeness_bed_dir = os.path.join(paths.observability_dir, "completeness_beds")
    write_expected_presence_bed(
        non_busco_high_confidence_losses,
        os.path.join(completeness_bed_dir, "non_busco_high_confidence_losses.bed"),
        pad_bp=args.bed_pad_bp,
    )
    if expected_presence.empty:
        high_conf_missing = expected_presence
    else:
        high_conf_missing = expected_presence[
            expected_presence.apply(_is_high_confidence, axis=1)
            & expected_presence.apply(_is_missing_actionable, axis=1)
        ]
    write_expected_presence_bed(
        high_conf_missing,
        os.path.join(completeness_bed_dir, "high_confidence_actionable_losses.bed"),
        pad_bp=args.bed_pad_bp,
    )
    protein_review = (
        reference_protein_audit[reference_protein_audit["review_priority"].isin(["P1", "P2"])]
        if not reference_protein_audit.empty
        else reference_protein_audit
    )
    protein_bed = protein_review.copy()
    if not protein_bed.empty:
        protein_bed["presence_class"] = protein_bed["protein_hit_class"]
    write_expected_presence_bed(
        protein_bed,
        os.path.join(completeness_bed_dir, "reference_protein_issues.bed"),
        pad_bp=args.bed_pad_bp,
    )
    write_actionable_summary(
        os.path.join(paths.observability_dir, "actionable_summary.md"),
        evidence_fate,
        expected_presence,
        review_loci,
        feature_profile,
        recommendations,
        release_readiness,
    )
    write_manifest(
        {
            "phase": "observability",
            "run_id": paths.run_id,
            "format": args.format,
            "locus_gap_bp": args.locus_gap_bp,
            "expected_genes": args.expected_genes,
            "expected_projections": args.expected_projections,
            "expected_gff3": getattr(args, "expected_gff3", None),
            "expected_proteins": getattr(args, "expected_proteins", None),
            "reference_protein_hits": getattr(args, "reference_protein_hits", None),
            "gffcompare_tmap": getattr(args, "gffcompare_tmap", None),
            "busco_complete_percent": getattr(args, "busco_complete_percent", None),
            "generated_expected_inputs": generated_expected,
            "outputs": {
                "audit_loci": len(audit_loci),
                "gene_to_locus": len(gene_to_locus),
                "evidence_fate": len(evidence_fate),
                "expected_gene_presence": len(expected_presence),
                "source_profile": len(source_profile),
                "expected_source_profile": len(expected_source_profile),
                "busco_expected_crosswalk": len(busco_crosswalk),
                "completeness_profile": len(completeness_profile),
                "non_busco_high_confidence_losses": len(non_busco_high_confidence_losses),
                "busco_proxy_calibration": len(busco_proxy_calibration),
                "copy_number_audit": len(copy_number_audit),
                "biotype_transition": len(biotype_transition),
                "recommendations": len(recommendations),
                "same_assembly_structure": len(same_assembly_structure),
                "reference_protein_audit": len(reference_protein_audit),
                "release_readiness": len(release_readiness),
                "feature_profile": len(feature_profile),
                "failure_mode_summary": len(failure_summary),
                "review_loci": len(review_loci),
                "review_beds": len(review_bed_paths),
            },
        },
        os.path.join(paths.root, "manifest.observability.json"),
    )
    print(f"[observability] wrote {len(review_loci)} review loci -> {paths.observability_dir}")
    print(f"[observability] summary: {os.path.join(paths.observability_dir, 'actionable_summary.md')}")
    return 0


def run_end_to_end(args: argparse.Namespace) -> int:
    if not args.run_id:
        args.run_id = make_run_id(args.core_db, args.layer_db)
    paths = _paths_from_args(args)
    ensure_dir(paths.extract_dir)
    ensure_dir(paths.loci_dir)

    core_params = DBParams(args.host, args.port, args.user, args.password, args.core_db)
    layer_params = DBParams(args.host, args.port, args.user, args.password, args.layer_db)

    with connect(core_params) as core_conn, connect(layer_params) as layer_conn:
        core_regions = set(list_seq_regions(core_conn, args.coord_system_name))
        layer_regions = set(list_seq_regions(layer_conn, args.coord_system_name))
        seq_regions = sorted(core_regions.union(layer_regions))

    print(f"[observability-run] run_id={paths.run_id}")
    print(f"[observability-run] extracting {len(seq_regions)} seq_regions")

    with connect(core_params) as core_conn:
        core_genes = extract_all_genes(core_conn, seq_regions, args.coord_system_name)
        core_tx = extract_all_transcripts(core_conn, seq_regions, args.coord_system_name)
        core_tr = extract_all_translations(core_conn)
    with connect(layer_params) as layer_conn:
        layer_genes = extract_all_genes(layer_conn, seq_regions, args.coord_system_name)
        layer_tx = extract_all_transcripts(layer_conn, seq_regions, args.coord_system_name)
        layer_tr = extract_all_translations(layer_conn)

    _write_table(core_genes, os.path.join(paths.extract_dir, f"core_genes.{args.format}"), args.format)
    _write_table(layer_genes, os.path.join(paths.extract_dir, f"layer_genes.{args.format}"), args.format)
    _write_table(core_tx, os.path.join(paths.extract_dir, f"core_transcripts.{args.format}"), args.format)
    _write_table(layer_tx, os.path.join(paths.extract_dir, f"layer_transcripts.{args.format}"), args.format)
    _write_table(core_tr, os.path.join(paths.extract_dir, f"core_translations.{args.format}"), args.format)
    _write_table(layer_tr, os.path.join(paths.extract_dir, f"layer_translations.{args.format}"), args.format)

    loci, expanded, gene_to_locus = build_loci(core_genes, layer_genes, args.locus_gap_bp)
    _write_table(loci, os.path.join(paths.loci_dir, f"loci.strict.{args.format}"), args.format)
    _write_table(expanded, os.path.join(paths.loci_dir, f"loci.expanded.{args.format}"), args.format)
    _write_table(gene_to_locus, os.path.join(paths.loci_dir, f"gene_to_locus.{args.format}"), args.format)

    generated_expected = populate_expected_inputs(args, paths, seq_regions)

    write_manifest(
        {
            "phase": "observability_run_extract_loci",
            "run_id": paths.run_id,
            "core_db": args.core_db,
            "layer_db": args.layer_db,
            "expected_core_db": getattr(args, "expected_core_db", None),
            "expected_gff3": getattr(args, "expected_gff3", None),
            "coord_system_name": args.coord_system_name,
            "seq_regions": seq_regions,
            "format": args.format,
            "generated_expected_inputs": generated_expected,
            "rows": {
                "core_genes": len(core_genes),
                "layer_genes": len(layer_genes),
                "core_transcripts": len(core_tx),
                "layer_transcripts": len(layer_tx),
                "core_translations": len(core_tr),
                "layer_translations": len(layer_tr),
                "loci": len(loci),
                "gene_to_locus": len(gene_to_locus),
            },
        },
        os.path.join(paths.root, "manifest.observability.run.json"),
    )
    return run_audit(args)


def expected_tables_from_core_gene_table(
    genes: pd.DataFrame,
    expected_source: str,
    confidence: str,
    projection_mode: str = "none",
    target_seq_regions: Optional[Iterable[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert an Ensembl core gene table into expected-gene audit inputs.

    This intentionally does not infer orthology, BUSCO ids, symbols, or real
    assembly projection. In `same_coordinates` mode it copies coordinates only
    for genes whose seq_region name exists in the audited target set.
    """
    norm = _normalise_gene_table(genes, "expected")
    target_region_set = {str(x) for x in target_seq_regions or []}
    expected_rows: list[dict] = []
    projection_rows: list[dict] = []
    used_ids: dict[str, int] = {}

    for _, row in norm.iterrows():
        reference_stable_id = _string_value(row.get("stable_id"))
        gene_id = _string_value(row.get("gene_id"))
        base_id = reference_stable_id or f"{expected_source}:gene:{gene_id}"
        if not base_id or base_id == f"{expected_source}:gene:":
            base_id = f"{expected_source}:row:{len(expected_rows) + 1}"
        seen = used_ids.get(base_id, 0)
        used_ids[base_id] = seen + 1
        expected_gene_id = base_id if seen == 0 else f"{base_id}.{seen + 1}"

        expected_rows.append(
            {
                "expected_gene_id": expected_gene_id,
                "expected_source": expected_source,
                "reference_stable_id": reference_stable_id,
                "symbol": "",
                "biotype": _string_value(row.get("biotype"), "unknown"),
                "orthogroup_id": "",
                "busco_id": "",
                "expected_copy_number": 1,
                "confidence": confidence,
            }
        )

        if projection_mode == "same_coordinates":
            seq_region_name = _string_value(row.get("seq_region_name"))
            seq_region_is_target = not target_region_set or seq_region_name in target_region_set
            projection_rows.append(
                {
                    "expected_gene_id": expected_gene_id,
                    "seq_region_name": seq_region_name if seq_region_is_target else "",
                    "seq_region_start": _int_value(row.get("seq_region_start")) if seq_region_is_target else 0,
                    "seq_region_end": _int_value(row.get("seq_region_end")) if seq_region_is_target else 0,
                    "seq_region_strand": _normalise_strand(row.get("seq_region_strand"))
                    if seq_region_is_target
                    else 0,
                    "projection_status": "mapped" if seq_region_is_target else "unmapped",
                    "projection_identity": 1.0 if seq_region_is_target else 0.0,
                    "projection_coverage": 1.0 if seq_region_is_target else 0.0,
                    "assembly_gap_overlap_bp": 0,
                    "repeat_overlap_bp": 0,
                }
            )

    expected = pd.DataFrame(expected_rows, columns=EXPECTED_GENE_COLUMNS)
    projections = pd.DataFrame(projection_rows, columns=EXPECTED_PROJECTION_COLUMNS)
    return expected, projections


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_gff3_attributes(raw: object) -> dict[str, str]:
    text = _string_value(raw)
    attrs: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            key, value = item, ""
        attrs[unquote(key)] = unquote(value)
    return attrs


def read_gff3_features(path: str) -> pd.DataFrame:
    columns = [
        "seq_region_name",
        "source",
        "feature_type",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "attributes",
        "feature_id",
        "parent_id",
        "name",
        "gene_id",
        "symbol",
        "biotype",
        "orthogroup_id",
        "busco_id",
    ]
    rows: list[dict] = []
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seq, source, feature_type, start, end, _score, strand, _phase, raw_attrs = parts
            attrs = parse_gff3_attributes(raw_attrs)
            feature_id = (
                attrs.get("ID")
                or attrs.get("gene_id")
                or attrs.get("transcript_id")
                or attrs.get("Name")
                or f"gff3:{feature_type}:{seq}:{start}:{end}:{strand}:{line_number}"
            )
            gene_id = attrs.get("gene_id") or attrs.get("gene") or attrs.get("locus_tag") or feature_id
            symbol = attrs.get("gene_name") or attrs.get("Name") or attrs.get("symbol") or ""
            biotype = (
                attrs.get("biotype")
                or attrs.get("gene_biotype")
                or attrs.get("transcript_biotype")
                or attrs.get("gene_type")
                or attrs.get("transcript_type")
                or attrs.get("gbkey")
                or ("protein_coding" if feature_type in {"mRNA", "CDS"} else "unknown")
            )
            busco_id = attrs.get("busco_id") or attrs.get("BUSCO") or ""
            rows.append(
                {
                    "seq_region_name": seq,
                    "source": source,
                    "feature_type": feature_type,
                    "seq_region_start": _int_value(start),
                    "seq_region_end": _int_value(end),
                    "seq_region_strand": _normalise_strand(strand),
                    "attributes": raw_attrs,
                    "feature_id": feature_id,
                    "parent_id": attrs.get("Parent", ""),
                    "name": attrs.get("Name", ""),
                    "gene_id": gene_id,
                    "symbol": symbol,
                    "biotype": biotype,
                    "orthogroup_id": attrs.get("orthogroup_id", "") or attrs.get("orthogroup", ""),
                    "busco_id": busco_id,
                }
            )
    return pd.DataFrame(rows, columns=columns)


def _gff3_expected_source_row(group: pd.DataFrame) -> pd.Series:
    gene_rows = group[group["feature_type"] == "gene"]
    if not gene_rows.empty:
        return gene_rows.iloc[0]
    return group.sort_values(["seq_region_start", "seq_region_end"]).iloc[0]


def expected_tables_from_gff3(
    path: str,
    expected_source: str,
    confidence: str,
    projection_mode: str = "none",
    target_seq_regions: Optional[Iterable[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    features = read_gff3_features(path)
    target_region_set = {str(x) for x in target_seq_regions or []}
    expected_rows: list[dict] = []
    projection_rows: list[dict] = []

    gene_features = features[features["feature_type"] == "gene"].copy()
    if not gene_features.empty:
        groups = [(row["feature_id"], pd.DataFrame([row])) for _, row in gene_features.iterrows()]
    else:
        transcript_types = {"mRNA", "transcript", "lnc_RNA", "ncRNA", "pseudogenic_transcript", "rRNA", "tRNA"}
        transcript_features = features[features["feature_type"].isin(transcript_types)].copy()
        if not transcript_features.empty:
            transcript_features["_expected_group_id"] = transcript_features.apply(
                lambda r: _string_value(r.get("parent_id")) or _string_value(r.get("feature_id")),
                axis=1,
            )
            groups = list(transcript_features.groupby("_expected_group_id", dropna=False))
        else:
            cds_features = features[features["feature_type"] == "CDS"].copy()
            cds_features["_expected_group_id"] = cds_features.apply(
                lambda r: _string_value(r.get("parent_id")) or _string_value(r.get("feature_id")),
                axis=1,
            )
            groups = list(cds_features.groupby("_expected_group_id", dropna=False))

    used_ids: dict[str, int] = {}
    for raw_group_id, group in groups:
        if group.empty:
            continue
        row = _gff3_expected_source_row(group)
        group_id = _string_value(raw_group_id) or _string_value(row.get("feature_id"))
        reference_id = _string_value(row.get("gene_id")) or group_id
        base_id = reference_id or group_id or f"gff3:{len(expected_rows) + 1}"
        seen = used_ids.get(base_id, 0)
        used_ids[base_id] = seen + 1
        expected_gene_id = base_id if seen == 0 else f"{base_id}.{seen + 1}"

        start = int(group["seq_region_start"].astype(int).min())
        end = int(group["seq_region_end"].astype(int).max())
        seq_region_name = _string_value(row.get("seq_region_name"))
        seq_region_is_target = not target_region_set or seq_region_name in target_region_set
        expected_rows.append(
            {
                "expected_gene_id": expected_gene_id,
                "expected_source": expected_source,
                "reference_stable_id": reference_id,
                "symbol": _string_value(row.get("symbol")),
                "biotype": _string_value(row.get("biotype"), "unknown"),
                "orthogroup_id": _string_value(row.get("orthogroup_id")),
                "busco_id": _string_value(row.get("busco_id")),
                "expected_copy_number": 1,
                "confidence": confidence,
            }
        )
        if projection_mode == "same_coordinates":
            projection_rows.append(
                {
                    "expected_gene_id": expected_gene_id,
                    "seq_region_name": seq_region_name if seq_region_is_target else "",
                    "seq_region_start": start if seq_region_is_target else 0,
                    "seq_region_end": end if seq_region_is_target else 0,
                    "seq_region_strand": _normalise_strand(row.get("seq_region_strand"))
                    if seq_region_is_target
                    else 0,
                    "projection_status": "mapped" if seq_region_is_target else "unmapped",
                    "projection_identity": 1.0 if seq_region_is_target else 0.0,
                    "projection_coverage": 1.0 if seq_region_is_target else 0.0,
                    "assembly_gap_overlap_bp": 0,
                    "repeat_overlap_bp": 0,
                }
            )

    expected = pd.DataFrame(expected_rows, columns=EXPECTED_GENE_COLUMNS)
    projections = pd.DataFrame(projection_rows, columns=EXPECTED_PROJECTION_COLUMNS)
    return expected, projections


def _expected_db_params(args: argparse.Namespace) -> DBParams:
    return DBParams(
        getattr(args, "expected_host", None) or args.host,
        getattr(args, "expected_port", None) or args.port,
        getattr(args, "expected_user", None) or args.user,
        args.password if getattr(args, "expected_password", None) is None else args.expected_password,
        args.expected_core_db,
    )


def populate_expected_inputs_from_core_db(
    args: argparse.Namespace,
    paths: AuditPaths,
    target_seq_regions: Iterable[str],
) -> dict:
    if not getattr(args, "expected_core_db", None):
        return {"mode": "not_requested"}

    expected_dir = os.path.join(paths.root, "expected")
    ensure_dir(expected_dir)
    generated = {
        "mode": "from_core_db",
        "expected_core_db": args.expected_core_db,
        "expected_source_name": args.expected_source_name,
        "expected_confidence": args.expected_confidence,
        "expected_projection_mode": args.expected_projection_mode,
        "expected_genes": args.expected_genes,
        "expected_projections": args.expected_projections,
        "rows": {},
    }

    if args.expected_genes:
        generated["mode"] = "not_generated_user_supplied_expected_genes"
        print(
            "[observability-run] using supplied --expected_genes; "
            "--expected_core_db will not overwrite it"
        )
        return generated

    expected_coord_system = args.expected_coord_system_name or args.coord_system_name
    expected_params = _expected_db_params(args)
    with connect(expected_params) as expected_conn:
        expected_regions = list_seq_regions(expected_conn, expected_coord_system)
        expected_genes_raw = extract_all_genes(expected_conn, expected_regions, expected_coord_system)

    expected_genes, expected_projections = expected_tables_from_core_gene_table(
        expected_genes_raw,
        expected_source=args.expected_source_name,
        confidence=args.expected_confidence,
        projection_mode=args.expected_projection_mode,
        target_seq_regions=target_seq_regions,
    )

    suffix = args.format
    expected_genes_path = os.path.join(expected_dir, f"expected_genes.{suffix}")
    _write_table(expected_genes, expected_genes_path, args.format)
    args.expected_genes = expected_genes_path
    generated["expected_genes"] = expected_genes_path
    generated["rows"]["expected_genes"] = len(expected_genes)
    print(
        f"[observability-run] generated {len(expected_genes)} expected genes "
        f"from {args.expected_core_db} -> {expected_genes_path}"
    )

    if args.expected_projection_mode == "same_coordinates":
        if args.expected_projections:
            generated["mode"] = "from_core_db_with_user_supplied_projections"
        else:
            expected_projections_path = os.path.join(expected_dir, f"expected_projections.{suffix}")
            _write_table(expected_projections, expected_projections_path, args.format)
            args.expected_projections = expected_projections_path
            generated["expected_projections"] = expected_projections_path
            generated["rows"]["expected_projections"] = len(expected_projections)
            unmapped_count = int((expected_projections["projection_status"] == "unmapped").sum())
            print(
                f"[observability-run] generated {len(expected_projections)} same-coordinate projections "
                f"({unmapped_count} unmapped seq_region-name mismatches) -> {expected_projections_path}"
            )
    elif not args.expected_projections:
        generated["rows"]["expected_projections"] = 0
        print(
            "[observability-run] no expected projections generated; expected genes will be "
            "classified as projection_unmapped until a projection file is supplied"
        )

    return generated


def populate_expected_inputs_from_gff3(
    args: argparse.Namespace,
    paths: AuditPaths,
    target_seq_regions: Iterable[str],
) -> dict:
    if not getattr(args, "expected_gff3", None):
        return {"mode": "not_requested"}

    expected_dir = os.path.join(paths.root, "expected")
    ensure_dir(expected_dir)
    generated = {
        "mode": "from_gff3",
        "expected_gff3": args.expected_gff3,
        "expected_gff3_source_name": args.expected_gff3_source_name,
        "expected_gff3_confidence": args.expected_gff3_confidence,
        "expected_gff3_projection_mode": args.expected_gff3_projection_mode,
        "expected_genes": args.expected_genes,
        "expected_projections": args.expected_projections,
        "rows": {},
    }

    if args.expected_genes:
        generated["mode"] = "not_generated_user_supplied_expected_genes"
        print(
            "[observability-run] using supplied --expected_genes; "
            "--expected_gff3 will not overwrite it"
        )
        return generated

    expected_genes, expected_projections = expected_tables_from_gff3(
        args.expected_gff3,
        expected_source=args.expected_gff3_source_name,
        confidence=args.expected_gff3_confidence,
        projection_mode=args.expected_gff3_projection_mode,
        target_seq_regions=target_seq_regions,
    )
    suffix = args.format
    expected_genes_path = os.path.join(expected_dir, f"expected_genes.{suffix}")
    _write_table(expected_genes, expected_genes_path, args.format)
    args.expected_genes = expected_genes_path
    generated["expected_genes"] = expected_genes_path
    generated["rows"]["expected_genes"] = len(expected_genes)
    print(
        f"[observability-run] generated {len(expected_genes)} expected genes "
        f"from {args.expected_gff3} -> {expected_genes_path}"
    )

    if args.expected_gff3_projection_mode == "same_coordinates":
        if args.expected_projections:
            generated["mode"] = "from_gff3_with_user_supplied_projections"
        else:
            expected_projections_path = os.path.join(expected_dir, f"expected_projections.{suffix}")
            _write_table(expected_projections, expected_projections_path, args.format)
            args.expected_projections = expected_projections_path
            generated["expected_projections"] = expected_projections_path
            generated["rows"]["expected_projections"] = len(expected_projections)
            unmapped_count = int((expected_projections["projection_status"] == "unmapped").sum())
            print(
                f"[observability-run] generated {len(expected_projections)} GFF3 same-coordinate projections "
                f"({unmapped_count} unmapped seq_region-name mismatches) -> {expected_projections_path}"
            )
    elif not args.expected_projections:
        generated["rows"]["expected_projections"] = 0
        print(
            "[observability-run] no GFF3 projections generated; expected genes will be "
            "classified as projection_unmapped until a projection file is supplied"
        )
    return generated


def populate_expected_inputs(
    args: argparse.Namespace,
    paths: AuditPaths,
    target_seq_regions: Iterable[str],
) -> dict:
    if getattr(args, "expected_genes", None):
        return {"mode": "user_supplied_expected_genes", "expected_genes": args.expected_genes}
    if getattr(args, "expected_gff3", None):
        return populate_expected_inputs_from_gff3(args, paths, target_seq_regions)
    if getattr(args, "expected_core_db", None):
        return populate_expected_inputs_from_core_db(args, paths, target_seq_regions)
    return {"mode": "not_requested"}


EXPECTED_GENE_COLUMNS = [
    "expected_gene_id",
    "expected_source",
    "reference_stable_id",
    "symbol",
    "biotype",
    "orthogroup_id",
    "busco_id",
    "expected_copy_number",
    "confidence",
]

EXPECTED_PROJECTION_COLUMNS = [
    "expected_gene_id",
    "seq_region_name",
    "seq_region_start",
    "seq_region_end",
    "seq_region_strand",
    "projection_status",
    "projection_identity",
    "projection_coverage",
    "assembly_gap_overlap_bp",
    "repeat_overlap_bp",
]


def init_expected_template(args: argparse.Namespace) -> int:
    ensure_dir(args.out_dir)
    genes_path = os.path.join(args.out_dir, "expected_genes.tsv")
    projections_path = os.path.join(args.out_dir, "expected_projections.tsv")
    proteins_path = os.path.join(args.out_dir, "expected_proteins.tsv")
    protein_hits_path = os.path.join(args.out_dir, "reference_protein_hits.tsv")
    pd.DataFrame(columns=EXPECTED_GENE_COLUMNS).to_csv(genes_path, sep="\t", index=False)
    pd.DataFrame(columns=EXPECTED_PROJECTION_COLUMNS).to_csv(projections_path, sep="\t", index=False)
    pd.DataFrame(columns=EXPECTED_PROTEIN_COLUMNS).to_csv(proteins_path, sep="\t", index=False)
    pd.DataFrame(columns=REFERENCE_PROTEIN_HIT_COLUMNS).to_csv(protein_hits_path, sep="\t", index=False)
    print(f"[observability] wrote {genes_path}")
    print(f"[observability] wrote {projections_path}")
    print(f"[observability] wrote {proteins_path}")
    print(f"[observability] wrote {protein_hits_path}")
    return 0


def _missing_columns(path: str, required_columns: list[str]) -> list[str]:
    df = _read_table(path)
    return [col for col in required_columns if col not in df.columns]


def validate_inputs(args: argparse.Namespace) -> int:
    errors = []
    try:
        paths = _paths_from_args(args)
    except ValueError as exc:
        errors.append(str(exc))
        paths = None

    if paths is not None:
        for base in ("core_genes", "layer_genes"):
            try:
                _load_named_table(paths.extract_dir, base, args.format)
            except FileNotFoundError as exc:
                errors.append(str(exc))
        for base in ("core_transcripts", "layer_transcripts", "core_translations", "layer_translations"):
            try:
                _load_named_table(paths.extract_dir, base, args.format, required=False)
            except FileNotFoundError:
                pass

    if args.expected_genes:
        missing = _missing_columns(args.expected_genes, ["expected_gene_id", "confidence"])
        if missing:
            errors.append(f"{args.expected_genes} missing columns: {missing}")
    elif args.expected_projections:
        errors.append("--expected_projections requires --expected_genes")
    if args.expected_projections:
        missing = _missing_columns(
            args.expected_projections,
            ["expected_gene_id", "seq_region_name", "seq_region_start", "seq_region_end"],
        )
        if missing:
            errors.append(f"{args.expected_projections} missing columns: {missing}")
    if getattr(args, "expected_proteins", None):
        missing = _missing_columns(args.expected_proteins, ["expected_gene_id", "query_protein_id", "confidence"])
        if missing:
            errors.append(f"{args.expected_proteins} missing columns: {missing}")
    if getattr(args, "reference_protein_hits", None):
        missing = _missing_columns(args.reference_protein_hits, ["query_protein_id"])
        if missing:
            errors.append(f"{args.reference_protein_hits} missing columns: {missing}")
    if getattr(args, "gffcompare_tmap", None) and not os.path.exists(args.gffcompare_tmap):
        errors.append(f"Could not find gffcompare tmap: {args.gffcompare_tmap}")
    if getattr(args, "expected_gff3", None) and not os.path.exists(args.expected_gff3):
        errors.append(f"Could not find expected GFF3: {args.expected_gff3}")

    if errors:
        for error in errors:
            print(f"[observability-validate] ERROR: {error}", file=sys.stderr)
        return 2
    print("[observability-validate] inputs look usable")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m ensembl.genes.rgb.observability",
        description="Implemented draft genebuild observability audit over RGB extract tables.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Extract from DBs, build loci, and run the observability audit")
    run.add_argument("--host", required=True)
    run.add_argument("--port", type=int, default=3306)
    run.add_argument("--user", required=True)
    run.add_argument("--password", default="")
    run.add_argument("--core_db", required=True)
    run.add_argument("--layer_db", required=True)
    run.add_argument("--coord_system_name", default=None)
    run.add_argument("--output_dir", required=True, help="Output root for the generated run directory")
    run.add_argument("--run_id", help="Optional run id; defaults to an auto-generated id")
    run.add_argument("--format", choices=["tsv", "csv", "parquet"], default="tsv")
    run.add_argument("--expected_genes", help="Optional expected gene catalogue TSV/CSV/Parquet")
    run.add_argument("--expected_projections", help="Optional expected gene projection TSV/CSV/Parquet")
    run.add_argument("--expected_gff3", help="Optional expected annotation GFF3/GFF3.gz")
    run.add_argument("--expected_proteins", help="Optional expected reference protein catalogue TSV/CSV/Parquet")
    run.add_argument("--reference_protein_hits", help="Optional BUSCO-like reference protein hit table")
    run.add_argument("--gffcompare_tmap", help="Optional GffCompare .tmap for same-assembly structure comparison")
    run.add_argument(
        "--expected_core_db",
        help=(
            "Optional reference core DB used to auto-generate expected_genes when "
            "--expected_genes is not supplied"
        ),
    )
    run.add_argument("--expected_host", help="Host for --expected_core_db; defaults to --host")
    run.add_argument("--expected_port", type=int, help="Port for --expected_core_db; defaults to --port")
    run.add_argument("--expected_user", help="User for --expected_core_db; defaults to --user")
    run.add_argument(
        "--expected_password",
        help="Password for --expected_core_db; defaults to --password",
    )
    run.add_argument(
        "--expected_coord_system_name",
        help="Coord-system name for --expected_core_db; defaults to --coord_system_name",
    )
    run.add_argument(
        "--expected_source_name",
        default="prior_ensembl",
        help="expected_source value for auto-generated expected genes",
    )
    run.add_argument(
        "--expected_confidence",
        choices=["high", "medium", "low", "unknown"],
        default="high",
        help="confidence value for auto-generated expected genes",
    )
    run.add_argument(
        "--expected_projection_mode",
        choices=["none", "same_coordinates"],
        default="none",
        help=(
            "How to auto-populate expected projections. Use same_coordinates only "
            "when --expected_core_db is already on the audited assembly coordinates"
        ),
    )
    run.add_argument(
        "--expected_gff3_source_name",
        default="gff3_annotation",
        help="expected_source value for auto-generated expected genes from --expected_gff3",
    )
    run.add_argument(
        "--expected_gff3_confidence",
        choices=["high", "medium", "low", "unknown"],
        default="high",
        help="confidence value for auto-generated expected genes from --expected_gff3",
    )
    run.add_argument(
        "--expected_gff3_projection_mode",
        choices=["none", "same_coordinates"],
        default="none",
        help=(
            "How to auto-populate projections from --expected_gff3. Use same_coordinates only "
            "when the GFF3 coordinates are already on the audited assembly"
        ),
    )
    run.add_argument("--locus_gap_bp", type=int, default=5000)
    run.add_argument("--protein_min_identity", type=float, default=0.30)
    run.add_argument("--protein_min_query_coverage", type=float, default=0.70)
    run.add_argument("--busco_complete_percent", type=float, help="Optional external BUSCO complete percentage")
    run.add_argument("--busco_floor_percent", type=float, default=95.0)
    run.add_argument("--max_high_confidence_actionable_loss_fraction", type=float, default=0.0)
    run.add_argument("--max_copy_number_issue_fraction", type=float, default=0.0)
    run.add_argument(
        "--require_expected_genes",
        action="store_true",
        help="Fail release readiness if no expected-gene catalogue was audited",
    )
    run.add_argument("--top_n", type=int, default=500)
    run.add_argument("--bed_pad_bp", type=int, default=0)
    run.set_defaults(func=run_end_to_end)

    audit = sub.add_parser("audit", help="Run audit over existing RGB extract/loci tables")
    audit.add_argument("--run_dir", help="Existing RGB run directory containing extract/ and optional loci/")
    audit.add_argument("--output_dir", help="RGB output root containing <run_id>/extract")
    audit.add_argument("--run_id", help="RGB run id")
    audit.add_argument("--format", choices=["tsv", "csv", "parquet"], default="tsv")
    audit.add_argument("--expected_genes", help="Optional expected gene catalogue TSV/CSV/Parquet")
    audit.add_argument("--expected_projections", help="Optional expected gene projection TSV/CSV/Parquet")
    audit.add_argument("--expected_gff3", help="Optional expected annotation GFF3/GFF3.gz")
    audit.add_argument(
        "--expected_gff3_source_name",
        default="gff3_annotation",
        help="expected_source value for auto-generated expected genes from --expected_gff3",
    )
    audit.add_argument(
        "--expected_gff3_confidence",
        choices=["high", "medium", "low", "unknown"],
        default="high",
        help="confidence value for auto-generated expected genes from --expected_gff3",
    )
    audit.add_argument(
        "--expected_gff3_projection_mode",
        choices=["none", "same_coordinates"],
        default="none",
        help=(
            "How to auto-populate projections from --expected_gff3. Use same_coordinates only "
            "when the GFF3 coordinates are already on the audited assembly"
        ),
    )
    audit.add_argument("--expected_proteins", help="Optional expected reference protein catalogue TSV/CSV/Parquet")
    audit.add_argument("--reference_protein_hits", help="Optional BUSCO-like reference protein hit table")
    audit.add_argument("--gffcompare_tmap", help="Optional GffCompare .tmap for same-assembly structure comparison")
    audit.add_argument("--locus_gap_bp", type=int, default=5000)
    audit.add_argument("--protein_min_identity", type=float, default=0.30)
    audit.add_argument("--protein_min_query_coverage", type=float, default=0.70)
    audit.add_argument("--busco_complete_percent", type=float, help="Optional external BUSCO complete percentage")
    audit.add_argument("--busco_floor_percent", type=float, default=95.0)
    audit.add_argument("--max_high_confidence_actionable_loss_fraction", type=float, default=0.0)
    audit.add_argument("--max_copy_number_issue_fraction", type=float, default=0.0)
    audit.add_argument(
        "--require_expected_genes",
        action="store_true",
        help="Fail release readiness if no expected-gene catalogue was audited",
    )
    audit.add_argument("--top_n", type=int, default=500)
    audit.add_argument("--bed_pad_bp", type=int, default=0)
    audit.set_defaults(func=run_audit)

    validate = sub.add_parser("validate-inputs", help="Check existing extract and expected-gene inputs")
    validate.add_argument("--run_dir", help="Existing RGB run directory containing extract/ and optional loci/")
    validate.add_argument("--output_dir", help="RGB output root containing <run_id>/extract")
    validate.add_argument("--run_id", help="RGB run id")
    validate.add_argument("--format", choices=["tsv", "csv", "parquet"], default="tsv")
    validate.add_argument("--expected_genes", help="Optional expected gene catalogue TSV/CSV/Parquet")
    validate.add_argument("--expected_projections", help="Optional expected gene projection TSV/CSV/Parquet")
    validate.add_argument("--expected_gff3", help="Optional expected annotation GFF3/GFF3.gz")
    validate.add_argument("--expected_proteins", help="Optional expected reference protein catalogue TSV/CSV/Parquet")
    validate.add_argument("--reference_protein_hits", help="Optional BUSCO-like reference protein hit table")
    validate.add_argument("--gffcompare_tmap", help="Optional GffCompare .tmap for same-assembly structure comparison")
    validate.set_defaults(func=validate_inputs)

    template = sub.add_parser("init-expected-template", help="Write empty expected gene/projection TSV templates")
    template.add_argument("--out_dir", required=True)
    template.set_defaults(func=init_expected_template)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

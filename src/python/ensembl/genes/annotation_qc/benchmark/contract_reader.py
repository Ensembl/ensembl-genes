"""
Read the files of one comparator result written under output contract 1.0.

Nothing here computes a metric. Values are read as written: counts as integers,
``value`` as the float the comparator wrote (unrounded; the text is kept as well),
``NA`` as None, and every status and reason verbatim. Files are read by column
name; a missing column required by the dashboard is an error naming the file.
Checksums are verified before import (results.verify_outputs), not here.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise.contract import METRIC_COLUMNS

METRICS = "metrics.tsv"
DETAILS = "comparison_details.tsv"
PAIRS = "cds_exact_one_to_one_pairs.tsv"
EVALUATED = "evaluated_transcripts.tsv.gz"
PER_SEQUENCE = "per_sequence.tsv"
SPLICE_POSITION = "splice_sites_by_position.tsv"
INTRON_LENGTH = "cds_introns_by_length.tsv"
STOP_AUDIT = "stop_codon_harmonisation.tsv"
MANIFEST = "comparison_manifest.json"

DETAIL_COLUMNS = (
    "source",
    "gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "classification",
    "classification_cds",
    "matched_id",
    "best_match_transcript_id",
    "exon_overlap",
    "intron_chain_match",
    "cds_overlap",
    "cds_intron_chain_match",
    "original_gene_id",
    "exon_coordinate_exact",
    "cds_coordinate_exact",
    "cds_matched_id",
    "best_cds_match_transcript_id",
    "best_match_own_transcript_id",
    "best_cds_match_own_transcript_id",
    "strand_mismatch_basis",
    "novel_category",
    "counterpart_count",
    "counterpart_ids",
    "gene_biotype",
    "evaluated_transcript_count",
    "max_exon_count",
    "max_cds_count",
    "has_cds",
    "cds_shared_bases",
    "cds_locus_recovered",
    "cds_outcome",
    "one_to_one_partner",
    "cds_intron_chain_any_pair",
    "intron_chain_any_pair",
    "cds_start_diff_bp",
    "cds_stop_diff_bp",
    "cds_introns_missing",
    "cds_introns_extra",
    "scope_ids",
)


class ContractError(ValueError):
    """A contract file lacks a required column or cannot be read."""


def read_tsv(path: Path, compressed: bool = False) -> pd.DataFrame:
    opener = gzip.open if compressed else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
        return pd.read_csv(handle, sep="\t", dtype=str, keep_default_na=False)


def _require(frame: pd.DataFrame, columns, name: str) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ContractError(f"{name}: missing columns {missing}")


def int_or_none(text):
    if text in ("", "NA", None):
        return None
    return int(float(text)) if "." in str(text) else int(text)


def float_or_none(text):
    if text in ("", "NA", None):
        return None
    return float(text)


def manifest(out_dir: Path) -> dict:
    return json.loads((Path(out_dir) / MANIFEST).read_text())


def metrics(out_dir: Path) -> list[dict]:
    """Rows of metrics.tsv with typed numbers; ``value_text`` keeps the written value."""
    frame = read_tsv(Path(out_dir) / METRICS)
    _require(frame, METRIC_COLUMNS, METRICS)
    rows = []
    for r in frame.to_dict("records"):
        rows.append(
            {
                "scope_id": r["scope_id"],
                "metric_id": r["metric_id"],
                "numerator": int_or_none(r["numerator"]),
                "denominator": int_or_none(r["denominator"]),
                "value": float_or_none(r["value"]),
                "value_text": r["value"],
                "status": r["status"],
                "reason": r["reason"],
                "numerator_id": r["numerator_id"],
                "denominator_id": r["denominator_id"],
                "unit": r["unit"],
                "side": r["side"],
                "source": "contract",
            }
        )
    return rows


def details(out_dir: Path) -> pd.DataFrame:
    frame = read_tsv(Path(out_dir) / DETAILS)
    _require(frame, DETAIL_COLUMNS, DETAILS)
    return frame


def pairs(out_dir: Path) -> pd.DataFrame:
    frame = read_tsv(Path(out_dir) / PAIRS)
    _require(
        frame,
        (
            "reference_gene_id",
            "query_gene_id",
            "chrom",
            "strand",
            "cds_start",
            "cds_end",
        ),
        PAIRS,
    )
    return frame


def evaluated(out_dir: Path) -> pd.DataFrame:
    frame = read_tsv(Path(out_dir) / EVALUATED, compressed=True)
    _require(
        frame,
        (
            "role",
            "gene_id",
            "transcript_id",
            "original_gene_id",
            "original_transcript_id",
            "chrom",
            "strand",
            "start",
            "end",
            "status",
            "exclusion_reason",
            "selected_by",
            "exon_count",
            "cds_segment_count",
            "cds_length",
            "transformed",
            "input_sha256",
            "canonical_tag",
            "transcript_biotype",
        ),
        EVALUATED,
    )
    return frame


def optional(out_dir: Path, name: str) -> pd.DataFrame | None:
    path = Path(out_dir) / name
    return read_tsv(path) if path.exists() else None

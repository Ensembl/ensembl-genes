"""
Explicit importer for comparator results written before output contract 1.0.

Pre-contract results (e.g. the human/mouse benchmark of 2026-10-05 and the pig
comparisons of 2026-10-06) have no ``metrics.tsv``. This module maps only fields
of ``comparison_summary.json`` whose definitions are identical to a contract
metric (``LEGACY_FIELDS``) onto that contract metric id, with the integer counts
as numerator and denominator. Every other contract metric is ``not_available``
with reason ``pre_contract_result``. Nothing is derived from per-gene rows into a
contract metric: CDS outcome categories, partial buckets, multi/single-segment
splits, any-pair intron chains, analysis scopes and splice tables need the
contract and stay not available.

The earlier dashboard estimate of one-to-one exact TP (min of unique genes among
per-gene best-pair exact rows; not a maximum matching) is kept under separate
``historical_*`` ids. It never fills ``exact_cds_one_to_one_*``.

Each row: scope_id ("all" only), metric_id, numerator, denominator, value,
status, reason, numerator_id, denominator_id, unit, side, source.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise import contract

SOURCE = "legacy_summary"
REASON = "pre_contract_result"
HISTORICAL_IDS = (
    "historical_exact_tp_best_pair",
    "historical_exact_recall_best_pair",
    "historical_exact_precision_best_pair",
    "historical_exact_f1_best_pair",
)
NOVEL_CATEGORIES = (
    "Overlaps_reference_pseudogene",
    "Overlaps_other_non_protein_coding_reference",
    "Overlaps_unevaluated_protein_coding_reference",
    "Overlaps_evaluated_reference_opposite_strand",
    "Novel_no_reference_overlap",
)


def _get(summary: dict, path: str):
    node = summary
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


# contract metric id -> (numerator field, denominator: "R", "Q", "R_plus_Q",
# another summary field, or None for a count). Only identical definitions.
LEGACY_FIELDS: dict[str, tuple[str, str | None]] = {
    "reference_genes": ("total_reference_genes", None),
    "query_genes": ("total_consensus_genes", None),
    "reference_multi_cds_genes": (
        "sensitivity_cds.multi_segment_cds_reference_genes",
        None,
    ),
    "locus_recovery_cds_overlap": ("cds_overlap_locus_recovery.recovered_count", "R"),
    "cds_coordinate_exact_reference": (
        "sensitivity_cds.cds_coordinate_exact_count",
        "R",
    ),
    "cds_coordinate_exact_query": (
        "specificity.consensus_cds_coordinate_exact_count",
        "Q",
    ),
    "exact_cds_one_to_one_tp": ("cds_exact_one_to_one.true_positives", None),
    "exact_cds_one_to_one_recall": ("cds_exact_one_to_one.true_positives", "R"),
    "exact_cds_one_to_one_precision": ("cds_exact_one_to_one.true_positives", "Q"),
    "exact_cds_one_to_one_f1": ("cds_exact_one_to_one.true_positives", "R_plus_Q"),
    "exact_cds_gene_pairs": ("cds_exact_one_to_one.exact_gene_pairs", None),
    "cds_intron_chain_best_pair": (
        "sensitivity_cds.cds_intron_chain_recovered",
        "sensitivity_cds.multi_segment_cds_reference_genes",
    ),
    "ref_unique_cds_introns_recovered": (
        "intron_support.cds.shared_introns",
        "intron_support.cds.reference_introns",
    ),
    "query_unique_cds_introns_supported": (
        "intron_support.cds.shared_introns",
        "intron_support.cds.query_introns",
    ),
    "ref_unique_exon_introns_recovered": (
        "intron_support.exon.shared_introns",
        "intron_support.exon.reference_introns",
    ),
    "query_unique_exon_introns_supported": (
        "intron_support.exon.shared_introns",
        "intron_support.exon.query_introns",
    ),
    "cds_structural_exact_match": ("sensitivity_cds.cds_exact_match_count", "R"),
    "exon_coordinate_exact_reference": ("sensitivity.exon_coordinate_exact_count", "R"),
    "exon_structural_exact_match": ("sensitivity.exact_match_count", "R"),
    "locus_detection_gene_span": ("sensitivity.locus_detected_count", "R"),
    "span_only_detection": ("locus_detection_exonic.span_only_detected_count", "R"),
    "missed_reference_genes": ("sensitivity.missed_count", "R"),
    "strand_mismatch_reference": ("sensitivity.strand_mismatch_count", "R"),
    "query_novel": ("specificity.novel_consensus_count", "Q"),
    "strand_mismatch_query": ("specificity.strand_mismatch_consensus_count", "Q"),
    "split_reference_genes": ("split_merge.gene_split_count", "R"),
    "merge_query_genes": ("split_merge.gene_merge_count", "Q"),
    "one_counterpart_reference_genes": ("split_merge.one_to_one_reference_genes", "R"),
    **{f"query_novel_{c}": (f"novel_categories.{c}", "Q") for c in NOVEL_CATEGORIES},
    "stop_codon_transcripts_with_cds": (
        "query_stop_codon_harmonisation.transcripts_with_cds",
        None,
    ),
    "stop_codon_eligible": (
        "query_stop_codon_harmonisation.eligible",
        "query_stop_codon_harmonisation.transcripts_with_cds",
    ),
    "stop_codon_extended": (
        "query_stop_codon_harmonisation.extended",
        "query_stop_codon_harmonisation.transcripts_with_cds",
    ),
    "stop_codon_unchanged": (
        "query_stop_codon_harmonisation.unchanged",
        "query_stop_codon_harmonisation.transcripts_with_cds",
    ),
    "stop_codon_skipped": (
        "query_stop_codon_harmonisation.skipped",
        "query_stop_codon_harmonisation.transcripts_with_cds",
    ),
}


def _row(
    mid,
    numerator=None,
    denominator=None,
    value=None,
    status="ok",
    reason="",
    source=SOURCE,
):
    spec = contract.METRICS.get(mid)
    return {
        "scope_id": contract.ALL_SCOPE,
        "metric_id": mid,
        "numerator": numerator,
        "denominator": denominator,
        "value": value,
        "status": status,
        "reason": reason,
        "numerator_id": mid if spec else "",
        "denominator_id": (spec[3] if spec else "") or "",
        "unit": spec[1] if spec else "",
        "side": spec[0] if spec else "",
        "source": source,
    }


def legacy_metric_rows(summary: dict) -> list[dict]:
    """Contract metric rows for a pre-contract summary (see module docstring)."""
    R, Q = summary.get("total_reference_genes"), summary.get("total_consensus_genes")
    stop_applied = bool(
        (summary.get("query_stop_codon_harmonisation") or {}).get("applied")
    )
    rows = []
    for mid in contract.METRICS:
        field = LEGACY_FIELDS.get(mid)
        if field is None:
            rows.append(_row(mid, status="not_available", reason=REASON))
            continue
        numerator = _get(summary, field[0])
        if (
            mid.startswith("stop_codon_")
            and summary.get("query_stop_codon_harmonisation") is not None
            and not stop_applied
        ):
            rows.append(
                _row(mid, status="not_applicable", reason="transformation_not_applied")
            )
            continue
        if numerator is None:
            rows.append(
                _row(mid, status="not_available", reason=f"{REASON}:{field[0]}_absent")
            )
            continue
        den_spec = field[1]
        if den_spec is None:
            rows.append(_row(mid, numerator, None, numerator))
            continue
        denominator = {
            "R": R,
            "Q": Q,
            "R_plus_Q": (R + Q) if R is not None and Q is not None else None,
        }.get(
            den_spec,
            _get(summary, den_spec) if den_spec not in ("R", "Q", "R_plus_Q") else None,
        )
        if denominator is None:
            rows.append(
                _row(mid, status="not_available", reason=f"{REASON}:{den_spec}_absent")
            )
            continue
        if mid == "exact_cds_one_to_one_f1":
            numerator = 2 * numerator
        if denominator == 0:
            rows.append(_row(mid, numerator, 0, None, "undefined", "zero_denominator"))
        else:
            rows.append(_row(mid, numerator, denominator, numerator / denominator))
    return rows


def historical_rows(details: pd.DataFrame | None, R: int, Q: int) -> list[dict]:
    """The earlier best-pair estimate of one-to-one exact TP (historical ids only)."""
    note = (
        "historical best-pair estimate: min(unique reference genes, unique query genes) "
        "among per-gene coordinate-exact rows; not a maximum matching"
    )
    if details is None:
        return [
            _row(
                mid,
                status="not_available",
                reason="comparison_details.tsv absent",
                source="historical",
            )
            for mid in HISTORICAL_IDS
        ]
    ref = details[details.source == "reference"]
    exact = ref[ref.cds_coordinate_exact == "True"]
    tp = min(exact.gene_id.nunique(), exact.cds_matched_id.nunique())
    recall = tp / R if R else None
    precision = tp / Q if Q else None
    f1 = 2 * recall * precision / (recall + precision) if recall and precision else None

    def status(denominator_ok: bool) -> dict:
        return (
            {}
            if denominator_ok
            else {"status": "undefined", "reason": "zero_denominator"}
        )

    rows = [
        _row(HISTORICAL_IDS[0], tp, None, tp, source="historical"),
        _row(HISTORICAL_IDS[1], tp, R, recall, source="historical", **status(bool(R))),
        _row(
            HISTORICAL_IDS[2], tp, Q, precision, source="historical", **status(bool(Q))
        ),
        _row(
            HISTORICAL_IDS[3],
            None,
            None,
            f1,
            source="historical",
            **status(f1 is not None),
        ),
    ]
    for row in rows:
        row["reason"] = row["reason"] or note
        row["unit"] = "rate" if row["metric_id"] != HISTORICAL_IDS[0] else "gene_pairs"
    return rows


def read_details(out_dir: Path) -> pd.DataFrame | None:
    path = Path(out_dir) / "comparison_details.tsv"
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def legacy_result(out_dir: Path) -> tuple[list[dict], pd.DataFrame | None]:
    """(metric rows, per-gene details or None) of one pre-contract result."""
    summary = json.loads((Path(out_dir) / "comparison_summary.json").read_text())
    details = read_details(out_dir)
    rows = legacy_metric_rows(summary)
    rows += historical_rows(
        details,
        summary.get("total_reference_genes"),
        summary.get("total_consensus_genes"),
    )
    return rows, details


# Headline table column -> (metric_id, field) for `benchmark verify-headline` against
# the 2026-10-05 benchmark's headline_metrics.tsv (wide). Only identical-definition
# columns are mapped; the four one-to-one columns of that table were computed with
# the historical best-pair derivation and are checked against the historical ids.
HEADLINE_COLUMNS = {
    "reference_coding_genes": ("reference_genes", "numerator"),
    "predicted_coding_genes": ("query_genes", "numerator"),
    "locus_recovery_same_strand_cds_overlap": (
        "locus_recovery_cds_overlap",
        "numerator",
    ),
    "cds_coordinate_exact_ref_genes": ("cds_coordinate_exact_reference", "numerator"),
    "cds_coordinate_exact_rate": ("cds_coordinate_exact_reference", "value"),
    "cds_intron_chain_best_pair": ("cds_intron_chain_best_pair", "numerator"),
    "cds_intron_chain_denominator_multi_segment": (
        "cds_intron_chain_best_pair",
        "denominator",
    ),
    "cds_intron_chain_best_pair_rate": ("cds_intron_chain_best_pair", "value"),
    "ref_unique_cds_introns": ("ref_unique_cds_introns_recovered", "denominator"),
    "ref_unique_cds_introns_recovered": (
        "ref_unique_cds_introns_recovered",
        "numerator",
    ),
    "ref_unique_cds_introns_recovered_rate": (
        "ref_unique_cds_introns_recovered",
        "value",
    ),
    "query_unique_cds_introns": ("query_unique_cds_introns_supported", "denominator"),
    "query_cds_introns_supported_rate": ("query_unique_cds_introns_supported", "value"),
    "query_cds_coordinate_exact_genes": ("cds_coordinate_exact_query", "numerator"),
    "query_cds_coordinate_exact_rate": ("cds_coordinate_exact_query", "value"),
    "exon_coordinate_exact_ref_genes": ("exon_coordinate_exact_reference", "numerator"),
    "structural_exact_match_with_utr": ("exon_structural_exact_match", "numerator"),
    "gene_span_locus_detection": ("locus_detection_gene_span", "numerator"),
    "gene_span_only_no_exon_overlap": ("span_only_detection", "numerator"),
    "missed_reference_genes": ("missed_reference_genes", "numerator"),
    "missed_rate": ("missed_reference_genes", "value"),
    "unmatched_predictions_novel": ("query_novel", "numerator"),
    "unmatched_predictions_novel_rate": ("query_novel", "value"),
    **{f"novel_{c}": (f"query_novel_{c}", "numerator") for c in NOVEL_CATEGORIES},
    "splits_ref_genes": ("split_reference_genes", "numerator"),
    "merges_query_genes": ("merge_query_genes", "numerator"),
    "one_to_one_ref_genes": ("one_counterpart_reference_genes", "numerator"),
    "strand_mismatch_ref": ("strand_mismatch_reference", "numerator"),
    "strand_mismatch_query": ("strand_mismatch_query", "numerator"),
    "exact_tp_one_to_one": ("historical_exact_tp_best_pair", "numerator"),
    "gene_recall_exact_cds": ("historical_exact_recall_best_pair", "value"),
    "gene_precision_exact_cds": ("historical_exact_precision_best_pair", "value"),
    "gene_f1_exact_cds": ("historical_exact_f1_best_pair", "value"),
}
# The benchmark rounded rates to 4 dp, and computed F1 from rounded recall and precision.
HEADLINE_TOLERANCE = {"gene_f1_exact_cds": 1.5e-4}


def compare_headline_value(column: str, expected, derived) -> bool:
    if expected is None or (isinstance(expected, float) and math.isnan(expected)):
        return derived is None
    if derived is None:
        return False
    return abs(float(expected) - float(derived)) <= HEADLINE_TOLERANCE.get(
        column, 5e-9 if float(expected).is_integer() else 5e-5
    )

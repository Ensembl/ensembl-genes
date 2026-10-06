"""
Per-gene outcomes and metric rows from one pairwise-compare output directory.

Packaged from the benchmark's ``build_headline.py``. Every metric row carries its
numerator, denominator, value, status and source:

    comparator          read from comparison_summary.json
    comparator_per_gene counted from comparison_details.tsv rows (same run)
    reference_models    needs per-gene CDS segment counts, computed with the
                        comparator's own parser/selection on the reference
    independent         from an independent-analysis table, used only after it has
                        been verified gene-by-gene against the comparator output

Status values: ``ok``, ``not_available`` (input needed is absent; ``note`` says
which), never a silent zero. One-to-one exact-CDS precision/recall/F1 use pairs
from the per-gene table (each query gene used once); they are not built from
summary counts with different denominators.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

DETECTED = {"Exact_Match", "Partial_Match", "Structural_Mismatch"}
NOVEL_CATEGORIES = (
    "Overlaps_reference_pseudogene",
    "Overlaps_other_non_protein_coding_reference",
    "Overlaps_unevaluated_protein_coding_reference",
    "Overlaps_evaluated_reference_opposite_strand",
    "Novel_no_reference_overlap",
)

# cds_status values for reference genes (mutually exclusive)
CDS_STATUS = {
    "exact": "coordinate-exact CDS",
    "terminal_diff": "same CDS intron chain, start/stop differs",
    "diff_chain": "≥ 0.8 reciprocal CDS overlap, different intron chain",
    "partial_lt08": "shares CDS, < 0.8 reciprocal overlap",
    "span_only": "same-strand gene-span partner, no CDS overlap",
    "strand_mismatch": "opposite-strand CDS (or exon) overlap only",
    "missed": "no same-strand partner",
}


def _rate(num, den):
    return round(num / den, 4) if den else None


def _metric(
    mid, num=None, den=None, value=None, source="comparator", status="ok", note=None
):
    if status == "ok" and value is None and num is not None:
        value = _rate(num, den) if den is not None else num
    return {
        "metric_id": mid,
        "numerator": num,
        "denominator": den,
        "value": value,
        "status": status,
        "source": source,
        "note": note,
    }


def _na(mid, source, note, den=None):
    return _metric(mid, den=den, source=source, status="not_available", note=note)


def read_details(out_dir: Path) -> pd.DataFrame | None:
    path = Path(out_dir) / "comparison_details.tsv"
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def reference_outcomes(det: pd.DataFrame) -> pd.DataFrame:
    ref = det[det.source == "reference"].copy()
    ovl = ref.cds_overlap.astype(float)
    exact = ref.cds_coordinate_exact == "True"
    status = pd.Series("missed", index=ref.index)
    status[ref.classification == "Strand_Mismatch"] = "strand_mismatch"
    detected = ref.classification.isin(DETECTED)
    status[detected & (ovl == 0)] = "span_only"
    status[(ovl > 0) & (ref.classification_cds == "Partial_Match")] = "partial_lt08"
    status[(ovl > 0) & (ref.classification_cds == "Structural_Mismatch")] = "diff_chain"
    status[(ovl > 0) & (ref.classification_cds == "Exact_Match")] = "terminal_diff"
    status[exact] = "exact"
    return pd.DataFrame(
        {
            "gene_id": ref.gene_id,
            "chrom": ref.chrom,
            "start": ref.start.astype(int),
            "end": ref.end.astype(int),
            "strand": ref.strand,
            "classification": ref.classification,
            "classification_cds": ref.classification_cds,
            "cds_status": status,
            "cds_exact": exact.astype(int),
            "cds_overlap": ovl,
            "exon_overlap": ref.exon_overlap.astype(float),
            "cds_intron_chain_match": ref.cds_intron_chain_match,
            "intron_chain_match": ref.intron_chain_match,
            "exon_coordinate_exact": (ref.exon_coordinate_exact == "True").astype(int),
            "cds_matched_id": ref.cds_matched_id,
            "best_cds_ref_tx": ref.best_cds_match_own_transcript_id,
            "best_cds_query_tx": ref.best_cds_match_transcript_id,
            "matched_id": ref.matched_id,
            "best_match_ref_tx": ref.best_match_own_transcript_id,
            "best_match_query_tx": ref.best_match_transcript_id,
            "counterpart_count": ref.counterpart_count.astype(int),
            "counterpart_ids": ref.counterpart_ids,
            "strand_mismatch_basis": ref.strand_mismatch_basis,
        }
    ).reset_index(drop=True)


def query_outcomes(det: pd.DataFrame) -> pd.DataFrame:
    q = det[det.source == "consensus"].copy()
    return pd.DataFrame(
        {
            "gene_id": q.gene_id,
            "original_gene_id": q.original_gene_id,
            "chrom": q.chrom,
            "start": q.start.astype(int),
            "end": q.end.astype(int),
            "strand": q.strand,
            "classification": q.classification,
            "classification_cds": q.classification_cds,
            "cds_exact": (q.cds_coordinate_exact == "True").astype(int),
            "cds_overlap": q.cds_overlap.astype(float),
            "cds_matched_id": q.cds_matched_id,
            "best_cds_ref_tx": q.best_cds_match_transcript_id,
            "best_cds_query_tx": q.best_cds_match_own_transcript_id,
            "novel_category": q.novel_category,
            "counterpart_count": q.counterpart_count.astype(int),
            "counterpart_ids": q.counterpart_ids,
            "strand_mismatch_basis": q.strand_mismatch_basis,
        }
    ).reset_index(drop=True)


def load_independent(
    path: Path | None, ref_out: pd.DataFrame
) -> tuple[pd.DataFrame | None, dict]:
    """Load and verify an independent reference-gene table against comparator outcomes."""
    if path is None:
        return None, {"status": "absent"}
    path = Path(path)
    if not path.exists():
        return None, {"status": "missing", "path": str(path)}
    ind = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    a = set(ref_out.gene_id)
    b = set(ind.gene_id)
    check = {"path": str(path), "genes_comparator": len(a), "genes_independent": len(b)}
    if a != b:
        return None, {
            **check,
            "status": "rejected",
            "reason": f"gene sets differ ({len(a ^ b)} genes)",
        }
    joined = ref_out.set_index("gene_id").join(
        ind.set_index("gene_id")[["cds_exact_any"]]
    )
    disagree = int(((joined.cds_exact == 1) != (joined.cds_exact_any == "True")).sum())
    if disagree:
        return None, {
            **check,
            "status": "rejected",
            "reason": f"{disagree} genes disagree on coordinate-exact CDS",
        }
    return ind, {**check, "status": "verified", "exact_disagreements": 0}


def derive_metrics(
    out_dir: Path,
    multi_segment: dict | None = None,
    independent: pd.DataFrame | None = None,
    independent_note: str | None = None,
) -> tuple[list[dict], pd.DataFrame | None, pd.DataFrame | None]:
    """
    Args:
            out_dir: pairwise-compare output directory
            multi_segment: reference gene_id -> True if any evaluated isoform has > 1 CDS segment
            independent: verified independent reference-gene table (load_independent)
    Returns:
            (metric rows, reference outcomes or None, query outcomes or None)
    """
    out_dir = Path(out_dir)
    s = json.loads((out_dir / "comparison_summary.json").read_text())
    R, Q = s["total_reference_genes"], s["total_consensus_genes"]
    sc, sens, spec = s["sensitivity_cds"], s["sensitivity"], s["specificity"]
    isup = s.get("intron_support", {}).get("cds")
    rows = [
        _metric("reference_coding_genes", R),
        _metric(
            "reference_eligible_transcripts",
            s["filter_info"]["post_filter_transcripts"],
        ),
        _metric("predicted_coding_genes", Q),
        _metric("cds_exact_ref", sc["cds_coordinate_exact_count"], R),
        _metric("cds_exact_query", spec["consensus_cds_coordinate_exact_count"], Q),
        _metric(
            "cds_intron_chain_best_pair",
            sc["cds_intron_chain_recovered"],
            sc["multi_segment_cds_reference_genes"],
        ),
        _metric("gene_span_locus_detection", sens["locus_detected_count"], R),
        _metric(
            "gene_span_only_no_exon_overlap",
            s["locus_detection_exonic"]["span_only_detected_count"],
            R,
        ),
        _metric("missed_reference_genes", sens["missed_count"], R),
        _metric("novel_predictions", spec["novel_consensus_count"], Q),
        _metric("splits", s["split_merge"]["gene_split_count"], R),
        _metric("merges", s["split_merge"]["gene_merge_count"], Q),
        _metric(
            "one_to_one_ref_genes", s["split_merge"]["one_to_one_reference_genes"], R
        ),
        _metric("strand_mismatch_ref", sens["strand_mismatch_count"], R),
        _metric("strand_mismatch_query", spec["strand_mismatch_consensus_count"], Q),
        _metric("exon_coordinate_exact_ref", sens["exon_coordinate_exact_count"], R),
        _metric("structural_exact_match_with_utr", sens["exact_match_count"], R),
        _metric("cds_structural_exact_match", sc["cds_exact_match_count"], R),
    ]
    rows += [
        _metric(f"novel_{c}", s["novel_categories"].get(c, 0), Q)
        for c in NOVEL_CATEGORIES
    ]
    if isup:
        rows += [
            _metric(
                "ref_cds_introns_recovered",
                isup["shared_introns"],
                isup["reference_introns"],
            ),
            _metric(
                "query_cds_introns_supported",
                isup["shared_introns"],
                isup["query_introns"],
            ),
        ]
    else:
        rows += [
            _na(
                "ref_cds_introns_recovered",
                "comparator",
                "intron_support missing from summary",
            ),
            _na(
                "query_cds_introns_supported",
                "comparator",
                "intron_support missing from summary",
            ),
        ]

    det = read_details(out_dir)
    if det is None:
        note = "comparison_details.tsv not available (summary-only result)"
        for mid in (
            "predicted_transcripts",
            "locus_recovery_cds_overlap",
            "detected_without_cds_overlap",
            "partial_cds",
            "partial_same_chain_terminal_diff",
            "partial_ge08_diff_chain",
            "partial_lt08",
            "ref_genes_without_cds_overlap",
            "query_matched_without_cds_overlap",
            "exact_pairs_query_gene_reused",
            "exact_tp_one_to_one",
            "exact_recall_one_to_one",
            "exact_precision_one_to_one",
            "exact_f1_one_to_one",
            "cds_exact_ref_multi",
            "cds_exact_ref_single",
            "cds_intron_chain_any_pair",
        ):
            rows.append(_na(mid, "comparator_per_gene", note))
        return rows, None, None

    ref = reference_outcomes(det)
    qry = query_outcomes(det)
    assert (
        len(ref) == R and len(qry) == Q
    ), "details and summary disagree on gene counts"
    labels = out_dir / "consensus_transcript_labels.tsv"
    n_labels = sum(1 for _ in open(labels)) - 1 if labels.exists() else None
    rows.append(
        _metric("predicted_transcripts", n_labels, source="comparator_per_gene")
        if n_labels is not None
        else _na(
            "predicted_transcripts",
            "comparator_per_gene",
            "consensus_transcript_labels.tsv missing",
        )
    )
    coding = ref.cds_overlap > 0
    detected = ref.classification.isin(DETECTED)
    exact = ref.cds_exact == 1
    partial = coding & ~exact
    pg = "comparator_per_gene"
    rows += [
        _metric("locus_recovery_cds_overlap", int(coding.sum()), R, source=pg),
        _metric(
            "detected_without_cds_overlap",
            int((detected & ~coding).sum()),
            R,
            source=pg,
        ),
        _metric("partial_cds", int(partial.sum()), R, source=pg),
        _metric(
            "partial_same_chain_terminal_diff",
            int((partial & (ref.classification_cds == "Exact_Match")).sum()),
            R,
            source=pg,
        ),
        _metric(
            "partial_ge08_diff_chain",
            int((partial & (ref.classification_cds == "Structural_Mismatch")).sum()),
            R,
            source=pg,
        ),
        _metric(
            "partial_lt08",
            int((partial & (ref.classification_cds == "Partial_Match")).sum()),
            R,
            source=pg,
        ),
        _metric("ref_genes_without_cds_overlap", int(R - coding.sum()), R, source=pg),
        _metric(
            "query_matched_without_cds_overlap",
            int(
                (
                    qry.classification.isin(DETECTED | {"Matched"})
                    & ~(qry.cds_overlap > 0)
                ).sum()
            ),
            Q,
            source=pg,
        ),
    ]
    pairs = ref.loc[exact, ["gene_id", "cds_matched_id"]]
    tp = min(pairs.gene_id.nunique(), pairs.cds_matched_id.nunique())
    recall, precision = (tp / R if R else None), (tp / Q if Q else None)
    f1 = 2 * recall * precision / (recall + precision) if recall and precision else None
    rows += [
        _metric(
            "exact_pairs_query_gene_reused",
            int(pairs.cds_matched_id.duplicated().sum()),
            source=pg,
        ),
        _metric("exact_tp_one_to_one", tp, source=pg),
        _metric("exact_recall_one_to_one", tp, R, source=pg),
        _metric("exact_precision_one_to_one", tp, Q, source=pg),
        _metric(
            "exact_f1_one_to_one",
            tp,
            None,
            value=round(f1, 4) if f1 is not None else None,
            source=pg,
            note="harmonic mean of one-to-one recall (TP/R) and precision (TP/Q)",
        ),
    ]
    if multi_segment is not None:
        is_multi = ref.gene_id.map(multi_segment)
        if is_multi.isna().any():
            note = f"{int(is_multi.isna().sum())} reference genes missing from reference models"
            rows += [
                _na("cds_exact_ref_multi", "reference_models", note),
                _na("cds_exact_ref_single", "reference_models", note),
            ]
        else:
            is_multi = is_multi.astype(bool)
            n_multi = int(is_multi.sum())
            agree = n_multi == sc["multi_segment_cds_reference_genes"]
            note = (
                None
                if agree
                else f"multi-segment count {n_multi} differs from comparator {sc['multi_segment_cds_reference_genes']}"
            )
            rows += [
                _metric(
                    "cds_exact_ref_multi",
                    int((exact & is_multi).sum()),
                    n_multi,
                    source="reference_models",
                    note=note,
                ),
                _metric(
                    "cds_exact_ref_single",
                    int((exact & ~is_multi).sum()),
                    int((~is_multi).sum()),
                    source="reference_models",
                    note=note,
                ),
            ]
            ref["multi_segment"] = is_multi.astype(int).values
    else:
        note = "needs the reference annotation (per-gene CDS segment counts); build the dashboard with reference models"
        rows += [
            _na("cds_exact_ref_multi", "reference_models", note),
            _na("cds_exact_ref_single", "reference_models", note),
        ]
    if independent is not None:
        ind = independent.set_index("gene_id")
        anyp = ind.cds_intron_chain_any_pair == "True"
        n_multi = int((ind.multi_segment == "True").sum())
        rows.append(
            _metric(
                "cds_intron_chain_any_pair",
                int(anyp.sum()),
                n_multi,
                source="independent",
            )
        )
        ref["any_pair_chain"] = (
            ref.gene_id.map(ind.cds_intron_chain_any_pair)
            .replace({"True": 1, "False": 0, "": None})
            .values
        )
    else:
        rows.append(
            _na(
                "cds_intron_chain_any_pair",
                "independent",
                independent_note
                or "independent analysis not run for this result (`annotation-qc benchmark independent`)",
                den=sc["multi_segment_cds_reference_genes"],
            )
        )
    return rows, ref, qry


# Headline table column -> (metric_id, field). Used to verify reproduction of the
# benchmark's headline_metrics.tsv (wide). The long table is the same data melted
# and must not be imported in addition.
HEADLINE_COLUMNS = {
    "reference_coding_genes": ("reference_coding_genes", "numerator"),
    "reference_eligible_transcripts": ("reference_eligible_transcripts", "numerator"),
    "predicted_coding_genes": ("predicted_coding_genes", "numerator"),
    "predicted_transcripts": ("predicted_transcripts", "numerator"),
    "locus_recovery_same_strand_cds_overlap": (
        "locus_recovery_cds_overlap",
        "numerator",
    ),
    "locus_recovery_same_strand_cds_overlap_rate": (
        "locus_recovery_cds_overlap",
        "value",
    ),
    "gene_span_locus_detection": ("gene_span_locus_detection", "numerator"),
    "gene_span_only_no_exon_overlap": ("gene_span_only_no_exon_overlap", "numerator"),
    "detected_without_cds_overlap": ("detected_without_cds_overlap", "numerator"),
    "cds_coordinate_exact_ref_genes": ("cds_exact_ref", "numerator"),
    "cds_coordinate_exact_rate": ("cds_exact_ref", "value"),
    "multi_segment_ref_genes": ("cds_exact_ref_multi", "denominator"),
    "cds_coordinate_exact_multi_segment": ("cds_exact_ref_multi", "numerator"),
    "cds_coordinate_exact_multi_segment_rate": ("cds_exact_ref_multi", "value"),
    "single_segment_ref_genes": ("cds_exact_ref_single", "denominator"),
    "cds_coordinate_exact_single_segment": ("cds_exact_ref_single", "numerator"),
    "cds_coordinate_exact_single_segment_rate": ("cds_exact_ref_single", "value"),
    "partial_cds_recovery": ("partial_cds", "numerator"),
    "partial_cds_recovery_rate": ("partial_cds", "value"),
    "partial_same_chain_terminal_diff": (
        "partial_same_chain_terminal_diff",
        "numerator",
    ),
    "partial_ge0.8_diff_chain": ("partial_ge08_diff_chain", "numerator"),
    "partial_lt0.8": ("partial_lt08", "numerator"),
    "cds_intron_chain_best_pair": ("cds_intron_chain_best_pair", "numerator"),
    "cds_intron_chain_denominator_multi_segment": (
        "cds_intron_chain_best_pair",
        "denominator",
    ),
    "cds_intron_chain_best_pair_rate": ("cds_intron_chain_best_pair", "value"),
    "cds_intron_chain_any_pair_independent": ("cds_intron_chain_any_pair", "numerator"),
    "cds_intron_chain_any_pair_rate": ("cds_intron_chain_any_pair", "value"),
    "ref_unique_cds_introns": ("ref_cds_introns_recovered", "denominator"),
    "ref_unique_cds_introns_recovered": ("ref_cds_introns_recovered", "numerator"),
    "ref_unique_cds_introns_recovered_rate": ("ref_cds_introns_recovered", "value"),
    "query_unique_cds_introns": ("query_cds_introns_supported", "denominator"),
    "query_cds_introns_supported_rate": ("query_cds_introns_supported", "value"),
    "query_cds_coordinate_exact_genes": ("cds_exact_query", "numerator"),
    "query_cds_coordinate_exact_rate": ("cds_exact_query", "value"),
    "exon_coordinate_exact_ref_genes": ("exon_coordinate_exact_ref", "numerator"),
    "structural_exact_match_with_utr": ("structural_exact_match_with_utr", "numerator"),
    "missed_reference_genes": ("missed_reference_genes", "numerator"),
    "missed_rate": ("missed_reference_genes", "value"),
    "ref_genes_without_cds_overlap": ("ref_genes_without_cds_overlap", "numerator"),
    "unmatched_predictions_novel": ("novel_predictions", "numerator"),
    "unmatched_predictions_novel_rate": ("novel_predictions", "value"),
    **{f"novel_{c}": (f"novel_{c}", "numerator") for c in NOVEL_CATEGORIES},
    "query_matched_without_cds_overlap": (
        "query_matched_without_cds_overlap",
        "numerator",
    ),
    "splits_ref_genes": ("splits", "numerator"),
    "merges_query_genes": ("merges", "numerator"),
    "one_to_one_ref_genes": ("one_to_one_ref_genes", "numerator"),
    "strand_mismatch_ref": ("strand_mismatch_ref", "numerator"),
    "strand_mismatch_query": ("strand_mismatch_query", "numerator"),
    "exact_pairs_query_gene_reused": ("exact_pairs_query_gene_reused", "numerator"),
    "exact_tp_one_to_one": ("exact_tp_one_to_one", "numerator"),
    "gene_recall_exact_cds": ("exact_recall_one_to_one", "value"),
    "gene_precision_exact_cds": ("exact_precision_one_to_one", "value"),
    "gene_f1_exact_cds": ("exact_f1_one_to_one", "value"),
}
# The benchmark computed F1 from recall and precision already rounded to 4 dp, so F1
# may differ from the unrounded harmonic mean in the 4th decimal.
HEADLINE_TOLERANCE = {"gene_f1_exact_cds": 1.5e-4}


def compare_headline_value(column: str, expected, derived) -> bool:
    if expected is None or (isinstance(expected, float) and math.isnan(expected)):
        return derived is None
    if derived is None:
        return False
    return abs(float(expected) - float(derived)) <= HEADLINE_TOLERANCE.get(
        column, 5e-9 if float(expected).is_integer() else 5e-5
    )

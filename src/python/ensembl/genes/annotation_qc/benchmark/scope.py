"""
Headline and structural metrics restricted to a sequence scope.

A scope is a configured set of sequence names (``references[].analysis_scopes``),
for example the assembled autosomes and X. It is applied to BOTH sides before
any value is calculated: reference genes and predicted genes on other sequences
are removed, and so are their one-to-one pairs and introns. The whole compared
assembly stays the primary analysis; its values are the comparator's own and are
never replaced by the values computed here.

Derivation, from one comparator result:

    R, Q, locus recovery, per-gene exact CDS, splits, merges, CDS intron chain
        per-gene rows of comparison_details.tsv on the scope's sequences
    one-to-one TP, precision, recall, F1
        rows of cds_exact_one_to_one_pairs.tsv on the scope's sequences, over the
        scoped R and Q. Exact-CDS edges only join genes on one sequence, so the
        maximum matching decomposes by sequence and its restriction to the scope
        is a maximum matching of the scoped genes.
    unique CDS introns
        benchmark/splice.py on the evaluated models, deduplicated as the
        comparator does

``check_whole_assembly`` repeats the derivation without a scope and compares it
with the comparator's summary values; any difference marks the scoped rows
inconsistent instead of showing them as verified values.
"""

from __future__ import annotations

import pandas as pd

from ensembl.genes.annotation_qc.benchmark.derive import _metric, _na

SOURCE = "scope_filter"
SCOPE_METRIC_IDS = (
    "reference_coding_genes",
    "predicted_coding_genes",
    "locus_recovery_cds_overlap",
    "cds_exact_ref",
    "cds_exact_query",
    "exact_tp_one_to_one",
    "exact_recall_one_to_one",
    "exact_precision_one_to_one",
    "exact_f1_one_to_one",
    "cds_intron_chain_best_pair",
    "splits",
    "merges",
    "ref_cds_introns_recovered",
    "query_cds_introns_supported",
)
WHOLE_ASSEMBLY = "whole_assembly"


def _within(frame: pd.DataFrame, sequences: set | None) -> pd.DataFrame:
    return frame if sequences is None else frame[frame.chrom.isin(sequences)]


def scoped_metrics(
    ref_out: pd.DataFrame,
    qry_out: pd.DataFrame,
    pairs: pd.DataFrame | None,
    sequences: set | None,
    introns: dict | None = None,
    source: str = SOURCE,
) -> list[dict]:
    """
    Metric rows (derive.py layout) for one result restricted to ``sequences``.

    Args:
            ref_out: reference outcomes (gene_id, chrom, cds_overlap, cds_exact,
                    cds_intron_chain_match, counterpart_count, multi_segment)
            qry_out: query outcomes (gene_id, chrom, cds_exact, counterpart_count)
            pairs: the comparator's one-to-one exact-CDS pairs (chrom column), or None
            sequences: scope sequence names; None = no restriction
            introns: splice.analyse() totals for this scope, or None
    """
    ref = _within(ref_out, sequences)
    qry = _within(qry_out, sequences)
    R, Q = len(ref), len(qry)
    rows = [
        _metric("reference_coding_genes", R, source=source),
        _metric("predicted_coding_genes", Q, source=source),
        _metric("locus_recovery_cds_overlap", int((ref.cds_overlap > 0).sum()), R, source=source),
        _metric("cds_exact_ref", int((ref.cds_exact == 1).sum()), R, source=source),
        _metric("cds_exact_query", int((qry.cds_exact == 1).sum()), Q, source=source),
        _metric("splits", int((ref.counterpart_count >= 2).sum()), R, source=source),
        _metric("merges", int((qry.counterpart_count >= 2).sum()), Q, source=source),
    ]
    if pairs is None:
        note = "cds_exact_one_to_one_pairs.tsv not available for this result"
        rows += [
            _na(mid, source, note)
            for mid in (
                "exact_tp_one_to_one",
                "exact_recall_one_to_one",
                "exact_precision_one_to_one",
                "exact_f1_one_to_one",
            )
        ]
    else:
        tp = len(_within(pairs, sequences))

        def rate(mid, num, den, label):
            if not den:
                return _na(mid, source, f"{label} undefined: denominator is 0", den=den, status="undefined")
            return _metric(mid, num, den, source=source)

        rows += [
            _metric("exact_tp_one_to_one", tp, source=source),
            rate("exact_recall_one_to_one", tp, R, "recall"),
            rate("exact_precision_one_to_one", tp, Q, "precision"),
            rate("exact_f1_one_to_one", 2 * tp, R + Q, "F1"),
        ]
    multi = ref.get("multi_segment")
    if multi is None or multi.isna().any():
        rows.append(
            _na(
                "cds_intron_chain_best_pair",
                source,
                "needs per-gene CDS segment counts from the reference models",
            )
        )
    else:
        chain = int((ref.cds_intron_chain_match == "True").sum())
        rows.append(
            _metric("cds_intron_chain_best_pair", chain, int((multi == 1).sum()), source=source)
        )
    if introns is None:
        note = "unique-intron analysis not available for this result"
        rows += [
            _na("ref_cds_introns_recovered", source, note),
            _na("query_cds_introns_supported", source, note),
        ]
    else:
        rows += [
            _metric(
                "ref_cds_introns_recovered",
                introns["shared_introns"],
                introns["reference_introns"],
                source=source,
            ),
            _metric(
                "query_cds_introns_supported",
                introns["shared_introns"],
                introns["query_introns"],
                source=source,
            ),
        ]
    return rows


def check_whole_assembly(derived: list[dict], comparator: list[dict]) -> dict:
    """Compare an unscoped derivation with the comparator's metric rows."""
    by_id = {m["metric_id"]: m for m in comparator}
    differences, compared = [], 0
    for row in derived:
        ref = by_id.get(row["metric_id"])
        if ref is None or ref["status"] not in ("ok", "inconsistent") or row["status"] != "ok":
            continue
        compared += 1
        for key in ("numerator", "denominator"):
            a, b = row[key], ref[key]
            if a is None or b is None:
                continue
            if float(a) != float(b):
                differences.append(f"{row['metric_id']}.{key}: derived {a}, comparator {b}")
        if row["value"] is not None and ref["value"] is not None and abs(row["value"] - ref["value"]) > 5e-5:
            differences.append(f"{row['metric_id']}.value: derived {row['value']}, comparator {ref['value']}")
    return {
        "status": "verified" if not differences else "disagrees",
        "metrics_compared": compared,
        "differences": differences,
    }


def mark_inconsistent(rows: list[dict], check: dict) -> list[dict]:
    if check["status"] == "verified":
        return rows
    note = "whole-assembly reproduction failed: " + "; ".join(check["differences"][:5])
    return [
        {**r, "status": "inconsistent", "note": note} if r["status"] == "ok" else r
        for r in rows
    ]


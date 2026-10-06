"""
Headline metrics and per-transcript labels for a pairwise comparison.

Key names follow gmb-compare so existing consumers keep working; "consensus"
in a key means the query annotation. summary["rate_denominators"] records the
denominator of every rate.

Reference-based (sensitivity-like) rates divide by the number of reference genes
after filtering (R), unless stated:

    locus detection     reference genes classified Exact_Match, Partial_Match or
                        Structural_Mismatch (a same-strand gene-span overlap) / R
    CDS-overlap locus   reference genes whose best same-strand CDS pair shares at
    recovery            least one CDS base (cds_overlap > 0) / R. A pair shares a
                        base exactly when its reciprocal CDS overlap is > 0, and
                        the best pair has the largest overlap, so this equals
                        "some same-strand locus partner shares CDS bases". Rate
                        is None when R = 0.
    exonic detection    the subset of those with some exonic overlap / R; the
                        difference ("span-only") is detection that rests on gene
                        spans alone
    Exact_Match         structural: identical intron chain and >= 0.8 reciprocal
                        overlap, terminal coordinates may differ / R
    coordinate exact    identical exon (or CDS) intervals, terminal coordinates
                        included / R
    CDS exact/any       classification_cds Exact_Match / any of the three match
                        classes / R
    intron chain        only reference genes with an intron count. Recovered =
                        detected gene whose best pair has an identical, non-empty
                        intron chain. exon_intron_chain_rate divides by detected
                        genes whose compared reference transcript has an intron;
                        exon_intron_chain_sensitivity divides by all reference
                        genes with a multi-exon transcript. The CDS equivalents
                        use CDS segments.
    split               reference genes with >= 2 query counterparts (count)

Query-based (precision-like) rates divide by the number of query genes (Q):

    matched             query genes with a same-strand reference partner
                        (Exact/Partial/Structural/Matched) / Q
    novel               query genes with no reference partner on either strand
                        and no opposite-strand feature overlap / Q, broken down
                        by novel_category
    exact, coordinate   query genes whose best pair is Exact_Match / coordinate
    exact               exact / Q
    merge               query genes with >= 2 reference counterparts (count)

One-to-one coordinate-exact CDS precision/recall/F1 are computed in
metrics/pairwise/matching.py.

summarise_intron_support compares unique introns: reference introns found in the
query / reference introns, and query introns found in the reference / query
introns.
"""

from collections import Counter, defaultdict

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    BASIS_CDS,
    BASIS_EXON,
    BASIS_NO_FEATURE,
    EXACT,
    MATCHED,
    MISSED,
    NO_CDS,
    NOVEL,
    NOVEL_CATEGORIES,
    PARTIAL,
    SPLIT_MERGE_MIN_SHARED_FRACTION,
    STRAND_MISMATCH,
    STRUCTURAL,
    GeneModels,
    intron_chain,
)

DETECTED = (EXACT, PARTIAL, STRUCTURAL)

RATE_DENOMINATORS = {
    "sensitivity.*_rate": "total_reference_genes",
    "sensitivity_cds.*_rate (except cds_intron_chain_*)": "total_reference_genes",
    "locus_detection_exonic.locus_detection_exonic_rate": "total_reference_genes",
    "intron_chain.exon_intron_chain_rate": (
        "intron_chain.exon_intron_chain_evaluated (detected reference genes whose "
        "compared transcript has an intron)"
    ),
    "intron_chain.exon_intron_chain_sensitivity": (
        "intron_chain.multi_exon_reference_genes"
    ),
    "sensitivity_cds.cds_intron_chain_rate": (
        "sensitivity_cds.cds_intron_chain_evaluated (CDS-detected reference genes "
        "whose compared CDS has an intron)"
    ),
    "sensitivity_cds.cds_intron_chain_sensitivity": (
        "sensitivity_cds.multi_segment_cds_reference_genes"
    ),
    "specificity.*_rate": "total_consensus_genes (query genes)",
    "cds_overlap_locus_recovery.rate": "total_reference_genes",
    "cds_exact_one_to_one.recall": "total_reference_genes",
    "cds_exact_one_to_one.precision": "total_consensus_genes (query genes)",
    "cds_exact_one_to_one.f1": "total_reference_genes + total_consensus_genes",
    "intron_support.*.reference_introns_recovered_rate": "reference_introns",
    "intron_support.*.query_introns_supported_rate": "query_introns",
}


def _rate(count: int, total: int) -> float:
    return round(count / total, 4) if total else 0.0


def _is_true(values: pd.Series) -> pd.Series:
    return values.astype("boolean").fillna(False).astype(bool)


def _applicable(values: pd.Series) -> pd.Series:
    return values.astype("boolean").notna()


def _stratified_counts(ref: pd.DataFrame) -> dict:
    strata = {
        k: Counter() for k in ("single_exon", "multi_exon", "cds_containing", "non_cds")
    }
    for row in ref.itertuples(index=False):
        # A gene is multi-exon when any of its transcripts has more than one exon.
        single = row.max_exon_count <= 1
        strata["single_exon" if single else "multi_exon"][row.classification] += 1
        has_cds = row.classification_cds not in (NO_CDS, MISSED, "")
        cds_key = "cds_containing" if has_cds or row.cds_overlap > 0 else "non_cds"
        strata[cds_key][row.classification] += 1
    return {k: dict(v) for k, v in strata.items()}


def _intron_chain_block(
    detected: pd.Series, chain: pd.Series, multi: pd.Series, prefix: str
) -> dict:
    recovered = int((detected & _is_true(chain)).sum())
    evaluated = int((detected & _applicable(chain)).sum())
    return {
        f"{prefix}_recovered": recovered,
        # gmb-compare name for the rate denominator; now multi-exon genes only.
        f"{prefix}_matched": evaluated,
        f"{prefix}_evaluated": evaluated,
        f"{prefix}_not_applicable_single_exon": int(
            (detected & ~_applicable(chain)).sum()
        ),
        f"{prefix}_rate": _rate(recovered, evaluated),
        f"{prefix}_sensitivity": _rate(recovered, int(multi.sum())),
    }


def summarise_comparison(ref: pd.DataFrame, query: pd.DataFrame) -> dict:
    """
    Compute summary metrics from classified reference and query genes.
    Args:
            ref: Reference results from classify_loci
            query: Query results from classify_loci
    Returns:
            dict with the gmb-compare summary keys plus the refinements listed
            in the module docstring
    """
    ref_counts = Counter(ref["classification"])
    query_counts = Counter(query["classification"])
    cds_counts = Counter(ref["classification_cds"])
    total_ref, total_query = len(ref), len(query)

    per_chromosome: dict[str, Counter] = defaultdict(Counter)
    for chrom, cls in zip(ref["chrom"], ref["classification"]):
        per_chromosome[chrom][cls] += 1

    detected = ref["classification"].isin(DETECTED)
    cds_matched = ref["classification_cds"].isin(DETECTED)
    exonic = detected & (ref["exon_overlap"] > 0)
    cds_recovered = ref["cds_overlap"].astype(float) > 0
    exon_coord = _is_true(ref["exon_coordinate_exact"])
    cds_coord = _is_true(ref["cds_coordinate_exact"])
    multi_exon = ref["max_exon_count"] > 1
    multi_cds = ref["max_cds_count"] > 1
    basis = Counter(ref["strand_mismatch_basis"])

    exon_chain = _intron_chain_block(
        detected, ref["intron_chain_match"], multi_exon, "exon_intron_chain"
    )
    exon_chain["multi_exon_reference_genes"] = int(multi_exon.sum())
    cds_chain = _intron_chain_block(
        cds_matched, ref["cds_intron_chain_match"], multi_cds, "cds_intron_chain"
    )
    cds_chain["multi_segment_cds_reference_genes"] = int(multi_cds.sum())

    # Query genes with a same-strand reference partner; each is counted once.
    query_matched = query["classification"].isin(DETECTED + (MATCHED,))
    query_novel = query["classification"] == NOVEL
    novel_breakdown = Counter(query.loc[query_novel, "novel_category"])

    summary = {
        "total_reference_genes": total_ref,
        "total_consensus_genes": total_query,
        "reference_classification": dict(ref_counts),
        "reference_classification_cds": dict(cds_counts),
        "consensus_classification": dict(query_counts),
        "per_chromosome": {c: dict(v) for c, v in sorted(per_chromosome.items())},
        "sensitivity": {
            "exact_match_count": ref_counts.get(EXACT, 0),
            "partial_match_count": ref_counts.get(PARTIAL, 0),
            "any_match_count": sum(ref_counts.get(c, 0) for c in DETECTED),
            "missed_count": ref_counts.get(MISSED, 0),
            "strand_mismatch_count": ref_counts.get(STRAND_MISMATCH, 0),
            "strand_mismatch_cds_count": basis.get(BASIS_CDS, 0),
            "strand_mismatch_exon_count": basis.get(BASIS_EXON, 0),
            # Missed genes that gmb-compare called Strand_Mismatch from gene-span contact.
            "opposite_strand_no_feature_overlap_count": basis.get(BASIS_NO_FEATURE, 0),
            "locus_detected_count": int(detected.sum()),
            "locus_detection_rate": _rate(int(detected.sum()), total_ref),
            "exon_coordinate_exact_count": int(exon_coord.sum()),
            "exon_coordinate_exact_rate": _rate(int(exon_coord.sum()), total_ref),
        },
        "sensitivity_cds": {
            "cds_exact_match_count": cds_counts.get(EXACT, 0),
            "cds_any_match_count": int(cds_matched.sum()),
            "cds_missed_count": cds_counts.get(MISSED, 0),
            "cds_exact_but_exon_differs": int(
                (
                    (ref["classification_cds"] == EXACT)
                    & (ref["classification"] != EXACT)
                ).sum()
            ),
            "cds_coordinate_exact_count": int(cds_coord.sum()),
            "cds_coordinate_exact_rate": _rate(int(cds_coord.sum()), total_ref),
            "cds_exact_match_not_coordinate_exact": int(
                ((ref["classification_cds"] == EXACT) & ~cds_coord).sum()
            ),
            # Reference genes whose best CDS pair (reference isoform, query
            # transcript) is not their best exon pair.
            "cds_pair_differs_from_exon_pair": int(
                (
                    (ref["best_cds_match_transcript_id"] != "")
                    & (
                        (
                            ref["best_cds_match_transcript_id"]
                            != ref["best_match_transcript_id"]
                        )
                        | (
                            ref["best_cds_match_own_transcript_id"]
                            != ref["best_match_own_transcript_id"]
                        )
                    )
                ).sum()
            ),
            **cds_chain,
        },
        "intron_chain": exon_chain,
        "specificity": {
            "novel_consensus_count": query_counts.get(NOVEL, 0),
            "matched_consensus_count": int(query_matched.sum()),
            "strand_mismatch_consensus_count": query_counts.get(STRAND_MISMATCH, 0),
            "novel_consensus_rate": _rate(query_counts.get(NOVEL, 0), total_query),
            "matched_consensus_rate": _rate(int(query_matched.sum()), total_query),
            "consensus_exact_match_count": query_counts.get(EXACT, 0),
            "consensus_exact_match_rate": _rate(
                query_counts.get(EXACT, 0), total_query
            ),
            "consensus_cds_exact_match_count": int(
                (query["classification_cds"] == EXACT).sum()
            ),
            "consensus_cds_exact_match_rate": _rate(
                int((query["classification_cds"] == EXACT).sum()), total_query
            ),
            "consensus_cds_coordinate_exact_count": int(
                _is_true(query["cds_coordinate_exact"]).sum()
            ),
            "consensus_cds_coordinate_exact_rate": _rate(
                int(_is_true(query["cds_coordinate_exact"]).sum()), total_query
            ),
            "consensus_exon_coordinate_exact_count": int(
                _is_true(query["exon_coordinate_exact"]).sum()
            ),
            "consensus_exon_coordinate_exact_rate": _rate(
                int(_is_true(query["exon_coordinate_exact"]).sum()), total_query
            ),
        },
        "novel_categories": {c: novel_breakdown.get(c, 0) for c in NOVEL_CATEGORIES},
        "split_merge": {
            "gene_split_count": int((ref["counterpart_count"] >= 2).sum()),
            "gene_split_query_genes": int(
                len(
                    {
                        gid
                        for ids in ref.loc[
                            ref["counterpart_count"] >= 2, "counterpart_ids"
                        ]
                        for gid in ids.split(",")
                    }
                )
            ),
            "gene_merge_count": int((query["counterpart_count"] >= 2).sum()),
            "gene_merge_reference_genes": int(
                len(
                    {
                        gid
                        for ids in query.loc[
                            query["counterpart_count"] >= 2, "counterpart_ids"
                        ]
                        for gid in ids.split(",")
                    }
                )
            ),
            "one_to_one_reference_genes": int((ref["counterpart_count"] == 1).sum()),
            "min_shared_fraction": SPLIT_MERGE_MIN_SHARED_FRACTION,
        },
        "stratified": _stratified_counts(ref),
        "locus_detection_exonic": {
            "locus_detected_exonic_count": int(exonic.sum()),
            "locus_detection_exonic_rate": _rate(int(exonic.sum()), total_ref),
            "span_only_detected_count": int((detected & ~exonic).sum()),
        },
        "cds_overlap_locus_recovery": {
            "recovered_count": int(cds_recovered.sum()),
            "reference_genes": total_ref,
            "rate": (
                round(int(cds_recovered.sum()) / total_ref, 4) if total_ref else None
            ),
            "gene_span_detected_without_cds_overlap": int(
                (detected & ~cds_recovered).sum()
            ),
            "definition": (
                "reference genes with a same-strand query locus partner whose "
                "transcripts share at least one CDS base (cds_overlap > 0, "
                "unrounded) / total_reference_genes"
            ),
        },
        "rate_denominators": RATE_DENOMINATORS,
    }

    if total_ref:
        sens, sens_cds = summary["sensitivity"], summary["sensitivity_cds"]
        sens["exact_match_rate"] = round(ref_counts.get(EXACT, 0) / total_ref, 4)
        sens["any_match_rate"] = round(sens["any_match_count"] / total_ref, 4)
        sens["missed_rate"] = round(ref_counts.get(MISSED, 0) / total_ref, 4)
        if sens_cds["cds_any_match_count"] + sens_cds["cds_missed_count"] > 0:
            sens_cds["cds_exact_match_rate"] = round(
                sens_cds["cds_exact_match_count"] / total_ref, 4
            )
            sens_cds["cds_any_match_rate"] = round(
                sens_cds["cds_any_match_count"] / total_ref, 4
            )
    return summary


def _unique_introns(models: GeneModels, kind: str) -> set:
    source = models.exons if kind == "exon" else models.cds
    introns = set()
    for chrom, strand, tids in zip(
        models.genes["Chromosome"],
        models.genes["Strand"],
        models.genes["transcript_ids"],
    ):
        for tid in tids:
            for start, end in intron_chain(source.get(tid, ())):
                introns.add((chrom, strand, start, end))
    return introns


def summarise_intron_support(reference: GeneModels, query: GeneModels) -> dict:
    """
    Compare unique introns (seqname, strand, start, end) between annotations.

    Exon introns include introns inside UTRs; CDS introns are those between CDS
    segments. Single-exon transcripts contribute nothing.
    Args:
            reference: GeneModels for the filtered reference
            query: GeneModels for the query
    Returns:
            {"exon": {...}, "cds": {...}} with counts, the reference-based recovery
            rate and the query-based support rate
    """
    result = {}
    for kind in ("exon", "cds"):
        ref_introns = _unique_introns(reference, kind)
        query_introns = _unique_introns(query, kind)
        shared = len(ref_introns & query_introns)
        result[kind] = {
            "reference_introns": len(ref_introns),
            "query_introns": len(query_introns),
            "shared_introns": shared,
            "reference_introns_recovered_rate": _rate(shared, len(ref_introns)),
            "query_introns_supported_rate": _rate(shared, len(query_introns)),
        }
    return result


def summarise_span_excess(excess: pd.DataFrame) -> dict:
    """
    Summarise gene_span_excess output.
    Args:
            excess: DataFrame with gene_id and excess_bp
    Returns:
            dict with gene counts and the largest excess
    """
    wider = excess[excess["excess_bp"] > 0]
    return {
        "genes_checked": len(excess),
        "genes_wider_than_transcripts": len(wider),
        "max_excess_bp": int(wider["excess_bp"].max()) if len(wider) else 0,
        "examples": wider.nlargest(5, "excess_bp")["gene_id"].tolist(),
    }


def query_transcript_labels(query: pd.DataFrame) -> pd.DataFrame:
    """
    One row per query transcript with its gene's comparison label.

    Exact_Match, Partial_Match, Structural_Mismatch and Matched become Matched;
    Strand_Mismatch is kept; everything else is Novel, with novel_category
    saying what reference feature, if any, the gene overlaps.
    Args:
            query: Query results from classify_loci
    Returns:
            DataFrame with transcript_id, gene_id, classification,
            best_ref_gene_id, best_ref_transcript_id, best_overlap,
            best_cds_overlap, best_cds_ref_gene_id, best_cds_ref_transcript_id,
            novel_category. best_overlap comes from the best exon pair and
            best_cds_overlap from the best CDS pair.
    """

    def label(cls: str) -> str:
        if cls in DETECTED or cls == MATCHED:
            return MATCHED
        return STRAND_MISMATCH if cls == STRAND_MISMATCH else NOVEL

    records = [
        (
            tid,
            row.gene_id,
            label(row.classification),
            row.matched_id,
            row.best_match_transcript_id,
            row.exon_overlap,
            row.cds_overlap,
            row.cds_matched_id,
            row.best_cds_match_transcript_id,
            row.novel_category,
        )
        for row in query.itertuples(index=False)
        for tid in row.transcript_ids
    ]
    return pd.DataFrame(
        records,
        columns=[
            "transcript_id",
            "gene_id",
            "classification",
            "best_ref_gene_id",
            "best_ref_transcript_id",
            "best_overlap",
            "best_cds_overlap",
            "best_cds_ref_gene_id",
            "best_cds_ref_transcript_id",
            "novel_category",
        ],
    )

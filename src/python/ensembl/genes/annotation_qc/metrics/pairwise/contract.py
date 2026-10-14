"""
Output contract of annotation-qc pairwise-compare and its canonical metric table.

CONTRACT_ID / CONTRACT_VERSION identify the set of files a comparison writes and
their columns (documented in the package README, "Output contract"). Minor
versions only add files, columns or metric ids; a changed definition, denominator
or column meaning needs a new major version. DEFINITION_VERSION changes whenever
the scientific definition of any metric changes.

metrics.tsv has one row per (scope, metric). Every value is a count or an
unrounded ratio of two counts written in the same row, with a status:

    ok              value computed
    undefined       the denominator is 0 (reason zero_denominator); no value
    not_applicable  the metric does not apply to this comparison (reason says why)

A metric that is missing is never written as 0. Rows are computed from the
classified genes (classify.classify_loci), the one-to-one exact-CDS pairs
(matching.py), the unique intron sets (summary.unique_introns) and the
stop-codon audit, restricted per scope to genes and introns on the scope's
sequences. Locus pairing, exact-CDS edges, split/merge counterparts and intron
identities never join two sequences, so a scope's values equal those of a
comparison run on those sequences alone.

Populations used as denominators (denominator_id):

    R              evaluated reference genes in the scope (after evaluation mode,
                   biotype filters, transcript selection, genome-mismatch
                   exclusion and regions)
    R_cds          members of R with at least one evaluated CDS transcript
    R_multi_cds    members of R with an evaluated transcript of > 1 CDS segment
    R_single_cds   R_cds minus R_multi_cds
    R_multi_exon   members of R with an evaluated transcript of > 1 exon
    Q, Q_cds       the same for evaluated query genes
    R_plus_Q       R + Q (one-to-one F1 = 2 TP / (R + Q))
    *_introns      unique (sequence, strand, start, end) introns
    query_transcripts_with_cds   query transcripts considered by --add-stop-codon
"""

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise import splice
from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    EXACT,
    MATCHED,
    MISSED,
    NOVEL,
    NOVEL_CATEGORIES,
    PARTIAL,
    STRAND_MISMATCH,
    STRUCTURAL,
    GeneModels,
)
from ensembl.genes.annotation_qc.metrics.pairwise.matching import exact_cds_edges
from ensembl.genes.annotation_qc.metrics.pairwise.summary import (
    CDS_OUTCOMES,
    reference_cds_outcomes,
    unique_introns,
)

CONTRACT_ID = "annotation-qc.pairwise-compare"
CONTRACT_VERSION = "1.0"
DEFINITION_VERSION = "2026-10"
ALL_SCOPE = "all"

STATUS_OK = "ok"
STATUS_UNDEFINED = "undefined"
STATUS_NOT_APPLICABLE = "not_applicable"
REASONS = {
    "zero_denominator": "the denominator is 0",
    "transformation_not_applied": "no coordinate transformation was requested",
    "selection_not_canonical": "the transcript selection is not canonical",
    "requires_one_reference_transcript_per_gene": (
        "needs a reference selection with one transcript per gene "
        "(longest_cds or canonical)"
    ),
}

METRIC_COLUMNS = [
    "scope_id",
    "metric_id",
    "numerator",
    "denominator",
    "value",
    "status",
    "reason",
    "numerator_id",
    "denominator_id",
    "unit",
    "side",
]

DETECTED = (EXACT, PARTIAL, STRUCTURAL)

# metric_id -> (side, unit, numerator description, denominator_id or "")
METRICS: dict[str, tuple[str, str, str, str]] = {
    "reference_genes": ("reference", "genes", "R", ""),
    "reference_genes_with_cds": ("reference", "genes", "R_cds", ""),
    "reference_multi_cds_genes": ("reference", "genes", "R_multi_cds", ""),
    "reference_multi_exon_genes": ("reference", "genes", "R_multi_exon", ""),
    "query_genes": ("query", "genes", "Q", ""),
    "query_genes_with_cds": ("query", "genes", "Q_cds", ""),
    "reference_transcripts_evaluated": ("reference", "transcripts", "count", ""),
    "query_transcripts_evaluated": ("query", "transcripts", "count", ""),
    "reference_canonical_fallback_genes": (
        "reference",
        "genes",
        "reference genes without an Ensembl_canonical transcript, scored on the "
        "longest-CDS fallback",
        "R",
    ),
    "locus_recovery_cds_overlap": (
        "reference",
        "genes",
        "reference genes sharing >= 1 CDS base with an evaluated transcript of a "
        "same-strand gene-span partner (cds_locus_recovered)",
        "R",
    ),
    "cds_coordinate_exact_reference": (
        "reference",
        "genes",
        "reference genes whose best CDS pair has identical CDS intervals",
        "R",
    ),
    "cds_coordinate_exact_query": (
        "query",
        "genes",
        "query genes whose best CDS pair has identical CDS intervals",
        "Q",
    ),
    "exact_cds_one_to_one_tp": (
        "both",
        "gene_pairs",
        "TP: pairs in the maximum one-to-one matching of identical-CDS genes",
        "",
    ),
    "exact_cds_one_to_one_recall": ("reference", "rate", "TP", "R"),
    "exact_cds_one_to_one_precision": ("query", "rate", "TP", "Q"),
    "exact_cds_one_to_one_f1": ("both", "score", "2 TP", "R_plus_Q"),
    "exact_cds_gene_pairs": (
        "both",
        "gene_pairs",
        "(reference gene, query gene) pairs sharing an identical CDS",
        "",
    ),
    "exact_cds_reference_with_candidate_unmatched": (
        "reference",
        "genes",
        "reference genes with an identical-CDS partner (any transcript pair) but "
        "not in the matching",
        "R",
    ),
    "exact_cds_query_with_candidate_unmatched": (
        "query",
        "genes",
        "query genes with an identical-CDS partner (any transcript pair) but not "
        "in the matching",
        "Q",
    ),
    "cds_intron_chain_best_pair": (
        "reference",
        "genes",
        "R_multi_cds genes whose best CDS pair (detected CDS class) has an "
        "identical CDS intron chain",
        "R_multi_cds",
    ),
    "cds_intron_chain_any_pair": (
        "reference",
        "genes",
        "R_multi_cds genes with an identical CDS intron chain between any "
        "evaluated transcript and any transcript of a same-strand partner",
        "R_multi_cds",
    ),
    "exon_intron_chain_best_pair": (
        "reference",
        "genes",
        "R_multi_exon genes whose best exon pair (detected class) has an "
        "identical intron chain",
        "R_multi_exon",
    ),
    "exon_intron_chain_any_pair": (
        "reference",
        "genes",
        "R_multi_exon genes with an identical intron chain in any transcript pair",
        "R_multi_exon",
    ),
    "ref_unique_cds_introns_recovered": (
        "reference",
        "unique_introns",
        "unique reference CDS introns also in the query",
        "ref_unique_cds_introns",
    ),
    "query_unique_cds_introns_supported": (
        "query",
        "unique_introns",
        "unique query CDS introns also in the reference",
        "query_unique_cds_introns",
    ),
    "ref_unique_exon_introns_recovered": (
        "reference",
        "unique_introns",
        "unique reference exon introns also in the query",
        "ref_unique_exon_introns",
    ),
    "query_unique_exon_introns_supported": (
        "query",
        "unique_introns",
        "unique query exon introns also in the reference",
        "query_unique_exon_introns",
    ),
    **{
        f"cds_outcome_{key}": ("reference", "genes", text, "R")
        for key, text in CDS_OUTCOMES.items()
    },
    "partial_cds": (
        "reference",
        "genes",
        "reference genes sharing CDS bases whose CDS is not coordinate-exact",
        "R",
    ),
    "cds_coordinate_exact_reference_multi_cds": (
        "reference",
        "genes",
        "coordinate-exact R_multi_cds genes",
        "R_multi_cds",
    ),
    "cds_coordinate_exact_reference_single_cds": (
        "reference",
        "genes",
        "coordinate-exact R_single_cds genes",
        "R_single_cds",
    ),
    "cds_structural_exact_match": (
        "reference",
        "genes",
        "CDS class Exact_Match (identical CDS intron chain, >= 0.8 reciprocal "
        "overlap; ends may differ)",
        "R",
    ),
    "exon_coordinate_exact_reference": (
        "reference",
        "genes",
        "reference genes whose best exon pair has identical exon intervals",
        "R",
    ),
    "exon_structural_exact_match": (
        "reference",
        "genes",
        "exon class Exact_Match (structural, includes UTRs)",
        "R",
    ),
    "locus_detection_gene_span": (
        "reference",
        "genes",
        "reference genes with a same-strand gene-span partner",
        "R",
    ),
    "locus_detection_exonic": (
        "reference",
        "genes",
        "span-detected reference genes with exonic overlap",
        "R",
    ),
    "span_only_detection": (
        "reference",
        "genes",
        "span-detected reference genes without exonic overlap",
        "R",
    ),
    "detected_without_cds_overlap": (
        "reference",
        "genes",
        "span-detected reference genes sharing no CDS base",
        "R",
    ),
    "ref_genes_without_cds_overlap": (
        "reference",
        "genes",
        "R minus locus_recovery_cds_overlap",
        "R",
    ),
    "missed_reference_genes": ("reference", "genes", "class Missed", "R"),
    "strand_mismatch_reference": ("reference", "genes", "class Strand_Mismatch", "R"),
    "query_matched": (
        "query",
        "genes",
        "query genes with a same-strand reference partner",
        "Q",
    ),
    "query_matched_without_cds_overlap": (
        "query",
        "genes",
        "matched query genes sharing no CDS base with a partner",
        "Q",
    ),
    "query_novel": ("query", "genes", "class Novel", "Q"),
    **{
        f"query_novel_{category}": ("query", "genes", f"Novel, {category}", "Q")
        for category in NOVEL_CATEGORIES
    },
    "strand_mismatch_query": ("query", "genes", "class Strand_Mismatch", "Q"),
    "split_reference_genes": (
        "reference",
        "genes",
        "reference genes with >= 2 query counterparts",
        "R",
    ),
    "split_query_genes_involved": (
        "query",
        "genes",
        "query genes that are counterparts of a split reference gene",
        "Q",
    ),
    "merge_query_genes": (
        "query",
        "genes",
        "query genes with >= 2 reference counterparts",
        "Q",
    ),
    "merge_reference_genes_involved": (
        "reference",
        "genes",
        "reference genes that are counterparts of a merging query gene",
        "R",
    ),
    "one_counterpart_reference_genes": (
        "reference",
        "genes",
        "reference genes with exactly one query counterpart",
        "R",
    ),
    "stop_codon_transcripts_with_cds": (
        "query",
        "transcripts",
        "query_transcripts_with_cds",
        "",
    ),
    "stop_codon_eligible": (
        "query",
        "transcripts",
        "frame established from the sequence",
        "query_transcripts_with_cds",
    ),
    "stop_codon_extended": (
        "query",
        "transcripts",
        "terminal CDS extended by an adjacent stop codon",
        "query_transcripts_with_cds",
    ),
    "stop_codon_unchanged": (
        "query",
        "transcripts",
        "unchanged (stop already in CDS, or no adjacent stop)",
        "query_transcripts_with_cds",
    ),
    "stop_codon_skipped": (
        "query",
        "transcripts",
        "skipped (frame not established or extension not possible)",
        "query_transcripts_with_cds",
    ),
}


def metric_definitions() -> dict[str, dict]:
    """The registry as JSON-ready records (written to the manifest)."""
    return {
        mid: {
            "side": side,
            "unit": unit,
            "numerator": numerator,
            "denominator_id": denominator or None,
        }
        for mid, (side, unit, numerator, denominator) in METRICS.items()
    }


def _true(values: pd.Series) -> pd.Series:
    return values.astype("boolean").fillna(False).astype(bool)


def _row(scope: str, mid: str, numerator, denominator=None, **status) -> dict:
    side, unit, _, denominator_id = METRICS[mid]
    numerator_id = "TP" if mid.startswith("exact_cds_one_to_one_") else mid
    if mid == "exact_cds_one_to_one_f1":
        numerator_id = "2TP"
    row = {
        "scope_id": scope,
        "metric_id": mid,
        "numerator": numerator,
        "denominator": denominator,
        "value": None,
        "status": STATUS_OK,
        "reason": "",
        "numerator_id": numerator_id,
        "denominator_id": denominator_id,
        "unit": unit,
        "side": side,
    }
    if status:
        row.update(numerator=None, denominator=None, **status)
    elif not denominator_id:
        row["value"] = numerator
    elif denominator:
        row["value"] = numerator / denominator
    else:
        row.update(status=STATUS_UNDEFINED, reason="zero_denominator")
    return row


def _within(frame: pd.DataFrame, column: str, sequences: set | None) -> pd.DataFrame:
    return frame if sequences is None else frame[frame[column].isin(sequences)]


def _introns(models: GeneModels, kind: str, sequences: set | None) -> set:
    introns = unique_introns(models, kind)
    if sequences is None:
        return introns
    return {intron for intron in introns if intron[0] in sequences}


def _counterpart_ids(frame: pd.DataFrame) -> set:
    return {
        gid
        for ids in frame.loc[frame["counterpart_count"] >= 2, "counterpart_ids"]
        for gid in ids.split(",")
        if gid
    }


def _stop_rows(scope: str, stop_audit: pd.DataFrame | None, sequences) -> list:
    ids = [mid for mid in METRICS if mid.startswith("stop_codon_")]
    if stop_audit is None:
        return [
            _row(
                scope,
                mid,
                None,
                status=STATUS_NOT_APPLICABLE,
                reason="transformation_not_applied",
            )
            for mid in ids
        ]
    audit = _within(stop_audit, "chrom", sequences)
    total = len(audit)
    eligible = int(audit["eligible"].astype(bool).sum())
    counts = audit["status"].value_counts()
    return [
        _row(scope, "stop_codon_transcripts_with_cds", total),
        _row(scope, "stop_codon_eligible", eligible, total),
        *(
            _row(scope, f"stop_codon_{status}", int(counts.get(status, 0)), total)
            for status in ("extended", "unchanged", "skipped")
        ),
    ]


def scope_metric_rows(  # pylint: disable=too-many-arguments,too-many-locals
    scope: str,
    sequences: set | None,
    ref: pd.DataFrame,
    query: pd.DataFrame,
    reference_models: GeneModels,
    query_models: GeneModels,
    pairs: pd.DataFrame,
    edges: list[tuple[str, str, str]],
    transcripts: pd.DataFrame,
    stop_audit: pd.DataFrame | None,
    canonical_reference: bool,
) -> list[dict]:
    """
    Every metric of METRICS for one scope.
    Args:
            scope: Scope id ("all" or a declared analysis scope)
            sequences: The scope's sequence names; None for every compared sequence
            ref, query: classify_loci results with cds_locus_recovered,
                    cds_outcome and one_to_one_partner added by the runner
            reference_models, query_models: The evaluated GeneModels
            pairs: One-to-one pair table (matching.one_to_one_exact_cds)
            edges: (reference gene_id, query gene_id, sequence) of every
                    identical-CDS gene pair (exact_edges)
            transcripts: evaluated_transcripts table (both roles)
            stop_audit: stop-codon audit table, or None when not applied
            canonical_reference: The reference selection is canonical
    Returns:
            list of metrics.tsv rows (dicts with METRIC_COLUMNS)
    """
    ref = _within(ref, "chrom", sequences)
    query = _within(query, "chrom", sequences)
    r, q = len(ref), len(query)
    r_multi_cds = ref["max_cds_count"] > 1
    r_multi_exon = ref["max_exon_count"] > 1
    r_cds = ref["has_cds"].astype(bool)
    recovered = ref["cds_shared_bases"].astype("int64") > 0
    exact = _true(ref["cds_coordinate_exact"])
    q_exact = _true(query["cds_coordinate_exact"])
    detected = ref["classification"].isin(DETECTED)
    cds_detected = ref["classification_cds"].isin(DETECTED)
    tp = len(_within(pairs, "chrom", sequences))
    evaluated = transcripts[transcripts["status"] == "evaluated"]
    evaluated = _within(evaluated, "chrom", sequences)
    q_matched = query["classification"].isin(DETECTED + (MATCHED,))
    novel = query[query["classification"] == NOVEL]

    rows = [
        _row(scope, "reference_genes", r),
        _row(scope, "reference_genes_with_cds", int(r_cds.sum())),
        _row(scope, "reference_multi_cds_genes", int(r_multi_cds.sum())),
        _row(scope, "reference_multi_exon_genes", int(r_multi_exon.sum())),
        _row(scope, "query_genes", q),
        _row(scope, "query_genes_with_cds", int(query["has_cds"].astype(bool).sum())),
        _row(
            scope,
            "reference_transcripts_evaluated",
            int((evaluated["role"] == "reference").sum()),
        ),
        _row(
            scope,
            "query_transcripts_evaluated",
            int((evaluated["role"] == "query").sum()),
        ),
    ]
    if canonical_reference:
        fallback = evaluated[
            (evaluated["role"] == "reference")
            & evaluated["selected_by"].str.startswith("canonical_fallback")
        ]
        rows.append(
            _row(
                scope,
                "reference_canonical_fallback_genes",
                int(fallback["gene_id"].nunique()),
                r,
            )
        )
    else:
        rows.append(
            _row(
                scope,
                "reference_canonical_fallback_genes",
                None,
                status=STATUS_NOT_APPLICABLE,
                reason="selection_not_canonical",
            )
        )
    # Candidates come from the identical-CDS graph that the matching uses. They
    # equal the per-gene coordinate-exact flags unless a transcript lies outside
    # its gene record, so that its gene-span partners miss the identical gene.
    scoped_edges = [e for e in edges if sequences is None or e[2] in sequences]
    candidates_ref = len({e[0] for e in scoped_edges})
    candidates_query = len({e[1] for e in scoped_edges})
    rows += [
        _row(scope, "locus_recovery_cds_overlap", int(recovered.sum()), r),
        _row(scope, "cds_coordinate_exact_reference", int(exact.sum()), r),
        _row(scope, "cds_coordinate_exact_query", int(q_exact.sum()), q),
        _row(scope, "exact_cds_one_to_one_tp", tp),
        _row(scope, "exact_cds_one_to_one_recall", tp, r),
        _row(scope, "exact_cds_one_to_one_precision", tp, q),
        _row(scope, "exact_cds_one_to_one_f1", 2 * tp, r + q),
        _row(
            scope,
            "exact_cds_gene_pairs",
            len(scoped_edges),
        ),
        _row(
            scope,
            "exact_cds_reference_with_candidate_unmatched",
            candidates_ref - tp,
            r,
        ),
        _row(
            scope,
            "exact_cds_query_with_candidate_unmatched",
            candidates_query - tp,
            q,
        ),
        _row(
            scope,
            "cds_intron_chain_best_pair",
            int((cds_detected & _true(ref["cds_intron_chain_match"])).sum()),
            int(r_multi_cds.sum()),
        ),
        _row(
            scope,
            "cds_intron_chain_any_pair",
            int(_true(ref["cds_intron_chain_any_pair"]).sum()),
            int(r_multi_cds.sum()),
        ),
        _row(
            scope,
            "exon_intron_chain_best_pair",
            int((detected & _true(ref["intron_chain_match"])).sum()),
            int(r_multi_exon.sum()),
        ),
        _row(
            scope,
            "exon_intron_chain_any_pair",
            int(_true(ref["intron_chain_any_pair"]).sum()),
            int(r_multi_exon.sum()),
        ),
    ]
    for kind in ("cds", "exon"):
        ref_introns = _introns(reference_models, kind, sequences)
        query_introns = _introns(query_models, kind, sequences)
        shared = len(ref_introns & query_introns)
        rows += [
            _row(
                scope, f"ref_unique_{kind}_introns_recovered", shared, len(ref_introns)
            ),
            _row(
                scope,
                f"query_unique_{kind}_introns_supported",
                shared,
                len(query_introns),
            ),
        ]
    outcomes = ref["cds_outcome"].value_counts()
    rows += [
        _row(scope, f"cds_outcome_{key}", int(outcomes.get(key, 0)), r)
        for key in CDS_OUTCOMES
    ]
    r_single_cds = r_cds & ~r_multi_cds
    rows += [
        _row(scope, "partial_cds", int((recovered & ~exact).sum()), r),
        _row(
            scope,
            "cds_coordinate_exact_reference_multi_cds",
            int((exact & r_multi_cds).sum()),
            int(r_multi_cds.sum()),
        ),
        _row(
            scope,
            "cds_coordinate_exact_reference_single_cds",
            int((exact & r_single_cds).sum()),
            int(r_single_cds.sum()),
        ),
        _row(
            scope,
            "cds_structural_exact_match",
            int((ref["classification_cds"] == EXACT).sum()),
            r,
        ),
        _row(
            scope,
            "exon_coordinate_exact_reference",
            int(_true(ref["exon_coordinate_exact"]).sum()),
            r,
        ),
        _row(
            scope,
            "exon_structural_exact_match",
            int((ref["classification"] == EXACT).sum()),
            r,
        ),
        _row(scope, "locus_detection_gene_span", int(detected.sum()), r),
        _row(
            scope,
            "locus_detection_exonic",
            int((detected & (ref["exon_overlap"] > 0)).sum()),
            r,
        ),
        _row(
            scope,
            "span_only_detection",
            int((detected & ~(ref["exon_overlap"] > 0)).sum()),
            r,
        ),
        _row(
            scope, "detected_without_cds_overlap", int((detected & ~recovered).sum()), r
        ),
        _row(scope, "ref_genes_without_cds_overlap", int(r - recovered.sum()), r),
        _row(
            scope,
            "missed_reference_genes",
            int((ref["classification"] == MISSED).sum()),
            r,
        ),
        _row(
            scope,
            "strand_mismatch_reference",
            int((ref["classification"] == STRAND_MISMATCH).sum()),
            r,
        ),
        _row(scope, "query_matched", int(q_matched.sum()), q),
        _row(
            scope,
            "query_matched_without_cds_overlap",
            int((q_matched & ~(query["cds_shared_bases"].astype("int64") > 0)).sum()),
            q,
        ),
        _row(scope, "query_novel", len(novel), q),
    ]
    categories = novel["novel_category"].value_counts()
    rows += [
        _row(scope, f"query_novel_{c}", int(categories.get(c, 0)), q)
        for c in NOVEL_CATEGORIES
    ]
    rows += [
        _row(
            scope,
            "strand_mismatch_query",
            int((query["classification"] == STRAND_MISMATCH).sum()),
            q,
        ),
        _row(
            scope,
            "split_reference_genes",
            int((ref["counterpart_count"] >= 2).sum()),
            r,
        ),
        _row(scope, "split_query_genes_involved", len(_counterpart_ids(ref)), q),
        _row(
            scope, "merge_query_genes", int((query["counterpart_count"] >= 2).sum()), q
        ),
        _row(scope, "merge_reference_genes_involved", len(_counterpart_ids(query)), r),
        _row(
            scope,
            "one_counterpart_reference_genes",
            int((ref["counterpart_count"] == 1).sum()),
            r,
        ),
    ]
    rows += _stop_rows(scope, stop_audit, sequences)
    return rows


def _one_based(intervals) -> tuple:
    return tuple((start + 1, end) for start, end in intervals)


def splice_analysis(
    ref_models: GeneModels, query_models: GeneModels, ref, scopes: dict
) -> dict:
    """
    splice.analyse on the evaluated models, for scope "all" and every scope.
    Args:
            ref_models, query_models: The evaluated GeneModels (one reference
                    transcript per gene)
            ref: Reference results with cds_locus_recovered
            scopes: analysis scope id -> sequence names
    Returns:
            scope id -> splice.analyse result
    """
    genes = [
        splice.ReferenceGene(
            gene_id=gene.gene_id,
            chrom=gene.Chromosome,
            strand=gene.Strand,
            transcripts=[
                _one_based(ref_models.cds.get(tid, ())) for tid in gene.transcript_ids
            ],
        )
        for gene in ref_models.genes.itertuples(index=False)
    ]
    support = splice.SupportIndex.build(
        (gene.Chromosome, gene.Strand, _one_based(query_models.cds.get(tid, ())))
        for gene in query_models.genes.itertuples(index=False)
        for tid in gene.transcript_ids
    )
    detected = dict(zip(ref["gene_id"], ref["cds_locus_recovered"]))
    return splice.analyse(
        genes,
        support,
        detected,
        {ALL_SCOPE: None, **scopes},
        one_per_gene=True,
    )


def add_contract_columns(ref, query, pairs, scopes: dict) -> None:
    """
    Add the per-gene contract columns derived from classification and matching.

    cds_locus_recovered (shared CDS bases > 0), scope_ids, cds_outcome (reference
    rows) and one_to_one_partner (from the one-to-one pair table).
    """
    for frame in (ref, query):
        frame["cds_locus_recovered"] = frame["cds_shared_bases"].astype("int64") > 0
        frame["scope_ids"] = [
            ",".join(sid for sid, seqs in scopes.items() if chrom in seqs)
            for chrom in frame["chrom"]
        ]
    ref["cds_outcome"] = reference_cds_outcomes(ref)
    query["cds_outcome"] = ""
    ref["one_to_one_partner"] = (
        ref["gene_id"]
        .map(dict(zip(pairs["reference_gene_id"], pairs["query_gene_id"])))
        .fillna("")
    )
    query["one_to_one_partner"] = (
        query["gene_id"]
        .map(dict(zip(pairs["query_gene_id"], pairs["reference_gene_id"])))
        .fillna("")
    )


def exact_edges(reference: GeneModels, query: GeneModels) -> list[tuple[str, str, str]]:
    """(reference gene_id, query gene_id, sequence) of every identical-CDS gene pair."""
    ref_ids = reference.genes["gene_id"].tolist()
    query_ids = query.genes["gene_id"].tolist()
    chroms = reference.genes["Chromosome"].tolist()
    return [
        (ref_ids[ref_row], query_ids[query_row], chroms[ref_row])
        for ref_row, query_row in exact_cds_edges(reference, query)
    ]


def metric_table(rows: list[dict]) -> pd.DataFrame:
    """metrics.tsv as a DataFrame (integer counts, unrounded values)."""
    table = pd.DataFrame(rows, columns=METRIC_COLUMNS)
    for column in ("numerator", "denominator"):
        table[column] = table[column].astype("Int64")
    # Counts stay integers and ratios floats (no common float dtype).
    table["value"] = pd.Series(
        [row["value"] for row in rows], index=table.index, dtype=object
    )
    return table

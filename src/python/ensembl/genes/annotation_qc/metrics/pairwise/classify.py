"""
Locus pairing and transcript-structure classification.

Definitions (gmb-compare pairing and thresholds; reporting refined as noted):

    pairing      genes are paired when their spans overlap by at least 1 bp,
                 ignoring strand. A gene with a same-strand partner is compared
                 with its same-strand partners only. A gene whose only partners
                 are on the other strand is Strand_Mismatch when it shares
                 feature bases with one of them: CDS bases when both genes have
                 CDS, otherwise exon bases (strand_mismatch_basis "CDS" or
                 "exon"). Opposite-strand contact through gene spans, introns or
                 UTR alone does not count; such a gene is Missed (reference) or
                 Novel (query) with strand_mismatch_basis "no_feature_overlap". A
                 gene with no partner is Missed or Novel.
    exon pair    for every same-strand partner, every own transcript is compared
                 with every partner transcript that has exons:
                   reciprocal overlap = min(covered fraction of A by B,
                                            covered fraction of B by A)
                   intron chain match = identical ordered intron coordinates
                 The best pair ranks by identical exon coordinates, then priority
                 (3 = intron chain match and reciprocal overlap >= 0.8, 2 =
                 overlap >= 0.8, 1 = overlap > 0, 0 = no exonic overlap), then
                 reciprocal overlap; the first pair found wins ties.
    exon class   Exact_Match (priority 3), Structural_Mismatch (priority 2),
                 otherwise Partial_Match. A reference gene with a same-strand
                 partner but no comparable transcripts is Partial_Match; a query
                 gene in that situation is Matched. Exact_Match is a structural
                 class: terminal exon ends may differ. exon_coordinate_exact is
                 True only when the best pair has identical exon intervals.
    pair ids     matched_id / best_match_transcript_id name the partner gene and
                 transcript of the best exon pair, cds_matched_id /
                 best_cds_match_transcript_id those of the best CDS pair, and
                 best_match_own_transcript_id / best_cds_match_own_transcript_id
                 the gene's own transcript in each pair.
    CDS pair     chosen independently of the exon pair: every own transcript
                 with CDS is compared with every partner transcript with CDS and
                 the best pair is ranked as above on CDS intervals. A query CDS
                 that is identical to one reference isoform is therefore found
                 even when another isoform shares more UTR.
    CDS class    the exon-class rules applied to the best CDS pair; CDS
                 Exact_Match tolerates different start/stop coordinates.
                 cds_coordinate_exact is True only for identical CDS intervals,
                 terminal coordinates included. No_CDS when none of the gene's
                 own transcripts has CDS; Missed when it has CDS but no partner
                 transcript does.
    intron chain intron_chain_match / cds_intron_chain_match are reported for
                 the best exon / CDS pair and are None (not applicable) when the
                 gene's own compared structure has a single exon / CDS segment.
                 Two single-exon structures still rank as an intron-chain match
                 for pairing, as in gmb-compare.
    counterparts reference gene R and query gene Q are counterparts when they are
                 same-strand locus partners and some transcript pair shares at
                 least SPLIT_MERGE_MIN_SHARED_FRACTION of the shorter structure:
                 CDS when both genes have CDS, otherwise exons. Gene-span, intron
                 or UTR-only contact does not qualify. A reference gene with >= 2
                 query counterparts is a split; a query gene with >= 2 reference
                 counterparts is a merge.
    novel        a Novel query gene is further labelled (novel_category) by its
                 exon overlap, on either strand, with the unfiltered reference;
                 see novel_categories.

Coordinates are zero-based half-open, so lengths are End - Start and two
features overlap when each starts before the other ends. Partners are visited in
order of partner gene Start (ascending), End (descending), then file order;
transcripts in order of their first exon in the file.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pyranges1 as pr

EXACT, PARTIAL, STRUCTURAL = "Exact_Match", "Partial_Match", "Structural_Mismatch"
STRAND_MISMATCH, MISSED, NOVEL, MATCHED, NO_CDS = (
    "Strand_Mismatch",
    "Missed",
    "Novel",
    "Matched",
    "No_CDS",
)
RECIPROCAL_THRESHOLD = 0.8
SPLIT_MERGE_MIN_SHARED_FRACTION = 0.1

# strand_mismatch_basis values
BASIS_CDS, BASIS_EXON, BASIS_NO_FEATURE = "CDS", "exon", "no_feature_overlap"

# novel_category values, in order of precedence
NOVEL_PSEUDOGENE = "Overlaps_reference_pseudogene"
NOVEL_OTHER_NON_CODING = "Overlaps_other_non_protein_coding_reference"
NOVEL_UNEVALUATED_CODING = "Overlaps_unevaluated_protein_coding_reference"
NOVEL_OPPOSITE_STRAND = "Overlaps_evaluated_reference_opposite_strand"
NOVEL_NO_OVERLAP = "Novel_no_reference_overlap"
NOVEL_CATEGORIES = (
    NOVEL_PSEUDOGENE,
    NOVEL_OTHER_NON_CODING,
    NOVEL_UNEVALUATED_CODING,
    NOVEL_OPPOSITE_STRAND,
    NOVEL_NO_OVERLAP,
)

RESULT_COLUMNS = [
    "gene_id",
    "original_gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "transcript_ids",
    "max_exon_count",
    "max_cds_count",
    "classification",
    "classification_cds",
    "matched_id",
    "best_match_transcript_id",
    "exon_overlap",
    "intron_chain_match",
    "cds_overlap",
    "cds_intron_chain_match",
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
]

Intervals = tuple[tuple[int, int], ...]


@dataclass
class GeneModels:
    """Gene spans plus exon and CDS intervals per transcript."""

    genes: pd.DataFrame
    exons: dict[str, Intervals]
    cds: dict[str, Intervals]


def _intervals_by_transcript(df: pd.DataFrame, feature: str) -> dict[str, Intervals]:
    rows = df.loc[df["Feature"] == feature, ["transcript_id", "Start", "End"]]
    rows = rows.sort_values(["transcript_id", "Start", "End"], kind="stable")
    tids = rows["transcript_id"].to_numpy()
    if len(tids) == 0:
        return {}
    pairs = list(zip(rows["Start"].tolist(), rows["End"].tolist()))
    cuts = np.flatnonzero(tids[1:] != tids[:-1]) + 1
    bounds = zip(np.r_[0, cuts], np.r_[cuts, len(tids)])
    return {tids[a]: tuple(pairs[a:b]) for a, b in bounds}


def build_gene_models(df: pd.DataFrame) -> GeneModels:
    """
    Collect gene spans and per-transcript exon/CDS intervals.

    A gene's transcripts are those with at least one exon, ordered by their
    first exon in the file.
    Args:
            df: Comparison-schema DataFrame
    Returns:
            GeneModels
    """
    exons = _intervals_by_transcript(df, "exon")
    cds = _intervals_by_transcript(df, "CDS")
    first_exons = df.loc[
        df["Feature"] == "exon", ["gene_id", "transcript_id"]
    ].drop_duplicates("transcript_id")
    tids_by_gene = first_exons.groupby("gene_id", sort=False)["transcript_id"].agg(
        tuple
    )

    genes = df.loc[
        df["Feature"] == "gene",
        [
            "gene_id",
            "original_gene_id",
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "gene_biotype",
            "source_feature",
        ],
    ].reset_index(drop=True)
    genes["transcript_ids"] = [tids_by_gene.get(g, ()) for g in genes["gene_id"]]
    genes["max_exon_count"] = [
        max((len(exons[t]) for t in tids), default=0)
        for tids in genes["transcript_ids"]
    ]
    genes["max_cds_count"] = [
        max((len(cds.get(t, ())) for t in tids), default=0)
        for tids in genes["transcript_ids"]
    ]
    return GeneModels(genes, exons, cds)


def shared_bases(a: Intervals, b: Intervals) -> int:
    """Bases of a covered by intervals in b."""
    return sum(
        max(0, min(end_a, end_b) - max(start_a, start_b))
        for start_a, end_a in a
        for start_b, end_b in b
    )


def overlap_fraction(a: Intervals, b: Intervals) -> float:
    """Fraction of the bases in a that are covered by intervals in b."""
    total = sum(end - start for start, end in a)
    if total == 0:
        return 0.0
    return shared_bases(a, b) / total


def intron_chain(intervals: Intervals) -> tuple[tuple[int, int], ...]:
    """Ordered (intron start, intron end) pairs; empty for a single interval."""
    return tuple(
        (intervals[i][1], intervals[i + 1][0]) for i in range(len(intervals) - 1)
    )


def intron_chain_status(own: Intervals, other: Intervals) -> bool | None:
    """
    Reported intron-chain match of own against other.
    Returns:
            None when own has no intron (not applicable), otherwise whether the
            two intron chains are identical
    """
    if len(own) < 2:
        return None
    return intron_chain(own) == intron_chain(other)


def compare_structures(a: Intervals, b: Intervals) -> tuple[float, bool, str]:
    """
    Compare two interval chains.
    Args:
            a, b: Sorted intervals
    Returns:
            (reciprocal overlap, intron chain match, Exact/Structural/Partial class).
            The chain match is True for two single-interval chains; use
            intron_chain_status for reporting.
    """
    reciprocal = min(overlap_fraction(a, b), overlap_fraction(b, a))
    chain_match = intron_chain(a) == intron_chain(b)
    if reciprocal >= RECIPROCAL_THRESHOLD:
        return reciprocal, chain_match, EXACT if chain_match else STRUCTURAL
    return reciprocal, chain_match, PARTIAL


def _priority(reciprocal: float, chain_match: bool) -> int:
    if reciprocal >= RECIPROCAL_THRESHOLD:
        return 3 if chain_match else 2
    return 1 if reciprocal > 0 else 0


def pair_loci(
    own: pd.DataFrame, other: pd.DataFrame
) -> dict[int, list[tuple[int, bool]]]:
    """
    Find overlapping gene spans, ignoring strand.
    Args:
            own, other: Gene tables from build_gene_models
    Returns:
            own row -> [(other row, same strand)], ordered by other Start, End
            (descending), row
    """
    if own.empty or other.empty:
        return {}

    def ranges(genes: pd.DataFrame) -> pr.PyRanges:
        return pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": genes["Chromosome"].to_numpy(),
                    "Start": genes["Start"].to_numpy(),
                    "End": genes["End"].to_numpy(),
                    "Strand": genes["Strand"].to_numpy(),
                    "row": np.arange(len(genes)),
                }
            )
        )

    joined = pd.DataFrame(
        ranges(own).join_overlaps(
            ranges(other), strand_behavior="ignore", suffix="_other"
        )
    )
    if joined.empty:
        return {}
    # Partner order decides ties: Start ascending, End descending (as gmb-compare's
    # pyranges 0.x join), then file order where gmb-compare's order was arbitrary.
    joined = joined.sort_values(
        ["row", "Start_other", "End_other", "row_other"],
        ascending=[True, True, False, True],
    )
    same = (joined["Strand"].astype(str) == joined["Strand_other"].astype(str)).tolist()
    pairs: dict[int, list[tuple[int, bool]]] = {}
    for row, other_row, same_strand in zip(
        joined["row"].tolist(), joined["row_other"].tolist(), same
    ):
        pairs.setdefault(row, []).append((other_row, same_strand))
    return pairs


def _gene_intervals(models: GeneModels, tids: tuple, kind: str) -> Intervals:
    """All exon or CDS intervals of a gene's transcripts (may repeat)."""
    source = models.exons if kind == "exon" else models.cds
    return tuple(iv for tid in tids for iv in source.get(tid, ()))


def _unpaired(label: str, max_exons: int, max_cds: int) -> dict:
    return {
        "classification": label,
        "classification_cds": label,
        "matched_id": "",
        "best_match_transcript_id": "",
        "exon_overlap": 0.0,
        "intron_chain_match": None if max_exons <= 1 else False,
        "cds_overlap": 0.0,
        "cds_intron_chain_match": None if max_cds <= 1 else False,
        "exon_coordinate_exact": False,
        "cds_coordinate_exact": False,
        "cds_matched_id": "",
        "best_cds_match_transcript_id": "",
        "best_match_own_transcript_id": "",
        "best_cds_match_own_transcript_id": "",
        "strand_mismatch_basis": "",
    }


def _best_match(
    own_tids: tuple,
    partners: list[int],
    own: GeneModels,
    other: GeneModels,
    default: str,
    max_exons: int,
    max_cds: int,
) -> dict:
    """Best exon pair and, independently, best CDS pair over same-strand partners."""
    other_ids = other.genes["gene_id"].tolist()
    other_tids = other.genes["transcript_ids"].tolist()
    best = _unpaired(PARTIAL, max_exons, max_cds)
    best["matched_id"] = other_ids[partners[0]]
    exon_key = cds_key = None
    for partner in partners:
        for own_tid in own_tids:
            own_exons = own.exons.get(own_tid)
            if not own_exons:
                continue
            own_cds = own.cds.get(own_tid, ())
            for other_tid in other_tids[partner]:
                other_exons = other.exons.get(other_tid)
                if not other_exons:
                    continue
                reciprocal, chain_match, exon_class = compare_structures(
                    own_exons, other_exons
                )
                identical = own_exons == other_exons
                key = (identical, _priority(reciprocal, chain_match), reciprocal)
                if exon_key is None or key > exon_key:
                    exon_key = key
                    best.update(
                        classification=exon_class,
                        matched_id=other_ids[partner],
                        best_match_transcript_id=other_tid,
                        best_match_own_transcript_id=own_tid,
                        exon_overlap=reciprocal,
                        intron_chain_match=intron_chain_status(own_exons, other_exons),
                        exon_coordinate_exact=identical,
                    )
                other_cds = other.cds.get(other_tid, ())
                if not (own_cds and other_cds):
                    continue
                cds_overlap, cds_chain, cds_class = compare_structures(
                    own_cds, other_cds
                )
                identical = own_cds == other_cds
                key = (identical, _priority(cds_overlap, cds_chain), cds_overlap)
                if cds_key is None or key > cds_key:
                    cds_key = key
                    best.update(
                        classification_cds=cds_class,
                        cds_matched_id=other_ids[partner],
                        best_cds_match_transcript_id=other_tid,
                        best_cds_match_own_transcript_id=own_tid,
                        cds_overlap=cds_overlap,
                        cds_intron_chain_match=intron_chain_status(own_cds, other_cds),
                        cds_coordinate_exact=identical,
                    )
    if exon_key is None:
        # No comparable transcripts: CDS stays Partial_Match, as in gmb-compare.
        best["classification"] = default
    elif cds_key is None:
        best["classification_cds"] = MISSED if max_cds else NO_CDS
    return best


def _opposite_strand_overlap(
    own_tids: tuple, partners: list[int], own: GeneModels, other: GeneModels
) -> tuple[int, str] | None:
    """
    First opposite-strand partner sharing feature bases with the own gene.

    CDS is used when both genes have CDS, exons otherwise; a CDS overlap with
    any partner is preferred over an exon overlap with an earlier one.
    Returns:
            (partner row, "CDS" or "exon"), or None
    """
    other_tids = other.genes["transcript_ids"].tolist()
    own_cds = _gene_intervals(own, own_tids, "CDS")
    own_exons = _gene_intervals(own, own_tids, "exon")
    exon_hit = None
    for partner in partners:
        other_cds = _gene_intervals(other, other_tids[partner], "CDS")
        if own_cds and other_cds:
            if shared_bases(own_cds, other_cds) > 0:
                return partner, BASIS_CDS
        elif exon_hit is None and shared_bases(
            own_exons, _gene_intervals(other, other_tids[partner], "exon")
        ):
            exon_hit = (partner, BASIS_EXON)
    return exon_hit


def _classify_side(
    own: GeneModels,
    other: GeneModels,
    pairs: dict[int, list[tuple[int, bool]]],
    unpaired_label: str,
    no_transcript_label: str,
) -> pd.DataFrame:
    other_ids = other.genes["gene_id"].tolist()
    records = []
    for row, (own_tids, max_exons, max_cds) in enumerate(
        zip(
            own.genes["transcript_ids"].tolist(),
            own.genes["max_exon_count"].tolist(),
            own.genes["max_cds_count"].tolist(),
        )
    ):
        candidates = pairs.get(row, [])
        same = [partner for partner, same_strand in candidates if same_strand]
        if not candidates:
            result = _unpaired(unpaired_label, max_exons, max_cds)
        elif same:
            result = _best_match(
                own_tids, same, own, other, no_transcript_label, max_exons, max_cds
            )
        else:
            opposite = [partner for partner, _ in candidates]
            hit = _opposite_strand_overlap(own_tids, opposite, own, other)
            if hit:
                result = _unpaired(STRAND_MISMATCH, max_exons, max_cds)
                result["matched_id"] = other_ids[hit[0]]
                result["strand_mismatch_basis"] = hit[1]
            else:
                result = _unpaired(unpaired_label, max_exons, max_cds)
                result["strand_mismatch_basis"] = BASIS_NO_FEATURE
        records.append(result)

    genes = own.genes.rename(
        columns={
            "Chromosome": "chrom",
            "Start": "start",
            "End": "end",
            "Strand": "strand",
        }
    )
    classified = pd.DataFrame.from_records(
        records, columns=list(_unpaired(MISSED, 0, 0))
    )
    for column in ("intron_chain_match", "cds_intron_chain_match"):
        classified[column] = classified[column].astype("boolean")
    classified["novel_category"] = ""
    classified["counterpart_count"] = 0
    classified["counterpart_ids"] = ""
    return pd.concat([genes.reset_index(drop=True), classified], axis=1)[RESULT_COLUMNS]


def _best_shared_fraction(
    ref_tids: tuple, query_tids: tuple, reference: GeneModels, query: GeneModels
) -> float:
    """Largest shared fraction of the shorter structure over transcript pairs."""
    use_cds = bool(
        _gene_intervals(reference, ref_tids, "CDS")
        and _gene_intervals(query, query_tids, "CDS")
    )
    source_ref = reference.cds if use_cds else reference.exons
    source_query = query.cds if use_cds else query.exons
    best = 0.0
    for ref_tid in ref_tids:
        a = source_ref.get(ref_tid)
        if not a:
            continue
        len_a = sum(end - start for start, end in a)
        for query_tid in query_tids:
            b = source_query.get(query_tid)
            if not b:
                continue
            shorter = min(len_a, sum(end - start for start, end in b))
            if shorter:
                best = max(best, shared_bases(a, b) / shorter)
    return best


def find_counterparts(
    reference: GeneModels,
    query: GeneModels,
    ref_pairs: dict[int, list[tuple[int, bool]]] | None = None,
) -> list[tuple[int, int]]:
    """
    Same-strand locus partners that share real coding (or exonic) sequence.
    Args:
            reference, query: GeneModels
            ref_pairs: pair_loci(reference.genes, query.genes), if already computed
    Returns:
            [(reference row, query row)] in reference row, partner order
    """
    if ref_pairs is None:
        ref_pairs = pair_loci(reference.genes, query.genes)
    ref_tids = reference.genes["transcript_ids"].tolist()
    query_tids = query.genes["transcript_ids"].tolist()
    links = []
    for ref_row in sorted(ref_pairs):
        for query_row, same_strand in ref_pairs[ref_row]:
            if not same_strand:
                continue
            fraction = _best_shared_fraction(
                ref_tids[ref_row], query_tids[query_row], reference, query
            )
            if fraction >= SPLIT_MERGE_MIN_SHARED_FRACTION:
                links.append((ref_row, query_row))
    return links


def _add_counterparts(
    results: pd.DataFrame, other_ids: list, links: list[tuple[int, int]]
) -> None:
    by_row: dict[int, list[str]] = {}
    for own_row, other_row in links:
        by_row.setdefault(own_row, []).append(other_ids[other_row])
    results["counterpart_count"] = [len(by_row.get(i, ())) for i in range(len(results))]
    results["counterpart_ids"] = [
        ",".join(by_row.get(i, ())) for i in range(len(results))
    ]


def _is_pseudogene(biotype: str, source_feature: str) -> bool:
    return "pseudogene" in biotype.lower() or source_feature == "pseudogene"


def novel_categories(
    query: GeneModels,
    novel_rows: list[int],
    context: GeneModels,
    evaluated_gene_ids: set,
) -> dict[int, str]:
    """
    Say what, if anything, each Novel query gene overlaps in the reference.

    Overlap is shared exon bases on either strand (gene spans when a gene has no
    exons) with any reference gene in context, which should be the reference
    before biotype/mode filtering. Categories, by precedence:

        Overlaps_reference_pseudogene                  filtered-out pseudogene
        Overlaps_other_non_protein_coding_reference    filtered-out non-coding gene
        Overlaps_unevaluated_protein_coding_reference  filtered-out coding gene
        Overlaps_evaluated_reference_opposite_strand   evaluated gene, other strand,
                                                       without qualifying CDS overlap
        Novel_no_reference_overlap                     nothing
    Args:
            query: Query GeneModels
            novel_rows: Query gene rows classified Novel
            context: Reference GeneModels including filtered-out genes
            evaluated_gene_ids: gene_ids of the reference genes that were compared
    Returns:
            query row -> novel_category
    """
    if not novel_rows:
        return {}
    pairs = pair_loci(query.genes, context.genes)
    query_tids = query.genes["transcript_ids"].tolist()
    ctx = context.genes
    ctx_ids = ctx["gene_id"].tolist()
    ctx_tids = ctx["transcript_ids"].tolist()
    ctx_biotype = ctx["gene_biotype"].astype(str).tolist()
    ctx_feature = ctx["source_feature"].astype(str).tolist()
    rank = {category: i for i, category in enumerate(NOVEL_CATEGORIES)}
    categories = {}
    for row in novel_rows:
        own_exons = _gene_intervals(query, query_tids[row], "exon")
        found = NOVEL_NO_OVERLAP
        for ctx_row, _ in pairs.get(row, []):
            ctx_exons = _gene_intervals(context, ctx_tids[ctx_row], "exon")
            if own_exons and ctx_exons and not shared_bases(own_exons, ctx_exons):
                continue
            if ctx_ids[ctx_row] in evaluated_gene_ids:
                category = NOVEL_OPPOSITE_STRAND
            elif _is_pseudogene(ctx_biotype[ctx_row], ctx_feature[ctx_row]):
                category = NOVEL_PSEUDOGENE
            elif ctx_biotype[ctx_row] == "protein_coding":
                category = NOVEL_UNEVALUATED_CODING
            else:
                category = NOVEL_OTHER_NON_CODING
            if rank[category] < rank[found]:
                found = category
        categories[row] = found
    return categories


def classify_loci(
    reference: GeneModels,
    query: GeneModels,
    reference_context: GeneModels | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Classify every reference gene (sensitivity) and every query gene (novelty).
    Args:
            reference: GeneModels for the filtered reference
            query: GeneModels for the query
            reference_context: GeneModels for the reference before biotype/mode
                    filtering, used only to label Novel query genes; defaults
                    to reference
    Returns:
            (reference results, query results), one row per gene, RESULT_COLUMNS
    """
    ref_pairs = pair_loci(reference.genes, query.genes)
    query_pairs = pair_loci(query.genes, reference.genes)
    ref_results = _classify_side(reference, query, ref_pairs, MISSED, PARTIAL)
    query_results = _classify_side(query, reference, query_pairs, NOVEL, MATCHED)

    links = find_counterparts(reference, query, ref_pairs)
    _add_counterparts(ref_results, query.genes["gene_id"].tolist(), links)
    _add_counterparts(
        query_results,
        reference.genes["gene_id"].tolist(),
        [(query_row, ref_row) for ref_row, query_row in links],
    )

    novel_rows = query_results.index[query_results["classification"] == NOVEL].tolist()
    context = reference if reference_context is None else reference_context
    for row, category in novel_categories(
        query, novel_rows, context, set(reference.genes["gene_id"])
    ).items():
        query_results.at[row, "novel_category"] = category
    return ref_results, query_results


def gene_span_excess(models: GeneModels) -> pd.DataFrame:
    """
    Measure how far each gene span extends beyond its transcripts' exons.

    Locus pairing uses gene spans, so a gene wider than its transcripts can
    pair with reference genes it shares no sequence with.
    Args:
            models: GeneModels
    Returns:
            DataFrame with gene_id and excess_bp (0 when the span is exact;
            genes without exons are omitted)
    """
    records = []
    for gene in models.genes.itertuples(index=False):
        exons = [iv for tid in gene.transcript_ids for iv in models.exons.get(tid, ())]
        if not exons:
            continue
        first = min(start for start, _ in exons)
        last = max(end for _, end in exons)
        excess = max(0, first - gene.Start) + max(0, gene.End - last)
        records.append((gene.gene_id, excess))
    return pd.DataFrame(records, columns=["gene_id", "excess_bp"]).astype(
        {"excess_bp": "int64"}
    )

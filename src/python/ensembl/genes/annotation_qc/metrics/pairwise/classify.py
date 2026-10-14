"""
Locus pairing and transcript-structure classification.

Definitions (gmb-compare pairing and thresholds; reporting refined as noted):

    pairing      genes are paired when their pairing spans overlap by at least
                 1 bp, ignoring strand. A gene's pairing span is its gene record
                 extended to cover the exons and CDS of its evaluated
                 transcripts; it equals the gene record whenever those
                 transcripts lie within it, as GFF3 and GTF require. (Before
                 definition 2026-10.1 the gene record alone was used, so a
                 transcript beyond its gene record was never compared with the
                 genes it overlaps.) Reported gene coordinates are always the
                 gene record's. A gene with a same-strand partner is compared
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
    any pair     cds_intron_chain_any_pair / intron_chain_any_pair: True when
                 ANY own transcript with > 1 CDS segment (exon) and ANY
                 transcript of a same-strand partner have identical, non-empty
                 CDS (exon) intron chains; None when no own transcript has an
                 intron. Unlike the best-pair columns this does not depend on
                 pair ranking, so any-pair >= best-pair recovery.
    best CDS pair descriptors (reported from the query's point of view on both
                 sides): cds_shared_bases (bases shared by the two CDS),
                 cds_start_diff_bp / cds_stop_diff_bp (> 0: the query CDS
                 extends further than the reference CDS at the start / stop end,
                 in the direction of transcription), cds_introns_missing
                 (reference CDS introns absent from the query transcript) and
                 cds_introns_extra (query CDS introns absent from the reference
                 transcript). None when the gene has no best CDS pair.

Coordinates are zero-based half-open, so lengths are End - Start and two
features overlap when each starts before the other ends. Partners are visited in
order of partner gene record Start (ascending), End (descending), then file order;
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
    # Appended for output contract 1.0 (see contract.py); earlier columns unchanged.
    "gene_biotype",
    "evaluated_transcript_count",
    "has_cds",
    "cds_shared_bases",
    "cds_intron_chain_any_pair",
    "intron_chain_any_pair",
    "cds_start_diff_bp",
    "cds_stop_diff_bp",
    "cds_introns_missing",
    "cds_introns_extra",
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
    genes["pair_start"], genes["pair_end"] = _pairing_spans(
        genes, first_exons, exons, cds
    )
    return GeneModels(genes, exons, cds)


def _pairing_spans(
    genes: pd.DataFrame,
    first_exons: pd.DataFrame,
    exons: dict[str, Intervals],
    cds: dict[str, Intervals],
) -> tuple[np.ndarray, np.ndarray]:
    """Gene record spans extended to the exons and CDS of the gene's transcripts."""
    starts = genes["Start"].to_numpy(copy=True)
    ends = genes["End"].to_numpy(copy=True)
    row_of = dict(zip(genes["gene_id"], range(len(genes))))
    for gene_id, tid in zip(first_exons["gene_id"], first_exons["transcript_id"]):
        row = row_of.get(gene_id)
        if row is None:
            continue
        for intervals in (exons.get(tid, ()), cds.get(tid, ())):
            if intervals:
                starts[row] = min(starts[row], intervals[0][0])
                ends[row] = max(ends[row], max(end for _, end in intervals))
    return starts, ends


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
    Find overlapping pairing spans (gene record spans for tables without
    pair_start/pair_end), ignoring strand.
    Args:
            own, other: Gene tables from build_gene_models
    Returns:
            own row -> [(other row, same strand)], ordered by other gene
            record Start, End (descending), row
    """
    if own.empty or other.empty:
        return {}

    def ranges(genes: pd.DataFrame) -> pr.PyRanges:
        start, end = (
            ("pair_start", "pair_end") if "pair_start" in genes else ("Start", "End")
        )
        return pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": genes["Chromosome"].to_numpy(),
                    "Start": genes[start].to_numpy(),
                    "End": genes[end].to_numpy(),
                    "Strand": genes["Strand"].to_numpy(),
                    "row": np.arange(len(genes)),
                    "record_start": genes["Start"].to_numpy(),
                    "record_end": genes["End"].to_numpy(),
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
    # Partner order decides ties: gene record Start ascending, End descending (as
    # gmb-compare's pyranges 0.x join), then file order where gmb-compare's order
    # was arbitrary. Pairing spans only decide which genes overlap.
    joined = joined.sort_values(
        ["row", "record_start_other", "record_end_other", "row_other"],
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
        "cds_shared_bases": 0,
        "cds_intron_chain_any_pair": None if max_cds <= 1 else False,
        "intron_chain_any_pair": None if max_exons <= 1 else False,
        "cds_start_diff_bp": None,
        "cds_stop_diff_bp": None,
        "cds_introns_missing": None,
        "cds_introns_extra": None,
    }


def cds_differences(
    reference_cds: Intervals, query_cds: Intervals, strand: str
) -> tuple[int, int, int, int]:
    """
    Describe how a query CDS differs from a reference CDS.
    Args:
            reference_cds, query_cds: Sorted CDS intervals of one transcript each
            strand: Strand of the reference gene ("-" reverses start and stop)
    Returns:
            (start_diff_bp, stop_diff_bp, reference introns missing from the
            query, query introns absent from the reference). A positive
            difference means the query CDS extends further at that end, in the
            direction of transcription.
    """
    left = reference_cds[0][0] - query_cds[0][0]
    right = query_cds[-1][1] - reference_cds[-1][1]
    start_diff, stop_diff = (right, left) if strand == "-" else (left, right)
    ref_introns = set(intron_chain(reference_cds))
    query_introns = set(intron_chain(query_cds))
    return (
        start_diff,
        stop_diff,
        len(ref_introns - query_introns),
        len(query_introns - ref_introns),
    )


def _best_match(  # pylint: disable=too-many-locals,too-many-branches
    own_tids: tuple,
    partners: list[int],
    own: GeneModels,
    other: GeneModels,
    default: str,
    max_exons: int,
    max_cds: int,
    own_is_reference: bool = True,
    strand: str = "+",
) -> dict:
    """Best exon pair and, independently, best CDS pair over same-strand partners."""
    other_ids = other.genes["gene_id"].tolist()
    other_tids = other.genes["transcript_ids"].tolist()
    best = _unpaired(PARTIAL, max_exons, max_cds)
    best["matched_id"] = other_ids[partners[0]]
    exon_key = cds_key = None
    best_cds_pair: tuple = ()
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
                if len(own_exons) > 1 and intron_chain(own_exons) == intron_chain(
                    other_exons
                ):
                    best["intron_chain_any_pair"] = True
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
                if len(own_cds) > 1 and intron_chain(own_cds) == intron_chain(
                    other_cds
                ):
                    best["cds_intron_chain_any_pair"] = True
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
                    best_cds_pair = (own_cds, other_cds)
    if exon_key is None:
        # No comparable transcripts: CDS stays Partial_Match, as in gmb-compare.
        best["classification"] = default
    elif cds_key is None:
        best["classification_cds"] = MISSED if max_cds else NO_CDS
    if best_cds_pair:
        own_cds, other_cds = best_cds_pair
        best["cds_shared_bases"] = shared_bases(own_cds, other_cds)
        ref_cds, query_cds = (
            (own_cds, other_cds) if own_is_reference else (other_cds, own_cds)
        )
        (
            best["cds_start_diff_bp"],
            best["cds_stop_diff_bp"],
            best["cds_introns_missing"],
            best["cds_introns_extra"],
        ) = cds_differences(ref_cds, query_cds, strand)
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


def _classify_side(  # pylint: disable=too-many-locals
    own: GeneModels,
    other: GeneModels,
    pairs: dict[int, list[tuple[int, bool]]],
    unpaired_label: str,
    no_transcript_label: str,
) -> pd.DataFrame:
    other_ids = other.genes["gene_id"].tolist()
    own_is_reference = unpaired_label == MISSED
    records = []
    for row, (own_tids, max_exons, max_cds, strand) in enumerate(
        zip(
            own.genes["transcript_ids"].tolist(),
            own.genes["max_exon_count"].tolist(),
            own.genes["max_cds_count"].tolist(),
            own.genes["Strand"].astype(str).tolist(),
        )
    ):
        candidates = pairs.get(row, [])
        same = [partner for partner, same_strand in candidates if same_strand]
        if not candidates:
            result = _unpaired(unpaired_label, max_exons, max_cds)
        elif same:
            result = _best_match(
                own_tids,
                same,
                own,
                other,
                no_transcript_label,
                max_exons,
                max_cds,
                own_is_reference,
                strand,
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
    for column in (
        "intron_chain_match",
        "cds_intron_chain_match",
        "intron_chain_any_pair",
        "cds_intron_chain_any_pair",
    ):
        classified[column] = classified[column].astype("boolean")
    for column in (
        "cds_start_diff_bp",
        "cds_stop_diff_bp",
        "cds_introns_missing",
        "cds_introns_extra",
    ):
        classified[column] = classified[column].astype("Int64")
    classified["novel_category"] = ""
    classified["counterpart_count"] = 0
    classified["counterpart_ids"] = ""
    genes = genes.reset_index(drop=True)
    genes["evaluated_transcript_count"] = [len(t) for t in genes["transcript_ids"]]
    genes["has_cds"] = genes["max_cds_count"] > 0
    return pd.concat([genes, classified], axis=1)[RESULT_COLUMNS]


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


def transcripts_outside_gene(models: GeneModels) -> pd.DataFrame:
    """
    Transcripts whose exons extend beyond their gene record.

    Input diagnostic: GFF3/GTF require transcripts within their gene. Locus
    pairing covers such transcripts through the gene's pairing span (see
    pairing_span_beyond_record).
    Args:
            models: GeneModels
    Returns:
            DataFrame with gene_id, transcript_id and the bases outside the gene
            span (outside_bp), one row per such transcript
    """
    records = []
    for gene in models.genes.itertuples(index=False):
        for tid in gene.transcript_ids:
            exons = models.exons.get(tid, ())
            if not exons:
                continue
            outside = max(0, gene.Start - exons[0][0]) + max(0, exons[-1][1] - gene.End)
            if outside:
                records.append((gene.gene_id, tid, outside))
    return pd.DataFrame(
        records, columns=["gene_id", "transcript_id", "outside_bp"]
    ).astype({"outside_bp": "int64"})


def pairing_span_beyond_record(models: GeneModels) -> pd.DataFrame:
    """
    Genes whose pairing span extends beyond their gene record.
    Args:
            models: GeneModels
    Returns:
            DataFrame with gene_id and beyond_bp (bases of the pairing span
            outside the gene record), one row per such gene
    """
    genes = models.genes
    beyond = (genes["Start"] - genes["pair_start"]) + (genes["pair_end"] - genes["End"])
    wider = beyond > 0
    return pd.DataFrame(
        {"gene_id": genes.loc[wider, "gene_id"], "beyond_bp": beyond[wider]}
    ).astype({"beyond_bp": "int64"})


def gene_span_excess(models: GeneModels) -> pd.DataFrame:
    """
    Measure how far each gene span extends beyond its transcripts' exons.

    Locus pairing uses pairing spans, never narrower than the gene record, so a
    gene wider than its transcripts can pair with reference genes it shares no
    sequence with.
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

"""
One-to-one coordinate-exact CDS matching (precision/recall/F1) and CDS-overlap
locus recovery.

Fixture coordinates are one-based inclusive GFF coordinates.
"""

import csv
import json

import pytest

from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    build_gene_models,
    classify_loci,
)
from ensembl.genes.annotation_qc.metrics.pairwise.matching import (
    PAIR_COLUMNS,
    maximum_matching,
    one_to_one_exact_cds,
)
from ensembl.genes.annotation_qc.metrics.pairwise.selection import (
    apply_evaluation_mode,
    select_transcripts,
    subset_to_regions,
)
from ensembl.genes.annotation_qc.metrics.pairwise.summary import summarise_comparison
from ensembl.genes.annotation_qc.parsers.annotation import (
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.regions import parse_region

from .conftest import gene_records
from .test_pairwise_runner import _run

X = [(100, 200), (300, 400)]
Y = [(1100, 1200), (1300, 1400)]


@pytest.fixture
def match(write_gff3):
    """Parse, optionally transform, classify and match; return the pieces."""

    def _match(ref_records, query_records, ref_transform=None, query_transform=None):
        ref = parse_annotation_for_comparison(write_gff3("ref.gff3", ref_records))
        query = parse_annotation_for_comparison(write_gff3("query.gff3", query_records))
        if ref_transform:
            ref = ref_transform(ref)
        if query_transform:
            query = query_transform(query)
        ref_models, query_models = build_gene_models(ref), build_gene_models(query)
        ref_results, query_results = classify_loci(ref_models, query_models)
        summary = summarise_comparison(ref_results, query_results)
        one_to_one, pairs = one_to_one_exact_cds(ref_models, query_models)
        return one_to_one, pairs, summary, ref_results.set_index("gene_id")

    return _match


def _pairs(pairs):
    return sorted(zip(pairs["reference_gene_id"], pairs["query_gene_id"]))


def _greedy(edges, order):
    """First-free-candidate pairing, for contrast with the maximum matching."""
    used, matched = set(), {}
    for left in order:
        for right in edges.get(left, ()):
            if right not in used:
                used.add(right)
                matched[left] = right
                break
    return matched


def test_maximum_matching_beats_greedy():
    edges = {"R1": ["Q1", "Q2"], "R2": ["Q1"]}
    assert len(_greedy(edges, ["R1", "R2"])) == 1
    assert maximum_matching(edges, ["R1", "R2"]) == {"R1": "Q2", "R2": "Q1"}
    # Longer augmenting path: R3 displaces R2, which displaces R1.
    chain = {"R1": ["Q1", "Q2"], "R2": ["Q1", "Q3"], "R3": ["Q3"]}
    assert len(maximum_matching(chain, ["R1", "R2", "R3"])) == 3


def test_competing_exact_matches_where_greedy_fails(match):
    # R1 has isoforms with CDS X and Y; R2 has CDS X only (its gene ends later,
    # so R1 is processed first and tries Q1 first). Greedy: R1-Q1, R2 unmatched.
    ref = gene_records(
        "1", "R1", "+", {"R1.x": X, "R1.y": Y}, with_cds={"R1.x": X, "R1.y": Y}
    ) + gene_records("1", "R2", "+", {"R2.1": X + [(500, 1500)]}, with_cds={"R2.1": X})
    query = gene_records("1", "Q1", "+", {"Q1.1": X}) + gene_records(
        "1", "Q2", "+", {"Q2.1": Y}
    )
    one, pairs, summary, ref_results = match(ref, query)
    assert one["exact_gene_pairs"] == 3
    assert one["true_positives"] == 2
    assert _pairs(pairs) == [("R1", "Q2"), ("R2", "Q1")]
    assert (one["precision"], one["recall"], one["f1"]) == (1.0, 1.0, 1.0)
    # The per-gene best CDS pairs both name Q1: a deduplication of those best
    # pairs would give TP = 1.
    assert set(ref_results["cds_matched_id"]) == {"Q1"}
    assert one["exact_components_with_several_edges"] == 1
    assert set(pairs["component_reference_genes"]) == {2}
    assert set(pairs["component_query_genes"]) == {2}


def test_several_matching_isoforms_count_once(match):
    # R.a and R.b differ only in UTR; Q has two transcripts matching R.a/R.b and R.c.
    utr = [(50, 200), (300, 400)]
    c = [(100, 200), (300, 450)]
    ref = gene_records(
        "1",
        "R",
        "+",
        {"R.a": X, "R.b": utr, "R.c": [(100, 200), (300, 500)]},
        with_cds={"R.a": X, "R.b": X, "R.c": c},
    )
    query = gene_records("1", "Q", "+", {"Q.1": X, "Q.2": [(100, 200), (300, 450)]})
    one, pairs, _, _ = match(ref, query)
    assert one["true_positives"] == 1 and one["exact_gene_pairs"] == 1
    row = pairs.iloc[0]
    assert row.shared_cds_signatures == 2
    assert row.supporting_reference_transcript_ids == "R.a,R.b,R.c"
    assert row.supporting_query_transcript_ids == "Q.1,Q.2"
    # Lowest signature reported, one-based conversion happens in the report.
    assert (row.cds_start, row.cds_end, row.cds_segments) == (99, 400, 2)


def test_duplicated_cds_across_genes(match):
    ref = gene_records("1", "R1", "+", {"R1.1": X}) + gene_records(
        "1", "R2", "+", {"R2.1": X}
    )
    query = gene_records("1", "Q", "+", {"Q.1": X})
    one, pairs, summary, _ = match(ref, query)
    # Both reference genes are coordinate-exact per gene, but one query gene
    # can only be one true positive.
    assert summary["sensitivity_cds"]["cds_coordinate_exact_count"] == 2
    assert summary["specificity"]["consensus_cds_coordinate_exact_count"] == 1
    assert one["reference_genes_with_exact_candidate"] == 2
    assert one["true_positives"] == 1
    assert one["unmatched_reference_genes"] == 1
    assert one["unmatched_query_genes"] == 0
    assert one["reference_genes_with_exact_candidate_unmatched"] == 1
    assert (one["recall"], one["precision"]) == (0.5, 1.0)
    assert one["f1"] == round(2 * 1 / (2 + 1), 4)
    assert _pairs(pairs) == [("R1", "Q")]  # deterministic: R1 sorts first
    assert pairs.iloc[0].query_exact_candidates == 2

    # The reverse: duplicated query models.
    one, _, summary, _ = match(
        query,
        ref,
    )
    assert summary["specificity"]["consensus_cds_coordinate_exact_count"] == 2
    assert one["true_positives"] == 1 and one["unmatched_query_genes"] == 1


def test_strands_and_sequences(match):
    minus = [(5100, 5200), (5300, 5400)]
    ref = (
        gene_records("1", "P", "+", {"P.1": X})
        + gene_records("1", "M", "-", {"M.1": minus})
        + gene_records("1", "S", "+", {"S.1": [(9000, 9100)]})
        + gene_records("2", "G1", "+", {"G1.1": X})
    )
    query = (
        gene_records("1", "QP", "-", {"QP.1": X})  # same coordinates, wrong strand
        + gene_records("1", "QM", "-", {"QM.1": minus})
        + gene_records("1", "QS", "+", {"QS.1": [(9000, 9100)]})
        + gene_records("3", "G1", "+", {"G1.1": X})  # same ID and CDS, other sequence
    )
    one, pairs, _, _ = match(ref, query)
    assert _pairs(pairs) == [("M", "QM"), ("S", "QS")]
    assert one["true_positives"] == 2
    assert set(pairs["strand"]) == {"+", "-"}
    single = pairs.set_index("reference_gene_id").loc["S"]
    assert single.cds_segments == 1


def test_sequence_scoped_ids(match):
    # The same identifiers on two sequences are namespaced by the parser and
    # must pair within their own sequence only.
    ref = gene_records("1", "g1", "+", {"g1.t1": X}) + gene_records(
        "2", "g1", "+", {"g1.t1": Y}
    )
    query = gene_records("1", "g1", "+", {"g1.t1": Y}) + gene_records(
        "2", "g1", "+", {"g1.t1": Y}
    )
    one, pairs, _, _ = match(ref, query)
    assert one["true_positives"] == 1
    row = pairs.iloc[0]
    assert (row.reference_gene_id, row.query_gene_id, row.chrom) == (
        "2:g1",
        "2:g1",
        "2",
    )
    assert (row.reference_original_gene_id, row.query_original_gene_id) == ("g1", "g1")
    assert row.supporting_query_transcript_ids == "2:g1.t1"


def test_single_segment_and_stop_codon_semantics(match):
    ref = gene_records("1", "R", "+", {"R.1": [(100, 400)]}) + gene_records(
        "1", "T", "+", {"T.1": [(1000, 1300)]}
    )
    # T': stop codon excluded (3 bp shorter) is not coordinate-exact.
    query = gene_records("1", "Q", "+", {"Q.1": [(100, 400)]}) + gene_records(
        "1", "QT", "+", {"QT.1": [(1000, 1297)]}
    )
    one, pairs, summary, _ = match(ref, query)
    assert _pairs(pairs) == [("R", "Q")]
    assert summary["sensitivity_cds"]["cds_exact_match_count"] == 2
    assert one["true_positives"] == 1


def test_genes_without_cds_stay_in_denominators(match):
    ref = gene_records("1", "R", "+", {"R.1": X})
    query = gene_records("1", "Q", "+", {"Q.1": X}) + gene_records(
        "1", "N", "+", {"N.1": [(3000, 3500)]}, with_cds=False
    )
    one, _, _, _ = match(ref, query)
    assert (one["query_genes"], one["query_genes_with_cds"]) == (2, 1)
    assert one["true_positives"] == 1 and one["unmatched_query_genes"] == 1
    assert one["precision"] == 0.5


def test_filtering_and_selection_change_the_population(match):
    longer = [(100, 200), (300, 600)]
    ref = (
        gene_records("1", "R", "+", {"R.short": X, "R.long": longer})
        + gene_records(
            "1", "L", "+", {"L.1": [(2000, 2100)]}, biotype="lncRNA", with_cds=False
        )
        + gene_records("2", "R2", "+", {"R2.1": X})
    )
    query = gene_records("1", "Q", "+", {"Q.1": X}) + gene_records(
        "2", "Q2", "+", {"Q2.1": X}
    )
    one, _, _, _ = match(ref, query)
    assert (one["reference_genes"], one["true_positives"]) == (3, 2)
    one, _, _, _ = match(
        ref, query, ref_transform=lambda df: apply_evaluation_mode(df, "protein_coding")
    )
    assert (one["reference_genes"], one["true_positives"]) == (2, 2)
    one, _, _, _ = match(
        ref,
        query,
        ref_transform=lambda df: select_transcripts(
            apply_evaluation_mode(df, "protein_coding"), "longest_cds"
        ),
    )
    # Only R.long is evaluated for R, so Q no longer matches it.
    assert (one["reference_genes"], one["true_positives"]) == (2, 1)
    region = [parse_region("2")]
    one, _, _, _ = match(
        ref,
        query,
        ref_transform=lambda df: subset_to_regions(df, region),
        query_transform=lambda df: subset_to_regions(df, region),
    )
    assert (one["reference_genes"], one["query_genes"], one["true_positives"]) == (
        1,
        1,
        1,
    )


def test_empty_populations_are_undefined_not_zero(match):
    records = gene_records("1", "R", "+", {"R.1": X})
    empty = lambda df: df.iloc[0:0]  # noqa: E731
    one, pairs, summary, _ = match(records, records, query_transform=empty)
    assert (one["query_genes"], one["true_positives"]) == (0, 0)
    assert one["precision"] is None
    assert one["recall"] == 0.0 and one["f1"] == 0.0
    assert list(pairs.columns) == PAIR_COLUMNS and pairs.empty
    one, _, summary, _ = match(records, records, ref_transform=empty)
    assert one["recall"] is None and one["precision"] == 0.0
    assert summary["cds_overlap_locus_recovery"]["rate"] is None
    one, _, _, _ = match(records, records, ref_transform=empty, query_transform=empty)
    assert (one["precision"], one["recall"], one["f1"]) == (None, None, None)


def test_cds_overlap_locus_recovery(match):
    long_cds = [(100, 30099)]
    ref = (
        gene_records("1", "A", "+", {"A.1": long_cds})  # shares 1 CDS base with QA
        + gene_records("1", "B", "+", {"B.1": [(40000, 40300)]})  # span-only partner
        + gene_records("1", "C", "+", {"C.1": [(50000, 50300)]})  # missed
    )
    query = gene_records("1", "QA", "+", {"QA.1": [(30099, 30200)]}) + gene_records(
        "1",
        "QB",
        "+",
        {"QB.1": [(39000, 39100), (40400, 40500)]},
    )
    _, _, summary, ref_results = match(ref, query)
    recovery = summary["cds_overlap_locus_recovery"]
    assert recovery["recovered_count"] == 1
    assert recovery["reference_genes"] == 3
    assert recovery["rate"] == round(1 / 3, 4)
    assert recovery["gene_span_detected_without_cds_overlap"] == 1
    assert recovery["recovered_count"] == int((ref_results["cds_overlap"] > 0).sum())
    # The tiny overlap is real but rounds to 0 at 4 decimals (details TSV).
    assert 0 < ref_results.loc["A", "cds_overlap"] < 0.00005
    assert summary["sensitivity"]["locus_detected_count"] == 2


def test_runner_writes_pairs_and_summary(write_gff3, tmp_path):
    ref = write_gff3(
        "ref.gff3",
        gene_records(
            "1", "R1", "+", {"R1.x": X, "R1.y": Y}, with_cds={"R1.x": X, "R1.y": Y}
        )
        + gene_records(
            "1", "R2", "+", {"R2.1": X + [(500, 1500)]}, with_cds={"R2.1": X}
        ),
    )
    query = write_gff3(
        "query.gff3",
        gene_records("1", "Q1", "+", {"Q1.1": X})
        + gene_records("1", "Q2", "+", {"Q2.1": Y})
        + gene_records("1", "Q3", "-", {"Q3.1": [(5000, 5100)]}),
    )
    out = tmp_path / "out"
    _run(["--query", query, "--reference", ref, "--outdir", str(out)])
    summary = json.loads((out / "comparison_summary.json").read_text())
    one = summary["cds_exact_one_to_one"]
    assert (one["true_positives"], one["reference_genes"], one["query_genes"]) == (
        2,
        2,
        3,
    )
    assert one["precision"] == round(2 / 3, 4) and one["f1"] == 0.8
    assert one["definition"]["f1"] == "2 TP / (R + Q)"
    assert summary["cds_overlap_locus_recovery"]["recovered_count"] == 2
    assert "cds_exact_one_to_one.f1" in summary["rate_denominators"]
    # Historical keys are unchanged.
    assert summary["sensitivity_cds"]["cds_coordinate_exact_count"] == 2
    with open(out / "cds_exact_one_to_one_pairs.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["reference_gene_id"] + "-" + r["query_gene_id"] for r in rows] == [
        "R1-Q2",
        "R2-Q1",
    ]
    r2 = rows[1]
    assert (
        r2["cds_start"],
        r2["cds_end"],
        r2["reference_start"],
        r2["reference_end"],
    ) == (
        "100",
        "400",
        "100",
        "1500",
    )
    tsv = dict(
        line.rstrip("\n").split("\t")
        for line in (out / "comparison_summary.tsv").read_text().splitlines()[1:]
    )
    assert tsv["cds_exact_one_to_one_true_positives"] == "2"
    assert tsv["cds_overlap_locus_recovery_rate"] == "1.0"

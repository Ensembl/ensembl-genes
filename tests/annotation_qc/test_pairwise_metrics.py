"""Locus pairing, classification, selection policies and summary denominators."""

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    build_gene_models,
    classify_loci,
    gene_span_excess,
    overlap_fraction,
)
from ensembl.genes.annotation_qc.metrics.pairwise.selection import (
    apply_evaluation_mode,
    select_transcripts,
    subset_to_regions,
)
from ensembl.genes.annotation_qc.metrics.pairwise.summary import (
    query_transcript_labels,
    summarise_comparison,
)
from ensembl.genes.annotation_qc.parsers.annotation import (
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.regions import Region

from .conftest import gene_records

MULTI = [(100, 200), (300, 400), (500, 600)]


@pytest.fixture
def compare(write_gff3):
    """Parse reference and query record lists, classify, and index results by gene."""

    def _compare(ref_records, query_records):
        ref = parse_annotation_for_comparison(write_gff3("ref.gff3", ref_records))
        query = parse_annotation_for_comparison(write_gff3("query.gff3", query_records))
        ref_results, query_results = classify_loci(
            build_gene_models(ref), build_gene_models(query)
        )
        return ref_results.set_index("gene_id"), query_results.set_index("gene_id")

    return _compare


def _shift(exons, by):
    return [(s + by, e + by) for s, e in exons]


class TestClassification:
    def test_exact_match(self, compare):
        ref, query = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}),
            gene_records("1", "Q", "+", {"Q.1": MULTI}),
        )
        row = ref.loc["R"]
        assert (row.classification, row.classification_cds) == (
            "Exact_Match",
            "Exact_Match",
        )
        assert row.exon_overlap == 1.0 and row.intron_chain_match
        assert (row.matched_id, row.best_match_transcript_id) == ("Q", "Q.1")
        assert query.loc["Q"].classification == "Exact_Match"

    def test_structural_mismatch_high_overlap_different_introns(self, compare):
        query_exons = [(100, 200), (300, 410), (500, 600)]
        query_exons[1] = (310, 400)  # moves an intron boundary, keeps >= 0.8 overlap
        ref, _ = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}),
            gene_records("1", "Q", "+", {"Q.1": query_exons}),
        )
        assert ref.loc["R"].classification == "Structural_Mismatch"
        assert not ref.loc["R"].intron_chain_match

    def test_partial_match_low_reciprocal_overlap(self, compare):
        ref, _ = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}),
            gene_records("1", "Q", "+", {"Q.1": [(100, 200)]}),
        )
        assert ref.loc["R"].classification == "Partial_Match"
        assert 0 < ref.loc["R"].exon_overlap < 0.8

    def test_missed_and_novel(self, compare):
        ref, query = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}),
            gene_records("1", "Q", "+", {"Q.1": _shift(MULTI, 10000)}),
        )
        assert ref.loc["R"].classification == "Missed"
        assert query.loc["Q"].classification == "Novel"

    def test_strand_mismatch_both_sides(self, compare):
        ref, query = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}),
            gene_records("1", "Q", "-", {"Q.1": MULTI}),
        )
        assert ref.loc["R"].classification == "Strand_Mismatch"
        assert ref.loc["R"].matched_id == "Q"
        assert query.loc["Q"].classification == "Strand_Mismatch"

    def test_reverse_strand_exact(self, compare):
        ref, _ = compare(
            gene_records("1", "R", "-", {"R.1": MULTI}),
            gene_records("1", "Q", "-", {"Q.1": MULTI}),
        )
        assert ref.loc["R"].classification_cds == "Exact_Match"

    def test_single_exon_exact_and_single_vs_multi(self, compare):
        ref, _ = compare(
            gene_records("1", "R1", "+", {"R1.1": [(100, 900)]})
            + gene_records("1", "R2", "+", {"R2.1": [(2000, 2900)]}),
            gene_records("1", "Q1", "+", {"Q1.1": [(100, 900)]})
            + gene_records("1", "Q2", "+", {"Q2.1": [(2000, 2400), (2450, 2900)]}),
        )
        # Two single-exon structures are still Exact_Match, but have no intron
        # chain to recover: the intron-chain result is not applicable (NA).
        assert ref.loc["R1"].classification == "Exact_Match"
        assert ref.loc["R1"].intron_chain_match is pd.NA
        assert ref.loc["R2"].classification == "Structural_Mismatch"

    def test_multi_isoform_best_pair_chosen(self, compare):
        ref, _ = compare(
            gene_records(
                "1", "R", "+", {"R.1": [(100, 200), (500, 600)], "R.2": MULTI}
            ),
            gene_records("1", "Q", "+", {"Q.1": [(100, 250)], "Q.2": MULTI}),
        )
        assert ref.loc["R"].classification == "Exact_Match"
        assert ref.loc["R"].best_match_transcript_id == "Q.2"

    def test_cds_compared_on_best_exon_pair(self, compare):
        cds_ref = {"R.1": [(150, 200), (300, 400), (500, 550)]}
        cds_query = {"Q.1": [(150, 200), (300, 400), (500, 520)]}
        ref, _ = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}, with_cds=cds_ref),
            gene_records("1", "Q", "+", {"Q.1": MULTI}, with_cds=cds_query),
        )
        row = ref.loc["R"]
        assert row.classification == "Exact_Match"
        assert row.classification_cds == "Exact_Match"  # same introns, >= 0.8 overlap
        # GFF lengths: reference 51 + 101 + 51 = 203 bp, query 51 + 101 + 21 = 173 bp.
        assert row.cds_overlap == pytest.approx(173 / 203)

    def test_no_cds_labels(self, compare):
        ref, query = compare(
            gene_records("1", "R1", "+", {"R1.1": MULTI}, with_cds=False)
            + gene_records("1", "R2", "+", {"R2.1": _shift(MULTI, 5000)}),
            gene_records("1", "Q1", "+", {"Q1.1": MULTI})
            + gene_records(
                "1", "Q2", "+", {"Q2.1": _shift(MULTI, 5000)}, with_cds=False
            ),
        )
        assert ref.loc["R1"].classification_cds == "No_CDS"  # reference has no CDS
        assert ref.loc["R2"].classification_cds == "Missed"  # query lacks CDS
        assert query.loc["Q2"].classification_cds == "No_CDS"

    def test_merge_one_query_gene_over_two_reference_genes(self, compare):
        ref, query = compare(
            gene_records("1", "R1", "+", {"R1.1": [(100, 200), (300, 400)]})
            + gene_records("1", "R2", "+", {"R2.1": [(600, 700), (800, 900)]}),
            gene_records(
                "1", "Q", "+", {"Q.1": [(100, 200), (300, 400), (600, 700), (800, 900)]}
            ),
        )
        assert set(ref["matched_id"]) == {"Q"}
        assert set(ref["classification"]) == {"Partial_Match"}
        assert query.loc["Q"].matched_id == "R1"  # first of two equal-overlap partners

    def test_split_two_query_genes_over_one_reference_gene(self, compare):
        ref, query = compare(
            gene_records(
                "1", "R", "+", {"R.1": [(100, 200), (300, 400), (600, 700), (800, 900)]}
            ),
            gene_records("1", "Q1", "+", {"Q1.1": [(100, 200), (300, 400)]})
            + gene_records("1", "Q2", "+", {"Q2.1": [(600, 700), (800, 900)]}),
        )
        assert ref.loc["R"].classification == "Partial_Match"
        assert set(query["matched_id"]) == {"R"}

    def test_one_base_shared_at_boundary_pairs(self, compare):
        # GFF 100..200 and 200..300 share base 200.
        ref, _ = compare(
            gene_records("1", "R", "+", {"R.1": [(100, 200)]}),
            gene_records("1", "Q", "-", {"Q.1": [(200, 300)]}),
        )
        assert ref.loc["R"].classification == "Strand_Mismatch"

    def test_adjacent_genes_do_not_pair(self, compare):
        ref, _ = compare(
            gene_records("1", "R", "+", {"R.1": [(100, 199)]}),
            gene_records("1", "Q", "+", {"Q.1": [(200, 300)]}),
        )
        assert ref.loc["R"].classification == "Missed"

    def test_duplicate_tiberius_ids_do_not_pair_across_chromosomes(self, compare):
        ref, query = compare(
            gene_records("1", "R1", "+", {"R1.1": MULTI}),
            gene_records("1", "g1", "+", {"g1.t1": MULTI})
            + gene_records("2", "g1", "+", {"g1.t1": MULTI}),
        )
        assert ref.loc["R1"].matched_id == "1:g1"
        assert query.loc["2:g1"].classification == "Novel"
        assert query.loc["2:g1"].original_gene_id == "g1"

    def test_span_only_detection_is_reported(self, compare, write_gff3):
        # Query gene record spans the reference gene but its exons are elsewhere.
        query = gene_records("1", "Q", "+", {"Q.1": [(5000, 5100)]})
        query[0] = ("1", "gene", 50, 5100, "+", "ID=Q;biotype=protein_coding")
        ref, query_results = compare(gene_records("1", "R", "+", {"R.1": MULTI}), query)
        assert ref.loc["R"].classification == "Partial_Match"
        assert ref.loc["R"].exon_overlap == 0.0
        summary = summarise_comparison(ref.reset_index(), query_results.reset_index())
        assert summary["sensitivity"]["locus_detected_count"] == 1
        assert summary["locus_detection_exonic"]["span_only_detected_count"] == 1
        parsed = parse_annotation_for_comparison(write_gff3("q.gff3", query))
        excess = gene_span_excess(build_gene_models(parsed)).set_index("gene_id")
        assert excess.loc["Q", "excess_bp"] == 4950


def test_overlap_fraction_half_open():
    assert overlap_fraction(((0, 10),), ((5, 20),)) == 0.5
    assert overlap_fraction(((0, 10),), ((10, 20),)) == 0.0


class TestSelection:
    @pytest.fixture
    def reference(self, write_gff3):
        records = gene_records(
            "1",
            "PC",
            "+",
            {"PC.1": [(100, 200)], "PC.2": [(100, 300)]},
            with_cds={"PC.1": [(100, 200)], "PC.2": [(100, 150)]},
        )
        records[2] = (*records[2][:5], records[2][5] + ";tag=Ensembl_canonical")
        records += gene_records(
            "1", "NC", "+", {"NC.1": [(1000, 1200)]}, biotype="lncRNA", with_cds=False
        )
        records += gene_records(
            "1", "TIE", "+", {"TIE.1": [(2000, 2100)], "TIE.2": [(2000, 2100)]}
        )
        return parse_annotation_for_comparison(write_gff3("ref.gff3", records))

    def _transcripts(self, df):
        return sorted(df.loc[df["Feature"] == "transcript", "transcript_id"])

    def test_protein_coding_keeps_whole_coding_genes(self, reference):
        kept = apply_evaluation_mode(reference, "protein_coding")
        assert set(kept.loc[kept["Feature"] == "gene", "gene_id"]) == {"PC", "TIE"}
        assert "PC.2" in self._transcripts(kept)

    def test_cds_only(self, reference):
        kept = apply_evaluation_mode(reference, "cds_only")
        assert "NC" not in set(kept["gene_id"])

    def test_canonical_tag_then_longest_cds(self, reference):
        kept = apply_evaluation_mode(reference, "canonical")
        assert self._transcripts(kept) == ["NC.1", "PC.1", "TIE.1"]

    def test_longest_cds_ignores_tag_and_breaks_ties_by_file_order(self, reference):
        kept = select_transcripts(reference, "longest_cds")
        # PC.1 CDS 101 bp > PC.2 51 bp; NC has no CDS -> longest exons; TIE -> first.
        assert self._transcripts(kept) == ["NC.1", "PC.1", "TIE.1"]

    def test_regions_keep_whole_genes(self, reference):
        kept = subset_to_regions(reference, [Region("1", 250, 260)])
        assert set(kept["gene_id"]) == {"PC"}
        assert len(kept[kept["Feature"] == "exon"]) == 2


def test_summary_denominators_and_labels(compare):
    ref, query = compare(
        gene_records("1", "R1", "+", {"R1.1": MULTI})
        + gene_records("1", "R2", "+", {"R2.1": _shift(MULTI, 5000)})
        + gene_records("1", "R3", "+", {"R3.1": _shift(MULTI, 10000)})
        + gene_records("1", "R4", "+", {"R4.1": _shift(MULTI, 20000)}),
        gene_records("1", "Q1", "+", {"Q1.1": MULTI, "Q1.2": MULTI[:2]})
        + gene_records("1", "Q2", "+", {"Q2.1": _shift(MULTI, 5000)[:1]})
        + gene_records("1", "Q3", "-", {"Q3.1": _shift(MULTI, 10000)})
        + gene_records("1", "Q4", "+", {"Q4.1": _shift(MULTI, 30000)}),
    )
    summary = summarise_comparison(ref.reset_index(), query.reset_index())
    assert summary["total_reference_genes"] == 4
    assert summary["reference_classification"] == {
        "Exact_Match": 1,
        "Partial_Match": 1,
        "Strand_Mismatch": 1,
        "Missed": 1,
    }
    sens = summary["sensitivity"]
    assert (sens["locus_detected_count"], sens["locus_detection_rate"]) == (2, 0.5)
    assert sens["exact_match_rate"] == 0.25 and sens["missed_rate"] == 0.25
    # CDS rates use all reference genes as the denominator.
    assert summary["sensitivity_cds"]["cds_exact_match_rate"] == 0.25
    # Intron-chain rate excludes Missed and Strand_Mismatch genes; the
    # sensitivity divides by every multi-exon reference gene.
    assert summary["intron_chain"] == {
        "exon_intron_chain_recovered": 1,
        "exon_intron_chain_matched": 2,
        "exon_intron_chain_evaluated": 2,
        "exon_intron_chain_not_applicable_single_exon": 0,
        "exon_intron_chain_rate": 0.5,
        "exon_intron_chain_sensitivity": 0.25,
        "multi_exon_reference_genes": 4,
    }
    spec = summary["specificity"]
    assert (spec["novel_consensus_count"], spec["matched_consensus_count"]) == (1, 2)
    assert summary["stratified"]["multi_exon"]["Exact_Match"] == 1

    labels = query_transcript_labels(query.reset_index()).set_index("transcript_id")
    assert labels["classification"].to_dict() == {
        "Q1.1": "Matched",
        "Q1.2": "Matched",
        "Q2.1": "Matched",
        "Q3.1": "Strand_Mismatch",
        "Q4.1": "Novel",
    }

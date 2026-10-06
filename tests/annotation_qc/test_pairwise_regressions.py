"""
Regression tests for the comparator fixes found in the P. falciparum
Vipsania/Tiberius evaluation: matched_consensus_count, coordinate-exact metrics,
single-exon intron chains, independent CDS pairing, feature-aware strand
mismatch, gene splits/merges, Novel context and rate denominators.

Fixture coordinates are one-based inclusive GFF coordinates.
"""

import csv
import json

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    build_gene_models,
    classify_loci,
)
from ensembl.genes.annotation_qc.metrics.pairwise.selection import (
    apply_evaluation_mode,
)
from ensembl.genes.annotation_qc.metrics.pairwise.summary import (
    query_transcript_labels,
    summarise_comparison,
    summarise_intron_support,
)
from ensembl.genes.annotation_qc.parsers.annotation import (
    parse_annotation_for_comparison,
)

from .conftest import gene_records
from .test_pairwise_runner import _run

MULTI = [(100, 200), (300, 400), (500, 600)]


@pytest.fixture
def compare(write_gff3):
    """
    Parse, optionally filter the reference to protein_coding (keeping the
    unfiltered reference as Novel context), classify and summarise.
    """

    def _compare(ref_records, query_records, protein_coding=False):
        ref_all = parse_annotation_for_comparison(write_gff3("ref.gff3", ref_records))
        ref = (
            apply_evaluation_mode(ref_all, "protein_coding")
            if protein_coding
            else ref_all
        )
        query = parse_annotation_for_comparison(write_gff3("query.gff3", query_records))
        ref_models, query_models = build_gene_models(ref), build_gene_models(query)
        ref_results, query_results = classify_loci(
            ref_models, query_models, reference_context=build_gene_models(ref_all)
        )
        summary = summarise_comparison(ref_results, query_results)
        summary["intron_support"] = summarise_intron_support(ref_models, query_models)
        return (
            ref_results.set_index("gene_id"),
            query_results.set_index("gene_id"),
            summary,
        )

    return _compare


def _shift(exons, by):
    return [(s + by, e + by) for s, e in exons]


def _gene_without_structure(seqid, gene, strand, start, end):
    """A transcript without exon or CDS rows: the gene has no comparable transcript."""
    return [
        (seqid, "gene", start, end, strand, f"ID={gene};biotype=protein_coding"),
        (seqid, "mRNA", start, end, strand, f"ID={gene}.1;Parent={gene}"),
    ]


# 1. matched_consensus_count -------------------------------------------------


class TestMatchedConsensusCount:
    def test_matched_class_counted_once(self, compare):
        # Q pairs with R on the same strand but R has no exon/CDS rows, so Q is
        # "Matched". gmb-compare's formula, Matched + total - Novel -
        # Strand_Mismatch, gave 2 for this single query gene.
        _, query, summary = compare(
            _gene_without_structure("1", "R", "+", 100, 600),
            gene_records("1", "Q", "+", {"Q.1": MULTI}),
        )
        assert query.loc["Q"].classification == "Matched"
        assert summary["total_consensus_genes"] == 1
        assert summary["specificity"]["matched_consensus_count"] == 1

    def test_every_query_class(self, compare):
        ref, query, summary = compare(
            gene_records("1", "R1", "+", {"R1.1": MULTI})
            + gene_records("1", "R2", "+", {"R2.1": _shift(MULTI, 2000)})
            + gene_records("1", "R3", "+", {"R3.1": _shift(MULTI, 4000)})
            + _gene_without_structure("1", "R4", "+", 6100, 6600),
            gene_records("1", "Q1", "+", {"Q1.1": MULTI})  # Exact_Match
            + gene_records("1", "Q2", "+", {"Q2.1": _shift(MULTI, 2000)[:1]})  # Partial
            + gene_records("1", "Q3", "-", {"Q3.1": _shift(MULTI, 4000)})  # Strand
            + gene_records("1", "Q4", "+", {"Q4.1": _shift(MULTI, 6000)})  # Matched
            + gene_records("1", "Q5", "+", {"Q5.1": _shift(MULTI, 9000)}),  # Novel
        )
        assert query["classification"].to_dict() == {
            "Q1": "Exact_Match",
            "Q2": "Partial_Match",
            "Q3": "Strand_Mismatch",
            "Q4": "Matched",
            "Q5": "Novel",
        }
        spec = summary["specificity"]
        assert spec["matched_consensus_count"] == 3  # Q1, Q2, Q4
        assert spec["matched_consensus_count"] <= summary["total_consensus_genes"]
        assert (
            spec["matched_consensus_count"]
            + spec["novel_consensus_count"]
            + spec["strand_mismatch_consensus_count"]
            == summary["total_consensus_genes"]
        )
        assert spec["matched_consensus_rate"] == 0.6


# 2-3. coordinate-exact CDS / exons -------------------------------------------


class TestCoordinateExact:
    CDS = [(150, 200), (300, 400), (500, 550)]

    def _run(self, compare, query_cds, query_exons=MULTI):
        ref, query, summary = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}, with_cds={"R.1": self.CDS}),
            gene_records(
                "1", "Q", "+", {"Q.1": query_exons}, with_cds={"Q.1": query_cds}
            ),
        )
        return ref.loc["R"], summary

    def test_identical_coordinates(self, compare):
        row, summary = self._run(compare, self.CDS)
        assert row.classification_cds == "Exact_Match"
        assert row.cds_coordinate_exact and row.exon_coordinate_exact
        assert summary["sensitivity_cds"]["cds_coordinate_exact_count"] == 1
        assert summary["sensitivity"]["exon_coordinate_exact_count"] == 1

    def test_same_intron_chain_different_start(self, compare):
        # Start codon 10 bp downstream: identical CDS intron chain, >= 0.8 overlap.
        row, summary = self._run(compare, [(160, 200), (300, 400), (500, 550)])
        assert row.classification_cds == "Exact_Match"  # structural class kept
        assert row.cds_intron_chain_match
        assert not row.cds_coordinate_exact
        assert row.exon_coordinate_exact  # exons are still identical
        sens_cds = summary["sensitivity_cds"]
        assert sens_cds["cds_exact_match_count"] == 1
        assert sens_cds["cds_coordinate_exact_count"] == 0
        assert sens_cds["cds_exact_match_not_coordinate_exact"] == 1

    def test_different_terminal_exon_end(self, compare):
        query_exons = [(100, 200), (300, 400), (500, 580)]
        row, _ = self._run(compare, self.CDS, query_exons)
        assert row.classification == "Exact_Match"  # intron chain + >= 0.8 overlap
        assert not row.exon_coordinate_exact
        assert row.cds_coordinate_exact

    def test_partial_overlap(self, compare):
        row, _ = self._run(compare, [(150, 200)], [(100, 200)])
        assert row.classification_cds == "Partial_Match"
        assert not row.cds_coordinate_exact and not row.exon_coordinate_exact


# 4-5. single-exon intron-chain semantics -------------------------------------


class TestSingleExonIntronChain:
    def test_single_exon_cases_are_not_applicable(self, compare):
        ref, _, summary = compare(
            gene_records("1", "EXACT", "+", {"EXACT.1": [(100, 900)]})
            + gene_records("1", "PART", "+", {"PART.1": [(2000, 2900)]})
            + gene_records("1", "NONE", "+", {"NONE.1": [(5000, 5900)]}),
            gene_records("1", "Q1", "+", {"Q1.1": [(100, 900)]})
            + gene_records("1", "Q2", "+", {"Q2.1": [(2000, 2300)]}),
        )
        assert ref.loc["EXACT"].classification == "Exact_Match"
        assert ref.loc["PART"].classification == "Partial_Match"
        assert ref.loc["NONE"].classification == "Missed"
        for gene in ("EXACT", "PART", "NONE"):
            assert ref.loc[gene].intron_chain_match is pd.NA
            assert ref.loc[gene].cds_intron_chain_match is pd.NA
        chain, cds = summary["intron_chain"], summary["sensitivity_cds"]
        assert chain["exon_intron_chain_recovered"] == 0
        assert chain["exon_intron_chain_evaluated"] == 0
        assert chain["exon_intron_chain_not_applicable_single_exon"] == 2
        assert chain["multi_exon_reference_genes"] == 0
        assert chain["exon_intron_chain_rate"] == 0.0
        assert cds["cds_intron_chain_recovered"] == 0
        assert cds["multi_segment_cds_reference_genes"] == 0

    def test_multi_exon_exact_and_single_exon_query_on_multi_exon_reference(
        self, compare
    ):
        ref, query, summary = compare(
            gene_records("1", "M", "+", {"M.1": MULTI})
            + gene_records("1", "N", "+", {"N.1": _shift(MULTI, 2000)}),
            gene_records("1", "Q1", "+", {"Q1.1": MULTI})
            + gene_records("1", "Q2", "+", {"Q2.1": [(2100, 2600)]}),
        )
        assert ref.loc["M"].intron_chain_match
        assert ref.loc["M"].cds_intron_chain_match
        # Reference introns not reproduced by a single-exon query: False, not NA.
        assert ref.loc["N"].intron_chain_match is not pd.NA
        assert not ref.loc["N"].intron_chain_match
        # The single-exon query has no introns of its own: NA on the query side.
        assert query.loc["Q2"].intron_chain_match is pd.NA
        chain = summary["intron_chain"]
        assert (
            chain["exon_intron_chain_recovered"],
            chain["exon_intron_chain_evaluated"],
        ) == (1, 2)
        assert chain["exon_intron_chain_sensitivity"] == 0.5


# 6. CDS pair chosen independently of the exon pair ---------------------------


def test_cds_pair_independent_of_exon_pair(compare):
    query_exons = MULTI
    query_cds = [(150, 200), (300, 400), (500, 550)]
    # Isoform A: identical exons (best exon/UTR match), CDS starts in exon 2.
    # Isoform B: longer UTRs, CDS identical to the query's.
    ref, query, summary = compare(
        gene_records(
            "1",
            "R",
            "+",
            {"R.A": MULTI, "R.B": [(50, 200), (300, 400), (500, 700)]},
            with_cds={"R.A": [(320, 400), (500, 550)], "R.B": query_cds},
        ),
        gene_records("1", "Q", "+", {"Q.1": query_exons}, with_cds={"Q.1": query_cds}),
    )
    q = query.loc["Q"]
    assert q.best_match_transcript_id == "R.A"
    assert q.exon_coordinate_exact
    assert q.best_cds_match_transcript_id == "R.B"
    assert q.cds_matched_id == "R"
    assert q.classification_cds == "Exact_Match" and q.cds_coordinate_exact
    assert q.cds_overlap == 1.0
    r = ref.loc["R"]
    assert r.best_match_transcript_id == "Q.1" and r.classification == "Exact_Match"
    assert r.cds_coordinate_exact and r.classification_cds == "Exact_Match"
    # From the reference side the switch is in the reference's own isoform.
    assert r.best_match_own_transcript_id == "R.A"
    assert r.best_cds_match_own_transcript_id == "R.B"
    assert summary["sensitivity_cds"]["cds_pair_differs_from_exon_pair"] == 1
    assert summary["sensitivity_cds"]["cds_coordinate_exact_count"] == 1
    assert summary["specificity"]["consensus_cds_coordinate_exact_count"] == 1
    assert summary["specificity"]["consensus_exon_coordinate_exact_count"] == 1

    labels = query_transcript_labels(query.reset_index()).set_index("transcript_id")
    assert labels.loc["Q.1", "best_ref_transcript_id"] == "R.A"
    assert labels.loc["Q.1", "best_cds_ref_transcript_id"] == "R.B"


def test_no_cds_uses_any_coding_isoform(compare):
    # The best exon pair is the reference's non-coding isoform; the coding
    # isoform still gives the CDS comparison (gmb-compare reported No_CDS).
    ref, _, _ = compare(
        gene_records(
            "1",
            "R",
            "+",
            {"R.nc": MULTI, "R.pc": [(150, 200), (300, 400)]},
            with_cds={"R.pc": [(150, 200), (300, 400)]},
        ),
        gene_records(
            "1", "Q", "+", {"Q.1": MULTI}, with_cds={"Q.1": [(150, 200), (300, 400)]}
        ),
    )
    assert ref.loc["R"].best_match_transcript_id == "Q.1"
    assert ref.loc["R"].classification_cds == "Exact_Match"
    assert ref.loc["R"].cds_coordinate_exact


def _cds_only_records(seqid, gene, strand, cds):
    """gene + mRNA + CDS rows, no exon rows (as in CDS-only prediction files)."""
    start, end = min(s for s, _ in cds), max(e for _, e in cds)
    records = [
        (seqid, "gene", start, end, strand, f"ID={gene}"),
        (seqid, "mRNA", start, end, strand, f"ID={gene}.1;Parent={gene}"),
    ]
    return records + [(seqid, "CDS", s, e, strand, f"Parent={gene}.1") for s, e in cds]


def test_cds_only_query_is_compared(compare):
    # Before exons were inferred for CDS-only transcripts, such a query gene had
    # no comparable transcript: its CDS was never compared and every reference
    # gene it covered was reported without a CDS partner.
    cds = [(150, 200), (300, 400), (500, 551)]
    ref, query, summary = compare(
        gene_records("1", "R", "+", {"R.1": MULTI}, with_cds={"R.1": cds})
        + gene_records(
            "1", "S", "-", {"S.1": _shift(MULTI, 2000)}, with_cds={"S.1": []}
        )
        + gene_records(
            "1",
            "T",
            "-",
            {"T.1": _shift(MULTI, 4000)},
            with_cds={"T.1": [(4147, 4200), (4300, 4400), (4500, 4550)]},
        ),
        _cds_only_records("1", "Q", "+", cds)
        # stop codon excluded: 3 bp shorter at the 3' end (minus strand: Start)
        + _cds_only_records("1", "QT", "-", [(4150, 4200), (4300, 4400), (4500, 4550)]),
    )
    q, r = query.loc["Q"], ref.loc["R"]
    assert q.classification_cds == "Exact_Match" and q.cds_coordinate_exact
    assert r.cds_coordinate_exact and r.cds_intron_chain_match
    assert r.cds_overlap == 1.0
    # Inferred exons equal the CDS, so the UTR-inclusive comparison sees no UTR.
    assert not r.exon_coordinate_exact
    t, qt = ref.loc["T"], query.loc["QT"]
    assert t.classification_cds == "Exact_Match" and t.cds_intron_chain_match
    assert not t.cds_coordinate_exact and not qt.cds_coordinate_exact
    assert summary["sensitivity_cds"]["cds_coordinate_exact_count"] == 1
    assert summary["sensitivity_cds"]["cds_exact_match_not_coordinate_exact"] == 1


# 7-8. feature-aware strand mismatch -----------------------------------------


class TestStrandMismatch:
    REF = gene_records(
        "1",
        "R",
        "+",
        {"R.1": [(100, 200), (300, 600)]},
        with_cds={"R.1": [(150, 200), (300, 400)]},  # 401-600 is 3' UTR
    )

    def test_opposite_strand_cds_overlap(self, compare):
        ref, query, summary = compare(
            self.REF, gene_records("1", "Q", "-", {"Q.1": [(150, 200), (300, 400)]})
        )
        assert ref.loc["R"].classification == "Strand_Mismatch"
        assert ref.loc["R"].strand_mismatch_basis == "CDS"
        assert ref.loc["R"].matched_id == "Q"
        assert query.loc["Q"].classification == "Strand_Mismatch"
        assert summary["sensitivity"]["strand_mismatch_cds_count"] == 1

    def test_opposite_strand_utr_only_overlap(self, compare):
        # Query CDS lies only in the reference 3' UTR: gmb-compare said
        # Strand_Mismatch; there is no shared coding sequence.
        ref, query, summary = compare(
            self.REF, gene_records("1", "Q", "-", {"Q.1": [(450, 700)]})
        )
        assert ref.loc["R"].classification == "Missed"
        assert ref.loc["R"].strand_mismatch_basis == "no_feature_overlap"
        assert query.loc["Q"].classification == "Novel"
        assert (
            query.loc["Q"].novel_category
            == "Overlaps_evaluated_reference_opposite_strand"
        )
        sens = summary["sensitivity"]
        assert sens["strand_mismatch_count"] == 0
        assert sens["opposite_strand_no_feature_overlap_count"] == 1

    def test_neighbouring_opposite_strand_gene_span_overlap(self, compare):
        query = gene_records("1", "Q", "-", {"Q.1": [(700, 900)]})
        query[0] = ("1", "gene", 550, 900, "-", "ID=Q;biotype=protein_coding")
        ref, query_results, _ = compare(self.REF, query)
        assert ref.loc["R"].classification == "Missed"
        assert query_results.loc["Q"].classification == "Novel"
        assert query_results.loc["Q"].novel_category == "Novel_no_reference_overlap"

    def test_exon_overlap_when_one_side_has_no_cds(self, compare):
        ref, _, _ = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}, with_cds=False),
            gene_records("1", "Q", "-", {"Q.1": [(350, 450)]}),
        )
        assert ref.loc["R"].classification == "Strand_Mismatch"
        assert ref.loc["R"].strand_mismatch_basis == "exon"

    def test_same_strand_cds_overlap(self, compare):
        ref, _, summary = compare(
            self.REF, gene_records("1", "Q", "+", {"Q.1": [(150, 200), (300, 400)]})
        )
        assert ref.loc["R"].classification != "Strand_Mismatch"
        assert ref.loc["R"].cds_coordinate_exact
        assert summary["sensitivity"]["strand_mismatch_count"] == 0


# 9-10. gene splits and merges ------------------------------------------------

FOUR = [(100, 200), (300, 400), (600, 700), (800, 900)]


class TestSplitMerge:
    def test_one_to_one(self, compare):
        ref, query, summary = compare(
            gene_records("1", "R", "+", {"R.1": MULTI}),
            gene_records("1", "Q", "+", {"Q.1": MULTI}),
        )
        assert (
            ref.loc["R"].counterpart_count == 1 and ref.loc["R"].counterpart_ids == "Q"
        )
        assert query.loc["Q"].counterpart_count == 1
        split_merge = summary["split_merge"]
        assert (split_merge["gene_split_count"], split_merge["gene_merge_count"]) == (
            0,
            0,
        )
        assert split_merge["one_to_one_reference_genes"] == 1

    def test_split(self, compare):
        ref, query, summary = compare(
            gene_records("1", "R", "+", {"R.1": FOUR}),
            gene_records("1", "Q1", "+", {"Q1.1": FOUR[:2]})
            + gene_records("1", "Q2", "+", {"Q2.1": FOUR[2:]}),
        )
        assert ref.loc["R"].counterpart_count == 2
        assert ref.loc["R"].counterpart_ids == "Q1,Q2"
        assert summary["split_merge"]["gene_split_count"] == 1
        assert summary["split_merge"]["gene_split_query_genes"] == 2
        assert summary["split_merge"]["gene_merge_count"] == 0

    def test_merge(self, compare):
        ref, query, summary = compare(
            gene_records("1", "R1", "+", {"R1.1": FOUR[:2]})
            + gene_records("1", "R2", "+", {"R2.1": FOUR[2:]}),
            gene_records("1", "Q", "+", {"Q.1": FOUR}),
        )
        assert query.loc["Q"].counterpart_count == 2
        assert query.loc["Q"].counterpart_ids == "R1,R2"
        assert summary["split_merge"]["gene_merge_count"] == 1
        assert summary["split_merge"]["gene_merge_reference_genes"] == 2
        assert summary["split_merge"]["gene_split_count"] == 0

    def test_non_qualifying_neighbours(self, compare):
        # Q2's span overlaps R through R's 3' UTR only; Q3 shares 5 bp of CDS
        # with R (< 10% of the shorter CDS). Neither is a counterpart.
        ref, query, summary = compare(
            gene_records(
                "1",
                "R",
                "+",
                {"R.1": [(100, 200), (300, 600)]},
                with_cds={"R.1": [(150, 200), (300, 400)]},
            ),
            gene_records("1", "Q1", "+", {"Q1.1": [(150, 200), (300, 400)]})
            + gene_records("1", "Q2", "+", {"Q2.1": [(550, 900)]})
            + gene_records("1", "Q3", "+", {"Q3.1": [(396, 595)]}),
        )
        assert ref.loc["R"].counterpart_ids == "Q1"
        assert summary["split_merge"]["gene_split_count"] == 0
        assert query.loc["Q2"].counterpart_count == 0
        assert query.loc["Q3"].counterpart_count == 0
        # Both are still locus partners for the existing classification.
        assert query.loc["Q2"].classification == "Partial_Match"

    def test_opposite_strand_is_not_a_counterpart(self, compare):
        ref, _, summary = compare(
            gene_records("1", "R", "+", {"R.1": FOUR}),
            gene_records("1", "Q1", "+", {"Q1.1": FOUR[:2]})
            + gene_records("1", "Q2", "-", {"Q2.1": FOUR[2:]}),
        )
        assert ref.loc["R"].counterpart_ids == "Q1"
        assert summary["split_merge"]["gene_split_count"] == 0


# 11-13. Novel context ---------------------------------------------------------


def _pseudogene(seqid, gene, strand, start, end):
    """Ensembl-style pseudogene / pseudogenic_transcript / exon records."""
    return [
        (seqid, "pseudogene", start, end, strand, f"ID=gene:{gene};biotype=pseudogene"),
        (
            seqid,
            "pseudogenic_transcript",
            start,
            end,
            strand,
            f"ID=transcript:{gene}.1;Parent=gene:{gene};biotype=pseudogene",
        ),
        (seqid, "exon", start, end, strand, f"Parent=transcript:{gene}.1"),
    ]


class TestNovelContext:
    @pytest.fixture
    def results(self, compare):
        reference = (
            gene_records("1", "PC", "+", {"PC.1": MULTI})
            + _pseudogene("1", "PSEUDO", "+", 2000, 2600)
            + gene_records(
                "1",
                "LNC",
                "+",
                {"LNC.1": [(4000, 4600)]},
                biotype="lncRNA",
                with_cds=False,
            )
        )
        query = (
            gene_records("1", "Q_PC", "+", {"Q_PC.1": MULTI})
            + gene_records("1", "Q_PSEUDO", "+", {"Q_PSEUDO.1": [(2100, 2500)]})
            + gene_records("1", "Q_LNC", "-", {"Q_LNC.1": [(4100, 4400)]})
            + gene_records("1", "Q_NEW", "+", {"Q_NEW.1": [(8000, 8600)]})
        )
        return compare(reference, query, protein_coding=True)

    def test_categories(self, results):
        ref, query, _ = results
        assert list(ref.index) == ["PC"]  # pseudogene and lncRNA are filtered out
        assert query.loc["Q_PC"].novel_category == ""
        assert query.loc["Q_NEW"].classification == "Novel"
        assert query.loc["Q_NEW"].novel_category == "Novel_no_reference_overlap"
        assert query.loc["Q_PSEUDO"].classification == "Novel"
        assert query.loc["Q_PSEUDO"].novel_category == "Overlaps_reference_pseudogene"
        assert query.loc["Q_LNC"].classification == "Novel"
        assert (
            query.loc["Q_LNC"].novel_category
            == "Overlaps_other_non_protein_coding_reference"
        )

    def test_not_counted_as_matches(self, results):
        _, _, summary = results
        assert summary["specificity"]["matched_consensus_count"] == 1
        assert summary["specificity"]["novel_consensus_count"] == 3
        assert summary["novel_categories"] == {
            "Overlaps_reference_pseudogene": 1,
            "Overlaps_other_non_protein_coding_reference": 1,
            "Overlaps_unevaluated_protein_coding_reference": 0,
            "Overlaps_evaluated_reference_opposite_strand": 0,
            "Novel_no_reference_overlap": 1,
        }


# 14. denominators -------------------------------------------------------------


def test_denominators(compare):
    ref, query, summary = compare(
        # R1 multi-exon, matched exactly; R2 multi-exon, missed; R3 single-exon,
        # matched; R4 single-exon, missed.
        gene_records("1", "R1", "+", {"R1.1": MULTI})
        + gene_records("1", "R2", "+", {"R2.1": _shift(MULTI, 2000)})
        + gene_records("1", "R3", "+", {"R3.1": [(4000, 4600)]})
        + gene_records("1", "R4", "+", {"R4.1": [(6000, 6600)]}),
        gene_records("1", "Q1", "+", {"Q1.1": MULTI})
        + gene_records("1", "Q3", "+", {"Q3.1": [(4000, 4600)]})
        + gene_records("1", "Q5", "+", {"Q5.1": [(9000, 9300), (9400, 9600)]}),
    )
    total_ref, total_query = (
        summary["total_reference_genes"],
        summary["total_consensus_genes"],
    )
    assert (total_ref, total_query) == (4, 3)

    sens, cds = summary["sensitivity"], summary["sensitivity_cds"]
    # Reference-based: divide by R = 4.
    assert sens["locus_detection_rate"] == 0.5
    assert sens["exon_coordinate_exact_rate"] == 0.5
    assert cds["cds_coordinate_exact_rate"] == 0.5
    assert cds["cds_exact_match_rate"] == 0.5
    # Intron chain: only multi-exon reference genes (R1, R2).
    chain = summary["intron_chain"]
    assert chain["multi_exon_reference_genes"] == 2
    assert chain["exon_intron_chain_recovered"] == 1
    assert chain["exon_intron_chain_evaluated"] == 1  # detected, multi-exon: R1
    assert chain["exon_intron_chain_rate"] == 1.0
    assert chain["exon_intron_chain_sensitivity"] == 0.5
    assert cds["cds_intron_chain_sensitivity"] == 0.5

    # Query-based: divide by Q = 3.
    spec = summary["specificity"]
    assert spec["matched_consensus_count"] == 2
    assert spec["matched_consensus_rate"] == round(2 / 3, 4)
    assert spec["novel_consensus_rate"] == round(1 / 3, 4)
    assert spec["consensus_cds_coordinate_exact_rate"] == round(2 / 3, 4)

    # Intron support: reference introns recovered vs query introns supported.
    introns = summary["intron_support"]["exon"]
    assert (introns["reference_introns"], introns["query_introns"]) == (4, 3)
    assert introns["shared_introns"] == 2
    assert introns["reference_introns_recovered_rate"] == 0.5
    assert introns["query_introns_supported_rate"] == round(2 / 3, 4)

    for block in ("sensitivity", "sensitivity_cds", "intron_chain", "specificity"):
        for name, value in summary[block].items():
            if name.endswith(("_rate", "_sensitivity")):
                assert 0.0 <= value <= 1.0, (block, name)
    assert set(summary["rate_denominators"].values()) >= {
        "total_reference_genes",
        "total_consensus_genes (query genes)",
    }


# End-to-end ---------------------------------------------------------------------


def test_runner_writes_refined_outputs(write_gff3, tmp_path):
    reference = write_gff3(
        "ref.gff3",
        gene_records("1", "R", "+", {"R.1": FOUR})
        + _pseudogene("1", "PSEUDO", "+", 3000, 3600),
    )
    query = write_gff3(
        "query.gff3",
        gene_records("1", "Q1", "+", {"Q1.1": FOUR[:2]})
        + gene_records("1", "Q2", "+", {"Q2.1": FOUR[2:]})
        + gene_records("1", "Q3", "+", {"Q3.1": [(3100, 3500)]}),
    )
    outdir = tmp_path / "out"
    summary = _run(
        [
            "--query",
            query,
            "--reference",
            reference,
            "--evaluation-mode",
            "protein_coding",
            "--outdir",
            str(outdir),
        ]
    )
    assert summary["split_merge"]["gene_split_count"] == 1
    assert summary["novel_categories"]["Overlaps_reference_pseudogene"] == 1
    assert summary["intron_support"]["cds"]["reference_introns"] == 3

    with open(outdir / "gene_splits.tsv") as handle:
        splits = list(csv.DictReader(handle, delimiter="\t"))
    assert len(splits) == 1
    assert (splits[0]["gene_id"], splits[0]["counterpart_count"]) == ("R", "2")
    assert splits[0]["counterpart_loci"] == "1:100-400,1:600-900"
    with open(outdir / "gene_merges.tsv") as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []

    with open(outdir / "comparison_details.tsv") as handle:
        rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["Q3"]["novel_category"] == "Overlaps_reference_pseudogene"
    assert rows["R"]["counterpart_ids"] == "Q1,Q2"
    assert rows["Q3"]["intron_chain_match"] == "NA"

    with open(outdir / "consensus_transcript_labels.tsv") as handle:
        labels = {
            row["transcript_id"]: row for row in csv.DictReader(handle, delimiter="\t")
        }
    assert labels["Q3.1"]["classification"] == "Novel"
    assert labels["Q3.1"]["novel_category"] == "Overlaps_reference_pseudogene"

    tsv = (outdir / "comparison_summary.tsv").read_text()
    assert "spec_matched_consensus_count\t2" in tsv
    assert "split_merge_gene_split_count\t1" in tsv
    assert json.loads((outdir / "comparison_summary.json").read_text())[
        "rate_denominators"
    ]

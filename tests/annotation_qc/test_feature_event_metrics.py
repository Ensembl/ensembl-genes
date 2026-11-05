"""Tests for the current feature and genome-sequence event metrics."""

from pathlib import Path

import pandas as pd
import pyranges1 as pr
from pyfaidx import Fasta

from ensembl.genes.annotation_qc.metrics.events import (
    compute_polyA_motif_metrics,
    compute_promoter_proximity_metrics,
    compute_splice_junction_metrics,
    compute_translation_metrics,
    derive_introns,
)
from ensembl.genes.annotation_qc.metrics.features import compute_feature_metrics
from ensembl.genes.annotation_qc.parsers.annotation import parse_annotation
from ensembl.genes.annotation_qc.runners.assess_feature_events import run


def _annotation() -> pr.PyRanges:
    rows = [
        ["chr1", 0, 110, ".", "gene", "gene:g1", "", "g1", "", ""],
        ["chr1", 10, 100, "+", "mRNA", "transcript:tx1", "gene:g1", "g1", "tx1", ""],
        ["chr1", 10, 50, "+", "exon", "ex1", "transcript:tx1", "g1", "tx1", ""],
        ["chr1", 60, 100, "+", "exon", "ex2", "transcript:tx1", "g1", "tx1", ""],
        ["chr1", 20, 32, "+", "CDS", "cds1", "transcript:tx1", "g1", "tx1", "0"],
        ["chr1", 60, 66, "+", "CDS", "cds2", "transcript:tx1", "g1", "tx1", "0"],
        [
            "chr1",
            10,
            20,
            "+",
            "five_prime_UTR",
            "u5",
            "transcript:tx1",
            "g1",
            "tx1",
            "",
        ],
        [
            "chr1",
            66,
            100,
            "+",
            "three_prime_UTR",
            "u3",
            "transcript:tx1",
            "g1",
            "tx1",
            "",
        ],
        ["chr2", 120, 220, ".", "gene", "gene:g2", "", "g2", "", ""],
        ["chr2", 120, 220, "-", "mRNA", "transcript:tx2", "gene:g2", "g2", "tx2", ""],
        ["chr2", 180, 220, "-", "exon", "ex3", "transcript:tx2", "g2", "tx2", ""],
        ["chr2", 120, 160, "-", "exon", "ex4", "transcript:tx2", "g2", "tx2", ""],
    ]
    columns = [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "Feature",
        "ID",
        "Parent",
        "gene_id",
        "transcript_id",
        "phase",
    ]
    return pr.PyRanges(pd.DataFrame(rows, columns=columns))


def _fasta(tmp_path: Path) -> Fasta:
    chr1 = list("C" * 240)
    chr1[14:25] = list("CCGACCATGGG")
    chr1[20:32] = list("ATGGCCCCCCCC")
    chr1[50:60] = list("gtccccccag")
    chr1[60:66] = list("CCCTAA")
    chr1[3:9] = list("TATAAA")
    chr1[82:88] = list("AATAAA")

    chr2 = list("C" * 240)
    chr2[160:162] = list("CT")
    chr2[178:180] = list("AC")
    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(
        ">chr1\n" + "".join(chr1) + "\n>chr2\n" + "".join(chr2) + "\n"
    )
    return Fasta(str(fasta_path))


def test_feature_metrics_are_transcript_rows():
    metrics = compute_feature_metrics(_annotation())
    tx1 = metrics.loc[metrics["feature_id"] == "tx1"].iloc[0]

    assert tx1["exon_count"] == 2
    assert tx1["exonic_length_bp"] == 80
    assert tx1["intron_count"] == 1
    assert tx1["cds_segment_count"] == 2
    assert tx1["cds_length_nt"] == 18
    assert bool(tx1["cds_frame_ok"])


def test_introns_and_splice_junctions_are_transcript_oriented(tmp_path):
    annotation = _annotation()
    introns = derive_introns(annotation)
    assert set(introns["transcript_id"]) == {"tx1", "tx2"}
    assert (
        introns.loc[introns["transcript_id"] == "tx1", "intron_length_bp"].iloc[0] == 10
    )

    junctions = compute_splice_junction_metrics(annotation, _fasta(tmp_path))
    tx1 = junctions.loc[junctions["transcript_id"] == "tx1"].iloc[0]
    tx2 = junctions.loc[junctions["transcript_id"] == "tx2"].iloc[0]
    assert tx1["junction_type"] == "GT-AG"
    assert bool(tx1["is_canonical"])
    assert tx2["junction_type"] == "GT-AG"


def test_translation_metrics_include_kozak_context(tmp_path):
    metrics = compute_translation_metrics(_annotation(), _fasta(tmp_path))
    tx1 = metrics.loc[metrics["transcript_id"] == "tx1"].iloc[0]

    assert tx1["start_codon"] == "ATG"
    assert tx1["stop_codon"] == "TAA"
    assert tx1["start_status"] == "canonical"
    assert tx1["stop_status"] == "canonical"
    assert tx1["internal_stop_count"] == 0
    assert tx1["kozak_window"] == "CCGACCATGG"
    assert tx1["kozak_minus3_base"] == "A"
    assert tx1["kozak_plus4_base"] == "G"
    assert tx1["kozak_context_class"] == "strong"


def test_boundary_motif_metrics_are_configurable_and_conservative(tmp_path):
    annotation = _annotation()
    fasta = _fasta(tmp_path)
    profile = {
        "promoter": {"motifs": {"tata_like": "TATAAA"}, "upstream": 10},
        "polya": {"motifs": {"canonical": "AATAAA"}, "upstream": 50},
    }

    promoter = compute_promoter_proximity_metrics(annotation, fasta, profile)
    tx1_promoter = promoter.loc[promoter["transcript_id"] == "tx1"].iloc[0]
    assert tx1_promoter["nearest_motif"] == "tata_like"
    assert tx1_promoter["nearest_motif_distance_bp"] == 1
    assert tx1_promoter["status"] == "sequence_supported_candidate"

    polya = compute_polyA_motif_metrics(annotation, fasta, profile)
    tx1_polya = polya.loc[polya["transcript_id"] == "tx1"].iloc[0]
    assert tx1_polya["nearest_motif"] == "canonical"
    assert tx1_polya["nearest_motif_distance_bp"] == 12

    unconfigured = compute_promoter_proximity_metrics(annotation, fasta)
    assert set(unconfigured["status"]) == {"not_configured"}


def test_gff3_parser_and_feature_event_runner(tmp_path):
    annotation_path = tmp_path / "annotation.gff3"
    annotation_path.write_text(
        "##gff-version 3\n"
        "chr1\tsource\tgene\t1\t110\t.\t+\t.\tID=g1\n"
        "chr1\tsource\tmRNA\t11\t100\t.\t+\t.\tID=tx1;Parent=g1\n"
        "chr1\tsource\texon\t11\t50\t.\t+\t.\tID=ex1;Parent=tx1\n"
        "chr1\tsource\texon\t61\t100\t.\t+\t.\tID=ex2;Parent=tx1\n"
        "chr1\tsource\tCDS\t21\t32\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\tsource\tCDS\t61\t66\t.\t+\t0\tID=cds2;Parent=tx1\n"
    )
    parsed = parse_annotation(str(annotation_path))
    assert set(parsed.loc[parsed["Feature"] == "CDS", "transcript_id"]) == {"tx1"}

    fasta = _fasta(tmp_path)
    fasta.close()
    fasta_path = tmp_path / "genome.fa"
    output_dir = tmp_path / "results"
    run(
        str(annotation_path),
        str(fasta_path),
        output_dir,
        promoter_motifs=["tata_like=TATAAA"],
        promoter_upstream=10,
    )
    assert (output_dir / "feature_metrics.tsv").is_file()
    assert (output_dir / "translation_events.tsv").is_file()
    assert (output_dir / "polya_events.tsv").is_file()

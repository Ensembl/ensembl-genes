"""
Runner to assess CDS translation validity.

Loads a GFF3/GTF annotation into a PyRanges object, extracts CDS sequences
from a genomic FASTA using pyfaidx, and runs CDS validity QC checks.

Usage:
        python assess_translation_validity.py \
            --annotation_file annotations.gff3 --genome_fasta genome.fa
"""

import argparse
import os

import pandas as pd
from pyfaidx import Fasta

from ensembl.genes.annotation_qc.metrics.events import (
    compute_translation_metrics,
    extract_cds_sequences,
    get_genetic_code,
    STANDARD_CODE,
)
from ensembl.genes.annotation_qc.parsers.annotation import parse_annotation
from ensembl.genes.annotation_qc.parsers.diamond import parse_diamond_hits


def pair_cds_with_sequences(annotation, fasta) -> dict:
    """
    Extract and concatenate CDS intervals from the FASTA for each transcript.

    CDS intervals are sorted by genomic position, with reverse-strand transcripts
    processed in descending order with reverse-complement applied to each interval.

    Args:
            annotation: PyRanges object containing annotation features
            fasta: pyfaidx Fasta object for the genomic sequence
    Returns:
            dict mapping transcript_id to concatenated CDS nucleotide sequence
    """
    return extract_cds_sequences(annotation, fasta)


def run_cds_qc(annotation, fasta, genetic_code: dict = STANDARD_CODE) -> pd.DataFrame:
    """
    Run CDS validity checks for all transcripts in the annotation.
    Args:
            annotation: PyRanges object containing annotation features
            fasta: pyfaidx Fasta object for the genomic sequence
            genetic_code: GeneticCode to use for stop codon evaluation (default: standard)
    Returns:
            DataFrame with one row per transcript and CDSMetrics columns
    """
    return compute_translation_metrics(annotation, fasta, genetic_code)


def _add_args(parser):
    parser.add_argument(
        "--annotation", required=True, help="Path to GFF3 or GTF annotation file."
    )
    parser.add_argument("--genome", required=True, help="Path to genomic FASTA file.")
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory to save results CSV (default: current directory).",
    )
    parser.add_argument(
        "--genetic-code",
        type=int,
        default=1,
        dest="genetic_code",
        metavar="TABLE_ID",
        help=(
            "NCBI genetic code table ID (default: 1 = standard). "
            "See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"
        ),
    )
    parser.add_argument(
        "--diamond-hits",
        default=None,
        dest="diamond_hits",
        help=(
            "Path to DIAMOND tabular hits TSV (optional). "
            "Best hit per transcript is merged into results."
        ),
    )


# CLI integration for shared annotation-qc framework
def register(subparsers):
    """Register this runner as a subcommand on the shared CLI subparsers object."""
    p = subparsers.add_parser(
        "translation-validity",
        help="Assess CDS translation validity (start/stop codons, frame, internal stops).",
    )
    _add_args(p)
    p.set_defaults(func=_run)


def _run(args):
    os.makedirs(args.outdir, exist_ok=True)
    genetic_code = get_genetic_code(args.genetic_code)
    print(f"Using genetic code: {genetic_code['table_id']} ({genetic_code['name']})")
    annotation = parse_annotation(args.annotation)
    fasta = Fasta(args.genome)
    results = run_cds_qc(annotation, fasta, genetic_code)
    if args.diamond_hits:
        diamond_df = parse_diamond_hits(args.diamond_hits)
        best_hits = (
            diamond_df.sort_values("bitscore", ascending=False)
            .drop_duplicates(subset="qseqid")
            .rename(columns={"qseqid": "transcript_id"})
        )
        results = results.merge(best_hits, on="transcript_id", how="left")
    print("Writing results to CSV...")
    results.to_csv(f"{args.outdir}/translation_validity_results.csv", index=False)
    return results


def main():
    """Entry point for standalone execution of the translation validity runner."""
    parser = argparse.ArgumentParser(
        description="Assess the validity of translations in gene annotations."
    )
    _add_args(parser)
    _run(parser.parse_args())


if __name__ == "__main__":
    main()

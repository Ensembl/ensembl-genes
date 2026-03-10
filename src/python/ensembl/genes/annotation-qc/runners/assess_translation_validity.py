"""
Runner to assess CDS translation validity.

Loads a GFF3/GTF annotation into a PyRanges object, extracts CDS sequences
from a genomic FASTA using pyfaidx, and runs CDS validity QC checks.

Usage:
	python assess_translation_validity.py --annotation_file annotations.gff3 --genome_fasta genome.fa
"""

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

import pandas as pd
import os

sys.path.insert(0, str(Path(__file__).parent.parent))

from parsers.annotation import parse_annotation
from parsers.sequence import parse_fasta
from metrics.events import compute_cds_metrics






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
	cds_df = annotation[annotation["Feature"] == "CDS"]

	if "transcript_id" in cds_df.columns and cds_df["transcript_id"].notna().any():
		group_col = "transcript_id"
	elif "Parent" in cds_df.columns:
		group_col = "Parent"
		cds_df = cds_df.assign(Parent=cds_df["Parent"].str.replace(r"^transcript:", "", regex=True))
	else:
		raise ValueError("Annotation has no 'transcript_id' or 'Parent' column — cannot group CDS by transcript.")

	transcript_sequences = {}
	for transcript_id, group in cds_df.groupby(group_col):
		strand = group["Strand"].iloc[0]
		group = group.sort_values("Start", ascending=(strand == "+"))

		seq_parts = []
		for _, row in group.iterrows():
			chrom = row["Chromosome"]
			start = int(row["Start"])
			end = int(row["End"])
			try:
				interval = fasta[chrom][start:end]
				seq_parts.append(
					interval.seq if strand == "+" else interval.reverse.complement.seq
				)
			except KeyError:
				seq_parts.append("")

		transcript_sequences[transcript_id] = "".join(seq_parts)

	return transcript_sequences


def run_cds_qc(annotation, fasta) -> pd.DataFrame:
	"""
	Run CDS validity checks for all transcripts in the annotation.
	Args:
		annotation: PyRanges object containing annotation features
		fasta: pyfaidx Fasta object for the genomic sequence
	Returns:
		DataFrame with one row per transcript and CDSMetrics columns
	"""
	cds_sequences = pair_cds_with_sequences(annotation, fasta)
	records = [
		{"transcript_id": tid, **asdict(compute_cds_metrics(seq))}
		for tid, seq in cds_sequences.items()
	]
	return pd.DataFrame(records)

# CLI integration for shared annotation-qc framework
def register(subparsers):
	"""Register this runner as a subcommand on the shared CLI subparsers object."""
	p = subparsers.add_parser(
		"translation-validity",
		help="Assess CDS translation validity (start/stop codons, frame, internal stops).",
	)
	p.add_argument("--annotation", required=True, help="Path to GFF3 or GTF annotation file.")
	p.add_argument("--genome", required=True, help="Path to genomic FASTA file.")
	p.add_argument("--outdir", default=".", help="Directory to save results CSV (default: current directory).")
	p.set_defaults(func=_run)


def _run(args):
	os.makedirs(args.outdir, exist_ok=True)
	annotation = parse_annotation(args.annotation)
	fasta = parse_fasta(args.genome)
	results = run_cds_qc(annotation, fasta)
	print("Writing results to CSV...")
	results.to_csv(f"{args.outdir}/translation_validity_results.csv", index=False)
	return results

# Standalone CLI entry point for direct execution
def parse_args():
	parser = argparse.ArgumentParser(description="Assess the validity of translations in gene annotations.")
	parser.add_argument("--annotation", required=True, help="Path to the gene annotation file (GFF3 or GTF).")
	parser.add_argument("--genome", required=True, help="Path to the genome FASTA file.")
	return parser.parse_args()


def main():
	args = parse_args()
	annotation = parse_annotation(args.annotation)
	fasta = parse_fasta(args.genome)
	results = run_cds_qc(annotation, fasta)
	print("Writing results to CSV...")
	results.to_csv(f"{args.outdir}/translation_validity_results.csv", index=False)
	return results


if __name__ == "__main__":
	main()

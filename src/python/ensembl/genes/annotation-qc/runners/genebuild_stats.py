"""
Runner to compute genebuild feature statistics from a GFF3/GTF annotation.

Calculates summary statistics for coding genes, non-coding genes, and
pseudogenes (gene counts, spans, transcript counts, exon/intron lengths,
etc.) and per-transcript 5'/3' UTR lengths and junction counts.

When --genome is supplied, also runs sequence-based event QC (no alignment
evidence required): per-transcript CDS translation metrics and per-intron
splice junction canonicality.

Usage:
	python genebuild_stats.py --annotation annotations.gff3 --outdir results/
	python genebuild_stats.py --annotation annotations.gff3 --genome genome.fa --outdir results/
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from parsers.annotation import parse_annotation
from metrics.features import (
	compute_coding_stats,
	compute_noncoding_stats,
	compute_pseudogene_stats,
	compute_transcript_utr_stats,
	compute_cds_utr5_overlap,
	compute_cds_utr5_overlap_summary,
	compute_translation_metrics,
	compute_translation_summary_stats,
	compute_splice_junction_metrics,
	compute_splice_junction_summary_stats,
)


def _to_long(stats: dict, category: str) -> pd.DataFrame:
	"""Pivot a flat stats dict to long format: category | metric | value."""
	return pd.DataFrame(
		[{"category": category, "metric": k, "value": v} for k, v in stats.items()]
	)


def _run(args):
	os.makedirs(args.outdir, exist_ok=True)

	print(f"Parsing annotation: {args.annotation}")
	gff = parse_annotation(args.annotation)

	print("Computing coding gene stats...")
	coding_stats, total_coding_seq_len = compute_coding_stats(gff)
	coding_stats["Total coding sequence length"] = total_coding_seq_len

	print("Computing non-coding gene stats...")
	noncoding_stats = compute_noncoding_stats(gff)

	print("Computing pseudogene stats...")
	pseudogene_stats = compute_pseudogene_stats(gff)

	print("Computing per-transcript UTR stats...")
	utr_df = compute_transcript_utr_stats(gff)

	print("Computing CDS / 5' UTR exon overlaps...")
	cds_utr5_df = compute_cds_utr5_overlap(gff)

	summary_parts = [
		_to_long(coding_stats, "coding"),
		_to_long(noncoding_stats, "noncoding"),
		_to_long(pseudogene_stats, "pseudogene"),
		_to_long(compute_cds_utr5_overlap_summary(cds_utr5_df), "cds_utr5_overlap"),
	]

	utr_path = os.path.join(args.outdir, "utr_stats.tsv")
	cds_utr5_path = os.path.join(args.outdir, "cds_utr5_overlap.tsv")

	utr_df.to_csv(utr_path, index=False, sep="\t")
	cds_utr5_df.to_csv(cds_utr5_path, index=False, sep="\t")

	written = [utr_path, cds_utr5_path]

	if args.genome:
		from parsers.sequence import parse_fasta
		from metrics.events import get_genetic_code

		fasta = parse_fasta(args.genome)
		genetic_code = get_genetic_code(args.genetic_code)
		print(f"Using genetic code: {genetic_code.table_id} ({genetic_code.name})")

		print("Computing translation metrics...")
		translation_df = compute_translation_metrics(gff, fasta, genetic_code)
		translation_path = os.path.join(args.outdir, "translation_metrics.tsv")
		translation_df.to_csv(translation_path, index=False, sep="\t")
		written.append(translation_path)
		summary_parts.append(_to_long(compute_translation_summary_stats(translation_df), "translation"))

		print("Computing splice junction metrics...")
		splice_df = compute_splice_junction_metrics(gff, fasta)
		splice_path = os.path.join(args.outdir, "splice_junction_metrics.tsv")
		splice_df.to_csv(splice_path, index=False, sep="\t")
		written.append(splice_path)
		summary_parts.append(_to_long(compute_splice_junction_summary_stats(splice_df), "splicing"))

	summary_path = os.path.join(args.outdir, "feature_stats.tsv")
	pd.concat(summary_parts, ignore_index=True).to_csv(summary_path, index=False, sep="\t")
	written.insert(0, summary_path)

	print("Done.\n" + "\n".join(f"  {p}" for p in written))


def _add_args(parser):
	parser.add_argument(
		"--annotation",
		required=True,
		help="Path to GFF3 or GTF annotation file (plain or .gz).",
	)
	parser.add_argument(
		"--genome",
		default=None,
		help=(
			"Path to genomic FASTA file. When provided, also runs per-transcript "
			"translation QC and per-intron splice junction assessment "
			"(no alignment evidence required)."
		),
	)
	parser.add_argument(
		"--genetic-code",
		type=int,
		default=1,
		dest="genetic_code",
		metavar="TABLE_ID",
		help="NCBI genetic code table ID used for translation QC (default: 1 = standard).",
	)
	parser.add_argument(
		"--outdir",
		default=".",
		help="Directory to write output TSVs (default: current directory).",
	)


def register(subparsers):
	"""Register this runner as a subcommand on the shared CLI subparsers object."""
	p = subparsers.add_parser(
		"feature-stats",
		help=(
			"Compute genebuild feature statistics: coding/non-coding/pseudogene "
			"summaries, per-transcript UTR metrics, and (with --genome) "
			"translation and splice junction QC."
		),
	)
	_add_args(p)
	p.set_defaults(func=_run)


def main():
	parser = argparse.ArgumentParser(
		description="Compute genebuild feature statistics from a GFF3/GTF annotation."
	)
	_add_args(parser)
	_run(parser.parse_args())


if __name__ == "__main__":
	main()

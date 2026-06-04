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

import pandas as pd

from ensembl.genes.annotation_qc.parsers.annotation import (
    CDS_FEATURE_TYPES,
    EXON_FEATURE_TYPES,
    FIVE_PRIME_UTR_FEATURE_TYPES,
    GENE_FEATURE_TYPES,
    THREE_PRIME_UTR_FEATURE_TYPES,
    TRANSCRIPT_FEATURE_TYPES,
    parse_annotation,
)
from ensembl.genes.annotation_qc.metrics.features import (
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
from ensembl.genes.annotation_qc.metadata.ncbi import ncbi_assembly_stats
from ensembl.genes.annotation_qc.reports.meta_sql import (
    write_json_from_tsv,
    write_meta_sql,
    write_tsv,
)


def _to_long(stats: dict, category: str) -> pd.DataFrame:
    """Pivot a flat stats dict to long format: category | metric | value."""
    return pd.DataFrame(
        [{"category": category, "metric": k, "value": v} for k, v in stats.items()]
    )


def _run(args):
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Parsing annotation: {args.annotation}")
    annotation = parse_annotation(args.annotation)
    gene_df = annotation[annotation["Feature"].isin(GENE_FEATURE_TYPES)].reset_index(
        drop=True
    )
    tx_df = annotation[
        annotation["Feature"].isin(TRANSCRIPT_FEATURE_TYPES)
    ].reset_index(drop=True)
    exon_df = annotation[annotation["Feature"].isin(EXON_FEATURE_TYPES)].reset_index(
        drop=True
    )
    cds_df = annotation[annotation["Feature"].isin(CDS_FEATURE_TYPES)].reset_index(
        drop=True
    )
    utr5_df = annotation[
        annotation["Feature"].isin(FIVE_PRIME_UTR_FEATURE_TYPES)
    ].reset_index(drop=True)
    utr3_df = annotation[
        annotation["Feature"].isin(THREE_PRIME_UTR_FEATURE_TYPES)
    ].reset_index(drop=True)

    print("Computing coding gene stats...")
    coding_stats, total_coding_seq_len = compute_coding_stats(
        gene_df, tx_df, exon_df, cds_df
    )
    coding_stats["Total coding sequence length"] = total_coding_seq_len

    print("Computing non-coding gene stats...")
    noncoding_stats = compute_noncoding_stats(gene_df, tx_df, exon_df, cds_df)

    print("Computing pseudogene stats...")
    pseudogene_stats = compute_pseudogene_stats(gene_df, tx_df, exon_df, cds_df)

    assembly_row = None
    if args.assembly_accession:
        asm = ncbi_assembly_stats(args.assembly_accession, outdir=args.outdir)
        scientific_name = str(args.scientific_name or asm.get("scientific_name") or "")
        taxon_id = str(args.taxon_id or asm.get("taxon_id") or "")
        assembly_row = {
            "Scientific name": scientific_name,
            "Sex": args.sex or "",
            "Breed/Cultivar/Strain": args.strain or "",
            "Taxonomy id": taxon_id,
            "Assembly name": asm.get("assembly_name") or "",
            "Assembly accession": asm.get("assembly_accession")
            or args.assembly_accession,
            "Assembly date": asm.get("assembly_date") or "",
            "Contig N50": asm.get("contig_n50") or "",
            "Total genome length": asm.get("total_length") or "",
            "Total coding sequence length": total_coding_seq_len,
            "Total gap length": asm.get("total_gap_length") or "",
            "Spanned gaps": asm.get("spanned_gaps") or "",
            "Chromosomes": asm.get("molecule_count") or "",
            "Toplevel sequences": asm.get("toplevel_sequences") or "",
            "Component sequences": asm.get("component_count") or "",
            "% GC": asm.get("gc_percent") or "",
        }
        for stats in (coding_stats, noncoding_stats, pseudogene_stats):
            stats["Scientific name"] = scientific_name

    print("Computing per-transcript UTR stats...")
    utr_df = compute_transcript_utr_stats(tx_df, gene_df, utr5_df, utr3_df)

    print("Computing CDS / 5' UTR exon overlaps...")
    cds_utr5_df = compute_cds_utr5_overlap(cds_df, utr5_df)

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

    if assembly_row:
        compat_outputs = _write_compat_reports(
            args,
            coding_stats,
            noncoding_stats,
            pseudogene_stats,
            assembly_row,
        )
        written.extend(compat_outputs)

    if args.genome:
        from ensembl.genes.annotation_qc.parsers.sequence import parse_fasta
        from ensembl.genes.annotation_qc.metrics.events import get_genetic_code

        fasta = parse_fasta(args.genome)
        genetic_code = get_genetic_code(args.genetic_code)
        print(f"Using genetic code: {genetic_code.table_id} ({genetic_code.name})")

        print("Computing translation metrics...")
        translation_df = compute_translation_metrics(cds_df, fasta, genetic_code)
        translation_path = os.path.join(args.outdir, "translation_metrics.tsv")
        translation_df.to_csv(translation_path, index=False, sep="\t")
        written.append(translation_path)
        summary_parts.append(
            _to_long(compute_translation_summary_stats(translation_df), "translation")
        )

        print("Computing splice junction metrics...")
        splice_df = compute_splice_junction_metrics(exon_df, fasta)
        splice_path = os.path.join(args.outdir, "splice_junction_metrics.tsv")
        splice_df.to_csv(splice_path, index=False, sep="\t")
        written.append(splice_path)
        summary_parts.append(
            _to_long(compute_splice_junction_summary_stats(splice_df), "splicing")
        )

    summary_path = os.path.join(args.outdir, "feature_stats.tsv")
    pd.concat(summary_parts, ignore_index=True).to_csv(
        summary_path, index=False, sep="\t"
    )
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
    parser.add_argument(
        "--assembly-accession",
        "--assembly_accession",
        dest="assembly_accession",
        default="",
        help="Optional GCA_/GCF_ accession for NCBI assembly stats.",
    )
    parser.add_argument(
        "--scientific-name",
        "--scientific_name",
        dest="scientific_name",
        default="",
        help="Override scientific name for compatibility reports.",
    )
    parser.add_argument(
        "--taxon-id",
        "--taxon_id",
        dest="taxon_id",
        default="",
        help="Override taxonomy ID for compatibility reports.",
    )
    parser.add_argument("--strain", default="", help="Strain/Breed/Cultivar.")
    parser.add_argument("--sex", default="", help="Sex metadata override.")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Also write JSON versions of compatibility TSV outputs.",
    )
    parser.add_argument(
        "--sql",
        action="store_true",
        help="Also write Ensembl meta table INSERT IGNORE SQL.",
    )
    parser.add_argument("--dbname", default="", help="Database name for SQL output.")
    parser.add_argument(
        "--species-id",
        "--species_id",
        dest="species_id",
        type=int,
        default=1,
        help="Species ID for meta inserts.",
    )
    parser.add_argument(
        "--sql-filename",
        "--sql_filename",
        dest="sql_filename",
        default="",
        help="Optional SQL output path.",
    )


def _write_compat_reports(
    args,
    coding_stats: dict,
    noncoding_stats: dict,
    pseudogene_stats: dict,
    assembly_row: dict,
) -> list[str]:
    """Write the standalone script's TSV/JSON/SQL compatibility outputs."""
    paths = []
    coding_headers = [
        "Scientific name",
        "Coding genes",
        "Average genomic span",
        "Average sequence length",
        "Average CDS length",
        "Shortest gene",
        "Longest gene",
        "Total transcripts",
        "Coding transcripts",
        "Transcripts per gene",
        "Coding transcripts per gene",
        "Total exons",
        "Total coding exons",
        "Average exon length",
        "Average coding exon length",
        "Average exons per transcript",
        "Average coding exons per coding transcript",
        "Total introns",
        "Average intron length",
        "Total coding sequence length",
    ]
    noncoding_headers = [
        "Scientific name",
        "Non-coding genes",
        "Small non-coding genes",
        "Long non-coding genes",
        "Misc non-coding genes",
        "Average genomic span",
        "Average sequence length",
        "Shortest gene",
        "Longest gene",
        "Total transcripts",
        "Transcripts per gene",
        "Total exons",
        "Average exon length",
        "Average exons per transcript",
        "Total introns",
        "Average intron length",
    ]
    pseudogene_headers = [
        "Scientific name",
        "Pseudogenes",
        "Average genomic span",
        "Average sequence length",
        "Shortest gene",
        "Longest gene",
        "Total transcripts",
        "Transcripts per gene",
        "Total exons",
        "Average exon length",
        "Average exons per transcript",
        "Total introns",
        "Average intron length",
    ]
    assembly_headers = [
        "Scientific name",
        "Sex",
        "Breed/Cultivar/Strain",
        "Taxonomy id",
        "Assembly name",
        "Assembly accession",
        "Assembly date",
        "Contig N50",
        "Total genome length",
        "Total coding sequence length",
        "Total gap length",
        "Spanned gaps",
        "Chromosomes",
        "Toplevel sequences",
        "Component sequences",
        "% GC",
    ]
    for name, headers, row in (
        ("coding_stats", coding_headers, coding_stats),
        ("noncoding_stats", noncoding_headers, noncoding_stats),
        ("pseudogene_stats", pseudogene_headers, pseudogene_stats),
        ("assembly_stats", assembly_headers, assembly_row),
    ):
        tsv_path = os.path.join(args.outdir, f"{name}.tsv")
        write_tsv(tsv_path, headers, row)
        paths.append(tsv_path)
        if args.json:
            json_path = os.path.join(args.outdir, f"{name}.json")
            write_json_from_tsv(tsv_path, json_path)
            paths.append(json_path)
    if args.sql:
        if not args.species_id:
            raise SystemExit("--species-id is required when using --sql")
        sql_path = args.sql_filename or os.path.join(
            args.outdir, f"stats_{args.assembly_accession}.sql"
        )
        write_meta_sql(
            out_path=sql_path,
            dbname=args.dbname,
            species_id=str(args.species_id),
            coding_stats=coding_stats,
            noncoding_stats=noncoding_stats,
            pseudogene_stats=pseudogene_stats,
            assembly_row=assembly_row,
        )
        paths.append(sql_path)
    return paths


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

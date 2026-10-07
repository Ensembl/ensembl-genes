"""
Write pairwise-comparison results.

Output files (names and columns kept from gmb-compare; "consensus" means the
query annotation):

    comparison_summary.json / .tsv    headline metrics
    comparison_details.tsv            one row per reference and query gene
    consensus_transcript_labels.tsv   one row per query transcript
    gene_splits.tsv                   reference genes with >= 2 query counterparts
    gene_merges.tsv                   query genes with >= 2 reference counterparts
    cds_exact_one_to_one_pairs.tsv    gene pairs of the one-to-one coordinate-exact
                                      CDS matching (one-based coordinates)
    stop_codon_harmonisation.tsv      only with --add-stop-codon: one row per query
                                      transcript with CDS, status and reason
    query_evaluated.stop_harmonised.gtf
                                      only with --add-stop-codon: the query exactly
                                      as compared (after harmonisation, transcript
                                      selection and regions), one-based GTF
    reference_filter_audit.json / .tsv
    evidence_attribution_labeled.tsv  only with --evidence-attribution
    comparison_manifest.json          inputs, checksums, options and versions

Gene coordinates in comparison_details.tsv are one-based inclusive, as in GFF.
Columns after original_gene_id were added after gmb-compare; intron-chain columns
are NA where the compared structure has no intron (single exon / CDS segment).
"""

import csv
import json
import os

import pandas as pd

HARMONISED_QUERY_GTF = "query_evaluated.stop_harmonised.gtf"
STOP_AUDIT_TSV = "stop_codon_harmonisation.tsv"

DETAIL_COLUMNS = [
    "source",
    "gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "classification",
    "classification_cds",
    "matched_id",
    "best_match_transcript_id",
    "exon_overlap",
    "intron_chain_match",
    "cds_overlap",
    "cds_intron_chain_match",
    "match_basis",
    "original_gene_id",
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
SPLIT_MERGE_COLUMNS = [
    "gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "classification",
    "counterpart_count",
    "counterpart_ids",
    "counterpart_loci",
]


def _value(value):
    """Summary TSV value: NA for an undefined (None) rate."""
    return "NA" if value is None else value


def _flag(value) -> str | bool:
    """True/False, or NA for a not-applicable (missing) value."""
    return "NA" if pd.isna(value) else bool(value)


def _write_json(data: dict, path: str) -> None:
    with open(path, "w") as handle:
        json.dump(data, handle, indent=2)


def write_summary(summary: dict, outdir: str) -> None:
    """
    Write comparison_summary.json and comparison_summary.tsv.
    Args:
            summary: Output of summarise_comparison (plus filter_info etc.)
            outdir: Output directory
    """
    _write_json(summary, os.path.join(outdir, "comparison_summary.json"))
    with open(
        os.path.join(outdir, "comparison_summary.tsv"), "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["total_reference_genes", summary["total_reference_genes"]])
        writer.writerow(["total_consensus_genes", summary["total_consensus_genes"]])
        for cls, count in sorted(summary["reference_classification"].items()):
            writer.writerow([f"ref_{cls}", count])
        for cls, count in sorted(summary["consensus_classification"].items()):
            writer.writerow([f"cons_{cls}", count])
        for key in ("sensitivity", "sensitivity_cds"):
            for name, value in summary[key].items():
                writer.writerow([f"sens_{name}", value])
        for name, value in summary["specificity"].items():
            writer.writerow([f"spec_{name}", value])
        for name, value in summary["locus_detection_exonic"].items():
            writer.writerow([f"exonic_{name}", value])
        for name, value in summary["intron_chain"].items():
            writer.writerow([f"intron_chain_{name}", value])
        for name, value in summary["novel_categories"].items():
            writer.writerow([f"novel_{name}", value])
        for name, value in summary["split_merge"].items():
            writer.writerow([f"split_merge_{name}", value])
        for kind, values in summary.get("intron_support", {}).items():
            for name, value in values.items():
                writer.writerow([f"introns_{kind}_{name}", value])
        for block, prefix in (
            ("cds_overlap_locus_recovery", "cds_overlap_locus_recovery"),
            ("cds_exact_one_to_one", "cds_exact_one_to_one"),
            ("query_stop_codon_harmonisation", "query_stop_codon_harmonisation"),
        ):
            for name, value in summary.get(block, {}).items():
                if isinstance(value, dict):
                    if name == "reasons":
                        for reason, count in value.items():
                            writer.writerow([f"{prefix}_{reason}", count])
                elif name != "definition":
                    writer.writerow([f"{prefix}_{name}", _value(value)])


def write_stop_codon_audit(audit: pd.DataFrame, outdir: str) -> str:
    """
    Write stop_codon_harmonisation.tsv (terminal CDS coordinates one-based).
    Args:
            audit: Audit table from stop_codons.harmonise_stop_codons
            outdir: Output directory
    Returns:
            Path of the written file
    """
    out = audit.copy()
    for column in ("terminal_cds_start", "new_terminal_cds_start"):
        out[column] = out[column] + 1
    path = os.path.join(outdir, STOP_AUDIT_TSV)
    out.to_csv(path, sep="\t", index=False)
    return path


def write_annotation_gtf(df: pd.DataFrame, outdir: str, filename: str) -> str:
    """
    Write a comparison-schema annotation as GTF (one-based, file IDs).

    Gene, transcript, exon and CDS rows are written in table order with the IDs
    as they appear in the source file (original_gene_id / original_transcript_id).
    Records the comparator inferred carry annotation_qc_record "<source_feature>".
    Args:
            df: Comparison-schema DataFrame
            outdir: Output directory
            filename: Output file name
    Returns:
            Path of the written file
    """
    path = os.path.join(outdir, filename)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("##gtf written by annotation-qc pairwise-compare\n")
        for row in df.itertuples(index=False):
            attributes = f'gene_id "{row.original_gene_id}";'
            if row.Feature != "gene":
                attributes += f' transcript_id "{row.original_transcript_id}";'
            if row.gene_biotype:
                attributes += f' gene_biotype "{row.gene_biotype}";'
            if str(row.source_feature).startswith("inferred_"):
                attributes += f' annotation_qc_record "{row.source_feature}";'
            handle.write(
                f"{row.Chromosome}\tannotation-qc\t{row.Feature}\t{int(row.Start) + 1}"
                f"\t{int(row.End)}\t.\t{row.Strand}\t.\t{attributes}\n"
            )
    return path


def write_exact_cds_pairs(pairs: pd.DataFrame, outdir: str) -> None:
    """
    Write cds_exact_one_to_one_pairs.tsv: one row per matched gene pair.
    Gene and CDS coordinates are converted to one-based inclusive.
    Args:
            pairs: Pair table from matching.one_to_one_exact_cds
            outdir: Output directory
    """
    out = pairs.copy()
    for column in ("reference_start", "query_start", "cds_start"):
        out[column] = out[column] + 1
    out.to_csv(
        os.path.join(outdir, "cds_exact_one_to_one_pairs.tsv"), sep="\t", index=False
    )


def write_details(ref: pd.DataFrame, query: pd.DataFrame, outdir: str) -> None:
    """
    Write comparison_details.tsv (reference genes first, then query genes).
    Args:
            ref, query: Results from classify_loci (zero-based coordinates)
            outdir: Output directory
    """
    with open(
        os.path.join(outdir, "comparison_details.tsv"), "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(DETAIL_COLUMNS)
        for source, results in (("reference", ref), ("consensus", query)):
            for row in results.itertuples(index=False):
                writer.writerow(
                    [
                        source,
                        row.gene_id,
                        row.chrom,
                        int(row.start) + 1,
                        int(row.end),
                        row.strand,
                        row.classification,
                        row.classification_cds,
                        row.matched_id,
                        row.best_match_transcript_id,
                        round(float(row.exon_overlap), 4),
                        _flag(row.intron_chain_match),
                        round(float(row.cds_overlap), 4),
                        _flag(row.cds_intron_chain_match),
                        "exon",
                        row.original_gene_id,
                        bool(row.exon_coordinate_exact),
                        bool(row.cds_coordinate_exact),
                        row.cds_matched_id,
                        row.best_cds_match_transcript_id,
                        row.best_match_own_transcript_id,
                        row.best_cds_match_own_transcript_id,
                        row.strand_mismatch_basis,
                        row.novel_category,
                        int(row.counterpart_count),
                        row.counterpart_ids,
                    ]
                )


def write_split_merge(ref: pd.DataFrame, query: pd.DataFrame, outdir: str) -> None:
    """
    Write gene_splits.tsv (reference genes with >= 2 query counterparts) and
    gene_merges.tsv (query genes with >= 2 reference counterparts).
    Args:
            ref, query: Results from classify_loci (zero-based coordinates)
            outdir: Output directory
    """
    for name, own, other in (
        ("gene_splits.tsv", ref, query),
        ("gene_merges.tsv", query, ref),
    ):
        loci = {
            row.gene_id: f"{row.chrom}:{int(row.start) + 1}-{int(row.end)}"
            for row in other.itertuples(index=False)
        }
        with open(os.path.join(outdir, name), "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(SPLIT_MERGE_COLUMNS)
            for row in own[own["counterpart_count"] >= 2].itertuples(index=False):
                ids = row.counterpart_ids.split(",")
                writer.writerow(
                    [
                        row.gene_id,
                        row.chrom,
                        int(row.start) + 1,
                        int(row.end),
                        row.strand,
                        row.classification,
                        int(row.counterpart_count),
                        row.counterpart_ids,
                        ",".join(loci[gid] for gid in ids),
                    ]
                )


def write_transcript_labels(labels: pd.DataFrame, outdir: str) -> None:
    """
    Write consensus_transcript_labels.tsv.
    Args:
            labels: Output of query_transcript_labels
            outdir: Output directory
    """
    with open(
        os.path.join(outdir, "consensus_transcript_labels.tsv"), "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(labels.columns.tolist())
        for row in labels.itertuples(index=False):
            writer.writerow(
                [
                    *row[:5],
                    round(float(row.best_overlap), 4),
                    round(float(row.best_cds_overlap), 4),
                    *row[7:],
                ]
            )


def write_filter_audit(audit: dict, outdir: str) -> None:
    """
    Write reference_filter_audit.json and reference_filter_audit.tsv.
    Args:
            audit: Filter settings and counts assembled by the runner
            outdir: Output directory
    """
    _write_json(audit, os.path.join(outdir, "reference_filter_audit.json"))
    with open(
        os.path.join(outdir, "reference_filter_audit.tsv"), "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "key", "value"])
        for key in (
            "evaluation_mode",
            "transcript_selection",
            "gene_biotypes_filter",
            "transcript_biotypes_filter",
        ):
            writer.writerow(["config", key, audit[key] or "none"])
        for key, value in audit.items():
            if isinstance(value, int) and not isinstance(value, bool):
                writer.writerow(["counts", key, value])
        for key in ("gene_biotype_counts", "transcript_biotype_counts"):
            for biotype, count in sorted(audit[key].items(), key=lambda item: -item[1]):
                writer.writerow([key.removesuffix("_counts"), biotype, count])
        for seqname, genes in audit.get("excluded_seqnames", {}).items():
            writer.writerow(["excluded_seqname", seqname, genes])


def write_evidence_attribution_labels(
    attribution: pd.DataFrame, labels: pd.DataFrame, outdir: str
) -> int:
    """
    Add each transcript's comparison label to a GMB evidence attribution table.
    Args:
            attribution: Evidence attribution table (transcript_id column)
            labels: Output of query_transcript_labels
            outdir: Output directory
    Returns:
            Number of attribution rows without a comparison label
    """
    merged = attribution.merge(
        labels[["transcript_id", "classification"]].rename(
            columns={"classification": "comparison_label"}
        ),
        on="transcript_id",
        how="left",
    )
    merged.to_csv(
        os.path.join(outdir, "evidence_attribution_labeled.tsv"), sep="\t", index=False
    )
    return int(merged["comparison_label"].isna().sum())


def write_manifest(manifest: dict, outdir: str) -> None:
    """
    Write comparison_manifest.json.
    Args:
            manifest: Provenance assembled by the runner
            outdir: Output directory
    """
    _write_json(manifest, os.path.join(outdir, "comparison_manifest.json"))


def _percent(value) -> str:
    """A rate as a percentage, or n/a when it is undefined (None)."""
    return "n/a" if value is None else f"{value:.1%}"


def _stop_codon_lines(summary: dict) -> list[str]:
    stop = summary.get("query_stop_codon_harmonisation", {})
    if not stop.get("applied"):
        return []
    return [
        "Query stop codons added (--add-stop-codon; harmonised coordinates "
        f"compared): {stop['extended']} extended, {stop['unchanged']} unchanged, "
        f"{stop['skipped']} skipped of {stop['transcripts_with_cds']}"
    ]


def format_console_summary(summary: dict) -> str:
    """
    Render the headline metrics for the terminal.
    Args:
            summary: Output of summarise_comparison
    Returns:
            Multi-line string
    """
    total = summary["total_reference_genes"]
    total_query = summary["total_consensus_genes"]
    sens, cds = summary["sensitivity"], summary["sensitivity_cds"]
    chain, exonic = summary["intron_chain"], summary["locus_detection_exonic"]
    spec, split_merge = summary["specificity"], summary["split_merge"]

    def pct(count: int, of: int = total) -> str:
        return f"{count / of:.1%}" if of else "n/a"

    def of(count: int, denominator: int, label: str) -> str:
        return f"{count:>6}  of {denominator} {label} ({pct(count, denominator)})"

    cds_chain = cds["cds_intron_chain_recovered"]
    recovery = summary.get("cds_overlap_locus_recovery", {})
    one_to_one = summary.get("cds_exact_one_to_one", {})

    lines = [
        f"Reference genes (R): {total}    Query genes (Q): {total_query}",
        "Reference-based (denominator R unless stated):",
        f"  CDS coordinate-exact:   {cds['cds_coordinate_exact_count']:>6}  ({pct(cds['cds_coordinate_exact_count'])})",
        f"  CDS exact (structural): {cds['cds_exact_match_count']:>6}  ({pct(cds['cds_exact_match_count'])})",
        f"  CDS any match:          {cds['cds_any_match_count']:>6}  ({pct(cds['cds_any_match_count'])})",
        f"  CDS intron chain:       {of(cds_chain, cds['multi_segment_cds_reference_genes'], 'multi-CDS-segment R')}",
        f"  Exon coordinate-exact:  {sens['exon_coordinate_exact_count']:>6}  ({pct(sens['exon_coordinate_exact_count'])})",
        f"  Exact match (with UTR): {sens['exact_match_count']:>6}  ({pct(sens['exact_match_count'])})",
        f"  Intron chain:           {of(chain['exon_intron_chain_recovered'], chain['multi_exon_reference_genes'], 'multi-exon R')}",
        f"  Locus, CDS overlap:     {recovery.get('recovered_count', 0):>6}  ({_percent(recovery.get('rate'))})",
        f"  Locus detected (span):  {sens['locus_detected_count']:>6}  ({pct(sens['locus_detected_count'])})",
        f"    with exonic overlap:  {exonic['locus_detected_exonic_count']:>6}  ({pct(exonic['locus_detected_exonic_count'])})",
        f"  Missed:                 {sens['missed_count']:>6}  ({pct(sens['missed_count'])})",
        f"  Strand mismatch:        {sens['strand_mismatch_count']:>6}  (CDS {sens['strand_mismatch_cds_count']}, exon {sens['strand_mismatch_exon_count']})",
        f"  Gene splits:            {split_merge['gene_split_count']:>6}",
        "Query-based (denominator Q):",
        f"  Matched query genes:    {spec['matched_consensus_count']:>6}  ({pct(spec['matched_consensus_count'], total_query)})",
        f"  CDS coordinate-exact:   {spec['consensus_cds_coordinate_exact_count']:>6}  ({pct(spec['consensus_cds_coordinate_exact_count'], total_query)})",
        f"  Novel query genes:      {spec['novel_consensus_count']:>6}  ({pct(spec['novel_consensus_count'], total_query)})",
    ]
    lines += [
        f"    {category}: {count}"
        for category, count in summary["novel_categories"].items()
        if count
    ]
    lines.append(f"  Gene merges:            {split_merge['gene_merge_count']:>6}")
    lines += _stop_codon_lines(summary)
    if one_to_one:
        lines += [
            "One-to-one coordinate-exact CDS (maximum matching):",
            f"  TP: {one_to_one['true_positives']}   recall (TP/R): {_percent(one_to_one['recall'])}"
            f"   precision (TP/Q): {_percent(one_to_one['precision'])}"
            f"   F1: {'n/a' if one_to_one['f1'] is None else format(one_to_one['f1'], '.4f')}",
        ]
    return "\n".join(lines)

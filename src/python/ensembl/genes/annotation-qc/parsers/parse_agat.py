"""
Select and rename AGAT metrics to Ensembl Genebuild stat names.

Reads a CSV of AGAT-format metrics, filters to the subset of metrics that
have a known Genebuild equivalent, renames them, and writes the result to a
new CSV file.

Usage:
    python select_and_rename_metrics.py --input <input.csv> --output <output.csv>
"""

import argparse
import pandas as pd
from pathlib import Path


def safe_agg(df: pd.DataFrame, columns: list[str], method: str = "sum") -> pd.Series:
    """
    Aggregate only the columns that exist in the DataFrame.

    Args:
        df: DataFrame with metrics as columns (after pivot if long format).
        columns: List of column names to aggregate.
        method: Aggregation method: 'sum', 'mean', 'min', 'max'.

    Returns:
        A Series with a single value of the aggregation.
    """
    existing = [c for c in columns if c in df.columns]
    if not existing:
        # Return 0 for sums, NaN for mean/min/max
        if method == "sum":
            return pd.Series([0], index=[0])
        else:
            return pd.Series([float("nan")], index=[0])

    if method == "sum":
        return df[existing].sum(axis=1)
    elif method == "mean":
        return df[existing].mean(axis=1)
    elif method == "min":
        return df[existing].min(axis=1)
    elif method == "max":
        return df[existing].max(axis=1)
    else:
        raise ValueError(f"Unsupported aggregation method: {method}")


def select_and_rename_metrics(input_csv: Path, output_csv: Path) -> None:
    """Filter and rename AGAT metric rows to their Genebuild equivalents.

    Reads a CSV whose rows each represent a single metric (with at least a
    ``metric`` column containing the AGAT metric name), drops any rows whose
    metric name is not present in the hard-coded mapping, renames the
    remaining metrics to their corresponding Genebuild stat names, and writes
    the resulting DataFrame to *output_csv*.

    Args:
        input_csv: Path to the input CSV file containing AGAT metric data.
        output_csv: Path where the filtered and renamed CSV will be written.

    Returns:
        None.  Results are written directly to *output_csv*.
    """

    original = pd.read_csv(input_csv)
    df = original.set_index("metric").T

    # metrics to keep and their new names
    nc_types = [
        "lnc_rna",
        "mirna",
        "ncrna",
        "rrna",
        "scrna",
        "snorna",
        "snrna",
        "trna",
        "y_rna",
    ]

    small_nc_types = [
        "snrna",
        "snorna",
        "trna",
        "mirna",
        "rrna",
        "scrna",
        "y_rna",
    ]

    df["average_coding_sequence_length"] = (
        df["mrna_total_exon_length_(bp)"] / df["mrna_number_of_mrna"]
    )

    df["nc_small_non_coding_genes"] = df[
        [f"{t}_number_of_ncrna_gene" for t in small_nc_types]
    ].sum(axis=1)

    df["nc_total_transcripts"] = df[[f"{t}_number_of_{t}" for t in nc_types]].sum(
        axis=1
    )

    df["nc_total_exons"] = safe_agg(
        df, [f"{t}_number_of_exon" for t in nc_types], method="sum"
    )

    df["nc_total_introns"] = safe_agg(
        df, [f"{t}_number_of_intron_in_exon" for t in nc_types], method="sum"
    )

    df["nc_non_coding_genes"] = df[[f"{t}_number_of_ncrna_gene" for t in nc_types]].sum(
        axis=1
    )

    df["nc_transcripts_per_gene"] = (
        df["nc_total_transcripts"] / df["nc_non_coding_genes"]
    )

    df["nc_total_exon_length"] = safe_agg(
        df, [f"{t}_total_exon_length_(bp)" for t in nc_types], method="sum"
    )

    df["nc_average_exons_per_transcript"] = (
        df["nc_total_exons"] / df["nc_total_transcripts"]
    )

    df["nc_total_exon_length"] = df[
        [f"{t}_total_exon_length_(bp)" for t in nc_types]
    ].sum(axis=1)

    df["nc_average_exon_length"] = df["nc_total_exon_length"] / df["nc_total_exons"]

    df["nc_total_intron_length"] = safe_agg(
        df, [f"{t}_total_intron_length_per_exon_(bp)" for t in nc_types], method="sum"
    )

    df["nc_average_intron_length"] = (
        df["nc_total_intron_length"] / df["nc_total_introns"]
    )

    df["nc_total_transcript_length"] = df[
        [f"{t}_total_{t}_length_(bp)" for t in nc_types]
    ].sum(axis=1)

    df["nc_average_sequence_length"] = (
        df["nc_total_transcript_length"] / df["nc_non_coding_genes"]
    )

    df["ps_average_sequence_length"] = (
        df["pseudogenic_transcript_total_exon_length_(bp)"]
        / df["pseudogenic_transcript_number_of_pseudogenic_transcript"]
    )

    df["total_transcripts"] = (
        df["mrna_number_of_mrna"]
        + df["nc_total_transcripts"]
        + df["pseudogenic_transcript_number_of_pseudogenic_transcript"]
    )

    df["total_genes"] = (
        df["mrna_number_of_gene"]
        + df["nc_non_coding_genes"]
        + df["pseudogenic_transcript_number_of_pseudogene"]
    )

    df["transcripts_per_gene"] = df["total_transcripts"] / df["total_genes"]

    new_map = {
        # coding
        "genebuild.stats.coding_genes": "mrna_number_of_gene",
        "genebuild.stats.coding_transcripts": "mrna_number_of_mrna",
        "genebuild.stats.coding_transcripts_per_gene": "mrna_mean_mrnas_per_gene",
        "genebuild.stats.total_coding_exons": "mrna_number_of_exon_in_cds",
        "genebuild.stats.total_transcript_exons": "mrna_number_of_exon",
        "genebuild.stats.total_coding_introns": "mrna_number_of_intron_in_exon",
        "genebuild.stats.longest_coding_gene_length": "mrna_longest_gene_(bp)",
        "genebuild.stats.shortest_coding_gene_length": "mrna_shortest_gene_(bp)",
        "genebuild.stats.total_transcripts": "total_transcripts",
        # coding averages
        "genebuild.stats.average_cds_length": "mrna_mean_cds_length_(bp)",
        "genebuild.stats.average_coding_exons_per_coding_tr": "mrna_mean_exons_per_cds",
        "genebuild.stats.average_coding_exons_per_transcrip": "mrna_mean_exons_per_mrna",
        "genebuild.stats.average_coding_exon_length": "mrna_mean_cds_piece_length_(bp)",
        "genebuild.stats.average_coding_exon_length": "mrna_mean_transcipt_exon_length_(bp)",
        "genebuild.stats.average_coding_genomic_span": "mrna_mean_gene_length_(bp)",
        "genebuild.stats.average_coding_intron_length": "mrna_mean_intron_in_exon_length_(bp)",
        "genebuild.stats.average_coding_sequence_length": "average_coding_sequence_length",
        "genebuild.stats.transcripts_per_gene": "transcripts_per_gene",
        # coding new
        "genebuild.stats.single_exon_coding_genes": "mrna_number_of_single_exon_gene",
        "genebuild.stats.single_exon_coding_transcripts": "mrna_number_of_single_exon_mrna",
        "genebuild.stats.overlapping_coding_genes": "mrna_number_gene_overlapping",
        "genebuild.stats.coding_transcripts_with_both_utrs": "mrna_number_of_mrnas_with_utr_both_sides",
        "genebuild.stats.coding_transcripts_with_utr": "mrna_number_of_mrnas_with_at_least_one_utr",
        "genebuild.stats.average_cds_intron_length": "mrna_mean_intron_in_cds_length_(bp)",
        "genebuild.stats.longest_transcript_length": "mrna_longest_mrna_(bp)",
        "genebuild.stats.longest_cds_length": "mrna_longest_cds_(bp)",
        "genebuild.stats.average_five_prime_utr_length": "mrna_mean_five_prime_utr_length_(bp)",
        "genebuild.stats.average_three_prime_utr_length": "mrna_mean_three_prime_utr_length_(bp)",
        # non coding
        "genebuild.stats.nc_long_non_coding_genes": "lnc_rna_number_of_ncrna_gene",
        "genebuild.stats.nc_misc_non_coding_genes": "ncrna_number_of_ncrna_gene",
        "genebuild.stats.nc_small_non_coding_genes": "nc_small_non_coding_genes",
        "genebuild.stats.nc_total_transcripts": "nc_total_transcripts",
        "genebuild.stats.nc_total_exons": "nc_total_exons",
        "genebuild.stats.nc_total_introns": "nc_total_introns",
        "genebuild.stats.nc_transcripts_per_gene": "nc_transcripts_per_gene",
        "genebuild.stats.nc_non_coding_genes": "nc_non_coding_genes",
        "genebuild.stats.nc_longest_gene_length": "lnc_rna_longest_ncrna_gene_(bp)",
        "genebuild.stats.nc_shortest_gene_length": "nc_shortest_gene_length",
        # non coding averages
        "genebuild.stats.nc_average_exons_per_transcript": "nc_average_exons_per_transcript",
        "genebuild.stats.nc_average_exon_length": "nc_average_exon_length",
        "genebuild.stats.nc_average_intron_length": "nc_average_intron_length",
        "genebuild.stats.nc_average_sequence_length": "nc_average_sequence_length",
        # non coding new
        "genebuild.stats.nc_single_exon_long_non_coding_genes": "lnc_rna_number_of_single_exon_ncrna_gene",
        "genebuild.stats.nc_overlapping_long_non_coding_genes": "lnc_rna_number_ncrna_gene_overlapping",
        "genebuild.stats.nc_mirna_genes": "mirna_number_of_ncrna_gene",
        "genebuild.stats.nc_unclassified_genes": "ncrna_number_of_ncrna_gene",
        "genebuild.stats.nc_rrna_genes": "rrna_number_of_ncrna_gene",
        "genebuild.stats.nc_scrna_genes": "scrna_number_of_ncrna_gene",
        "genebuild.stats.nc_snorna_genes": "snorna_number_of_ncrna_gene",
        "genebuild.stats.nc_snrna_genes": "snrna_number_of_ncrna_gene",
        "genebuild.stats.nc_trna_genes": "trna_number_of_ncrna_gene",
        "genebuild.stats.nc_y_rna_genes": "y_rna_number_of_ncrna_gene",
        # pseudogenes
        "genebuild.stats.ps_pseudogenes": "pseudogenic_transcript_number_of_pseudogene",
        "genebuild.stats.ps_total_transcripts": "pseudogenic_transcript_number_of_pseudogenic_transcript",
        "genebuild.stats.ps_total_exons": "pseudogenic_transcript_number_of_exon",
        "genebuild.stats.ps_total_introns": "pseudogenic_transcript_number_of_intron_in_exon",
        "genebuild.stats.ps_transcripts_per_gene": "pseudogenic_transcript_mean_pseudogenic_transcripts_per_pseudogene",
        "genebuild.stats.ps_longest_gene_length": "pseudogenic_transcript_longest_pseudogene_(bp)",
        "genebuild.stats.ps_shortest_gene_length": "pseudogenic_transcript_shortest_pseudogene_(bp)",
        # pseudogenes averages
        "genebuild.stats.ps_average_exons_per_transcript": "pseudogenic_transcript_mean_exons_per_pseudogenic_transcript",
        "genebuild.stats.ps_average_exon_length": "pseudogenic_transcript_mean_exon_length_(bp)",
        "genebuild.stats.ps_average_genomic_span": "pseudogenic_transcript_mean_pseudogene_length_(bp)",
        "genebuild.stats.ps_average_intron_length": "pseudogenic_transcript_mean_intron_in_exon_length_(bp)",
        "genebuild.stats.ps_average_sequence_length": "ps_average_sequence_length",
    }

    df_final = df.T.reset_index()
    df_final.columns = ["metric", "value"]

    agat_to_genebuild = {v: k for k, v in new_map.items()}

    df_final = df_final[df_final["metric"].isin(agat_to_genebuild.keys())].copy()
    df_final["metric"] = df_final["metric"].map(agat_to_genebuild)

    # Convert numeric columns
    def format_value(x):
        x = round(float(x), 2)
        if x.is_integer():
            return int(x)
        return x

    df_final["value"] = df_final["value"].apply(format_value)

    df_final.to_csv(output_csv, index=False)


def main() -> None:
    """Parse command-line arguments and run the metric selection/renaming step."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input CSV")
    parser.add_argument("--output", required=True, help="Output CSV")

    args = parser.parse_args()

    select_and_rename_metrics(Path(args.input), Path(args.output))


if __name__ == "__main__":
    main()

"""
Select and rename AGAT metrics to Ensembl Genebuild stat names.

Reads an AGAT stats *text* file, parses it to a table of metrics,
computes derived metrics, maps them to Genebuild stat names using a JSON
config, and writes the result to a CSV.

Usage:
    python parse_agat.py \
        --input_txt <agat_stats.txt> \
        [--mapping_json <value_map.json>] \
        --output <genebuild_stats.csv>
"""

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def numeric_series(series: pd.Series) -> pd.Series:
    """Coerce a Series to numeric values, using 0 for non-numeric entries."""
    return pd.to_numeric(series, errors="coerce").fillna(0)


def safe_agg(
    dataframe: pd.DataFrame, columns: List[str], method: str = "sum"
) -> pd.Series:
    """Aggregate only the columns that exist in the DataFrame."""
    existing = [column for column in columns if column in dataframe.columns]

    if not existing:
        if method == "sum":
            return pd.Series(0, index=dataframe.index)
        return pd.Series(float("nan"), index=dataframe.index)

    numeric_df = dataframe[existing].apply(pd.to_numeric, errors="coerce")

    if method == "sum":
        return numeric_df.sum(axis=1)
    if method == "mean":
        return numeric_df.mean(axis=1)
    if method == "min":
        return numeric_df.min(axis=1)
    if method == "max":
        return numeric_df.max(axis=1)
    raise ValueError(f"Unsupported aggregation method: {method}")


def safe_col(dataframe: pd.DataFrame, col: str) -> pd.Series:
    """Return column if present, otherwise a Series of zeros aligned to df."""
    if col in dataframe.columns:
        return numeric_series(dataframe[col])
    return pd.Series(0, index=dataframe.index)


def parse_agat_txt_to_df(stats_file: Path) -> pd.DataFrame:
    """
    Parse an AGAT statistics text report into a DataFrame with columns
    ['metric', 'value'].
    """
    rows: List[Tuple[str, str]] = []
    current_section: Optional[str] = None
    minus_isoform = False

    section_pattern = re.compile(r"-+\s*(.+?)\s*-+")
    value_pattern = re.compile(r"(.+?)\s+([0-9.]+)$")

    with stats_file.open(encoding="utf-8") as stats_handle:
        for line in stats_handle:
            line = line.strip()

            if not line:
                continue

            # Detect isoform message
            if "shortest isoforms excluded" in line.lower():
                minus_isoform = True
                continue

            # Detect section header
            section_match = section_pattern.match(line)
            if section_match:
                current_section = section_match.group(1).strip()
                minus_isoform = False
                continue

            match = value_pattern.match(line)
            if match:
                metric = match.group(1).strip()
                value = match.group(2)

                if current_section:
                    metric = f"{current_section}_{metric}"

                metric = metric.replace(" ", "_").lower()

                if minus_isoform:
                    metric = f"{metric}_minus_shortest_isoform"

                rows.append((metric, value))

    metrics_df = pd.DataFrame(rows, columns=["metric", "value"])
    return metrics_df


def select_and_rename_metrics(
    input_txt: Path, mapping_json: Optional[Path], output_csv: Path
) -> None:
    """
    Parse AGAT text, compute derived metrics, and map to Genebuild equivalents.

    If mapping_json is None, load 'value_map.json' from the same folder
    as this script.
    """
    # 0. Resolve mapping path
    if mapping_json is None:
        script_dir = Path(__file__).resolve().parent
        mapping_json = script_dir.parent / "config" / "value_map.json"

    # 1. Parse txt -> DataFrame of metrics
    original = parse_agat_txt_to_df(input_txt)

    # 2. Pivot to one-row-wide DataFrame with metrics as columns
    metrics_df = original.set_index("metric").T

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

    metrics_df["average_coding_sequence_length"] = safe_col(
        metrics_df, "mrna_total_exon_length_(bp)"
    ) / safe_col(metrics_df, "mrna_number_of_mrna")

    metrics_df["nc_small_non_coding_genes"] = safe_agg(
        metrics_df, [f"{nc_type}_number_of_ncrna_gene" for nc_type in small_nc_types]
    )

    metrics_df["nc_total_transcripts"] = safe_agg(
        metrics_df, [f"{nc_type}_number_of_{nc_type}" for nc_type in nc_types]
    )

    metrics_df["nc_total_exons"] = safe_agg(
        metrics_df, [f"{nc_type}_number_of_exon" for nc_type in nc_types], method="sum"
    )

    metrics_df["nc_total_introns"] = safe_agg(
        metrics_df,
        [f"{nc_type}_number_of_intron_in_exon" for nc_type in nc_types],
        method="sum",
    )

    metrics_df["nc_non_coding_genes"] = safe_agg(
        metrics_df, [f"{nc_type}_number_of_ncrna_gene" for nc_type in nc_types]
    )

    metrics_df["nc_transcripts_per_gene"] = (
        metrics_df["nc_total_transcripts"] / metrics_df["nc_non_coding_genes"]
    )

    metrics_df["nc_total_exon_length"] = safe_agg(
        metrics_df, [f"{nc_type}_total_exon_length_(bp)" for nc_type in nc_types]
    )

    metrics_df["nc_average_exons_per_transcript"] = (
        metrics_df["nc_total_exons"] / metrics_df["nc_total_transcripts"]
    )

    metrics_df["nc_average_exon_length"] = (
        metrics_df["nc_total_exon_length"] / metrics_df["nc_total_exons"]
    )

    metrics_df["nc_total_intron_length"] = safe_agg(
        metrics_df,
        [f"{nc_type}_total_intron_length_per_exon_(bp)" for nc_type in nc_types],
        method="sum",
    )

    metrics_df["nc_average_intron_length"] = (
        metrics_df["nc_total_intron_length"] / metrics_df["nc_total_introns"]
    )

    metrics_df["nc_total_transcript_length"] = safe_agg(
        metrics_df, [f"{nc_type}_total_{nc_type}_length_(bp)" for nc_type in nc_types]
    )

    metrics_df["nc_average_sequence_length"] = (
        metrics_df["nc_total_transcript_length"] / metrics_df["nc_non_coding_genes"]
    )

    # Only create ps_average_sequence_length if both columns exist
    if (
        "pseudogenic_transcript_total_exon_length_(bp)" in metrics_df.columns
        and "pseudogenic_transcript_number_of_pseudogenic_transcript"
        in metrics_df.columns
    ):
        metrics_df["ps_average_sequence_length"] = safe_col(
            metrics_df, "pseudogenic_transcript_total_exon_length_(bp)"
        ) / safe_col(
            metrics_df, "pseudogenic_transcript_number_of_pseudogenic_transcript"
        )

    metrics_df["total_transcripts"] = (
        safe_col(metrics_df, "mrna_number_of_mrna")
        + safe_col(metrics_df, "nc_total_transcripts")
        + safe_col(
            metrics_df, "pseudogenic_transcript_number_of_pseudogenic_transcript"
        )
    )

    metrics_df["total_genes"] = (
        safe_col(metrics_df, "mrna_number_of_gene")
        + safe_col(metrics_df, "nc_non_coding_genes")
        + safe_col(metrics_df, "pseudogenic_transcript_number_of_pseudogene")
    )

    metrics_df["transcripts_per_gene"] = (
        metrics_df["total_transcripts"] / metrics_df["total_genes"]
    )

    # 3. Load mapping from JSON (genebuild -> AGAT metric)
    with Path(mapping_json).open(encoding="utf-8") as mapping_handle:
        new_map: Dict[str, str] = json.load(mapping_handle)

    # 4. Apply mapping
    df_final = metrics_df.T.reset_index()
    df_final.columns = ["metric", "value"]

    # We need AGAT -> Genebuild for filtering & renaming
    agat_to_genebuild = {v: k for k, v in new_map.items()}

    df_final = df_final[df_final["metric"].isin(agat_to_genebuild.keys())].copy()
    df_final["metric"] = df_final["metric"].map(agat_to_genebuild)

    # 5. Convert numeric columns
    def format_value(value):
        rounded_value = round(float(value), 2)
        if rounded_value.is_integer():
            return int(rounded_value)
        return rounded_value

    df_final["value"] = df_final["value"].apply(format_value)

    df_final.to_csv(output_csv, index=False)


def main() -> None:
    """Parse command-line arguments and run the parsing/mapping step."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_txt",
        required=True,
        help="AGAT stats text file (output of agat_sp_statistics.pl)",
    )
    parser.add_argument(
        "--mapping_json",
        help="JSON mapping {genebuild.stat: agat_metric_name} "
        "(defaults to value_map.json in this folder)",
    )
    parser.add_argument("--output", required=True, help="Output CSV")

    args = parser.parse_args()

    mapping_path = Path(args.mapping_json) if args.mapping_json else None

    select_and_rename_metrics(
        Path(args.input_txt),
        mapping_path,
        Path(args.output),
    )


if __name__ == "__main__":
    main()

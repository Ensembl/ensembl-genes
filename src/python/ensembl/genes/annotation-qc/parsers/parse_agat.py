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
import csv
import json
import re
from pathlib import Path
import pandas as pd


def safe_agg(df: pd.DataFrame, columns: list[str], method: str = "sum") -> pd.Series:
    """Aggregate only the columns that exist in the DataFrame."""
    existing = [c for c in columns if c in df.columns]

    if not existing:
        if method == "sum":
            return pd.Series(0, index=df.index)
        else:
            return pd.Series(float("nan"), index=df.index)

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


def safe_col(df: pd.DataFrame, col: str) -> pd.Series:
    """Return column if present, otherwise a Series of zeros aligned to df."""
    if col in df.columns:
        return df[col]
    return pd.Series(0, index=df.index)


def parse_agat_txt_to_df(stats_file: Path) -> pd.DataFrame:
    """
    Parse an AGAT statistics text report into a DataFrame with columns
    ['metric', 'value'].
    """
    rows: list[tuple[str, str]] = []
    current_section: str | None = None
    minus_isoform = False

    section_pattern = re.compile(r"-+\s*(.+?)\s*-+")
    value_pattern = re.compile(r"(.+?)\s+([0-9.]+)$")

    with stats_file.open() as f:
        for line in f:
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

    df = pd.DataFrame(rows, columns=["metric", "value"])
    return df


def select_and_rename_metrics(
    input_txt: Path, mapping_json: Path | None, output_csv: Path
) -> None:
    """
    Parse AGAT text, compute derived metrics, and map to Genebuild equivalents.

    If mapping_json is None, load 'value_map.json' from the same folder
    as this script.
    """
    # 0. Resolve mapping path
    if mapping_json is None:
        script_dir = Path(__file__).resolve().parent
        mapping_json = script_dir / "value_map.json"

    # 1. Parse txt -> DataFrame of metrics
    original = parse_agat_txt_to_df(input_txt)

    # 2. Pivot to one-row-wide DataFrame with metrics as columns
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

    df["average_coding_sequence_length"] = safe_col(
        df, "mrna_total_exon_length_(bp)"
    ) / safe_col(df, "mrna_number_of_mrna")

    df["nc_small_non_coding_genes"] = safe_agg(
        df, [f"{t}_number_of_ncrna_gene" for t in small_nc_types]
    )

    df["nc_total_transcripts"] = safe_agg(df, [f"{t}_number_of_{t}" for t in nc_types])

    df["nc_total_exons"] = safe_agg(
        df, [f"{t}_number_of_exon" for t in nc_types], method="sum"
    )

    df["nc_total_introns"] = safe_agg(
        df, [f"{t}_number_of_intron_in_exon" for t in nc_types], method="sum"
    )

    df["nc_non_coding_genes"] = safe_agg(
        df, [f"{t}_number_of_ncrna_gene" for t in nc_types]
    )

    df["nc_transcripts_per_gene"] = (
        df["nc_total_transcripts"] / df["nc_non_coding_genes"]
    )

    df["nc_total_exon_length"] = safe_agg(
        df, [f"{t}_total_exon_length_(bp)" for t in nc_types]
    )

    df["nc_average_exons_per_transcript"] = (
        df["nc_total_exons"] / df["nc_total_transcripts"]
    )

    df["nc_average_exon_length"] = df["nc_total_exon_length"] / df["nc_total_exons"]

    df["nc_total_intron_length"] = safe_agg(
        df, [f"{t}_total_intron_length_per_exon_(bp)" for t in nc_types], method="sum"
    )

    df["nc_average_intron_length"] = (
        df["nc_total_intron_length"] / df["nc_total_introns"]
    )

    df["nc_total_transcript_length"] = safe_agg(
        df, [f"{t}_total_{t}_length_(bp)" for t in nc_types]
    )

    df["nc_average_sequence_length"] = (
        df["nc_total_transcript_length"] / df["nc_non_coding_genes"]
    )

    # Only create ps_average_sequence_length if both columns exist
    if (
        "pseudogenic_transcript_total_exon_length_(bp)" in df.columns
        and "pseudogenic_transcript_number_of_pseudogenic_transcript" in df.columns
    ):
        df["ps_average_sequence_length"] = (
            df["pseudogenic_transcript_total_exon_length_(bp)"]
            / df["pseudogenic_transcript_number_of_pseudogenic_transcript"]
        )

    df["total_transcripts"] = (
        safe_col(df, "mrna_number_of_mrna")
        + safe_col(df, "nc_total_transcripts")
        + safe_col(df, "pseudogenic_transcript_number_of_pseudogenic_transcript")
    )

    df["total_genes"] = (
        safe_col(df, "mrna_number_of_gene")
        + safe_col(df, "nc_non_coding_genes")
        + safe_col(df, "pseudogenic_transcript_number_of_pseudogene")
    )

    df["transcripts_per_gene"] = df["total_transcripts"] / df["total_genes"]

    # 3. Load mapping from JSON (genebuild -> AGAT metric)
    with Path(mapping_json).open() as fh:
        new_map: dict[str, str] = json.load(fh)

    # 4. Apply mapping
    df_final = df.T.reset_index()
    df_final.columns = ["metric", "value"]

    # We need AGAT -> Genebuild for filtering & renaming
    agat_to_genebuild = {v: k for k, v in new_map.items()}

    df_final = df_final[df_final["metric"].isin(agat_to_genebuild.keys())].copy()
    df_final["metric"] = df_final["metric"].map(agat_to_genebuild)

    # 5. Convert numeric columns
    def format_value(x):
        x = round(float(x), 2)
        if x.is_integer():
            return int(x)
        return x

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

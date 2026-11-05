"""
Parse annotation files using PyRanges1.

The main entry point is parse_annotation, which detects whether the input file
is GFF3 or GTF and calls the appropriate parser.
"""

from pathlib import Path

import pandas as pd
import pyranges1 as pr

GENE_FEATURE_TYPES = frozenset(
    {"gene", "ncRNA_gene", "pseudogene", "transposable_element_gene"}
)
TRANSCRIPT_FEATURE_TYPES = frozenset(
    {
        "mRNA",
        "transcript",
        "snoRNA",
        "snRNA",
        "miRNA",
        "rRNA",
        "tRNA",
        "lnc_RNA",
        "lncRNA",
        "ncRNA",
        "antisense_RNA",
        "sRNA",
        "scaRNA",
        "pseudogenic_transcript",
        "scRNA",
        "Y_RNA",
    }
)


def _as_dataframe(annotation) -> pd.DataFrame:
    """Return a copy of a PyRanges or DataFrame annotation table."""
    if isinstance(annotation, pd.DataFrame):
        return annotation.copy()
    if hasattr(annotation, "df"):
        return annotation.df.copy()
    if hasattr(annotation, "as_df"):
        return annotation.as_df().copy()
    raise TypeError("annotation must be a PyRanges or pandas DataFrame")


def _clean_id(value: object) -> str:
    """Remove common GFF/GTF identifier prefixes from one value."""
    value = "" if pd.isna(value) else str(value)
    value = value.split(",", 1)[0].strip()
    for prefix in ("gene:", "transcript:", "mRNA:"):
        if value.startswith(prefix):
            return value[len(prefix) :]
    return value


def standardize_annotation(annotation) -> pd.DataFrame:
    """Normalize common GFF3/GTF columns used by annotation-QC metrics."""
    df = _as_dataframe(annotation)
    required = {"Feature", "Chromosome", "Start", "End", "Strand"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Annotation is missing required columns: {sorted(missing)}")

    for column in (
        "ID",
        "Parent",
        "gene_id",
        "transcript_id",
        "biotype",
        "tag",
        "phase",
    ):
        if column not in df.columns:
            df[column] = pd.NA

    feature = df["Feature"].astype(str)
    gene_mask = feature.isin(GENE_FEATURE_TYPES)
    transcript_mask = feature.isin(TRANSCRIPT_FEATURE_TYPES)
    child_mask = ~(gene_mask | transcript_mask)

    for index in df.index:
        if gene_mask.loc[index] and pd.isna(df.at[index, "gene_id"]):
            df.at[index, "gene_id"] = _clean_id(df.at[index, "ID"])
        if transcript_mask.loc[index] and pd.isna(df.at[index, "transcript_id"]):
            df.at[index, "transcript_id"] = _clean_id(df.at[index, "ID"])
        if transcript_mask.loc[index] and pd.isna(df.at[index, "gene_id"]):
            df.at[index, "gene_id"] = _clean_id(df.at[index, "Parent"])
        if child_mask.loc[index] and pd.isna(df.at[index, "transcript_id"]):
            df.at[index, "transcript_id"] = _clean_id(df.at[index, "Parent"])

    transcript_to_gene = (
        df.loc[transcript_mask, ["transcript_id", "gene_id"]]
        .dropna(subset=["transcript_id", "gene_id"])
        .drop_duplicates("transcript_id")
        .set_index("transcript_id")["gene_id"]
        .to_dict()
    )
    missing_gene = child_mask & df["gene_id"].isna()
    df.loc[missing_gene, "gene_id"] = df.loc[missing_gene, "transcript_id"].map(
        transcript_to_gene
    )

    if "gene_type" in df.columns:
        df["biotype"] = df["biotype"].fillna(df["gene_type"])
    if "transcript_type" in df.columns:
        df["biotype"] = df["biotype"].fillna(df["transcript_type"])
    if "Frame" in df.columns:
        df["phase"] = df["phase"].fillna(df["Frame"])

    df["feature_length"] = df["End"].astype(int) - df["Start"].astype(int)
    return df


def parse_gff3(file_path: str):
    """
    Parse a GFF3 file using PyRanges1.
    Args:
            file_path: Path to the GFF3 file
    Returns:
            PyRanges object
    """
    gff3_file = pr.read_gff3(file_path)
    return gff3_file


def parse_gtf(file_path: str):
    """
    Parse a GTF file using PyRanges1.
    Args:
            file_path: Path to the GTF file
    Returns:
            PyRanges object
    """

    gtf_file = pr.read_gtf(file_path)
    return gtf_file


def parse_annotation(file_path: str):
    """
    Parse a GTF or GFF3 file using PyRanges1. Main function and entry to the script
    Args:
            file_path: Path to the GTF file
    Returns:
            PyRanges object
    """

    path = Path(file_path)
    suffixes = [s.lower() for s in path.suffixes]
    real_suffix = (
        suffixes[-2]
        if suffixes and suffixes[-1] == ".gz" and len(suffixes) >= 2
        else suffixes[-1]
    )

    if real_suffix in {".gff3", ".gff"}:
        print(f"Parsing GFF3/GFF... ({file_path})")
        data = parse_gff3(file_path)
    elif real_suffix == ".gtf":
        print(f"Parsing GTF... ({file_path})")
        data = parse_gtf(file_path)
    else:
        raise ValueError(
            "Unsupported file type. Use .gff, .gff3, .gtf, or their .gz variants"
        )

    return standardize_annotation(data)

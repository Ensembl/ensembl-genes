import pyranges1 as pr
from pathlib import Path

import pandas as pd


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

FIVE_PRIME_UTR_FEATURE_TYPES = frozenset({"five_prime_UTR"})
THREE_PRIME_UTR_FEATURE_TYPES = frozenset({"three_prime_UTR"})
CDS_FEATURE_TYPES = frozenset({"CDS"})
EXON_FEATURE_TYPES = frozenset({"exon"})


def _as_dataframe(annotation) -> pd.DataFrame:
    """Return annotation data as a pandas DataFrame."""
    if isinstance(annotation, pd.DataFrame):
        return annotation.copy()
    if hasattr(annotation, "df"):
        return annotation.df.copy()
    if hasattr(annotation, "as_df"):
        return annotation.as_df().copy()
    return annotation.copy()


def _strip_prefix(values: pd.Series, prefix: str) -> pd.Series:
    return values.astype("string").str.removeprefix(prefix)


def _first_parent(values: pd.Series) -> pd.Series:
    return values.astype("string").str.split(",").str[0]


def standardize_annotation(annotation) -> pd.DataFrame:
    """
    Normalize parser output to the annotation table shape used by QC metrics.

    The returned DataFrame keeps the original parser columns and fills common
    structural columns where possible: ``gene_id``, ``transcript_id``,
    ``biotype``, and ``feature_length``. Missing optional metadata stays empty
    rather than being guessed.
    """
    df = _as_dataframe(annotation)

    for col in (
        "Feature",
        "ID",
        "Parent",
        "gene_id",
        "transcript_id",
        "biotype",
        "tag",
        "phase",
    ):
        if col not in df.columns:
            df[col] = pd.NA

    gene_mask = df["Feature"].isin(GENE_FEATURE_TYPES)
    tx_mask = df["Feature"].isin(TRANSCRIPT_FEATURE_TYPES)
    child_mask = ~(gene_mask | tx_mask)

    missing_gene_id = gene_mask & df["gene_id"].isna() & df["ID"].notna()
    df.loc[missing_gene_id, "gene_id"] = _strip_prefix(
        df.loc[missing_gene_id, "ID"], "gene:"
    )

    missing_tx_id = tx_mask & df["transcript_id"].isna() & df["ID"].notna()
    df.loc[missing_tx_id, "transcript_id"] = _strip_prefix(
        df.loc[missing_tx_id, "ID"], "transcript:"
    )

    missing_tx_gene = tx_mask & df["gene_id"].isna() & df["Parent"].notna()
    df.loc[missing_tx_gene, "gene_id"] = _strip_prefix(
        _first_parent(df.loc[missing_tx_gene, "Parent"]), "gene:"
    )

    missing_child_tx = child_mask & df["transcript_id"].isna() & df["Parent"].notna()
    df.loc[missing_child_tx, "transcript_id"] = _strip_prefix(
        _first_parent(df.loc[missing_child_tx, "Parent"]), "transcript:"
    )

    tx_to_gene = (
        df.loc[tx_mask & df["transcript_id"].notna(), ["transcript_id", "gene_id"]]
        .dropna(subset=["gene_id"])
        .drop_duplicates(subset=["transcript_id"])
        .set_index("transcript_id")["gene_id"]
        .to_dict()
    )
    missing_child_gene = child_mask & df["gene_id"].isna()
    df.loc[missing_child_gene, "gene_id"] = df.loc[
        missing_child_gene, "transcript_id"
    ].map(tx_to_gene)

    if "gene_type" in df.columns:
        missing_biotype = df["biotype"].isna()
        df.loc[missing_biotype, "biotype"] = df.loc[missing_biotype, "gene_type"]

    if "Frame" in df.columns:
        missing_phase = df["phase"].isna()
        df.loc[missing_phase, "phase"] = df.loc[missing_phase, "Frame"]

    if "Start" in df.columns and "End" in df.columns:
        df["feature_length"] = df["End"] - df["Start"]

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
    Parse a GTF/GFF3 file and return a standardized annotation DataFrame.

    Args:
            file_path: Path to the GTF or GFF3 file
    Returns:
            pandas DataFrame with normalized structural columns
    """

    path = Path(file_path)
    suffixes = [s.lower() for s in path.suffixes]
    real_suffix = (
        suffixes[-2]
        if suffixes and suffixes[-1] == ".gz" and len(suffixes) >= 2
        else suffixes[-1]
    )

    if real_suffix == ".gff3":
        print(f"Parsing GFF3... ({file_path})")
        data = parse_gff3(file_path)
    elif real_suffix == ".gtf":
        print(f"Parsing GTF... ({file_path})")
        data = parse_gtf(file_path)
    else:
        raise ValueError("Unsupported file type. Use .gff3, .gtf, .gff3.gz, or .gtf.gz")

    return standardize_annotation(data)

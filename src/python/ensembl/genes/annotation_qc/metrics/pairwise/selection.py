"""
Reference filtering and transcript-selection policies for pairwise comparison.

Inputs and outputs are comparison-schema DataFrames (see
parsers/annotation_normalise.py). Nothing here reads or writes files.

Policies (unchanged from gmb-compare):

    evaluation mode
        all             no filtering
        protein_coding  genes whose gene_biotype is protein_coding, with all of
                        their transcripts (transcript biotype is not checked)
        cds_only        transcripts that have CDS, and genes with such a transcript
        canonical       one transcript per gene, as transcript selection 'canonical'

    transcript selection (per gene, among its transcripts in file order)
        all             keep every transcript
        canonical       first transcript tagged Ensembl_canonical, otherwise as
                        longest_cds
        longest_cds     longest total CDS; if no transcript has CDS, longest total
                        exon length; ties go to the first transcript
"""

import pandas as pd

from ensembl.genes.annotation_qc.parsers.regions import Region

EVALUATION_MODES = ("all", "protein_coding", "cds_only", "canonical")
TRANSCRIPT_SELECTIONS = ("all", "longest_cds", "canonical")


def _keep(
    df: pd.DataFrame, gene_ids: set, transcript_ids: set | None = None
) -> pd.DataFrame:
    """Keep genes in gene_ids; keep other rows whose transcript (or gene) is kept."""
    is_gene = df["Feature"] == "gene"
    if transcript_ids is None:
        keep = df["gene_id"].isin(gene_ids)
    else:
        keep = (is_gene & df["gene_id"].isin(gene_ids)) | (
            ~is_gene & df["transcript_id"].isin(transcript_ids)
        )
    return df[keep].copy()


def filter_by_biotype(
    df: pd.DataFrame,
    gene_biotypes: list[str] | None = None,
    transcript_biotypes: list[str] | None = None,
) -> pd.DataFrame:
    """
    Keep genes and/or transcripts with the given biotypes.

    A gene filter keeps whole genes. A transcript filter keeps matching
    transcripts and the genes that still have one.
    Args:
            df: Comparison-schema DataFrame
            gene_biotypes: gene_biotype values to keep
            transcript_biotypes: transcript_biotype values to keep
    Returns:
            Filtered DataFrame
    """
    if gene_biotypes:
        genes = df[(df["Feature"] == "gene") & df["gene_biotype"].isin(gene_biotypes)]
        df = _keep(df, set(genes["gene_id"]))
    if transcript_biotypes:
        transcripts = df[
            (df["Feature"] == "transcript")
            & df["transcript_biotype"].isin(transcript_biotypes)
        ]
        df = _keep(df, set(transcripts["gene_id"]), set(transcripts["transcript_id"]))
    return df


def select_transcripts(df: pd.DataFrame, mode: str = "all") -> pd.DataFrame:
    """
    Keep one representative transcript per gene (see module docstring).

    Genes are kept even when they have no transcripts.
    Args:
            df: Comparison-schema DataFrame
            mode: "all", "longest_cds" or "canonical"
    Returns:
            Filtered DataFrame
    """
    if mode not in TRANSCRIPT_SELECTIONS:
        raise ValueError(f"Unknown transcript selection '{mode}'")
    transcripts = df[df["Feature"] == "transcript"]
    if mode == "all" or transcripts.empty:
        return df

    lengths = (df["End"] - df["Start"]).rename("length")
    cds_length = lengths[df["Feature"] == "CDS"].groupby(df["transcript_id"]).sum()
    exon_length = lengths[df["Feature"] == "exon"].groupby(df["transcript_id"]).sum()
    ranked = pd.DataFrame(
        {
            "gene_id": transcripts["gene_id"].to_numpy(),
            "transcript_id": transcripts["transcript_id"].to_numpy(),
            "canonical": transcripts["tags"]
            .str.contains("Ensembl_canonical", regex=False)
            .to_numpy(),
        }
    ).drop_duplicates("transcript_id")
    ranked["cds"] = ranked["transcript_id"].map(cds_length).fillna(0)
    ranked["exons"] = ranked["transcript_id"].map(exon_length).fillna(0)
    ranked["order"] = range(len(ranked))

    def choose(group: pd.DataFrame) -> str:
        if mode == "canonical" and group["canonical"].any():
            return group.loc[group["canonical"], "transcript_id"].iloc[0]
        column = "cds" if group["cds"].max() > 0 else "exons"
        # idxmax returns the first maximum, i.e. the earliest transcript on ties.
        return group.loc[group[column].idxmax(), "transcript_id"]

    chosen = {choose(group) for _, group in ranked.groupby("gene_id", sort=False)}
    is_gene = df["Feature"] == "gene"
    return df[is_gene | df["transcript_id"].isin(chosen)].copy()


def apply_evaluation_mode(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    """
    Apply an evaluation-mode preset to a reference annotation.
    Args:
            df: Comparison-schema DataFrame
            mode: One of EVALUATION_MODES
    Returns:
            Filtered DataFrame
    """
    if mode == "all":
        return df
    if mode == "protein_coding":
        return filter_by_biotype(df, gene_biotypes=["protein_coding"])
    if mode == "cds_only":
        with_cds = set(df.loc[df["Feature"] == "CDS", "transcript_id"])
        transcripts = df[
            (df["Feature"] == "transcript") & df["transcript_id"].isin(with_cds)
        ]
        return _keep(df, set(transcripts["gene_id"]), set(transcripts["transcript_id"]))
    if mode == "canonical":
        return select_transcripts(df, "canonical")
    raise ValueError(f"Unknown evaluation mode '{mode}'")


def drop_seqnames(df: pd.DataFrame, seqnames: set[str]) -> pd.DataFrame:
    """
    Remove every feature on the given sequences.
    Args:
            df: Comparison-schema DataFrame
            seqnames: Sequence names to drop
    Returns:
            Filtered DataFrame
    """
    return df[~df["Chromosome"].isin(seqnames)].copy()


def subset_to_regions(df: pd.DataFrame, regions: list[Region]) -> pd.DataFrame:
    """
    Keep whole genes whose span overlaps any region, with all of their features.
    Args:
            df: Comparison-schema DataFrame
            regions: Regions with zero-based half-open coordinates
    Returns:
            Filtered DataFrame
    """
    genes = df[df["Feature"] == "gene"]
    hit = pd.Series(False, index=genes.index)
    for region in regions:
        on_seq = genes["Chromosome"] == region.seqname
        if region.start is None:
            hit |= on_seq
        else:
            hit |= (
                on_seq & (genes["Start"] < region.end) & (genes["End"] > region.start)
            )
    return _keep(df, set(genes.loc[hit, "gene_id"]))


def summarise_annotation(df: pd.DataFrame) -> dict:
    """
    Count genes, transcripts, exons, CDS, biotypes and canonical tags.
    Args:
            df: Comparison-schema DataFrame
    Returns:
            dict of counts used by the reference filter audit
    """
    genes = df[df["Feature"] == "gene"]
    transcripts = df[df["Feature"] == "transcript"]
    return {
        "genes": len(genes),
        "transcripts": int(transcripts["transcript_id"].nunique()),
        "exons": int((df["Feature"] == "exon").sum()),
        "cds": int((df["Feature"] == "CDS").sum()),
        "cds_containing_transcripts": int(
            df.loc[df["Feature"] == "CDS", "transcript_id"].nunique()
        ),
        "canonical_transcripts": int(
            transcripts["tags"].str.contains("Ensembl_canonical", regex=False).sum()
        ),
        "gene_biotype_counts": {
            k: int(v) for k, v in genes["gene_biotype"].value_counts().items()
        },
        "transcript_biotype_counts": {
            k: int(v)
            for k, v in transcripts["transcript_biotype"].value_counts().items()
        },
    }


from __future__ import annotations

import pandas as pd
from typing import Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_GENE_FEATURE_TYPES = frozenset(
    {"gene", "ncRNA_gene", "pseudogene", "transposable_element_gene"}
)

_TRANSCRIPT_FEATURE_TYPES = frozenset(
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

# Mirrors leannes_stats.py biotype sets
_SMALL_NC_BIOTYPES = frozenset(
    {"miRNA", "rRNA", "snRNA", "snoRNA", "tRNA", "scaRNA", "sRNA", "ncRNA"}
)
_LONG_NC_BIOTYPES = frozenset({"lncRNA", "lincRNA", "lnc_RNA"})


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _biotype_group(biotype: str) -> str:
    """
    Classify a gene biotype string following leannes_stats.py conventions.

    Returns one of: ``coding``, ``pseudogene``, ``lnoncoding``,
    ``snoncoding``, ``mnoncoding``.
    """
    bt = (biotype or "").strip()
    if "pseudogene" in bt.lower():
        return "pseudogene"
    if bt == "protein_coding":
        return "coding"
    if bt in _LONG_NC_BIOTYPES:
        return "lnoncoding"
    if bt in _SMALL_NC_BIOTYPES:
        return "snoncoding"
    return "mnoncoding"


def _compute_introns(
    exon_df: pd.DataFrame, tx_col: str = "transcript_id"
) -> pd.DataFrame:
    """
    Derive intron lengths from exon coordinates using a vectorised approach.

    Exons are sorted by start position within each transcript; a running
    cumulative-max of exon end positions is used so that overlapping or
    adjacent exons produce no intron.  Coordinates are assumed to follow
    the pyranges convention (0-based, half-open), so intron length equals
    ``next_exon_start - prev_cummax_end``.

    Parameters
    ----------
    exon_df:
        DataFrame with at least *tx_col*, ``Start``, and ``End`` columns.
    tx_col:
        Column name used to group exons by transcript.

    Returns
    -------
    pd.DataFrame
        One row per intron with columns *tx_col* and ``intron_len``.
    """
    if exon_df.empty:
        return pd.DataFrame(columns=[tx_col, "intron_len"])

    df = exon_df[[tx_col, "Start", "End"]].copy().sort_values([tx_col, "Start"])
    df["cum_max_end"] = df.groupby(tx_col)["End"].cummax()
    df["prev_cum_max_end"] = df.groupby(tx_col)["cum_max_end"].shift(1)

    gap_mask = df["prev_cum_max_end"].notna() & (df["prev_cum_max_end"] < df["Start"])
    introns = df[gap_mask].copy()
    introns["intron_len"] = introns["Start"] - introns["prev_cum_max_end"]
    return introns[[tx_col, "intron_len"]]


def _canonical_info(
    tx_df: pd.DataFrame,
    exon_df: pd.DataFrame,
    cds_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Return the spliced exon and CDS lengths of the canonical transcript per gene.

    Selection mirrors leannes_stats.py / Ensembl conventions:

    1. Prefer transcripts tagged with ``Ensembl_canonical``.
    2. Among eligible transcripts rank by CDS spliced length DESC then
       exon spliced length DESC.
    3. If no transcript is tagged, all transcripts for that gene are eligible.

    Parameters
    ----------
    tx_df:
        Transcript rows with columns ``transcript_id``, ``gene_id``, ``tag``.
    exon_df:
        Exon rows with columns ``transcript_id``, ``exon_len``.
    cds_df:
        CDS rows with columns ``transcript_id``, ``cds_len``.

    Returns
    -------
    pd.DataFrame
        Indexed by ``gene_id`` with columns ``exon_len`` and ``cds_len``
        for the canonical transcript of each gene.
    """
    exon_len_per_tx = (
        exon_df.groupby("transcript_id")["exon_len"].sum()
        if not exon_df.empty
        else pd.Series(dtype=float)
    )
    cds_len_per_tx = (
        cds_df.groupby("transcript_id")["cds_len"].sum()
        if not cds_df.empty
        else pd.Series(dtype=float)
    )

    tx = tx_df[["transcript_id", "gene_id", "tag"]].copy()
    tx["exon_len"] = tx["transcript_id"].map(exon_len_per_tx).fillna(0.0)
    tx["cds_len"] = tx["transcript_id"].map(cds_len_per_tx).fillna(0.0)
    tx["is_canonical"] = tx["tag"].str.contains("Ensembl_canonical", na=False)

    gene_has_canonical = tx.groupby("gene_id")["is_canonical"].any()
    tx["gene_has_canonical"] = tx["gene_id"].map(gene_has_canonical)

    eligible = tx[tx["is_canonical"] | ~tx["gene_has_canonical"]]
    best = (
        eligible.sort_values(["cds_len", "exon_len"], ascending=False)
        .groupby("gene_id")[["exon_len", "cds_len"]]
        .first()
    )
    return best


def _prepare_dataframes(gff) -> tuple:
    """
    Split a pyranges GFF3 object into typed DataFrames with resolved IDs.

    Parses the ``Parent`` attribute to derive ``gene_id`` for transcript
    rows and ``transcript_id`` for exon/CDS rows, then propagates
    ``gene_id`` down to exon and CDS rows via a transcript→gene mapping.

    Genes with CDS features (regardless of biotype) are classified as
    ``coding`` unless already classified as ``pseudogene``.

    Returns
    -------
    tuple of (gene_df, tx_df, exon_df, cds_df)
        pandas DataFrames with resolved ``gene_id`` / ``transcript_id``
        columns and pre-computed length columns (``gene_len``,
        ``exon_len``, ``cds_len``).  ``gene_df`` has a ``group`` column.
    """
    # --- Gene features ---
    gene_df = gff[gff["Feature"].isin(_GENE_FEATURE_TYPES)].reset_index(drop=True).copy()
    gene_df["gene_len"] = gene_df["End"] - gene_df["Start"]
    # pyranges1 already parses gene_id; fall back to stripping the ID prefix
    missing = gene_df["gene_id"].isna()
    if missing.any():
        gene_df.loc[missing, "gene_id"] = (
            gene_df.loc[missing, "ID"].str.split(":", n=1).str[-1]
        )
    gene_df["group"] = gene_df["biotype"].fillna("").apply(_biotype_group)

    # --- Transcript features ---
    tx_df = (
        gff[gff["Feature"].isin(_TRANSCRIPT_FEATURE_TYPES)]
        .reset_index(drop=True)
        .copy()
    )
    # Derive gene_id from Parent (e.g. "gene:ENSGID" → "ENSGID")
    tx_df["gene_id"] = (
        tx_df["Parent"].str.split(",").str[0].str.removeprefix("gene:")
    )

    # --- Exon features ---
    exon_df = gff[gff["Feature"] == "exon"].reset_index(drop=True).copy()
    exon_df["transcript_id"] = (
        exon_df["Parent"].str.split(",").str[0].str.removeprefix("transcript:")
    )
    exon_df["exon_len"] = exon_df["End"] - exon_df["Start"]

    # --- CDS features ---
    cds_df = gff[gff["Feature"] == "CDS"].reset_index(drop=True).copy()
    cds_df["transcript_id"] = (
        cds_df["Parent"].str.split(",").str[0].str.removeprefix("transcript:")
    )
    cds_df["cds_len"] = cds_df["End"] - cds_df["Start"]

    # Propagate gene_id to exon/CDS via transcript→gene lookup
    tx_to_gene = tx_df.set_index("transcript_id")["gene_id"].to_dict()
    exon_df["gene_id"] = exon_df["transcript_id"].map(tx_to_gene)
    cds_df["gene_id"] = cds_df["transcript_id"].map(tx_to_gene)

    # Genes that have CDS (and aren't pseudogenes) are coding
    genes_with_cds = set(cds_df["gene_id"].dropna())
    coding_override = (
        gene_df["gene_id"].isin(genes_with_cds) & (gene_df["group"] != "pseudogene")
    )
    gene_df.loc[coding_override, "group"] = "coding"

    return gene_df, tx_df, exon_df, cds_df


def _stats_for_group(
    gene_df: pd.DataFrame,
    tx_df: pd.DataFrame,
    exon_df: pd.DataFrame,
    cds_df: pd.DataFrame,
    group: str,
    include_cds: bool = False,
) -> tuple:
    """
    Compute gene/transcript/exon/intron summary statistics for a gene group.

    Used as a shared backend by the three public ``compute_*`` functions.

    Returns
    -------
    tuple of (metrics_dict, canonical_info_df)
        *metrics_dict* contains counts and averages.
        *canonical_info_df* is the per-gene canonical transcript table
        (columns ``exon_len``, ``cds_len``) needed by callers that also
        want ``total_coding_sequence_length``.
    """
    g = gene_df[gene_df["group"] == group]
    gene_count = len(g)

    gene_ids = set(g["gene_id"])
    txs = tx_df[tx_df["gene_id"].isin(gene_ids)]
    exons = exon_df[exon_df["gene_id"].isin(gene_ids)]
    cds = cds_df[cds_df["gene_id"].isin(gene_ids)] if include_cds else cds_df.iloc[:0]

    canon = _canonical_info(txs, exons, cds)

    total_tx = len(txs)
    total_exons = len(exons)
    total_exon_len = int(exons["exon_len"].sum())

    introns = _compute_introns(exons)
    total_introns = len(introns)
    total_intron_len = int(introns["intron_len"].sum()) if not introns.empty else 0

    metrics: dict = {
        "gene_count": gene_count,
        "avg_genomic_span": round(g["gene_len"].mean(), 2) if gene_count else "",
        "avg_sequence_length": (
            round(canon["exon_len"].sum() / gene_count, 2) if gene_count else ""
        ),
        "shortest_gene": int(g["gene_len"].min()) if gene_count else "",
        "longest_gene": int(g["gene_len"].max()) if gene_count else "",
        "total_tx": total_tx,
        "tx_per_gene": (
            round(total_tx / gene_count, 2) if gene_count else ""
        ),
        "total_exons": total_exons,
        "avg_exon_len": (
            round(total_exon_len / total_exons, 2) if total_exons else ""
        ),
        "avg_exons_per_tx": (
            round(total_exons / total_tx, 2) if total_tx else ""
        ),
        "total_introns": total_introns,
        "avg_intron_len": (
            round(total_intron_len / total_introns, 2) if total_introns else ""
        ),
    }

    if include_cds:
        tx_with_cds = set(cds["transcript_id"].dropna())
        coding_txs = txs[txs["transcript_id"].isin(tx_with_cds)]
        coding_tx_count = len(coding_txs)
        total_cds_segments = len(cds)
        total_cds_len = int(cds["cds_len"].sum())

        metrics.update(
            {
                "coding_tx": coding_tx_count,
                "coding_tx_per_gene": (
                    round(coding_tx_count / gene_count, 2) if gene_count else ""
                ),
                "total_cds_segments": total_cds_segments,
                "avg_cds_len": (
                    round(total_cds_len / coding_tx_count, 2)
                    if coding_tx_count
                    else ""
                ),
                "avg_cds_segment_len": (
                    round(total_cds_len / total_cds_segments, 2)
                    if total_cds_segments
                    else ""
                ),
                "avg_cds_segments_per_coding_tx": (
                    round(total_cds_segments / coding_tx_count, 2)
                    if coding_tx_count
                    else ""
                ),
            }
        )

    return metrics, canon


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_coding_stats(gff) -> Tuple[dict, int]:
    """
    Compute summary statistics for coding genes from a pyranges GFF3 object.

    Mirrors the metrics produced by ``leannes_stats.process_coding_genes``.
    Coding genes are defined as genes with ``biotype=protein_coding`` or
    genes that carry at least one CDS feature (excluding pseudogenes).

    The canonical transcript per gene is selected by preferring the
    ``Ensembl_canonical``-tagged transcript and breaking ties by CDS then
    exon spliced length.

    Parameters
    ----------
    gff:
        PyRanges object returned by ``parsers.annotation.parse_annotation``.

    Returns
    -------
    tuple of (stats_dict, total_coding_sequence_length)
        *stats_dict* contains the metrics listed below.
        *total_coding_sequence_length* is the sum of CDS spliced lengths
        across the canonical coding transcript of each coding gene (0 if
        the canonical transcript has no CDS).

    Stats keys
    ----------
    Coding genes, Average genomic span, Average sequence length,
    Average CDS length, Shortest gene, Longest gene, Total transcripts,
    Coding transcripts, Transcripts per gene, Coding transcripts per gene,
    Total exons, Total coding exons, Average exon length,
    Average coding exon length, Average exons per transcript,
    Average coding exons per coding transcript, Total introns,
    Average intron length.
    """
    gene_df, tx_df, exon_df, cds_df = _prepare_dataframes(gff)
    m, canon = _stats_for_group(gene_df, tx_df, exon_df, cds_df, "coding", include_cds=True)

    gene_count = m["gene_count"]

    if gene_count == 0:
        return (
            {
                "Coding genes": 0,
                "Average genomic span": "",
                "Average sequence length": "",
                "Average CDS length": "",
                "Shortest gene": "",
                "Longest gene": "",
                "Total transcripts": "",
                "Coding transcripts": "",
                "Transcripts per gene": "",
                "Coding transcripts per gene": "",
                "Total exons": "",
                "Total coding exons": "",
                "Average exon length": "",
                "Average coding exon length": "",
                "Average exons per transcript": "",
                "Average coding exons per coding transcript": "",
                "Total introns": "",
                "Average intron length": "",
            },
            0,
        )

    total_coding_sequence_length = int(canon["cds_len"].sum())

    stats = {
        "Coding genes": gene_count,
        "Average genomic span": m["avg_genomic_span"],
        "Average sequence length": m["avg_sequence_length"],
        "Average CDS length": m["avg_cds_len"],
        "Shortest gene": m["shortest_gene"],
        "Longest gene": m["longest_gene"],
        "Total transcripts": m["total_tx"],
        "Coding transcripts": m["coding_tx"],
        "Transcripts per gene": m["tx_per_gene"],
        "Coding transcripts per gene": m["coding_tx_per_gene"],
        "Total exons": m["total_exons"],
        "Total coding exons": m["total_cds_segments"],
        "Average exon length": m["avg_exon_len"],
        "Average coding exon length": m["avg_cds_segment_len"],
        "Average exons per transcript": m["avg_exons_per_tx"],
        "Average coding exons per coding transcript": m["avg_cds_segments_per_coding_tx"],
        "Total introns": m["total_introns"],
        "Average intron length": m["avg_intron_len"],
    }
    return stats, total_coding_sequence_length


def compute_noncoding_stats(gff) -> dict:
    """
    Compute summary statistics for non-coding genes from a pyranges GFF3 object.

    Mirrors the metrics produced by ``leannes_stats.process_non_coding_genes``.
    Non-coding genes are sub-classified into small non-coding, long non-coding,
    and miscellaneous non-coding categories based on gene biotype.

    Parameters
    ----------
    gff:
        PyRanges object returned by ``parsers.annotation.parse_annotation``.

    Returns
    -------
    dict
        Stats keys: Non-coding genes, Small non-coding genes,
        Long non-coding genes, Misc non-coding genes,
        Average genomic span, Average sequence length, Shortest gene,
        Longest gene, Total transcripts, Transcripts per gene,
        Total exons, Average exon length, Average exons per transcript,
        Total introns, Average intron length.
    """
    gene_df, tx_df, exon_df, cds_df = _prepare_dataframes(gff)

    nc_groups = {"lnoncoding", "snoncoding", "mnoncoding"}
    nc_genes = gene_df[gene_df["group"].isin(nc_groups)]
    gene_count = len(nc_genes)

    if gene_count == 0:
        return {
            "Non-coding genes": 0,
            "Small non-coding genes": "",
            "Long non-coding genes": "",
            "Misc non-coding genes": "",
            "Average genomic span": "",
            "Average sequence length": "",
            "Shortest gene": "",
            "Longest gene": "",
            "Total transcripts": "",
            "Transcripts per gene": "",
            "Total exons": "",
            "Average exon length": "",
            "Average exons per transcript": "",
            "Total introns": "",
            "Average intron length": "",
        }

    sn = int((nc_genes["group"] == "snoncoding").sum())
    ln = int((nc_genes["group"] == "lnoncoding").sum())
    mn = int((nc_genes["group"] == "mnoncoding").sum())

    nc_gene_ids = set(nc_genes["gene_id"])
    txs = tx_df[tx_df["gene_id"].isin(nc_gene_ids)]
    exons = exon_df[exon_df["gene_id"].isin(nc_gene_ids)]
    canon = _canonical_info(txs, exons, cds_df.iloc[:0])

    total_tx = len(txs)
    total_exons = len(exons)
    total_exon_len = int(exons["exon_len"].sum())

    introns = _compute_introns(exons)
    total_introns = len(introns)
    total_intron_len = int(introns["intron_len"].sum()) if not introns.empty else 0

    return {
        "Non-coding genes": gene_count,
        "Small non-coding genes": sn,
        "Long non-coding genes": ln,
        "Misc non-coding genes": mn,
        "Average genomic span": round(nc_genes["gene_len"].mean(), 2),
        "Average sequence length": (
            round(canon["exon_len"].sum() / gene_count, 2) if gene_count else ""
        ),
        "Shortest gene": int(nc_genes["gene_len"].min()),
        "Longest gene": int(nc_genes["gene_len"].max()),
        "Total transcripts": total_tx,
        "Transcripts per gene": round(total_tx / gene_count, 2) if gene_count else "",
        "Total exons": total_exons,
        "Average exon length": (
            round(total_exon_len / total_exons, 2) if total_exons else ""
        ),
        "Average exons per transcript": (
            round(total_exons / total_tx, 2) if total_tx else ""
        ),
        "Total introns": total_introns,
        "Average intron length": (
            round(total_intron_len / total_introns, 2) if total_introns else ""
        ),
    }


def compute_transcript_utr_stats(gff) -> pd.DataFrame:
    """
    Compute per-transcript 5' and 3' UTR metrics from a pyranges GFF3 object.

    UTR segments are taken directly from ``five_prime_UTR`` and
    ``three_prime_UTR`` features in the GFF3.  Junction counts equal the
    number of introns within each UTR region, i.e. the number of gaps
    between consecutive UTR segments for a transcript (a transcript with
    *n* non-adjacent UTR segments has *n-1* UTR junctions).

    Parameters
    ----------
    gff:
        PyRanges object returned by ``parsers.annotation.parse_annotation``.

    Returns
    -------
    pd.DataFrame
        One row per transcript with columns:

        * ``transcript_id``
        * ``gene_id``
        * ``five_prime_utr_length``  – total spliced 5' UTR length (bp)
        * ``three_prime_utr_length`` – total spliced 3' UTR length (bp)
        * ``five_prime_utr_junctions``  – intron count within 5' UTR
        * ``three_prime_utr_junctions`` – intron count within 3' UTR

        Transcripts with no UTR annotation receive 0 for the relevant
        columns.
    """
    tx_df = (
        gff[gff["Feature"].isin(_TRANSCRIPT_FEATURE_TYPES)]
        .reset_index(drop=True)
        .copy()
    )
    tx_df["gene_id"] = (
        tx_df["Parent"].str.split(",").str[0].str.removeprefix("gene:")
    )
    base = tx_df[["transcript_id", "gene_id"]].copy()

    def _utr_metrics(feature_name: str):
        utr = gff[gff["Feature"] == feature_name].reset_index(drop=True).copy()
        if utr.empty:
            return pd.Series(dtype="int64"), pd.Series(dtype="int64")

        utr["transcript_id"] = (
            utr["Parent"].str.split(",").str[0].str.removeprefix("transcript:")
        )
        utr["seg_len"] = utr["End"] - utr["Start"]

        lengths = utr.groupby("transcript_id")["seg_len"].sum()
        introns = _compute_introns(utr)
        junctions = (
            introns.groupby("transcript_id").size()
            if not introns.empty
            else pd.Series(dtype="int64")
        )
        return lengths, junctions

    utr5_len, utr5_junc = _utr_metrics("five_prime_UTR")
    utr3_len, utr3_junc = _utr_metrics("three_prime_UTR")

    base["five_prime_utr_length"] = (
        base["transcript_id"].map(utr5_len).fillna(0).astype(int)
    )
    base["three_prime_utr_length"] = (
        base["transcript_id"].map(utr3_len).fillna(0).astype(int)
    )
    base["five_prime_utr_junctions"] = (
        base["transcript_id"].map(utr5_junc).fillna(0).astype(int)
    )
    base["three_prime_utr_junctions"] = (
        base["transcript_id"].map(utr3_junc).fillna(0).astype(int)
    )
    return base.reset_index(drop=True)


def compute_pseudogene_stats(gff) -> dict:
    """
    Compute summary statistics for pseudogenes from a pyranges GFF3 object.

    Mirrors the metrics produced by ``leannes_stats.process_pseudogenes``.
    Pseudogenes are identified by a biotype containing the string
    ``pseudogene`` (case-insensitive).

    Parameters
    ----------
    gff:
        PyRanges object returned by ``parsers.annotation.parse_annotation``.

    Returns
    -------
    dict
        Stats keys: Pseudogenes, Average genomic span,
        Average sequence length, Shortest gene, Longest gene,
        Total transcripts, Transcripts per gene, Total exons,
        Average exon length, Average exons per transcript,
        Total introns, Average intron length.
    """
    gene_df, tx_df, exon_df, cds_df = _prepare_dataframes(gff)
    m, _ = _stats_for_group(gene_df, tx_df, exon_df, cds_df, "pseudogene")

    gene_count = m["gene_count"]

    if gene_count == 0:
        return {
            "Pseudogenes": 0,
            "Average genomic span": "",
            "Average sequence length": "",
            "Shortest gene": "",
            "Longest gene": "",
            "Total transcripts": "",
            "Transcripts per gene": "",
            "Total exons": "",
            "Average exon length": "",
            "Average exons per transcript": "",
            "Total introns": "",
            "Average intron length": "",
        }

    return {
        "Pseudogenes": gene_count,
        "Average genomic span": m["avg_genomic_span"],
        "Average sequence length": m["avg_sequence_length"],
        "Shortest gene": m["shortest_gene"],
        "Longest gene": m["longest_gene"],
        "Total transcripts": m["total_tx"],
        "Transcripts per gene": m["tx_per_gene"],
        "Total exons": m["total_exons"],
        "Average exon length": m["avg_exon_len"],
        "Average exons per transcript": m["avg_exons_per_tx"],
        "Total introns": m["total_introns"],
        "Average intron length": m["avg_intron_len"],
    }

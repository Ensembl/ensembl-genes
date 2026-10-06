"""
Normalise PyRanges1 GFF3/GTF tables into the pairwise-comparison schema.

Schema (one row per gene, transcript, exon or CDS; COMPARISON_COLUMNS):

    Chromosome            sequence name as written in the file (seqname maps are
                          applied later, by the runner)
    Start, End            zero-based half-open (Start = GFF start - 1, End = GFF end)
    Strand                "+", "-" or "."
    Feature               "gene", "transcript", "exon" or "CDS"
    source_feature        column 3 as written (e.g. mRNA, ncRNA_gene, pseudogene)
    gene_id               comparison key of the owning gene, unique genome-wide
    transcript_id         comparison key of the owning transcript ("" on genes)
    ID, Parent            source identifiers with Ensembl type prefixes
                          ("gene:", "transcript:", ...) removed. Children with
                          several parents are split into one row per parent, so
                          Parent always holds a single identifier.
    gene_biotype          propagated from the gene to its transcripts and children
    transcript_biotype    propagated from the transcript to its children
    tags                  transcript tags, comma-separated (e.g. Ensembl_canonical)
    original_gene_id      gene_id before cross-sequence namespacing
    original_transcript_id
    record_index          0-based order of the source record in the file

Identifiers reused on more than one sequence (Tiberius numbers genes g1, g2, ...
per chromosome) are namespaced as "<seqname>:<id>" in gene_id/transcript_id;
the file values are kept in original_gene_id/original_transcript_id.

Relationships that cannot be resolved (child without a known parent, transcript
without a gene, duplicate gene/transcript identifiers on one sequence) raise
ValueError rather than being guessed from overlapping coordinates.
"""

from collections import Counter

import numpy as np
import pandas as pd

GENE_FEATURES = frozenset({"gene", "ncRNA_gene", "pseudogene"})
TRANSCRIPT_FEATURES = frozenset(
    {
        "mRNA",
        "transcript",
        "lnc_RNA",
        "ncRNA",
        "rRNA",
        "tRNA",
        "snRNA",
        "snoRNA",
        "miRNA",
        "pre_miRNA",
        "SRP_RNA",
        "RNase_P_RNA",
        "RNase_MRP_RNA",
        "telomerase_RNA",
        "scRNA",
        "processed_transcript",
        "V_gene_segment",
        "D_gene_segment",
        "J_gene_segment",
        "C_gene_segment",
        "pseudogenic_transcript",
    }
)
CHILD_FEATURES = frozenset({"exon", "CDS"})

COMPARISON_COLUMNS = [
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Feature",
    "source_feature",
    "gene_id",
    "transcript_id",
    "ID",
    "Parent",
    "gene_biotype",
    "transcript_biotype",
    "tags",
    "original_gene_id",
    "original_transcript_id",
    "record_index",
]

_ENSEMBL_PREFIX = r"^(?:gene|transcript|chromosome|mRNA):"
_MAX_EXAMPLES = 5


def _attr(df: pd.DataFrame, name: str) -> pd.Series:
    """Return an attribute column as strings, with "" for absent values."""
    if name not in df.columns:
        return pd.Series("", index=df.index, dtype=object)
    return df[name].astype(object).where(df[name].notna(), "").astype(str).str.strip()


def _first_non_empty(*series: pd.Series) -> pd.Series:
    result = series[0]
    for other in series[1:]:
        result = result.where(result != "", other)
    return result


def _fail(message: str, examples) -> None:
    examples = list(examples)
    shown = ", ".join(str(e) for e in examples[:_MAX_EXAMPLES])
    more = (
        f" (+{len(examples) - _MAX_EXAMPLES} more)"
        if len(examples) > _MAX_EXAMPLES
        else ""
    )
    raise ValueError(f"{message}: {shown}{more}")


def _base_table(raw: pd.DataFrame) -> tuple[pd.DataFrame, Counter]:
    """Keep gene/transcript/exon/CDS records and add the fixed columns."""
    source_feature = raw["Feature"].astype(str)
    kind = pd.Series("", index=raw.index, dtype=object)
    kind[source_feature.isin(GENE_FEATURES)] = "gene"
    kind[source_feature.isin(TRANSCRIPT_FEATURES)] = "transcript"
    kind[source_feature.isin(CHILD_FEATURES)] = source_feature
    dropped = Counter(source_feature[kind == ""])

    base = pd.DataFrame(
        {
            "Chromosome": raw["Chromosome"].astype(str),
            "Start": raw["Start"].astype("int64"),
            "End": raw["End"].astype("int64"),
            "Strand": raw["Strand"].astype(str),
            "Feature": kind,
            "source_feature": source_feature,
            "record_index": np.arange(len(raw), dtype="int64"),
        },
        index=raw.index,
    )
    keep = kind != ""
    return base[keep].copy(), dropped


def _check_unique(df: pd.DataFrame, key: str, what: str) -> None:
    dup = df.duplicated(["Chromosome", key], keep=False)
    if dup.any():
        _fail(
            f"Duplicate {what} identifiers on the same sequence",
            sorted(
                {
                    f"{c}:{k}"
                    for c, k in df.loc[dup, ["Chromosome", key]].itertuples(index=False)
                }
            ),
        )


def _namespace_duplicates(df: pd.DataFrame, diagnostics: dict) -> pd.DataFrame:
    """Namespace gene_id/transcript_id values that occur on more than one sequence."""
    df["original_gene_id"] = df["gene_id"]
    df["original_transcript_id"] = df["transcript_id"]
    for column, label in (("gene_id", "gene"), ("transcript_id", "transcript")):
        present = df[df[column] != ""]
        per_id = present.groupby(column)["Chromosome"].nunique()
        reused = set(per_id.index[per_id > 1])
        diagnostics[f"namespaced_{label}_ids"] = len(reused)
        if reused:
            mask = df[column].isin(reused)
            df.loc[mask, column] = (
                df.loc[mask, "Chromosome"] + ":" + df.loc[mask, column]
            )
    return df


def _finish(df: pd.DataFrame, diagnostics: dict) -> pd.DataFrame:
    """Namespace IDs, propagate biotypes, order rows and record counts."""
    df = _namespace_duplicates(df, diagnostics)
    genes = df[df["Feature"] == "gene"]
    transcripts = df[df["Feature"] == "transcript"]
    gene_biotype = genes.set_index("gene_id")["gene_biotype"]
    tx_biotype = transcripts.set_index("transcript_id")["transcript_biotype"]
    not_gene = df["Feature"] != "gene"
    df.loc[not_gene, "gene_biotype"] = (
        df.loc[not_gene, "gene_id"].map(gene_biotype).fillna("")
    )
    child = df["Feature"].isin(CHILD_FEATURES)
    df.loc[child, "transcript_biotype"] = (
        df.loc[child, "transcript_id"].map(tx_biotype).fillna("")
    )

    df = df.sort_values("record_index", kind="stable").reset_index(drop=True)
    # Records created for GTF files without gene/transcript lines have no source line.
    synthetic = df["record_index"] % 1 != 0
    df["record_index"] = df["record_index"].where(~synthetic, -1).astype("int64")
    df = df[COMPARISON_COLUMNS]

    counts = df["Feature"].value_counts()
    tx_with_exons = set(df.loc[df["Feature"] == "exon", "transcript_id"])
    genes_with_tx = set(df.loc[df["Feature"] == "transcript", "gene_id"])
    all_tx = df.loc[df["Feature"] == "transcript", "transcript_id"]
    all_genes = df.loc[df["Feature"] == "gene", "gene_id"]
    diagnostics.update(
        {
            "genes": int(counts.get("gene", 0)),
            "transcripts": int(counts.get("transcript", 0)),
            "exon_rows": int(counts.get("exon", 0)),
            "cds_rows": int(counts.get("CDS", 0)),
            "transcripts_without_exons": int((~all_tx.isin(tx_with_exons)).sum()),
            "genes_without_transcripts": int((~all_genes.isin(genes_with_tx)).sum()),
        }
    )
    df.attrs["parse_diagnostics"] = diagnostics
    return df


def normalise_gff3(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise a PyRanges1 GFF3 table (ID=/Parent= hierarchy).
    Args:
            raw: DataFrame from pyranges1.read_gff3, rows in file order
    Returns:
            DataFrame with COMPARISON_COLUMNS
    """
    df, dropped = _base_table(raw)
    raw = raw.loc[df.index]
    diagnostics = {
        "records_read": len(df) + sum(dropped.values()),
        "dropped_feature_types": dict(dropped),
    }

    df["ID"] = _attr(raw, "ID").str.replace(_ENSEMBL_PREFIX, "", regex=True)
    parent_raw = _attr(raw, "Parent")
    biotype = _attr(raw, "biotype")
    is_gene = df["Feature"] == "gene"
    is_tx = df["Feature"] == "transcript"
    df["gene_biotype"] = _first_non_empty(_attr(raw, "gene_biotype"), biotype).where(
        is_gene, ""
    )
    df["transcript_biotype"] = _first_non_empty(
        _attr(raw, "transcript_biotype"), biotype
    ).where(is_tx, "")
    df["tags"] = _attr(raw, "tag").where(is_tx, "")

    # Genes: gene_id attribute if present, otherwise ID.
    genes = df[is_gene].copy()
    genes["gene_id"] = _first_non_empty(_attr(raw, "gene_id")[is_gene], genes["ID"])
    if (genes["gene_id"] == "").any():
        _fail(
            "Gene records without ID",
            genes.loc[genes["gene_id"] == "", "record_index"] + 1,
        )
    _check_unique(genes, "ID", "gene")
    genes["transcript_id"] = ""
    genes["Parent"] = ""
    gene_key = genes.set_index(["Chromosome", "ID"])["gene_id"]

    # Transcripts: exactly one parent gene.
    tx = df[is_tx].copy()
    tx_parents = parent_raw[is_tx].str.split(",")
    multi = tx_parents.str.len() > 1
    if multi.any():
        _fail("Transcripts with more than one parent gene", tx.loc[multi, "ID"])
    tx["Parent"] = tx_parents.str[0].str.replace(_ENSEMBL_PREFIX, "", regex=True)
    tx["transcript_id"] = _first_non_empty(_attr(raw, "transcript_id")[is_tx], tx["ID"])
    if (tx["transcript_id"] == "").any():
        _fail(
            "Transcript records without ID",
            tx.loc[tx["transcript_id"] == "", "record_index"] + 1,
        )
    _check_unique(tx, "ID", "transcript")
    tx["gene_id"] = pd.Series(
        gene_key.reindex(
            pd.MultiIndex.from_arrays([tx["Chromosome"], tx["Parent"]])
        ).to_numpy(),
        index=tx.index,
    )
    orphan = tx["gene_id"].isna()
    if orphan.any():
        _fail(
            "Transcripts whose Parent is not a gene record on the same sequence",
            [
                f"{t} (Parent={p})"
                for t, p in tx.loc[orphan, ["ID", "Parent"]].itertuples(index=False)
            ],
        )
    tx_key = tx.set_index(["Chromosome", "ID"])[["transcript_id", "gene_id"]]

    # Exons/CDS: one row per parent transcript.
    children = df[df["Feature"].isin(CHILD_FEATURES)].copy()
    children["Parent"] = parent_raw[children.index].str.split(",")
    diagnostics["multi_parent_children"] = int((children["Parent"].str.len() > 1).sum())
    children = children.explode("Parent")
    children["Parent"] = (
        children["Parent"].str.strip().str.replace(_ENSEMBL_PREFIX, "", regex=True)
    )
    no_parent = children["Parent"] == ""
    if no_parent.any():
        _fail(
            "Exon/CDS records without Parent",
            children.loc[no_parent, "record_index"] + 1,
        )
    resolved = tx_key.reindex(
        pd.MultiIndex.from_arrays([children["Chromosome"], children["Parent"]])
    )
    children["transcript_id"] = resolved["transcript_id"].to_numpy()
    children["gene_id"] = resolved["gene_id"].to_numpy()
    orphan = children["transcript_id"].isna()
    if orphan.any():
        _fail(
            "Exon/CDS records whose Parent is not a transcript record on the same sequence",
            sorted(set(children.loc[orphan, "Parent"])),
        )

    return _finish(pd.concat([genes, tx, children]), diagnostics)


def normalise_gtf(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise a PyRanges1 GTF table (gene_id/transcript_id attributes).

    Gene and transcript records are optional in GTF; missing ones are created from
    the span of their children and counted in the parse diagnostics.
    Args:
            raw: DataFrame from pyranges1.read_gtf, rows in file order
    Returns:
            DataFrame with COMPARISON_COLUMNS
    """
    df, dropped = _base_table(raw)
    raw = raw.loc[df.index]
    diagnostics = {
        "records_read": len(df) + sum(dropped.values()),
        "dropped_feature_types": dict(dropped),
    }

    df["gene_id"] = _attr(raw, "gene_id")
    df["transcript_id"] = _attr(raw, "transcript_id")
    is_gene = df["Feature"] == "gene"
    is_tx = df["Feature"] == "transcript"
    df.loc[is_gene, "transcript_id"] = ""
    df["gene_biotype"] = _first_non_empty(
        _attr(raw, "gene_biotype"),
        _attr(raw, "gene_type"),
        _attr(raw, "biotype").where(is_gene, ""),
    )
    df["transcript_biotype"] = _first_non_empty(
        _attr(raw, "transcript_biotype"),
        _attr(raw, "transcript_type"),
        _attr(raw, "biotype").where(is_tx, ""),
    )
    df["tags"] = _attr(raw, "tag").where(is_tx, "")

    missing = (df["gene_id"] == "") & (is_gene | is_tx)
    missing |= (df["transcript_id"] == "") & ~is_gene
    if missing.any():
        _fail(
            "GTF records without gene_id/transcript_id (line order)",
            df.loc[missing, "record_index"] + 1,
        )

    genes = df[is_gene].copy()
    tx = df[is_tx].copy()
    children = df[df["Feature"].isin(CHILD_FEATURES)].copy()
    _check_unique(genes, "gene_id", "gene")
    _check_unique(tx, "transcript_id", "transcript")

    # Child gene_id may be omitted; otherwise it must agree with its transcript.
    tx_gene = tx.set_index(["Chromosome", "transcript_id"])["gene_id"]
    child_key = pd.MultiIndex.from_arrays(
        [children["Chromosome"], children["transcript_id"]]
    )
    known_gene = pd.Series(tx_gene.reindex(child_key).to_numpy(), index=children.index)
    conflict = (
        known_gene.notna()
        & (children["gene_id"] != "")
        & (children["gene_id"] != known_gene)
    )
    if conflict.any():
        _fail(
            "Exon/CDS gene_id disagrees with its transcript record",
            children.loc[conflict, "transcript_id"],
        )
    children["gene_id"] = children["gene_id"].where(
        children["gene_id"] != "", known_gene
    )
    if children["gene_id"].isna().any() or (children["gene_id"] == "").any():
        _fail(
            "Exon/CDS records without gene_id and without a transcript record",
            children.loc[
                children["gene_id"].isna() | (children["gene_id"] == ""),
                "transcript_id",
            ],
        )

    # Create transcript/gene records that the GTF omits.
    new_tx = _span_records(
        children[~pd.Series(child_key.isin(tx_gene.index), index=children.index)],
        ["Chromosome", "transcript_id", "gene_id"],
        "transcript",
    )
    tx = pd.concat([tx, new_tx])
    gene_ids = genes.set_index(["Chromosome", "gene_id"]).index
    tx_key = pd.MultiIndex.from_arrays([tx["Chromosome"], tx["gene_id"]])
    new_genes = _span_records(
        tx[~tx_key.isin(gene_ids)], ["Chromosome", "gene_id"], "gene"
    )
    new_genes["transcript_id"] = ""
    genes = pd.concat([genes, new_genes])
    diagnostics["synthesised_transcripts"] = len(new_tx)
    diagnostics["synthesised_genes"] = len(new_genes)

    tx_gene = tx.set_index(["Chromosome", "transcript_id"])["gene_id"]
    if tx_gene.index.duplicated().any():
        _fail(
            "Transcript identifiers assigned to more than one gene",
            tx_gene.index[tx_gene.index.duplicated()],
        )

    genes["ID"] = genes["gene_id"]
    genes["Parent"] = ""
    tx["ID"] = tx["transcript_id"]
    tx["Parent"] = tx["gene_id"]
    children["ID"] = ""
    children["Parent"] = children["transcript_id"]
    diagnostics["multi_parent_children"] = 0

    # GTF often repeats gene_biotype on every line; take the gene's first value.
    combined = pd.concat([genes, tx, children])
    first_gene_biotype = (
        combined[combined["gene_biotype"] != ""]
        .groupby(["Chromosome", "gene_id"])["gene_biotype"]
        .first()
    )
    gene_key = pd.MultiIndex.from_arrays([combined["Chromosome"], combined["gene_id"]])
    combined["gene_biotype"] = (
        pd.Series(
            first_gene_biotype.reindex(gene_key).to_numpy(), index=combined.index
        ).fillna("")
    ).where(combined["Feature"] == "gene", combined["gene_biotype"])
    return _finish(combined, diagnostics)


def _span_records(rows: pd.DataFrame, keys: list[str], feature: str) -> pd.DataFrame:
    """Create one record per key group spanning its rows, placed at its first row."""
    if rows.empty:
        return rows.iloc[0:0].assign(Feature=feature)
    spans = (
        rows.groupby(keys, sort=False)
        .agg(
            Start=("Start", "min"),
            End=("End", "max"),
            Strand=("Strand", "first"),
            record_index=("record_index", "min"),
            gene_biotype=("gene_biotype", "first"),
            transcript_biotype=("transcript_biotype", "first"),
        )
        .reset_index()
    )
    spans["Feature"] = feature
    spans["source_feature"] = f"inferred_{feature}"
    spans["tags"] = ""
    # Place the synthetic parent just before its first child.
    spans["record_index"] = spans["record_index"].astype("float64") - 0.5
    return spans


def normalise_tiberius(raw: pd.DataFrame, raw_columns: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise the Tiberius hybrid: bare IDs on gene/transcript rows, GTF children.

    pyranges1.read_gtf cannot read bare column-9 identifiers, so they are taken
    from ``raw_columns`` (the file's records read verbatim). Rows are aligned by
    position and every aligned row must agree on sequence, type and coordinates.
    A transcript's gene is the text before its last "." (g1.t1 -> g1).
    Args:
            raw: DataFrame from pyranges1.read_gtf
            raw_columns: DataFrame with Chromosome, Feature, gff_start, End, attributes
    Returns:
            DataFrame with COMPARISON_COLUMNS
    """
    raw = raw.reset_index(drop=True).copy()
    raw_columns = raw_columns.reset_index(drop=True)
    aligned = len(raw) == len(raw_columns) and (
        (raw["Chromosome"].astype(str) == raw_columns["Chromosome"]).all()
        and (raw["Feature"].astype(str) == raw_columns["Feature"]).all()
        and (raw["Start"] + 1 == raw_columns["gff_start"]).all()
        and (raw["End"] == raw_columns["End"]).all()
    )
    if not aligned:
        raise ValueError(
            "Could not align Tiberius records with the PyRanges1 GTF table"
        )

    bare = raw_columns["attributes"].str.strip().str.rstrip(";").str.strip()
    is_bare = ~bare.str.contains('"', regex=False) & ~bare.str.contains(
        "=", regex=False
    )
    for column in ("gene_id", "transcript_id"):
        if column not in raw.columns:
            raw[column] = pd.NA
        raw[column] = raw[column].astype(object)

    gene_rows = (raw["Feature"].astype(str) == "gene") & is_bare
    raw.loc[gene_rows, "gene_id"] = bare[gene_rows]
    tx_rows = (raw["Feature"].astype(str) == "transcript") & is_bare
    tids = bare[tx_rows]
    dot = tids.str.rfind(".")
    raw.loc[tx_rows, "transcript_id"] = tids
    raw.loc[tx_rows, "gene_id"] = [t[:d] if d > 0 else t for t, d in zip(tids, dot)]
    return normalise_gtf(raw)

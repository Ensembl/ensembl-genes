"""Normalized source-model feature tables for expected-centric audits."""

# pylint: disable=missing-function-docstring,unnecessary-lambda,too-many-locals

from __future__ import annotations

import hashlib

import pandas as pd


def _string_value(value: object, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    text = str(value)
    if text.lower() == "nan":
        return default
    return text


def _int_value(value: object, default: int = 0) -> int:
    if value is None or pd.isna(value):
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _normalise_strand(value: object) -> int:
    if value in ("+", "+1", "1", 1):
        return 1
    if value in ("-", "-1", -1):
        return -1
    return _int_value(value, 0)


def _hash_key(parts: list[str]) -> str:
    if not parts:
        return ""
    return hashlib.sha1("|".join(parts).encode("utf-8")).hexdigest()


def normalise_source_genes(
    genes: pd.DataFrame,
    *,
    source_name: str,
    source_role: str,
    source_class: str = "",
) -> pd.DataFrame:
    columns = [
        "source_name",
        "source_role",
        "source_class",
        "gene_id",
        "stable_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "biotype",
        "logic_name",
    ]
    out = genes.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = pd.NA
    out["source_name"] = source_name
    out["source_role"] = source_role
    out["source_class"] = source_class or source_role
    if not out.empty:
        out["gene_id"] = out["gene_id"].map(lambda x: _string_value(x))
        out["stable_id"] = out["stable_id"].map(lambda x: _string_value(x))
        out["seq_region_name"] = out["seq_region_name"].map(lambda x: _string_value(x))
        out["seq_region_start"] = out["seq_region_start"].map(lambda x: _int_value(x))
        out["seq_region_end"] = out["seq_region_end"].map(lambda x: _int_value(x))
        out["seq_region_strand"] = out["seq_region_strand"].map(_normalise_strand)
        out["biotype"] = out["biotype"].map(lambda x: _string_value(x, "unknown"))
        out["logic_name"] = out["logic_name"].map(lambda x: _string_value(x, ""))
    return out[columns]


def normalise_source_transcripts(
    transcripts: pd.DataFrame,
    translations: pd.DataFrame | None = None,
    *,
    source_name: str,
    source_role: str,
    source_class: str = "",
) -> pd.DataFrame:
    columns = [
        "source_name",
        "source_role",
        "source_class",
        "transcript_id",
        "gene_id",
        "stable_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "biotype",
        "logic_name",
        "has_translation",
    ]
    out = transcripts.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = pd.NA
    tr_ids = set()
    if (
        translations is not None
        and not translations.empty
        and "transcript_id" in translations
    ):
        tr_ids = set(
            translations["transcript_id"].dropna().map(lambda x: _string_value(x))
        )
    out["source_name"] = source_name
    out["source_role"] = source_role
    out["source_class"] = source_class or source_role
    if not out.empty:
        out["transcript_id"] = out["transcript_id"].map(lambda x: _string_value(x))
        out["gene_id"] = out["gene_id"].map(lambda x: _string_value(x))
        out["stable_id"] = out["stable_id"].map(lambda x: _string_value(x))
        out["seq_region_name"] = out["seq_region_name"].map(lambda x: _string_value(x))
        out["seq_region_start"] = out["seq_region_start"].map(lambda x: _int_value(x))
        out["seq_region_end"] = out["seq_region_end"].map(lambda x: _int_value(x))
        out["seq_region_strand"] = out["seq_region_strand"].map(_normalise_strand)
        out["biotype"] = out["biotype"].map(lambda x: _string_value(x, "unknown"))
        out["logic_name"] = out["logic_name"].map(lambda x: _string_value(x, ""))
        out["has_translation"] = out["transcript_id"].map(lambda x: int(x in tr_ids))
    return out[columns]


def normalise_source_features(
    features: pd.DataFrame,
    *,
    feature_kind: str,
    source_name: str,
    source_role: str,
    source_class: str = "",
) -> pd.DataFrame:
    columns = [
        "source_name",
        "source_role",
        "source_class",
        "feature_kind",
        "feature_id",
        "transcript_id",
        "gene_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "rank",
        "phase",
        "end_phase",
    ]
    rename = {
        "exon_id": "feature_id",
        "cds_id": "feature_id",
        "exon_rank": "rank",
        "cds_rank": "rank",
    }
    out = features.rename(
        columns={k: v for k, v in rename.items() if k in features}
    ).copy()
    for col in columns:
        if col not in out.columns:
            out[col] = pd.NA
    out["source_name"] = source_name
    out["source_role"] = source_role
    out["source_class"] = source_class or source_role
    out["feature_kind"] = feature_kind
    if not out.empty:
        out["feature_id"] = out.apply(
            lambda r: _string_value(r.get("feature_id")) or f"{feature_kind}:{r.name}",
            axis=1,
        )
        out["transcript_id"] = out["transcript_id"].map(lambda x: _string_value(x))
        out["gene_id"] = out["gene_id"].map(lambda x: _string_value(x))
        out["seq_region_name"] = out["seq_region_name"].map(lambda x: _string_value(x))
        out["seq_region_start"] = out["seq_region_start"].map(lambda x: _int_value(x))
        out["seq_region_end"] = out["seq_region_end"].map(lambda x: _int_value(x))
        out["seq_region_strand"] = out["seq_region_strand"].map(_normalise_strand)
        out["rank"] = out["rank"].map(lambda x: _int_value(x))
        out["phase"] = out["phase"].map(lambda x: _int_value(x, -1))
        out["end_phase"] = out["end_phase"].map(lambda x: _int_value(x, -1))
    return out[columns]


def build_introns(exons: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "source_name",
        "source_role",
        "source_class",
        "transcript_id",
        "gene_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "rank",
        "intron_key",
    ]
    if exons.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict] = []
    for tx_id, group in exons.groupby("transcript_id", dropna=False):
        ordered = group.sort_values("rank")
        previous = None
        rank = 0
        for _, exon in ordered.iterrows():
            if previous is None:
                previous = exon
                continue
            start = _int_value(previous.get("seq_region_end")) + 1
            end = _int_value(exon.get("seq_region_start")) - 1
            if end >= start:
                rank += 1
                seq = _string_value(exon.get("seq_region_name"))
                strand = _normalise_strand(exon.get("seq_region_strand"))
                key = f"{seq}:{strand}:{start}:{end}"
                rows.append(
                    {
                        "source_name": exon.get("source_name"),
                        "source_role": exon.get("source_role"),
                        "source_class": exon.get("source_class"),
                        "transcript_id": _string_value(tx_id),
                        "gene_id": _string_value(exon.get("gene_id")),
                        "seq_region_name": seq,
                        "seq_region_start": start,
                        "seq_region_end": end,
                        "seq_region_strand": strand,
                        "rank": rank,
                        "intron_key": key,
                    }
                )
            previous = exon
    return pd.DataFrame(rows, columns=columns)


def build_intron_chains(
    introns: pd.DataFrame, transcripts: pd.DataFrame
) -> pd.DataFrame:
    columns = [
        "source_name",
        "source_role",
        "source_class",
        "transcript_id",
        "gene_id",
        "intron_count",
        "ordered_intron_chain_key",
        "intron_chain_hash",
    ]
    rows: list[dict] = []
    intron_groups = {
        tx_id: group.sort_values("rank")
        for tx_id, group in introns.groupby("transcript_id", dropna=False)
    }
    for _, tx in transcripts.iterrows():
        tx_id = _string_value(tx.get("transcript_id"))
        group = intron_groups.get(tx_id)
        keys = [] if group is None else group["intron_key"].map(_string_value).tolist()
        rows.append(
            {
                "source_name": tx.get("source_name"),
                "source_role": tx.get("source_role"),
                "source_class": tx.get("source_class"),
                "transcript_id": tx_id,
                "gene_id": _string_value(tx.get("gene_id")),
                "intron_count": len(keys),
                "ordered_intron_chain_key": ",".join(keys),
                "intron_chain_hash": _hash_key(keys),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def build_model_quality(
    transcripts: pd.DataFrame,
    exons: pd.DataFrame,
    cds: pd.DataFrame,
    intron_chains: pd.DataFrame,
) -> pd.DataFrame:
    columns = [
        "source_name",
        "source_role",
        "source_class",
        "transcript_id",
        "gene_id",
        "exon_count",
        "intron_count",
        "cds_exon_count",
        "has_cds",
        "has_translation",
        "cds_span",
        "protein_length",
        "start_complete",
        "stop_complete",
        "ordered_intron_chain_key",
        "intron_chain_hash",
    ]
    exon_counts = (
        exons.groupby("transcript_id").size().to_dict() if not exons.empty else {}
    )
    cds_counts = cds.groupby("transcript_id").size().to_dict() if not cds.empty else {}
    cds_span = {}
    if not cds.empty:
        cds_span = (
            cds.assign(_span=cds["seq_region_end"] - cds["seq_region_start"] + 1)
            .groupby("transcript_id")["_span"]
            .sum()
            .to_dict()
        )
    chains = (
        intron_chains.set_index("transcript_id")
        if not intron_chains.empty
        else pd.DataFrame()
    )
    rows: list[dict] = []
    for _, tx in transcripts.iterrows():
        tx_id = _string_value(tx.get("transcript_id"))
        chain = chains.loc[tx_id] if not chains.empty and tx_id in chains.index else {}
        rows.append(
            {
                "source_name": tx.get("source_name"),
                "source_role": tx.get("source_role"),
                "source_class": tx.get("source_class"),
                "transcript_id": tx_id,
                "gene_id": _string_value(tx.get("gene_id")),
                "exon_count": int(exon_counts.get(tx_id, 0)),
                "intron_count": int(
                    chain.get("intron_count", max(0, exon_counts.get(tx_id, 0) - 1))
                    if isinstance(chain, pd.Series)
                    else 0
                ),
                "cds_exon_count": int(cds_counts.get(tx_id, 0)),
                "has_cds": int(cds_counts.get(tx_id, 0) > 0),
                "has_translation": _int_value(tx.get("has_translation")),
                "cds_span": int(cds_span.get(tx_id, 0)),
                "protein_length": pd.NA,
                "start_complete": pd.NA,
                "stop_complete": pd.NA,
                "ordered_intron_chain_key": (
                    chain.get("ordered_intron_chain_key", "")
                    if isinstance(chain, pd.Series)
                    else ""
                ),
                "intron_chain_hash": (
                    chain.get("intron_chain_hash", "")
                    if isinstance(chain, pd.Series)
                    else ""
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def combine_source_models(
    sources: list[dict[str, pd.DataFrame | str]],
) -> dict[str, pd.DataFrame]:
    gene_tables = []
    transcript_tables = []
    exon_tables = []
    cds_tables = []
    for source in sources:
        name = str(source.get("name", "source"))
        role = str(source.get("role", "candidate"))
        source_class = str(source.get("class", role))
        genes = source.get("genes")
        transcripts = source.get("transcripts")
        exons = source.get("exons")
        cds = source.get("cds")
        translations = source.get("translations")
        if isinstance(genes, pd.DataFrame):
            gene_tables.append(
                normalise_source_genes(
                    genes, source_name=name, source_role=role, source_class=source_class
                )
            )
        if isinstance(transcripts, pd.DataFrame):
            transcript_tables.append(
                normalise_source_transcripts(
                    transcripts,
                    translations if isinstance(translations, pd.DataFrame) else None,
                    source_name=name,
                    source_role=role,
                    source_class=source_class,
                )
            )
        if isinstance(exons, pd.DataFrame):
            exon_tables.append(
                normalise_source_features(
                    exons,
                    feature_kind="exon",
                    source_name=name,
                    source_role=role,
                    source_class=source_class,
                )
            )
        if isinstance(cds, pd.DataFrame):
            cds_tables.append(
                normalise_source_features(
                    cds,
                    feature_kind="CDS",
                    source_name=name,
                    source_role=role,
                    source_class=source_class,
                )
            )

    genes_df = (
        pd.concat(gene_tables, ignore_index=True) if gene_tables else pd.DataFrame()
    )
    tx_df = (
        pd.concat(transcript_tables, ignore_index=True)
        if transcript_tables
        else pd.DataFrame()
    )
    exon_df = (
        pd.concat(exon_tables, ignore_index=True) if exon_tables else pd.DataFrame()
    )
    cds_df = pd.concat(cds_tables, ignore_index=True) if cds_tables else pd.DataFrame()
    intron_df = build_introns(exon_df)
    chain_df = build_intron_chains(intron_df, tx_df)
    quality_df = build_model_quality(tx_df, exon_df, cds_df, chain_df)
    return {
        "source_gene": genes_df,
        "source_transcript": tx_df,
        "source_exon": exon_df,
        "source_cds": cds_df,
        "source_intron": intron_df,
        "source_intron_chain": chain_df,
        "source_model_quality": quality_df,
    }

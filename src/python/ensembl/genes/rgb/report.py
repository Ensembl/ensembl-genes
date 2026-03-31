from __future__ import annotations

from typing import Dict, Tuple

import pandas as pd


def _interval_overlap(a_start: pd.Series, a_end: pd.Series, b_start: pd.Series, b_end: pd.Series) -> pd.Series:
    return (a_start <= b_end) & (a_end >= b_start)


def compute_retention_by_biotype(
    gene_map: pd.DataFrame,
    core_tx: pd.DataFrame,
    layer_tx: pd.DataFrame,
    layer_tr: pd.DataFrame,
    biotype_to_tier: Dict[str, Dict[str, object]] | Dict[str, Tuple[str, int]]
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (retention_by_biotype, retention_by_tier, crosswalk, by_logic) dataframes.

    A layer transcript is 'built' if it overlaps any core transcript on the same strand
    within the same strict locus (we map transcripts to locus via gene_id).
    """
    # Normalize mapping structure
    b2t: Dict[str, Tuple[str, int]] = {}
    for k, v in biotype_to_tier.items():
        if isinstance(v, dict):
            tier = str(v.get("tier") or v.get("layer") or v.get("id") or "")
            idx = int(v.get("tier_index") or 0)
            b2t[k] = (tier, idx)
        else:
            tier, idx = v
            b2t[k] = (str(tier), int(idx))

    # Map transcripts to strict locus via gene_id
    layer_tx2 = layer_tx.merge(
        gene_map[gene_map.db_kind == "layer"][['gene_id', 'locus_id_strict']],
        on="gene_id",
        how="left",
    )
    core_tx2 = core_tx.merge(
        gene_map[gene_map.db_kind == "core"][['gene_id', 'locus_id_strict']],
        on="gene_id",
        how="left",
    )

    # Build an index of core transcripts per locus for fast lookup
    core_by_loc = {k: g for k, g in core_tx2.groupby("locus_id_strict")}

    # Basic features
    layer_tx2["biotype"] = layer_tx2["biotype"].fillna("unknown").astype(str)
    layer_tx2["tx_len_bp"] = (layer_tx2["seq_region_end"].astype(int) - layer_tx2["seq_region_start"].astype(int) + 1)
    has_tr = set(layer_tr["transcript_id"].astype(int).tolist()) if not layer_tr.empty else set()
    layer_tx2["has_cds"] = layer_tx2["transcript_id"].astype(int).isin(has_tr).astype(int)

    built_flags = []
    core_biotypes = []
    for loc_id, grp in layer_tx2.groupby("locus_id_strict"):
        core_grp = core_by_loc.get(loc_id)
        if core_grp is None or core_grp.empty:
            built_flags.extend([0] * len(grp))
            core_biotypes.extend([None] * len(grp))
            continue
        for _, r in grp.iterrows():
            ov = _interval_overlap(
                pd.Series([int(r.seq_region_start)]*len(core_grp)),
                pd.Series([int(r.seq_region_end)]*len(core_grp)),
                core_grp["seq_region_start"].astype(int),
                core_grp["seq_region_end"].astype(int),
            )
            is_built = bool(ov.any())
            built_flags.append(1 if is_built else 0)
            core_biotypes.append(core_grp.loc[ov].iloc[0]["biotype"] if is_built else None)

    layer_tx2["built"] = built_flags
    layer_tx2["core_biotype"] = core_biotypes

    # Aggregation by layer biotype
    def agg(df: pd.DataFrame) -> pd.Series:
        n = len(df)
        nb = int(df["built"].sum())
        nl = n - nb
        built_pct = (nb / n) if n else 0.0
        left_cds_pct = (df.loc[df["built"] == 0, "has_cds"].mean() if nl else 0.0)
        med_len = float(df["tx_len_bp"].median()) if n else 0.0
        bt = str(df.name)
        tier, idx = b2t.get(bt, ("", 0))
        return pd.Series(
            {
                "layer_biotype": bt,
                "tier": tier,
                "tier_index": int(idx),
                "n_layer_tx": int(n),
                "n_built_tx": int(nb),
                "built_pct": float(built_pct),
                "n_leftout_tx": int(nl),
                "leftout_cds_pct": float(left_cds_pct),
                "med_tx_len_bp": float(med_len),
            }
        )

    by_bt = layer_tx2.groupby("biotype").apply(agg).reset_index(drop=True)
    by_bt = by_bt.sort_values(["tier_index", "layer_biotype"]).reset_index(drop=True)

    # Aggregation by tier
    tmp = by_bt.groupby(["tier", "tier_index"]).agg(
        n_layer_tx=("n_layer_tx", "sum"),
        n_built_tx=("n_built_tx", "sum"),
    ).reset_index()
    tmp["built_pct"] = tmp["n_built_tx"] / tmp["n_layer_tx"].clip(lower=1)
    by_tier = tmp.sort_values("tier_index").reset_index(drop=True)

    # Crosswalk: layer biotype -> core biotype counts (only built)
    cw = (
        layer_tx2.loc[layer_tx2["built"] == 1, ["biotype", "core_biotype"]]
        .value_counts()
        .reset_index(name="n_pairs")
        .rename(columns={"biotype": "layer_biotype"})
    )
    total_by_bt = layer_tx2.loc[layer_tx2["built"] == 1].groupby("biotype").size().to_dict()
    cw["pct_of_layer_biotype"] = cw.apply(lambda r: r["n_pairs"] / max(1, total_by_bt.get(r["layer_biotype"], 0)), axis=1)

    # Per-logic_name retention (useful when biotype mapping is sparse)
    def agg_ln(df: pd.DataFrame) -> pd.Series:
        n = len(df)
        nb = int(df["built"].sum())
        med_len = float(df["tx_len_bp"].median()) if n else 0.0
        return pd.Series({
            "logic_name": str(df.name),
            "n_layer_tx": int(n),
            "n_built_tx": int(nb),
            "built_pct": float((nb / n) if n else 0.0),
        })

    by_logic = layer_tx2.groupby("logic_name").apply(agg_ln).reset_index(drop=True)

    return by_bt, by_tier, cw, by_logic

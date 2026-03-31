from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


def _tier_of_biotype(biotype: str, b2t: Dict[str, dict | tuple]) -> Tuple[str, int]:
    v = b2t.get(biotype)
    if v is None:
        return ("unknown", 999)
    if isinstance(v, dict):
        return (str(v.get("tier", "unknown")), int(v.get("tier_index", 999)))
    tier, idx = v
    return (str(tier), int(idx))


def _quality_score(has_cds: int, length_bp: int) -> float:
    # Simple proxy: CDS dominates, then diminishing returns with length (cap at 2kb)
    return 0.6 * (1.0 if has_cds else 0.0) + 0.4 * min(max(length_bp, 0) / 2000.0, 1.0)


def analyze_competition(
    loci_df: pd.DataFrame,
    gene_map_df: pd.DataFrame,
    core_tx: pd.DataFrame,
    layer_tx: pd.DataFrame,
    layer_tr: pd.DataFrame,
    biotype_to_tier: Dict[str, dict | tuple],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Prep mapping
    b2t = biotype_to_tier or {}

    # Attach locus to layer/core transcripts via gene_id
    layer_tx2 = layer_tx.merge(
        gene_map_df[gene_map_df.db_kind == "layer"][['gene_id', 'locus_id_strict']],
        on="gene_id",
        how="left",
    )
    core_tx2 = core_tx.merge(
        gene_map_df[gene_map_df.db_kind == "core"][['gene_id', 'locus_id_strict']],
        on="gene_id",
        how="left",
    )
    # Basic fields
    layer_tx2["biotype"] = layer_tx2["biotype"].fillna("unknown").astype(str)
    layer_tx2["tx_len_bp"] = (layer_tx2["seq_region_end"].astype(int) - layer_tx2["seq_region_start"].astype(int) + 1)
    has_tr = set(layer_tr["transcript_id"].astype(int).tolist()) if not layer_tr.empty else set()
    layer_tx2["has_cds"] = layer_tx2["transcript_id"].astype(int).isin(has_tr).astype(int)
    layer_tx2["qscore"] = [
        _quality_score(int(h), int(l)) for h, l in zip(layer_tx2["has_cds"], layer_tx2["tx_len_bp"])
    ]
    # Strand harmonization for later
    # Use gene_map_df strand for locus; core/layer transcripts have their own strands, but loci are split by strand already

    # Build per-locus aggregates
    loci_df = loci_df.copy()
    if "locus_id" not in loci_df.columns:
        # reconstruct if necessary
        loci_df["locus_id"] = (
            loci_df.seq_region_name.astype(str)
            + ":" + loci_df.seq_region_strand.astype(str)
            + ":" + loci_df.locus_start.astype(int).astype(str)
            + ":" + loci_df.locus_end.astype(int).astype(str)
            + ":" + loci_df.groupby(["seq_region_name","seq_region_strand"]).cumcount().astype(str)
        )

    layer_by_loc = {k: g for k, g in layer_tx2.groupby("locus_id_strict")}
    core_by_loc = {k: g for k, g in core_tx2.groupby("locus_id_strict")}

    rows: List[dict] = []
    for _, loc in loci_df.iterrows():
        loc_id = loc.get("locus_id")
        seq = loc["seq_region_name"]
        strand = int(loc["seq_region_strand"])
        start = int(loc["locus_start"]) if "locus_start" in loc else None
        end = int(loc["locus_end"]) if "locus_end" in loc else None

        L = layer_by_loc.get(loc_id, pd.DataFrame(columns=layer_tx2.columns))
        C = core_by_loc.get(loc_id, pd.DataFrame(columns=core_tx2.columns))

        n_layer_tx = int(len(L))
        n_logic = int(L["logic_name"].nunique()) if n_layer_tx else 0
        cds_frac = float(L["has_cds"].mean()) if n_layer_tx else 0.0
        med_len = float(L["tx_len_bp"].median()) if n_layer_tx else 0.0

        # Per-logic/tier scores
        tier_scores: Dict[int, float] = {}
        tier_labels: Dict[int, str] = {}
        logic_scores = (
            L.groupby("logic_name")["qscore"].sum().sort_values(ascending=False) if n_layer_tx else pd.Series(dtype=float)
        )
        top_logic = str(logic_scores.index[0]) if len(logic_scores) else None
        top_logic_score = float(logic_scores.iloc[0]) if len(logic_scores) else 0.0

        for _, r in L.iterrows():
            tier_label, tier_idx = _tier_of_biotype(str(r["biotype"]), b2t)
            tier_scores[tier_idx] = tier_scores.get(tier_idx, 0.0) + float(r["qscore"])
            tier_labels[tier_idx] = tier_label
        # pick best tier present: minimal index with score >= 2.0, else minimal with >0
        present_tiers = sorted([(idx, sc) for idx, sc in tier_scores.items() if sc > 0.0], key=lambda x: x[0])
        best_tier_idx = None
        best_tier_label = None
        if present_tiers:
            cand = [t for t in present_tiers if t[1] >= 2.0]
            chosen = (cand[0] if cand else present_tiers[0])
            best_tier_idx = int(chosen[0])
            best_tier_label = tier_labels.get(best_tier_idx, None)

        # Core support tier (from overlapping layer tx)
        support_tier_idx = None
        support_tier_label = None
        n_core_tx = int(len(C))
        if n_core_tx and n_layer_tx:
            # For each layer tx, if overlaps any core tx (same locus), count towards support
            # (We assume locus already enforces strand; otherwise, check strands equal)
            support_scores: Dict[int, float] = {}
            for _, lr in L.iterrows():
                s, e = int(lr["seq_region_start"]), int(lr["seq_region_end"])
                ov = ((C["seq_region_start"].astype(int) <= e) & (C["seq_region_end"].astype(int) >= s))
                if ov.any():
                    tier_label, tier_idx = _tier_of_biotype(str(lr["biotype"]), b2t)
                    support_scores[tier_idx] = support_scores.get(tier_idx, 0.0) + float(lr["qscore"])
            if support_scores:
                chosen_idx = sorted([(idx, sc) for idx, sc in support_scores.items()], key=lambda x: x[0])[0][0]
                support_tier_idx = int(chosen_idx)
                support_tier_label = tier_labels.get(support_tier_idx, None)

        delta_tier = None
        if best_tier_idx is not None and support_tier_idx is not None:
            delta_tier = int(support_tier_idx) - int(best_tier_idx)

        # No-core reason
        no_core_category = None
        if n_core_tx == 0 and n_layer_tx > 0:
            # crude heuristics
            if best_tier_idx is not None and best_tier_idx <= 2 and cds_frac >= 0.5 and med_len >= 1200:
                no_core_category = "EXPECTED_TO_BUILD_BUT_MISSING"
            elif cds_frac < 0.2:
                no_core_category = "LACKS_CDS"
            elif med_len < 800:
                no_core_category = "SHORT_TRANSCRIPTS"
            elif best_tier_idx is None or best_tier_idx >= 10:
                no_core_category = "LOW_TIER_ONLY"
            else:
                no_core_category = "OTHER_NO_BUILD"

        rows.append(
            {
                "seq_region_name": seq,
                "seq_region_strand": strand,
                "locus_start": start,
                "locus_end": end,
                "locus_id": loc_id,
                "layer_transcript_count": n_layer_tx,
                "layer_logic_name_count": n_logic,
                "layer_cds_frac": cds_frac,
                "layer_median_tx_len_bp": med_len,
                "top_logic_name": top_logic,
                "top_logic_score": top_logic_score,
                "best_tier_label_present": best_tier_label,
                "best_tier_index_present": best_tier_idx,
                "core_transcript_count": n_core_tx,
                "support_best_tier_label": support_tier_label,
                "support_best_tier_index": support_tier_idx,
                "delta_tier_index": delta_tier,
                "no_core_category": no_core_category,
            }
        )

    per_locus = pd.DataFrame(rows)

    # Aggregated summary
    agg_rows = []
    # Distribution of no-core categories
    if not per_locus.empty:
        nc = per_locus[per_locus["no_core_category"].notna()]["no_core_category"].value_counts()
        for k, v in nc.items():
            agg_rows.append({"metric": "no_core_category", "key": k, "value": int(v)})
        # Delta tier stats among core-present loci
        dt = per_locus[per_locus["delta_tier_index"].notna()]["delta_tier_index"].astype(int)
        if not dt.empty:
            for q in (0.1, 0.5, 0.9):
                agg_rows.append({"metric": "delta_tier_quantile", "key": str(q), "value": float(dt.quantile(q))})
    agg = pd.DataFrame(agg_rows)

    return per_locus, agg


from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _merge_intervals(starts: Sequence[int], ends: Sequence[int]) -> List[Tuple[int, int]]:
    if len(starts) == 0:
        return []
    ivs = sorted(zip(map(int, starts), map(int, ends)))
    out: List[Tuple[int, int]] = []
    s, e = ivs[0]
    for s2, e2 in ivs[1:]:
        if s2 <= e + 0:  # overlap/touch
            e = max(e, e2)
        else:
            out.append((s, e))
            s, e = s2, e2
    out.append((s, e))
    return out


def _total_span(ivs: List[Tuple[int, int]]) -> int:
    return int(sum(e - s + 1 for s, e in ivs))


def _intersect_bp(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> int:
    i = j = 0
    cov = 0
    while i < len(a) and j < len(b):
        s1, e1 = a[i]
        s2, e2 = b[j]
        s = max(s1, s2)
        e = min(e1, e2)
        if s <= e:
            cov += (e - s + 1)
        if e1 < e2:
            i += 1
        else:
            j += 1
    return int(cov)


def _quantile_cap(x: pd.Series, q: float, default: float) -> float:
    x = x.dropna()
    if x.empty:
        return default
    return float(np.quantile(x.values, q))


def _safe_ratio(n: int, d: int) -> float:
    return float(n) / float(d if d != 0 else 1)


def load_mapping(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    try:
        if path.endswith(".yml") or path.endswith(".yaml"):
            try:
                import yaml  # type: ignore

                with open(path, "r", encoding="utf-8") as fh:
                    data = yaml.safe_load(fh) or {}
                return {str(k): str(v) for k, v in data.items()}
            except Exception:
                pass
        # TSV fallback: logic_name\tclass
        m: Dict[str, str] = {}
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = [p.strip() for p in line.split("\t")]
                if len(parts) >= 2:
                    m[parts[0]] = parts[1]
        return m
    except FileNotFoundError:
        return {}


def summarize_loci(
    loci_df: pd.DataFrame,
    gene_map_df: pd.DataFrame,
    core_genes: pd.DataFrame,
    layer_genes: pd.DataFrame,
    core_tx: pd.DataFrame,
    layer_tx: pd.DataFrame,
    evidence_map: Optional[Dict[str, str]] = None,
    locus_gap_bp: int = 5000,
) -> pd.DataFrame:
    evidence_map = evidence_map or {}
    # Pre-index for joins
    gm = gene_map_df[[
        "db_kind",
        "gene_id",
        "seq_region_name",
        "seq_region_strand",
        "locus_id_strict",
    ]].rename(columns={"locus_id_strict": "locus_id"})

    # Attach locus_id to genes
    core_genes2 = core_genes.merge(
        gm[gm.db_kind == "core"][["gene_id", "locus_id"]], on="gene_id", how="left"
    )
    layer_genes2 = layer_genes.merge(
        gm[gm.db_kind == "layer"][["gene_id", "locus_id"]], on="gene_id", how="left"
    )

    # Map transcripts to locus via gene_id
    core_tx2 = core_tx.merge(
        gm[gm.db_kind == "core"][["gene_id", "locus_id"]], on="gene_id", how="left"
    )
    layer_tx2 = layer_tx.merge(
        gm[gm.db_kind == "layer"][["gene_id", "locus_id"]], on="gene_id", how="left"
    )

    # Build aggregation per locus
    loci = loci_df.copy()
    if "locus_id" not in loci.columns:
        # Fallback: derive deterministic order per (seq,strand) but may not match mapping
        loci = loci.sort_values(["seq_region_name", "seq_region_strand", "locus_start", "locus_end"]).copy()
        loci["locus_id"] = (
            loci.seq_region_name.astype(str)
            + ":"
            + loci.seq_region_strand.astype(str)
            + ":"
            + loci.locus_start.astype(int).astype(str)
            + ":"
            + loci.locus_end.astype(int).astype(str)
            + ":"
            + loci.groupby(["seq_region_name", "seq_region_strand"]).cumcount().astype(str)
        )

    # Interval unions per locus
    def union_bp(df: pd.DataFrame) -> int:
        ivs = _merge_intervals(df.seq_region_start.values, df.seq_region_end.values)
        return _total_span(ivs)

    def cover_bp(layer_df: pd.DataFrame, core_df: pd.DataFrame) -> Tuple[int, int, int]:
        la = _merge_intervals(layer_df.seq_region_start.values, layer_df.seq_region_end.values)
        co = _merge_intervals(core_df.seq_region_start.values, core_df.seq_region_end.values)
        layer_span = _total_span(la)
        core_span = _total_span(co)
        inter = _intersect_bp(la, co)
        return layer_span, core_span, inter

    # Grouped stats per locus_id
    stats = []
    # Build fast lookup dicts
    core_by_locus = {k: v for k, v in core_genes2.groupby("locus_id")}
    layer_by_locus = {k: v for k, v in layer_genes2.groupby("locus_id")}
    core_tx_by_locus = {k: v for k, v in core_tx2.groupby("locus_id")}
    layer_tx_by_locus = {k: v for k, v in layer_tx2.groupby("locus_id")}

    for _, loc in loci.iterrows():
        cg = core_by_locus.get(loc.locus_id, pd.DataFrame(columns=core_genes2.columns))
        lg = layer_by_locus.get(loc.locus_id, pd.DataFrame(columns=layer_genes2.columns))
        ctx = core_tx_by_locus.get(loc.locus_id, pd.DataFrame(columns=core_tx2.columns))
        ltx = layer_tx_by_locus.get(loc.locus_id, pd.DataFrame(columns=layer_tx2.columns))

        # Counts
        core_gene_count = len(cg)
        layer_gene_count = len(lg)
        core_tx_count = len(ctx)
        layer_tx_count = len(ltx)

        # Span and coverage
        layer_span_bp, core_span_bp, covered_bp = cover_bp(lg, cg)
        uncovered_bp = max(layer_span_bp - covered_bp, 0)
        coverage_fraction = (_safe_ratio(covered_bp, layer_span_bp) if layer_span_bp > 0 else 0.0)

        # Diversity
        logic_names = set(ltx.logic_name.dropna().astype(str).tolist())
        logic_name_count = len(logic_names)
        classes = set((evidence_map.get(ln, "unknown") for ln in logic_names))
        class_count = len(classes)

        # Compression
        gene_ratio = _safe_ratio(layer_gene_count, core_gene_count)
        tx_ratio = _safe_ratio(layer_tx_count, core_tx_count)

        stats.append(
            {
                "seq_region_name": loc.seq_region_name,
                "seq_region_strand": int(loc.seq_region_strand),
                "locus_start": int(loc.locus_start),
                "locus_end": int(loc.locus_end),
                "locus_length": int(loc.locus_end) - int(loc.locus_start) + 1,
                "core_gene_count": int(core_gene_count),
                "core_transcript_count": int(core_tx_count),
                "layer_gene_count": int(layer_gene_count),
                "layer_transcript_count": int(layer_tx_count),
                "layer_logic_name_count": int(logic_name_count),
                "layer_evidence_class_count": int(class_count),
                "layer_span_bp": int(layer_span_bp),
                "core_span_bp": int(core_span_bp),
                "layer_bp_covered_by_core": int(covered_bp),
                "layer_bp_uncovered_by_core": int(uncovered_bp),
                "coverage_fraction": float(coverage_fraction),
                "layer_to_core_gene_ratio": float(gene_ratio),
                "layer_to_core_tx_ratio": float(tx_ratio),
            }
        )

    out = pd.DataFrame(stats)

    # Scoring for no-core loci
    no_core = out[(out.core_gene_count == 0) & (out.layer_gene_count > 0)].copy()
    # Robust caps (P95)
    cap_g = _quantile_cap(no_core.layer_gene_count, 0.95, 5.0) or 5.0
    cap_tx = _quantile_cap(no_core.layer_transcript_count, 0.95, 10.0) or 10.0
    cap_l = _quantile_cap(no_core.layer_logic_name_count, 0.95, 5.0) or 5.0
    cap_c = _quantile_cap(no_core.layer_evidence_class_count, 0.95, 4.0) or 4.0
    cap_bp = _quantile_cap(no_core.layer_span_bp, 0.95, 10000.0) or 10000.0

    def _norm(x, cap):
        return min(float(x), float(cap)) / float(cap if cap else 1.0)

    S = []
    for _, r in out.iterrows():
        if r.core_gene_count == 0 and r.layer_gene_count > 0:
            n_obj = _norm(r.layer_gene_count, cap_g)
            n_tx = _norm(r.layer_transcript_count, cap_tx)
            n_logic = _norm(r.layer_logic_name_count, cap_l)
            n_class = _norm(r.layer_evidence_class_count, cap_c)
            span = _norm(r.layer_span_bp, cap_bp)
            # strand consistency (we assume locus strand is the grouping strand; set to 1.0 as placeholder)
            strand = 1.0
            # coherence needs subcluster count; placeholder 1.0 for Phase 2
            coherence = 1.0
            score = 0.25 * n_obj + 0.2 * n_tx + 0.15 * n_logic + 0.15 * n_class + 0.15 * span + 0.05 * strand + 0.05 * coherence
            S.append(score)
        else:
            S.append(0.0)
    out["no_core_score"] = S
    out["evidence_rich_no_core_flag"] = ((out.core_gene_count == 0) & (out.layer_gene_count > 0) & (out.no_core_score >= 0.6)).astype(int)

    # Diagnostic categories (coarse)
    cat = []
    for _, r in out.iterrows():
        if r.core_gene_count == 0 and r.layer_gene_count > 0:
            cat.append("LAYER_ONLY_NO_CORE")
        elif r.core_gene_count > 0 and r.layer_gene_count == 0:
            cat.append("CORE_ONLY")
        elif r.coverage_fraction >= 0.8 and r.layer_to_core_tx_ratio <= 1.5:
            cat.append("CORE_AND_LAYER_CONCORDANT")
        elif (r.coverage_fraction < 0.8) or (r.layer_to_core_tx_ratio > 2.0):
            cat.append("CORE_PRESENT_LAYER_EXCESS")
        else:
            cat.append("OTHER")
    out["diagnostic_category"] = cat

    # Simple underbuilt flag
    out["underbuilt_locus_flag"] = (
        (out.core_gene_count > 0)
        & ((out.coverage_fraction < 0.6) | (out.layer_to_core_tx_ratio >= 2.0))
        & (out.layer_span_bp - out.layer_bp_covered_by_core >= 2000)
    ).astype(int)

    return out

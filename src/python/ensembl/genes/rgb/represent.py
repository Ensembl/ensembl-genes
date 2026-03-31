from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


@dataclass
class MatchResult:
    class_label: str
    core_transcript_id: Optional[int]
    core_gene_id: Optional[int]
    matched_introns: int
    total_layer_introns: int
    max_splice_delta_bp: int
    exon_overlap_frac: float


def _build_exon_index(exons: pd.DataFrame) -> Dict[Tuple[str, int], Dict[int, List[Tuple[int, int]]]]:
    idx: Dict[Tuple[str, int], Dict[int, List[Tuple[int, int]]]] = {}
    if exons.empty:
        return idx
    for (seq, strand), g in exons.groupby(["seq_region_name", "seq_region_strand"]):
        d: Dict[int, List[Tuple[int, int]]] = {}
        for tid, gt in g.sort_values(["transcript_id", "exon_rank"]).groupby("transcript_id"):
            ex = list(zip(gt["exon_start"].astype(int).to_list(), gt["exon_end"].astype(int).to_list()))
            d[int(tid)] = ex
        idx[(str(seq), int(strand))] = d
    return idx


def _bounds_from_exons(ex_list: List[Tuple[int, int]]) -> Tuple[int, int]:
    if not ex_list:
        return (0, -1)
    s = min(s for s, _ in ex_list)
    e = max(e for _, e in ex_list)
    return (s, e)


def _introns_from_exons(ex_list: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if len(ex_list) < 2:
        return []
    ex = sorted(ex_list)
    return [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]


def _match_introns(layer_introns: List[Tuple[int, int]], core_introns: List[Tuple[int, int]], tol_bp: int) -> Tuple[int, int]:
    """Return (matched_count, max_delta). Delta is max(|d_end|,|d_start|) across matched pairs."""
    if not layer_introns:
        return (0, 0)
    matched = 0
    max_delta = 0
    if not core_introns:
        return (0, 0)
    # For speed, use a set of rounded intron coordinates as an approximation, but also compute exact deltas
    for (lend, lstart) in layer_introns:
        best = None
        for (cend, cstart) in core_introns:
            if abs(lend - cend) <= tol_bp and abs(lstart - cstart) <= tol_bp:
                d = max(abs(lend - cend), abs(lstart - cstart))
                if best is None or d < best:
                    best = d
        if best is not None:
            matched += 1
            if best > max_delta:
                max_delta = best
    return (matched, max_delta)


def _exon_overlap_frac_single(layer_exon: Tuple[int, int], core_exons: List[Tuple[int, int]]) -> float:
    if not core_exons:
        return 0.0
    s, e = layer_exon
    L = max(0, e - s + 1)
    if L == 0:
        return 0.0
    cov = 0
    for (xs, xe) in core_exons:
        a = max(s, xs)
        b = min(e, xe)
        if a <= b:
            cov += (b - a + 1)
    return cov / L


def classify_layer_tx(
    seq: str,
    strand: int,
    layer_tid: int,
    layer_gene_id: int,
    layer_exons: List[Tuple[int, int]],
    core_exon_index: Dict[int, List[Tuple[int, int]]],
    core_bounds_arr: np.ndarray,
    core_tid_arr: np.ndarray,
    core_gene_by_tid: Dict[int, int],
    *,
    tol_bp: int = 5,
    intron_subset_frac: float = 0.8,
    single_exon_overlap_frac: float = 0.8,
) -> MatchResult:
    # Determine layer introns and span
    layer_exons = sorted(layer_exons)
    lstart, lend = _bounds_from_exons(layer_exons)
    lintr = _introns_from_exons(layer_exons)
    # Preselect candidate core transcripts by span overlap
    # core_bounds_arr is shape (N,2)
    if core_bounds_arr.size == 0:
        return MatchResult("UNREPRESENTED", None, None, 0, len(lintr), 0, 0.0)
    ov = (core_bounds_arr[:, 0] <= lend) & (core_bounds_arr[:, 1] >= lstart)
    cand_idx = np.where(ov)[0]
    if cand_idx.size == 0:
        return MatchResult("UNREPRESENTED", None, None, 0, len(lintr), 0, 0.0)

    best = None  # (score_tuple, tid, gid, matched, total, max_delta, exon_frac)

    if len(lintr) >= 1:
        total = len(lintr)
        for i in cand_idx:
            tid = int(core_tid_arr[i])
            cex = core_exon_index.get(tid, [])
            cintr = _introns_from_exons(cex)
            matched, maxd = _match_introns(lintr, cintr, tol_bp)
            frac = matched / total if total else 0.0
            # score tuple prioritizes exact subset, then near subset, then matched count, then smaller delta
            exact = 1 if matched == total else 0
            near = 1 if frac >= intron_subset_frac else 0
            score = (exact, near, matched, -maxd)
            if best is None or score > best[0]:
                best = (score, tid, int(core_gene_by_tid.get(tid, -1)), matched, total, maxd, 0.0)
        if best is None:
            return MatchResult("UNREPRESENTED", None, None, 0, total, 0, 0.0)
        exact, near, matched, maxd = best[0][0], best[0][1], best[3], best[5]
        if exact:
            return MatchResult("EXACT_SUBSET", best[1], best[2], matched, total, maxd, 0.0)
        if near:
            return MatchResult("NEAR_SUBSET", best[1], best[2], matched, total, maxd, 0.0)
        # fallthrough to overlap only if any span overlap existed
        return MatchResult("OVERLAP_ONLY", best[1], best[2], matched, total, maxd, 0.0)

    # single-exon fallback
    if len(layer_exons) == 1:
        se = layer_exons[0]
        best_frac = 0.0
        best_tid = None
        for i in cand_idx:
            tid = int(core_tid_arr[i])
            cex = core_exon_index.get(tid, [])
            frac = _exon_overlap_frac_single(se, cex)
            if frac > best_frac:
                best_frac = frac
                best_tid = tid
        if best_tid is not None and best_frac >= single_exon_overlap_frac:
            return MatchResult("SINGLE_EXON_COVERED", best_tid, int(core_gene_by_tid.get(best_tid, -1)), 0, 0, 0, best_frac)
        if best_tid is not None and best_frac > 0:
            return MatchResult("OVERLAP_ONLY", best_tid, int(core_gene_by_tid.get(best_tid, -1)), 0, 0, 0, best_frac)
        return MatchResult("UNREPRESENTED", None, None, 0, 0, 0, 0.0)

    # No exons info for this transcript
    return MatchResult("UNREPRESENTED", None, None, 0, 0, 0, 0.0)


def represent_layer_vs_core(
    layer_exons: pd.DataFrame,
    core_exons: pd.DataFrame,
    layer_tx: pd.DataFrame,
    core_tx: pd.DataFrame,
    *,
    tol_bp: int = 5,
    intron_subset_frac: float = 0.8,
    single_exon_overlap_frac: float = 0.8,
) -> pd.DataFrame:
    # Build exon indexes
    lex_idx = _build_exon_index(layer_exons)
    cex_idx = _build_exon_index(core_exons)

    # Index core by (seq,strand) with bounds and ids
    results: List[dict] = []
    # maps for quick gene lookup
    core_gene_by_tid = {int(r.transcript_id): int(r.gene_id) for _, r in core_tx.iterrows()} if not core_tx.empty else {}
    # layer minimal fields
    if not layer_tx.empty:
        layer_tx = layer_tx[[
            "transcript_id",
            "gene_id",
            "seq_region_name",
            "seq_region_strand",
            "logic_name",
            "biotype",
        ] + ([c for c in ("stable_id",) if c in layer_tx.columns])].copy()

    for key, group in layer_tx.groupby(["seq_region_name", "seq_region_strand"], dropna=False):
        seq = str(key[0]); strand = int(key[1])
        # Build core arrays for this key
        cidx = cex_idx.get((seq, strand), {})
        if not cidx:
            core_bounds = np.empty((0, 2), dtype=np.int64)
            core_tids = np.empty((0,), dtype=np.int64)
        else:
            tids = np.array(list(cidx.keys()), dtype=np.int64)
            bounds = np.array([_bounds_from_exons(cidx[int(t)]) for t in tids], dtype=np.int64)
            core_bounds = bounds
            core_tids = tids

        for _, r in group.iterrows():
            tid = int(r.transcript_id)
            gid = int(r.gene_id)
            lex = lex_idx.get((seq, strand), {}).get(tid, [])
            m = classify_layer_tx(
                seq, strand, tid, gid, lex, cidx,
                core_bounds, core_tids, core_gene_by_tid,
                tol_bp=tol_bp, intron_subset_frac=intron_subset_frac,
                single_exon_overlap_frac=single_exon_overlap_frac,
            )
            results.append({
                "layer_transcript_id": tid,
                "layer_gene_id": gid,
                "seq_region_name": seq,
                "seq_region_strand": strand,
                "logic_name": r.get("logic_name", None),
                "biotype": r.get("biotype", None),
                "representation_class": m.class_label,
                "matched_core_transcript_id": m.core_transcript_id,
                "matched_core_gene_id": m.core_gene_id,
                "matched_introns": m.matched_introns,
                "total_layer_introns": m.total_layer_introns,
                "max_splice_delta_bp": m.max_splice_delta_bp,
                "single_exon_overlap_frac": m.exon_overlap_frac,
            })

    return pd.DataFrame(results)


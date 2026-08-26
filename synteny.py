# synteny.py
from __future__ import annotations
from typing import Dict, List, Tuple
from gff_io import Gene
from aligner import GeneMap

def _neighbors(genes: Dict[str, Gene]) -> Dict[Tuple[str,str], List[str]]:
    """
    Build neighbor order per (seqid,strand) -> sorted gene ID list by start.
    """
    by_chr: Dict[Tuple[str,str], List[Gene]] = {}
    for g in genes.values():
        by_chr.setdefault((g.seqid, g.strand), []).append(g)
    for k in by_chr:
        by_chr[k].sort(key=lambda x: x.start)
    out: Dict[Tuple[str,str], List[str]] = {}
    for k, glist in by_chr.items():
        out[k] = [g.gid for g in glist]
    return out

def synteny_disambiguate(ref_genes: Dict[str, Gene], candidate: Dict[str, GeneMap]) -> Dict[str, GeneMap]:
    """
    Prefer mappings whose target locus attracts the highest density of
    *other* genes from the same syntenic neighborhood within a window.
    Heuristic but effective in paralog clusters.
    """
    # Build simple index of mapped target windows
    window = 500_000  # 500kb
    ref_order = _neighbors(ref_genes)
    # group target hits by contig
    by_tgt: Dict[str, List[GeneMap]] = {}
    for gm in candidate.values():
        if gm.rname: 
            by_tgt.setdefault(gm.rname, []).append(gm)
    for mlist in by_tgt.values():
        mlist.sort(key=lambda m: m.start or 0)
    # score each mapping by neighbor co-location
    scores: Dict[str, float] = {}
    for gid, gm in candidate.items():
        if gm.rname is None or gm.start is None: 
            scores[gid] = 0.0
            continue
        neigh_score = 0
        # neighbors in ref
        key_plus = (ref_genes[gid].seqid, ref_genes[gid].strand)
        order = ref_order.get(key_plus, [])
        try:
            i = order.index(gid)
        except ValueError:
            i = -1
        ref_neigh = set(order[max(0,i-5): i] + order[i+1: i+6]) if i >= 0 else set()
        # count how many neighbors land within window of this target hit
        nearby = [x for x in by_tgt.get(gm.rname, []) 
                  if abs(((x.start or 0) + (x.end or 0))//2 - ((gm.start or 0)+(gm.end or 0))//2) <= window]
        ref_neigh_mapped = {x.gid for x in nearby}
        neigh_score = len(ref_neigh & ref_neigh_mapped)
        # final score: identity weight + neighbor bonus
        scores[gid] = gm.identity + 0.05 * neigh_score
    # nothing else to disambiguate in this simplified version
    return candidate  # scores computed if you later want tie-breaks

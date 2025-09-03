# cigar_map.py
from __future__ import annotations
from typing import List, Tuple

def parse_cigar(cg: str) -> List[Tuple[int, str]]:
    num = ""
    out: List[Tuple[int, str]] = []
    for ch in cg:
        if ch.isdigit():
            num += ch
        else:
            if not num: raise ValueError(f"Bad CIGAR: {cg}")
            out.append((int(num), ch))
            num = ""
    if num: raise ValueError(f"Trailing length in CIGAR: {cg}")
    return out

def project_ref_interval(cigar: str, aln_pos_1: int, q_start_1: int, q_end_1: int) -> Tuple[int,int]:
    """
    Map a query (reference gene) interval [q_start_1, q_end_1] (1-based, inclusive)
    through an alignment that starts on target at aln_pos_1 with CIGAR 'cigar'.
    Returns target 1-based inclusive coordinates. Raises if region not fully covered.
    """
    ops = parse_cigar(cigar)
    q = 1
    t = aln_pos_1
    t_start = t_end = None
    q_s, q_e = q_start_1, q_end_1
    for ln, op in ops:
        if op in "MX=":      # consumes both
            for _ in range(ln):
                if q == q_s: t_start = t
                if q == q_e: 
                    t_end = t
                    break
                q += 1; t += 1
            if t_end is not None: break
        elif op == "I":      # query insertion: consume query only
            q += ln
        elif op == "D" or op == "N":  # deletion/skip: consume target only
            t += ln
        elif op == "S" or op == "H":  # clipping: consume query
            q += ln
        elif op == "P":               # padding: ignore
            continue
        else:
            raise ValueError(f"Unhandled CIGAR op {op}")
    if t_start is None or t_end is None:
        raise RuntimeError("Interval not fully mapped by alignment")
    if t_start > t_end: t_start, t_end = t_end, t_start
    return t_start, t_end

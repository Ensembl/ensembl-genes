#!/usr/bin/env python3
"""
Lifton ↔ Reference ID Mapper
=============================

Maps transcripts/genes from a LiftOn GFF3 onto a reference (e.g., Ensembl) GFF3
on the *same target assembly*, with emphasis on exon/intron structure and local
context. Produces:

1) transcript_pairs.tsv — one-to-one transcript mapping with component scores
2) gene_pairs.tsv        — gene mapping aggregated from transcript pairs
3) (optional) mapped_reference.gff3 — reference GFF3 annotated with LiftOn IDs

Speed/accuracy balance:
- Candidate search by locus (same contig/strand; overlap ±window)
- Fast features: intron-chain similarity, internal-exon Jaccard, exon Jaccard,
  exon-count similarity, boundary similarity; optionally uses LiftOn-supplied
  protein_identity/dna_identity attributes as a small prior when present.
- Greedy maximum-weight matching to enforce one-to-one transcript mapping.

No external dependencies (stdlib only).

Example
-------
python lifton_id_mapper.py \
  --lifton liftOn.gff3 \
  --reference ensembl.gff3 \
  --out-prefix out/mapping \
  --window 100000 \
  --min-score 0.75 \
  --rename-mode alias  # or: rename

Outputs
-------
- out/mapping.transcript_pairs.tsv
- out/mapping.gene_pairs.tsv
- out/mapping.mapped_reference.gff3  (if --rewrite-reference provided)

Notes
-----
- Designed for GFF3 with gene/mRNA(or transcript)/exon (CDS optional).
- Biotype is ignored for scoring (captured as metadata only).
- UTR differences are expected; internal exons and intron chains are weighted
  more than boundaries/UTRs.
"""
from __future__ import annotations
import argparse
import collections
import math
import os
import sys
import gzip
from bisect import bisect_left
from typing import Dict, List, Tuple, Iterable, Optional, Set

###############################################################################
# Utility: GFF3 parsing/writing
###############################################################################

def _open_maybe_gzip(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')

def parse_gff3_attributes(attr_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    if not attr_str or attr_str == '.':
        return d
    parts = attr_str.strip().split(';')
    for p in parts:
        if not p:
            continue
        if '=' in p:
            k, v = p.split('=', 1)
            d[k] = v
        else:
            # tolerate bare keys
            d[p] = ''
    return d

def format_gff3_attributes(attrs: Dict[str, str]) -> str:
    # keep order stable for readability
    return ';'.join(f"{k}={v}" for k, v in attrs.items()) if attrs else '.'

###############################################################################
# Data models
###############################################################################

class Transcript:
    __slots__ = (
        'id', 'gene_id', 'contig', 'strand', 'exons', 'attrs',
        '_sorted', '_merged_exons', '_merged_internal_exons', '_intron_pairs',
        '_tss', '_tes'
    )
    def __init__(self, tid: str, gene_id: str, contig: str, strand: str):
        self.id = tid
        self.gene_id = gene_id
        self.contig = contig
        self.strand = strand
        self.exons: List[Tuple[int,int]] = []  # [ (start,end) ] inclusive 1-based
        self.attrs: Dict[str, str] = {}
        self._sorted = False
        self._merged_exons: Optional[List[Tuple[int,int]]] = None
        self._merged_internal_exons: Optional[List[Tuple[int,int]]] = None
        self._intron_pairs: Optional[List[Tuple[int,int]]] = None
        self._tss: Optional[int] = None
        self._tes: Optional[int] = None

    def add_exon(self, start: int, end: int):
        self.exons.append((min(start, end), max(start, end)))
        self._sorted = False
        self._merged_exons = None
        self._merged_internal_exons = None
        self._intron_pairs = None
        self._tss = None
        self._tes = None

    def _ensure_sorted(self):
        if not self._sorted:
            self.exons.sort()
            self._sorted = True

    def exon_count(self) -> int:
        return len(self.exons)

    def merged_exons(self) -> List[Tuple[int,int]]:
        if self._merged_exons is not None:
            return self._merged_exons
        self._ensure_sorted()
        out: List[Tuple[int,int]] = []
        for s,e in self.exons:
            if not out or s > out[-1][1] + 1:
                out.append([s,e])
            else:
                out[-1][1] = max(out[-1][1], e)
        self._merged_exons = [(int(s), int(e)) for s,e in out]
        return self._merged_exons

    def merged_internal_exons(self) -> List[Tuple[int,int]]:
        if self._merged_internal_exons is not None:
            return self._merged_internal_exons
        self._ensure_sorted()
        n = len(self.exons)
        if n <= 2:
            # no internal exons; fall back to merged full exons (weighting will handle it)
            self._merged_internal_exons = self.merged_exons()
            return self._merged_internal_exons
        # internal = drop first and last in transcriptional order
        # transcriptional order = genomic ascending on '+'; descending on '-'
        if self.strand == '+':
            internals = self.exons[1:-1]
        else:
            # exons are stored ascending; for '-' strand, internal exons are also 1..n
            internals = self.exons[1:-1]
        # merge internal intervals
        out: List[Tuple[int,int]] = []
        for s,e in internals:
            if not out or s > out[-1][1] + 1:
                out.append([s,e])
            else:
                out[-1][1] = max(out[-1][1], e)
        self._merged_internal_exons = [(int(s), int(e)) for s,e in out] if out else []
        return self._merged_internal_exons

    def intron_pairs(self) -> List[Tuple[int,int]]:
        if self._intron_pairs is not None:
            return self._intron_pairs
        self._ensure_sorted()
        pairs: List[Tuple[int,int]] = []
        if len(self.exons) >= 2:
            for i in range(len(self.exons) - 1):
                left = self.exons[i][1]
                right = self.exons[i+1][0]
                pairs.append((left, right))  # splice sites in genomic coords
        self._intron_pairs = pairs
        return pairs

    def span(self) -> Tuple[int,int]:
        self._ensure_sorted()
        if not self.exons:
            return (0,0)
        return (self.exons[0][0], self.exons[-1][1])

    def tss(self) -> int:
        if self._tss is not None:
            return self._tss
        s, e = self.span()
        self._tss = s if self.strand == '+' else e
        return self._tss

    def tes(self) -> int:
        if self._tes is not None:
            return self._tes
        s, e = self.span()
        self._tes = e if self.strand == '+' else s
        return self._tes

class Gene:
    __slots__ = ('id', 'contig', 'strand', 'transcripts', 'span_start', 'span_end', 'attrs')
    def __init__(self, gid: str, contig: str, strand: str):
        self.id = gid
        self.contig = contig
        self.strand = strand
        self.transcripts: Dict[str, Transcript] = {}
        self.span_start = 10**18
        self.span_end = -1
        self.attrs: Dict[str, str] = {}

    def add_tx(self, tx: Transcript):
        self.transcripts[tx.id] = tx
        s,e = tx.span()
        if s and e:
            self.span_start = min(self.span_start, s)
            self.span_end = max(self.span_end, e)

    def span(self) -> Tuple[int,int]:
        if self.span_end < self.span_start:
            return (0,0)
        return (self.span_start, self.span_end)

class Annotation:
    def __init__(self, name: str):
        self.name = name
        self.genes: Dict[str, Gene] = {}
        self.tx_index: Dict[str, Transcript] = {}

    def ensure_gene(self, gid: str, contig: str, strand: str) -> Gene:
        g = self.genes.get(gid)
        if g is None:
            g = Gene(gid, contig, strand)
            self.genes[gid] = g
        return g

###############################################################################
# Parsing GFF3 into Annotation
###############################################################################

def load_gff3_as_annotation(path: str, label: str) -> Annotation:
    ann = Annotation(label)
    # Temporary storage: transcript attributes before exons
    tx_meta: Dict[str, Tuple[str,str,Dict[str,str]]] = {}

    with _open_maybe_gzip(path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_s = parts
            start_i = int(start)
            end_i = int(end)
            attrs = parse_gff3_attributes(attrs_s)
            ftype_l = ftype.lower()

            if ftype_l in ('mrna', 'transcript'):
                tid = attrs.get('ID') or attrs.get('transcript_id') or attrs.get('Name')
                parent = attrs.get('Parent') or attrs.get('gene_id')
                if not tid or not parent:
                    continue
                # Some IDs come as transcript:ENST..., gene:ENSG...
                tid = tid
                gid = parent
                tx_meta[tid] = (seqid, strand, attrs)
                # Create placeholders (we'll add exons later)
                if gid not in ann.genes:
                    ann.genes[gid] = Gene(gid, seqid, strand)
                if tid not in ann.tx_index:
                    t = Transcript(tid, gid, seqid, strand)
                    t.attrs.update(attrs)
                    ann.tx_index[tid] = t

            elif ftype_l == 'exon':
                parent = attrs.get('Parent')
                if not parent:
                    continue
                # Parent may be comma-separated; we assign to each
                for tid in parent.split(','):
                    t = ann.tx_index.get(tid)
                    if t is None:
                        # exon precedes transcript line: create minimal transcript
                        seqid2 = seqid
                        strand2 = strand
                        t = Transcript(tid, attrs.get('gene_id',''), seqid2, strand2)
                        ann.tx_index[tid] = t
                    t.add_exon(start_i, end_i)

            elif ftype_l == 'gene':
                gid = attrs.get('ID') or attrs.get('gene_id') or attrs.get('Name')
                if not gid:
                    continue
                g = ann.genes.get(gid)
                if g is None:
                    g = Gene(gid, seqid, strand)
                    g.attrs.update(attrs)
                    ann.genes[gid] = g

    # Link transcripts to genes & compute spans
    for tid, t in list(ann.tx_index.items()):
        # fill basic metadata from tx_meta if present
        if tid in tx_meta:
            seqid, strand, a = tx_meta[tid]
            t.contig = seqid
            t.strand = strand
            # keep attrs union (prefer earlier keys)
            for k,v in a.items():
                ann.tx_index[tid].attrs.setdefault(k, v)
        gid = t.gene_id or t.attrs.get('Parent') or t.attrs.get('gene_id')
        if not gid:
            # create synthetic gene ID based on transcript if needed
            gid = f"gene_of:{tid}"
            t.gene_id = gid
        g = ann.genes.get(gid)
        if g is None:
            g = Gene(gid, t.contig, t.strand)
            ann.genes[gid] = g
        t.gene_id = gid
        g.add_tx(t)

    return ann

###############################################################################
# Interval index (per contig/strand for genes)
###############################################################################

class GeneIntervalIndex:
    def __init__(self, ann: Annotation):
        # contig -> strand -> list of (start, end, gene_id)
        self.idx: Dict[str, Dict[str, List[Tuple[int,int,str]]]] = collections.defaultdict(lambda: collections.defaultdict(list))
        for gid, g in ann.genes.items():
            s,e = g.span()
            if s == 0 and e == 0:
                # derive from transcripts if needed
                ss = 10**18; ee = -1
                for tx in g.transcripts.values():
                    ts,te = tx.span()
                    if ts:
                        ss = min(ss, ts); ee = max(ee, te)
                if ss <= ee:
                    s,e = ss,ee
            if s == 0 and e == 0:
                continue
            self.idx[g.contig][g.strand].append((s,e,gid))
        # sort by start for binary search
        for contig in self.idx:
            for strand in self.idx[contig]:
                self.idx[contig][strand].sort(key=lambda x: x[0])

    def query(self, contig: str, strand: str, qs: int, qe: int, window: int = 100000) -> List[str]:
        out: List[str] = []
        lst = self.idx.get(contig, {}).get(strand, [])
        if not lst:
            return out
        starts = [x[0] for x in lst]
        left = max(0, bisect_left(starts, qs - window) - 5)  # small cushion
        max_start = qe + window
        for i in range(left, len(lst)):
            s,e,gid = lst[i]
            if s > max_start:
                break
            if e >= qs - window:
                out.append(gid)
        return out

###############################################################################
# Similarity metrics
###############################################################################

def _interval_len(intervals: List[Tuple[int,int]]) -> int:
    return sum(e - s + 1 for s,e in intervals)

def _interval_intersection_len(a: List[Tuple[int,int]], b: List[Tuple[int,int]]) -> int:
    i = j = 0
    total = 0
    while i < len(a) and j < len(b):
        s1,e1 = a[i]
        s2,e2 = b[j]
        if e1 < s2:
            i += 1
            continue
        if e2 < s1:
            j += 1
            continue
        overlap_s = max(s1, s2)
        overlap_e = min(e1, e2)
        if overlap_s <= overlap_e:
            total += (overlap_e - overlap_s + 1)
        if e1 <= e2:
            i += 1
        else:
            j += 1
    return total

def jaccard_len(a: List[Tuple[int,int]], b: List[Tuple[int,int]]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = _interval_intersection_len(a,b)
    la = _interval_len(a)
    lb = _interval_len(b)
    denom = la + lb - inter
    return (inter / denom) if denom > 0 else 0.0

def intron_chain_similarity(q: Transcript, r: Transcript) -> float:
    iq = q.intron_pairs()
    ir = r.intron_pairs()
    if not iq and not ir:
        # single-exon both: use exon Jaccard as proxy (closer to structure)
        return jaccard_len(q.merged_exons(), r.merged_exons())
    if not iq or not ir:
        return 0.0
    set_q = set(iq)
    set_r = set(ir)
    inter = len(set_q & set_r)
    union = len(set_q | set_r)
    return inter / union if union else 0.0

def exon_count_similarity(q: Transcript, r: Transcript) -> float:
    nq = q.exon_count(); nr = r.exon_count()
    if nq == 0 and nr == 0:
        return 1.0
    if nq == 0 or nr == 0:
        return 0.0
    return 1.0 - (abs(nq - nr) / max(nq, nr))

def boundary_similarity(q: Transcript, r: Transcript, scale: int = 1000) -> float:
    # Penalize large shifts in TSS/TES; small shifts tolerated (UTR changes)
    dtss = abs(q.tss() - r.tss())
    dtes = abs(q.tes() - r.tes())
    # Map distance to 0..1 with a soft decay. e.g., 0 at >= 5*scale
    def sim(d):
        x = d / float(scale)
        return 1.0 / (1.0 + x)
    return 0.5 * (sim(dtss) + sim(dtes))

###############################################################################
# Scoring function
###############################################################################

def lifton_identity_prior(t: Transcript) -> float:
    # Use LiftOn-provided identities (if present) as a light boost.
    # Expect attributes like protein_identity=0.998; dna_identity=0.999
    vals: List[float] = []
    for k in ('protein_identity', 'dna_identity'):
        v = t.attrs.get(k)
        if v is None:
            continue
        try:
            vals.append(float(v))
        except Exception:
            pass
    if not vals:
        return 0.0
    # map [0,1] to [0,1] directly; we'll weight it lightly in score
    return max(0.0, min(1.0, sum(vals) / len(vals)))

def score_pair(q: Transcript, r: Transcript) -> Tuple[float, Dict[str, float]]:
    # Components
    intron_sim = intron_chain_similarity(q, r)
    jacc_all = jaccard_len(q.merged_exons(), r.merged_exons())
    jacc_internal = jaccard_len(q.merged_internal_exons(), r.merged_internal_exons())
    ex_ct_sim = exon_count_similarity(q, r)
    bnd_sim = boundary_similarity(q, r)
    prior = lifton_identity_prior(q)  # Only from LiftOn side

    # Weighted sum — internal structure is emphasized
    score = (
        0.50 * intron_sim +
        0.25 * jacc_internal +
        0.10 * jacc_all +
        0.05 * ex_ct_sim +
        0.05 * bnd_sim +
        0.05 * prior
    )

    details = {
        'intron_sim': intron_sim,
        'jacc_internal': jacc_internal,
        'jacc_all': jacc_all,
        'exon_count_sim': ex_ct_sim,
        'boundary_sim': bnd_sim,
        'lifton_identity_prior': prior,
    }
    return score, details

###############################################################################
# Matching logic
###############################################################################

Pair = collections.namedtuple('Pair', 'q_id r_id score details')


def candidate_transcripts_for(q: Transcript, ref_ann: Annotation, idx: GeneIntervalIndex, window: int) -> List[Transcript]:
    qs, qe = q.span()
    gids = idx.query(q.contig, q.strand, qs, qe, window)
    cands: List[Transcript] = []
    for gid in gids:
        g = ref_ann.genes.get(gid)
        if not g:
            continue
        cands.extend(g.transcripts.values())
    return cands


def compute_pairs(lifton: Annotation, reference: Annotation, window: int, min_candidate_score: float, topk: int) -> List[Pair]:
    idx = GeneIntervalIndex(reference)
    pairs: List[Pair] = []

    # Pre-score all candidates per LiftOn transcript
    per_q_candidates: Dict[str, List[Pair]] = {}
    for qid, q in lifton.tx_index.items():
        if not q.exons:
            continue
        cands = candidate_transcripts_for(q, reference, idx, window)
        best: List[Tuple[float, Transcript, Dict[str,float]]] = []
        for r in cands:
            if r.contig != q.contig or r.strand != q.strand:
                continue
            sc, det = score_pair(q, r)
            if sc >= min_candidate_score:
                best.append((sc, r, det))
        if not best:
            continue
        best.sort(key=lambda x: x[0], reverse=True)
        best = best[:topk]
        per_q_candidates[qid] = [Pair(qid, r.id, sc, det) for sc, r, det in best]

    # Flatten and greedy assign for one-to-one
    all_pairs = [p for lst in per_q_candidates.values() for p in lst]
    all_pairs.sort(key=lambda p: p.score, reverse=True)

    assigned_q: Set[str] = set()
    assigned_r: Set[str] = set()
    chosen: List[Pair] = []

    for p in all_pairs:
        if p.q_id in assigned_q or p.r_id in assigned_r:
            continue
        assigned_q.add(p.q_id)
        assigned_r.add(p.r_id)
        chosen.append(p)

    # Second pass: try to assign remaining q to their next available candidate
    remaining_q = [qid for qid in per_q_candidates.keys() if qid not in assigned_q]
    for qid in remaining_q:
        for p in per_q_candidates[qid]:
            if p.r_id not in assigned_r:
                assigned_q.add(qid)
                assigned_r.add(p.r_id)
                chosen.append(p)
                break

    chosen.sort(key=lambda p: p.score, reverse=True)
    return chosen

###############################################################################
# Gene aggregation
###############################################################################

def aggregate_gene_pairs(chosen_pairs: List[Pair], lifton: Annotation, reference: Annotation, min_gene_fraction: float = 0.6) -> List[Tuple[str,str,float,float,int]]:
    # For each LiftOn gene, tally weights to reference genes
    tallies: Dict[str, Dict[str, float]] = collections.defaultdict(lambda: collections.defaultdict(float))
    counts: Dict[str, Dict[str, int]] = collections.defaultdict(lambda: collections.defaultdict(int))
    for p in chosen_pairs:
        q_tx = lifton.tx_index[p.q_id]
        r_tx = reference.tx_index[p.r_id]
        tallies[q_tx.gene_id][r_tx.gene_id] += p.score
        counts[q_tx.gene_id][r_tx.gene_id] += 1

    out: List[Tuple[str,str,float,float,int]] = []
    for q_gid, bucket in tallies.items():
        total = sum(bucket.values())
        if total <= 0:
            continue
        r_gid, w = max(bucket.items(), key=lambda kv: kv[1])
        frac = w / total
        n = counts[q_gid][r_gid]
        out.append((q_gid, r_gid, w, frac, n))
    # filter by fraction
    out = [row for row in out if row[3] >= min_gene_fraction]
    # sort by weight descending
    out.sort(key=lambda x: x[2], reverse=True)
    return out

###############################################################################
# Writing outputs
###############################################################################

def write_transcript_pairs(pairs: List[Pair], lifton: Annotation, reference: Annotation, path: str, conf_thresh: float, good_thresh: float):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as out:
        header = [
            'lifton_tx','ref_tx','score','status','intron_sim','jacc_internal','jacc_all',
            'exon_count_sim','boundary_sim','lifton_identity_prior','lifton_gene','ref_gene',
            'contig','strand','lifton_exons','ref_exons'
        ]
        out.write('\t'.join(header) + '\n')
        for p in pairs:
            status = 'confident' if p.score >= conf_thresh else ('good' if p.score >= good_thresh else 'low')
            q = lifton.tx_index[p.q_id]; r = reference.tx_index[p.r_id]
            row = [
                p.q_id, p.r_id, f"{p.score:.6f}", status,
                f"{p.details['intron_sim']:.6f}", f"{p.details['jacc_internal']:.6f}", f"{p.details['jacc_all']:.6f}",
                f"{p.details['exon_count_sim']:.6f}", f"{p.details['boundary_sim']:.6f}", f"{p.details['lifton_identity_prior']:.6f}",
                q.gene_id, r.gene_id, q.contig, q.strand, str(q.exon_count()), str(r.exon_count())
            ]
            out.write('\t'.join(map(str,row)) + '\n')


def write_gene_pairs(gpairs: List[Tuple[str,str,float,float,int]], path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as out:
        header = ['lifton_gene','ref_gene','weighted_score','fraction_of_total','n_transcripts']
        out.write('\t'.join(header) + '\n')
        for qg, rg, w, frac, n in gpairs:
            out.write('\t'.join([qg, rg, f"{w:.6f}", f"{frac:.6f}", str(n)]) + '\n')

###############################################################################
# Reference GFF3 rewriting
###############################################################################

def rewrite_reference_gff3(reference_gff3: str, out_path: str, tx_map: Dict[str,str], gene_map: Dict[str,str], mode: str = 'alias'):
    """
    mode = 'alias'  -> keep original IDs; add lifton_* attributes on gene/mRNA
    mode = 'rename' -> replace ID with LiftOn id; update Parent on children
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Precompute ID rewrite mapping for 'rename' mode
    id_rewrite: Dict[str,str] = {}
    if mode == 'rename':
        id_rewrite.update(gene_map)
        id_rewrite.update(tx_map)

    def _rewrite_parents(val: str) -> str:
        # Parent can be comma-separated
        parts = val.split(',')
        for i,p in enumerate(parts):
            parts[i] = id_rewrite.get(p, p)
        return ','.join(parts)

    with _open_maybe_gzip(reference_gff3) as inp, open(out_path, 'w') as out:
        for line in inp:
            if line.startswith('#'):
                out.write(line)
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                out.write(line)
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_s = parts
            attrs = parse_gff3_attributes(attrs_s)
            ftype_l = ftype.lower()
            changed = False

            if ftype_l in ('gene', 'mrna', 'transcript', 'tRNA', 'rRNA'):
                idv = attrs.get('ID')
                if idv:
                    if ftype_l == 'gene' and idv in gene_map:
                        lifton_id = gene_map[idv]
                        if mode == 'alias':
                            attrs['lifton_gene_id'] = lifton_id
                            changed = True
                        else:
                            attrs['ID'] = lifton_id
                            id_rewrite[idv] = lifton_id
                            changed = True
                    elif ftype_l in ('mrna','transcript') and idv in tx_map:
                        lifton_id = tx_map[idv]
                        if mode == 'alias':
                            attrs['lifton_transcript_id'] = lifton_id
                            changed = True
                        else:
                            attrs['ID'] = lifton_id
                            id_rewrite[idv] = lifton_id
                            changed = True

                # Update Parent if rename mode and mapping exists
                if mode == 'rename' and 'Parent' in attrs:
                    new_parent = _rewrite_parents(attrs['Parent'])
                    if new_parent != attrs['Parent']:
                        attrs['Parent'] = new_parent
                        changed = True

            else:
                # child features — update Parent only in rename mode
                if mode == 'rename' and 'Parent' in attrs:
                    new_parent = _rewrite_parents(attrs['Parent'])
                    if new_parent != attrs['Parent']:
                        attrs['Parent'] = new_parent
                        changed = True

            if changed:
                parts[8] = format_gff3_attributes(attrs)
                out.write('\t'.join(parts) + '\n')
            else:
                out.write(line)

###############################################################################
# CLI
###############################################################################

def main():
    ap = argparse.ArgumentParser(description='Map LiftOn transcripts/genes onto a reference GFF3 and transfer IDs.')
    ap.add_argument('--lifton', required=True, help='LiftOn GFF3 (query)')
    ap.add_argument('--reference', required=True, help='Reference/Ensembl GFF3 (target)')
    ap.add_argument('--out-prefix', required=True, help='Prefix for output files')

    ap.add_argument('--window', type=int, default=100000, help='Window (bp) around LiftOn transcript span for candidate reference genes [default: 100000]')
    ap.add_argument('--topk', type=int, default=5, help='Keep top-K candidates per transcript before matching [default: 5]')
    ap.add_argument('--min-score', type=float, default=0.60, help='Minimum candidate score to consider [default: 0.60]')
    ap.add_argument('--good', type=float, default=0.75, help='Score threshold for status=good [default: 0.75]')
    ap.add_argument('--confident', type=float, default=0.85, help='Score threshold for status=confident [default: 0.85]')
    ap.add_argument('--gene-fraction', type=float, default=0.60, help='Min fraction of transcript weight to assign LiftOn gene to a reference gene [default: 0.60]')

    ap.add_argument('--rewrite-reference', action='store_true', help='Write mapped_reference.gff3 with LiftOn IDs (alias or rename mode)')
    ap.add_argument('--rename-mode', choices=['alias','rename'], default='alias', help='When rewriting reference GFF3: add aliases or actually rename IDs [default: alias]')

    args = ap.parse_args()

    lifton = load_gff3_as_annotation(args.lifton, 'LiftOn')
    reference = load_gff3_as_annotation(args.reference, 'Reference')

    pairs = compute_pairs(lifton, reference, window=args.window, min_candidate_score=args.min_score, topk=args.topk)

    tx_pairs_path = f"{args.out_prefix}.transcript_pairs.tsv"
    write_transcript_pairs(pairs, lifton, reference, tx_pairs_path, conf_thresh=args.confident, good_thresh=args.good)

    gpairs = aggregate_gene_pairs(pairs, lifton, reference, min_gene_fraction=args.gene_fraction)
    gene_pairs_path = f"{args.out_prefix}.gene_pairs.tsv"
    write_gene_pairs(gpairs, gene_pairs_path)

    if args.rewrite_reference:
        # Build maps restricted to pairs with adequate confidence
        tx_map: Dict[str,str] = {}
        for p in pairs:
            if p.score >= args.good:  # allow good+ to propagate
                tx_map[p.r_id] = p.q_id
        # Gene map from aggregated pairs
        gene_map: Dict[str,str] = {}
        for qg, rg, w, frac, n in gpairs:
            gene_map[rg] = qg

        out_gff = f"{args.out_prefix}.mapped_reference.gff3"
        rewrite_reference_gff3(args.reference, out_gff, tx_map, gene_map, mode=args.rename_mode)

    # Summary to stderr
    total_q = len([t for t in lifton.tx_index.values() if t.exons])
    mapped = len(pairs)
    sys.stderr.write(f"LiftOn transcripts with exons: {total_q}\n")
    sys.stderr.write(f"Mapped transcripts (one-to-one): {mapped}\n")
    if pairs:
        sys.stderr.write(f"Top score: {pairs[0].score:.3f}, median: {sorted(p.score for p in pairs)[len(pairs)//2]:.3f}\n")
    sys.stderr.write(f"Gene pairs: {len(gpairs)} (fraction threshold {args.gene_fraction})\n")

if __name__ == '__main__':
    main()

# aligner.py
from __future__ import annotations
import os, shutil, subprocess, tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from fasta_io import FastaReader, revcomp
from gff_io import Gene, Transcript

@dataclass
class Aln:
    qname: str
    rname: str
    pos: int        # 1-based start on target
    cigar: str
    mapq: int
    is_rev: bool
    nm: Optional[int]  # edit distance if present
    length_q: int

@dataclass
class GeneMap:
    gid: str
    rname: Optional[str]
    start: Optional[int]
    end: Optional[int]
    strand: Optional[str]
    identity: float
    coverage: float
    cigar: Optional[str]

def _parse_sam_line(ln: str) -> Optional[Aln]:
    if ln.startswith("@"): return None
    p = ln.rstrip("\n").split("\t")
    if len(p) < 11: return None
    qname, flag, rname, pos, mapq, cigar = p[0], int(p[1]), p[2], int(p[3]), int(p[4]), p[5]
    if rname == "*": return None
    tags = { t.split(":",2)[0]: t.split(":",2)[2] for t in p[11:] if ":" in t }
    nm = int(tags["NM"]) if "NM" in tags else None
    is_rev = bool(flag & 16)
    lq = int(tags["ql"]) if "ql" in tags else 0
    return Aln(qname, rname, pos, cigar, mapq, is_rev, nm, lq)

def _write_gene_fasta(tmpfa: Path, ref_fa: FastaReader, genes: Dict[str, Gene]) -> Dict[str, Tuple[str,int,int,bool]]:
    """
    Write per-gene sequences (always forward in reference genomic orientation).
    Returns map gid -> (seqid, start, end, gene_on_minus_strand)
    """
    m: Dict[str, Tuple[str,int,int,bool]] = {}
    with tmpfa.open("w") as fh:
        for g in genes.values():
            start, end = (g.start, g.end)
            seq = ref_fa.slice(g.seqid, start, end)
            if g.strand == "-":
                # keep query in forward reference orientation: DO NOT revcomp here
                pass
            fh.write(f">{g.gid}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
            m[g.gid] = (g.seqid, start, end, g.strand == "-")
    return m

def run_minimap2(query_fa: Path, target_fa: Path, threads: int = 4, extra: Optional[List[str]] = None) -> List[Aln]:
    if shutil.which("minimap2") is None:
        return []  # handled by fallback
    cmd = ["minimap2", "-a", "--eqx", "--end-bonus", "5", "-N", "50", "-p", "0.5", "-t", str(threads)]
    if extra: cmd.extend(extra)
    cmd.extend([str(target_fa), str(query_fa)])
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
    alns: List[Aln] = []
    for ln in proc.stdout.splitlines():
        al = _parse_sam_line(ln)
        if al: alns.append(al)
    return alns

def fallback_exact_align(ref_fa: FastaReader, tgt_fa: FastaReader, genes: Dict[str, Gene]) -> List[Aln]:
    out: List[Aln] = []
    for g in genes.values():
        seq = ref_fa.slice(g.seqid, g.start, g.end)
        for chrom in tgt_fa.chromosomes:
            s = tgt_fa.get(chrom)
            i = s.find(seq)
            is_rev = False
            if i == -1:
                rc = revcomp(seq)
                i = s.find(rc)
                is_rev = i != -1
            if i != -1:
                out.append(Aln(g.gid, chrom, i+1, f"{len(seq)}M", 60, is_rev, 0, len(seq)))
                break
    return out

def align_genes(ref_fasta: Path, ref_genes: Dict[str, Gene], target_fasta: Path, threads: int = 4) -> Dict[str, GeneMap]:
    ref = FastaReader(ref_fasta)
    tgt = FastaReader(target_fasta)
    with tempfile.TemporaryDirectory() as td:
        qfa = Path(td) / "genes.fa"
        meta = _write_gene_fasta(qfa, ref, ref_genes)
        alns = run_minimap2(qfa, Path(target_fasta), threads=threads)
        if not alns:
            alns = fallback_exact_align(ref, tgt, ref_genes)
    # choose best per gene (highest mapq; then shortest NM; then coverage)
    best: Dict[str, Aln] = {}
    for a in alns:
        cur = best.get(a.qname)
        if cur is None or (a.mapq, -(cur.nm or 1<<30)) < (a.mapq, -(a.nm or 1<<30)):
            best[a.qname] = a
    # convert to GeneMap (identity approx)
    gmaps: Dict[str, GeneMap] = {}
    for gid, aln in best.items():
        qlen = ref_genes[gid].end - ref_genes[gid].start + 1
        mm = aln.nm or 0
        ident = max(0.0, 1.0 - (mm / max(1, qlen)))
        strand = "-" if aln.is_rev else "+"
        gmaps[gid] = GeneMap(gid, aln.rname, aln.pos, aln.pos + qlen - 1, strand, ident, 1.0, aln.cigar)
    return gmaps

# gff_io.py
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Iterable, Tuple

def parse_attrs(s: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for kv in filter(None, s.split(";")):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out

def fmt_attrs(d: Dict[str, str]) -> str:
    return ";".join(f"{k}={v}" for k, v in d.items())

@dataclass
class Exon:
    start: int
    end: int

@dataclass
class CDS:
    start: int
    end: int
    phase: int = 0

@dataclass
class Transcript:
    tid: str
    seqid: str
    start: int
    end: int
    strand: str
    attrs: Dict[str, str] = field(default_factory=dict)
    exons: List[Exon] = field(default_factory=list)
    cdss: List[CDS] = field(default_factory=list)

@dataclass
class Gene:
    gid: str
    seqid: str
    start: int
    end: int
    strand: str
    attrs: Dict[str, str] = field(default_factory=dict)
    transcripts: Dict[str, Transcript] = field(default_factory=dict)

def load_gff(path: str | Path) -> Dict[str, Gene]:
    """
    Minimal, Ensembl-style GFF3 parser (gene/mRNA/exon/CDS). Keeps hierarchy.
    """
    genes: Dict[str, Gene] = {}
    path = Path(path)
    with path.open() as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9: 
                continue
            seqid, _src, ftype, s, e, _score, strand, phase, attr = p
            s, e = int(s), int(e)
            a = parse_attrs(attr)
            fid = a.get("ID")
            parent = a.get("Parent")
            if ftype == "gene" and fid:
                genes[fid] = Gene(fid, seqid, s, e, strand, a, {})
            elif ftype in ("mRNA", "transcript") and fid and parent and parent in genes:
                genes[parent].transcripts[fid] = Transcript(fid, seqid, s, e, strand, a)
            elif ftype == "exon" and parent:
                for pid in parent.split(","):
                    for g in genes.values():
                        tx = g.transcripts.get(pid)
                        if tx:
                            tx.exons.append(Exon(s, e))
                            break
            elif ftype == "CDS" and parent:
                for pid in parent.split(","):
                    for g in genes.values():
                        tx = g.transcripts.get(pid)
                        if tx:
                            cd = CDS(s, e, 0 if phase == "." else int(phase))
                            tx.cdss.append(cd)
                            break
    # canonical ordering
    for g in genes.values():
        for tx in g.transcripts.values():
            tx.exons.sort(key=lambda x: x.start)
            tx.cdss.sort(key=lambda x: x.start)
    return genes

def write_gff(genes: Dict[str, Gene], out: str | Path) -> None:
    out = Path(out)
    rows: List[Tuple] = []
    for g in genes.values():
        rows.append((g.seqid, "idmapper", "gene", g.start, g.end, ".", g.strand, ".", fmt_attrs(g.attrs)))
        for tx in g.transcripts.values():
            rows.append((tx.seqid, "idmapper", "mRNA", tx.start, tx.end, ".", tx.strand, ".", fmt_attrs(tx.attrs)))
            for i, ex in enumerate(tx.exons, 1):
                attrs = {"ID": f"{tx.tid}.exon{i}", "Parent": tx.tid}
                rows.append((tx.seqid, "idmapper", "exon", ex.start, ex.end, ".", tx.strand, ".", fmt_attrs(attrs)))
            for i, cd in enumerate(tx.cdss, 1):
                attrs = {"ID": f"{tx.tid}.cds{i}", "Parent": tx.tid}
                rows.append((tx.seqid, "idmapper", "CDS", cd.start, cd.end, ".", tx.strand, str(cd.phase), fmt_attrs(attrs)))
    rows.sort(key=lambda r: (r[0], r[3], r[2]))
    with out.open("w") as fh:
        fh.write("##gff-version 3\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

"""Independent GFF3 / GTF / Tiberius-hybrid reader.

Packaged unchanged in substance from the 2026-10-05 benchmark (scripts/indep_gff.py).
It deliberately shares no code with the comparator parsers so that it can cross-check
them. Coordinates stay one-based inclusive; genes and transcripts are keyed by
(seqname, identity field), so identity reuse on different sequences stays distinct.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter, defaultdict

GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
CHILD_TYPES = {"exon", "CDS"}
# Everything with a gene Parent is a transcript in GFF3 (independent of any SO list).

_GTF_KV = re.compile(r'(\S+)\s+"([^"]*)"')
_GFF_PREFIX = re.compile(r"^(gene|transcript|mRNA):")


def open_text(path):
    with open(path, "rb") as h:
        gz = h.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path, "rt")


def gff3_attrs(text: str) -> dict:
    out = {}
    for part in text.strip().strip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def gtf_attrs(text: str) -> dict:
    return {k: v for k, v in _GTF_KV.findall(text)}


def detect(path: str, n: int = 300) -> str:
    gff = gtf = bare = 0
    with open_text(path) as h:
        seen = 0
        for line in h:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            seen += 1
            a = p[8].strip()
            if "ID=" in a or "Parent=" in a:
                gff += 1
            if 'gene_id "' in a or 'transcript_id "' in a:
                gtf += 1
            if (
                p[2] in ("gene", "transcript")
                and re.fullmatch(r"\S+", a)
                and "=" not in a
                and '"' not in a
            ):
                bare += 1
            if seen >= n:
                break
    if bare and gtf:
        return "tiberius"
    return "gff3" if gff >= gtf else "gtf"


class Annotation:
    def __init__(self):
        self.genes: dict = {}
        self.tx: dict = {}
        self.exons: dict = defaultdict(list)
        self.cds: dict = defaultdict(list)
        self.feature_counts: Counter = Counter()
        self.seq_genes: Counter = Counter()
        self.seq_max_end: dict = {}
        self.problems: Counter = Counter()
        self.examples: dict = defaultdict(list)
        self.phase_values: Counter = Counter()
        self.child_phase_on_exon: int = 0
        self.multi_parent_children = 0
        self.duplicate_rows = 0
        self.identity_reuse_same_seq: Counter = Counter()  # level -> ids
        self.identity_reuse_cross_seq: Counter = Counter()
        self.name_reuse: Counter = Counter()
        self.ids_by_level: dict = {
            "gene": defaultdict(set),
            "transcript": defaultdict(set),
        }
        self.names_by_level: dict = {"gene": Counter(), "transcript": Counter()}
        self.format = ""

    def note(self, key, example):
        self.problems[key] += 1
        if len(self.examples[key]) < 5:
            self.examples[key].append(example)


def read(path: str, fmt: str | None = None, keep_gene=None) -> Annotation:
    """Read an annotation. keep_gene(attrs, ftype) -> bool limits stored genes (references)."""
    fmt = fmt or detect(path)
    ann = Annotation()
    ann.format = fmt
    if fmt == "gff3":
        _read_gff3(path, ann, keep_gene)
    else:
        _read_gtf(path, ann, fmt == "tiberius")
    _finalise(ann)
    return ann


def _track_seq(ann, chrom, end):
    if end > ann.seq_max_end.get(chrom, 0):
        ann.seq_max_end[chrom] = end


def _read_gff3(path, ann, keep_gene):
    # Pass 1: genes and transcripts (anything whose Parent is a gene).
    raw_tx = []
    seen_lines = set()
    with open_text(path) as h:
        for line in h:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, ftype = p[0], p[2]
            start, end = int(p[3]), int(p[4])
            ann.feature_counts[ftype] += 1
            _track_seq(ann, chrom, end)
            if (
                ftype in GENE_TYPES
                or ftype in ("mRNA", "transcript")
                or "Parent=gene:" in p[8]
                or ftype.endswith("RNA")
                or ftype.endswith("_segment")
                or ftype
                in (
                    "pseudogenic_transcript",
                    "unconfirmed_transcript",
                    "processed_transcript",
                )
            ):
                a = gff3_attrs(p[8])
                gid = _GFF_PREFIX.sub("", a.get("ID", ""))
                if ftype in GENE_TYPES:
                    if keep_gene is not None and not keep_gene(a, ftype):
                        ann.genes[(chrom, gid)] = None  # known but not stored
                        continue
                    key = (chrom, gid)
                    if key in ann.genes and ann.genes[key] is not None:
                        ann.identity_reuse_same_seq["gene"] += 1
                        ann.note("gene_id_reused_same_seq", f"{chrom}:{gid}")
                    ann.genes[key] = dict(
                        chrom=chrom,
                        strand=p[6],
                        start=start,
                        end=end,
                        id=gid,
                        name=a.get("Name", ""),
                        biotype=a.get("biotype", a.get("gene_biotype", "")),
                        ftype=ftype,
                    )
                    ann.ids_by_level["gene"][gid].add(chrom)
                    if a.get("Name"):
                        ann.names_by_level["gene"][a["Name"]] += 1
                    ann.seq_genes[chrom] += 1
                elif "Parent" in a:
                    raw_tx.append((chrom, p[6], start, end, ftype, a))
    for chrom, strand, start, end, ftype, a in raw_tx:
        tid = _GFF_PREFIX.sub("", a.get("ID", ""))
        parents = [_GFF_PREFIX.sub("", x) for x in a["Parent"].split(",")]
        if (chrom, parents[0]) not in ann.genes:
            if ftype not in CHILD_TYPES:
                ann.note(
                    "transcript_parent_not_gene_on_same_seq",
                    f"{chrom}:{a.get('ID','')}->{parents[0]}",
                )
            continue
        if ann.genes[(chrom, parents[0])] is None:
            ann.tx[(chrom, tid)] = None
            continue
        if len(parents) > 1:
            ann.note("transcript_multi_parent", f"{chrom}:{tid}")
        key = (chrom, tid)
        if key in ann.tx and ann.tx[key] is not None:
            ann.identity_reuse_same_seq["transcript"] += 1
            ann.note("transcript_id_reused_same_seq", f"{chrom}:{tid}")
        ann.tx[key] = dict(
            gene=(chrom, parents[0]),
            chrom=chrom,
            strand=strand,
            start=start,
            end=end,
            id=tid,
            ftype=ftype,
            biotype=a.get("biotype", a.get("transcript_biotype", "")),
            tags=a.get("tag", ""),
        )
        ann.ids_by_level["transcript"][tid].add(chrom)
        if a.get("Name"):
            ann.names_by_level["transcript"][a["Name"]] += 1
    # Pass 2: exons and CDS.
    with open_text(path) as h:
        for line in h:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in CHILD_TYPES:
                continue
            hl = hash(line)
            if hl in seen_lines:
                ann.duplicate_rows += 1
                ann.note("identical_duplicate_child_row", line.strip()[:120])
            seen_lines.add(hl)
            chrom, start, end = p[0], int(p[3]), int(p[4])
            a = gff3_attrs(p[8])
            parents = [
                _GFF_PREFIX.sub("", x) for x in a.get("Parent", "").split(",") if x
            ]
            if len(parents) > 1:
                ann.multi_parent_children += 1
            if p[2] == "exon" and p[7] not in (".", ""):
                ann.child_phase_on_exon += 1
            for parent in parents:
                key = (chrom, parent)
                if key not in ann.tx:
                    ann.note("orphan_child", f"{chrom}:{parent}")
                    continue
                if ann.tx[key] is None:
                    continue
                if p[2] == "exon":
                    ann.exons[key].append((start, end))
                else:
                    ann.cds[key].append((start, end, p[7]))
                    ann.phase_values[p[7]] += 1
    _dup_child_rows(ann)


def _dup_child_rows(ann):
    """Identical repeated exon/CDS intervals within one transcript."""
    for store, label in ((ann.exons, "exon"), (ann.cds, "CDS")):
        for key, ivs in store.items():
            plain = [(s, e) for s, e, *_ in ivs]
            if len(plain) != len(set(plain)):
                ann.note(
                    f"repeated_{label}_interval_in_transcript", f"{key[0]}:{key[1]}"
                )


def _read_gtf(path, ann, tiberius):
    gene_rows = {}
    with open_text(path) as h:
        for line in h:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, ftype, start, end, strand = p[0], p[2], int(p[3]), int(p[4]), p[6]
            ann.feature_counts[ftype] += 1
            _track_seq(ann, chrom, end)
            raw = p[8].strip()
            a = gtf_attrs(raw)
            if tiberius and ftype in ("gene", "transcript") and not a:
                ident = raw.rstrip(";").strip()
                if ftype == "gene":
                    a = {"gene_id": ident}
                else:
                    a = {"transcript_id": ident, "gene_id": ident.rsplit(".", 1)[0]}
            gid, tid = a.get("gene_id", ""), a.get("transcript_id", "")
            if ftype == "gene":
                key = (chrom, gid)
                if key in gene_rows:
                    ann.identity_reuse_same_seq["gene"] += 1
                    ann.note("gene_id_reused_same_seq", f"{chrom}:{gid}")
                gene_rows[key] = dict(
                    chrom=chrom,
                    strand=strand,
                    start=start,
                    end=end,
                    id=gid,
                    name=a.get("gene_name", ""),
                    biotype=a.get("gene_biotype", a.get("gene_type", "")),
                    ftype=ftype,
                )
                if a.get("gene_name"):
                    ann.names_by_level["gene"][a["gene_name"]] += 1
            elif ftype == "transcript":
                key = (chrom, tid)
                if key in ann.tx:
                    ann.identity_reuse_same_seq["transcript"] += 1
                    ann.note("transcript_id_reused_same_seq", f"{chrom}:{tid}")
                ann.tx[key] = dict(
                    gene=(chrom, gid),
                    chrom=chrom,
                    strand=strand,
                    start=start,
                    end=end,
                    id=tid,
                    ftype=ftype,
                    biotype=a.get("transcript_biotype", a.get("transcript_type", "")),
                    tags=a.get("tag", ""),
                    explicit=True,
                )
            elif ftype in CHILD_TYPES:
                key = (chrom, tid)
                if key not in ann.tx:
                    ann.tx[key] = dict(
                        gene=(chrom, gid),
                        chrom=chrom,
                        strand=strand,
                        start=start,
                        end=end,
                        id=tid,
                        ftype="inferred",
                        biotype="",
                        tags="",
                        explicit=False,
                    )
                t = ann.tx[key]
                if t["gene"] != (chrom, gid) and gid:
                    ann.note("child_gene_id_conflict", f"{chrom}:{tid}")
                if t["strand"] != strand:
                    ann.note("child_strand_conflict", f"{chrom}:{tid}")
                if not t.get("explicit"):
                    t["start"], t["end"] = min(t["start"], start), max(t["end"], end)
                if ftype == "exon":
                    ann.exons[key].append((start, end))
                    if p[7] not in (".", ""):
                        ann.child_phase_on_exon += 1
                else:
                    ann.cds[key].append((start, end, p[7]))
                    ann.phase_values[p[7]] += 1
    # genes: explicit rows, else inferred from transcripts
    for key, t in ann.tx.items():
        g = t["gene"]
        if g not in gene_rows:
            gene_rows[g] = dict(
                chrom=t["chrom"],
                strand=t["strand"],
                start=t["start"],
                end=t["end"],
                id=g[1],
                name="",
                biotype="",
                ftype="inferred",
            )
        else:
            gr = gene_rows[g]
            if gr["ftype"] == "inferred":
                gr["start"], gr["end"] = min(gr["start"], t["start"]), max(
                    gr["end"], t["end"]
                )
    ann.genes = gene_rows
    for chrom, gid in gene_rows:
        ann.ids_by_level["gene"][gid].add(chrom)
        ann.seq_genes[chrom] += 1
    for chrom, tid in ann.tx:
        ann.ids_by_level["transcript"][tid].add(chrom)
    _dup_child_rows(ann)


def _finalise(ann):
    for level in ("gene", "transcript"):
        ann.identity_reuse_cross_seq[level] = sum(
            1 for chroms in ann.ids_by_level[level].values() if len(chroms) > 1
        )
        ann.name_reuse[level] = sum(
            1 for c in ann.names_by_level[level].values() if c > 1
        )
    for key in list(ann.exons):
        ann.exons[key] = sorted(ann.exons[key])
    for key in list(ann.cds):
        ann.cds[key] = sorted(ann.cds[key])


def cds_chain(ivs) -> tuple:
    return tuple((s, e) for s, e, *_ in ivs)


def intron_chain(ivs) -> tuple:
    ivs = [(s, e) for s, e, *_ in ivs]
    return tuple((ivs[i][1] + 1, ivs[i + 1][0] - 1) for i in range(len(ivs) - 1))


def genes_transcripts(ann):
    by_gene = defaultdict(list)
    for tkey, t in ann.tx.items():
        if t is None:
            continue
        by_gene[t["gene"]].append(tkey)
    return by_gene

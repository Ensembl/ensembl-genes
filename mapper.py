# mapper.py
from __future__ import annotations
from typing import Dict, Tuple, Optional
from collections import Counter
import re

from gff_io import Gene, Transcript, Exon, CDS, load_gff, write_gff
from id_manager import IDManager
from aligner import align_genes, GeneMap
from cigar_map import project_ref_interval
from fasta_io import FastaReader, revcomp
from synteny import synteny_disambiguate


def _flip_strand(s: str) -> str:
    return "+" if s == "-" else "-"

# --- helpers for namespace and species-specific core prefix inference ---

_CORE_RE = re.compile(r"^[A-Za-z]+(?=\d+)")  # letters before the first digit

def _extract_core_prefix(stable_id: str) -> Optional[str]:
    core = stable_id.split(":", 1)[-1]  # strip namespace if present
    m = _CORE_RE.match(core)
    return m.group(0) if m else None

def _derive_type_prefix(base_prefix: str, kind: str) -> str:
    """
    Given a species base prefix (e.g. ENSCGRG), derive a type-specific prefix:
      gene -> keep as-is, transcript -> ...T, exon -> ...E, protein -> ...P
    If base ends with one of G/T/E/P, replace that last letter; otherwise append.
    """
    if not base_prefix:
        return base_prefix
    target_letter = {"gene":"G", "transcript":"T", "exon":"E", "protein":"P"}[kind]
    if base_prefix[-1] in "GTEP":
        return base_prefix[:-1] + target_letter
    return base_prefix + target_letter

def _infer_core_prefix_for(kind: str, tgt_genes: dict, fallback_from: Optional[str] = None) -> str:
    """
    Infer a species-specific core prefix for 'gene'/'transcript'/'exon'/'protein'.
    Safely handles Exon/CDS that lack attrs. If nothing is discoverable and
    fallback_from is provided, derive from that (e.g., exon/protein from gene).
    Otherwise raise.
    """
    seen: list[str] = []

    if kind == "gene":
        for g in tgt_genes.values():
            sid = g.attrs.get("ID") or g.gid
            p = _extract_core_prefix(sid) if isinstance(sid, str) else None
            if p: seen.append(p)

    elif kind == "transcript":
        for g in tgt_genes.values():
            for t in g.transcripts.values():
                sid = t.attrs.get("ID") or t.tid
                p = _extract_core_prefix(sid) if isinstance(sid, str) else None
                if p: seen.append(p)

    elif kind == "exon":
        # Try to learn from exon attributes when available
        for g in tgt_genes.values():
            for t in g.transcripts.values():
                for ex in getattr(t, "exons", []):
                    ex_attrs = getattr(ex, "attrs", None)
                    if isinstance(ex_attrs, dict):
                        sid = ex_attrs.get("exon_id") or ex_attrs.get("Name") or ex_attrs.get("ID")
                        if isinstance(sid, str):
                            p = _extract_core_prefix(sid)
                            if p: seen.append(p)

    elif kind == "protein":
        # Try to learn from CDS/protein attributes when available
        for g in tgt_genes.values():
            for t in g.transcripts.values():
                for cd in getattr(t, "cdss", []):
                    cd_attrs = getattr(cd, "attrs", None)
                    if isinstance(cd_attrs, dict):
                        sid = cd_attrs.get("protein_id") or cd_attrs.get("ID")
                        sid_core = sid.split(":",1)[-1] if isinstance(sid, str) else None
                        if isinstance(sid_core, str):
                            p = _extract_core_prefix(sid_core)
                            if p: seen.append(p)

    if seen:
        return Counter(seen).most_common(1)[0][0]

    # Nothing discoverable: derive from fallback base if provided
    if fallback_from:
        return _derive_type_prefix(fallback_from, kind)

    raise ValueError(f"Could not infer a species-specific core prefix for {kind} IDs from the target GFF3.")


def _infer_namespace_for(kind: str, genes: dict, default_token: str) -> str:
    """
    Infer namespace token (bit before ':') from existing IDs, safely handling
    cases where Exon/CDS objects have no attrs. Falls back to default_token.
    """
    if kind == "gene":
        for g in genes.values():
            sid = g.attrs.get("ID") or g.gid
            if isinstance(sid, str) and ":" in sid:
                return sid.split(":", 1)[0]
        return default_token

    if kind == "transcript":
        for g in genes.values():
            for t in g.transcripts.values():
                sid = t.attrs.get("ID") or t.tid
                if isinstance(sid, str) and ":" in sid:
                    return sid.split(":", 1)[0]
        return default_token

    if kind == "cds":
        # Many Ensembl GFF3s use ID=CDS:ENSP... on CDS rows; if not present in your
        # parsed objects, just fall back to the default.
        for g in genes.values():
            for t in g.transcripts.values():
                for cd in getattr(t, "cdss", []):
                    sid = getattr(cd, "id", None)
                    if isinstance(sid, str) and ":" in sid:
                        return sid.split(":", 1)[0]
                    cd_attrs = getattr(cd, "attrs", None)
                    if isinstance(cd_attrs, dict):
                        sid2 = cd_attrs.get("ID")
                        if isinstance(sid2, str) and ":" in sid2:
                            return sid2.split(":", 1)[0]
        return default_token

    if kind == "exon":
        # Ensembl exons usually don't have namespaced IDs; they use exon_id/Name.
        # We only attempt to discover a namespace if an exon carries ID=exon:...
        for g in genes.values():
            for t in g.transcripts.values():
                for ex in getattr(t, "exons", []):
                    ex_attrs = getattr(ex, "attrs", None)
                    if isinstance(ex_attrs, dict):
                        sid = ex_attrs.get("ID")
                        if isinstance(sid, str) and ":" in sid:
                            return sid.split(":", 1)[0]
        return default_token  # typically unused for exons

    return default_token


def _get_version(attrs: Dict[str, str]) -> int:
    v = attrs.get("version")
    try:
        return int(v) if v is not None else 1
    except ValueError:
        return 1


def _set_version(attrs: Dict[str, str], v: int) -> None:
    attrs["version"] = str(v)


# --- robust projection helpers (CIGAR + fallback exact search) ---

def _search_interval_in_window(
    exon_seq: str,
    target_seq: str,
) -> Optional[Tuple[int, int, str]]:
    """
    Try to find exon_seq or its reverse-complement inside target_seq.
    Returns (start_1_based, end_1_based, strand) relative to the *window string*.
    """
    idx = target_seq.find(exon_seq)
    if idx != -1:
        s = idx + 1
        e = s + len(exon_seq) - 1
        return s, e, "+"
    rc = revcomp(exon_seq)
    idx = target_seq.find(rc)
    if idx != -1:
        s = idx + 1
        e = s + len(exon_seq) - 1
        return s, e, "-"
    return None


def _project_gene_and_children(
    ref_gene: Gene,
    gm: GeneMap,
    ref_reader: FastaReader,
    tgt_reader: FastaReader,
    pad: int = 10000
) -> Optional[Gene]:
    """
    Use the gene-level alignment to project transcript/exon/CDS coordinates.
    Prefer fast CIGAR-based projection; if an interval is not fully covered,
    fall back to exact search of the exon/CDS sequence within a padded window
    around the mapped gene on the target contig.

    Returns a new Gene in target space, or None if projection fails.
    """
    if gm.rname is None or gm.cigar is None or gm.start is None or gm.end is None:
        return None

    # Build a new Gene in target space; copy attrs so we preserve Ensembl fields
    t_strand = ref_gene.strand if gm.strand == "+" else _flip_strand(ref_gene.strand)
    out = Gene(
        gid=ref_gene.gid,
        seqid=gm.rname,
        start=gm.start,
        end=gm.end,
        strand=t_strand,
        attrs=dict(ref_gene.attrs),
        transcripts={}
    )

    # Prepare target window for fallback searches
    full_tgt = tgt_reader.get(gm.rname)
    win_s = max(1, gm.start - pad)
    win_e = min(len(full_tgt), gm.end + pad)
    tgt_window = full_tgt[win_s - 1: win_e]  # 1-based to 0-based

    # For each transcript, project exon/CDS
    g_start = ref_gene.start
    any_fail = False

    for tid, tx in ref_gene.transcripts.items():
        # copy transcript attrs
        ttx = Transcript(
            tid=tx.tid,
            seqid=gm.rname,
            start=gm.start,
            end=gm.end,
            strand=t_strand,
            attrs=dict(tx.attrs),
        )
        min_s, max_e = gm.end, gm.start

        # exons (preserve/clone attrs if available)
        for ex in tx.exons:
            ex_attrs = dict(getattr(ex, "attrs", {}))  # may be empty if gff_io doesn't store them
            q_s = ex.start - g_start + 1
            q_e = ex.end   - g_start + 1
            try:
                t_s, t_e = project_ref_interval(gm.cigar, gm.start, q_s, q_e)
                if t_s > t_e:
                    t_s, t_e = t_e, t_s
            except RuntimeError:
                # Fallback: search exon sequence in target window
                exon_seq = ref_reader.slice(ref_gene.seqid, ex.start, ex.end)
                hit = _search_interval_in_window(exon_seq, tgt_window)
                if not hit:
                    any_fail = True
                    break
                w_s, w_e, _str = hit
                t_s = win_s + (w_s - 1)
                t_e = win_s + (w_e - 1)
                if t_s > t_e:
                    t_s, t_e = t_e, t_s
            new_ex = Exon(t_s, t_e)
            # attach attrs back if Exon supports it
            if hasattr(new_ex, "attrs"):
                new_ex.attrs.update(ex_attrs)
            ttx.exons.append(new_ex)
            min_s, max_e = min(min_s, t_s), max(max_e, t_e)
        if any_fail:
            break

        # CDS blocks
        for cd in tx.cdss:
            cd_attrs = dict(getattr(cd, "attrs", {}))
            q_s = cd.start - g_start + 1
            q_e = cd.end   - g_start + 1
            try:
                t_s, t_e = project_ref_interval(gm.cigar, gm.start, q_s, q_e)
                if t_s > t_e:
                    t_s, t_e = t_e, t_s
            except RuntimeError:
                cds_seq = ref_reader.slice(ref_gene.seqid, cd.start, cd.end)
                hit = _search_interval_in_window(cds_seq, tgt_window)
                if not hit:
                    any_fail = True
                    break
                w_s, w_e, _str = hit
                t_s = win_s + (w_s - 1)
                t_e = win_s + (w_e - 1)
                if t_s > t_e:
                    t_s, t_e = t_e, t_s
            new_cd = CDS(t_s, t_e, cd.phase if t_strand == ref_gene.strand else 0)
            if hasattr(new_cd, "attrs"):
                new_cd.attrs.update(cd_attrs)
            ttx.cdss.append(new_cd)
            min_s, max_e = min(min_s, t_s), max(max_e, t_e)

        if any_fail:
            break

        ttx.start, ttx.end = min_s, max_e
        out.transcripts[tid] = ttx

    if any_fail or not out.transcripts:
        return None

    out.start = min(min(tx.start for tx in out.transcripts.values()), out.start)
    out.end   = max(max(tx.end for tx in out.transcripts.values()), out.end)
    return out


def map_ids(
    ref_fasta: str, ref_gff: str, tgt_fasta: str, tgt_gff: str,
    out_gff: str, report_txt: str, threads: int = 8,
    identity_min: float = 0.80
) -> None:
    ref = load_gff(ref_gff)
    tgt = load_gff(tgt_gff)
    idm = IDManager(ref, tgt)

    # readers for fallback search
    ref_reader = FastaReader(ref_fasta)
    tgt_reader = FastaReader(tgt_fasta)

    # infer species-specific prefixes and namespaces from TARGET
    ns_gene = _infer_namespace_for("gene", tgt, "gene")
    ns_tx   = _infer_namespace_for("transcript", tgt, "transcript")
    ns_cds  = _infer_namespace_for("cds", tgt, "CDS")

    gene_core = _infer_core_prefix_for("gene", tgt)                          # e.g., ENSCGRG
    tx_core   = _infer_core_prefix_for("transcript", tgt, gene_core)         # derive ...T if needed
    ex_core   = _infer_core_prefix_for("exon", tgt, gene_core)               # derive ...E if needed
    prot_core = _infer_core_prefix_for("protein", tgt, gene_core)            # derive ...P if needed

    # 1) Align genes (minimap2, fallback exact), then synteny disambiguate
    gmaps = align_genes(ref_fasta, ref, tgt_fasta, threads=threads)
    gmaps = synteny_disambiguate(ref, gmaps)

    # 2) Project features into target coordinates; decide mapped/missing/new
    updated: Dict[str, Gene] = {}
    mapped_genes: set[str] = set()

    for gid, refgene in ref.items():
        gm = gmaps.get(gid)
        if gm and gm.identity >= identity_min and gm.rname:
            # project robustly
            projected = _project_gene_and_children(refgene, gm, ref_reader, tgt_reader)
            if projected is None:
                # cannot project all parts => treat as missing
                idm.set(refgene.gid, refgene.gid, "missing")
                continue

            # --- MAPPED: keep IDs; bump version= attribute where present ---
            # Gene
            g_attrs = projected.attrs
            g_old_v = _get_version(g_attrs)
            _set_version(g_attrs, g_old_v + 1)
            # Ensure ID (namespaced) and gene_id (core) are consistent
            g_id_ns = refgene.attrs.get("ID", refgene.gid)
            g_core  = g_id_ns.split(":", 1)[-1]
            g_attrs["ID"] = g_id_ns
            g_attrs["gene_id"] = g_core

            # Transcripts
            for tid, tx in projected.transcripts.items():
                t_attrs = tx.attrs
                t_old_v = _get_version(t_attrs)
                _set_version(t_attrs, t_old_v + 1)
                t_id_ns = refgene.transcripts[tid].attrs.get("ID", tid)
                t_core  = t_id_ns.split(":", 1)[-1]
                t_attrs["ID"] = t_id_ns
                t_attrs["transcript_id"] = t_core
                t_attrs["Parent"] = g_id_ns

                # Exons: increment version if present; preserve exon_id/Name if available
                for ex in getattr(tx, "exons", []):
                    if hasattr(ex, "attrs"):
                        ex_old_v = _get_version(ex.attrs)
                        if "version" in ex.attrs:
                            _set_version(ex.attrs, ex_old_v + 1)
                        # Parent should be transcript namespaced ID
                        ex.attrs["Parent"] = t_id_ns

                # CDS: keep protein_id and ID (CDS namespace) as-is (Ensembl CDS lines usually no version)
                for cd in getattr(tx, "cdss", []):
                    if hasattr(cd, "attrs"):
                        cd.attrs["Parent"] = t_id_ns

            updated[g_id_ns] = projected
            idm.set(refgene.gid, g_id_ns, "mapped")
            mapped_genes.add(gid)
        else:
            idm.set(refgene.gid, refgene.gid, "missing")

    # 3) For target genes without a mapped ID, create NEW stable IDs (version=1)
    for gid, tgene in tgt.items():
        # skip ones we already replaced via mapping (by overlap on seqid,start)
        overlaps = False
        for pg in updated.values():
            if pg.seqid == tgene.seqid and not (pg.end < tgene.start or tgene.end < pg.start):
                overlaps = True
                break
        if overlaps:
            continue

        # new core IDs (species-specific prefixes)
        new_gene_core = idm.new_id(gene_core).split(":", 1)[-1]  # get core generated by IDManager
        # But IDManager.new_id returns ns:core if ns provided; we want just the core here
        new_gene_core = new_gene_core if "." not in new_gene_core else new_gene_core.split(".",1)[0]

        ng = Gene(
            gid=tgene.gid, seqid=tgene.seqid, start=tgene.start, end=tgene.end,
            strand=tgene.strand, attrs=dict(tgene.attrs), transcripts={}
        )

        # Set gene attrs: ID namespaced, version=1, gene_id core
        ng.attrs["ID"] = f"{ns_gene}:{new_gene_core}"
        ng.attrs["gene_id"] = new_gene_core
        ng.attrs["version"] = "1"

        for tid, ttx in tgene.transcripts.items():
            # new transcript core
            new_tx_core = idm.new_id(tx_core).split(":", 1)[-1]
            new_tx_core = new_tx_core if "." not in new_tx_core else new_tx_core.split(".",1)[0]

            ntx = Transcript(
                tid=ttx.tid, seqid=ttx.seqid, start=ttx.start, end=ttx.end,
                strand=ttx.strand, attrs=dict(ttx.attrs),
                exons=list(getattr(ttx, "exons", [])), cdss=list(getattr(ttx, "cdss", []))
            )
            ntx.attrs["ID"] = f"{ns_tx}:{new_tx_core}"
            ntx.attrs["transcript_id"] = new_tx_core
            ntx.attrs["Parent"] = ng.attrs["ID"]
            ntx.attrs["version"] = "1"

            # exons: give them version=1 and update Parent; if exon_id/Name exist, keep or regenerate using ex_core
            for ex in getattr(ntx, "exons", []):
                if hasattr(ex, "attrs"):
                    ex.attrs["Parent"] = ntx.attrs["ID"]
                    if "version" in ex.attrs:
                        ex.attrs["version"] = "1"
                    # if no exon_id provided, you can generate one:
                    if "exon_id" not in ex.attrs and ex_core:
                        # allocate a new exon core via IDManager to keep width
                        new_ex_core = idm.new_id(ex_core).split(":",1)[-1]
                        new_ex_core = new_ex_core if "." not in new_ex_core else new_ex_core.split(".",1)[0]
                        ex.attrs["exon_id"] = new_ex_core
                        ex.attrs.setdefault("Name", new_ex_core)

            # CDS: keep or create protein IDs if needed (no version by Ensembl convention)
            for cd in getattr(ntx, "cdss", []):
                if hasattr(cd, "attrs"):
                    cd.attrs["Parent"] = ntx.attrs["ID"]
                    if "protein_id" not in cd.attrs and prot_core:
                        new_p_core = idm.new_id(prot_core).split(":",1)[-1]
                        new_p_core = new_p_core if "." not in new_p_core else new_p_core.split(".",1)[0]
                        cd.attrs["protein_id"] = new_p_core
                    if "ID" not in cd.attrs and prot_core:
                        cd.attrs["ID"] = f"{ns_cds}:{cd.attrs['protein_id']}"

            ng.transcripts[ntx.attrs["ID"]] = ntx

        updated[ng.attrs["ID"]] = ng
        idm.set(gid, ng.attrs["ID"], "new")

    # 4) Write outputs
    write_gff(updated, out_gff)

    # 5) Report
    mapped = sum(1 for _, st in idm.mapping.values() if st == "mapped")
    missing = sum(1 for _, st in idm.mapping.values() if st == "missing")
    new = sum(1 for _, st in idm.mapping.values() if st == "new")
    total = len(idm.mapping)
    with open(report_txt, "w") as fh:
        fh.write("ID Mapping Report\n")
        fh.write("=================\n\n")
        fh.write(f"Genes in reference: {len(ref)}\n")
        fh.write(f"Mapped: {mapped} ({mapped*100/max(1,total):.1f}%)\n")
        fh.write(f"Missing: {missing} ({missing*100/max(1,total):.1f}%)\n")
        fh.write(f"New: {new}\n\n")
        fh.write("Missing gene IDs:\n")
        for old, (_newid, status) in idm.mapping.items():
            if status == "missing":
                fh.write(f"  {old}\n")

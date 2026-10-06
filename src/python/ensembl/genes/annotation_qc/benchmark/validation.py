"""
Input and configuration validation.

Every annotation is scanned once (results cached by sha256) for: format,
compression vs file name, feature types, sequence names and maximum coordinates,
transcript-level types the comparator parser does not recognise, identity
(``ID``/``gene_id``/``transcript_id``) reuse, display-name reuse, verbatim
duplicated exon/CDS rows, GTF transcript_id reuse without transcript lines and
canonical tags. Genomes are checked for sequence names and lengths.

Validation reports problems; it never repairs biological structures. Known
comparator parser limitations are probed against the *current* code each time, so
their reported status reflects this checkout.
"""

from __future__ import annotations

import gzip
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from ensembl.genes.annotation_qc.benchmark.config import Experiment
from ensembl.genes.annotation_qc.benchmark.preparation import is_gzip
from ensembl.genes.annotation_qc.benchmark.provenance import now
from ensembl.genes.annotation_qc.benchmark.workspace import (
    Workspace,
    read_json,
    write_json,
)
from ensembl.genes.annotation_qc.parsers import annotation_normalise as normalise
from ensembl.genes.annotation_qc.parsers.annotation import (
    detect_annotation_format,
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.seqnames import build_seqname_mapping
from ensembl.genes.annotation_qc.parsers.sequence import parse_sequence_lengths

SCAN_VERSION = "2"
GTF_SPAN_WARNING = (
    2_000_000  # inferred GTF transcripts longer than this are flagged for review
)
_PREFIX = re.compile(r"^(?:gene|transcript|chromosome|mRNA):")
_GTF_KV = re.compile(r'(\S+)\s+"([^"]*)"')


def _attrs_gff(text: str) -> dict:
    out = {}
    for part in text.strip().strip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def scan_annotation(path: Path) -> dict:
    """Single streaming pass over an annotation; see module docstring."""
    path = Path(path)
    fmt = detect_annotation_format(str(path))
    compressed = is_gzip(path)
    feature_counts: Counter = Counter()
    seq_genes: Counter = Counter()
    seq_max_end: dict[str, int] = {}
    gene_ids: dict[tuple, int] = Counter()
    tx_ids: dict[tuple, int] = Counter()
    id_seqs: dict[str, set] = defaultdict(set)
    names: Counter = Counter()
    parented: list[tuple] = []  # (seq, type, parent) for GFF3 non-gene, non-child rows
    gff_gene_keys: set = set()
    child_hashes: set = set()
    duplicate_child_rows = 0
    duplicate_examples: list[str] = []
    canonical_rows = 0
    gtf_tx: dict = {}
    gtf_has_tx_rows = False
    nonchild_ids: set = set()
    child_parents: Counter = Counter()
    opener = gzip.open if compressed else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            seq, ftype = p[0], p[2]
            try:
                end = int(p[4])
            except ValueError:
                continue
            feature_counts[ftype] += 1
            if end > seq_max_end.get(seq, 0):
                seq_max_end[seq] = end
            if "Ensembl_canonical" in p[8] and ftype not in normalise.CHILD_FEATURES:
                canonical_rows += 1
            if ftype in normalise.CHILD_FEATURES:
                h = hash(line)
                if h in child_hashes:
                    duplicate_child_rows += 1
                    if len(duplicate_examples) < 5:
                        duplicate_examples.append(line.strip()[:160])
                child_hashes.add(h)
            if fmt == "gff3":
                if ftype in normalise.CHILD_FEATURES:
                    a = _attrs_gff(p[8])
                    for parent in a.get("Parent", "").split(","):
                        child_parents[(seq, _PREFIX.sub("", parent.strip()))] += 1
                    continue
                a = _attrs_gff(p[8])
                ident = _PREFIX.sub("", a.get("ID", ""))
                if ident:
                    nonchild_ids.add((seq, ident))
                if ftype in normalise.GENE_FEATURES:
                    gene_ids[(seq, ident)] += 1
                    gff_gene_keys.add((seq, ident))
                    id_seqs["gene:" + ident].add(seq)
                    seq_genes[seq] += 1
                    if a.get("Name"):
                        names[a["Name"]] += 1
                elif "Parent" in a:
                    parent = _PREFIX.sub("", a["Parent"].split(",")[0])
                    parented.append((seq, ftype, parent, ident))
            else:
                if ftype == "gene":
                    raw = p[8].strip()
                    gid = _GTF_KV.findall(raw)
                    ident = dict(gid).get("gene_id") if gid else raw.rstrip(";").strip()
                    gene_ids[(seq, ident)] += 1
                    id_seqs["gene:" + ident].add(seq)
                    seq_genes[seq] += 1
                    if dict(gid).get("gene_name"):
                        names[dict(gid)["gene_name"]] += 1
                elif ftype == "transcript":
                    gtf_has_tx_rows = True
                    raw = p[8].strip()
                    kv = dict(_GTF_KV.findall(raw))
                    ident = kv.get("transcript_id") or raw.rstrip(";").strip()
                    tx_ids[(seq, ident)] += 1
                    id_seqs["transcript:" + ident].add(seq)
                elif ftype in normalise.CHILD_FEATURES:
                    kv = dict(_GTF_KV.findall(p[8]))
                    tid = kv.get("transcript_id")
                    if tid:
                        rec = gtf_tx.setdefault(
                            (seq, tid),
                            {
                                "min": int(p[3]),
                                "max": end,
                                "strands": set(),
                                "genes": set(),
                            },
                        )
                        rec["min"] = min(rec["min"], int(p[3]))
                        rec["max"] = max(rec["max"], end)
                        rec["strands"].add(p[6])
                        if kv.get("gene_id"):
                            rec["genes"].add(kv["gene_id"])

    result = {
        "scan_version": SCAN_VERSION,
        "format": fmt,
        "compressed": compressed,
        "name_suffix": "".join(Path(path).suffixes[-2:]),
        "feature_counts": dict(feature_counts),
        "sequences": {
            s: {"max_end": seq_max_end[s], "genes": seq_genes.get(s, 0)}
            for s in sorted(seq_max_end)
        },
        "genes": sum(seq_genes.values()),
        "identity_reuse_same_sequence": {
            "gene": sorted(f"{s}:{i}" for (s, i), n in gene_ids.items() if n > 1)[:20],
            "transcript": sorted(f"{s}:{i}" for (s, i), n in tx_ids.items() if n > 1)[
                :20
            ],
        },
        "identity_reuse_across_sequences": sum(
            1 for seqs in id_seqs.values() if len(seqs) > 1
        ),
        "display_name_reuse": sum(1 for n in names.values() if n > 1),
        "duplicate_child_rows": duplicate_child_rows,
        "duplicate_child_examples": duplicate_examples,
        "canonical_tagged_records": canonical_rows,
    }
    if fmt == "gff3":
        transcript_like = Counter(
            t for s, t, parent, _ in parented if (s, parent) in gff_gene_keys
        )
        tx_counter = Counter(
            (s, i) for s, t, parent, i in parented if (s, parent) in gff_gene_keys
        )
        result["transcript_types"] = dict(transcript_like)
        result["unrecognised_transcript_types"] = {
            t: n
            for t, n in transcript_like.items()
            if t not in normalise.TRANSCRIPT_FEATURES
        }
        result["identity_reuse_same_sequence"]["transcript"] = sorted(
            f"{s}:{i}" for (s, i), n in tx_counter.items() if n > 1
        )[:20]
        for s, t, parent, i in parented:
            if (s, parent) in gff_gene_keys:
                id_seqs["transcript:" + i].add(s)
        result["identity_reuse_across_sequences"] = sum(
            1 for seqs in id_seqs.values() if len(seqs) > 1
        )
        orphans = {k: n for k, n in child_parents.items() if k not in nonchild_ids}
        result["orphan_child_rows"] = sum(orphans.values())
        result["orphan_child_examples"] = [
            f"{s}:{p or '(no Parent)'}" for s, p in list(orphans)[:5]
        ]
    else:
        result["gtf_has_transcript_rows"] = gtf_has_tx_rows
        if not gtf_has_tx_rows:
            conflicts = [
                f"{s}:{t}"
                for (s, t), r in gtf_tx.items()
                if len(r["strands"]) > 1 or len(r["genes"]) > 1
            ]
            long = [
                f"{s}:{t} ({r['max'] - r['min'] + 1:,} bp)"
                for (s, t), r in gtf_tx.items()
                if r["max"] - r["min"] + 1 > GTF_SPAN_WARNING
            ]
            result["gtf_inferred_transcript_conflicts"] = conflicts[:20]
            result["gtf_inferred_transcript_conflict_count"] = len(conflicts)
            result["gtf_inferred_transcripts_over_span_limit"] = long[:20]
            result["gtf_inferred_transcripts_over_span_limit_count"] = len(long)
    return result


def cached_scan(ws: Workspace, path: Path) -> dict:
    sha = ws.checksums.sha256(path)
    cache = ws.root / "cache" / f"scan_{sha}.json"
    data = read_json(cache)
    if data and data.get("scan_version") == SCAN_VERSION:
        return data
    data = scan_annotation(path)
    write_json(cache, data)
    return data


def cached_fasta_lengths(ws: Workspace, path: Path) -> dict:
    sha = ws.checksums.sha256(path)
    cache = ws.root / "cache" / f"fasta_{sha}.json"
    data = read_json(cache)
    if data is None:
        data = parse_sequence_lengths(str(path))
        write_json(cache, data)
    return data


# ---------------------------------------------------------------------------
# Known parser limitations, probed against the current comparator code.


def _probe(text: str, suffix: str, check) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / f"probe{suffix}"
        if suffix.endswith(".bgz"):
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        try:
            return check(parse_annotation_for_comparison(str(path)))
        except Exception as error:  # noqa: BLE001 - any failure is the observation
            return f"present ({type(error).__name__})"


def parser_limitations() -> list[dict]:
    gff = "##gff-version 3\n1\tt\tgene\t1\t90\t.\t+\t.\tID=g1\n1\tt\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n1\tt\texon\t1\t90\t.\t+\t.\tParent=t1\n1\tt\tCDS\t1\t90\t.\t+\t0\tParent=t1\n"
    unlisted = gff.replace("\tmRNA\t", "\tunconfirmed_transcript\t")
    dup = gff + "1\tt\tCDS\t1\t90\t.\t+\t0\tParent=t1\n"
    gtf = (
        '1\tt\tCDS\t1\t90\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        '1\tt\tCDS\t500001\t500090\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    )
    return [
        {
            "id": "bgz_suffix",
            "title": "gzip/BGZF annotation not named *.gz",
            "status": _probe(gff, ".gff3.bgz", lambda df: "resolved"),
            "handling": "use the 'decompress' preparation step or rename to *.gz; validation flags affected inputs",
        },
        {
            "id": "unlisted_transcript_types",
            "title": "Ensembl unconfirmed_transcript / gene_segment transcript types",
            "status": _probe(unlisted, ".gff3", lambda df: "resolved"),
            "handling": "use the 'retype_transcript_types' preparation step (audited, IDs/Parents unchanged)",
        },
        {
            "id": "duplicate_child_rows",
            "title": "verbatim duplicated exon/CDS rows counted twice",
            "status": _probe(
                dup,
                ".gff3",
                lambda df: (
                    "present (duplicates kept)"
                    if int((df["Feature"] == "CDS").sum()) > 1
                    else "resolved"
                ),
            ),
            "handling": "not repaired; validation reports affected files as errors",
        },
        {
            "id": "gtf_inferred_transcript_fusion",
            "title": "GTF without transcript lines: transcript_id reused at two loci is fused",
            "status": _probe(
                gtf,
                ".gtf",
                lambda df: (
                    "present (fused silently)"
                    if int((df["Feature"] == "transcript").sum()) == 1
                    else "resolved"
                ),
            ),
            "handling": "not repaired; validation flags strand/gene conflicts and very long inferred transcripts",
        },
    ]


# ---------------------------------------------------------------------------


def _issue(level: str, where: str, message: str) -> dict:
    return {"level": level, "where": where, "message": message}


def check_annotation(
    scan: dict, path: Path, where: str, prepared: bool, issues: list, will_prepare: bool
) -> None:
    if (
        scan["compressed"]
        and not str(path).lower().endswith(".gz")
        and not will_prepare
    ):
        issues.append(
            _issue(
                "error",
                where,
                f"{path.name} is gzip-compressed but not named *.gz; the comparator cannot read it "
                "(known limitation). Add a 'decompress' preparation step.",
            )
        )
    unlisted = scan.get("unrecognised_transcript_types") or {}
    if unlisted:
        level = "info" if prepared else "error"
        issues.append(
            _issue(
                level,
                where,
                f"transcript types not recognised by the parser: {unlisted}"
                + (
                    " (handled by retype_transcript_types)"
                    if prepared
                    else " — add a retype_transcript_types step"
                ),
            )
        )
    reuse = scan["identity_reuse_same_sequence"]
    if reuse["gene"] or reuse["transcript"]:
        issues.append(
            _issue(
                "error",
                where,
                f"identity reused on the same sequence (comparator rejects): {reuse}",
            )
        )
    if scan["identity_reuse_across_sequences"]:
        issues.append(
            _issue(
                "info",
                where,
                f"{scan['identity_reuse_across_sequences']} identifiers reused on different "
                "sequences; the comparator namespaces them as seqname:id",
            )
        )
    if scan["duplicate_child_rows"]:
        issues.append(
            _issue(
                "error",
                where,
                f"{scan['duplicate_child_rows']} verbatim duplicated exon/CDS rows; the comparator "
                f"counts them twice (known limitation). Examples: {scan['duplicate_child_examples'][:2]}",
            )
        )
    if scan.get("orphan_child_rows"):
        issues.append(
            _issue(
                "error",
                where,
                f"{scan['orphan_child_rows']} exon/CDS rows whose Parent is not a record on the same "
                f"sequence (the comparator stops): {scan['orphan_child_examples']}",
            )
        )
    if scan.get("gtf_inferred_transcript_conflict_count"):
        issues.append(
            _issue(
                "error",
                where,
                f"{scan['gtf_inferred_transcript_conflict_count']} transcript_ids without transcript "
                "lines have children on both strands or in several genes; they would be fused",
            )
        )
    if scan.get("gtf_inferred_transcripts_over_span_limit_count"):
        issues.append(
            _issue(
                "warning",
                where,
                f"{scan['gtf_inferred_transcripts_over_span_limit_count']} inferred GTF transcripts span "
                f"> {GTF_SPAN_WARNING:,} bp; check for transcript_id reuse at distinct loci",
            )
        )


def validate_experiment(
    experiment: Experiment,
    ws: Workspace | None = None,
    prepared_records: dict | None = None,
) -> dict:
    """
    Validate configuration inputs. ``prepared_records`` (reference_id -> record from
    preparation.prepare_reference) lets checks run on prepared files; without them the
    source files are checked.
    """
    ws = ws or Workspace(experiment)
    prepared_records = prepared_records or {}
    issues: list[dict] = []
    references = {}
    for ref in experiment.references.values():
        where = f"reference {ref.id}"
        entry = {"id": ref.id, "inputs": {}}
        for label, path in (
            ("annotation", ref.annotation),
            ("genome", ref.genome),
            ("seqname_map", ref.seqname_map),
        ):
            if path is None:
                continue
            if not Path(path).exists():
                issues.append(_issue("error", where, f"{label} not found: {path}"))
            else:
                entry["inputs"][label] = ws.checksums.record(path)
        if ref.genome is None:
            issues.append(
                _issue(
                    "warning",
                    where,
                    "no genome FASTA: sequence names and bounds are not checked against an assembly",
                )
            )
        if "annotation" not in entry["inputs"]:
            references[ref.id] = entry
            continue
        prep = prepared_records.get(ref.id)
        ann_path = (
            Path(prep["annotation"]["output"]["path"]) if prep else Path(ref.annotation)
        )
        will_prepare = bool(ref.annotation_preparation.steps)
        scan = cached_scan(ws, ann_path)
        entry["annotation_scan"] = scan
        entry["checked_file"] = str(ann_path)
        check_annotation(
            scan,
            ann_path,
            where,
            prepared=prep is not None and will_prepare,
            issues=issues,
            will_prepare=will_prepare,
        )
        if prep is None and will_prepare:
            src_scan = scan
            if src_scan.get("unrecognised_transcript_types"):
                types = {
                    t
                    for s in ref.annotation_preparation.steps
                    if s["name"] == "retype_transcript_types"
                    for t in s["types"]
                }
                left = set(src_scan["unrecognised_transcript_types"]) - types
                if left:
                    issues.append(
                        _issue(
                            "error",
                            where,
                            f"retype step does not cover: {sorted(left)}",
                        )
                    )
        entry["canonical_available"] = scan["canonical_tagged_records"] > 0
        lengths = None
        if ref.genome and Path(ref.genome).exists():
            genome_path = (
                Path(prep["genome"]["output"]["path"])
                if prep and prep.get("genome")
                else Path(ref.genome)
            )
            lengths = cached_fasta_lengths(ws, genome_path)
            entry["genome_sequences"] = len(lengths)
        mapping = (
            build_seqname_mapping(str(ref.seqname_map), None)
            if ref.seqname_map and Path(ref.seqname_map).exists()
            else {}
        )
        entry["sequence_check"] = _sequence_check(scan, lengths, where, issues, mapping)
        scope = (prep or {}).get("scope") or {}
        scope_seqs = scope.get("sequences")
        entry["scope"] = {
            "type": ref.scope.get("type"),
            "sequences": len(scope_seqs) if scope_seqs else None,
        }
        if scope_seqs and lengths:
            missing = [s for s in scope_seqs if s not in lengths]
            if missing:
                issues.append(
                    _issue(
                        "error",
                        where,
                        f"{len(missing)} scope sequences not in the genome FASTA: {missing[:5]}",
                    )
                )
        references[ref.id] = entry

    runs = {}
    for run in experiment.runs.values():
        where = f"run {run.run_id}"
        entry = {"run_id": run.run_id, "reference": run.reference}
        if not Path(run.annotation).exists():
            issues.append(
                _issue("error", where, f"annotation not found: {run.annotation}")
            )
            runs[run.run_id] = entry
            continue
        entry["input"] = ws.checksums.record(run.annotation)
        scan = cached_scan(ws, Path(run.annotation))
        entry["annotation_scan"] = scan
        check_annotation(
            scan,
            Path(run.annotation),
            where,
            prepared=False,
            issues=issues,
            will_prepare=bool(run.preparation.steps),
        )
        if run.format != "auto" and run.format != scan["format"]:
            issues.append(
                _issue(
                    "warning",
                    where,
                    f"configured format {run.format} but content looks like {scan['format']}",
                )
            )
        ref_entry = references.get(run.reference, {})
        ref_prep = prepared_records.get(run.reference)
        lengths = None
        ref = experiment.references[run.reference]
        if ref.genome and Path(ref.genome).exists():
            genome_path = (
                Path(ref_prep["genome"]["output"]["path"])
                if ref_prep and ref_prep.get("genome")
                else Path(ref.genome)
            )
            lengths = cached_fasta_lengths(ws, genome_path)
        mapping = (
            build_seqname_mapping(str(ref.seqname_map), None)
            if ref.seqname_map and Path(ref.seqname_map).exists()
            else {}
        )
        entry["sequence_check"] = _sequence_check(scan, lengths, where, issues, mapping)
        scope_seqs = ((ref_prep or {}).get("scope") or {}).get("sequences")
        if scope_seqs is not None:
            scope_set = set(scope_seqs)
            outside = {
                s: v["genes"]
                for s, v in scan["sequences"].items()
                if mapping.get(s, s) not in scope_set
            }
            entry["outside_scope"] = {
                "sequences": len(outside),
                "genes": sum(outside.values()),
            }
            ref_scope = {s for s in scope_set}
            mapped = {mapping.get(s, s) for s in scan["sequences"]}
            no_preds = [s for s in ref_scope if s not in mapped]
            entry["scope_sequences_without_predictions"] = len(no_preds)
        if not run.provenance:
            issues.append(
                _issue("warning", where, "no provenance recorded; shown as unknown")
            )
        for field_name in ("tool_version", "model"):
            if not getattr(run, field_name):
                entry.setdefault("unknown", []).append(field_name)
        if ref_entry and not ref_entry.get("canonical_available", True):
            entry["policy_notes"] = (
                "canonical policy unavailable for this reference (no Ensembl_canonical tags)"
            )
        runs[run.run_id] = entry

    report = {
        "experiment": experiment.id,
        "validated": now(),
        "ok": not any(i["level"] == "error" for i in issues),
        "issues": issues,
        "references": references,
        "runs": runs,
        "parser_limitations": parser_limitations(),
    }
    write_json(ws.validation_file(), report)
    return report


def _sequence_check(
    scan: dict,
    lengths: dict | None,
    where: str,
    issues: list,
    mapping: dict | None = None,
) -> dict:
    if lengths is None:
        return {"checked": False}
    mapping = mapping or {}
    seqs = {
        mapping.get(s, s): v for s, v in scan["sequences"].items()
    }  # comparison names, as the comparator applies them
    missing = [s for s in seqs if s not in lengths]
    beyond = [s for s, v in seqs.items() if s in lengths and v["max_end"] > lengths[s]]
    if missing:
        issues.append(
            _issue(
                "error",
                where,
                f"{len(missing)} sequences not in the genome FASTA: {missing[:5]} "
                "(check the assembly or supply a seqname_map)",
            )
        )
    if beyond:
        issues.append(
            _issue(
                "error",
                where,
                f"features beyond sequence ends on {beyond[:5]}; possibly a different assembly",
            )
        )
    return {
        "checked": True,
        "sequences": len(scan["sequences"]),
        "missing_from_fasta": missing[:20],
        "missing_count": len(missing),
        "out_of_bounds": beyond[:20],
        "out_of_bounds_count": len(beyond),
    }

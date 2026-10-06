"""
Build the dashboard dataset (SQLite) from prepared inputs and completed results.

Gene models are parsed once with the comparator's own parser and selection
functions (so evaluated transcripts and CDS segment counts are exactly the ones the
comparator used) and stored as compact per-transcript rows indexed by sequence and
start. The dashboard then answers bounded locus queries without re-reading
annotations.

The build is incremental: each component (reference models, run models, one result)
has a key derived from its inputs; unchanged components are skipped.

Identity: reference and query genes keep the comparator's sequence-scoped
``gene_id`` (``seqname:id`` when an ID is reused on several sequences) and the file
value in ``original_gene_id``. Display names are stored only as labels.
"""

from __future__ import annotations

import gzip
import json
import re
import sqlite3
import time
from pathlib import Path

import pandas as pd

from ensembl.genes.annotation_qc.benchmark.config import Experiment
from ensembl.genes.annotation_qc.benchmark.derive import (
    derive_metrics,
    load_independent,
)
from ensembl.genes.annotation_qc.benchmark.preparation import is_gzip
from ensembl.genes.annotation_qc.benchmark.provenance import now, stable_hash
from ensembl.genes.annotation_qc.benchmark.runs import STATUS_COMPLETE, plan
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace, read_json
from ensembl.genes.annotation_qc.metrics.pairwise.classify import build_gene_models
from ensembl.genes.annotation_qc.metrics.pairwise.selection import (
    apply_evaluation_mode,
    filter_by_biotype,
    select_transcripts,
    subset_to_regions,
)
from ensembl.genes.annotation_qc.parsers.annotation import (
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.regions import Region
from ensembl.genes.annotation_qc.parsers.seqnames import (
    apply_seqname_mapping,
    build_seqname_mapping,
)

STORE_VERSION = "1"
METRIC_DICTIONARY = Path(__file__).with_name("metric_dictionary.json")
PALETTE = [
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#4a3aa7",
    "#e34948",
]

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS components (component TEXT PRIMARY KEY, key TEXT, built TEXT, info TEXT);
CREATE TABLE IF NOT EXISTS ref_genes (reference_id TEXT, gene_id TEXT, original_gene_id TEXT, name TEXT, chrom TEXT,
    start INTEGER, end INTEGER, strand TEXT, biotype TEXT, source_feature TEXT, in_scope INTEGER, evaluated INTEGER, search TEXT);
CREATE TABLE IF NOT EXISTS ref_gene_policy (reference_id TEXT, policy_id TEXT, gene_id TEXT, evaluated_transcripts TEXT,
    multi_segment INTEGER, canonical_fallback INTEGER);
CREATE TABLE IF NOT EXISTS transcripts (source TEXT, gene_id TEXT, transcript_id TEXT, chrom TEXT, start INTEGER, end INTEGER,
    strand TEXT, biotype TEXT, tags TEXT, eligible INTEGER, exons TEXT, cds TEXT);
CREATE TABLE IF NOT EXISTS query_genes (run_id TEXT, gene_id TEXT, original_gene_id TEXT, name TEXT, chrom TEXT, start INTEGER,
    end INTEGER, strand TEXT, in_scope INTEGER, search TEXT);
CREATE TABLE IF NOT EXISTS seq_maxlen (source TEXT, chrom TEXT, max_len INTEGER);
CREATE TABLE IF NOT EXISTS results (reference_id TEXT, run_id TEXT, policy_id TEXT, status TEXT, reason TEXT, source TEXT, record TEXT);
CREATE TABLE IF NOT EXISTS metrics (reference_id TEXT, run_id TEXT, policy_id TEXT, metric_id TEXT, numerator REAL,
    denominator REAL, value REAL, status TEXT, source TEXT, note TEXT);
CREATE TABLE IF NOT EXISTS ref_outcomes (reference_id TEXT, policy_id TEXT, run_id TEXT, gene_id TEXT, classification TEXT,
    classification_cds TEXT, cds_status TEXT, cds_exact INTEGER, cds_overlap REAL, exon_overlap REAL, cds_intron_chain_match TEXT,
    intron_chain_match TEXT, exon_coordinate_exact INTEGER, cds_matched_id TEXT, best_cds_ref_tx TEXT, best_cds_query_tx TEXT,
    matched_id TEXT, best_match_ref_tx TEXT, best_match_query_tx TEXT, counterpart_count INTEGER, counterpart_ids TEXT,
    strand_mismatch_basis TEXT, multi_segment INTEGER, any_pair_chain INTEGER);
CREATE TABLE IF NOT EXISTS query_outcomes (reference_id TEXT, policy_id TEXT, run_id TEXT, gene_id TEXT, original_gene_id TEXT,
    chrom TEXT, start INTEGER, end INTEGER, strand TEXT, classification TEXT, classification_cds TEXT, cds_exact INTEGER,
    cds_overlap REAL, cds_matched_id TEXT, best_cds_ref_tx TEXT, best_cds_query_tx TEXT, novel_category TEXT,
    counterpart_count INTEGER, counterpart_ids TEXT, strand_mismatch_basis TEXT);
CREATE INDEX IF NOT EXISTS ix_ref_genes ON ref_genes (reference_id, chrom, start);
CREATE INDEX IF NOT EXISTS ix_ref_genes_id ON ref_genes (reference_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_rgp ON ref_gene_policy (reference_id, policy_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_tx ON transcripts (source, chrom, start);
CREATE INDEX IF NOT EXISTS ix_tx_gene ON transcripts (source, gene_id);
CREATE INDEX IF NOT EXISTS ix_tx_id ON transcripts (source, transcript_id);
CREATE INDEX IF NOT EXISTS ix_qg ON query_genes (run_id, chrom, start);
CREATE INDEX IF NOT EXISTS ix_qg_id ON query_genes (run_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_metrics ON metrics (reference_id, policy_id, run_id);
CREATE INDEX IF NOT EXISTS ix_ro ON ref_outcomes (reference_id, policy_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_ro_run ON ref_outcomes (reference_id, policy_id, run_id, cds_status);
CREATE INDEX IF NOT EXISTS ix_qo ON query_outcomes (reference_id, policy_id, run_id, chrom, start);
CREATE INDEX IF NOT EXISTS ix_qo_id ON query_outcomes (run_id, policy_id, gene_id);
"""

_PREFIX = re.compile(r"^(?:gene|transcript|chromosome|mRNA):")
_GTF_KV = re.compile(r'(\S+)\s+"([^"]*)"')


def connect(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("PRAGMA journal_mode=WAL")
    con.executescript(SCHEMA)
    return con


def _component_key(con, component) -> str | None:
    row = con.execute(
        "SELECT key FROM components WHERE component=?", (component,)
    ).fetchone()
    return row[0] if row else None


def _set_component(con, component, key, info=None):
    con.execute(
        "INSERT OR REPLACE INTO components VALUES (?,?,?,?)",
        (component, key, now(), json.dumps(info or {})),
    )


def gene_names(path: Path) -> dict[tuple[str, str], str]:
    """(seqname, file identity) -> display name from GFF3 Name / GTF gene_name. Labels only."""
    names = {}
    opener = gzip.open if is_gzip(path) else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if len(p) < 9 or p[2] not in ("gene", "ncRNA_gene", "pseudogene"):
                continue
            attrs = p[8]
            if "Name=" in attrs or "ID=" in attrs:
                kv = dict(
                    x.split("=", 1)
                    for x in attrs.strip().strip(";").split(";")
                    if "=" in x
                )
                ident = kv.get("gene_id") or _PREFIX.sub("", kv.get("ID", ""))
                if kv.get("Name"):
                    names[(p[0], ident)] = kv["Name"]
            else:
                kv = dict(_GTF_KV.findall(attrs))
                if kv.get("gene_name"):
                    names[(p[0], kv.get("gene_id", ""))] = kv["gene_name"]
    return names


def _intervals(df: pd.DataFrame) -> pd.DataFrame:
    """transcript_id -> 'start-end,...' strings (one-based inclusive) for exons and CDS."""
    feat = df.loc[
        df.Feature.isin(["exon", "CDS"]), ["transcript_id", "Feature", "Start", "End"]
    ]
    feat = feat.sort_values(["transcript_id", "Feature", "Start", "End"], kind="stable")
    iv = (feat.Start + 1).astype(str) + "-" + feat.End.astype(str)
    agg = (
        iv.groupby([feat.transcript_id, feat.Feature], sort=False)
        .agg(",".join)
        .unstack("Feature")
    )
    for col in ("exon", "CDS"):
        if col not in agg.columns:
            agg[col] = ""
    return agg.fillna("")


def _transcript_rows(
    source: str, df: pd.DataFrame, eligible: set | None
) -> list[tuple]:
    tx = df[df.Feature == "transcript"].drop_duplicates("transcript_id")
    ivs = _intervals(df)
    rows = []
    for r in tx.itertuples(index=False):
        exons = ivs.at[r.transcript_id, "exon"] if r.transcript_id in ivs.index else ""
        cds = ivs.at[r.transcript_id, "CDS"] if r.transcript_id in ivs.index else ""
        rows.append(
            (
                source,
                r.gene_id,
                r.transcript_id,
                r.Chromosome,
                int(r.Start) + 1,
                int(r.End),
                r.Strand,
                r.transcript_biotype,
                r.tags,
                int(eligible is None or r.transcript_id in eligible),
                exons,
                cds,
            )
        )
    return rows


def _scope_regions(prep: dict) -> list[Region] | None:
    seqs = (prep.get("scope") or {}).get("sequences")
    return [Region(s) for s in seqs] if seqs else None


def build_reference(
    con, experiment: Experiment, ws: Workspace, ref, prep: dict, log=print
) -> None:
    ev = experiment.evaluation
    key = stable_hash(
        {
            "v": STORE_VERSION,
            "annotation": prep["annotation"]["output"]["sha256"],
            "scope": prep["scope"]["sha256"],
            "seqmap": (prep.get("seqname_map") or {}).get("sha256"),
            "evaluation": ev,
            "policies": {
                p.id: p.reference_transcript_selection
                for p in experiment.policies.values()
            },
            "format": ref.annotation_format,
        }
    )
    component = f"reference:{ref.id}"
    if _component_key(con, component) == key:
        log(f"  reference {ref.id}: models up to date")
        return
    t0 = time.time()
    path = Path(prep["annotation"]["output"]["path"])
    log(f"  reference {ref.id}: parsing {path.name} with the comparator parser (once)")
    df = parse_annotation_for_comparison(str(path), ref.annotation_format)
    mapping = build_seqname_mapping(
        str(ref.seqname_map) if ref.seqname_map else None, None
    )
    raw_names = gene_names(path)
    names = {(mapping.get(c, c), i): n for (c, i), n in raw_names.items()}
    df, _ = apply_seqname_mapping(df, mapping)
    regions = _scope_regions(prep)
    scope = {r.seqname for r in regions} if regions else None

    base = apply_evaluation_mode(df, ev["evaluation_mode"])
    if ev["reference_gene_biotypes"]:
        base = filter_by_biotype(base, ev["reference_gene_biotypes"], None)
    eligible_df = (
        filter_by_biotype(base, None, ev["reference_transcript_biotypes"])
        if ev["reference_transcript_biotypes"]
        else base
    )
    eligible = set(
        eligible_df.loc[eligible_df.Feature == "transcript", "transcript_id"]
    )
    evaluated_genes = set()
    policy_rows = []
    for policy in experiment.policies.values():
        selected = select_transcripts(
            eligible_df, policy.reference_transcript_selection
        )
        if regions:
            selected = subset_to_regions(selected, regions)
        models = build_gene_models(selected)
        canonical = set()
        if policy.reference_transcript_selection == "canonical":
            tx = eligible_df[eligible_df.Feature == "transcript"]
            canonical = set(
                tx.loc[
                    tx.tags.str.contains("Ensembl_canonical", regex=False), "gene_id"
                ]
            )
        for g in models.genes.itertuples(index=False):
            evaluated_genes.add(g.gene_id)
            fallback = int(
                policy.reference_transcript_selection == "canonical"
                and g.gene_id not in canonical
            )
            policy_rows.append(
                (
                    ref.id,
                    policy.id,
                    g.gene_id,
                    ",".join(g.transcript_ids),
                    int(g.max_cds_count > 1),
                    fallback,
                )
            )

    genes = df[df.Feature == "gene"]
    tx_by_gene = (
        base[base.Feature == "transcript"]
        .groupby("gene_id")["transcript_id"]
        .agg(" ".join)
    )
    gene_rows = []
    for g in genes.itertuples(index=False):
        name = names.get((g.Chromosome, g.original_gene_id), "")
        in_scope = int(scope is None or g.Chromosome in scope)
        search = " ".join(
            filter(
                None,
                [g.gene_id, g.original_gene_id, name, tx_by_gene.get(g.gene_id, "")],
            )
        ).lower()
        gene_rows.append(
            (
                ref.id,
                g.gene_id,
                g.original_gene_id,
                name,
                g.Chromosome,
                int(g.Start) + 1,
                int(g.End),
                g.Strand,
                g.gene_biotype,
                g.source_feature,
                in_scope,
                int(g.gene_id in evaluated_genes),
                search,
            )
        )
    tx_rows = _transcript_rows(f"ref:{ref.id}", base, eligible)
    del df
    source = f"ref:{ref.id}"
    with con:
        con.execute("DELETE FROM ref_genes WHERE reference_id=?", (ref.id,))
        con.execute("DELETE FROM ref_gene_policy WHERE reference_id=?", (ref.id,))
        con.execute("DELETE FROM transcripts WHERE source=?", (source,))
        con.execute(
            "DELETE FROM seq_maxlen WHERE source IN (?,?)",
            (source, f"refgene:{ref.id}"),
        )
        con.executemany(
            "INSERT INTO ref_genes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", gene_rows
        )
        con.executemany("INSERT INTO ref_gene_policy VALUES (?,?,?,?,?,?)", policy_rows)
        con.executemany(
            "INSERT INTO transcripts VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", tx_rows
        )
        con.execute(
            "INSERT INTO seq_maxlen SELECT source, chrom, MAX(end-start+1) FROM transcripts WHERE source=? GROUP BY chrom",
            (source,),
        )
        con.execute(
            "INSERT INTO seq_maxlen SELECT ?, chrom, MAX(end-start+1) FROM ref_genes WHERE reference_id=? GROUP BY chrom",
            (f"refgene:{ref.id}", ref.id),
        )
        _set_component(
            con,
            component,
            key,
            {
                "genes": len(gene_rows),
                "transcripts": len(tx_rows),
                "seconds": round(time.time() - t0, 1),
            },
        )
    log(
        f"  reference {ref.id}: {len(gene_rows):,} genes, {len(tx_rows):,} transcripts stored in {time.time() - t0:.0f} s"
    )


def build_run_models(
    con, experiment: Experiment, ws: Workspace, run, prep: dict, log=print
) -> None:
    query = ws.run_preparation_record(run.run_id) or {
        "output": ws.checksums.record(run.annotation)
    }
    path = Path(query["output"]["path"])
    ref = experiment.references[run.reference]
    key = stable_hash(
        {
            "v": STORE_VERSION,
            "query": query["output"]["sha256"],
            "scope": prep["scope"]["sha256"],
            "seqmap": (prep.get("seqname_map") or {}).get("sha256"),
            "format": run.format,
        }
    )
    component = f"run:{run.run_id}"
    if _component_key(con, component) == key:
        return
    df = parse_annotation_for_comparison(str(path), run.format)
    mapping = build_seqname_mapping(
        str(ref.seqname_map) if ref.seqname_map else None, None
    )
    names = {(mapping.get(c, c), i): n for (c, i), n in gene_names(path).items()}
    df, _ = apply_seqname_mapping(df, mapping)
    seqs = (prep.get("scope") or {}).get("sequences")
    scope = set(seqs) if seqs else None
    genes = df[df.Feature == "gene"]
    tx_by_gene = (
        df[df.Feature == "transcript"].groupby("gene_id")["transcript_id"].agg(" ".join)
    )
    rows = []
    for g in genes.itertuples(index=False):
        name = names.get((g.Chromosome, g.original_gene_id), "")
        search = " ".join(
            filter(
                None,
                [g.gene_id, g.original_gene_id, name, tx_by_gene.get(g.gene_id, "")],
            )
        ).lower()
        rows.append(
            (
                run.run_id,
                g.gene_id,
                g.original_gene_id,
                name,
                g.Chromosome,
                int(g.Start) + 1,
                int(g.End),
                g.Strand,
                int(scope is None or g.Chromosome in scope),
                search,
            )
        )
    source = f"run:{run.run_id}"
    tx_rows = _transcript_rows(source, df, None)
    with con:
        con.execute("DELETE FROM query_genes WHERE run_id=?", (run.run_id,))
        con.execute("DELETE FROM transcripts WHERE source=?", (source,))
        con.execute("DELETE FROM seq_maxlen WHERE source=?", (source,))
        con.executemany("INSERT INTO query_genes VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
        con.executemany(
            "INSERT INTO transcripts VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", tx_rows
        )
        con.execute(
            "INSERT INTO seq_maxlen SELECT source, chrom, MAX(end-start+1) FROM transcripts WHERE source=? GROUP BY chrom",
            (source,),
        )
        _set_component(
            con, component, key, {"genes": len(rows), "transcripts": len(tx_rows)}
        )
    log(f"  run {run.run_id}: {len(rows):,} genes stored")


def _independent_table(
    ws: Workspace, run_id: str, policy_id: str
) -> tuple[Path | None, str]:
    computed = ws.independent_dir(run_id) / f"{policy_id}.reference_genes.tsv"
    if computed.exists():
        return computed, "computed by `annotation-qc benchmark independent`"
    imported = read_json(ws.independent_dir(run_id) / f"{policy_id}.import.json")
    if imported:
        return Path(imported["source_path"]), "imported from the benchmark directory"
    return None, "independent analysis not run for this result"


def build_result(con, experiment, ws, task, reference_ready: bool, log=print) -> None:
    run, policy, ref = task.run, task.policy, task.reference
    record = (
        ws.run_record(run.run_id, policy.id) if task.status == STATUS_COMPLETE else None
    )
    ind_path, ind_origin = _independent_table(ws, run.run_id, policy.id)
    ind_sha = ws.checksums.sha256(ind_path) if ind_path and ind_path.exists() else None
    ref_key = _component_key(con, f"reference:{ref.id}") if reference_ready else None
    key = stable_hash(
        {
            "v": STORE_VERSION,
            "status": task.status,
            "reason": task.reason,
            "cache_key": task.key,
            "record_key": (record or {}).get("cache_key"),
            "independent": ind_sha,
            "ref": ref_key,
        }
    )
    component = f"result:{run.run_id}:{policy.id}"
    if _component_key(con, component) == key:
        return
    with con:
        for table in ("results", "metrics", "ref_outcomes", "query_outcomes"):
            con.execute(
                f"DELETE FROM {table} WHERE run_id=? AND policy_id=?",
                (run.run_id, policy.id),
            )
        public = {k: v for k, v in (record or {}).items() if k not in ("environment",)}
        attempts = ws.attempts(run.run_id, policy.id)
        if attempts:
            public["attempts"] = [
                {
                    k: a.get(k)
                    for k in (
                        "status",
                        "started",
                        "finished",
                        "exit_code",
                        "source",
                        "partial_dir",
                    )
                }
                for a in attempts[-10:]
            ]
        public["independent"] = {
            "origin": ind_origin,
            "path": str(ind_path) if ind_path else None,
        }
        con.execute(
            "INSERT INTO results VALUES (?,?,?,?,?,?,?)",
            (
                ref.id,
                run.run_id,
                policy.id,
                task.status,
                task.reason,
                (record or {}).get("source"),
                json.dumps(public, default=str),
            ),
        )
        if record is None:
            _set_component(con, component, key)
            return
        out_dir = Path(record["output_dir"])
        multi = None
        if reference_ready:
            multi = {
                g: bool(m)
                for g, m in con.execute(
                    "SELECT gene_id, multi_segment FROM ref_gene_policy WHERE reference_id=? AND policy_id=?",
                    (ref.id, policy.id),
                )
            }
        independent, ind_check = None, {"status": "absent"}
        det_ref = None
        rows, det_ref, det_q = derive_metrics(out_dir, multi)
        if ind_path is not None and det_ref is not None:
            independent, ind_check = load_independent(ind_path, det_ref)
            if independent is not None:
                rows, det_ref, det_q = derive_metrics(out_dir, multi, independent)
            else:
                rows, det_ref, det_q = derive_metrics(
                    out_dir,
                    multi,
                    None,
                    independent_note=f"independent table {ind_check.get('status')}: {ind_check.get('reason', '')}",
                )
        public["independent"].update(ind_check)
        con.execute(
            "UPDATE results SET record=? WHERE run_id=? AND policy_id=?",
            (json.dumps(public, default=str), run.run_id, policy.id),
        )
        con.executemany(
            "INSERT INTO metrics VALUES (?,?,?,?,?,?,?,?,?,?)",
            [
                (
                    ref.id,
                    run.run_id,
                    policy.id,
                    r["metric_id"],
                    r["numerator"],
                    r["denominator"],
                    r["value"],
                    r["status"],
                    r["source"],
                    r["note"],
                )
                for r in rows
            ],
        )
        if det_ref is not None:
            for col in ("multi_segment", "any_pair_chain"):
                if col not in det_ref:
                    det_ref[col] = None
            cols = [
                "gene_id",
                "classification",
                "classification_cds",
                "cds_status",
                "cds_exact",
                "cds_overlap",
                "exon_overlap",
                "cds_intron_chain_match",
                "intron_chain_match",
                "exon_coordinate_exact",
                "cds_matched_id",
                "best_cds_ref_tx",
                "best_cds_query_tx",
                "matched_id",
                "best_match_ref_tx",
                "best_match_query_tx",
                "counterpart_count",
                "counterpart_ids",
                "strand_mismatch_basis",
                "multi_segment",
                "any_pair_chain",
            ]
            data = (
                det_ref[cols]
                .astype(object)
                .where(det_ref[cols].notna(), None)
                .values.tolist()
            )
            con.executemany(
                f"INSERT INTO ref_outcomes VALUES (?,?,?,{','.join('?' * len(cols))})",
                [(ref.id, policy.id, run.run_id, *row) for row in data],
            )
            qcols = [
                "gene_id",
                "original_gene_id",
                "chrom",
                "start",
                "end",
                "strand",
                "classification",
                "classification_cds",
                "cds_exact",
                "cds_overlap",
                "cds_matched_id",
                "best_cds_ref_tx",
                "best_cds_query_tx",
                "novel_category",
                "counterpart_count",
                "counterpart_ids",
                "strand_mismatch_basis",
            ]
            con.executemany(
                f"INSERT INTO query_outcomes VALUES (?,?,?,{','.join('?' * len(qcols))})",
                [
                    (ref.id, policy.id, run.run_id, *row)
                    for row in det_q[qcols].astype(object).values.tolist()
                ],
            )
        _set_component(con, component, key)
    log(f"  result {run.run_id}/{policy.id}: {task.status}")


def run_colors(experiment: Experiment) -> dict[str, str]:
    """Colour by tool (stable across references); later runs of the same tool on one reference get lighter tints."""
    tools = []
    for run in experiment.runs.values():
        if run.tool not in tools:
            tools.append(run.tool)
    seen: dict[tuple, int] = {}
    colors = {}
    for run in experiment.runs.values():
        if run.color:
            colors[run.run_id] = run.color
            continue
        idx = tools.index(run.tool)
        base = PALETTE[idx] if idx < len(PALETTE) else "#8a8a85"
        k = seen.get((run.reference, run.tool), 0)
        seen[(run.reference, run.tool)] = k + 1
        colors[run.run_id] = _tint(base, 0.38 * k) if k else base
    return colors


def _tint(hex_color: str, amount: float) -> str:
    amount = min(amount, 0.75)
    r, g, b = (int(hex_color[i : i + 2], 16) for i in (1, 3, 5))
    mix = lambda c: round(c + (255 - c) * amount)
    return f"#{mix(r):02x}{mix(g):02x}{mix(b):02x}"


def _read_text(path: Path) -> str | None:
    try:
        return Path(path).read_text()
    except OSError:
        return None


def build_dashboard(experiment: Experiment, ws: Workspace, log=print) -> Path:
    db = ws.dashboard_db()
    con = connect(db)
    tasks = plan(experiment, ws, prepare=False)
    prepared = {
        ref_id: ws.preparation_record(ref_id) for ref_id in experiment.references
    }
    ready = {}
    for ref in experiment.references.values():
        prep = prepared[ref.id]
        if prep is None or not Path(prep["annotation"]["output"]["path"]).exists():
            log(
                f"  reference {ref.id}: prepared annotation unavailable; locus browser and per-gene segment counts disabled"
            )
            ready[ref.id] = False
            continue
        build_reference(con, experiment, ws, ref, prep, log)
        ready[ref.id] = True
    model_errors = {}
    for run in experiment.runs.values():
        prep = prepared.get(run.reference)
        annotation = (ws.run_preparation_record(run.run_id) or {}).get(
            "output", {}
        ).get("path") or str(run.annotation)
        if not (prep and Path(annotation).exists()):
            model_errors[run.run_id] = "annotation or prepared reference unavailable"
            continue
        try:
            build_run_models(con, experiment, ws, run, prep, log)
        except (
            ValueError,
            OSError,
        ) as error:  # the comparator parser rejected the file; record, do not repair
            model_errors[run.run_id] = f"{type(error).__name__}: {error}"
            with con:
                con.execute("DELETE FROM query_genes WHERE run_id=?", (run.run_id,))
                con.execute(
                    "DELETE FROM transcripts WHERE source=?", (f"run:{run.run_id}",)
                )
                con.execute(
                    "DELETE FROM components WHERE component=?", (f"run:{run.run_id}",)
                )
            log(
                f"  run {run.run_id}: models unavailable ({model_errors[run.run_id][:120]})"
            )
    for task in tasks:
        build_result(
            con, experiment, ws, task, ready.get(task.reference.id, False), log
        )

    validation = read_json(ws.validation_file())
    meta = {
        "store_version": STORE_VERSION,
        "built": now(),
        "experiment": {
            "id": experiment.id,
            "title": experiment.title,
            "description": experiment.description,
            "synthetic": experiment.synthetic,
            "config_path": str(experiment.path),
        },
        "config": experiment.raw,
        "evaluation": experiment.evaluation,
        "policies": [p.__dict__ for p in experiment.policies.values()],
        "references": [
            {
                **{
                    k: (str(v) if isinstance(v, Path) else v)
                    for k, v in r.__dict__.items()
                    if k not in ("raw", "annotation_preparation", "genome_preparation")
                },
                "models_available": ready.get(r.id, False),
                "preparation": prepared.get(r.id),
            }
            for r in experiment.references.values()
        ],
        "runs": [
            {
                "run_id": r.run_id,
                "reference": r.reference,
                "tool": r.tool,
                "tool_version": r.tool_version,
                "model": r.model,
                "annotation": str(r.annotation),
                "format": r.format,
                "provenance": r.provenance,
                "caveats": r.caveats,
                "synthetic": r.synthetic,
                "color": c,
            }
            for r, c in (
                (r, run_colors(experiment)[r.run_id]) for r in experiment.runs.values()
            )
        ],
        "validation": validation,
        "model_errors": model_errors,
        "import_report": read_json(ws.import_report()),
        "metric_dictionary": json.loads(METRIC_DICTIONARY.read_text()),
        "context_documents": [
            {"path": str(p), "name": Path(p).name, "text": _read_text(p)}
            for p in experiment.imports.get("context_documents") or []
        ],
        "paper_analogous": [
            {"path": str(p), "name": Path(p).name, "data": read_json(p)}
            for p in experiment.imports.get("paper_analogous") or []
        ],
    }
    with con:
        con.executemany(
            "INSERT OR REPLACE INTO meta VALUES (?,?)",
            [(k, json.dumps(v, default=str)) for k, v in meta.items()],
        )
    con.execute("PRAGMA optimize")
    con.close()
    log(f"Dashboard data: {db}")
    return db

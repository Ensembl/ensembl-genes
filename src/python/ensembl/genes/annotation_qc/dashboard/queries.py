"""
Read-only queries over a dashboard dataset (benchmark/store.py, store version 4).

Every metric returned is a stored row of the comparator's metrics.tsv (or a
legacy-import row), selected by (reference, policy, scope, run, metric id); this
module filters, indexes and formats only. All locus queries are bounded: a window
is at most MAX_WINDOW bp, overlap searches use the per-sequence maximum feature
length so the (source, chrom, start) index applies, and each track returns at
most a fixed number of records with a ``truncated`` flag. Nothing here re-reads
annotation files.

Datasets of earlier store versions are refused with a message (they used
dashboard-side metric derivations); build a new workspace instead.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

METRIC_DICTIONARY = (
    Path(__file__).resolve().parent.parent / "benchmark" / "metric_dictionary.json"
)
SUPPORTED_STORE = "4"
MAX_WINDOW = 2_000_000
TRACK_LIMIT = 60
CONTEXT_LIMIT = 80
PAGE_LIMIT = 500
CSV_LIMIT = 1_000_000
ALL = "all"

CDS_OUTCOMES = (
    "exact",
    "terminal_diff",
    "diff_chain",
    "partial_lt08",
    "span_only",
    "strand_mismatch",
    "missed",
)
CHROM_ORDER = (
    "(CASE WHEN {c} GLOB '[0-9]*' THEN CAST({c} AS INTEGER) ELSE 100000 END), {c}"
)


class QueryError(ValueError):
    pass


def connect(db: Path) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True, check_same_thread=False)
    con.row_factory = sqlite3.Row
    return con


def store_version(con) -> str | None:
    try:
        row = con.execute("SELECT value FROM meta WHERE key='store_version'").fetchone()
    except sqlite3.OperationalError:
        return None
    return json.loads(row[0]) if row else None


def check_version(con) -> None:
    version = store_version(con)
    if version != SUPPORTED_STORE:
        raise QueryError(
            f"dataset store version {version!r} is not supported (needs {SUPPORTED_STORE}); "
            "rebuild it with `annotation-qc benchmark build` into a new workspace"
        )


def meta(con) -> dict:
    check_version(con)
    data = {
        row["key"]: json.loads(row["value"])
        for row in con.execute("SELECT key, value FROM meta")
    }
    data["metric_dictionary_at_build"] = data.get("metric_dictionary")
    data["metric_dictionary"] = json.loads(METRIC_DICTIONARY.read_text())
    return data


def _results(con, reference_id: str) -> list[dict]:
    out = []
    for r in con.execute(
        "SELECT run_id, policy_id, status, reason, source, record FROM results WHERE reference_id=?",
        (reference_id,),
    ):
        r = dict(r)
        record = json.loads(r.pop("record") or "{}")
        for key in (
            "identity",
            "outputs_check",
            "analyses",
            "transformations",
            "provenance",
            "pairs_check",
            "independent",
            "stop_codon_harmonisation",
            "limitations",
            "output_dir",
            "finished",
            "wall_seconds",
            "peak_rss_bytes",
            "difference",
            "attempts",
            "contract_version",
            "definition_version",
            "legacy",
            "options",
            "inputs",
            "scope",
            "thresholds",
            "gene_span_checks",
            "output_files",
            "assessment_reason",
            "imported",
        ):
            r[key] = record.get(key)
        out.append(r)
    return out


def overview(con, reference_id: str) -> dict:
    """Result status and every stored metric row (all scopes) of one reference."""
    metrics = [
        dict(r)
        for r in con.execute(
            "SELECT run_id, policy_id, scope_id, metric_id, numerator, denominator, value, value_text, "
            "status, reason, numerator_id, denominator_id, unit, side, source FROM metrics "
            "WHERE reference_id=?",
            (reference_id,),
        )
    ]
    return {"results": _results(con, reference_id), "metrics": metrics}


def metrics_rows(
    con, reference_id: str, policy: str | None, scope: str | None
) -> list[dict]:
    """Metric rows for download (optionally one policy and scope)."""
    where, args = ["reference_id=?"], [reference_id]
    for col, value in (("policy_id", policy), ("scope_id", scope)):
        if value:
            where.append(f"{col}=?")
            args.append(value)
    return [
        dict(r)
        for r in con.execute(
            "SELECT reference_id, run_id, policy_id, scope_id, metric_id, numerator, denominator, value_text AS value, "
            f"status, reason, numerator_id, denominator_id, unit, side, source FROM metrics WHERE {' AND '.join(where)} "
            "ORDER BY run_id, policy_id, scope_id, metric_id",
            args,
        )
    ]


def splice_sites(con, p: dict) -> dict:
    """The comparator's binned splice-site and intron-length tables for one policy and scope."""
    args = (p["reference"], p.get("scope") or ALL, p["policy"])
    where = "reference_id=? AND scope_id=? AND policy_id=?"
    sites = [
        dict(r)
        for r in con.execute(
            f"SELECT run_id, min_intron_length, coordinate, kind, bin, lo, hi, assessed, supported, missing_locus "
            f"FROM splice_sites WHERE {where} ORDER BY run_id, min_intron_length, coordinate, kind, bin",
            args,
        )
    ]
    lengths = [
        dict(r)
        for r in con.execute(
            f"SELECT run_id, bin, lo, hi, label, eligible, exact, donor, acceptor, missing_locus "
            f"FROM intron_lengths WHERE {where} ORDER BY run_id, bin",
            args,
        )
    ]
    return {"scope": args[1], "policy": args[2], "sites": sites, "lengths": lengths}


def per_sequence(con, p: dict) -> dict:
    """Per-sequence drill-down rows of one result (never summed into scope values)."""
    rows = [
        dict(r)
        for r in con.execute(
            f"SELECT * FROM per_sequence WHERE reference_id=? AND policy_id=? AND run_id=? "
            f"ORDER BY {CHROM_ORDER.format(c='chrom')}",
            (p["reference"], p["policy"], p["run"]),
        )
    ]
    return {"rows": rows}


def _int(value, name, default=None):
    if value in (None, ""):
        return default
    try:
        return int(value)
    except ValueError as error:
        raise QueryError(f"{name} must be an integer") from error


def _scope_clause(alias: str, scope: str | None) -> tuple[str, list]:
    if not scope or scope == ALL:
        return "", []
    return (
        f" AND (',' || COALESCE({alias}.scope_ids, '') || ',') LIKE ?",
        [f"%,{scope},%"],
    )


def genes(con, p: dict, csv=False) -> dict:
    """Evaluated reference genes of one policy (and scope) with per-run outcomes."""
    ref, policy = p["reference"], p["policy"]
    where = ["g.reference_id=?", "g.policy_id=?"]
    args: list = [ref, policy]
    clause, extra = _scope_clause("g", p.get("scope"))
    if clause:
        where.append(clause[5:])
        args += extra
    if p.get("q"):
        where.append("g.search LIKE ?")
        args.append(f"%{p['q'].strip().lower()}%")
    if p.get("chrom"):
        where.append("g.chrom=?")
        args.append(p["chrom"])
        start, end = _int(p.get("start"), "start"), _int(p.get("end"), "end")
        if start is not None and end is not None:
            where += ["g.start<=?", "g.end>=?"]
            args += [end, start]
    if p.get("multi") == "1":
        where.append("g.max_cds_count>1")
    elif p.get("multi") == "0":
        where.append("g.max_cds_count<=1")
    focus = p.get("focus_run")
    exists = (
        "EXISTS (SELECT 1 FROM ref_outcomes o WHERE o.reference_id=g.reference_id AND "
        "o.policy_id=g.policy_id AND o.gene_id=g.gene_id AND o.run_id=? AND {})"
    )
    if focus and p.get("status"):
        statuses = [s for s in p["status"].split(",") if s in CDS_OUTCOMES]
        if statuses:
            where.append(
                exists.format(f"o.cds_outcome IN ({','.join('?' * len(statuses))})")
            )
            args += [focus, *statuses]
    if focus and p.get("split") == "1":
        where.append(exists.format("o.counterpart_count>=2"))
        args.append(focus)
    if focus and p.get("one_to_one") in ("0", "1"):
        where.append(
            exists.format(
                "o.one_to_one_partner IS NOT NULL"
                if p["one_to_one"] == "1"
                else "o.one_to_one_partner IS NULL"
            )
        )
        args.append(focus)
    compare = p.get("compare_run")
    if focus and compare and p.get("contrast"):
        good, bad = (
            (focus, compare) if p["contrast"] == "exact_in_focus" else (compare, focus)
        )
        where.append(exists.format("o.cds_exact=1"))
        where.append(exists.format("o.cds_exact=0"))
        args += [good, bad]
    base = f"FROM ref_gene_policy g WHERE {' AND '.join(where)}"
    total = con.execute(f"SELECT COUNT(*) {base}", args).fetchone()[0]
    limit = CSV_LIMIT if csv else min(_int(p.get("limit"), "limit", 100), PAGE_LIMIT)
    offset = 0 if csv else _int(p.get("offset"), "offset", 0)
    rows = [
        dict(r)
        for r in con.execute(
            f"SELECT g.gene_id, g.chrom, g.start, g.end, g.strand, g.biotype, g.max_cds_count, "
            f"g.max_exon_count, g.evaluated_transcript_count, g.has_cds, g.scope_ids, "
            f"(SELECT name FROM ref_genes rg WHERE rg.reference_id=g.reference_id AND rg.gene_id=g.gene_id) AS name, "
            f"(SELECT original_gene_id FROM ref_genes rg WHERE rg.reference_id=g.reference_id AND rg.gene_id=g.gene_id) "
            f"AS original_gene_id {base} ORDER BY {CHROM_ORDER.format(c='g.chrom')}, g.start LIMIT ? OFFSET ?",
            [*args, limit, offset],
        )
    ]
    ids = [r["gene_id"] for r in rows]
    outcomes: dict = {}
    for i in range(0, len(ids), 900):
        chunk = ids[i : i + 900]
        for o in con.execute(
            f"SELECT run_id, gene_id, cds_outcome, classification, classification_cds, cds_exact, cds_overlap, "
            f"cds_shared_bases, cds_locus_recovered, cds_matched_id, one_to_one_partner, counterpart_count, "
            f"counterpart_ids, cds_intron_chain_match, cds_intron_chain_any_pair FROM ref_outcomes "
            f"WHERE reference_id=? AND policy_id=? AND gene_id IN ({','.join('?' * len(chunk))})",
            [ref, policy, *chunk],
        ):
            outcomes.setdefault(o["gene_id"], {})[o["run_id"]] = dict(o)
    for r in rows:
        r["outcomes"] = outcomes.get(r["gene_id"], {})
    return {"total": total, "offset": offset, "rows": rows}


def _stop_audit(con, run_id: str, transcript_ids: list[str]) -> dict:
    """transcript id -> stop-codon audit row (stop-harmonised runs only)."""
    found = {}
    ids = [t for t in transcript_ids if t]
    for i in range(0, len(ids), 900):
        chunk = ids[i : i + 900]
        for r in con.execute(
            f"SELECT * FROM stop_codon_audit WHERE run_id=? AND transcript_id IN ({','.join('?' * len(chunk))})",
            [run_id, *chunk],
        ):
            found[r["transcript_id"]] = dict(r)
    return found


def pairs(con, p: dict, csv=False) -> dict:
    """The comparator's one-to-one exact-CDS pairs for one result (optionally one scope)."""
    where = "reference_id=? AND policy_id=? AND run_id=?"
    args = [p["reference"], p["policy"], p["run"]]
    total = con.execute(
        f"SELECT COUNT(*) FROM exact_pairs WHERE {where}", args
    ).fetchone()[0]
    limit = CSV_LIMIT if csv else min(_int(p.get("limit"), "limit", 100), PAGE_LIMIT)
    offset = 0 if csv else _int(p.get("offset"), "offset", 0)
    rows = [
        dict(r)
        for r in con.execute(
            f"SELECT * FROM exact_pairs WHERE {where} ORDER BY {CHROM_ORDER.format(c='chrom')}, cds_start, "
            f"reference_gene_id LIMIT ? OFFSET ?",
            [*args, limit, offset],
        )
    ]
    return {"total": total, "offset": offset, "rows": rows}


def stop_audit(con, p: dict, csv=False) -> dict:
    """Per-transcript stop-codon audit of a stop-harmonised run (optionally by status/reason)."""
    where, args = ["run_id=?"], [p["run"]]
    for key in ("status", "reason"):
        if p.get(key):
            where.append(f"{key}=?")
            args.append(p[key])
    clause = " AND ".join(where)
    total = con.execute(
        f"SELECT COUNT(*) FROM stop_codon_audit WHERE {clause}", args
    ).fetchone()[0]
    limit = CSV_LIMIT if csv else min(_int(p.get("limit"), "limit", 100), PAGE_LIMIT)
    offset = 0 if csv else _int(p.get("offset"), "offset", 0)
    rows = [
        dict(r)
        for r in con.execute(
            f"SELECT * FROM stop_codon_audit WHERE {clause} ORDER BY {CHROM_ORDER.format(c='chrom')}, "
            f"terminal_cds_start LIMIT ? OFFSET ?",
            [*args, limit, offset],
        )
    ]
    counts = [
        dict(r)
        for r in con.execute(
            "SELECT status, reason, COUNT(*) AS transcripts FROM stop_codon_audit WHERE run_id=? "
            "GROUP BY status, reason ORDER BY status, reason",
            (p["run"],),
        )
    ]
    return {"total": total, "offset": offset, "rows": rows, "counts": counts}


def query_genes(con, p: dict, csv=False) -> dict:
    """Prediction table for one run and policy (Novel context, merges, strand mismatch)."""
    ref, policy, run = p["reference"], p["policy"], p.get("run")
    if not run:
        raise QueryError("run is required")
    where = ["o.reference_id=?", "o.policy_id=?", "o.run_id=?"]
    args: list = [ref, policy, run]
    clause, extra = _scope_clause("o", p.get("scope"))
    if clause:
        where.append(clause[5:])
        args += extra
    if p.get("q"):
        where.append("o.search LIKE ?")
        args.append(f"%{p['q'].strip().lower()}%")
    if p.get("chrom"):
        where.append("o.chrom=?")
        args.append(p["chrom"])
    if p.get("classification"):
        where.append("o.classification=?")
        args.append(p["classification"])
    if p.get("novel_category"):
        where.append("o.novel_category=?")
        args.append(p["novel_category"])
    if p.get("merge") == "1":
        where.append("o.counterpart_count>=2")
    if p.get("exact") in ("0", "1"):
        where.append("o.cds_exact=?")
        args.append(int(p["exact"]))
    if p.get("one_to_one") in ("0", "1"):
        where.append(
            "o.one_to_one_partner IS NOT NULL"
            if p["one_to_one"] == "1"
            else "o.one_to_one_partner IS NULL"
        )
    base = f"FROM query_outcomes o WHERE {' AND '.join(where)}"
    total = con.execute(f"SELECT COUNT(*) {base}", args).fetchone()[0]
    limit = CSV_LIMIT if csv else min(_int(p.get("limit"), "limit", 100), PAGE_LIMIT)
    offset = 0 if csv else _int(p.get("offset"), "offset", 0)
    rows = [
        dict(r)
        for r in con.execute(
            f"SELECT o.gene_id, o.original_gene_id, o.chrom, o.start, o.end, o.strand, o.classification, "
            f"o.classification_cds, o.cds_exact, o.cds_overlap, o.cds_matched_id, o.best_cds_ref_tx, "
            f"o.best_cds_query_tx, o.novel_category, o.counterpart_count, o.counterpart_ids, "
            f"o.strand_mismatch_basis, o.one_to_one_partner, o.evaluated_transcript_count, o.has_cds, "
            f"o.scope_ids {base} ORDER BY {CHROM_ORDER.format(c='o.chrom')}, o.start LIMIT ? OFFSET ?",
            [*args, limit, offset],
        )
    ]
    return {"total": total, "offset": offset, "rows": rows}


def _tx_status(con, source: str, policy: str, ids: list[str]) -> dict:
    found = {}
    for i in range(0, len(ids), 900):
        chunk = ids[i : i + 900]
        for r in con.execute(
            f"SELECT transcript_id, status, exclusion_reason, selected_by, exon_count, cds_segment_count, "
            f"cds_length, transformed FROM tx_status WHERE source=? AND policy_id=? "
            f"AND transcript_id IN ({','.join('?' * len(chunk))})",
            [source, policy, *chunk],
        ):
            found[r["transcript_id"]] = dict(r)
    return found


def gene_detail(con, reference: str, policy: str, gene_id: str) -> dict:
    """One reference gene: per-run comparator outcome, compared isoforms and CDS differences."""
    gene = con.execute(
        "SELECT * FROM ref_gene_policy WHERE reference_id=? AND policy_id=? AND gene_id=?",
        (reference, policy, gene_id),
    ).fetchone()
    if gene is None:
        raise QueryError(
            f"reference gene {gene_id} is not evaluated under policy {policy}"
        )
    gene = dict(gene)
    gene.pop("search", None)
    names = con.execute(
        "SELECT name, original_gene_id FROM ref_genes WHERE reference_id=? AND gene_id=?",
        (reference, gene_id),
    ).fetchone()
    gene["name"], gene["original_gene_id"] = (
        (names[0], names[1]) if names else (None, gene_id)
    )
    isoforms = [
        dict(r)
        for r in con.execute(
            "SELECT transcript_id, status, exclusion_reason, selected_by, exon_count, cds_segment_count, "
            "cds_length FROM tx_status WHERE source=? AND policy_id=? AND gene_id=? ORDER BY transcript_id",
            (f"ref:{reference}", policy, gene_id),
        )
    ]
    outcomes = {}
    for o in con.execute(
        "SELECT * FROM ref_outcomes WHERE reference_id=? AND policy_id=? AND gene_id=?",
        (reference, policy, gene_id),
    ):
        o = dict(o)
        o["query_stop_codon"] = (
            _stop_audit(con, o["run_id"], [o["best_cds_query_tx"]]).get(
                o["best_cds_query_tx"]
            )
            if o["best_cds_query_tx"]
            else None
        )
        outcomes[o["run_id"]] = o
    runs = {
        r["run_id"]
        for r in con.execute(
            "SELECT run_id FROM results WHERE reference_id=?", (reference,)
        )
    }
    return {
        "gene": gene,
        "isoforms": isoforms,
        "outcomes": outcomes,
        "runs_without_outcome": sorted(runs - set(outcomes)),
    }


def locus(con, p: dict) -> dict:  # pylint: disable=too-many-locals
    ref = p["reference"]
    chrom = p.get("chrom")
    start, end = _int(p.get("start"), "start"), _int(p.get("end"), "end")
    if not chrom or start is None or end is None:
        raise QueryError("chrom, start and end are required")
    if end < start:
        start, end = end, start
    clamped = False
    if end - start + 1 > MAX_WINDOW:
        mid = (start + end) // 2
        start, end = mid - MAX_WINDOW // 2, mid + MAX_WINDOW // 2
        clamped = True
    policy = p.get("policy")
    runs = [r for r in (p.get("runs") or "").split(",") if r]
    focus_gene = p.get("gene_id")
    show = p.get("isoforms", "partners")  # partners | evaluated | all

    def maxlen(source):
        row = con.execute(
            "SELECT max_len FROM seq_maxlen WHERE source=? AND chrom=?", (source, chrom)
        ).fetchone()
        return row[0] if row else 0

    def overlapping(table, source_col, source, limit, ml_source):
        ml = maxlen(ml_source)
        rows = [
            dict(r)
            for r in con.execute(
                f"SELECT * FROM {table} WHERE {source_col}=? AND chrom=? AND start<=? AND start>=? AND end>=? "
                f"ORDER BY start LIMIT ?",
                (source, chrom, end, start - ml, start, limit + 1),
            )
        ]
        return rows[:limit], len(rows) > limit

    ref_source = f"ref:{ref}"
    ref_genes, genes_trunc = overlapping(
        "ref_genes", "reference_id", ref, CONTEXT_LIMIT, f"genes:{ref_source}"
    )
    for g in ref_genes:
        g.pop("search", None)
    partners: dict = {}
    outcome_by_gene: dict = {}
    if policy and runs and ref_genes:
        ids = [g["gene_id"] for g in ref_genes]
        for o in con.execute(
            f"SELECT run_id, gene_id, best_cds_ref_tx, best_cds_query_tx, cds_outcome, cds_matched_id, "
            f"one_to_one_partner FROM ref_outcomes WHERE reference_id=? AND policy_id=? "
            f"AND gene_id IN ({','.join('?' * len(ids))}) AND run_id IN ({','.join('?' * len(runs))})",
            [ref, policy, *ids, *runs],
        ):
            if o["best_cds_ref_tx"]:
                partners.setdefault(o["best_cds_ref_tx"], []).append(o["run_id"])
            outcome_by_gene.setdefault(o["gene_id"], {})[o["run_id"]] = dict(o)
    ref_tx_all, ref_tx_trunc = overlapping(
        "transcripts", "source", ref_source, 4000, ref_source
    )
    status = (
        _tx_status(con, ref_source, policy, [t["transcript_id"] for t in ref_tx_all])
        if policy
        else {}
    )
    ref_tracks, hidden = [], 0
    for t in ref_tx_all:
        st = status.get(t["transcript_id"]) or {}
        t["evaluated"] = st.get("status") == "evaluated"
        t["exclusion_reason"] = st.get("exclusion_reason")
        t["selected_by"] = st.get("selected_by")
        t["partner_of_runs"] = partners.get(t["transcript_id"], [])
        t["canonical"] = "Ensembl_canonical" in (t["tags"] or "")
        keep = (
            show == "all"
            or (show == "evaluated" and t["evaluated"])
            or (
                show == "partners"
                and t["evaluated"]
                and (
                    t["partner_of_runs"] or t["canonical"] or t["selected_by"] != "all"
                )
            )
        )
        if keep:
            ref_tracks.append(t)
        else:
            hidden += 1
    ref_tracks.sort(
        key=lambda t: (t["gene_id"] != focus_gene, not t["partner_of_runs"], t["start"])
    )
    ref_trunc = len(ref_tracks) > TRACK_LIMIT
    ref_tracks = ref_tracks[:TRACK_LIMIT]
    run_tracks = {}
    for run in runs:
        source = f"run:{run}"
        txs, trunc = overlapping(
            "transcripts", "source", source, TRACK_LIMIT * 4, source
        )
        st = (
            _tx_status(con, source, policy, [t["transcript_id"] for t in txs])
            if policy
            else {}
        )
        q_out = {}
        ids = list({t["gene_id"] for t in txs})
        if policy and ids:
            for o in con.execute(
                f"SELECT gene_id, classification, classification_cds, cds_exact, novel_category, counterpart_count, "
                f"one_to_one_partner FROM query_outcomes WHERE run_id=? AND policy_id=? "
                f"AND gene_id IN ({','.join('?' * len(ids))})",
                [run, policy, *ids],
            ):
                q_out[o["gene_id"]] = dict(o)
        stops = _stop_audit(con, run, [t["transcript_id"] for t in txs])
        shown, run_hidden = [], 0
        for t in txs:
            s = st.get(t["transcript_id"]) or {}
            t["evaluated"] = s.get("status") == "evaluated"
            t["exclusion_reason"] = s.get("exclusion_reason")
            t["selected_by"] = s.get("selected_by")
            t["outcome"] = q_out.get(t["gene_id"])
            t["stop_codon"] = stops.get(t["transcript_id"])
            if show == "all" or t["evaluated"]:
                shown.append(t)
            else:
                run_hidden += 1
        run_tracks[run] = {
            "transcripts": shown[:TRACK_LIMIT],
            "truncated": trunc or len(shown) > TRACK_LIMIT,
            "hidden_excluded": run_hidden,
        }
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "window_clamped": clamped,
        "reference_genes": ref_genes,
        "reference_genes_truncated": genes_trunc,
        "reference_transcripts": ref_tracks,
        "reference_transcripts_truncated": ref_trunc or ref_tx_trunc,
        "reference_transcripts_hidden": hidden,
        "isoform_mode": show,
        "runs": run_tracks,
        "outcomes": outcome_by_gene,
        "limits": {
            "max_window": MAX_WINDOW,
            "track_limit": TRACK_LIMIT,
            "context_limit": CONTEXT_LIMIT,
        },
    }


def chromosomes(con, reference: str) -> list[str]:
    return [
        r[0]
        for r in con.execute(
            f"SELECT DISTINCT chrom FROM ref_gene_policy WHERE reference_id=? "
            f"ORDER BY {CHROM_ORDER.format(c='chrom')}",
            (reference,),
        )
    ]

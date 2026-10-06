"""
Read-only queries over a dashboard dataset (benchmark/store.py schema).

All locus queries are bounded: a window is at most MAX_WINDOW bp, overlap
searches use the per-sequence maximum feature length so the (source, chrom, start)
index applies, and each track returns at most a fixed number of records with a
``truncated`` flag. Nothing here re-reads annotation files.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

MAX_WINDOW = 2_000_000
TRACK_LIMIT = 60
CONTEXT_LIMIT = 80
PAGE_LIMIT = 500
CSV_LIMIT = 1_000_000

REF_STATUSES = (
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


def meta(con) -> dict:
    return {
        row["key"]: json.loads(row["value"])
        for row in con.execute("SELECT key, value FROM meta")
    }


def overview(con, reference_id: str) -> dict:
    """Result status and metric rows for every run and policy of one reference."""
    results = [
        dict(r)
        for r in con.execute(
            "SELECT run_id, policy_id, status, reason, source, record FROM results WHERE reference_id=?",
            (reference_id,),
        )
    ]
    for r in results:
        record = json.loads(r.pop("record") or "{}")
        r["independent"] = record.get("independent")
        r["wall_seconds"] = record.get("wall_seconds")
        r["peak_rss_bytes"] = record.get("peak_rss_bytes")
        r["output_dir"] = record.get("output_dir")
        r["finished"] = record.get("finished")
        r["difference"] = record.get("difference")
        r["attempts"] = record.get("attempts")
    metrics = [
        dict(r)
        for r in con.execute(
            "SELECT run_id, policy_id, metric_id, numerator, denominator, value, status, source, note FROM metrics WHERE reference_id=?",
            (reference_id,),
        )
    ]
    return {"results": results, "metrics": metrics}


def _int(value, name, default=None):
    if value in (None, ""):
        return default
    try:
        return int(value)
    except ValueError as error:
        raise QueryError(f"{name} must be an integer") from error


def genes(con, p: dict, csv=False) -> dict:
    """Reference-gene table for one reference and policy, with per-run outcomes."""
    ref, policy = p["reference"], p["policy"]
    where = ["g.reference_id=?", "rp.policy_id=?"]
    args: list = [ref, policy]
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
    if p.get("multi") in ("0", "1"):
        where.append("rp.multi_segment=?")
        args.append(int(p["multi"]))
    focus = p.get("focus_run")
    if focus and p.get("status"):
        statuses = [s for s in p["status"].split(",") if s in REF_STATUSES]
        if statuses:
            where.append(
                f"EXISTS (SELECT 1 FROM ref_outcomes o WHERE o.reference_id=g.reference_id AND o.policy_id=rp.policy_id "
                f"AND o.gene_id=g.gene_id AND o.run_id=? AND o.cds_status IN ({','.join('?' * len(statuses))}))"
            )
            args += [focus, *statuses]
    if focus and p.get("split") == "1":
        where.append(
            "EXISTS (SELECT 1 FROM ref_outcomes o WHERE o.reference_id=g.reference_id AND o.policy_id=rp.policy_id "
            "AND o.gene_id=g.gene_id AND o.run_id=? AND o.counterpart_count>=2)"
        )
        args.append(focus)
    compare = p.get("compare_run")
    if focus and compare and p.get("contrast"):
        # focus run exact, compare run not exact (or the reverse when contrast=missing_in_focus)
        good, bad = (
            (focus, compare) if p["contrast"] == "exact_in_focus" else (compare, focus)
        )
        where.append(
            "EXISTS (SELECT 1 FROM ref_outcomes o WHERE o.reference_id=g.reference_id AND o.policy_id=rp.policy_id "
            "AND o.gene_id=g.gene_id AND o.run_id=? AND o.cds_exact=1)"
        )
        where.append(
            "EXISTS (SELECT 1 FROM ref_outcomes o WHERE o.reference_id=g.reference_id AND o.policy_id=rp.policy_id "
            "AND o.gene_id=g.gene_id AND o.run_id=? AND o.cds_exact=0)"
        )
        args += [good, bad]
    clause = " AND ".join(where)
    base = f"FROM ref_genes g JOIN ref_gene_policy rp ON rp.reference_id=g.reference_id AND rp.gene_id=g.gene_id WHERE {clause}"
    total = con.execute(f"SELECT COUNT(*) {base}", args).fetchone()[0]
    limit = CSV_LIMIT if csv else min(_int(p.get("limit"), "limit", 100), PAGE_LIMIT)
    offset = 0 if csv else _int(p.get("offset"), "offset", 0)
    rows = [
        dict(r)
        for r in con.execute(
            f"SELECT g.gene_id, g.original_gene_id, g.name, g.chrom, g.start, g.end, g.strand, g.biotype, rp.multi_segment, "
            f"rp.canonical_fallback, rp.evaluated_transcripts {base} ORDER BY {CHROM_ORDER.format(c='g.chrom')}, g.start LIMIT ? OFFSET ?",
            [*args, limit, offset],
        )
    ]
    ids = [r["gene_id"] for r in rows]
    outcomes: dict = {}
    for i in range(0, len(ids), 900):
        chunk = ids[i : i + 900]
        for o in con.execute(
            f"SELECT run_id, gene_id, cds_status, classification, classification_cds, cds_exact, cds_overlap, cds_matched_id, "
            f"counterpart_count, counterpart_ids, any_pair_chain, cds_intron_chain_match FROM ref_outcomes "
            f"WHERE reference_id=? AND policy_id=? AND gene_id IN ({','.join('?' * len(chunk))})",
            [ref, policy, *chunk],
        ):
            outcomes.setdefault(o["gene_id"], {})[o["run_id"]] = dict(o)
    for r in rows:
        r["evaluated_transcript_count"] = (
            len(r.pop("evaluated_transcripts").split(","))
            if r.get("evaluated_transcripts")
            else 0
        )
        r["outcomes"] = outcomes.get(r["gene_id"], {})
    return {"total": total, "offset": offset, "rows": rows}


def query_genes(con, p: dict, csv=False) -> dict:
    """Prediction table for one run and policy (Novel context, merges, strand mismatch)."""
    ref, policy, run = p["reference"], p["policy"], p.get("run")
    if not run:
        raise QueryError("run is required")
    where = ["o.reference_id=?", "o.policy_id=?", "o.run_id=?"]
    args: list = [ref, policy, run]
    if p.get("q"):
        where.append("qg.search LIKE ?")
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
    base = (
        f"FROM query_outcomes o LEFT JOIN query_genes qg ON qg.run_id=o.run_id AND qg.gene_id=o.gene_id "
        f"WHERE {' AND '.join(where)}"
    )
    total = con.execute(f"SELECT COUNT(*) {base}", args).fetchone()[0]
    limit = CSV_LIMIT if csv else min(_int(p.get("limit"), "limit", 100), PAGE_LIMIT)
    offset = 0 if csv else _int(p.get("offset"), "offset", 0)
    rows = [
        dict(r)
        for r in con.execute(
            f"SELECT o.gene_id, o.original_gene_id, qg.name, o.chrom, o.start, o.end, o.strand, o.classification, o.classification_cds, "
            f"o.cds_exact, o.cds_overlap, o.cds_matched_id, o.best_cds_ref_tx, o.novel_category, o.counterpart_count, o.counterpart_ids, "
            f"o.strand_mismatch_basis {base} ORDER BY {CHROM_ORDER.format(c='o.chrom')}, o.start LIMIT ? OFFSET ?",
            [*args, limit, offset],
        )
    ]
    return {"total": total, "offset": offset, "rows": rows}


def _parse_iv(text: str) -> list[tuple[int, int]]:
    return (
        [tuple(int(x) for x in part.split("-")) for part in text.split(",") if part]
        if text
        else []
    )


def _tx(con, source: str, transcript_id: str):
    row = con.execute(
        "SELECT * FROM transcripts WHERE source=? AND transcript_id=?",
        (source, transcript_id),
    ).fetchone()
    return dict(row) if row else None


def explain_pair(ref_tx: dict | None, q_tx: dict | None) -> dict:
    """Describe CDS differences between two stored transcripts (one-based inclusive)."""
    if not ref_tx or not q_tx:
        return {"available": False}
    a, b = _parse_iv(ref_tx["cds"]), _parse_iv(q_tx["cds"])
    if not a or not b:
        return {"available": False}
    ia = [(a[i][1] + 1, a[i + 1][0] - 1) for i in range(len(a) - 1)]
    ib = [(b[i][1] + 1, b[i + 1][0] - 1) for i in range(len(b) - 1)]
    plus = ref_tx["strand"] != "-"
    # positive = query CDS extends further than the reference at that end
    start_diff = (a[0][0] - b[0][0]) if plus else (b[-1][1] - a[-1][1])
    stop_diff = (b[-1][1] - a[-1][1]) if plus else (a[0][0] - b[0][0])
    shared = sum(max(0, min(e1, e2) - max(s1, s2) + 1) for s1, e1 in a for s2, e2 in b)
    la, lb = sum(e - s + 1 for s, e in a), sum(e - s + 1 for s, e in b)
    out = {
        "available": True,
        "identical_cds": a == b,
        "same_cds_intron_chain": ia == ib,
        "reference_cds_segments": len(a),
        "query_cds_segments": len(b),
        "reference_cds_bp": la,
        "query_cds_bp": lb,
        "reciprocal_cds_overlap": (
            round(min(shared / la, shared / lb), 4) if la and lb else 0
        ),
        "start_side_diff_bp": start_diff,
        "stop_side_diff_bp": stop_diff,
        "reference_introns_missing_in_query": [
            f"{s}-{e}" for s, e in ia if (s, e) not in set(ib)
        ][:10],
        "query_introns_not_in_reference": [
            f"{s}-{e}" for s, e in ib if (s, e) not in set(ia)
        ][:10],
        "reference_introns_missing_count": sum(1 for x in ia if x not in set(ib)),
        "query_introns_extra_count": sum(1 for x in ib if x not in set(ia)),
    }
    notes = []
    if out["identical_cds"]:
        notes.append("CDS identical (start, every splice site and stop)")
    else:
        if out["same_cds_intron_chain"]:
            notes.append("same CDS intron chain")
        if start_diff:
            notes.append(
                f"start side: query CDS {'longer' if start_diff > 0 else 'shorter'} by {abs(start_diff)} bp"
            )
        if stop_diff:
            notes.append(
                f"stop side: query CDS {'longer' if stop_diff > 0 else 'shorter'} by {abs(stop_diff)} bp"
            )
        if out["reference_introns_missing_count"]:
            notes.append(
                f"{out['reference_introns_missing_count']} reference CDS intron(s) absent in query"
            )
        if out["query_introns_extra_count"]:
            notes.append(
                f"{out['query_introns_extra_count']} query CDS intron(s) not in reference"
            )
    out["summary"] = "; ".join(notes)
    return out


def gene_detail(con, reference: str, policy: str, gene_id: str) -> dict:
    gene = con.execute(
        "SELECT * FROM ref_genes WHERE reference_id=? AND gene_id=?",
        (reference, gene_id),
    ).fetchone()
    if gene is None:
        raise QueryError(f"reference gene {gene_id} not found")
    gene = dict(gene)
    gene.pop("search", None)
    rp = con.execute(
        "SELECT * FROM ref_gene_policy WHERE reference_id=? AND policy_id=? AND gene_id=?",
        (reference, policy, gene_id),
    ).fetchone()
    gene["policy"] = dict(rp) if rp else None
    runs = {
        r["run_id"]
        for r in con.execute(
            "SELECT run_id FROM results WHERE reference_id=?", (reference,)
        )
    }
    outcomes = {}
    for o in con.execute(
        "SELECT * FROM ref_outcomes WHERE reference_id=? AND policy_id=? AND gene_id=?",
        (reference, policy, gene_id),
    ):
        o = dict(o)
        ref_tx = (
            _tx(con, f"ref:{reference}", o["best_cds_ref_tx"])
            if o["best_cds_ref_tx"]
            else None
        )
        q_tx = (
            _tx(con, f"run:{o['run_id']}", o["best_cds_query_tx"])
            if o["best_cds_query_tx"]
            else None
        )
        o["explanation"] = explain_pair(ref_tx, q_tx)
        o["partner_query_gene"] = o["cds_matched_id"] or o["matched_id"]
        outcomes[o["run_id"]] = o
    missing = sorted(runs - set(outcomes))
    return {"gene": gene, "outcomes": outcomes, "runs_without_outcome": missing}


def locus(con, p: dict) -> dict:
    ref = p["reference"]
    chrom = p.get("chrom")
    start, end = _int(p.get("start"), "start"), _int(p.get("end"), "end")
    if not chrom or start is None or end is None:
        raise QueryError("chrom, start and end are required")
    if end < start:
        start, end = end, start
    truncated_window = False
    if end - start + 1 > MAX_WINDOW:
        mid = (start + end) // 2
        start, end = mid - MAX_WINDOW // 2, mid + MAX_WINDOW // 2
        truncated_window = True
    policy = p.get("policy")
    runs = [r for r in (p.get("runs") or "").split(",") if r]
    focus_gene = p.get("gene_id")
    show = p.get("isoforms", "partners")  # partners | evaluated | all

    def maxlen(source):
        row = con.execute(
            "SELECT max_len FROM seq_maxlen WHERE source=? AND chrom=?", (source, chrom)
        ).fetchone()
        return row[0] if row else 0

    def overlapping(table, source_col, source, limit, extra="", extra_args=()):
        ml = maxlen(source if table == "transcripts" else f"refgene:{ref}")
        sql = (
            f"SELECT * FROM {table} WHERE {source_col}=? AND chrom=? AND start<=? AND start>=? AND end>=? {extra} "
            f"ORDER BY start LIMIT ?"
        )
        rows = [
            dict(r)
            for r in con.execute(
                sql, (source, chrom, end, start - ml, start, *extra_args, limit + 1)
            )
        ]
        return rows[:limit], len(rows) > limit

    ref_genes, genes_trunc = overlapping(
        "ref_genes", "reference_id", ref, CONTEXT_LIMIT
    )
    for g in ref_genes:
        g.pop("search", None)
    evaluated = {}
    if policy:
        ids = [g["gene_id"] for g in ref_genes]
        for i in range(0, len(ids), 900):
            chunk = ids[i : i + 900]
            for r in con.execute(
                f"SELECT gene_id, evaluated_transcripts, multi_segment FROM ref_gene_policy WHERE reference_id=? "
                f"AND policy_id=? AND gene_id IN ({','.join('?' * len(chunk))})",
                [ref, policy, *chunk],
            ):
                evaluated[r["gene_id"]] = set(
                    (r["evaluated_transcripts"] or "").split(",")
                )
    partners = {}  # ref transcript -> [run_ids] ; query transcript -> ref gene
    outcome_by_gene = {}
    if policy and runs:
        ids = [g["gene_id"] for g in ref_genes]
        for i in range(0, len(ids), 900):
            chunk = ids[i : i + 900]
            for o in con.execute(
                f"SELECT run_id, gene_id, best_cds_ref_tx, best_cds_query_tx, cds_status, cds_matched_id FROM ref_outcomes "
                f"WHERE reference_id=? AND policy_id=? AND gene_id IN ({','.join('?' * len(chunk))}) "
                f"AND run_id IN ({','.join('?' * len(runs))})",
                [ref, policy, *chunk, *runs],
            ):
                if o["best_cds_ref_tx"]:
                    partners.setdefault(o["best_cds_ref_tx"], []).append(o["run_id"])
                outcome_by_gene.setdefault(o["gene_id"], {})[o["run_id"]] = dict(o)
    ref_tx_all, ref_tx_trunc = overlapping("transcripts", "source", f"ref:{ref}", 4000)
    ref_tracks = []
    hidden = 0
    for t in ref_tx_all:
        ev = t["transcript_id"] in evaluated.get(t["gene_id"], set())
        t["evaluated"] = ev
        t["partner_of_runs"] = partners.get(t["transcript_id"], [])
        t["canonical"] = "Ensembl_canonical" in (t["tags"] or "")
        keep = (
            show == "all"
            or (show == "evaluated" and ev)
            or (
                show == "partners"
                and (
                    t["partner_of_runs"]
                    or (ev and t["canonical"])
                    or (ev and len(evaluated.get(t["gene_id"], ())) == 1)
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
        txs, trunc = overlapping("transcripts", "source", f"run:{run}", TRACK_LIMIT)
        q_out = {}
        if policy:
            ids = list({t["gene_id"] for t in txs})
            if ids:
                for o in con.execute(
                    f"SELECT gene_id, classification, classification_cds, cds_exact, novel_category, counterpart_count "
                    f"FROM query_outcomes WHERE run_id=? AND policy_id=? AND gene_id IN ({','.join('?' * len(ids))})",
                    [run, policy, *ids],
                ):
                    q_out[o["gene_id"]] = dict(o)
        for t in txs:
            t["outcome"] = q_out.get(t["gene_id"])
        run_tracks[run] = {"transcripts": txs, "truncated": trunc}
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "window_clamped": truncated_window,
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
            f"SELECT DISTINCT chrom FROM ref_genes WHERE reference_id=? AND evaluated=1 ORDER BY {CHROM_ORDER.format(c='chrom')}",
            (reference,),
        )
    ]

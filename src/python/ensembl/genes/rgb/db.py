from __future__ import annotations

import contextlib
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import pymysql


@dataclass
class DBParams:
    host: str
    port: int
    user: str
    password: str
    database: str


@contextlib.contextmanager
def connect(params: DBParams):
    conn = pymysql.connect(
        host=params.host,
        port=params.port,
        user=params.user,
        password=params.password,
        database=params.database,
        cursorclass=pymysql.cursors.DictCursor,
        charset="utf8mb4",
        autocommit=True,
        read_timeout=3600,
        write_timeout=3600,
    )
    try:
        yield conn
    finally:
        conn.close()


def list_seq_regions(
    conn,
    coord_system_name: Optional[str] = None,
    allowlist: Optional[Iterable[str]] = None,
) -> List[str]:
    sql = [
        "SELECT sr.name AS name",
        "FROM seq_region sr",
    ]
    params: Tuple = ()
    if coord_system_name:
        sql.append("JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id")
        sql.append("WHERE cs.name = %s")
        params = (coord_system_name,)
    sql.append("ORDER BY name")
    with conn.cursor() as cur:
        cur.execute("\n".join(sql), params)
        rows = [r["name"] for r in cur.fetchall()]
    if allowlist is not None:
        allow = set(allowlist)
        rows = [r for r in rows if r in allow]
    return rows


def _gene_query(coord_system_name: Optional[str]) -> str:
    where_cs = (
        "JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id AND cs.name = %s"
        if coord_system_name
        else ""
    )
    return f"""
SELECT
  g.gene_id,
  g.stable_id,
  sr.name AS seq_region_name,
  g.seq_region_start AS seq_region_start,
  g.seq_region_end AS seq_region_end,
  g.seq_region_strand AS seq_region_strand,
  g.biotype,
  g.canonical_transcript_id,
  a.logic_name
FROM gene g
JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
{where_cs}
LEFT JOIN analysis a ON g.analysis_id = a.analysis_id
WHERE sr.name = %s
ORDER BY g.seq_region_start, g.seq_region_end;
"""


def fetch_genes_for_region(
    conn, seq_region_name: str, coord_system_name: Optional[str]
) -> pd.DataFrame:
    sql = _gene_query(coord_system_name)
    params = (
        (coord_system_name, seq_region_name)
        if coord_system_name
        else (seq_region_name,)
    )
    with conn.cursor() as cur:
        cur.execute(sql, params)
        rows = cur.fetchall()
    if not rows:
        return pd.DataFrame(
            columns=[
                "gene_id",
                "stable_id",
                "seq_region_name",
                "seq_region_start",
                "seq_region_end",
                "seq_region_strand",
                "biotype",
                "canonical_transcript_id",
                "logic_name",
            ]
        )
    return pd.DataFrame(rows)


def extract_all_genes(
    conn,
    seq_regions: Iterable[str],
    coord_system_name: Optional[str],
) -> pd.DataFrame:
    dfs = []
    for sr in seq_regions:
        df = fetch_genes_for_region(conn, sr, coord_system_name)
        if not df.empty:
            dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def _transcript_query(coord_system_name: Optional[str]) -> str:
    where_cs = (
        "JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id AND cs.name = %s"
        if coord_system_name
        else ""
    )
    return f"""
SELECT
  t.transcript_id,
  t.gene_id,
  sr.name AS seq_region_name,
  t.seq_region_start AS seq_region_start,
  t.seq_region_end AS seq_region_end,
  t.seq_region_strand AS seq_region_strand,
  t.biotype,
  a.logic_name
FROM transcript t
JOIN seq_region sr ON t.seq_region_id = sr.seq_region_id
{where_cs}
LEFT JOIN analysis a ON t.analysis_id = a.analysis_id
WHERE sr.name = %s
ORDER BY t.seq_region_start, t.seq_region_end;
"""


def fetch_transcripts_for_region(
    conn, seq_region_name: str, coord_system_name: Optional[str]
) -> pd.DataFrame:
    sql = _transcript_query(coord_system_name)
    params = (
        (coord_system_name, seq_region_name)
        if coord_system_name
        else (seq_region_name,)
    )
    with conn.cursor() as cur:
        cur.execute(sql, params)
        rows = cur.fetchall()
    if not rows:
        return pd.DataFrame(
            columns=[
                "transcript_id",
                "gene_id",
                "seq_region_name",
                "seq_region_start",
                "seq_region_end",
                "seq_region_strand",
                "biotype",
                "logic_name",
            ]
        )
    return pd.DataFrame(rows)


def extract_all_transcripts(
    conn,
    seq_regions: Iterable[str],
    coord_system_name: Optional[str],
) -> pd.DataFrame:
    dfs = []
    for sr in seq_regions:
        df = fetch_transcripts_for_region(conn, sr, coord_system_name)
        if not df.empty:
            dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def extract_all_translations(conn) -> pd.DataFrame:
    """Fetch translations table minimally (transcript_id presence suffices)."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT tl.translation_id, tl.transcript_id, tl.start_exon_id, tl.end_exon_id, tl.seq_start, tl.seq_end
            FROM translation tl
            """
        )
        rows = cur.fetchall()
    if not rows:
        return pd.DataFrame(
            columns=[
                "translation_id",
                "transcript_id",
                "start_exon_id",
                "end_exon_id",
                "seq_start",
                "seq_end",
            ]
        )
    return pd.DataFrame(rows)


def _exon_query(coord_system_name: Optional[str]) -> str:
    where_cs = (
        "JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id AND cs.name = %s"
        if coord_system_name
        else ""
    )
    return f"""
SELECT
  et.transcript_id,
  e.exon_id,
  e.stable_id AS exon_stable_id,
  sr.name AS seq_region_name,
  e.seq_region_start AS seq_region_start,
  e.seq_region_end AS seq_region_end,
  e.seq_region_strand AS seq_region_strand,
  et.rank AS exon_rank,
  e.phase,
  e.end_phase
FROM exon e
JOIN exon_transcript et ON e.exon_id = et.exon_id
JOIN seq_region sr ON e.seq_region_id = sr.seq_region_id
{where_cs}
WHERE sr.name = %s
ORDER BY et.transcript_id, et.rank;
"""


def fetch_exons_for_region(
    conn, seq_region_name: str, coord_system_name: Optional[str]
) -> pd.DataFrame:
    sql = _exon_query(coord_system_name)
    params = (
        (coord_system_name, seq_region_name)
        if coord_system_name
        else (seq_region_name,)
    )
    with conn.cursor() as cur:
        cur.execute(sql, params)
        rows = cur.fetchall()
    columns = [
        "transcript_id",
        "exon_id",
        "exon_stable_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "exon_rank",
        "phase",
        "end_phase",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns)


def extract_all_exons(
    conn,
    seq_regions: Iterable[str],
    coord_system_name: Optional[str],
) -> pd.DataFrame:
    dfs = []
    for sr in seq_regions:
        df = fetch_exons_for_region(conn, sr, coord_system_name)
        if not df.empty:
            dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def cds_from_exons_and_translations(
    exons: pd.DataFrame, translations: pd.DataFrame
) -> pd.DataFrame:
    columns = [
        "transcript_id",
        "translation_id",
        "exon_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "cds_rank",
    ]
    if exons.empty or translations.empty:
        return pd.DataFrame(columns=columns)

    exon_by_tx = {
        tx_id: group.sort_values("exon_rank")
        for tx_id, group in exons.groupby("transcript_id", dropna=False)
    }
    rows: list[dict] = []
    for _, tr in translations.iterrows():
        tx_id = tr.get("transcript_id")
        tx_exons = exon_by_tx.get(tx_id)
        if tx_exons is None or tx_exons.empty:
            continue
        start_exon_id = tr.get("start_exon_id")
        end_exon_id = tr.get("end_exon_id")
        started = False
        cds_rank = 0
        for _, exon in tx_exons.iterrows():
            exon_id = exon.get("exon_id")
            if exon_id == start_exon_id:
                started = True
            if not started:
                continue
            start = int(exon.get("seq_region_start"))
            end = int(exon.get("seq_region_end"))
            strand = int(exon.get("seq_region_strand") or 0)
            if exon_id == start_exon_id:
                offset = max(0, int(tr.get("seq_start") or 1) - 1)
                if strand == -1:
                    end -= offset
                else:
                    start += offset
            if exon_id == end_exon_id:
                offset = max(1, int(tr.get("seq_end") or (end - start + 1)))
                if strand == -1:
                    start = int(exon.get("seq_region_end")) - offset + 1
                else:
                    end = int(exon.get("seq_region_start")) + offset - 1
            cds_rank += 1
            rows.append(
                {
                    "transcript_id": tx_id,
                    "translation_id": tr.get("translation_id"),
                    "exon_id": exon_id,
                    "seq_region_name": exon.get("seq_region_name"),
                    "seq_region_start": start,
                    "seq_region_end": end,
                    "seq_region_strand": exon.get("seq_region_strand"),
                    "cds_rank": cds_rank,
                }
            )
            if exon_id == end_exon_id:
                break
    return pd.DataFrame(rows, columns=columns)

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
    conn, coord_system_name: Optional[str] = None, allowlist: Optional[Iterable[str]] = None
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
    return f'''
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
'''


def fetch_genes_for_region(
    conn, seq_region_name: str, coord_system_name: Optional[str]
) -> pd.DataFrame:
    sql = _gene_query(coord_system_name)
    params = (coord_system_name, seq_region_name) if coord_system_name else (seq_region_name,)
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
    return f'''
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
'''


def fetch_transcripts_for_region(
    conn, seq_region_name: str, coord_system_name: Optional[str]
) -> pd.DataFrame:
    sql = _transcript_query(coord_system_name)
    params = (coord_system_name, seq_region_name) if coord_system_name else (seq_region_name,)
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

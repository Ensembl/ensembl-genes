"""
Build the dashboard dataset (SQLite) from comparator results.

The dataset is an index of comparator outputs; it computes no metric.

Contract results (output contract 1.0, benchmark/contract_reader.py):

    metrics.tsv                    every row verbatim (scope, metric, numerator,
                                   denominator, value, status, reason, populations,
                                   unit, side) -> ``metrics``
    comparison_details.tsv         per-gene rows -> ``ref_outcomes`` / ``query_outcomes``
    cds_exact_one_to_one_pairs.tsv -> ``exact_pairs``
    splice / intron-length tables  -> ``splice_sites`` / ``intron_lengths`` (per scope)
    per_sequence.tsv               -> ``per_sequence`` (drill-down only)
    stop_codon_harmonisation.tsv   -> ``stop_codon_audit``
    evaluated_transcripts.tsv.gz   -> ``tx_status`` and, with the input files, the
                                      locus models (benchmark/models.py)

Pre-contract results go through the explicit legacy importer (legacy.py):
identical-definition summary fields only, every other contract metric
not_available, best-pair estimates under historical ids. They have no locus
models (no evaluated-transcript export) and no scoped or splice values.

The build is incremental: each component (one result, one reference's models,
one run's models) has a key derived from its inputs and this code; unchanged
components are skipped.
"""

from __future__ import annotations

import gzip
import json
import re
import sqlite3
import time
from pathlib import Path

import pandas as pd

from ensembl.genes.annotation_qc.benchmark import contract_reader as cr
from ensembl.genes.annotation_qc.benchmark import legacy
from ensembl.genes.annotation_qc.benchmark.config import Experiment
from ensembl.genes.annotation_qc.benchmark.models import (
    ModelsUnavailable,
    evaluated_status,
    parse_input,
    rebuild,
)
from ensembl.genes.annotation_qc.benchmark.preparation import is_gzip
from ensembl.genes.annotation_qc.benchmark.provenance import (
    now,
    sha256_file,
    stable_hash,
)
from ensembl.genes.annotation_qc.benchmark.results import (
    DISPLAYED,
    STATUS_LEGACY,
)
from ensembl.genes.annotation_qc.benchmark.runs import load_registry, plan
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace, read_json
from ensembl.genes.annotation_qc.metrics.pairwise import contract
from ensembl.genes.annotation_qc.parsers.seqnames import build_seqname_mapping

# 4: contract 1.0 ingestion; no metric derivation; exact evaluated models.
STORE_VERSION = "4"
# Content of the code that turns comparator outputs into dataset rows.
BUILD_CODE = stable_hash(
    {
        name: sha256_file(Path(__file__).with_name(name))
        for name in ("store.py", "legacy.py", "models.py", "contract_reader.py")
    }
)
METRIC_DICTIONARY = Path(__file__).with_name("metric_dictionary.json")
ALL = contract.ALL_SCOPE
# Reserved colours for runs with a display role (Okabe-Ito vermilion and blue: distinct
# for common colour-vision deficiencies); the palette for the other tools avoids both.
ROLE_COLORS = {"focus": "#d55e00", "primary_comparator": "#0072b2"}
PALETTE = [
    "#1d8a5a",
    "#b5407a",
    "#5b4bb7",
    "#a8780a",
    "#2e8a99",
    "#7a5230",
    "#6b7a86",
    "#8f3b3b",
]

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS components (component TEXT PRIMARY KEY, key TEXT, built TEXT, info TEXT);
CREATE TABLE IF NOT EXISTS results (reference_id TEXT, run_id TEXT, policy_id TEXT, status TEXT, reason TEXT,
    source TEXT, record TEXT);
CREATE TABLE IF NOT EXISTS metrics (reference_id TEXT, run_id TEXT, policy_id TEXT, scope_id TEXT, metric_id TEXT,
    numerator INTEGER, denominator INTEGER, value REAL, value_text TEXT, status TEXT, reason TEXT, numerator_id TEXT,
    denominator_id TEXT, unit TEXT, side TEXT, source TEXT);
CREATE TABLE IF NOT EXISTS ref_outcomes (reference_id TEXT, policy_id TEXT, run_id TEXT, gene_id TEXT, chrom TEXT,
    start INTEGER, end INTEGER, strand TEXT, biotype TEXT, classification TEXT, classification_cds TEXT,
    cds_outcome TEXT, cds_exact INTEGER, cds_overlap REAL, cds_shared_bases INTEGER, cds_locus_recovered INTEGER,
    exon_overlap REAL, cds_intron_chain_match TEXT, cds_intron_chain_any_pair TEXT, intron_chain_match TEXT,
    intron_chain_any_pair TEXT, exon_coordinate_exact INTEGER, cds_matched_id TEXT, best_cds_ref_tx TEXT,
    best_cds_query_tx TEXT, matched_id TEXT, best_match_ref_tx TEXT, best_match_query_tx TEXT,
    counterpart_count INTEGER, counterpart_ids TEXT, strand_mismatch_basis TEXT, evaluated_transcript_count INTEGER,
    max_exon_count INTEGER, max_cds_count INTEGER, has_cds INTEGER, one_to_one_partner TEXT,
    cds_start_diff_bp INTEGER, cds_stop_diff_bp INTEGER, cds_introns_missing INTEGER, cds_introns_extra INTEGER,
    scope_ids TEXT);
CREATE TABLE IF NOT EXISTS query_outcomes (reference_id TEXT, policy_id TEXT, run_id TEXT, gene_id TEXT,
    original_gene_id TEXT, chrom TEXT, start INTEGER, end INTEGER, strand TEXT, classification TEXT,
    classification_cds TEXT, cds_exact INTEGER, cds_overlap REAL, cds_shared_bases INTEGER, cds_matched_id TEXT,
    best_cds_ref_tx TEXT, best_cds_query_tx TEXT, novel_category TEXT, counterpart_count INTEGER,
    counterpart_ids TEXT, strand_mismatch_basis TEXT, evaluated_transcript_count INTEGER, max_cds_count INTEGER,
    has_cds INTEGER, one_to_one_partner TEXT, scope_ids TEXT, search TEXT);
CREATE TABLE IF NOT EXISTS ref_gene_policy (reference_id TEXT, policy_id TEXT, gene_id TEXT, chrom TEXT,
    start INTEGER, end INTEGER, strand TEXT, biotype TEXT, evaluated_transcript_count INTEGER,
    max_cds_count INTEGER, max_exon_count INTEGER, has_cds INTEGER, scope_ids TEXT, search TEXT);
CREATE TABLE IF NOT EXISTS exact_pairs (reference_id TEXT, policy_id TEXT, run_id TEXT, reference_gene_id TEXT,
    query_gene_id TEXT, query_original_gene_id TEXT, chrom TEXT, strand TEXT, cds_start INTEGER, cds_end INTEGER,
    cds_segments INTEGER, supporting_reference_transcript_ids TEXT, supporting_query_transcript_ids TEXT,
    reference_exact_candidates INTEGER, query_exact_candidates INTEGER, component_reference_genes INTEGER,
    component_query_genes INTEGER);
CREATE TABLE IF NOT EXISTS splice_sites (reference_id TEXT, scope_id TEXT, run_id TEXT, policy_id TEXT,
    min_intron_length INTEGER, coordinate TEXT, kind TEXT, bin INTEGER, lo REAL, hi REAL, assessed INTEGER,
    supported INTEGER, missing_locus INTEGER);
CREATE TABLE IF NOT EXISTS intron_lengths (reference_id TEXT, scope_id TEXT, run_id TEXT, policy_id TEXT,
    bin INTEGER, lo INTEGER, hi INTEGER, label TEXT, eligible INTEGER, exact INTEGER, donor INTEGER,
    acceptor INTEGER, missing_locus INTEGER);
CREATE TABLE IF NOT EXISTS per_sequence (reference_id TEXT, run_id TEXT, policy_id TEXT, chrom TEXT,
    reference_genes INTEGER, reference_cds_locus_recovered INTEGER, reference_cds_coordinate_exact INTEGER,
    query_genes INTEGER, query_cds_coordinate_exact INTEGER, one_to_one_tp INTEGER);
CREATE TABLE IF NOT EXISTS stop_codon_audit (run_id TEXT, transcript_id TEXT, original_transcript_id TEXT,
    gene_id TEXT, chrom TEXT, strand TEXT, genetic_code TEXT, genetic_code_source TEXT, status TEXT, reason TEXT,
    last_codon TEXT, next_codon TEXT, exon_representation TEXT, terminal_cds_start INTEGER,
    terminal_cds_end INTEGER, new_terminal_cds_start INTEGER, new_terminal_cds_end INTEGER);
CREATE TABLE IF NOT EXISTS ref_genes (reference_id TEXT, gene_id TEXT, original_gene_id TEXT, name TEXT,
    chrom TEXT, start INTEGER, end INTEGER, strand TEXT, biotype TEXT, source_feature TEXT, evaluated INTEGER,
    search TEXT);
CREATE TABLE IF NOT EXISTS query_genes (run_id TEXT, gene_id TEXT, original_gene_id TEXT, name TEXT, chrom TEXT,
    start INTEGER, end INTEGER, strand TEXT, search TEXT);
CREATE TABLE IF NOT EXISTS transcripts (source TEXT, gene_id TEXT, transcript_id TEXT, original_transcript_id TEXT,
    chrom TEXT, start INTEGER, end INTEGER, strand TEXT, biotype TEXT, tags TEXT, exons TEXT, cds TEXT,
    transformed TEXT);
CREATE TABLE IF NOT EXISTS tx_status (source TEXT, policy_id TEXT, transcript_id TEXT, gene_id TEXT, status TEXT,
    exclusion_reason TEXT, selected_by TEXT, exon_count INTEGER, cds_segment_count INTEGER, cds_length INTEGER,
    transformed TEXT);
CREATE TABLE IF NOT EXISTS seq_maxlen (source TEXT, chrom TEXT, max_len INTEGER);
CREATE INDEX IF NOT EXISTS ix_results ON results (reference_id, run_id, policy_id);
CREATE INDEX IF NOT EXISTS ix_metrics ON metrics (reference_id, policy_id, scope_id, run_id);
CREATE INDEX IF NOT EXISTS ix_ro ON ref_outcomes (reference_id, policy_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_ro_run ON ref_outcomes (reference_id, policy_id, run_id, cds_outcome);
CREATE INDEX IF NOT EXISTS ix_qo ON query_outcomes (reference_id, policy_id, run_id, chrom, start);
CREATE INDEX IF NOT EXISTS ix_qo_id ON query_outcomes (run_id, policy_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_rgp ON ref_gene_policy (reference_id, policy_id, chrom, start);
CREATE INDEX IF NOT EXISTS ix_rgp_id ON ref_gene_policy (reference_id, policy_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_pairs ON exact_pairs (reference_id, policy_id, run_id, reference_gene_id);
CREATE INDEX IF NOT EXISTS ix_splice ON splice_sites (reference_id, scope_id, policy_id, run_id);
CREATE INDEX IF NOT EXISTS ix_lengths ON intron_lengths (reference_id, scope_id, policy_id, run_id);
CREATE INDEX IF NOT EXISTS ix_perseq ON per_sequence (reference_id, policy_id, run_id);
CREATE INDEX IF NOT EXISTS ix_stop ON stop_codon_audit (run_id, transcript_id);
CREATE INDEX IF NOT EXISTS ix_ref_genes ON ref_genes (reference_id, chrom, start);
CREATE INDEX IF NOT EXISTS ix_ref_genes_id ON ref_genes (reference_id, gene_id);
CREATE INDEX IF NOT EXISTS ix_qg ON query_genes (run_id, chrom, start);
CREATE INDEX IF NOT EXISTS ix_tx ON transcripts (source, chrom, start);
CREATE INDEX IF NOT EXISTS ix_tx_id ON transcripts (source, transcript_id);
CREATE INDEX IF NOT EXISTS ix_txs ON tx_status (source, policy_id, transcript_id);
CREATE INDEX IF NOT EXISTS ix_txs_gene ON tx_status (source, policy_id, gene_id);
"""

_PREFIX = re.compile(r"^(?:gene|transcript|chromosome|mRNA):")
_GTF_KV = re.compile(r'(\S+)\s+"([^"]*)"')
_TRUE = {"True": 1, "False": 0, "true": 1, "false": 0}


def connect(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("PRAGMA journal_mode=WAL")
    row = None
    try:
        row = con.execute("SELECT value FROM meta WHERE key='store_version'").fetchone()
    except sqlite3.OperationalError:
        pass
    if row is not None and json.loads(row[0]) != STORE_VERSION:
        raise RuntimeError(
            f"{path} was built with store version {json.loads(row[0])}; this code writes version "
            f"{STORE_VERSION}. Build into a new workspace (output_dir) instead of overwriting it."
        )
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
        (component, key, now(), json.dumps(info or {}, default=str)),
    )


def gene_names(path: Path) -> dict[tuple[str, str], str]:
    """(seqname, file identity) -> display name from GFF3 Name / GTF gene_name. Labels only."""
    names = {}
    opener = gzip.open if is_gzip(path) else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
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


def _int(v):
    if v in ("", "NA", None):
        return None
    return int(float(v))


def _float(v):
    if v in ("", "NA", None):
        return None
    return float(v)


def _bool(v):
    return _TRUE.get(v) if v not in ("", "NA", None) else None


def _col(frame: pd.DataFrame, name: str) -> pd.Series:
    if name in frame.columns:
        return frame[name]
    return pd.Series([""] * len(frame), index=frame.index, dtype=object)


REF_DETAIL = (
    "gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "gene_biotype",
    "classification",
    "classification_cds",
    "cds_outcome",
    "cds_coordinate_exact",
    "cds_overlap",
    "cds_shared_bases",
    "cds_locus_recovered",
    "exon_overlap",
    "cds_intron_chain_match",
    "cds_intron_chain_any_pair",
    "intron_chain_match",
    "intron_chain_any_pair",
    "exon_coordinate_exact",
    "cds_matched_id",
    "best_cds_match_own_transcript_id",
    "best_cds_match_transcript_id",
    "matched_id",
    "best_match_own_transcript_id",
    "best_match_transcript_id",
    "counterpart_count",
    "counterpart_ids",
    "strand_mismatch_basis",
    "evaluated_transcript_count",
    "max_exon_count",
    "max_cds_count",
    "has_cds",
    "one_to_one_partner",
    "cds_start_diff_bp",
    "cds_stop_diff_bp",
    "cds_introns_missing",
    "cds_introns_extra",
    "scope_ids",
)
QUERY_DETAIL = (
    "gene_id",
    "original_gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "classification",
    "classification_cds",
    "cds_coordinate_exact",
    "cds_overlap",
    "cds_shared_bases",
    "cds_matched_id",
    "best_cds_match_transcript_id",
    "best_cds_match_own_transcript_id",
    "novel_category",
    "counterpart_count",
    "counterpart_ids",
    "strand_mismatch_basis",
    "evaluated_transcript_count",
    "max_cds_count",
    "has_cds",
    "one_to_one_partner",
    "scope_ids",
)
INT_FIELDS = {
    "start",
    "end",
    "cds_shared_bases",
    "counterpart_count",
    "evaluated_transcript_count",
    "max_exon_count",
    "max_cds_count",
    "cds_start_diff_bp",
    "cds_stop_diff_bp",
    "cds_introns_missing",
    "cds_introns_extra",
}
BOOL_FIELDS = {
    "cds_coordinate_exact",
    "cds_locus_recovered",
    "exon_coordinate_exact",
    "has_cds",
}
FLOAT_FIELDS = {"cds_overlap", "exon_overlap"}
KEEP_NA = {"cds_intron_chain_match", "intron_chain_match"}


def _typed(name: str, value):
    if name in INT_FIELDS:
        return _int(value)
    if name in BOOL_FIELDS:
        return _bool(value)
    if name in FLOAT_FIELDS:
        return _float(value)
    if name in KEEP_NA:
        return value or "NA"
    return value if value not in ("",) else None


def _outcome_rows(det: pd.DataFrame, where: tuple) -> tuple[list, list]:
    """Per-gene rows as written by the comparator (typed; absent columns -> NULL)."""
    ref = det[det.source == "reference"]
    qry = det[det.source == "consensus"]
    ref_cols = [_col(ref, c) for c in REF_DETAIL]
    q_cols = [_col(qry, c) for c in QUERY_DETAIL]
    ref_rows = [
        (*where, *(_typed(n, v) for n, v in zip(REF_DETAIL, values)))
        for values in zip(*ref_cols)
    ]
    q_rows = []
    for values in zip(*q_cols):
        row = [_typed(n, v) for n, v in zip(QUERY_DETAIL, values)]
        q_rows.append((*where, *row, " ".join(filter(None, values[:2])).lower()))
    return ref_rows, q_rows


PAIR_COLUMNS = [
    "reference_gene_id",
    "query_gene_id",
    "query_original_gene_id",
    "chrom",
    "strand",
    "cds_start",
    "cds_end",
    "cds_segments",
    "supporting_reference_transcript_ids",
    "supporting_query_transcript_ids",
    "reference_exact_candidates",
    "query_exact_candidates",
    "component_reference_genes",
    "component_query_genes",
]
AUDIT_COLUMNS = (
    "transcript_id",
    "original_transcript_id",
    "gene_id",
    "chrom",
    "strand",
    "genetic_code",
    "genetic_code_source",
    "status",
    "reason",
    "last_codon",
    "next_codon",
    "exon_representation",
    "terminal_cds_start",
    "terminal_cds_end",
    "new_terminal_cds_start",
    "new_terminal_cds_end",
)
RESULT_TABLES = (
    "results",
    "metrics",
    "ref_outcomes",
    "query_outcomes",
    "exact_pairs",
    "splice_sites",
    "intron_lengths",
    "per_sequence",
)


def _limitations(summary: dict, policy, run) -> list[dict]:
    """Limitations that the comparator's own outputs establish for this result."""
    notes = []
    outside = ((summary.get("gene_span_checks") or {}).get("query") or {}).get(
        "transcripts_outside"
    ) or {}
    if outside.get("transcripts_outside_gene_span"):
        notes.append(
            {
                "id": "query_transcripts_outside_gene_records",
                "affects": ["cds_coordinate_exact_query", "locus pairing"],
                "text": (
                    f"{outside['transcripts_outside_gene_span']} query transcripts in "
                    f"{outside['genes_with_transcripts_outside_span']} genes lie outside their own "
                    f"gene record (by up to {outside.get('max_outside_bp', 0):,} bp). Loci are paired "
                    "by gene-record spans, so such a transcript can share an identical CDS with a "
                    "reference gene its gene is never paired with: the span-paired per-gene exact "
                    "count can be lower than an independent identical-CDS count. One-to-one TP, "
                    "recall, precision and F1 use the identical-CDS matching graph and are not affected."
                ),
            }
        )
    if (
        run.variant == "stop_harmonised"
        and policy.reference_transcript_selection == "all"
    ):
        notes.append(
            {
                "id": "stop_harmonised_best_pair_chain",
                "affects": ["cds_intron_chain_best_pair"],
                "text": (
                    "Supplementary stop-harmonised variant under policy A: best-pair CDS "
                    "intron-chain counts can differ slightly from the original-coordinate run "
                    "because the extended terminal CDS changes best-pair ranking (0.8 reciprocal "
                    "overlap). Any-pair chains, unique introns, R and Q are unchanged."
                ),
            }
        )
    return notes


def _pairs_check(pairs: pd.DataFrame | None, metric_rows: list[dict]) -> dict:
    tp = next(
        (
            r["numerator"]
            for r in metric_rows
            if r["scope_id"] == ALL and r["metric_id"] == "exact_cds_one_to_one_tp"
        ),
        None,
    )
    if pairs is None:
        return {"status": "absent", "file": cr.PAIRS, "true_positives": tp}
    ok = (
        tp is not None
        and len(pairs) == tp
        and pairs.reference_gene_id.is_unique
        and pairs.query_gene_id.is_unique
    )
    return {
        "status": "verified" if ok else "disagrees",
        "file": cr.PAIRS,
        "rows": len(pairs),
        "true_positives": tp,
    }


def _independent_check(
    ws: Workspace, run_id: str, policy_id: str, det: pd.DataFrame | None
) -> dict:
    """Independent any-pair table against the comparator's per-gene rows (verification only)."""
    computed = ws.independent_dir(run_id) / f"{policy_id}.reference_genes.tsv"
    imported = read_json(ws.independent_dir(run_id) / f"{policy_id}.import.json")
    path = (
        computed
        if computed.exists()
        else Path(imported["source_path"]) if imported else None
    )
    if path is None or det is None:
        return {"status": "absent"}
    if not path.exists():
        return {"status": "missing", "path": str(path)}
    ind = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    ind = ind.set_index("gene_id")
    ref = det[det.source == "reference"].set_index("gene_id")
    if set(ind.index) != set(ref.index):
        return {"status": "gene_sets_differ", "path": str(path)}
    out: dict = {"status": "compared", "path": str(path), "genes": len(ref)}
    if (
        "cds_intron_chain_any_pair" in ref.columns
        and "cds_intron_chain_any_pair" in ind
    ):
        a = ref.cds_intron_chain_any_pair.replace({"NA": "False"})
        b = ind.cds_intron_chain_any_pair.reindex(ref.index).replace({"": "False"})
        out["any_pair_chain_disagreements"] = int((a != b).sum())
    return out


def _contract_tables(out_dir: Path, ids: tuple) -> tuple[list, list, list, list]:
    ref_id, run_id, policy_id = ids
    position = cr.optional(out_dir, cr.SPLICE_POSITION)
    length = cr.optional(out_dir, cr.INTRON_LENGTH)
    per_seq = cr.optional(out_dir, cr.PER_SEQUENCE)
    audit = cr.optional(out_dir, cr.STOP_AUDIT)
    site_rows = (
        []
        if position is None
        else [
            (
                ref_id,
                r.scope_id,
                run_id,
                policy_id,
                _int(r.min_intron_length),
                r.coordinate,
                r.kind,
                _int(r.bin),
                _float(r.lo),
                _float(r.hi),
                _int(r.assessed),
                _int(r.supported),
                _int(r.missing_locus),
            )
            for r in position.itertuples(index=False)
        ]
    )
    length_rows = (
        []
        if length is None
        else [
            (
                ref_id,
                r.scope_id,
                run_id,
                policy_id,
                _int(r.bin),
                _int(r.lo_bp),
                _int(r.hi_bp),
                r.label,
                _int(r.eligible),
                _int(r.exact),
                _int(r.donor),
                _int(r.acceptor),
                _int(r.missing_locus),
            )
            for r in length.itertuples(index=False)
        ]
    )
    seq_cols = [
        "chrom",
        "reference_genes",
        "reference_cds_locus_recovered",
        "reference_cds_coordinate_exact",
        "query_genes",
        "query_cds_coordinate_exact",
        "one_to_one_tp",
    ]
    seq_rows = (
        []
        if per_seq is None
        else [
            (ref_id, run_id, policy_id, row[0], *(_int(x) for x in row[1:]))
            for row in per_seq[seq_cols].values.tolist()
        ]
    )
    audit_rows = (
        []
        if audit is None
        else [
            (run_id, *row)
            for row in audit[list(AUDIT_COLUMNS)]
            .replace({"": None})
            .astype(object)
            .values.tolist()
        ]
    )
    return site_rows, length_rows, seq_rows, audit_rows


def _public_record(record, manifest, summary, task) -> dict:
    manifest = manifest or {}
    options = manifest.get("options") or {}
    public = {
        k: v
        for k, v in (record or {}).items()
        if k not in ("environment", "key_components")
    }
    public.update(
        assessment_reason=task.reason,
        contract_version=manifest.get("contract_version"),
        definition_version=manifest.get("definition_version"),
        analyses=manifest.get("analyses"),
        transformations=manifest.get("transformations"),
        provenance=manifest.get("provenance"),
        thresholds=manifest.get("thresholds"),
        scope=manifest.get("scope"),
        inputs={
            k: {kk: v.get(kk) for kk in ("path", "sha256", "compression")}
            for k, v in (manifest.get("inputs") or {}).items()
            if isinstance(v, dict)
        },
        options={
            k: options.get(k)
            for k in (
                "add_stop_codon",
                "stop_codon_genetic_code",
                "stop_codon_sequence_code",
                "reference_transcript_selection",
                "query_transcript_selection",
                "evaluation_mode",
                "reference_transcript_biotypes",
                "genome_mismatch",
                "seqname_map",
                "analysis_scope",
            )
        },
        stop_codon_harmonisation=summary.get("query_stop_codon_harmonisation"),
        gene_span_checks=summary.get("gene_span_checks"),
        limitations=_limitations(summary, task.policy, task.run),
        output_files=sorted((manifest.get("outputs") or {}).keys()),
    )
    if task.assessment:
        public["identity"] = task.assessment.get("identity")
        public["outputs_check"] = task.assessment.get("outputs_check")
    return public


def build_result(con, ws: Workspace, task, log=print) -> None:
    """Index one (run, policy) result: contract outputs verbatim or the legacy importer."""
    run, policy, ref = task.run, task.policy, task.reference
    record = ws.run_record(run.run_id, policy.id) if task.status in DISPLAYED else None
    out_dir = Path(record["output_dir"]) if record else None
    manifest = read_json(out_dir / cr.MANIFEST) if out_dir else None
    key = stable_hash(
        {
            "v": STORE_VERSION,
            "code": BUILD_CODE,
            "status": task.status,
            "reason": task.reason,
            "output_dir": str(out_dir),
            "outputs": (manifest or {}).get("outputs"),
            "created": (manifest or {}).get("created"),
            "identity": (task.assessment or {}).get("identity"),
        }
    )
    component = f"result:{run.run_id}:{policy.id}"
    if _component_key(con, component) == key:
        return
    where = (ref.id, run.run_id, policy.id)
    with con:
        for table in RESULT_TABLES:
            con.execute(
                f"DELETE FROM {table} WHERE reference_id=? AND run_id=? AND policy_id=?",
                where,
            )
    attempts = ws.attempts(run.run_id, policy.id)
    if record is None:
        public = {"assessment_reason": task.reason}
        if attempts:
            public["attempts"] = attempts[-5:]
        with con:
            con.execute(
                "INSERT INTO results VALUES (?,?,?,?,?,?,?)",
                (
                    *where,
                    task.status,
                    task.reason,
                    None,
                    json.dumps(public, default=str),
                ),
            )
            _set_component(con, component, key)
        return
    assert out_dir is not None  # a displayed result always has its output directory
    summary = read_json(out_dir / "comparison_summary.json") or {}
    pairs = cr.optional(out_dir, cr.PAIRS)
    public = _public_record(record, manifest, summary, task)
    site_rows: list = []
    length_rows: list = []
    seq_rows: list = []
    audit_rows: list = []
    if task.status == STATUS_LEGACY:
        rows, det = legacy.legacy_result(out_dir)
        public["legacy"] = True
    else:
        # contract results: metrics.tsv rows only (no dashboard-side estimates; the
        # historical best-pair rows belong to the legacy importer)
        rows = cr.metrics(out_dir)
        det = cr.details(out_dir)
        site_rows, length_rows, seq_rows, audit_rows = _contract_tables(out_dir, where)
    public["pairs_check"] = _pairs_check(pairs, rows)
    public["independent"] = _independent_check(ws, run.run_id, policy.id, det)
    keyed = (
        ref.id,
        policy.id,
        run.run_id,
    )  # column order of the per-gene and pair tables
    ref_rows, q_rows = _outcome_rows(det, keyed) if det is not None else ([], [])
    pair_rows = (
        []
        if pairs is None
        else [
            (*keyed, *row)
            for row in pairs.reindex(columns=PAIR_COLUMNS).fillna("").values.tolist()
        ]
    )
    with con:
        con.execute(
            "INSERT INTO results VALUES (?,?,?,?,?,?,?)",
            (
                *where,
                task.status,
                task.reason,
                record.get("source"),
                json.dumps(public, default=str),
            ),
        )
        con.executemany(
            f"INSERT INTO metrics VALUES ({','.join('?' * 16)})",
            [
                (
                    *where,
                    r["scope_id"],
                    r["metric_id"],
                    r["numerator"],
                    r["denominator"],
                    r["value"],
                    r.get(
                        "value_text", "NA" if r["value"] is None else repr(r["value"])
                    ),
                    r["status"],
                    r["reason"],
                    r["numerator_id"],
                    r["denominator_id"],
                    r["unit"],
                    r["side"],
                    r["source"],
                )
                for r in rows
            ],
        )
        con.executemany(
            f"INSERT INTO ref_outcomes VALUES ({','.join('?' * (3 + len(REF_DETAIL)))})",
            ref_rows,
        )
        con.executemany(
            f"INSERT INTO query_outcomes VALUES ({','.join('?' * (4 + len(QUERY_DETAIL)))})",
            q_rows,
        )
        con.executemany(
            f"INSERT INTO exact_pairs VALUES ({','.join('?' * 17)})", pair_rows
        )
        con.executemany(
            f"INSERT INTO splice_sites VALUES ({','.join('?' * 13)})", site_rows
        )
        con.executemany(
            f"INSERT INTO intron_lengths VALUES ({','.join('?' * 13)})", length_rows
        )
        con.executemany(
            f"INSERT INTO per_sequence VALUES ({','.join('?' * 10)})", seq_rows
        )
        if audit_rows:
            con.execute("DELETE FROM stop_codon_audit WHERE run_id=?", (run.run_id,))
            con.executemany(
                f"INSERT INTO stop_codon_audit VALUES ({','.join('?' * 17)})",
                audit_rows,
            )
        _set_component(con, component, key, {"metric_rows": len(rows)})
    log(f"  result {run.run_id}/{policy.id}: {task.status} ({len(rows)} metric rows)")


# ----------------------------------------------------------------------------- models


def _contract_results(ws: Workspace, tasks, select) -> list[tuple]:
    """(task, out_dir, manifest) of displayed contract results matching ``select``."""
    out = []
    for task in tasks:
        if task.status not in DISPLAYED or task.status == STATUS_LEGACY:
            continue
        if not select(task):
            continue
        record = ws.run_record(task.run.run_id, task.policy.id) or {}
        out_dir = Path(record.get("output_dir", ""))
        manifest = read_json(out_dir / cr.MANIFEST)
        if manifest:
            out.append((task, out_dir, manifest))
    return out


def _export(out_dir: Path, role: str) -> pd.DataFrame:
    frame = cr.evaluated(out_dir)
    return frame[frame.role == role].reset_index(drop=True)


def _clear_models(con, source: str, gene_table: str, gene_key: str, owner: str) -> None:
    con.execute("DELETE FROM transcripts WHERE source=?", (source,))
    con.execute("DELETE FROM tx_status WHERE source=?", (source,))
    con.execute(
        "DELETE FROM seq_maxlen WHERE source IN (?,?)", (source, f"genes:{source}")
    )
    con.execute(f"DELETE FROM {gene_table} WHERE {gene_key}=?", (owner,))


def _store_models(con, source, models: pd.DataFrame, statuses: dict) -> None:
    columns = [
        "gene_id",
        "transcript_id",
        "original_transcript_id",
        "chrom",
        "start",
        "end",
        "strand",
        "biotype",
        "tags",
        "exons",
        "cds",
        "transformed",
    ]
    con.executemany(
        f"INSERT INTO transcripts VALUES ({','.join('?' * 13)})",
        [(source, *row) for row in models[columns].astype(object).values.tolist()],
    )
    for policy_id, status in statuses.items():
        con.executemany(
            f"INSERT INTO tx_status VALUES ({','.join('?' * 11)})",
            [
                (
                    source,
                    policy_id,
                    t,
                    g,
                    s,
                    er or None,
                    sb or None,
                    int(ec),
                    int(cc),
                    int(cl),
                    tr or None,
                )
                for t, g, s, er, sb, ec, cc, cl, tr in status.values.tolist()
            ],
        )
    con.execute(
        "INSERT INTO seq_maxlen SELECT source, chrom, MAX(end-start+1) FROM transcripts "
        "WHERE source=? GROUP BY chrom",
        (source,),
    )


def _gene_rows(parsed, names, owner, evaluated, tx_ids) -> list[tuple]:
    genes = parsed[parsed.Feature == "gene"].drop_duplicates("gene_id")
    rows: list[tuple] = []
    for g in genes.itertuples(index=False):
        name = names.get((g.Chromosome, g.original_gene_id), "")
        search = " ".join(
            filter(
                None, [g.gene_id, g.original_gene_id, name, tx_ids.get(g.gene_id, "")]
            )
        ).lower()
        base = (
            owner,
            g.gene_id,
            g.original_gene_id,
            name,
            g.Chromosome,
            int(g.Start) + 1,
            int(g.End),
            g.Strand,
        )
        if evaluated is None:
            rows.append((*base, search))
        else:
            rows.append(
                (
                    *base,
                    g.gene_biotype,
                    g.source_feature,
                    int(g.gene_id in evaluated),
                    search,
                )
            )
    return rows


def _mapped_names(manifest: dict, role: str) -> dict:
    """Gene display names keyed by comparison sequence name (labels only)."""
    inputs = manifest.get("inputs") or {}
    mapping = build_seqname_mapping(
        (inputs.get("seqname_map") or {}).get("path"),
        (inputs.get("assembly_report") or {}).get("path"),
    )
    path = Path(inputs[role]["path"])
    return {(mapping.get(c, c), i): n for (c, i), n in gene_names(path).items()}


def _rebuild_all(parsed, results, role) -> tuple:
    """Models from the first result; every result's export checked against the same parse."""
    statuses, models, checked = {}, None, []
    for task, out_dir, manifest in results:
        export = _export(out_dir, role)
        if set(export.input_sha256) != {manifest["inputs"][role]["sha256"]}:
            raise ModelsUnavailable(
                f"{out_dir}: export rows name another input checksum than the manifest"
            )
        audit = cr.optional(out_dir, cr.STOP_AUDIT) if role == "query" else None
        rebuilt = rebuild(parsed, export, audit)
        if models is None:
            models = rebuilt
        else:
            a = rebuilt.set_index("transcript_id")[["exons", "cds"]].sort_index()
            b = models.set_index("transcript_id")[["exons", "cds"]].sort_index()
            if not a.equals(b):
                raise ModelsUnavailable(
                    f"policy {task.policy.id} evaluated other coordinates than policy "
                    f"{results[0][0].policy.id} for the same {role} input"
                )
        statuses[task.policy.id] = evaluated_status(export)
        checked.append(task.policy.id)
    return models, statuses, checked


def build_reference_models(
    con, ws: Workspace, ref, tasks, model_errors, log=print
) -> bool:
    results = _contract_results(ws, tasks, lambda t: t.reference.id == ref.id)
    source = f"ref:{ref.id}"
    first: dict = {}
    for (
        item
    ) in results:  # one result per policy: the reference side is the same for every run
        first.setdefault(item[0].policy.id, item)
    chosen = list(first.values())
    if not chosen:
        model_errors[source] = (
            "no contract result with an evaluated-transcript export for this reference "
            "(pre-contract results cannot be rebuilt without selecting transcripts again)"
        )
        with con:
            _clear_models(con, source, "ref_genes", "reference_id", ref.id)
            con.execute(
                "DELETE FROM components WHERE component=?", (f"reference:{ref.id}",)
            )
        return False
    manifest0 = chosen[0][2]
    key = stable_hash(
        {
            "v": STORE_VERSION,
            "code": BUILD_CODE,
            "input": manifest0["inputs"]["reference"]["sha256"],
            "seqmap": (manifest0["inputs"].get("seqname_map") or {}).get("sha256"),
            "exports": {
                t.policy.id: (m.get("outputs") or {})
                .get(cr.EVALUATED, {})
                .get("sha256")
                for t, _, m in chosen
            },
        }
    )
    component = f"reference:{ref.id}"
    if _component_key(con, component) == key:
        return True
    t0 = time.time()
    try:
        parsed = parse_input(manifest0, "reference", ws.checksums)
        models, statuses, checked = _rebuild_all(parsed, chosen, "reference")
        names = _mapped_names(manifest0, "reference")
    except (ModelsUnavailable, ValueError, OSError, KeyError) as error:
        model_errors[source] = f"{type(error).__name__}: {error}"
        log(f"  reference {ref.id}: models unavailable ({model_errors[source][:160]})")
        with con:
            _clear_models(con, source, "ref_genes", "reference_id", ref.id)
            con.execute("DELETE FROM components WHERE component=?", (component,))
        return False
    evaluated: set = set()
    for status in statuses.values():
        evaluated |= set(status.loc[status.status == "evaluated", "gene_id"])
    tx_ids = models.groupby("gene_id").transcript_id.agg(" ".join).to_dict()
    gene_rows = _gene_rows(parsed, names, ref.id, evaluated, tx_ids)
    with con:
        _clear_models(con, source, "ref_genes", "reference_id", ref.id)
        con.executemany(
            f"INSERT INTO ref_genes VALUES ({','.join('?' * 12)})", gene_rows
        )
        con.execute(
            "INSERT INTO seq_maxlen SELECT ?, chrom, MAX(end-start+1) FROM ref_genes "
            "WHERE reference_id=? GROUP BY chrom",
            (f"genes:{source}", ref.id),
        )
        _store_models(con, source, models, statuses)
        _set_component(
            con,
            component,
            key,
            {
                "transcripts": len(models),
                "genes": len(gene_rows),
                "policies_checked": checked,
                "seconds": round(time.time() - t0, 1),
            },
        )
    log(
        f"  reference {ref.id}: {len(models):,} transcripts rebuilt and checked for policies "
        f"{', '.join(checked)} in {time.time() - t0:.0f} s"
    )
    return True


def build_run_models(
    con, ws: Workspace, run, tasks, model_errors, parse_cache, log=print
):
    results = _contract_results(ws, tasks, lambda t: t.run.run_id == run.run_id)
    source = f"run:{run.run_id}"
    component = f"run:{run.run_id}"
    if not results:
        model_errors[source] = (
            "no contract result with an evaluated-transcript export for this run"
        )
        with con:
            _clear_models(con, source, "query_genes", "run_id", run.run_id)
            con.execute("DELETE FROM components WHERE component=?", (component,))
        return
    manifest0 = results[0][2]
    key = stable_hash(
        {
            "v": STORE_VERSION,
            "code": BUILD_CODE,
            "input": manifest0["inputs"]["query"]["sha256"],
            "seqmap": (manifest0["inputs"].get("seqname_map") or {}).get("sha256"),
            "outputs": {
                t.policy.id: {
                    name: (m.get("outputs") or {}).get(name, {}).get("sha256")
                    for name in (cr.EVALUATED, cr.STOP_AUDIT)
                }
                for t, _, m in results
            },
        }
    )
    if _component_key(con, component) == key:
        return
    t0 = time.time()
    try:
        cache_key = (
            manifest0["inputs"]["query"]["sha256"],
            (manifest0["inputs"].get("seqname_map") or {}).get("sha256"),
            (manifest0.get("options") or {}).get("query_format"),
        )
        if parse_cache.get("key") != cache_key:
            parse_cache.clear()
            parse_cache.update(
                key=cache_key, frame=parse_input(manifest0, "query", ws.checksums)
            )
        parsed = parse_cache["frame"]
        models, statuses, checked = _rebuild_all(parsed, results, "query")
        names = _mapped_names(manifest0, "query")
    except (ModelsUnavailable, ValueError, OSError, KeyError) as error:
        model_errors[source] = f"{type(error).__name__}: {error}"
        log(f"  run {run.run_id}: models unavailable ({model_errors[source][:160]})")
        with con:
            _clear_models(con, source, "query_genes", "run_id", run.run_id)
            con.execute("DELETE FROM components WHERE component=?", (component,))
        return
    tx_ids = models.groupby("gene_id").transcript_id.agg(" ".join).to_dict()
    gene_rows = _gene_rows(parsed, names, run.run_id, None, tx_ids)
    transformed = int((models.transformed != "").sum())
    with con:
        _clear_models(con, source, "query_genes", "run_id", run.run_id)
        con.executemany(
            f"INSERT INTO query_genes VALUES ({','.join('?' * 9)})", gene_rows
        )
        _store_models(con, source, models, statuses)
        _set_component(
            con,
            component,
            key,
            {
                "transcripts": len(models),
                "transformed": transformed,
                "policies_checked": checked,
                "seconds": round(time.time() - t0, 1),
            },
        )
    log(
        f"  run {run.run_id}: {len(models):,} transcripts rebuilt ({transformed:,} with an "
        f"audited transformation) in {time.time() - t0:.0f} s"
    )


def _model_component_info(con) -> dict:
    return {
        comp: json.loads(info or "{}")
        for comp, info in con.execute(
            "SELECT component, info FROM components "
            "WHERE component LIKE 'reference:%' OR component LIKE 'run:%'"
        )
    }


def refresh_gene_policy(con) -> None:
    """One row per (reference, policy, evaluated gene), from the comparator's per-gene rows."""
    with con:
        con.execute("DELETE FROM ref_gene_policy")
        con.execute(
            "INSERT INTO ref_gene_policy SELECT o.reference_id, o.policy_id, o.gene_id, MIN(o.chrom), "
            "MIN(o.start), MIN(o.end), MIN(o.strand), MIN(o.biotype), MAX(o.evaluated_transcript_count), "
            "MAX(o.max_cds_count), MAX(o.max_exon_count), MAX(o.has_cds), MIN(o.scope_ids), "
            "LOWER(o.gene_id || ' ' || COALESCE(MIN(g.search), '')) "
            "FROM ref_outcomes o LEFT JOIN ref_genes g "
            "ON g.reference_id=o.reference_id AND g.gene_id=o.gene_id "
            "GROUP BY o.reference_id, o.policy_id, o.gene_id"
        )


def run_colors(experiment: Experiment) -> dict[str, str]:
    """
    Colour by tool; later runs of a tool and supplementary companions get lighter tints.
    A tool with a display role (focus / primary comparator) takes the role's reserved
    colour; the other tools share PALETTE in configuration order.
    """
    role_of: dict[str, str] = {}
    for run in experiment.runs.values():
        if run.display_role and run.tool not in role_of:
            role_of[run.tool] = run.display_role
    tools: list[str] = []
    for run in experiment.runs.values():
        if run.tool not in tools and run.tool not in role_of:
            tools.append(run.tool)
    seen: dict[tuple, int] = {}
    colors = {}
    for run in experiment.runs.values():
        if run.color:
            colors[run.run_id] = run.color
            continue
        base = (
            ROLE_COLORS[role_of[run.tool]]
            if run.tool in role_of
            else PALETTE[tools.index(run.tool) % len(PALETTE)]
        )
        if run.companion_of:
            colors[run.run_id] = _tint(base, 0.45)
            continue
        k = seen.get((run.reference, run.tool), 0)
        seen[(run.reference, run.tool)] = k + 1
        colors[run.run_id] = _tint(base, 0.3 * k) if k else base
    return colors


def _tint(hex_color: str, amount: float) -> str:
    amount = min(amount, 0.75)
    r, g, b = (int(hex_color[i : i + 2], 16) for i in (1, 3, 5))

    def mix(c):
        return round(c + (255 - c) * amount)

    return f"#{mix(r):02x}{mix(g):02x}{mix(b):02x}"


def _read_text(path: Path) -> str | None:
    try:
        return Path(path).read_text()
    except OSError:
        return None


def _input_audit(path: Path | None) -> dict | None:
    """Processing-scope facts from a prediction input audit (shown, never used in a value)."""
    if path is None:
        return None
    data = read_json(path)
    if data is None:
        return {"path": str(path), "available": False}
    return {
        "path": str(path),
        "available": True,
        "sequences_with_records": data.get("sequences_with_records"),
        "processed_sequences_record": data.get("processed_sequences_record"),
        "reference_pc_genes_on_sequences_without_records": data.get(
            "reference_pc_genes_on_sequences_without_records"
        ),
        "sha256": data.get("sha256"),
        "detected_format": data.get("detected_format"),
        "source_check": data.get("source_check"),
    }


def _scopes(con, ref) -> list[dict]:
    present = {
        sid
        for (sid,) in con.execute(
            "SELECT DISTINCT scope_id FROM metrics WHERE reference_id=?", (ref.id,)
        )
    }
    out = [
        {
            "id": ALL,
            "label": ref.scope.get("label") or "All compared sequences",
            "role": "primary",
            "sequences": None,
            "available": True,
            "description": ref.scope.get("note")
            or "every sequence the comparison covered",
        }
    ]
    for scope in ref.analysis_scopes:
        out.append(
            {
                **scope,
                "sequence_count": len(scope["sequences"]),
                "available": scope["id"] in present,
                "reason": (
                    None
                    if scope["id"] in present
                    else "no result of this reference reports this scope in metrics.tsv"
                ),
            }
        )
    return out


def _metric_definitions(ws: Workspace, tasks) -> dict:
    definitions: dict = {}
    for task in tasks:
        if task.status not in DISPLAYED:
            continue
        record = ws.run_record(task.run.run_id, task.policy.id) or {}
        manifest = read_json(Path(record.get("output_dir", "")) / cr.MANIFEST)
        if manifest and manifest.get("metric_definitions"):
            definitions.setdefault(
                manifest["definition_version"], manifest["metric_definitions"]
            )
    return definitions


def _run_meta(experiment: Experiment) -> list[dict]:
    colors = run_colors(experiment)
    order = list(experiment.runs)
    return [
        {
            "run_id": r.run_id,
            "reference": r.reference,
            "tool": r.tool,
            "run_label": r.run_label,
            "label": r.display_label,
            "role": r.role,
            "variant": r.variant,
            "comparator_options": r.comparator_options,
            "companion_of": r.companion_of,
            "display_role": r.display_role,
            "config_order": order.index(r.run_id),
            "query_coordinates": r.variant,
            "tool_version": r.tool_version,
            "model": r.model,
            "annotation": str(r.annotation),
            "source": r.source,
            "legacy_ids": r.legacy_ids,
            "input_audit": _input_audit(r.input_audit),
            "format": r.format,
            "provenance": r.provenance,
            "caveats": r.caveats,
            "synthetic": r.synthetic,
            "color": colors[r.run_id],
        }
        for r in experiment.runs.values()
    ]


def build_dashboard(experiment: Experiment, ws: Workspace, log=print) -> Path:
    db = ws.dashboard_db()
    con = connect(db)
    tasks = plan(experiment, ws, prepare=False)
    registry = load_registry(experiment)
    for task in tasks:
        build_result(con, ws, task, log)
    model_errors: dict[str, str] = {}
    ready = {
        ref.id: build_reference_models(con, ws, ref, tasks, model_errors, log)
        for ref in experiment.references.values()
    }
    parse_cache: dict = {}
    for run in sorted(
        experiment.runs.values(), key=lambda r: (str(r.annotation), r.run_id)
    ):
        build_run_models(con, ws, run, tasks, model_errors, parse_cache, log)
    parse_cache.clear()
    refresh_gene_policy(con)
    meta = {
        "store_version": STORE_VERSION,
        "built": now(),
        "experiment": {
            "id": experiment.id,
            "title": experiment.title,
            "description": experiment.description,
            "synthetic": experiment.synthetic,
            "config_path": str(experiment.path),
            "config_version": experiment.config_version,
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
                "caveats": list(r.raw.get("caveats") or []),
                "models_available": ready.get(r.id, False),
                "models_error": model_errors.get(f"ref:{r.id}"),
                "preparation": ws.preparation_record(r.id),
                "scopes": _scopes(con, r),
            }
            for r in experiment.references.values()
        ],
        "runs": _run_meta(experiment),
        "display": experiment.display,
        "validation": read_json(ws.validation_file()),
        "model_errors": model_errors,
        "model_components": _model_component_info(con),
        "import_report": read_json(ws.import_report()),
        "registry": {"sources": registry.sources, "builds": registry.public()},
        "contract": {
            "contract_id": contract.CONTRACT_ID,
            "contract_version": contract.CONTRACT_VERSION,
            "definition_version": contract.DEFINITION_VERSION,
        },
        "metric_definitions": _metric_definitions(ws, tasks),
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

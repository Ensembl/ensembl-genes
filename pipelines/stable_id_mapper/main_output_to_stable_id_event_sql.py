#!/usr/bin/env python3
"""
Generate quick core-DB SQL from main.py stable-ID mapper output.

The SQL does four things:
  1. Creates backup copies of the core tables it will touch.
  2. Stages stable-ID update decisions for gene/transcript/translation.
  3. Updates core feature tables by internal DB IDs, using a two-phase temporary
     stable_id rename to avoid cascading/collision problems.
  4. Inserts audit/history rows into stable_id_event.

Use --dry-run to emit validation-only SQL that creates temporary staging tables
and reports counts, without persistent backups, updates, or event inserts.

Inputs are file-based. This script does not connect to MySQL.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Optional

FEATURE_ORDER = ("gene", "transcript", "translation")
PK_BY_TYPE = {
    "gene": "gene_id",
    "transcript": "transcript_id",
    "translation": "translation_id",
}
TABLE_BY_TYPE = {
    "gene": "gene",
    "transcript": "transcript",
    "translation": "translation",
}


@dataclass(frozen=True)
class Feature:
    feature_type: str
    stable_id: str
    version: int
    seqid: str
    start: int
    end: int
    strand: str
    parent_stable_id: Optional[str] = None

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class IdRange:
    prefix: str
    start: int
    end: int
    width: int


@dataclass
class IdAllocator:
    range: IdRange
    next_value: int
    reserved: set[str]

    def allocate(self) -> str:
        while self.next_value <= self.range.end:
            candidate = f"{self.range.prefix}{self.next_value:0{self.range.width}d}"
            self.next_value += 1
            if candidate in self.reserved:
                continue
            self.reserved.add(candidate)
            return candidate
        raise RuntimeError(
            f"Ran out of IDs for prefix {self.range.prefix} "
            f"in range {self.range.start}-{self.range.end}"
        )


@dataclass(frozen=True)
class Decision:
    feature_type: str
    action: str
    current_stable_id: Optional[str]
    current_version: int
    old_stable_id: Optional[str]
    old_version: int
    new_stable_id: Optional[str]
    new_version: int
    mapping_session_id: int
    score: float
    reason: str


@dataclass(frozen=True)
class GffRecord:
    stable_id: str
    version: int
    id_namespace: str
    feature_type_lc: str
    seqid: str
    start: int
    end: int
    strand: str
    parent_stable_id: Optional[str]
    translation_stable_id: Optional[str]
    translation_version: int


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if not attr_text or attr_text == ".":
        return attrs
    for item in attr_text.split(";"):
        if not item:
            continue
        if "=" not in item:
            attrs[item] = ""
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def split_stable_id(raw: Optional[str]) -> tuple[Optional[str], Optional[int]]:
    """Return stable ID without GFF3 namespace/version plus optional version."""
    if not raw:
        return None, None
    token = raw.split(",", 1)[0]
    if ":" in token:
        token = token.split(":", 1)[1]
    if "." in token:
        maybe_core, maybe_version = token.rsplit(".", 1)
        if maybe_version.isdigit():
            return maybe_core, int(maybe_version)
    return token, None


def id_namespace(raw: Optional[str]) -> str:
    if raw and ":" in raw:
        return raw.split(":", 1)[0].lower()
    return ""


def parse_version(attrs: Dict[str, str], fallback: Optional[int], default: int) -> int:
    value = attrs.get("version")
    if value is not None:
        try:
            return int(value)
        except ValueError:
            pass
    if fallback is not None:
        return fallback
    return default


def first_parent(attrs: Dict[str, str]) -> Optional[str]:
    parent, _version = split_stable_id(attrs.get("Parent"))
    return parent


def parse_missing_gene_ids(report_path: Path) -> set[str]:
    missing: set[str] = set()
    in_missing_section = False
    with report_path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped == "Missing gene IDs:":
                in_missing_section = True
                continue
            if in_missing_section and stripped.endswith("gene IDs:"):
                break
            if not in_missing_section or not stripped:
                continue
            stable_id, _version = split_stable_id(stripped)
            if stable_id:
                missing.add(stable_id)
    return missing


def default_report_path(mapped_gff: Path) -> Optional[Path]:
    candidate = mapped_gff.parent / "report.txt"
    if candidate.exists():
        return candidate
    return None


def overlap_len(a: Feature, b: Feature) -> int:
    if a.seqid != b.seqid or a.strand != b.strand:
        return 0
    return max(0, min(a.end, b.end) - max(a.start, b.start) + 1)


def overlap_score(a: Feature, b: Feature) -> float:
    overlap = overlap_len(a, b)
    if overlap == 0:
        return 0.0
    return overlap / max(a.length, b.length)


def load_features(
    path: Path,
    default_translation_version: int,
) -> Dict[str, Dict[str, Feature]]:
    features: Dict[str, Dict[str, Feature]] = {
        "gene": {},
        "transcript": {},
        "translation": {},
    }
    records: list[GffRecord] = []

    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            seqid, _source, feature_type, start, end, _score, strand, _phase, attr_text = cols
            feature_type_lc = feature_type.lower()
            attrs = parse_attrs(attr_text)
            start_i = int(start)
            end_i = int(end)
            parent = first_parent(attrs)

            raw_id = attrs.get("ID")
            stable_id, embedded_version = split_stable_id(raw_id)
            if not stable_id:
                continue

            translation_stable_id = None
            translation_version = default_translation_version
            if feature_type_lc == "cds":
                translation_stable_id, translation_embedded_version = split_stable_id(
                    attrs.get("protein_id") or attrs.get("ID")
                )
                translation_version = parse_version(
                    attrs, translation_embedded_version, default_translation_version
                )

            records.append(
                GffRecord(
                    stable_id=stable_id,
                    version=parse_version(attrs, embedded_version, default=1),
                    id_namespace=id_namespace(raw_id),
                    feature_type_lc=feature_type_lc,
                    seqid=seqid,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    parent_stable_id=parent,
                    translation_stable_id=translation_stable_id,
                    translation_version=translation_version,
                )
            )

    gene_ids = {record.stable_id for record in records if record.id_namespace == "gene"}
    if not gene_ids:
        gene_ids = {
            record.stable_id
            for record in records
            if record.parent_stable_id is None and record.feature_type_lc != "cds"
        }
    for record in records:
        if record.stable_id in gene_ids:
            features["gene"][record.stable_id] = Feature(
                feature_type="gene",
                stable_id=record.stable_id,
                version=record.version,
                seqid=record.seqid,
                start=record.start,
                end=record.end,
                strand=record.strand,
            )

    transcript_ids = {
        record.stable_id
        for record in records
        if record.parent_stable_id in gene_ids and record.feature_type_lc != "cds"
    }
    for record in records:
        if record.stable_id in transcript_ids:
            features["transcript"][record.stable_id] = Feature(
                feature_type="transcript",
                stable_id=record.stable_id,
                version=record.version,
                seqid=record.seqid,
                start=record.start,
                end=record.end,
                strand=record.strand,
                parent_stable_id=record.parent_stable_id,
            )

    for record in records:
        if record.feature_type_lc != "cds" or not record.translation_stable_id:
            continue
        if record.parent_stable_id not in transcript_ids:
            continue

        # Several CDS rows can belong to one translation. Keep the full span so
        # overlap-based matching remains sane if needed later.
        previous = features["translation"].get(record.translation_stable_id)
        if previous:
            features["translation"][record.translation_stable_id] = Feature(
                feature_type="translation",
                stable_id=record.translation_stable_id,
                version=record.translation_version,
                seqid=previous.seqid,
                start=min(previous.start, record.start),
                end=max(previous.end, record.end),
                strand=previous.strand,
                parent_stable_id=previous.parent_stable_id,
            )
        else:
            features["translation"][record.translation_stable_id] = Feature(
                feature_type="translation",
                stable_id=record.translation_stable_id,
                version=record.translation_version,
                seqid=record.seqid,
                start=record.start,
                end=record.end,
                strand=record.strand,
                parent_stable_id=record.parent_stable_id,
            )

    return features


def parse_range(value: str) -> IdRange:
    """
    Parse PREFIX:START-END.

    START and END may include leading zeroes; their widest width is preserved
    when rendering new stable IDs.
    """
    try:
        prefix, interval = value.split(":", 1)
        start_text, end_text = interval.split("-", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Bad range {value!r}; expected PREFIX:START-END"
        ) from exc
    if not prefix:
        raise argparse.ArgumentTypeError("Range prefix cannot be empty")
    if not start_text.isdigit() or not end_text.isdigit():
        raise argparse.ArgumentTypeError("Range START and END must be numeric")
    start = int(start_text)
    end = int(end_text)
    if end < start:
        raise argparse.ArgumentTypeError("Range END must be >= START")
    return IdRange(prefix=prefix, start=start, end=end, width=max(len(start_text), len(end_text)))


def reserve_existing_ids(
    *feature_sets: Dict[str, Dict[str, Feature]],
) -> Dict[str, set[str]]:
    reserved = {feature_type: set() for feature_type in FEATURE_ORDER}
    for features_by_type in feature_sets:
        for feature_type in FEATURE_ORDER:
            reserved[feature_type].update(features_by_type[feature_type])
    return reserved


def make_allocators(
    args: argparse.Namespace,
    reserved: Dict[str, set[str]],
) -> Dict[str, IdAllocator]:
    ranges = {
        "gene": args.gene_range,
        "transcript": args.transcript_range,
        "translation": args.translation_range,
    }
    return {
        feature_type: IdAllocator(
            range=id_range,
            next_value=id_range.start,
            reserved=set(reserved[feature_type]),
        )
        for feature_type, id_range in ranges.items()
    }


def best_unused_match(
    query: Feature,
    candidates: Iterable[Feature],
    used_target_ids: set[str],
    min_overlap: float,
    parent_filter: Optional[str] = None,
) -> Optional[Feature]:
    best: Optional[Feature] = None
    best_score = 0.0
    for candidate in candidates:
        if candidate.stable_id in used_target_ids:
            continue
        if parent_filter is not None and candidate.parent_stable_id != parent_filter:
            continue
        score = overlap_score(query, candidate)
        if score >= min_overlap and score > best_score:
            best = candidate
            best_score = score
    return best


def score_for_versions(old_version: int, new_version: int, args: argparse.Namespace) -> float:
    if old_version == new_version:
        return args.score_same_version
    return args.score_changed_version


def build_gene_decisions(
    old_features: Dict[str, Dict[str, Feature]],
    target_features: Dict[str, Dict[str, Feature]],
    mapped_features: Dict[str, Dict[str, Feature]],
    missing_gene_ids: set[str],
    allocators: Dict[str, IdAllocator],
    args: argparse.Namespace,
) -> tuple[list[Decision], dict[str, str], set[str]]:
    decisions: list[Decision] = []
    used_targets: set[str] = set()
    old_gene_to_target_gene: dict[str, str] = {}

    old_genes = old_features["gene"]
    target_genes = target_features["gene"]
    mapped_genes = mapped_features["gene"]

    for old_id in sorted(missing_gene_ids):
        old = old_genes.get(old_id)
        decisions.append(
            Decision(
                feature_type="gene",
                action="missing",
                current_stable_id=None,
                current_version=0,
                old_stable_id=old.stable_id if old else old_id,
                old_version=old.version if old else 1,
                new_stable_id=None,
                new_version=0,
                mapping_session_id=args.mapping_session_id,
                score=args.score_missing,
                reason="old gene listed as missing by main.py report",
            )
        )

    for mapped_id in sorted(mapped_genes):
        mapped = mapped_genes[mapped_id]
        old = old_genes.get(mapped_id)
        target = best_unused_match(
            mapped,
            target_genes.values(),
            used_targets,
            args.min_gene_overlap,
        )
        if target is None:
            decisions.append(
                Decision(
                    feature_type="gene",
                    action="missing",
                    current_stable_id=None,
                    current_version=0,
                    old_stable_id=old.stable_id if old else mapped.stable_id,
                    old_version=old.version if old else mapped.version,
                    new_stable_id=None,
                    new_version=0,
                    mapping_session_id=args.mapping_session_id,
                    score=args.score_missing,
                    reason="mapped gene did not match a target gene by coordinate overlap",
                )
            )
            continue

        used_targets.add(target.stable_id)
        current_stable_id = target.stable_id
        current_version = target.version

        if old is None:
            new_id = allocators["gene"].allocate()
            decisions.append(
                Decision(
                    feature_type="gene",
                    action="new",
                    current_stable_id=current_stable_id,
                    current_version=current_version,
                    old_stable_id=None,
                    old_version=0,
                    new_stable_id=new_id,
                    new_version=1,
                    mapping_session_id=args.mapping_session_id,
                    score=args.score_new,
                    reason="mapped output gene not in old annotation",
                )
            )
            continue

        old_gene_to_target_gene[old.stable_id] = current_stable_id
        decisions.append(
            Decision(
                feature_type="gene",
                action="mapped",
                current_stable_id=current_stable_id,
                current_version=current_version,
                old_stable_id=old.stable_id,
                old_version=old.version,
                new_stable_id=old.stable_id,
                new_version=mapped.version,
                mapping_session_id=args.mapping_session_id,
                score=score_for_versions(old.version, mapped.version, args),
                reason="old gene present in mapped output from main.py",
            )
        )

    return decisions, old_gene_to_target_gene, missing_gene_ids


def build_transcript_decisions(
    old_features: Dict[str, Dict[str, Feature]],
    target_features: Dict[str, Dict[str, Feature]],
    mapped_features: Dict[str, Dict[str, Feature]],
    old_gene_to_target_gene: dict[str, str],
    missing_gene_ids: set[str],
    allocators: Dict[str, IdAllocator],
    args: argparse.Namespace,
) -> tuple[list[Decision], dict[str, str], set[str]]:
    decisions: list[Decision] = []
    used_targets: set[str] = set()
    old_transcript_to_target_transcript: dict[str, str] = {}
    missing_transcript_ids: set[str] = set()

    old_transcripts = old_features["transcript"]
    target_transcripts = target_features["transcript"]
    mapped_transcripts = mapped_features["transcript"]

    for mapped_id in sorted(mapped_transcripts):
        mapped = mapped_transcripts[mapped_id]
        old = old_transcripts.get(mapped_id)
        target_parent = (
            old_gene_to_target_gene.get(old.parent_stable_id)
            if old and old.parent_stable_id
            else None
        )
        target = (
            best_unused_match(
                mapped,
                target_transcripts.values(),
                used_targets,
                args.min_transcript_overlap,
                parent_filter=target_parent,
            )
            if mapped and target_parent
            else None
        )
        if target is None and mapped:
            target = best_unused_match(
                mapped,
                target_transcripts.values(),
                used_targets,
                args.min_transcript_overlap,
            )

        if target is None:
            decisions.append(
                Decision(
                    feature_type="transcript",
                    action="missing",
                    current_stable_id=None,
                    current_version=0,
                    old_stable_id=old.stable_id if old else mapped.stable_id,
                    old_version=old.version if old else mapped.version,
                    new_stable_id=None,
                    new_version=0,
                    mapping_session_id=args.mapping_session_id,
                    score=args.score_missing,
                    reason="mapped transcript did not match a target transcript by coordinate overlap",
                )
            )
            continue

        used_targets.add(target.stable_id)
        current_stable_id = target.stable_id
        current_version = target.version

        if old is None:
            new_id = allocators["transcript"].allocate()
            decisions.append(
                Decision(
                    feature_type="transcript",
                    action="new",
                    current_stable_id=current_stable_id,
                    current_version=current_version,
                    old_stable_id=None,
                    old_version=0,
                    new_stable_id=new_id,
                    new_version=1,
                    mapping_session_id=args.mapping_session_id,
                    score=args.score_new,
                    reason="mapped output transcript not in old annotation",
                )
            )
            continue

        old_transcript_to_target_transcript[old.stable_id] = current_stable_id
        decisions.append(
            Decision(
                feature_type="transcript",
                action="mapped",
                current_stable_id=current_stable_id,
                current_version=current_version,
                old_stable_id=old.stable_id,
                old_version=old.version,
                new_stable_id=old.stable_id,
                new_version=mapped.version,
                mapping_session_id=args.mapping_session_id,
                score=score_for_versions(old.version, mapped.version, args),
                reason="old transcript present in mapped output from main.py",
            )
        )

    for old_id in sorted(old_transcripts):
        old = old_transcripts[old_id]
        if old.stable_id in mapped_transcripts:
            continue
        if old.parent_stable_id not in missing_gene_ids:
            continue
        missing_transcript_ids.add(old.stable_id)
        decisions.append(
            Decision(
                feature_type="transcript",
                action="missing",
                current_stable_id=None,
                current_version=0,
                old_stable_id=old.stable_id,
                old_version=old.version,
                new_stable_id=None,
                new_version=0,
                mapping_session_id=args.mapping_session_id,
                score=args.score_missing,
                reason="old transcript belongs to a gene listed as missing by main.py report",
            )
        )

    return decisions, old_transcript_to_target_transcript, missing_transcript_ids


def translations_by_parent(features: Dict[str, Feature]) -> Dict[str, list[Feature]]:
    by_parent: Dict[str, list[Feature]] = {}
    for feature in features.values():
        if not feature.parent_stable_id:
            continue
        by_parent.setdefault(feature.parent_stable_id, []).append(feature)
    return by_parent


def build_translation_decisions(
    old_features: Dict[str, Dict[str, Feature]],
    target_features: Dict[str, Dict[str, Feature]],
    old_transcript_to_target_transcript: dict[str, str],
    missing_transcript_ids: set[str],
    allocators: Dict[str, IdAllocator],
    args: argparse.Namespace,
) -> list[Decision]:
    decisions: list[Decision] = []
    used_targets: set[str] = set()

    old_by_parent = translations_by_parent(old_features["translation"])
    target_by_parent = translations_by_parent(target_features["translation"])

    mapped_or_missing_transcripts = set(old_transcript_to_target_transcript) | missing_transcript_ids
    for old_transcript_id in sorted(mapped_or_missing_transcripts):
        old_translations = old_by_parent.get(old_transcript_id, [])
        if not old_translations:
            continue
        target_transcript_id = old_transcript_to_target_transcript.get(old_transcript_id)
        target_translations = target_by_parent.get(target_transcript_id or "", [])

        if len(old_translations) == 1 and len(target_translations) == 1:
            old = old_translations[0]
            target = target_translations[0]
            used_targets.add(target.stable_id)
            decisions.append(
                Decision(
                    feature_type="translation",
                    action="mapped",
                    current_stable_id=target.stable_id,
                    current_version=target.version,
                    old_stable_id=old.stable_id,
                    old_version=old.version,
                    new_stable_id=old.stable_id,
                    new_version=old.version,
                    mapping_session_id=args.mapping_session_id,
                    score=args.score_translation_from_transcript,
                    reason="single old and target translation under a mapped transcript",
                )
            )
            continue

        for old in old_translations:
            decisions.append(
                Decision(
                    feature_type="translation",
                    action="missing",
                    current_stable_id=None,
                    current_version=0,
                    old_stable_id=old.stable_id,
                    old_version=old.version,
                    new_stable_id=None,
                    new_version=0,
                    mapping_session_id=args.mapping_session_id,
                    score=args.score_missing,
                    reason="old translation could not be paired under a mapped transcript",
                )
            )

    for target_id in sorted(target_features["translation"]):
        if target_id in used_targets:
            continue
        target = target_features["translation"][target_id]
        new_id = allocators["translation"].allocate()
        decisions.append(
            Decision(
                feature_type="translation",
                action="new",
                current_stable_id=target.stable_id,
                current_version=target.version,
                old_stable_id=None,
                old_version=0,
                new_stable_id=new_id,
                new_version=1,
                mapping_session_id=args.mapping_session_id,
                score=args.score_new,
                reason="target translation did not receive an old stable ID",
            )
        )

    return decisions


def build_decisions(
    old_features: Dict[str, Dict[str, Feature]],
    target_features: Dict[str, Dict[str, Feature]],
    mapped_features: Dict[str, Dict[str, Feature]],
    missing_gene_ids: set[str],
    allocators: Dict[str, IdAllocator],
    args: argparse.Namespace,
) -> list[Decision]:
    gene_decisions, old_gene_to_target_gene, missing_gene_ids = build_gene_decisions(
        old_features, target_features, mapped_features, missing_gene_ids, allocators, args
    )
    (
        transcript_decisions,
        old_transcript_to_target_transcript,
        missing_transcript_ids,
    ) = build_transcript_decisions(
        old_features,
        target_features,
        mapped_features,
        old_gene_to_target_gene,
        missing_gene_ids,
        allocators,
        args,
    )
    translation_decisions = (
        build_translation_decisions(
            old_features,
            target_features,
            old_transcript_to_target_transcript,
            missing_transcript_ids,
            allocators,
            args,
        )
        if args.include_translations
        else []
    )
    return gene_decisions + transcript_decisions + translation_decisions


def sql_string(value: Optional[str]) -> str:
    if value is None:
        return "NULL"
    return "'" + value.replace("\\", "\\\\").replace("'", "''") + "'"


def sql_number(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    if value.is_integer():
        return str(int(value))
    return repr(value)


def write_values_insert(
    handle,
    table_name: str,
    columns: list[str],
    rows: list[list[str]],
    batch_size: int,
) -> None:
    if not rows:
        return
    col_text = ", ".join(columns)
    for start in range(0, len(rows), batch_size):
        batch = rows[start : start + batch_size]
        handle.write(f"INSERT INTO {table_name} ({col_text}) VALUES\n")
        handle.write(",\n".join("  (" + ", ".join(row) + ")" for row in batch))
        handle.write(";\n\n")


def write_sql(decisions: list[Decision], path: Path, args: argparse.Namespace) -> None:
    dry_run = args.dry_run
    update_decisions = [
        decision
        for decision in decisions
        if decision.current_stable_id and decision.new_stable_id
    ]

    backup_tables = ["gene", "transcript", "translation", "stable_id_event"]
    if not args.include_translations:
        backup_tables.remove("translation")

    with path.open("w") as handle:
        handle.write("-- Generated by main_output_to_stable_id_event_sql.py\n")
        handle.write("-- Review before running against a core database.\n")
        if dry_run:
            handle.write(
                "-- DRY RUN: creates temporary staging tables and SELECT checks only.\n"
                "-- No backup tables, core updates, stable_id_event rows, or deletes are emitted.\n\n"
            )
        else:
            handle.write("-- Backup table creation is DDL and autocommits in MySQL.\n\n")

            for table in backup_tables:
                handle.write(
                    f"CREATE TABLE {args.backup_prefix}_{table} AS SELECT * FROM {table};\n"
                )
            handle.write("\n")

        handle.write("START TRANSACTION;\n\n")

        if dry_run:
            decision_table = "tmp_stable_id_mapper_decision"
            handle.write(
                f"CREATE TEMPORARY TABLE {decision_table} (\n"
                "  type VARCHAR(32) NOT NULL,\n"
                "  action VARCHAR(32) NOT NULL,\n"
                "  current_stable_id VARCHAR(128) NULL,\n"
                "  current_version INT NOT NULL,\n"
                "  old_stable_id VARCHAR(128) NULL,\n"
                "  old_version INT NOT NULL,\n"
                "  new_stable_id VARCHAR(128) NULL,\n"
                "  new_version INT NOT NULL,\n"
                "  mapping_session_id INT NOT NULL,\n"
                "  score DOUBLE NOT NULL,\n"
                "  reason VARCHAR(255) NOT NULL\n"
                ");\n\n"
            )
            decision_rows = [
                [
                    sql_string(decision.feature_type),
                    sql_string(decision.action),
                    sql_string(decision.current_stable_id),
                    str(decision.current_version),
                    sql_string(decision.old_stable_id),
                    str(decision.old_version),
                    sql_string(decision.new_stable_id),
                    str(decision.new_version),
                    str(decision.mapping_session_id),
                    sql_number(decision.score),
                    sql_string(decision.reason[:255]),
                ]
                for decision in decisions
            ]
            write_values_insert(
                handle,
                decision_table,
                [
                    "type",
                    "action",
                    "current_stable_id",
                    "current_version",
                    "old_stable_id",
                    "old_version",
                    "new_stable_id",
                    "new_version",
                    "mapping_session_id",
                    "score",
                    "reason",
                ],
                decision_rows,
                args.batch_size,
            )
            handle.write(
                "-- Decision summary expected from this mapper run.\n"
                f"SELECT type, action, COUNT(*) AS decisions\n"
                f"FROM {decision_table}\n"
                "GROUP BY type, action\n"
                "ORDER BY type, action;\n\n"
            )

        if args.replace_events_for_session and not dry_run:
            handle.write(
                "DELETE FROM stable_id_event "
                f"WHERE mapping_session_id = {args.mapping_session_id};\n\n"
            )
        elif args.replace_events_for_session:
            handle.write(
                "-- Dry-run: skipped DELETE FROM stable_id_event for this mapping session.\n\n"
            )

        for feature_type in FEATURE_ORDER:
            type_updates = [
                decision for decision in update_decisions if decision.feature_type == feature_type
            ]
            if not type_updates:
                continue

            update_table = f"tmp_stable_id_mapper_{feature_type}_update"
            pk_table = f"tmp_stable_id_mapper_{feature_type}_pk"
            db_table = TABLE_BY_TYPE[feature_type]
            pk_col = PK_BY_TYPE[feature_type]

            handle.write(
                f"CREATE TEMPORARY TABLE {update_table} (\n"
                "  current_stable_id VARCHAR(128) NOT NULL,\n"
                "  new_stable_id VARCHAR(128) NOT NULL,\n"
                "  new_version INT NOT NULL,\n"
                "  old_stable_id VARCHAR(128) NULL,\n"
                "  old_version INT NOT NULL,\n"
                "  score DOUBLE NOT NULL,\n"
                "  action VARCHAR(32) NOT NULL,\n"
                "  reason VARCHAR(255) NOT NULL,\n"
                "  PRIMARY KEY (current_stable_id)\n"
                ");\n\n"
            )

            rows = [
                [
                    sql_string(decision.current_stable_id),
                    sql_string(decision.new_stable_id),
                    str(decision.new_version),
                    sql_string(decision.old_stable_id),
                    str(decision.old_version),
                    sql_number(decision.score),
                    sql_string(decision.action),
                    sql_string(decision.reason[:255]),
                ]
                for decision in type_updates
            ]
            write_values_insert(
                handle,
                update_table,
                [
                    "current_stable_id",
                    "new_stable_id",
                    "new_version",
                    "old_stable_id",
                    "old_version",
                    "score",
                    "action",
                    "reason",
                ],
                rows,
                args.batch_size,
            )

            handle.write(
                f"CREATE TEMPORARY TABLE {pk_table} AS\n"
                f"SELECT f.{pk_col}, u.*\n"
                f"FROM {db_table} f\n"
                f"JOIN {update_table} u ON f.stable_id = u.current_stable_id;\n\n"
            )

            handle.write(
                f"-- Fail early if any staged {feature_type} stable IDs are not present in the DB.\n"
                f"SELECT '{feature_type}' AS type, COUNT(*) AS staged_rows FROM {update_table};\n"
                f"SELECT '{feature_type}' AS type, COUNT(*) AS matched_db_rows FROM {pk_table};\n\n"
            )
            handle.write(
                f"-- Inspect staged {feature_type} updates that did not match a DB row.\n"
                f"SELECT '{feature_type}' AS type, u.current_stable_id, u.new_stable_id, "
                "u.action, u.reason\n"
                f"FROM {update_table} u\n"
                f"LEFT JOIN {pk_table} p ON p.current_stable_id = u.current_stable_id\n"
                "WHERE p.current_stable_id IS NULL\n"
                "LIMIT 20;\n\n"
            )
            handle.write(
                f"-- Inspect duplicate DB matches for staged {feature_type} stable IDs.\n"
                f"SELECT '{feature_type}' AS type, current_stable_id, COUNT(*) AS db_matches\n"
                f"FROM {pk_table}\n"
                "GROUP BY current_stable_id\n"
                "HAVING COUNT(*) > 1\n"
                "LIMIT 20;\n\n"
            )

            if dry_run:
                handle.write(
                    f"-- Dry-run: skipped two-phase UPDATE of {db_table}.\n\n"
                )
                continue

            handle.write(
                f"-- Two-phase update avoids A->B, B->C cascading/collision issues.\n"
                f"UPDATE {db_table} f\n"
                f"JOIN {pk_table} p ON f.{pk_col} = p.{pk_col}\n"
                f"SET f.stable_id = CONCAT('__stable_id_mapper_tmp_{feature_type}_', f.{pk_col}),\n"
                "    f.version = 0;\n\n"
            )
            handle.write(
                f"UPDATE {db_table} f\n"
                f"JOIN {pk_table} p ON f.{pk_col} = p.{pk_col}\n"
                "SET f.stable_id = p.new_stable_id,\n"
                "    f.version = p.new_version;\n\n"
            )

        if dry_run:
            handle.write(
                "-- Dry-run: skipped INSERT INTO stable_id_event. "
                "Review the temporary decision-table summary and TSV output for all events.\n"
                "ROLLBACK;\n"
            )
            return

        event_rows = [
            [
                sql_string(decision.old_stable_id),
                str(decision.old_version),
                sql_string(decision.new_stable_id),
                str(decision.new_version),
                str(decision.mapping_session_id),
                sql_string(decision.feature_type),
                sql_number(decision.score),
            ]
            for decision in decisions
        ]
        write_values_insert(
            handle,
            "stable_id_event",
            [
                "old_stable_id",
                "old_version",
                "new_stable_id",
                "new_version",
                "mapping_session_id",
                "type",
                "score",
            ],
            event_rows,
            args.batch_size,
        )

        handle.write("COMMIT;\n")


def write_tsv(decisions: list[Decision], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "type",
                "action",
                "current_stable_id",
                "current_version",
                "old_stable_id",
                "old_version",
                "new_stable_id",
                "new_version",
                "mapping_session_id",
                "score",
                "reason",
            ]
        )
        for decision in decisions:
            writer.writerow(
                [
                    decision.feature_type,
                    decision.action,
                    decision.current_stable_id or "",
                    decision.current_version,
                    decision.old_stable_id or "",
                    decision.old_version,
                    decision.new_stable_id or "",
                    decision.new_version,
                    decision.mapping_session_id,
                    decision.score,
                    decision.reason,
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate core-table update SQL plus stable_id_event rows from "
            "main.py mapper output."
        )
    )
    parser.add_argument("--ref-gff", required=True, type=Path)
    parser.add_argument("--target-gff", required=True, type=Path)
    parser.add_argument("--mapped-gff", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--mapping-session-id", required=True, type=int)
    parser.add_argument("--gene-range", required=True, type=parse_range)
    parser.add_argument("--transcript-range", required=True, type=parse_range)
    parser.add_argument("--translation-range", required=True, type=parse_range)
    parser.add_argument("--output-sql", required=True, type=Path)
    parser.add_argument("--output-tsv", type=Path)
    parser.add_argument(
        "--backup-prefix",
        default=f"stable_id_mapper_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
    )
    parser.add_argument("--include-translations", action="store_true")
    parser.add_argument("--default-translation-version", type=int, default=1)
    parser.add_argument("--score-same-version", type=float, default=1.0)
    parser.add_argument("--score-changed-version", type=float, default=0.0)
    parser.add_argument("--score-missing", type=float, default=0.0)
    parser.add_argument("--score-new", type=float, default=0.0)
    parser.add_argument("--score-translation-from-transcript", type=float, default=0.0)
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--min-gene-overlap", type=float, default=0.10)
    parser.add_argument("--min-transcript-overlap", type=float, default=0.10)
    parser.add_argument("--replace-events-for-session", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Emit validation SQL only: temporary staging tables and SELECT "
            "checks, with no backups, core-table updates, event deletes, or "
            "stable_id_event inserts."
        ),
    )
    return parser.parse_args()


def print_summary(decisions: list[Decision]) -> None:
    counts: Dict[tuple[str, str], int] = {}
    for decision in decisions:
        key = (decision.feature_type, decision.action)
        counts[key] = counts.get(key, 0) + 1
    mapped_or_new = sum(1 for decision in decisions if decision.new_stable_id)
    missing = sum(1 for decision in decisions if decision.new_stable_id is None)
    parts = [
        f"{feature_type}:{action}={count}"
        for (feature_type, action), count in sorted(counts.items())
    ]
    sys.stderr.write(
        f"Generated {len(decisions)} decisions "
        f"({mapped_or_new} assigned, {missing} missing): "
        + ", ".join(parts)
        + "\n"
    )


def main() -> None:
    args = parse_args()
    old_features = load_features(args.ref_gff, args.default_translation_version)
    target_features = load_features(args.target_gff, args.default_translation_version)
    mapped_features = load_features(args.mapped_gff, args.default_translation_version)
    report_path = args.report or default_report_path(args.mapped_gff)
    if report_path:
        missing_gene_ids = parse_missing_gene_ids(report_path)
        sys.stderr.write(
            f"Using missing gene list from {report_path}: "
            f"{len(missing_gene_ids)} genes\n"
        )
    else:
        missing_gene_ids = set(old_features["gene"]) - set(mapped_features["gene"])
        sys.stderr.write(
            "Warning: no main.py report found; missing genes inferred from "
            "ref-gff minus mapped-gff.\n"
        )

    reserved = reserve_existing_ids(old_features, target_features, mapped_features)
    allocators = make_allocators(args, reserved)
    decisions = build_decisions(
        old_features, target_features, mapped_features, missing_gene_ids, allocators, args
    )

    write_sql(decisions, args.output_sql, args)
    if args.output_tsv:
        write_tsv(decisions, args.output_tsv)
    print_summary(decisions)

    if args.include_translations:
        sys.stderr.write(
            "Note: translation mappings are inferred only when one old and one "
            "target translation sit under a mapped transcript.\n"
        )


if __name__ == "__main__":
    main()

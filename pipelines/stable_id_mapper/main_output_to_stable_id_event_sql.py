#!/usr/bin/env python3
"""
Generate SQL stable-ID updates from main.py mapper output.

This is a file-based post-processor. It does not connect to MySQL.
Review the generated SQL before running it against a core database.
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
import textwrap
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional

FEATURE_ORDER = ("gene", "transcript", "translation")
ACTION_ORDER = ("mapped", "missing", "new")
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


@dataclass
class IdAllocator:
    prefix: str
    next_value: int
    end: int
    width: int
    reserved: set[str]

    def allocate(self) -> str:
        while self.next_value <= self.end:
            candidate = f"{self.prefix}{self.next_value:0{self.width}d}"
            self.next_value += 1
            if candidate in self.reserved:
                continue
            self.reserved.add(candidate)
            return candidate
        raise RuntimeError(
            f"Ran out of IDs for {self.prefix}:{self.next_value}-{self.end}"
        )


@dataclass(frozen=True)
class RawFeature:
    stable_id: str
    version: int
    feature_type: str
    seqid: str
    start: int
    end: int
    strand: str
    parent_ids: tuple[str, ...]


@dataclass(frozen=True)
class CdsRow:
    parent_stable_id: str
    parent_version: Optional[int]
    seqid: str
    start: int
    end: int
    strand: str


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if not attr_text or attr_text == ".":
        return attrs
    for item in attr_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
        else:
            attrs[item] = ""
    return attrs


def split_stable_id(raw: Optional[str]) -> tuple[Optional[str], Optional[int]]:
    if not raw:
        return None, None
    token = raw.strip()
    if not token:
        return None, None
    if ":" in token:
        token = token.split(":", 1)[1]
    if "." in token:
        stable_id, maybe_version = token.rsplit(".", 1)
        if maybe_version.isdigit():
            return stable_id, int(maybe_version)
    return token, None


def parent_ids(raw: Optional[str]) -> tuple[tuple[str, Optional[int]], ...]:
    if not raw:
        return ()
    out: list[tuple[str, Optional[int]]] = []
    for token in raw.split(","):
        stable_id, version = split_stable_id(token)
        if stable_id:
            out.append((stable_id, version))
    return tuple(out)


def parse_version(
    attrs: dict[str, str],
    embedded_version: Optional[int],
    default: int = 1,
) -> int:
    version_text = attrs.get("version")
    if version_text is not None:
        try:
            return int(version_text)
        except ValueError:
            pass
    if embedded_version is not None:
        return embedded_version
    return default


def parse_gff3(path: str | Path) -> dict[str, dict[str, Feature]]:
    features: dict[str, dict[str, Feature]] = {
        "gene": {},
        "transcript": {},
        "translation": {},
    }
    raw_features: list[RawFeature] = []
    cds_rows: list[CdsRow] = []

    with Path(path).open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            seqid, _source, feature_type, start, end, _score, strand, _phase, attrs_text = fields
            attrs = parse_attrs(attrs_text)
            start_i = int(start)
            end_i = int(end)
            feature_type_lc = feature_type.lower()

            if feature_type_lc == "exon":
                continue

            parents = parent_ids(attrs.get("Parent"))
            if feature_type_lc == "cds":
                for parent_stable_id, parent_version in parents:
                    cds_rows.append(
                        CdsRow(
                            parent_stable_id=parent_stable_id,
                            parent_version=parent_version,
                            seqid=seqid,
                            start=start_i,
                            end=end_i,
                            strand=strand,
                        )
                    )
                continue

            stable_id, embedded_version = split_stable_id(attrs.get("ID"))
            if not stable_id:
                continue
            raw_features.append(
                RawFeature(
                    stable_id=stable_id,
                    version=parse_version(attrs, embedded_version),
                    feature_type=feature_type,
                    seqid=seqid,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    parent_ids=tuple(parent for parent, _version in parents),
                )
            )

    gene_ids = {record.stable_id for record in raw_features if not record.parent_ids}
    for record in raw_features:
        if record.stable_id not in gene_ids:
            continue
        features["gene"][record.stable_id] = Feature(
            stable_id=record.stable_id,
            version=record.version,
            seqid=record.seqid,
            start=record.start,
            end=record.end,
            strand=record.strand,
        )

    transcript_ids: set[str] = set()
    for record in raw_features:
        parent_gene = next(
            (parent for parent in record.parent_ids if parent in gene_ids),
            None,
        )
        if parent_gene is None:
            continue
        transcript_ids.add(record.stable_id)
        features["transcript"][record.stable_id] = Feature(
            stable_id=record.stable_id,
            version=record.version,
            seqid=record.seqid,
            start=record.start,
            end=record.end,
            strand=record.strand,
            parent_stable_id=parent_gene,
        )

    for cds in cds_rows:
        if cds.parent_stable_id not in transcript_ids:
            continue
        transcript = features["transcript"][cds.parent_stable_id]
        version = cds.parent_version if cds.parent_version is not None else transcript.version
        previous = features["translation"].get(cds.parent_stable_id)
        if previous is None:
            features["translation"][cds.parent_stable_id] = Feature(
                stable_id=cds.parent_stable_id,
                version=version,
                seqid=cds.seqid,
                start=cds.start,
                end=cds.end,
                strand=cds.strand,
                parent_stable_id=cds.parent_stable_id,
            )
        else:
            features["translation"][cds.parent_stable_id] = Feature(
                stable_id=previous.stable_id,
                version=previous.version,
                seqid=previous.seqid,
                start=min(previous.start, cds.start),
                end=max(previous.end, cds.end),
                strand=previous.strand,
                parent_stable_id=previous.parent_stable_id,
            )

    return features


def parse_missing_gene_ids(report_path: str | Path) -> set[str]:
    missing: set[str] = set()
    in_missing_section = False
    with Path(report_path).open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped == "Missing gene IDs:":
                in_missing_section = True
                continue
            if not in_missing_section:
                continue
            if not stripped:
                continue
            if stripped.endswith(":") and not stripped.startswith("gene:"):
                break
            stable_id, _version = split_stable_id(stripped)
            if stable_id:
                missing.add(stable_id)
    return missing


def overlap_score(a: Feature, b: Feature) -> float:
    if a.seqid != b.seqid or a.strand != b.strand:
        return 0.0
    overlap = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    if overlap == 0:
        return 0.0
    return overlap / max(a.length, b.length)


def best_unused_match(
    query: Feature,
    candidates: Iterable[Feature],
    used_ids: set[str],
    min_overlap: float,
    parent_filter: Optional[str] = None,
) -> Optional[Feature]:
    best: Optional[Feature] = None
    best_score = -1.0
    for candidate in sorted(candidates, key=lambda feature: feature.stable_id):
        if candidate.stable_id in used_ids:
            continue
        if parent_filter is not None and candidate.parent_stable_id != parent_filter:
            continue
        score = overlap_score(query, candidate)
        if score >= min_overlap and score > best_score:
            best = candidate
            best_score = score
    return best


def parse_id_range(value: str) -> tuple[str, int, int, int]:
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
        raise argparse.ArgumentTypeError("Range START and END must be integers")
    start = int(start_text)
    end = int(end_text)
    if end < start:
        raise argparse.ArgumentTypeError("Range END must be >= START")
    return prefix, start, end, max(len(start_text), len(end_text))


def collect_reserved_ids(*feature_sets: dict[str, dict[str, Feature]]) -> set[str]:
    reserved: set[str] = set()
    for features_by_type in feature_sets:
        for features in features_by_type.values():
            reserved.update(features)
    return reserved


def make_allocator(
    id_range: tuple[str, int, int, int],
    reserved_ids: set[str],
) -> IdAllocator:
    prefix, start, end, width = id_range
    return IdAllocator(
        prefix=prefix,
        next_value=start,
        end=end,
        width=width,
        reserved=set(reserved_ids),
    )


def old_version_for(
    features: dict[str, Feature],
    stable_id: str,
    fallback: Feature,
) -> int:
    old = features.get(stable_id)
    return old.version if old else fallback.version


def build_gene_decisions(
    ref_features: dict[str, dict[str, Feature]],
    target_features: dict[str, dict[str, Feature]],
    mapped_features: dict[str, dict[str, Feature]],
    missing_gene_ids: set[str],
    allocators: dict[str, IdAllocator],
    mapping_session_id: int,
    min_overlap: float,
) -> tuple[list[Decision], dict[str, str]]:
    decisions: list[Decision] = []
    used_targets: set[str] = set()
    old_gene_to_target_gene: dict[str, str] = {}

    for old_id in sorted(missing_gene_ids):
        old = ref_features["gene"].get(old_id)
        decisions.append(
            Decision(
                feature_type="gene",
                action="missing",
                current_stable_id=None,
                current_version=0,
                old_stable_id=old_id,
                old_version=old.version if old else 1,
                new_stable_id=None,
                new_version=0,
                mapping_session_id=mapping_session_id,
                score=0.0,
                reason="gene listed as missing by main.py report",
            )
        )

    for old_id in sorted(mapped_features["gene"]):
        mapped = mapped_features["gene"][old_id]
        target = best_unused_match(
            mapped,
            target_features["gene"].values(),
            used_targets,
            min_overlap,
        )
        old_version = old_version_for(ref_features["gene"], old_id, mapped)
        if target is None:
            decisions.append(
                Decision(
                    feature_type="gene",
                    action="missing",
                    current_stable_id=None,
                    current_version=0,
                    old_stable_id=old_id,
                    old_version=old_version,
                    new_stable_id=None,
                    new_version=0,
                    mapping_session_id=mapping_session_id,
                    score=0.0,
                    reason="mapped gene did not match a target gene by coordinate overlap",
                )
            )
            continue

        used_targets.add(target.stable_id)
        old_gene_to_target_gene[old_id] = target.stable_id
        decisions.append(
            Decision(
                feature_type="gene",
                action="mapped",
                current_stable_id=target.stable_id,
                current_version=target.version,
                old_stable_id=old_id,
                old_version=old_version,
                new_stable_id=old_id,
                new_version=mapped.version,
                mapping_session_id=mapping_session_id,
                score=overlap_score(mapped, target),
                reason="mapped gene matched a target gene by coordinate overlap",
            )
        )

    for target_id in sorted(target_features["gene"]):
        if target_id in used_targets:
            continue
        target = target_features["gene"][target_id]
        decisions.append(
            Decision(
                feature_type="gene",
                action="new",
                current_stable_id=target.stable_id,
                current_version=target.version,
                old_stable_id=None,
                old_version=0,
                new_stable_id=allocators["gene"].allocate(),
                new_version=1,
                mapping_session_id=mapping_session_id,
                score=0.0,
                reason="target gene was not claimed by any mapped gene",
            )
        )

    return decisions, old_gene_to_target_gene


def build_transcript_decisions(
    ref_features: dict[str, dict[str, Feature]],
    target_features: dict[str, dict[str, Feature]],
    mapped_features: dict[str, dict[str, Feature]],
    missing_gene_ids: set[str],
    old_gene_to_target_gene: dict[str, str],
    allocators: dict[str, IdAllocator],
    mapping_session_id: int,
    min_overlap: float,
) -> tuple[list[Decision], dict[str, str], set[str]]:
    decisions: list[Decision] = []
    used_targets: set[str] = set()
    old_transcript_to_target_transcript: dict[str, str] = {}
    missing_transcript_ids: set[str] = set()
    mapped_transcript_ids = set(mapped_features["transcript"])

    for old_id in sorted(ref_features["transcript"]):
        old = ref_features["transcript"][old_id]
        if old.parent_stable_id not in missing_gene_ids:
            continue
        if old_id in mapped_transcript_ids:
            continue
        missing_transcript_ids.add(old_id)
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
                mapping_session_id=mapping_session_id,
                score=0.0,
                reason="transcript belongs to a gene listed as missing",
            )
        )

    for old_id in sorted(mapped_features["transcript"]):
        mapped = mapped_features["transcript"][old_id]
        old = ref_features["transcript"].get(old_id)
        old_parent = old.parent_stable_id if old else mapped.parent_stable_id
        target_parent = old_gene_to_target_gene.get(old_parent or "")

        if target_parent is None:
            target = best_unused_match(
                mapped,
                target_features["transcript"].values(),
                used_targets,
                min_overlap,
            )
        else:
            target = best_unused_match(
                mapped,
                target_features["transcript"].values(),
                used_targets,
                min_overlap,
                parent_filter=target_parent,
            )

        old_version = old.version if old else mapped.version
        if target is None:
            missing_transcript_ids.add(old_id)
            decisions.append(
                Decision(
                    feature_type="transcript",
                    action="missing",
                    current_stable_id=None,
                    current_version=0,
                    old_stable_id=old_id,
                    old_version=old_version,
                    new_stable_id=None,
                    new_version=0,
                    mapping_session_id=mapping_session_id,
                    score=0.0,
                    reason="mapped transcript did not match a target transcript by coordinate overlap",
                )
            )
            continue

        used_targets.add(target.stable_id)
        old_transcript_to_target_transcript[old_id] = target.stable_id
        decisions.append(
            Decision(
                feature_type="transcript",
                action="mapped",
                current_stable_id=target.stable_id,
                current_version=target.version,
                old_stable_id=old_id,
                old_version=old_version,
                new_stable_id=old_id,
                new_version=mapped.version,
                mapping_session_id=mapping_session_id,
                score=overlap_score(mapped, target),
                reason="mapped transcript matched a target transcript by coordinate overlap",
            )
        )

    for target_id in sorted(target_features["transcript"]):
        if target_id in used_targets:
            continue
        target = target_features["transcript"][target_id]
        decisions.append(
            Decision(
                feature_type="transcript",
                action="new",
                current_stable_id=target.stable_id,
                current_version=target.version,
                old_stable_id=None,
                old_version=0,
                new_stable_id=allocators["transcript"].allocate(),
                new_version=1,
                mapping_session_id=mapping_session_id,
                score=0.0,
                reason="target transcript was not claimed by any mapped transcript",
            )
        )

    return decisions, old_transcript_to_target_transcript, missing_transcript_ids


def translations_by_parent(
    translations: dict[str, Feature],
) -> dict[str, list[Feature]]:
    grouped: dict[str, list[Feature]] = defaultdict(list)
    for translation in translations.values():
        parent = translation.parent_stable_id or translation.stable_id
        grouped[parent].append(translation)
    for grouped_translations in grouped.values():
        grouped_translations.sort(key=lambda feature: feature.stable_id)
    return grouped


def projected_translation(
    old_translation: Feature,
    mapped_translations: dict[str, Feature],
) -> Feature:
    mapped = mapped_translations.get(old_translation.stable_id)
    if mapped is None:
        return old_translation
    return Feature(
        stable_id=old_translation.stable_id,
        version=old_translation.version,
        seqid=mapped.seqid,
        start=mapped.start,
        end=mapped.end,
        strand=mapped.strand,
        parent_stable_id=old_translation.parent_stable_id,
    )


def new_translation_decision(
    target: Feature,
    allocator: IdAllocator,
    mapping_session_id: int,
    reason: str,
) -> Decision:
    return Decision(
        feature_type="translation",
        action="new",
        current_stable_id=target.stable_id,
        current_version=target.version,
        old_stable_id=None,
        old_version=0,
        new_stable_id=allocator.allocate(),
        new_version=1,
        mapping_session_id=mapping_session_id,
        score=0.0,
        reason=reason,
    )


def missing_translation_decision(
    old: Feature,
    mapping_session_id: int,
    reason: str,
) -> Decision:
    return Decision(
        feature_type="translation",
        action="missing",
        current_stable_id=None,
        current_version=0,
        old_stable_id=old.stable_id,
        old_version=old.version,
        new_stable_id=None,
        new_version=0,
        mapping_session_id=mapping_session_id,
        score=0.0,
        reason=reason,
    )


def mapped_translation_decision(
    old: Feature,
    target: Feature,
    mapping_session_id: int,
    score: float,
    reason: str,
) -> Decision:
    return Decision(
        feature_type="translation",
        action="mapped",
        current_stable_id=target.stable_id,
        current_version=target.version,
        old_stable_id=old.stable_id,
        old_version=old.version,
        new_stable_id=old.stable_id,
        new_version=old.version,
        mapping_session_id=mapping_session_id,
        score=score,
        reason=reason,
    )


def build_translation_decisions(
    ref_features: dict[str, dict[str, Feature]],
    target_features: dict[str, dict[str, Feature]],
    mapped_features: dict[str, dict[str, Feature]],
    old_transcript_to_target_transcript: dict[str, str],
    missing_transcript_ids: set[str],
    allocators: dict[str, IdAllocator],
    mapping_session_id: int,
    min_overlap: float,
) -> list[Decision]:
    decisions: list[Decision] = []
    used_targets: set[str] = set()
    old_by_parent = translations_by_parent(ref_features["translation"])
    target_by_parent = translations_by_parent(target_features["translation"])
    mapped_translations = mapped_features["translation"]

    for old_transcript_id in sorted(old_transcript_to_target_transcript):
        target_transcript_id = old_transcript_to_target_transcript[old_transcript_id]
        old_translations = old_by_parent.get(old_transcript_id, [])
        target_translations = target_by_parent.get(target_transcript_id, [])
        old_count = len(old_translations)
        target_count = len(target_translations)

        if old_count == 0 and target_count == 0:
            continue

        if old_count == 1 and target_count == 1:
            old = old_translations[0]
            target = target_translations[0]
            used_targets.add(target.stable_id)
            decisions.append(
                mapped_translation_decision(
                    old,
                    target,
                    mapping_session_id,
                    1.0,
                    "single old and target translation under a mapped transcript",
                )
            )
            continue

        if old_count == 1 and target_count == 0:
            decisions.append(
                missing_translation_decision(
                    old_translations[0],
                    mapping_session_id,
                    "old translation has no target translation under mapped transcript",
                )
            )
            continue

        if old_count == 0 and target_count == 1:
            target = target_translations[0]
            used_targets.add(target.stable_id)
            decisions.append(
                new_translation_decision(
                    target,
                    allocators["translation"],
                    mapping_session_id,
                    "target translation has no old translation under mapped transcript",
                )
            )
            continue

        sys.stderr.write(
            "Warning: translation count mismatch for transcript "
            f"{old_transcript_id}: old={old_count}, target={target_count}\n"
        )
        used_in_pair: set[str] = set()
        projected_old = [
            projected_translation(old, mapped_translations)
            for old in old_translations
        ]
        old_by_stable_id = {old.stable_id: old for old in old_translations}
        for old_query in projected_old:
            target = best_unused_match(
                old_query,
                target_translations,
                used_in_pair,
                min_overlap,
            )
            old = old_by_stable_id[old_query.stable_id]
            if target is None:
                decisions.append(
                    missing_translation_decision(
                        old,
                        mapping_session_id,
                        "old translation did not match a target translation by coordinate overlap",
                    )
                )
                continue
            used_in_pair.add(target.stable_id)
            used_targets.add(target.stable_id)
            decisions.append(
                mapped_translation_decision(
                    old,
                    target,
                    mapping_session_id,
                    overlap_score(old_query, target),
                    "old translation matched a target translation by coordinate overlap",
                )
            )
        for target in target_translations:
            if target.stable_id in used_in_pair:
                continue
            used_targets.add(target.stable_id)
            decisions.append(
                new_translation_decision(
                    target,
                    allocators["translation"],
                    mapping_session_id,
                    "target translation was not matched by an old translation",
                )
            )

    for old_transcript_id in sorted(missing_transcript_ids):
        for old in old_by_parent.get(old_transcript_id, []):
            decisions.append(
                missing_translation_decision(
                    old,
                    mapping_session_id,
                    "old translation belongs to a transcript that was not mapped",
                )
            )

    for target_id in sorted(target_features["translation"]):
        target = target_features["translation"][target_id]
        if target.stable_id in used_targets:
            continue
        decisions.append(
            new_translation_decision(
                target,
                allocators["translation"],
                mapping_session_id,
                "target translation belongs to a transcript that was not mapped",
            )
        )

    return decisions


def build_decisions(
    ref_features: dict[str, dict[str, Feature]],
    target_features: dict[str, dict[str, Feature]],
    mapped_features: dict[str, dict[str, Feature]],
    missing_gene_ids: set[str],
    allocators: dict[str, IdAllocator],
    mapping_session_id: int,
    min_overlap: float,
    include_translations: bool,
) -> list[Decision]:
    gene_decisions, old_gene_to_target_gene = build_gene_decisions(
        ref_features,
        target_features,
        mapped_features,
        missing_gene_ids,
        allocators,
        mapping_session_id,
        min_overlap,
    )
    (
        transcript_decisions,
        old_transcript_to_target_transcript,
        missing_transcript_ids,
    ) = build_transcript_decisions(
        ref_features,
        target_features,
        mapped_features,
        missing_gene_ids,
        old_gene_to_target_gene,
        allocators,
        mapping_session_id,
        min_overlap,
    )
    translation_decisions = (
        build_translation_decisions(
            ref_features,
            target_features,
            mapped_features,
            old_transcript_to_target_transcript,
            missing_transcript_ids,
            allocators,
            mapping_session_id,
            min_overlap,
        )
        if include_translations
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
    column_text = ", ".join(columns)
    for start in range(0, len(rows), batch_size):
        batch = rows[start : start + batch_size]
        handle.write(f"INSERT INTO {table_name} ({column_text}) VALUES\n")
        handle.write(",\n".join("  (" + ", ".join(row) + ")" for row in batch))
        handle.write(";\n\n")


def decision_rows(decisions: list[Decision]) -> list[list[str]]:
    return [
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


def write_sql(decisions: list[Decision], path: str | Path, args: argparse.Namespace) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    feature_order = ["gene", "transcript"]
    if args.include_translations:
        feature_order.append("translation")

    with path.open("w") as handle:
        handle.write("-- Generated by main_output_to_stable_id_event_sql.py\n")
        handle.write("-- Review before running against a core database.\n")
        handle.write("-- This SQL was generated from main.py mapper output.\n\n")

        if not args.dry_run:
            for table in ["gene", "transcript", "stable_id_event"]:
                handle.write(
                    f"CREATE TABLE {args.backup_prefix}_{table} AS SELECT * FROM {table};\n"
                )
            if args.include_translations:
                handle.write(
                    f"CREATE TABLE {args.backup_prefix}_translation AS SELECT * FROM translation;\n"
                )
            handle.write("\n")

        handle.write("START TRANSACTION;\n\n")

        if args.dry_run:
            decision_table = "tmp_stable_id_mapper_decisions"
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
                decision_rows(decisions),
                args.batch_size,
            )
            handle.write(
                f"SELECT type, action, COUNT(*) AS decisions\n"
                f"FROM {decision_table}\n"
                "GROUP BY type, action\n"
                "ORDER BY type, action;\n\n"
            )

        if args.replace_events_for_session and not args.dry_run:
            handle.write(
                "DELETE FROM stable_id_event "
                f"WHERE mapping_session_id = {args.mapping_session_id};\n\n"
            )

        for feature_type in feature_order:
            updates = [
                decision
                for decision in decisions
                if decision.feature_type == feature_type
                and decision.current_stable_id
                and decision.new_stable_id
            ]
            if not updates:
                continue

            table_name = TABLE_BY_TYPE[feature_type]
            pk_col = PK_BY_TYPE[feature_type]
            update_table = f"tmp_stable_id_mapper_{feature_type}_update"
            pk_table = f"tmp_stable_id_mapper_{feature_type}_pk"

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
            update_rows = [
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
                for decision in updates
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
                update_rows,
                args.batch_size,
            )
            handle.write(
                f"CREATE TEMPORARY TABLE {pk_table} AS\n"
                f"SELECT f.{pk_col}, u.*\n"
                f"FROM {table_name} f\n"
                f"JOIN {update_table} u ON f.stable_id = u.current_stable_id;\n\n"
            )
            handle.write(
                f"SELECT '{feature_type}' AS type,\n"
                f"       (SELECT COUNT(*) FROM {update_table}) AS staged_count,\n"
                f"       (SELECT COUNT(*) FROM {pk_table}) AS matched_count;\n\n"
            )
            handle.write(
                f"SELECT '{feature_type}' AS type, u.current_stable_id, u.new_stable_id, "
                "u.action, u.reason\n"
                f"FROM {update_table} u\n"
                f"LEFT JOIN {pk_table} p ON p.current_stable_id = u.current_stable_id\n"
                "WHERE p.current_stable_id IS NULL\n"
                "LIMIT 20;\n\n"
            )
            handle.write(
                f"SELECT '{feature_type}' AS type, current_stable_id, COUNT(*) AS db_matches\n"
                f"FROM {pk_table}\n"
                "GROUP BY current_stable_id\n"
                "HAVING COUNT(*) > 1\n"
                "LIMIT 20;\n\n"
            )

            if args.dry_run:
                continue

            handle.write(
                f"UPDATE {table_name} f\n"
                f"JOIN {pk_table} p ON f.{pk_col} = p.{pk_col}\n"
                f"SET f.stable_id = CONCAT('__tmp_{feature_type}_', f.{pk_col}),\n"
                "    f.version = 0;\n\n"
            )
            handle.write(
                f"UPDATE {table_name} f\n"
                f"JOIN {pk_table} p ON f.{pk_col} = p.{pk_col}\n"
                "SET f.stable_id = p.new_stable_id,\n"
                "    f.version = p.new_version;\n\n"
            )

        if args.dry_run:
            handle.write("ROLLBACK;\n")
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


def write_tsv(decisions: list[Decision], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
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


def print_summary(decisions: list[Decision]) -> None:
    counts = Counter((decision.feature_type, decision.action) for decision in decisions)
    assigned = sum(1 for decision in decisions if decision.new_stable_id is not None)
    missing = sum(1 for decision in decisions if decision.action == "missing")
    parts: list[str] = []
    for feature_type in FEATURE_ORDER:
        for action in ACTION_ORDER:
            count = counts.get((feature_type, action))
            if count:
                parts.append(f"{feature_type}:{action}={count}")
    sys.stderr.write(
        f"Generated {len(decisions)} decisions "
        f"({assigned} assigned, {missing} missing): "
        + ", ".join(parts)
        + "\n"
    )


def make_allocators(
    args: argparse.Namespace,
    reserved_ids: set[str],
) -> dict[str, IdAllocator]:
    return {
        "gene": make_allocator(args.gene_range, reserved_ids),
        "transcript": make_allocator(args.transcript_range, reserved_ids),
        "translation": make_allocator(args.translation_range, reserved_ids),
    }


def run_pipeline(args: argparse.Namespace) -> list[Decision]:
    ref_features = parse_gff3(args.ref_gff)
    target_features = parse_gff3(args.target_gff)
    mapped_features = parse_gff3(args.mapped_gff)
    missing_gene_ids = parse_missing_gene_ids(args.report)
    reserved_ids = collect_reserved_ids(ref_features, target_features, mapped_features)
    allocators = make_allocators(args, reserved_ids)
    decisions = build_decisions(
        ref_features,
        target_features,
        mapped_features,
        missing_gene_ids,
        allocators,
        args.mapping_session_id,
        args.min_overlap,
        args.include_translations,
    )
    write_sql(decisions, args.output_sql, args)
    if args.output_tsv:
        write_tsv(decisions, args.output_tsv)
    print_summary(decisions)
    return decisions


def write_test_file(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip(), encoding="utf-8")


def assert_has_decision(
    decisions: list[Decision],
    feature_type: str,
    action: str,
    current_stable_id: Optional[str] = None,
    old_stable_id: Optional[str] = None,
    new_stable_id: Optional[str] = None,
) -> None:
    for decision in decisions:
        if decision.feature_type != feature_type or decision.action != action:
            continue
        if current_stable_id is not None and decision.current_stable_id != current_stable_id:
            continue
        if old_stable_id is not None and decision.old_stable_id != old_stable_id:
            continue
        if new_stable_id is not None and decision.new_stable_id != new_stable_id:
            continue
        return
    raise AssertionError(
        "Missing decision "
        f"type={feature_type} action={action} current={current_stable_id} "
        f"old={old_stable_id} new={new_stable_id}"
    )


def run_tests() -> None:
    old_g1 = "ENSNWIG00000000001"
    old_g2 = "ENSNWIG00000000002"
    old_g3 = "ENSNWIG00000000003"
    old_t1 = "ENSNWIT00000000001"
    old_t2 = "ENSNWIT00000000002"
    old_t3 = "ENSNWIT00000000003"
    target_g1 = "ENSNWIG00000090001"
    target_g2 = "ENSNWIG00000090002"
    target_g3 = "ENSNWIG00000090003"
    target_t1 = "ENSNWIT00000090001"
    target_t2 = "ENSNWIT00000090002"
    target_t3 = "ENSNWIT00000090003"

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        ref_gff = tmp / "ref.gff3"
        target_gff = tmp / "target.gff3"
        mapped_gff = tmp / "mapped.gff3"
        report = tmp / "report.txt"
        output_sql = tmp / "out.sql"
        output_tsv = tmp / "out.tsv"

        write_test_file(
            ref_gff,
            f"""
            ##gff-version 3
            chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:{old_g1}.1
            chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:{old_t1}.1;Parent=gene:{old_g1}.1
            chr1\ttest\tCDS\t150\t450\t.\t+\t0\tParent=transcript:{old_t1}.1
            chr1\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene:{old_g2}.1
            chr1\ttest\tmRNA\t1000\t1500\t.\t+\t.\tID=transcript:{old_t2}.1;Parent=gene:{old_g2}.1
            chr1\ttest\tgene\t2000\t2500\t.\t+\t.\tID=gene:{old_g3}.1
            chr1\ttest\tmRNA\t2000\t2500\t.\t+\t.\tID=transcript:{old_t3}.1;Parent=gene:{old_g3}.1
            chr1\ttest\tCDS\t2050\t2450\t.\t+\t0\tParent=transcript:{old_t3}.1
            """,
        )
        write_test_file(
            target_gff,
            f"""
            ##gff-version 3
            chrT\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:{target_g1}.1
            chrT\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:{target_t1}.1;Parent=gene:{target_g1}.1
            chrT\ttest\tCDS\t150\t450\t.\t+\t0\tParent=transcript:{target_t1}.1
            chrT\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene:{target_g2}.1
            chrT\ttest\tmRNA\t1000\t1500\t.\t+\t.\tID=transcript:{target_t2}.1;Parent=gene:{target_g2}.1
            chrT\ttest\tgene\t3000\t3500\t.\t+\t.\tID=gene:{target_g3}.1
            chrT\ttest\tmRNA\t3000\t3500\t.\t+\t.\tID=transcript:{target_t3}.1;Parent=gene:{target_g3}.1
            chrT\ttest\tCDS\t3050\t3450\t.\t+\t0\tParent=transcript:{target_t3}.1
            """,
        )
        write_test_file(
            mapped_gff,
            f"""
            ##gff-version 3
            chrT\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:{old_g1}.1
            chrT\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:{old_t1}.1;Parent=gene:{old_g1}.1
            chrT\ttest\texon\t100\t500\t.\t+\t.\tID=exon:{old_t1}.exon1;Parent=transcript:{old_t1}.1
            chrT\ttest\tCDS\t150\t450\t.\t+\t0\tParent=transcript:{old_t1}.1
            chrT\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene:{old_g3}.1
            chrT\ttest\tmRNA\t1000\t1500\t.\t+\t.\tID=transcript:{old_t3}.1;Parent=gene:{old_g3}.1
            chrT\ttest\tCDS\t1050\t1450\t.\t+\t0\tParent=transcript:{old_t3}.1
            """,
        )
        write_test_file(
            report,
            f"""
            Stable ID mapper report

            Missing gene IDs:
              gene:{old_g2}
            """,
        )

        args = argparse.Namespace(
            ref_gff=ref_gff,
            target_gff=target_gff,
            mapped_gff=mapped_gff,
            report=report,
            mapping_session_id=42,
            gene_range=parse_id_range("ENSNWIG:5000001-5000010"),
            transcript_range=parse_id_range("ENSNWIT:5000001-5000010"),
            translation_range=parse_id_range("ENSNWIP:5000001-5000010"),
            output_sql=output_sql,
            output_tsv=output_tsv,
            include_translations=True,
            dry_run=True,
            backup_prefix="stable_id_mapper_backup_test",
            batch_size=2,
            replace_events_for_session=False,
            min_overlap=0.10,
        )
        decisions = run_pipeline(args)

        expected_counts = {
            ("gene", "mapped"): 2,
            ("gene", "missing"): 1,
            ("gene", "new"): 1,
            ("transcript", "mapped"): 2,
            ("transcript", "missing"): 1,
            ("transcript", "new"): 1,
            ("translation", "mapped"): 1,
            ("translation", "missing"): 1,
            ("translation", "new"): 1,
        }
        observed_counts = Counter(
            (decision.feature_type, decision.action) for decision in decisions
        )
        if observed_counts != expected_counts:
            raise AssertionError(
                f"Unexpected decision counts: {observed_counts} != {expected_counts}"
            )

        assert_has_decision(
            decisions,
            "gene",
            "mapped",
            current_stable_id=target_g1,
            old_stable_id=old_g1,
            new_stable_id=old_g1,
        )
        assert_has_decision(
            decisions,
            "gene",
            "missing",
            old_stable_id=old_g2,
        )
        assert_has_decision(
            decisions,
            "gene",
            "new",
            current_stable_id=target_g3,
        )
        assert_has_decision(
            decisions,
            "transcript",
            "mapped",
            current_stable_id=target_t1,
            old_stable_id=old_t1,
            new_stable_id=old_t1,
        )
        assert_has_decision(
            decisions,
            "transcript",
            "missing",
            old_stable_id=old_t2,
        )
        assert_has_decision(
            decisions,
            "translation",
            "mapped",
            current_stable_id=target_t1,
            old_stable_id=old_t1,
            new_stable_id=old_t1,
        )
        assert_has_decision(
            decisions,
            "translation",
            "missing",
            old_stable_id=old_t3,
        )
        assert_has_decision(
            decisions,
            "translation",
            "new",
            current_stable_id=target_t3,
        )

        sql_text = output_sql.read_text(encoding="utf-8")
        if "CREATE TEMPORARY TABLE tmp_stable_id_mapper_decisions" not in sql_text:
            raise AssertionError("Dry-run decision table was not written")
        if "ROLLBACK;" not in sql_text:
            raise AssertionError("Dry-run SQL should end with ROLLBACK")
        tsv_header = output_tsv.read_text(encoding="utf-8").splitlines()[0]
        if not tsv_header.startswith("type\taction\tcurrent_stable_id"):
            raise AssertionError("TSV header is incorrect")

    sys.stderr.write("All tests passed.\n")


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Ensembl core SQL from main.py stable-ID mapper output."
    )
    parser.add_argument("--test", action="store_true", help="run built-in synthetic tests")
    parser.add_argument("--ref-gff", type=Path)
    parser.add_argument("--target-gff", type=Path)
    parser.add_argument("--mapped-gff", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--mapping-session-id", type=int)
    parser.add_argument("--gene-range", type=parse_id_range)
    parser.add_argument("--transcript-range", type=parse_id_range)
    parser.add_argument("--translation-range", type=parse_id_range)
    parser.add_argument("--output-sql", type=Path)
    parser.add_argument("--output-tsv", type=Path)
    parser.add_argument("--include-translations", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--backup-prefix",
        default=f"stable_id_mapper_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
    )
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--replace-events-for-session", action="store_true")
    parser.add_argument("--min-overlap", type=float, default=0.10)

    args = parser.parse_args(argv)
    if args.test:
        return args

    required = [
        ("ref_gff", "--ref-gff"),
        ("target_gff", "--target-gff"),
        ("mapped_gff", "--mapped-gff"),
        ("report", "--report"),
        ("mapping_session_id", "--mapping-session-id"),
        ("gene_range", "--gene-range"),
        ("transcript_range", "--transcript-range"),
        ("translation_range", "--translation-range"),
        ("output_sql", "--output-sql"),
    ]
    missing = [flag for attr, flag in required if getattr(args, attr) is None]
    if missing:
        parser.error("the following arguments are required: " + ", ".join(missing))
    if args.batch_size < 1:
        parser.error("--batch-size must be >= 1")
    if args.min_overlap < 0:
        parser.error("--min-overlap must be >= 0")
    return args


def main() -> None:
    args = parse_args()
    if args.test:
        run_tests()
        return
    run_pipeline(args)


if __name__ == "__main__":
    main()

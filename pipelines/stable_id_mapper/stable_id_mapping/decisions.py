"""Stable-ID decision rules."""

from __future__ import annotations

import sys
from collections import defaultdict
from typing import Iterable, Optional

from .ids import IdAllocator
from .models import Decision, Feature
from .scoring import ScoreEvidence, score_with_evidence


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


def old_version_for(
    features: dict[str, Feature],
    stable_id: str,
    fallback: Feature,
) -> int:
    old = features.get(stable_id)
    return old.version if old else fallback.version


def mapped_version(old_version: int) -> int:
    return old_version + 1


def build_gene_decisions(
    ref_features: dict[str, dict[str, Feature]],
    target_features: dict[str, dict[str, Feature]],
    mapped_features: dict[str, dict[str, Feature]],
    missing_gene_ids: set[str],
    allocators: dict[str, IdAllocator],
    mapping_session_id: int,
    min_overlap: float,
    score_evidence: ScoreEvidence,
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
                reason="gene listed as missing by mapping report",
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
        new_version = mapped_version(old_version)
        if target is None:
            if old_id in used_targets:
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
                        reason="mapped gene fallback target stable ID was already claimed",
                    )
                )
                continue

            used_targets.add(old_id)
            target_feature = target_features["gene"].get(old_id)
            if target_feature is not None:
                used_targets.add(target_feature.stable_id)
            old_gene_to_target_gene[old_id] = old_id
            decisions.append(
                Decision(
                    feature_type="gene",
                    action="mapped",
                    current_stable_id=old_id,
                    current_version=target_feature.version if target_feature else 0,
                    old_stable_id=old_id,
                    old_version=old_version,
                    new_stable_id=old_id,
                    new_version=new_version,
                    mapping_session_id=mapping_session_id,
                    score=0.0,
                    reason="mapped gene trusted after coordinate overlap failed",
                )
            )
            continue

        used_targets.add(target.stable_id)
        old_gene_to_target_gene[old_id] = target.stable_id
        score, score_source = score_with_evidence(
            score_evidence,
            "gene",
            old_id,
            target.stable_id,
            overlap_score(mapped, target),
            "coordinate_overlap",
        )
        decisions.append(
            Decision(
                feature_type="gene",
                action="mapped",
                current_stable_id=target.stable_id,
                current_version=target.version,
                old_stable_id=old_id,
                old_version=old_version,
                new_stable_id=old_id,
                new_version=new_version,
                mapping_session_id=mapping_session_id,
                score=score,
                reason=(
                    "mapped gene matched a target gene by coordinate overlap; "
                    f"score_source={score_source}"
                ),
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
    score_evidence: ScoreEvidence,
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
        new_version = mapped_version(old_version)
        if target is None:
            if old_id in used_targets:
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
                        reason="mapped transcript fallback target stable ID was already claimed",
                    )
                )
                continue

            used_targets.add(old_id)
            target_feature = target_features["transcript"].get(old_id)
            if target_feature is not None:
                used_targets.add(target_feature.stable_id)
            old_transcript_to_target_transcript[old_id] = old_id
            decisions.append(
                Decision(
                    feature_type="transcript",
                    action="mapped",
                    current_stable_id=old_id,
                    current_version=target_feature.version if target_feature else 0,
                    old_stable_id=old_id,
                    old_version=old_version,
                    new_stable_id=old_id,
                    new_version=new_version,
                    mapping_session_id=mapping_session_id,
                    score=0.0,
                    reason="mapped transcript trusted after coordinate overlap failed",
                )
            )
            continue

        used_targets.add(target.stable_id)
        old_transcript_to_target_transcript[old_id] = target.stable_id
        score, score_source = score_with_evidence(
            score_evidence,
            "transcript",
            old_id,
            target.stable_id,
            overlap_score(mapped, target),
            "coordinate_overlap",
        )
        decisions.append(
            Decision(
                feature_type="transcript",
                action="mapped",
                current_stable_id=target.stable_id,
                current_version=target.version,
                old_stable_id=old_id,
                old_version=old_version,
                new_stable_id=old_id,
                new_version=new_version,
                mapping_session_id=mapping_session_id,
                score=score,
                reason=(
                    "mapped transcript matched a target transcript by coordinate overlap; "
                    f"score_source={score_source}"
                ),
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
        new_version=mapped_version(old.version),
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
    score_evidence: ScoreEvidence,
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
            transcript_evidence = score_evidence.get(
                "transcript",
                old_transcript_id,
                target_transcript_id,
            )
            score = transcript_evidence.score if transcript_evidence else 1.0
            source = (
                "parent_transcript_lifton_score"
                if transcript_evidence
                else "single_translation_pair"
            )
            decisions.append(
                mapped_translation_decision(
                    old,
                    target,
                    mapping_session_id,
                    score,
                    "single old and target translation under a mapped transcript; "
                    f"score_source={source}",
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
            fallback_score = overlap_score(old_query, target)
            transcript_evidence = score_evidence.get(
                "transcript",
                old_transcript_id,
                target_transcript_id,
            )
            score = transcript_evidence.score if transcript_evidence else fallback_score
            source = (
                "parent_transcript_lifton_score"
                if transcript_evidence
                else "coordinate_overlap"
            )
            decisions.append(
                mapped_translation_decision(
                    old,
                    target,
                    mapping_session_id,
                    score,
                    "old translation matched a target translation by coordinate overlap; "
                    f"score_source={source}",
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
    score_evidence: Optional[ScoreEvidence] = None,
) -> list[Decision]:
    if score_evidence is None:
        score_evidence = ScoreEvidence()
    gene_decisions, old_gene_to_target_gene = build_gene_decisions(
        ref_features,
        target_features,
        mapped_features,
        missing_gene_ids,
        allocators,
        mapping_session_id,
        min_overlap,
        score_evidence,
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
        score_evidence,
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
            score_evidence,
        )
        if include_translations
        else []
    )
    return gene_decisions + transcript_decisions + translation_decisions

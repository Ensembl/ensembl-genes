"""File-based stable-ID event pipeline."""

from __future__ import annotations

import sys
from collections import Counter

from .config import StableIdEventConfig
from .decisions import build_decisions
from .gff3 import parse_gff3, parse_missing_gene_ids
from .ids import IdAllocator, collect_reserved_ids, make_allocator
from .models import ACTION_ORDER, FEATURE_ORDER, Decision
from .outputs import write_sql, write_tsv
from .scoring import write_score_evidence_tsv


def make_allocators(
    config: StableIdEventConfig,
    reserved_ids: set[str],
) -> dict[str, IdAllocator]:
    return {
        "gene": make_allocator(config.gene_range, reserved_ids),
        "transcript": make_allocator(config.transcript_range, reserved_ids),
        "translation": make_allocator(config.translation_range, reserved_ids),
    }


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
    evidence_counts = Counter(
        decision.feature_type
        for decision in decisions
        if decision.feature_type in ("gene", "transcript")
        and decision.action == "mapped"
        and "by lifton structural evidence" in decision.reason
    )
    overlap_counts = Counter(
        decision.feature_type
        for decision in decisions
        if decision.feature_type in ("gene", "transcript")
        and decision.action == "mapped"
        and "matched a target" in decision.reason
        and "by coordinate overlap" in decision.reason
    )
    sys.stderr.write(
        "Mapped by LiftOn structural evidence: "
        f"{evidence_counts['gene']} genes, {evidence_counts['transcript']} transcripts; "
        "mapped by coordinate overlap: "
        f"{overlap_counts['gene']} genes, {overlap_counts['transcript']} transcripts\n"
    )


def run_pipeline(config: StableIdEventConfig) -> list[Decision]:
    config.validate()
    ref_features = parse_gff3(config.ref_gff)
    target_features = parse_gff3(config.target_gff)
    mapped_features = parse_gff3(config.mapped_gff)
    missing_gene_ids = parse_missing_gene_ids(config.report)
    reserved_ids = collect_reserved_ids(ref_features, target_features, mapped_features)
    allocators = make_allocators(config, reserved_ids)
    decisions = build_decisions(
        ref_features,
        target_features,
        mapped_features,
        missing_gene_ids,
        allocators,
        config.mapping_session_id,
        config.min_overlap,
        config.include_translations,
        config.score_evidence,
    )
    write_sql(decisions, config.output_sql, config)
    if config.output_tsv:
        write_tsv(decisions, config.output_tsv)
    if config.output_score_evidence_tsv:
        write_score_evidence_tsv(config.score_evidence, config.output_score_evidence_tsv)
    print_summary(decisions)
    return decisions

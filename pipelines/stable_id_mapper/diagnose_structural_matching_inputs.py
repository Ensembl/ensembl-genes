#!/usr/bin/env python3
"""Diagnose GFF3 transcript-structure linkage for stable-ID structural matching."""

from __future__ import annotations

import argparse
import gzip
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from stable_id_mapping.gff3 import parse_attrs, split_stable_id

TRANSCRIPT_TYPES = {"mrna", "transcript"}
CHILD_TYPES = {"exon", "cds"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize whether exon/CDS Parent attributes can be linked to "
            "explicit mRNA/transcript IDs."
        )
    )
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--label", default="annotation")
    parser.add_argument("--example-limit", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = summarize_gff3_structure(args.gff)
    print_summary(summary, args.label, args.example_limit)


def summarize_gff3_structure(path: Path) -> dict[str, object]:
    feature_counts: Counter[str] = Counter()
    transcript_ids: set[str] = set()
    transcript_core_to_ids: dict[str, set[str]] = defaultdict(set)
    child_parent_links = 0
    child_rows_without_parent = 0
    exact_child_links: Counter[str] = Counter()
    core_child_links: Counter[str] = Counter()
    unmatched_child_parent_links: Counter[str] = Counter()
    parent_prefix_counts: Counter[str] = Counter()
    examples: dict[str, list[str]] = {
        "exact": [],
        "core_only": [],
        "unmatched": [],
        "no_parent": [],
    }

    rows = list(iter_gff3_rows(path))
    for fields in rows:
        _seqid, _source, feature_type, _start, _end, _score, _strand, _phase, attrs_text = fields
        feature_type_lc = feature_type.lower()
        feature_counts[feature_type_lc] += 1
        attrs = parse_attrs(attrs_text)
        if feature_type_lc not in TRANSCRIPT_TYPES:
            continue
        transcript_id = attrs.get("ID") or attrs.get("transcript_id") or attrs.get("Name")
        if not transcript_id:
            continue
        transcript_ids.add(transcript_id)
        core, _version = split_stable_id(transcript_id)
        if core:
            transcript_core_to_ids[core].add(transcript_id)

    transcripts_with_exact_children: set[str] = set()
    transcripts_with_core_children: set[str] = set()
    for fields in rows:
        seqid, _source, feature_type, start, end, _score, _strand, _phase, attrs_text = fields
        feature_type_lc = feature_type.lower()
        if feature_type_lc not in CHILD_TYPES:
            continue
        attrs = parse_attrs(attrs_text)
        parent_text = attrs.get("Parent")
        if not parent_text:
            child_rows_without_parent += 1
            add_example(
                examples["no_parent"],
                f"{feature_type}:{seqid}:{start}-{end}",
            )
            continue
        for parent_id in parent_text.split(","):
            child_parent_links += 1
            parent_id = parent_id.strip()
            parent_prefix_counts[parent_prefix(parent_id)] += 1
            if parent_id in transcript_ids:
                exact_child_links[feature_type_lc] += 1
                core_child_links[feature_type_lc] += 1
                transcripts_with_exact_children.add(parent_id)
                transcripts_with_core_children.add(parent_id)
                add_example(
                    examples["exact"],
                    f"{feature_type}:{seqid}:{start}-{end} Parent={parent_id}",
                )
                continue
            parent_core, _version = split_stable_id(parent_id)
            core_matches = transcript_core_to_ids.get(parent_core or "")
            if core_matches:
                core_child_links[feature_type_lc] += 1
                transcripts_with_core_children.update(core_matches)
                add_example(
                    examples["core_only"],
                    f"{feature_type}:{seqid}:{start}-{end} Parent={parent_id} "
                    f"matches_transcript_core={parent_core}",
                )
                continue
            unmatched_child_parent_links[feature_type_lc] += 1
            add_example(
                examples["unmatched"],
                f"{feature_type}:{seqid}:{start}-{end} Parent={parent_id}",
            )

    return {
        "feature_counts": feature_counts,
        "transcript_count": len(transcript_ids),
        "child_parent_links": child_parent_links,
        "child_rows_without_parent": child_rows_without_parent,
        "exact_child_links": exact_child_links,
        "core_child_links": core_child_links,
        "unmatched_child_parent_links": unmatched_child_parent_links,
        "parent_prefix_counts": parent_prefix_counts,
        "transcripts_with_exact_children": len(transcripts_with_exact_children),
        "transcripts_with_core_children": len(transcripts_with_core_children),
        "examples": examples,
    }


def iter_gff3_rows(path: Path) -> Iterable[list[str]]:
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9:
                yield fields


def add_example(examples: list[str], example: str, limit: int = 20) -> None:
    if len(examples) < limit:
        examples.append(example)


def parent_prefix(parent_id: str) -> str:
    if ":" not in parent_id:
        return "<none>"
    return parent_id.split(":", 1)[0]


def print_summary(summary: dict[str, object], label: str, example_limit: int) -> None:
    feature_counts = summary["feature_counts"]
    assert isinstance(feature_counts, Counter)
    print(f"{label} feature counts:")
    print(
        "  "
        + ", ".join(
            f"{feature_type}={count}"
            for feature_type, count in feature_counts.most_common(15)
        )
    )
    print(f"{label} explicit transcripts: {summary['transcript_count']}")
    print(f"{label} child Parent links: {summary['child_parent_links']}")
    print(f"{label} child rows without Parent: {summary['child_rows_without_parent']}")
    print(
        f"{label} transcripts with exact exon/CDS children: "
        f"{summary['transcripts_with_exact_children']}"
    )
    print(
        f"{label} transcripts with core-compatible exon/CDS children: "
        f"{summary['transcripts_with_core_children']}"
    )
    print_counter("Exact child links", summary["exact_child_links"])
    print_counter("Core-compatible child links", summary["core_child_links"])
    print_counter("Unmatched child Parent links", summary["unmatched_child_parent_links"])
    print_counter("Child Parent prefixes", summary["parent_prefix_counts"])
    examples = summary["examples"]
    assert isinstance(examples, dict)
    for key in ("exact", "core_only", "unmatched", "no_parent"):
        values = examples[key]
        print(f"{key} examples:")
        if not values:
            print("  none")
            continue
        for value in values[:example_limit]:
            print(f"  {value}")


def print_counter(label: str, counter_obj: object) -> None:
    assert isinstance(counter_obj, Counter)
    if not counter_obj:
        print(f"{label}: none")
        return
    print(
        f"{label}: "
        + ", ".join(f"{key}={value}" for key, value in counter_obj.most_common())
    )


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open(encoding="utf-8")


if __name__ == "__main__":
    main()

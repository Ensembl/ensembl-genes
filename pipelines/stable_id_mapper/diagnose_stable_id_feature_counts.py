#!/usr/bin/env python3
"""Report the feature universe seen by the stable-ID decision parser."""

from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, Optional

from stable_id_mapping.gff3 import parse_attrs, parse_gff3, split_stable_id

BIOTYPE_KEYS = (
    "biotype",
    "gene_biotype",
    "transcript_biotype",
    "gbkey",
    "gene_type",
    "transcript_type",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Show raw GFF3 counts and the gene/transcript/translation counts "
            "that the stable-ID mapper will use."
        )
    )
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--label", default="annotation")
    parser.add_argument("--example-limit", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = summarize_feature_counts(args.gff)
    print_summary(summary, args.label, args.example_limit)


def summarize_feature_counts(path: Path) -> dict[str, object]:
    parsed_features = parse_gff3(path)
    raw_feature_types: Counter[str] = Counter()
    parentless_gene_rows = 0
    parented_gene_rows = 0
    gene_rows_with_parseable_id = 0
    gene_id_prefixes: Counter[str] = Counter()
    gene_stable_prefixes: Counter[str] = Counter()
    gene_biotypes: Counter[str] = Counter()
    duplicate_ids: Counter[str] = Counter()
    seen_ids: Counter[str] = Counter()
    examples: dict[str, list[str]] = {
        "gene_rows": [],
        "parented_gene_rows": [],
        "unparseable_gene_ids": [],
    }

    for fields in iter_gff3_rows(path):
        seqid, _source, feature_type, start, end, _score, strand, _phase, attrs_text = fields
        feature_type_lc = feature_type.lower()
        raw_feature_types[feature_type_lc] += 1
        attrs = parse_attrs(attrs_text)
        feature_id = attrs.get("ID")
        if feature_id:
            seen_ids[feature_id] += 1
        if feature_type_lc != "gene":
            continue

        parent = attrs.get("Parent")
        if parent:
            parented_gene_rows += 1
            add_example(
                examples["parented_gene_rows"],
                f"{seqid}:{start}-{end}({strand}) ID={feature_id or '<missing>'} Parent={parent}",
            )
        else:
            parentless_gene_rows += 1
            add_example(
                examples["gene_rows"],
                f"{seqid}:{start}-{end}({strand}) ID={feature_id or '<missing>'} "
                f"biotype={first_attr(attrs, BIOTYPE_KEYS) or '<missing>'}",
            )

        if feature_id:
            gene_id_prefixes[prefix_before_colon(feature_id)] += 1
        stable_id, _version = split_stable_id(feature_id)
        if stable_id:
            gene_rows_with_parseable_id += 1
            gene_stable_prefixes[stable_prefix(stable_id)] += 1
        else:
            add_example(
                examples["unparseable_gene_ids"],
                f"{seqid}:{start}-{end}({strand}) ID={feature_id or '<missing>'}",
            )
        gene_biotypes[first_attr(attrs, BIOTYPE_KEYS) or "<missing>"] += 1

    duplicate_ids = Counter({feature_id: count for feature_id, count in seen_ids.items() if count > 1})
    return {
        "parsed_counts": {
            feature_type: len(features)
            for feature_type, features in parsed_features.items()
        },
        "raw_feature_types": raw_feature_types,
        "parentless_gene_rows": parentless_gene_rows,
        "parented_gene_rows": parented_gene_rows,
        "gene_rows_with_parseable_id": gene_rows_with_parseable_id,
        "gene_id_prefixes": gene_id_prefixes,
        "gene_stable_prefixes": gene_stable_prefixes,
        "gene_biotypes": gene_biotypes,
        "duplicate_ids": duplicate_ids,
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


def first_attr(attrs: dict[str, str], keys: Iterable[str]) -> Optional[str]:
    for key in keys:
        value = attrs.get(key)
        if value:
            return value
    return None


def prefix_before_colon(value: str) -> str:
    if ":" in value:
        return value.split(":", 1)[0]
    return "<none>"


def stable_prefix(stable_id: str) -> str:
    match = re.match(r"([A-Za-z_]+)", stable_id)
    return match.group(1) if match else "<none>"


def add_example(examples: list[str], example: str, limit: int = 20) -> None:
    if len(examples) < limit:
        examples.append(example)


def print_summary(summary: dict[str, object], label: str, example_limit: int) -> None:
    parsed_counts = summary["parsed_counts"]
    assert isinstance(parsed_counts, dict)
    print(f"{label} stable-ID parser counts:")
    for feature_type in ("gene", "transcript", "translation"):
        print(f"  {feature_type}: {parsed_counts.get(feature_type, 0)}")

    print_counter(f"{label} raw feature_type counts", summary["raw_feature_types"], limit=20)
    print(f"{label} parentless gene rows: {summary['parentless_gene_rows']}")
    print(f"{label} parented gene rows: {summary['parented_gene_rows']}")
    print(f"{label} gene rows with parseable ID: {summary['gene_rows_with_parseable_id']}")
    print_counter(f"{label} gene ID prefixes", summary["gene_id_prefixes"])
    print_counter(f"{label} gene stable-ID prefixes", summary["gene_stable_prefixes"])
    print_counter(f"{label} gene biotypes", summary["gene_biotypes"], limit=20)
    print_counter(f"{label} duplicate raw IDs", summary["duplicate_ids"], limit=20)

    examples = summary["examples"]
    assert isinstance(examples, dict)
    for key in ("gene_rows", "parented_gene_rows", "unparseable_gene_ids"):
        print(f"{label} {key} examples:")
        values = examples[key]
        if not values:
            print("  none")
            continue
        for value in values[:example_limit]:
            print(f"  {value}")


def print_counter(label: str, counter_obj: object, limit: int = 10) -> None:
    assert isinstance(counter_obj, Counter)
    if not counter_obj:
        print(f"{label}: none")
        return
    print(f"{label}:")
    for key, value in counter_obj.most_common(limit):
        print(f"  {key}: {value}")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open(encoding="utf-8")


if __name__ == "__main__":
    main()

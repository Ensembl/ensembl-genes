#!/usr/bin/env python3
"""Pick representative examples from a completed stable-ID mapping run."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from stable_id_mapping.gff3 import split_stable_id


@dataclass(frozen=True)
class Example:
    category: str
    feature_type: str
    old_stable_id: str
    target_stable_id: str
    score: str
    locus_score: str
    note: str
    source: str


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select representative examples from stable-ID mapper outputs."
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed stable-ID mapper output directory",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=5,
        help="Maximum examples per category",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output TSV path. Defaults to RUN_DIR/reports/real_example_candidates.tsv",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    output = args.output or args.run_dir / "reports" / "real_example_candidates.tsv"
    examples = pick_examples(args.run_dir, args.limit)
    write_examples(examples, output)
    print(f"Wrote {len(examples)} examples: {output}")
    for category in sorted({example.category for example in examples}):
        count = sum(1 for example in examples if example.category == category)
        print(f"  {category}: {count}")


def pick_examples(run_dir: Path, limit: int) -> list[Example]:
    decisions = _read_tsv(run_dir / "stable_id_decisions.tsv")
    locus_path = run_dir / "matching" / "lifton.gene_locus_comparison.tsv"
    locus_rows = _read_tsv(locus_path) if locus_path.exists() else []
    examples: list[Example] = []
    examples.extend(_decision_examples(decisions, "gene_structural_mapped", "gene", "lifton structural evidence", limit))
    examples.extend(_decision_examples(decisions, "gene_coordinate_mapped", "gene", "coordinate overlap", limit))
    examples.extend(_decision_examples(decisions, "transcript_structural_mapped", "transcript", "lifton structural evidence", limit))
    examples.extend(_locus_examples(locus_rows, "high_locus_no_accepted_structure", "no_accepted_structure", limit))
    examples.extend(_locus_examples(locus_rows, "locus_structure_disagree", "different", limit))
    examples.extend(_no_locus_examples(locus_rows, decisions, limit))
    return examples


def write_examples(examples: list[Example], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "category",
                "feature_type",
                "old_stable_id",
                "target_stable_id",
                "score",
                "locus_score",
                "note",
                "source",
            ]
        )
        for example in examples:
            writer.writerow(
                [
                    example.category,
                    example.feature_type,
                    example.old_stable_id,
                    example.target_stable_id,
                    example.score,
                    example.locus_score,
                    example.note,
                    example.source,
                ]
            )


def _decision_examples(
    decisions: list[dict[str, str]],
    category: str,
    feature_type: str,
    reason_text: str,
    limit: int,
) -> list[Example]:
    rows = [
        row
        for row in decisions
        if row.get("type") == feature_type
        and row.get("action") == "mapped"
        and reason_text in row.get("reason", "")
    ]
    rows.sort(key=lambda row: (-_float(row.get("score")), row.get("old_stable_id", "")))
    return [
        Example(
            category=category,
            feature_type=feature_type,
            old_stable_id=row.get("old_stable_id", ""),
            target_stable_id=row.get("current_stable_id", ""),
            score=row.get("score", ""),
            locus_score="",
            note=row.get("reason", ""),
            source="stable_id_decisions.tsv",
        )
        for row in rows[:limit]
    ]


def _locus_examples(
    locus_rows: list[dict[str, str]],
    category: str,
    locus_status: str,
    limit: int,
) -> list[Example]:
    rows = [
        row
        for row in locus_rows
        if row.get("locus_vs_structure") == locus_status
        and _float(row.get("locus_score")) > 0
    ]
    rows.sort(key=lambda row: (-_float(row.get("locus_score")), row.get("lifton_gene", "")))
    return [
        Example(
            category=category,
            feature_type="gene",
            old_stable_id=_core(row.get("lifton_gene", "")),
            target_stable_id=_core(
                row.get("target_gene_by_locus", "")
                or row.get("structure_accepted_target_gene", "")
            ),
            score=row.get("best_tx_structure_score", "") or row.get("structure_weighted_score", ""),
            locus_score=row.get("locus_score", ""),
            note=(
                f"locus_vs_structure={row.get('locus_vs_structure', '')};"
                f"old_cov={row.get('old_gene_coverage', '')};"
                f"target_cov={row.get('target_gene_coverage', '')};"
                f"structure_target={row.get('structure_accepted_target_gene', '')}"
            ),
            source="matching/lifton.gene_locus_comparison.tsv",
        )
        for row in rows[:limit]
    ]


def _no_locus_examples(
    locus_rows: list[dict[str, str]],
    decisions: list[dict[str, str]],
    limit: int,
) -> list[Example]:
    missing_old_genes = {
        row.get("old_stable_id", "")
        for row in decisions
        if row.get("type") == "gene" and row.get("action") == "missing"
    }
    rows = [
        row
        for row in locus_rows
        if row.get("locus_vs_structure") == "no_locus_candidate"
        and _core(row.get("lifton_gene", "")) in missing_old_genes
    ]
    rows.sort(key=lambda row: row.get("lifton_gene", ""))
    return [
        Example(
            category="missing_projected_gene_no_locus",
            feature_type="gene",
            old_stable_id=_core(row.get("lifton_gene", "")),
            target_stable_id="",
            score="",
            locus_score=row.get("locus_score", ""),
            note=(
                f"projected={row.get('lifton_contig', '')}:"
                f"{row.get('lifton_start', '')}-{row.get('lifton_end', '')}"
                f"({row.get('lifton_strand', '')})"
            ),
            source="matching/lifton.gene_locus_comparison.tsv",
        )
        for row in rows[:limit]
    ]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _float(value: Optional[str]) -> float:
    if value is None or value == "":
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def _core(value: str) -> str:
    if not value:
        return ""
    stable_id, _version = split_stable_id(value)
    return stable_id or value


if __name__ == "__main__":
    main()

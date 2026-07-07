"""Rules and thresholds for stable-ID mapping."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


DEFAULT_RULES_PATH = Path(__file__).resolve().parents[1] / "stable_id_mapping_rules.json"


@dataclass(frozen=True)
class CoordinateOverlapRules:
    min_overlap: float


@dataclass(frozen=True)
class StructuralMatchingRules:
    window: int
    topk: int
    min_score: float
    good_score: float
    confident_score: float
    gene_fraction: float
    score_weights: dict[str, float]


@dataclass(frozen=True)
class MappingRules:
    coordinate_overlap: CoordinateOverlapRules
    structural_matching: StructuralMatchingRules
    source_path: Path


def load_mapping_rules(path: Optional[Path] = None) -> MappingRules:
    source_path = path or DEFAULT_RULES_PATH
    with source_path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    coordinate = data.get("coordinate_overlap", {})
    structural = data.get("structural_matching", {})
    score_weights = structural.get("score_weights", {})

    rules = MappingRules(
        coordinate_overlap=CoordinateOverlapRules(
            min_overlap=_float_rule(coordinate, "min_overlap"),
        ),
        structural_matching=StructuralMatchingRules(
            window=_int_rule(structural, "window"),
            topk=_int_rule(structural, "topk"),
            min_score=_float_rule(structural, "min_score"),
            good_score=_float_rule(structural, "good_score"),
            confident_score=_float_rule(structural, "confident_score"),
            gene_fraction=_float_rule(structural, "gene_fraction"),
            score_weights={
                key: float(value)
                for key, value in score_weights.items()
            },
        ),
        source_path=source_path,
    )
    _validate_rules(rules)
    return rules


def default_score_weights() -> dict[str, float]:
    return dict(load_mapping_rules().structural_matching.score_weights)


def _float_rule(data: dict, key: str) -> float:
    if key not in data:
        raise ValueError(f"Rules config is missing {key!r}")
    return float(data[key])


def _int_rule(data: dict, key: str) -> int:
    if key not in data:
        raise ValueError(f"Rules config is missing {key!r}")
    return int(data[key])


def _validate_rules(rules: MappingRules) -> None:
    if rules.coordinate_overlap.min_overlap < 0:
        raise ValueError("coordinate_overlap.min_overlap must be >= 0")
    structural = rules.structural_matching
    if structural.window < 0:
        raise ValueError("structural_matching.window must be >= 0")
    if structural.topk < 1:
        raise ValueError("structural_matching.topk must be >= 1")
    if structural.min_score < 0:
        raise ValueError("structural_matching.min_score must be >= 0")
    if structural.good_score < 0:
        raise ValueError("structural_matching.good_score must be >= 0")
    if structural.confident_score < 0:
        raise ValueError("structural_matching.confident_score must be >= 0")
    if structural.gene_fraction < 0:
        raise ValueError("structural_matching.gene_fraction must be >= 0")
    if not structural.score_weights:
        raise ValueError("structural_matching.score_weights must not be empty")

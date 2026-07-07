from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stable_id_mapping.rules import load_mapping_rules


def test_load_mapping_rules_from_json(tmp_path: Path) -> None:
    rules_path = tmp_path / "rules.json"
    rules_path.write_text(
        json.dumps(
            {
                "coordinate_overlap": {"min_overlap": 0.2},
                "structural_matching": {
                    "window": 42,
                    "topk": 3,
                    "min_score": 0.4,
                    "good_score": 0.6,
                    "confident_score": 0.8,
                    "gene_fraction": 0.7,
                    "score_weights": {
                        "query_coverage": 2,
                        "span_containment": 1,
                    },
                },
            }
        ),
        encoding="utf-8",
    )

    rules = load_mapping_rules(rules_path)

    assert rules.coordinate_overlap.min_overlap == 0.2
    assert rules.structural_matching.window == 42
    assert rules.structural_matching.topk == 3
    assert rules.structural_matching.min_score == 0.4
    assert rules.structural_matching.good_score == 0.6
    assert rules.structural_matching.confident_score == 0.8
    assert rules.structural_matching.gene_fraction == 0.7
    assert rules.structural_matching.score_weights["query_coverage"] == 2
    assert rules.source_path == rules_path

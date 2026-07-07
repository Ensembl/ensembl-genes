from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stable_id_mapping.ids import parse_id_range
from stable_id_mapping.synthetic_cases import write_synthetic_cases
from stable_id_mapping.workflow import SingleSpeciesRunConfig, run_single_species_pipeline


def decision_counts(result) -> dict[str, dict[str, int]]:
    counts = Counter(
        (decision.feature_type, decision.action)
        for decision in result.decisions
    )
    return {
        feature_type: {
            action: counts[(feature_type, action)]
            for action in ("mapped", "missing", "new")
        }
        for feature_type in ("gene", "transcript", "translation")
    }


def run_synthetic_case(case):
    return run_single_species_pipeline(
        SingleSpeciesRunConfig(
            ref_fasta=case.ref_fasta,
            ref_gff=case.ref_gff,
            target_fasta=case.target_fasta,
            target_gff=case.target_gff,
            db_name="synthetic_core_test",
            mapping_session_id=1,
            gene_range=parse_id_range("ENSFAKEG:990001-990999"),
            transcript_range=parse_id_range("ENSFAKET:990001-990999"),
            translation_range=parse_id_range("ENSFAKEP:990001-990999"),
            output_dir=case.case_dir / "run",
            existing_lifton_gff=case.lifton_gff,
            include_translations=True,
            dry_run_sql=True,
        )
    )


def test_high_identity_synthetic_case_maps_all_features(tmp_path: Path) -> None:
    case = write_synthetic_cases(tmp_path, ["high_identity"])[0]

    result = run_synthetic_case(case)

    assert decision_counts(result) == case.expected_decisions
    assert result.match_summary.gene_pairs == 3
    assert result.match_summary.transcript_pairs == 3


def test_isoform_shift_synthetic_case_maps_gene_only(tmp_path: Path) -> None:
    case = write_synthetic_cases(tmp_path, ["isoform_shift_gene_only"])[0]

    result = run_synthetic_case(case)

    assert decision_counts(result) == case.expected_decisions
    comparison_text = result.match_summary.gene_locus_comparison_path.read_text(
        encoding="utf-8"
    )
    assert "no_accepted_structure" in comparison_text
    assert "\t1.000000\t" in comparison_text


def test_unrelated_empty_lifton_synthetic_case_maps_nothing(tmp_path: Path) -> None:
    case = write_synthetic_cases(tmp_path, ["unrelated_empty_lifton"])[0]

    result = run_synthetic_case(case)

    assert decision_counts(result) == case.expected_decisions
    assert result.match_summary.gene_pairs == 0
    assert result.match_summary.transcript_pairs == 0

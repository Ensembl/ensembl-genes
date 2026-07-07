from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stable_id_mapping.scoring import (
    confidence_label,
    load_lifton_score_evidence,
    score_with_evidence,
)


def write_text(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def test_load_lifton_score_evidence_normalizes_ids_and_gene_score(tmp_path: Path) -> None:
    transcript_pairs = tmp_path / "pairs.transcript_pairs.tsv"
    gene_pairs = tmp_path / "pairs.gene_pairs.tsv"
    write_text(
        transcript_pairs,
        """
lifton_tx	ref_tx	score	status	intron_sim	jacc_internal	jacc_all	exon_count_sim	boundary_sim	lifton_identity_prior	lifton_gene	ref_gene	contig	strand	lifton_exons	ref_exons
transcript:OLDT.7	transcript:TARGETT.2	0.910000	confident	1	1	1	1	1	0	gene:OLDG.2	gene:TARGETG.1	chr1	+	3	3
        """,
    )
    write_text(
        gene_pairs,
        """
lifton_gene	ref_gene	weighted_score	fraction_of_total	n_transcripts
gene:OLDG.2	gene:TARGETG.1	1.800000	0.500000	2
        """,
    )

    evidence = load_lifton_score_evidence(transcript_pairs, gene_pairs)

    assert evidence.get("transcript", "OLDT", "TARGETT").score == 0.91
    assert evidence.get("gene", "OLDG", "TARGETG").score == 0.45
    score, source = score_with_evidence(
        evidence,
        "transcript",
        "OLDT",
        "TARGETT",
        fallback_score=0.1,
        fallback_source="fallback",
    )
    assert score == 0.91
    assert "lifton_transcript_structure" in source


def test_score_with_evidence_uses_fallback_when_pair_missing(tmp_path: Path) -> None:
    evidence = load_lifton_score_evidence(
        empty_tsv(tmp_path / "transcripts.tsv", "lifton_tx\tref_tx\tscore\n"),
        empty_tsv(
            tmp_path / "genes.tsv",
            "lifton_gene\tref_gene\tweighted_score\tfraction_of_total\tn_transcripts\n",
        ),
    )

    score, source = score_with_evidence(
        evidence,
        "gene",
        "missing_old",
        "missing_target",
        fallback_score=0.42,
        fallback_source="coordinate_overlap",
    )

    assert score == 0.42
    assert source == "coordinate_overlap"


def test_confidence_label_boundaries() -> None:
    assert confidence_label(0.90) == "high"
    assert confidence_label(0.80) == "medium"
    assert confidence_label(0.65) == "low"
    assert confidence_label(0.20) == "review"


def empty_tsv(path: Path, header: str) -> Path:
    path.write_text(header, encoding="utf-8")
    return path


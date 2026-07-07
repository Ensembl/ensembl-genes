"""Score evidence parsing and confidence labels."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .gff3 import split_stable_id


@dataclass(frozen=True)
class MappingEvidence:
    feature_type: str
    old_stable_id: str
    target_stable_id: str
    score: float
    source: str
    confidence: str
    detail: str = ""


@dataclass(frozen=True)
class ScoreEvidence:
    by_feature_type: dict[str, dict[tuple[str, str], MappingEvidence]] = field(
        default_factory=lambda: {"gene": {}, "transcript": {}, "translation": {}}
    )

    def get(
        self,
        feature_type: str,
        old_stable_id: Optional[str],
        target_stable_id: Optional[str],
    ) -> Optional[MappingEvidence]:
        if old_stable_id is None or target_stable_id is None:
            return None
        return self.by_feature_type.get(feature_type, {}).get(
            (old_stable_id, target_stable_id)
        )

    def add(self, evidence: MappingEvidence) -> None:
        self.by_feature_type.setdefault(evidence.feature_type, {})[
            (evidence.old_stable_id, evidence.target_stable_id)
        ] = evidence

    def all(self) -> list[MappingEvidence]:
        out: list[MappingEvidence] = []
        for feature_scores in self.by_feature_type.values():
            out.extend(feature_scores.values())
        return sorted(
            out,
            key=lambda evidence: (
                evidence.feature_type,
                evidence.old_stable_id,
                evidence.target_stable_id,
            ),
        )


def confidence_label(score: float) -> str:
    if score >= 0.85:
        return "high"
    if score >= 0.75:
        return "medium"
    if score >= 0.60:
        return "low"
    return "review"


def clamp_score(score: float) -> float:
    return max(0.0, min(1.0, score))


def core_stable_id(value: str) -> str:
    stable_id, _version = split_stable_id(value)
    if stable_id is None:
        raise ValueError(f"Could not parse stable ID from {value!r}")
    return stable_id


def load_lifton_score_evidence(
    transcript_pairs_path: str | Path,
    gene_pairs_path: str | Path,
) -> ScoreEvidence:
    evidence = ScoreEvidence()
    load_transcript_scores(transcript_pairs_path, evidence)
    load_gene_scores(gene_pairs_path, evidence)
    return evidence


def load_transcript_scores(
    transcript_pairs_path: str | Path,
    evidence: ScoreEvidence,
) -> None:
    with Path(transcript_pairs_path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            old_id = core_stable_id(row["lifton_tx"])
            target_id = core_stable_id(row["ref_tx"])
            score = clamp_score(float(row["score"]))
            status = row.get("status") or confidence_label(score)
            detail_parts = [
                f"{key}={row[key]}"
                for key in (
                    "intron_sim",
                    "jacc_internal",
                    "jacc_all",
                    "exon_count_sim",
                    "boundary_sim",
                    "lifton_identity_prior",
                )
                if key in row
            ]
            evidence.add(
                MappingEvidence(
                    feature_type="transcript",
                    old_stable_id=old_id,
                    target_stable_id=target_id,
                    score=score,
                    source="lifton_transcript_structure",
                    confidence=status,
                    detail=";".join(detail_parts),
                )
            )


def load_gene_scores(
    gene_pairs_path: str | Path,
    evidence: ScoreEvidence,
) -> None:
    with Path(gene_pairs_path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            old_id = core_stable_id(row["lifton_gene"])
            target_id = core_stable_id(row["ref_gene"])
            weighted_score = float(row["weighted_score"])
            fraction = float(row["fraction_of_total"])
            n_transcripts = int(row["n_transcripts"])
            average_score = weighted_score / n_transcripts if n_transcripts else 0.0
            score = clamp_score(average_score * fraction)
            evidence.add(
                MappingEvidence(
                    feature_type="gene",
                    old_stable_id=old_id,
                    target_stable_id=target_id,
                    score=score,
                    source="lifton_gene_aggregate",
                    confidence=confidence_label(score),
                    detail=(
                        f"weighted_score={weighted_score:.6f};"
                        f"fraction_of_total={fraction:.6f};"
                        f"n_transcripts={n_transcripts}"
                    ),
                )
            )


def score_with_evidence(
    evidence: ScoreEvidence,
    feature_type: str,
    old_stable_id: str,
    target_stable_id: str,
    fallback_score: float,
    fallback_source: str,
) -> tuple[float, str]:
    matched_evidence = evidence.get(feature_type, old_stable_id, target_stable_id)
    if matched_evidence is None:
        return fallback_score, fallback_source
    return matched_evidence.score, (
        f"{matched_evidence.source} confidence={matched_evidence.confidence}"
    )


def write_score_evidence_tsv(
    evidence: ScoreEvidence,
    path: str | Path,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "type",
                "old_stable_id",
                "target_stable_id",
                "score",
                "confidence",
                "source",
                "detail",
            ]
        )
        for item in evidence.all():
            writer.writerow(
                [
                    item.feature_type,
                    item.old_stable_id,
                    item.target_stable_id,
                    item.score,
                    item.confidence,
                    item.source,
                    item.detail,
                ]
            )


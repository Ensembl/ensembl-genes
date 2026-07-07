"""Run structural matching between LiftOn output and target annotation."""

from __future__ import annotations

import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path

import lifton_id_mapper


@dataclass(frozen=True)
class LiftonMatchConfig:
    lifton_gff: Path
    target_gff: Path
    out_prefix: Path
    window: int = 100000
    topk: int = 5
    min_score: float = 0.60
    good: float = 0.75
    confident: float = 0.85
    gene_fraction: float = 0.60

    @property
    def transcript_pairs_path(self) -> Path:
        return Path(f"{self.out_prefix}.transcript_pairs.tsv")

    @property
    def gene_pairs_path(self) -> Path:
        return Path(f"{self.out_prefix}.gene_pairs.tsv")

    def validate(self) -> None:
        for path in (self.lifton_gff, self.target_gff):
            if not path.exists():
                raise FileNotFoundError(path)
        if self.window < 0:
            raise ValueError("window must be >= 0")
        if self.topk < 1:
            raise ValueError("topk must be >= 1")
        if self.min_score < 0:
            raise ValueError("min_score must be >= 0")


@dataclass(frozen=True)
class LiftonMatchSummary:
    transcript_pairs: int
    gene_pairs: int
    transcript_pairs_path: Path
    gene_pairs_path: Path


@dataclass(frozen=True)
class MatchDiagnostics:
    query_transcripts: int
    queries_with_window_candidates: int
    scored_candidates: int
    queries_above_min_score: int
    best_scores: tuple[float, ...] = field(default_factory=tuple)
    no_candidate_examples: tuple[str, ...] = field(default_factory=tuple)
    below_threshold_examples: tuple[str, ...] = field(default_factory=tuple)

    @property
    def max_score(self) -> float:
        return max(self.best_scores) if self.best_scores else 0.0

    @property
    def median_score(self) -> float:
        return statistics.median(self.best_scores) if self.best_scores else 0.0

    @property
    def p90_score(self) -> float:
        if not self.best_scores:
            return 0.0
        ordered = sorted(self.best_scores)
        index = min(len(ordered) - 1, int(0.9 * (len(ordered) - 1)))
        return ordered[index]


def _transcripts_with_structure(annotation) -> int:
    return sum(1 for transcript in annotation.tx_index.values() if transcript.exons)


def _matching_diagnostics(
    lifton,
    target,
    window: int,
    min_score: float,
) -> MatchDiagnostics:
    index = lifton_id_mapper.GeneIntervalIndex(target)
    query_transcripts = 0
    queries_with_window_candidates = 0
    scored_candidates = 0
    queries_above_min_score = 0
    best_scores: list[float] = []
    no_candidate_examples: list[str] = []
    below_threshold_examples: list[str] = []

    for query_id, query in lifton.tx_index.items():
        if not query.exons:
            continue
        query_transcripts += 1
        candidates = lifton_id_mapper.candidate_transcripts_for(
            query,
            target,
            index,
            window,
        )
        if not candidates:
            if len(no_candidate_examples) < 3:
                no_candidate_examples.append(query_id)
            continue
        queries_with_window_candidates += 1

        best_score = 0.0
        for candidate in candidates:
            if candidate.contig != query.contig or candidate.strand != query.strand:
                continue
            score, _details = lifton_id_mapper.score_pair(query, candidate)
            scored_candidates += 1
            best_score = max(best_score, score)

        if best_score >= min_score:
            queries_above_min_score += 1
        elif len(below_threshold_examples) < 3:
            below_threshold_examples.append(f"{query_id}:{best_score:.3f}")
        best_scores.append(best_score)

    return MatchDiagnostics(
        query_transcripts=query_transcripts,
        queries_with_window_candidates=queries_with_window_candidates,
        scored_candidates=scored_candidates,
        queries_above_min_score=queries_above_min_score,
        best_scores=tuple(best_scores),
        no_candidate_examples=tuple(no_candidate_examples),
        below_threshold_examples=tuple(below_threshold_examples),
    )


def run_lifton_matching(config: LiftonMatchConfig) -> LiftonMatchSummary:
    config.validate()
    lifton = lifton_id_mapper.load_gff3_as_annotation(str(config.lifton_gff), "LiftOn")
    target = lifton_id_mapper.load_gff3_as_annotation(str(config.target_gff), "Target")
    sys.stderr.write(
        "Structural matching inputs: "
        f"lifton_genes={len(lifton.genes)}, "
        f"lifton_transcripts={len(lifton.tx_index)}, "
        f"lifton_transcripts_with_structure={_transcripts_with_structure(lifton)}, "
        f"target_genes={len(target.genes)}, "
        f"target_transcripts={len(target.tx_index)}, "
        f"target_transcripts_with_structure={_transcripts_with_structure(target)}\n"
    )
    diagnostics = _matching_diagnostics(
        lifton,
        target,
        config.window,
        config.min_score,
    )
    sys.stderr.write(
        "Structural matching candidate diagnostics: "
        f"query_transcripts={diagnostics.query_transcripts}, "
        f"queries_with_window_candidates={diagnostics.queries_with_window_candidates}, "
        f"scored_candidates={diagnostics.scored_candidates}, "
        f"queries_above_min_score={diagnostics.queries_above_min_score}, "
        f"best_score_max={diagnostics.max_score:.3f}, "
        f"best_score_median={diagnostics.median_score:.3f}, "
        f"best_score_p90={diagnostics.p90_score:.3f}, "
        f"no_candidate_examples={','.join(diagnostics.no_candidate_examples) or 'none'}, "
        "below_threshold_examples="
        f"{','.join(diagnostics.below_threshold_examples) or 'none'}\n"
    )
    pairs = lifton_id_mapper.compute_pairs(
        lifton,
        target,
        window=config.window,
        min_candidate_score=config.min_score,
        topk=config.topk,
    )
    lifton_id_mapper.write_transcript_pairs(
        pairs,
        lifton,
        target,
        str(config.transcript_pairs_path),
        conf_thresh=config.confident,
        good_thresh=config.good,
    )
    gene_pairs = lifton_id_mapper.aggregate_gene_pairs(
        pairs,
        lifton,
        target,
        min_gene_fraction=config.gene_fraction,
    )
    lifton_id_mapper.write_gene_pairs(gene_pairs, str(config.gene_pairs_path))
    if not pairs:
        sys.stderr.write(
            "Warning: structural matching produced no transcript pairs. "
            "Check contig names, strands, transcript structure rows, and match thresholds.\n"
        )
    sys.stderr.write(
        "Structural matching outputs: "
        f"transcript_pairs={len(pairs)}, gene_pairs={len(gene_pairs)}\n"
    )
    return LiftonMatchSummary(
        transcript_pairs=len(pairs),
        gene_pairs=len(gene_pairs),
        transcript_pairs_path=config.transcript_pairs_path,
        gene_pairs_path=config.gene_pairs_path,
    )

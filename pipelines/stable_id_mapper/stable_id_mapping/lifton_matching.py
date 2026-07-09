"""Run structural matching between LiftOn output and target annotation."""

from __future__ import annotations

import statistics
import sys
from bisect import bisect_left
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import lifton_id_mapper
from .rules import default_score_weights


@dataclass(frozen=True)
class LiftonMatchConfig:
    lifton_gff: Path
    target_gff: Path
    out_prefix: Path
    window: int = 100000
    topk: int = 5
    min_score: float = 0.30
    good: float = 0.45
    confident: float = 0.60
    gene_fraction: float = 0.60
    score_weights: dict[str, float] = field(default_factory=default_score_weights)

    @property
    def transcript_pairs_path(self) -> Path:
        return Path(f"{self.out_prefix}.transcript_pairs.tsv")

    @property
    def gene_pairs_path(self) -> Path:
        return Path(f"{self.out_prefix}.gene_pairs.tsv")

    @property
    def gene_locus_comparison_path(self) -> Path:
        return Path(f"{self.out_prefix}.gene_locus_comparison.tsv")

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
    gene_locus_comparison_path: Optional[Path] = None


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


def _gene_feature_span(gene) -> tuple[int, int]:
    start = getattr(gene, "feature_start", 0)
    end = getattr(gene, "feature_end", 0)
    if start and end:
        return min(start, end), max(start, end)
    return gene.span()


def _span_len(start: int, end: int) -> int:
    if start <= 0 or end < start:
        return 0
    return end - start + 1


def _span_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def _gene_locus_metrics(query_gene, target_gene) -> dict[str, float | int]:
    q_start, q_end = _gene_feature_span(query_gene)
    t_start, t_end = _gene_feature_span(target_gene)
    q_len = _span_len(q_start, q_end)
    t_len = _span_len(t_start, t_end)
    overlap_bp = _span_overlap(q_start, q_end, t_start, t_end)
    if q_len == 0 or t_len == 0 or overlap_bp == 0:
        return {
            "overlap_bp": overlap_bp,
            "old_gene_coverage": 0.0,
            "target_gene_coverage": 0.0,
            "reciprocal_overlap": 0.0,
            "span_containment": 0.0,
            "length_ratio": 0.0,
            "center_distance": abs(((q_start + q_end) / 2) - ((t_start + t_end) / 2)),
            "locus_score": 0.0,
        }
    old_coverage = overlap_bp / q_len
    target_coverage = overlap_bp / t_len
    reciprocal_overlap = overlap_bp / max(q_len, t_len)
    span_containment = overlap_bp / min(q_len, t_len)
    length_ratio = min(q_len, t_len) / max(q_len, t_len)
    center_distance = abs(((q_start + q_end) / 2) - ((t_start + t_end) / 2))
    return {
        "overlap_bp": overlap_bp,
        "old_gene_coverage": old_coverage,
        "target_gene_coverage": target_coverage,
        "reciprocal_overlap": reciprocal_overlap,
        "span_containment": span_containment,
        "length_ratio": length_ratio,
        "center_distance": center_distance,
        "locus_score": reciprocal_overlap,
    }


def _gene_locus_index(annotation) -> dict[str, dict[str, list[tuple[int, int, str]]]]:
    index: dict[str, dict[str, list[tuple[int, int, str]]]] = {}
    for gene_id, gene in annotation.genes.items():
        start, end = _gene_feature_span(gene)
        if _span_len(start, end) == 0:
            continue
        by_strand = index.setdefault(gene.contig, {})
        by_strand.setdefault(gene.strand, []).append((start, end, gene_id))
    for by_strand in index.values():
        for rows in by_strand.values():
            rows.sort(key=lambda item: item[0])
    return index


def _query_gene_locus_index(
    index: dict[str, dict[str, list[tuple[int, int, str]]]],
    contig: str,
    strand: str,
    start: int,
    end: int,
    window: int,
) -> list[str]:
    rows = index.get(contig, {}).get(strand, [])
    if not rows:
        return []
    starts = [row[0] for row in rows]
    left = max(0, bisect_left(starts, start - window) - 5)
    max_start = end + window
    out: list[str] = []
    for row_start, row_end, gene_id in rows[left:]:
        if row_start > max_start:
            break
        if row_end >= start - window:
            out.append(gene_id)
    return out


def _best_locus_gene(query_gene, target, index, window: int):
    q_start, q_end = _gene_feature_span(query_gene)
    if _span_len(q_start, q_end) == 0:
        return None, None
    candidates = _query_gene_locus_index(
        index,
        query_gene.contig,
        query_gene.strand,
        q_start,
        q_end,
        window,
    )
    best_gene = None
    best_metrics = None
    best_key = None
    for candidate_id in candidates:
        candidate = target.genes.get(candidate_id)
        if candidate is None:
            continue
        metrics = _gene_locus_metrics(query_gene, candidate)
        if metrics["overlap_bp"] == 0:
            continue
        key = (
            metrics["locus_score"],
            metrics["span_containment"],
            metrics["old_gene_coverage"],
            -metrics["center_distance"],
            candidate_id,
        )
        if best_key is None or key > best_key:
            best_gene = candidate
            best_metrics = metrics
            best_key = key
    return best_gene, best_metrics


def _structural_gene_lookup(gene_pairs) -> dict[str, tuple[str, float, float, int]]:
    return {
        lifton_gene: (target_gene, weighted_score, fraction, n_transcripts)
        for lifton_gene, target_gene, weighted_score, fraction, n_transcripts in gene_pairs
    }


def _best_transcript_structure_by_gene(pairs, lifton, target) -> dict[str, tuple[str, str, float]]:
    best: dict[str, tuple[str, str, float]] = {}
    for pair in pairs:
        query_tx = lifton.tx_index.get(pair.q_id)
        target_tx = target.tx_index.get(pair.r_id)
        if query_tx is None or target_tx is None:
            continue
        previous = best.get(query_tx.gene_id)
        if previous is None or pair.score > previous[2]:
            best[query_tx.gene_id] = (target_tx.gene_id, pair.r_id, pair.score)
    return best


def write_gene_locus_comparison(
    lifton,
    target,
    pairs,
    gene_pairs,
    path: Path,
    window: int,
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    index = _gene_locus_index(target)
    structural_genes = _structural_gene_lookup(gene_pairs)
    best_transcripts = _best_transcript_structure_by_gene(pairs, lifton, target)
    rows_written = 0
    header = [
        "lifton_gene",
        "target_gene_by_locus",
        "locus_score",
        "overlap_bp",
        "old_gene_coverage",
        "target_gene_coverage",
        "reciprocal_overlap",
        "span_containment",
        "length_ratio",
        "center_distance",
        "lifton_contig",
        "lifton_start",
        "lifton_end",
        "lifton_strand",
        "target_contig",
        "target_start",
        "target_end",
        "target_strand",
        "structure_accepted_target_gene",
        "structure_weighted_score",
        "structure_fraction_of_total",
        "structure_n_transcripts",
        "best_tx_structure_target_gene",
        "best_tx_structure_target_tx",
        "best_tx_structure_score",
        "locus_vs_structure",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for lifton_gene_id in sorted(lifton.genes):
            lifton_gene = lifton.genes[lifton_gene_id]
            q_start, q_end = _gene_feature_span(lifton_gene)
            locus_gene, metrics = _best_locus_gene(lifton_gene, target, index, window)
            target_gene_id = ""
            target_contig = ""
            target_start = ""
            target_end = ""
            target_strand = ""
            if locus_gene is not None and metrics is not None:
                t_start, t_end = _gene_feature_span(locus_gene)
                target_gene_id = locus_gene.id
                target_contig = locus_gene.contig
                target_start = str(t_start)
                target_end = str(t_end)
                target_strand = locus_gene.strand
            else:
                metrics = {
                    "locus_score": 0.0,
                    "overlap_bp": 0,
                    "old_gene_coverage": 0.0,
                    "target_gene_coverage": 0.0,
                    "reciprocal_overlap": 0.0,
                    "span_containment": 0.0,
                    "length_ratio": 0.0,
                    "center_distance": 0.0,
                }

            structural = structural_genes.get(lifton_gene_id)
            if structural is None:
                structural_target_gene = ""
                structural_weighted_score = ""
                structural_fraction = ""
                structural_n_transcripts = ""
            else:
                (
                    structural_target_gene,
                    weighted_score,
                    fraction,
                    n_transcripts,
                ) = structural
                structural_weighted_score = f"{weighted_score:.6f}"
                structural_fraction = f"{fraction:.6f}"
                structural_n_transcripts = str(n_transcripts)

            best_tx = best_transcripts.get(lifton_gene_id)
            if best_tx is None:
                best_tx_target_gene = ""
                best_tx_target_tx = ""
                best_tx_score = ""
            else:
                best_tx_target_gene, best_tx_target_tx, score = best_tx
                best_tx_score = f"{score:.6f}"

            if not target_gene_id:
                comparison = "no_locus_candidate"
            elif not structural_target_gene:
                comparison = "no_accepted_structure"
            elif target_gene_id == structural_target_gene:
                comparison = "same"
            else:
                comparison = "different"

            row = [
                lifton_gene_id,
                target_gene_id,
                f"{metrics['locus_score']:.6f}",
                str(metrics["overlap_bp"]),
                f"{metrics['old_gene_coverage']:.6f}",
                f"{metrics['target_gene_coverage']:.6f}",
                f"{metrics['reciprocal_overlap']:.6f}",
                f"{metrics['span_containment']:.6f}",
                f"{metrics['length_ratio']:.6f}",
                f"{metrics['center_distance']:.1f}",
                lifton_gene.contig,
                str(q_start),
                str(q_end),
                lifton_gene.strand,
                target_contig,
                target_start,
                target_end,
                target_strand,
                structural_target_gene,
                structural_weighted_score,
                structural_fraction,
                structural_n_transcripts,
                best_tx_target_gene,
                best_tx_target_tx,
                best_tx_score,
                comparison,
            ]
            handle.write("\t".join(row) + "\n")
            rows_written += 1
    return rows_written


def _matching_diagnostics(
    lifton,
    target,
    window: int,
    min_score: float,
    score_weights: dict[str, float],
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
            score, _details = lifton_id_mapper.score_pair(
                query,
                candidate,
                score_weights=score_weights,
            )
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
        config.score_weights,
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
        score_weights=config.score_weights,
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
    locus_comparison_rows = write_gene_locus_comparison(
        lifton,
        target,
        pairs,
        gene_pairs,
        config.gene_locus_comparison_path,
        config.window,
    )
    if not pairs:
        sys.stderr.write(
            "Warning: structural matching produced no transcript pairs. "
            "Check contig names, strands, transcript structure rows, and match thresholds.\n"
        )
    sys.stderr.write(
        "Structural matching outputs: "
        f"transcript_pairs={len(pairs)}, gene_pairs={len(gene_pairs)}, "
        f"gene_locus_comparison_rows={locus_comparison_rows}\n"
    )
    return LiftonMatchSummary(
        transcript_pairs=len(pairs),
        gene_pairs=len(gene_pairs),
        transcript_pairs_path=config.transcript_pairs_path,
        gene_pairs_path=config.gene_pairs_path,
        gene_locus_comparison_path=config.gene_locus_comparison_path,
    )

"""Run structural matching between LiftOn output and target annotation."""

from __future__ import annotations

import sys
from dataclasses import dataclass
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


def _transcripts_with_structure(annotation) -> int:
    return sum(1 for transcript in annotation.tx_index.values() if transcript.exons)


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

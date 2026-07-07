"""Single-species stable-ID mapping workflow."""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .config import StableIdEventConfig, default_backup_prefix
from .ids import StableIdRange
from .lifton import DEFAULT_LIFTON_FEATURE_TYPES, LiftonRunConfig, run_lifton
from .lifton_matching import LiftonMatchConfig, LiftonMatchSummary, run_lifton_matching
from .models import Decision
from .pipeline import run_pipeline
from .reports import write_missing_gene_report
from .scoring import load_lifton_score_evidence


@dataclass(frozen=True)
class SingleSpeciesRunConfig:
    ref_fasta: Path
    ref_gff: Path
    target_fasta: Path
    target_gff: Path
    db_name: str
    mapping_session_id: int
    gene_range: StableIdRange
    transcript_range: StableIdRange
    translation_range: StableIdRange
    output_dir: Path
    existing_lifton_gff: Optional[Path] = None
    existing_transcript_pairs: Optional[Path] = None
    existing_gene_pairs: Optional[Path] = None
    lifton_threads: int = 8
    lifton_executable: str = "lifton"
    lifton_feature_types: tuple[str, ...] = DEFAULT_LIFTON_FEATURE_TYPES
    lifton_feature_types_file: Optional[Path] = None
    lifton_extra_args: tuple[str, ...] = field(default_factory=tuple)
    include_translations: bool = True
    dry_run_sql: bool = True
    replace_events_for_session: bool = False
    backup_prefix: str = field(default_factory=default_backup_prefix)
    batch_size: int = 500
    min_overlap: float = 0.10
    match_window: int = 100000
    match_topk: int = 5
    match_min_score: float = 0.60
    match_good: float = 0.75
    match_confident: float = 0.85
    match_gene_fraction: float = 0.60

    @property
    def lifton_output_gff(self) -> Path:
        if self.existing_lifton_gff is not None:
            return self.existing_lifton_gff
        return self.output_dir / "lifton" / "projected_ref_on_target.gff3"

    @property
    def missing_report(self) -> Path:
        return self.output_dir / "reports" / "missing_genes.txt"

    @property
    def match_out_prefix(self) -> Path:
        return self.output_dir / "matching" / "lifton"

    @property
    def output_sql(self) -> Path:
        suffix = ".dry_run.sql" if self.dry_run_sql else ".sql"
        return self.output_dir / "sql" / f"stable_id_updates{suffix}"

    @property
    def output_tsv(self) -> Path:
        return self.output_dir / "stable_id_decisions.tsv"

    @property
    def output_score_evidence_tsv(self) -> Path:
        return self.output_dir / "score_evidence.tsv"

    def validate_inputs(self) -> None:
        for path in (self.ref_fasta, self.ref_gff, self.target_fasta, self.target_gff):
            if not path.exists():
                raise FileNotFoundError(path)
        for path in (
            self.existing_lifton_gff,
            self.existing_transcript_pairs,
            self.existing_gene_pairs,
        ):
            if path is not None and not path.exists():
                raise FileNotFoundError(path)
        if (self.existing_transcript_pairs is None) != (self.existing_gene_pairs is None):
            raise ValueError(
                "existing_transcript_pairs and existing_gene_pairs must be provided together"
            )
        if not self.db_name:
            raise ValueError("db_name is required")
        if self.mapping_session_id < 1:
            raise ValueError("mapping_session_id must be >= 1")


@dataclass(frozen=True)
class SingleSpeciesRunResult:
    lifton_gff: Path
    missing_report: Path
    match_summary: LiftonMatchSummary
    output_sql: Path
    output_tsv: Path
    output_score_evidence_tsv: Path
    decisions: list[Decision]


def run_single_species_pipeline(
    config: SingleSpeciesRunConfig,
) -> SingleSpeciesRunResult:
    config.validate_inputs()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    if config.existing_lifton_gff is None:
        lifton_config = LiftonRunConfig(
            ref_gff=config.ref_gff,
            ref_fasta=config.ref_fasta,
            target_fasta=config.target_fasta,
            output_gff=config.lifton_output_gff,
            threads=config.lifton_threads,
            executable=config.lifton_executable,
            feature_types=config.lifton_feature_types,
            feature_types_file=config.lifton_feature_types_file,
            extra_args=config.lifton_extra_args,
        )
        run_lifton(lifton_config)

    missing_gene_report = config.missing_report
    write_missing_gene_report(
        config.ref_gff,
        config.lifton_output_gff,
        missing_gene_report,
    )

    if config.existing_transcript_pairs is None or config.existing_gene_pairs is None:
        match_summary = run_lifton_matching(
            LiftonMatchConfig(
                lifton_gff=config.lifton_output_gff,
                target_gff=config.target_gff,
                out_prefix=config.match_out_prefix,
                window=config.match_window,
                topk=config.match_topk,
                min_score=config.match_min_score,
                good=config.match_good,
                confident=config.match_confident,
                gene_fraction=config.match_gene_fraction,
            )
        )
    else:
        transcript_pair_count = _count_data_rows(config.existing_transcript_pairs)
        gene_pair_count = _count_data_rows(config.existing_gene_pairs)
        if transcript_pair_count == 0 and gene_pair_count == 0:
            sys.stderr.write(
                "Warning: reusing empty structural pair tables; "
                "run with only --existing-lifton-gff to regenerate matching.\n"
            )
        match_summary = LiftonMatchSummary(
            transcript_pairs=transcript_pair_count,
            gene_pairs=gene_pair_count,
            transcript_pairs_path=config.existing_transcript_pairs,
            gene_pairs_path=config.existing_gene_pairs,
        )
    score_evidence = load_lifton_score_evidence(
        match_summary.transcript_pairs_path,
        match_summary.gene_pairs_path,
    )
    decisions = run_pipeline(
        StableIdEventConfig(
            ref_gff=config.ref_gff,
            target_gff=config.target_gff,
            mapped_gff=config.lifton_output_gff,
            report=missing_gene_report,
            mapping_session_id=config.mapping_session_id,
            gene_range=config.gene_range,
            transcript_range=config.transcript_range,
            translation_range=config.translation_range,
            output_sql=config.output_sql,
            output_tsv=config.output_tsv,
            db_name=config.db_name,
            include_translations=config.include_translations,
            dry_run=config.dry_run_sql,
            backup_prefix=config.backup_prefix,
            batch_size=config.batch_size,
            replace_events_for_session=config.replace_events_for_session,
            min_overlap=config.min_overlap,
            score_evidence=score_evidence,
            output_score_evidence_tsv=config.output_score_evidence_tsv,
        )
    )

    return SingleSpeciesRunResult(
        lifton_gff=config.lifton_output_gff,
        missing_report=missing_gene_report,
        match_summary=match_summary,
        output_sql=config.output_sql,
        output_tsv=config.output_tsv,
        output_score_evidence_tsv=config.output_score_evidence_tsv,
        decisions=decisions,
    )


def _count_data_rows(path: Path) -> int:
    with path.open() as handle:
        return max(0, sum(1 for _line in handle) - 1)

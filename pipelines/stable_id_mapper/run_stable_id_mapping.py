#!/usr/bin/env python3
"""Run one species stable-ID mapping workflow."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

from run_lifton_projection import parse_feature_types
from stable_id_mapping.config import default_backup_prefix
from stable_id_mapping.ids import parse_id_range
from stable_id_mapping.lifton import (
    LiftonRunConfig,
    build_lifton_command,
    format_command,
    prepare_lifton_feature_types_file,
)
from stable_id_mapping.rules import DEFAULT_RULES_PATH, load_mapping_rules
from stable_id_mapping.workflow import SingleSpeciesRunConfig, run_single_species_pipeline


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a single-species stable-ID mapping workflow."
    )
    parser.add_argument("--ref-fasta", type=Path, required=True)
    parser.add_argument("--ref-gff", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--target-gff", type=Path, required=True)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--mapping-session-id", type=int, required=True)
    parser.add_argument("--gene-range", type=parse_id_range, required=True)
    parser.add_argument("--transcript-range", type=parse_id_range, required=True)
    parser.add_argument("--translation-range", type=parse_id_range, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--rules-config",
        type=Path,
        default=DEFAULT_RULES_PATH,
        help="JSON file with stable-ID mapping rules and thresholds",
    )

    parser.add_argument(
        "--existing-lifton-gff",
        type=Path,
        help="Reuse an existing LiftOn-projected GFF3 instead of running LiftOn",
    )
    parser.add_argument(
        "--existing-transcript-pairs",
        type=Path,
        help="Reuse an existing lifton.transcript_pairs.tsv file",
    )
    parser.add_argument(
        "--existing-gene-pairs",
        type=Path,
        help="Reuse an existing lifton.gene_pairs.tsv file",
    )

    parser.add_argument("--lifton-threads", type=int, default=8)
    parser.add_argument("--lifton-executable", default="lifton")
    parser.add_argument(
        "--lifton-feature-types",
        type=parse_feature_types,
        default=None,
        help="Comma-separated parent feature types written to the lifton -f file",
    )
    parser.add_argument(
        "--lifton-feature-types-file",
        type=Path,
        help="Feature-type file passed to lifton -f",
    )
    parser.add_argument(
        "--extra-lifton-arg",
        action="append",
        default=[],
        help="Additional argument passed to lifton before FASTA positional arguments",
    )

    parser.add_argument("--match-window", type=int)
    parser.add_argument("--match-topk", type=int)
    parser.add_argument("--match-min-score", type=float)
    parser.add_argument("--match-good", type=float)
    parser.add_argument("--match-confident", type=float)
    parser.add_argument("--match-gene-fraction", type=float)

    parser.add_argument("--no-translations", action="store_true")
    parser.add_argument(
        "--write-executable-sql",
        action="store_true",
        help="Write COMMIT SQL instead of dry-run ROLLBACK SQL",
    )
    parser.add_argument("--replace-events-for-session", action="store_true")
    parser.add_argument("--backup-prefix", default=default_backup_prefix())
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--min-overlap", type=float)
    parser.add_argument(
        "--dry-run-lifton-command",
        action="store_true",
        help="Print the LiftOn command for this run and exit",
    )
    return parser.parse_args(argv)


def config_from_args(args: argparse.Namespace) -> SingleSpeciesRunConfig:
    rules = load_mapping_rules(args.rules_config)
    structural = rules.structural_matching
    kwargs = {}
    if args.lifton_feature_types is not None:
        kwargs["lifton_feature_types"] = args.lifton_feature_types
    return SingleSpeciesRunConfig(
        ref_fasta=args.ref_fasta,
        ref_gff=args.ref_gff,
        target_fasta=args.target_fasta,
        target_gff=args.target_gff,
        db_name=args.db_name,
        mapping_session_id=args.mapping_session_id,
        gene_range=args.gene_range,
        transcript_range=args.transcript_range,
        translation_range=args.translation_range,
        output_dir=args.output_dir,
        existing_lifton_gff=args.existing_lifton_gff,
        existing_transcript_pairs=args.existing_transcript_pairs,
        existing_gene_pairs=args.existing_gene_pairs,
        lifton_threads=args.lifton_threads,
        lifton_executable=args.lifton_executable,
        lifton_feature_types_file=args.lifton_feature_types_file,
        lifton_extra_args=tuple(args.extra_lifton_arg),
        include_translations=not args.no_translations,
        dry_run_sql=not args.write_executable_sql,
        replace_events_for_session=args.replace_events_for_session,
        backup_prefix=args.backup_prefix,
        batch_size=args.batch_size,
        min_overlap=(
            args.min_overlap
            if args.min_overlap is not None
            else rules.coordinate_overlap.min_overlap
        ),
        match_window=(
            args.match_window if args.match_window is not None else structural.window
        ),
        match_topk=args.match_topk if args.match_topk is not None else structural.topk,
        match_min_score=(
            args.match_min_score
            if args.match_min_score is not None
            else structural.min_score
        ),
        match_good=args.match_good if args.match_good is not None else structural.good_score,
        match_confident=(
            args.match_confident
            if args.match_confident is not None
            else structural.confident_score
        ),
        match_gene_fraction=(
            args.match_gene_fraction
            if args.match_gene_fraction is not None
            else structural.gene_fraction
        ),
        match_score_weights=dict(structural.score_weights),
        rules_config=args.rules_config,
        **kwargs,
    )


def lifton_command_for_config(
    config: SingleSpeciesRunConfig,
    prepare_feature_types: bool = False,
) -> list[str]:
    if config.existing_lifton_gff is not None:
        return []
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
    if prepare_feature_types:
        prepare_lifton_feature_types_file(lifton_config)
    return build_lifton_command(lifton_config)


def main() -> None:
    args = parse_args()
    config = config_from_args(args)
    if config.output_dir.exists():
        sys.stderr.write(f"Using existing output directory {config.output_dir}\n")
    if args.dry_run_lifton_command:
        command = lifton_command_for_config(config, prepare_feature_types=True)
        if command:
            print(format_command(command))
        else:
            print(f"Reusing existing LiftOn GFF3: {config.existing_lifton_gff}")
        return

    result = run_single_species_pipeline(config)
    sys.stderr.write("Stable-ID mapping run complete.\n")
    sys.stderr.write(f"LiftOn GFF3: {result.lifton_gff}\n")
    sys.stderr.write(f"Transcript scores: {result.match_summary.transcript_pairs_path}\n")
    sys.stderr.write(f"Gene scores: {result.match_summary.gene_pairs_path}\n")
    sys.stderr.write(f"Score evidence: {result.output_score_evidence_tsv}\n")
    sys.stderr.write(f"Stable-ID TSV: {result.output_tsv}\n")
    sys.stderr.write(f"SQL: {result.output_sql}\n")


if __name__ == "__main__":
    main()

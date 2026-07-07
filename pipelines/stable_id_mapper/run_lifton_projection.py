#!/usr/bin/env python3
"""Run LiftOn projection for stable-ID mapping."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

from stable_id_mapping.lifton import (
    DEFAULT_LIFTON_FEATURE_TYPES,
    LiftonRunConfig,
    build_lifton_command,
    format_command,
    prepare_lifton_feature_types_file,
    run_lifton,
)


def parse_feature_types(value: str) -> tuple[str, ...]:
    feature_types = tuple(item.strip() for item in value.split(",") if item.strip())
    if not feature_types:
        raise argparse.ArgumentTypeError("At least one feature type is required")
    return feature_types


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run LiftOn to project a reference GFF3 onto a target FASTA."
    )
    parser.add_argument("--ref-gff", type=Path, required=True)
    parser.add_argument("--ref-fasta", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--output-gff", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--lifton-executable", default="lifton")
    parser.add_argument(
        "--feature-types",
        type=parse_feature_types,
        default=DEFAULT_LIFTON_FEATURE_TYPES,
        help="Comma-separated feature types passed to lifton -f",
    )
    parser.add_argument(
        "--feature-types-file",
        type=Path,
        help="Feature-type file passed to lifton -f; overrides --feature-types",
    )
    parser.add_argument(
        "--extra-lifton-arg",
        action="append",
        default=[],
        help="Additional argument passed to lifton before FASTA positional arguments",
    )
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument(
        "--dry-run-command",
        action="store_true",
        help="Print the lifton command without executing it",
    )
    return parser.parse_args(argv)


def config_from_args(args: argparse.Namespace) -> LiftonRunConfig:
    return LiftonRunConfig(
        ref_gff=args.ref_gff,
        ref_fasta=args.ref_fasta,
        target_fasta=args.target_fasta,
        output_gff=args.output_gff,
        threads=args.threads,
        executable=args.lifton_executable,
        feature_types=args.feature_types,
        feature_types_file=args.feature_types_file,
        extra_args=tuple(args.extra_lifton_arg),
        work_dir=args.work_dir,
    )


def main() -> None:
    args = parse_args()
    config = config_from_args(args)
    if config.output_gff.exists():
        sys.stderr.write(f"Warning: overwriting existing output {config.output_gff}\n")
    if args.dry_run_command:
        prepare_lifton_feature_types_file(config)
        command = build_lifton_command(config)
        print(format_command(command))
        return
    command = build_lifton_command(config)
    sys.stderr.write(format_command(command) + "\n")
    run_lifton(config)


if __name__ == "__main__":
    main()

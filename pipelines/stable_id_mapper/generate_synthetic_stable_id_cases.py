#!/usr/bin/env python3
"""Generate controlled synthetic stable-ID mapping examples."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from stable_id_mapping.synthetic_cases import (
    available_case_names,
    format_shell_command,
    stable_id_mapping_command,
    write_synthetic_cases,
)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write synthetic stable-ID mapping fixtures with known outcomes."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("pipelines/stable_id_mapper/out/synthetic_cases"),
        help="Directory where case folders will be written",
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=available_case_names(),
        help="Case to generate. May be provided more than once. Defaults to all cases.",
    )
    parser.add_argument(
        "--rules-config",
        default="pipelines/stable_id_mapper/stable_id_mapping_rules_locus_calibration.json",
        help="Rules config path to include in printed commands",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    cases = write_synthetic_cases(args.output_dir, args.case)
    for case in cases:
        print(f"\n{case.name}: {case.description}")
        print(f"  files: {case.case_dir}")
        print(f"  expected: {case.expected_decisions}")
        command = stable_id_mapping_command(
            case,
            rules_config=args.rules_config,
        )
        print("  command:")
        print("    " + format_shell_command(command))


if __name__ == "__main__":
    main()

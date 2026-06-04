"""annotation-qc CLI entry point."""

import argparse

from ensembl.genes.annotation_qc.runners import (
    assess_translation_validity,
    genebuild_stats,
)

RUNNERS = [
    assess_translation_validity,
    genebuild_stats,
]


def main() -> None:
    """Parse arguments and dispatch to the selected runner."""
    parser = argparse.ArgumentParser(
        prog="annotation-qc",
        description="Run quality checks on genome annotations.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for runner in RUNNERS:
        runner.register(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

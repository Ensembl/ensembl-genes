"""
annotation-qc CLI entry point.

Each runner registers its own subcommand via a register(subparsers) function.
To add a new runner: import it here and append it to RUNNERS.

Usage:
        python cli.py translation-validity \
            --annotation_file annotations.gff3 --genome_fasta genome.fa
"""

import argparse

from ensembl.genes.annotation_qc.parsers import interpro, parse_agat
from ensembl.genes.annotation_qc.runners import (
    assess_feature_events,
    assess_translation_validity,
)

RUNNERS = [
    assess_translation_validity,
    assess_feature_events,
]

PARSERS = [
    parse_agat,
    interpro,
]


def main():
    """Parse arguments and dispatch to the appropriate runner subcommand."""
    parser = argparse.ArgumentParser(
        prog="annotation-qc",
        description="Run quality checks on genome annotations.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for runner in RUNNERS:
        runner.register(subparsers)

    for parser_module in PARSERS:
        parser_module.register(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

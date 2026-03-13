"""
annotation-qc CLI entry point.

Each runner registers its own subcommand via a register(subparsers) function.
To add a new runner: import it here and append it to RUNNERS.

Usage:
	python cli.py translation-validity --annotation_file annotations.gff3 --genome_fasta genome.fa
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from runners import assess_translation_validity, genebuild_stats

RUNNERS = [
	assess_translation_validity,
	genebuild_stats,
]


def main():
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

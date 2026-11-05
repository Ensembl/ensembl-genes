#!/usr/bin/env python3

"""
Summarise InterProScan protein hits.

Reads an InterProScan TSV file and reports the percentage of query proteins
that have at least one InterProScan hit.

Usage:
    python interpro.py \
        --input_tsv interproscan.tsv \
        --query_protein proteins.fa \
        --output results/
"""

import argparse
import csv
from pathlib import Path


def count_fasta(fasta: Path) -> int:
    """Count protein sequences in a FASTA file."""
    with fasta.open() as fh:
        return sum(1 for line in fh if line.startswith(">"))


def count_hits(tsv: Path) -> int:
    """Count unique query proteins with at least one InterProScan hit."""
    with tsv.open() as fh:
        return len(
            {
                line.split("\t", 1)[0]
                for line in fh
                if line.strip() and not line.startswith("#")
            }
        )


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add InterProScan parser arguments to an argparse parser."""
    parser.add_argument(
        "--input_tsv",
        type=Path,
        required=True,
        help="InterProScan TSV output.",
    )

    parser.add_argument(
        "--query_protein",
        type=Path,
        required=True,
        help="Protein FASTA used as InterProScan input.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory.",
    )

    parser.add_argument(
        "--sample_name",
        help="Sample/accession name. Defaults to FASTA filename.",
    )


def run(args: argparse.Namespace) -> None:
    """Run the InterProScan summary calculation."""
    if not args.input_tsv.is_file():
        raise FileNotFoundError(f"Input TSV not found: {args.input_tsv}")

    if not args.query_protein.is_file():
        raise FileNotFoundError(f"Protein FASTA not found: {args.query_protein}")

    args.output.mkdir(parents=True, exist_ok=True)

    sample = args.sample_name or args.query_protein.stem

    total = count_fasta(args.query_protein)
    hits = count_hits(args.input_tsv)

    if hits > total:
        raise ValueError(
            f"Found {hits} proteins with hits but FASTA contains only {total} proteins."
        )

    percent = 100 * hits / total if total else 0.0

    outfile = args.output / "interpro_hit_summary.tsv"

    with outfile.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        writer.writerow(
            [
                "accession",
                "protein_count",
                "protein_with_hit_count",
                "protein_with_hit_percent",
            ]
        )

        writer.writerow(
            [
                sample,
                total,
                hits,
                f"{percent:.2f}",
            ]
        )

    print(f"Wrote {outfile}")


def register(subparsers) -> None:
    """Register InterProScan parsing as an annotation QC subcommand."""
    parser = subparsers.add_parser(
        "parse-interpro",
        help="Summarise InterProScan protein hits.",
    )
    add_arguments(parser)
    parser.set_defaults(func=run)


def main() -> None:
    """Parse standalone InterProScan arguments and run the parser."""
    parser = argparse.ArgumentParser(description="Summarise InterProScan protein hits.")
    add_arguments(parser)
    run(parser.parse_args())


if __name__ == "__main__":
    main()

# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Create the final EAM gene confidence CSV from a scored evidence TSV."""
import argparse
from pathlib import Path

from ensembl.genes.annotation_qc.parsers.high_confidence import parse_tsv_frame
from ensembl.genes.annotation_qc.metrics.gene_confidence_summary import (
    summarize_gene_confidence,
)
from ensembl.genes.annotation_qc.reports.gene_confidence_summary import (
    write_gene_confidence_summary,
)


def add_arguments(parser):
    """Declare the scored evidence input and final CSV destination."""
    parser.add_argument(
        "--evidence",
        required=True,
        type=Path,
        help="Scored .evidence.tsv containing gene_id and confidence",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Destination metric,value CSV"
    )


def register(subparsers):
    """Register the final gene confidence summary command."""
    parser = subparsers.add_parser("high-confidence-gene-summary", help=__doc__)
    add_arguments(parser)
    parser.set_defaults(func=run_cli)


def run(args):
    """Read evidence, calculate unique-gene metrics and export the final CSV."""
    if args.evidence.resolve() == args.output.resolve():
        raise ValueError("Output must not overwrite the input evidence file")
    metrics = summarize_gene_confidence(parse_tsv_frame(args.evidence))
    write_gene_confidence_summary(metrics, args.output)
    return metrics


def run_cli(args):
    """Return a clear CLI failure for malformed input or file errors."""
    try:
        metrics = run(args)
    except (ValueError, OSError) as error:
        raise SystemExit(f"error: {error}") from error
    print(f"Summarized {metrics['eam.total_genes']} scored genes: {args.output}")
    return metrics


def main():
    """Support execution as a standalone installed module."""
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    run_cli(parser.parse_args())


if __name__ == "__main__":
    main()

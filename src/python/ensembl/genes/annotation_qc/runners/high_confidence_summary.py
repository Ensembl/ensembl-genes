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
"""Assemble evidence-summary and comparative QC reports."""
from __future__ import annotations
import argparse
from pathlib import Path
from ensembl.genes.annotation_qc.parsers.high_confidence import (
    parse_confidence_run,
    parse_confidence_models,
)
from ensembl.genes.annotation_qc.metrics.high_confidence_summary import (
    summarize_confidence,
    prepare_confidence_plots,
)
from ensembl.genes.annotation_qc.reports.high_confidence_summary import (
    write_summary,
    write_comparison,
    write_figures,
)


def add_arguments(parser):
    """Declare single-run and model-comparison report inputs."""
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--evidence", type=Path)
    source.add_argument(
        "--models",
        type=Path,
        help="TSV with model/evidence and optional manifest/introns paths",
    )
    parser.add_argument("--label", default="Prediction method")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--introns", type=Path)
    parser.add_argument("--species")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Also render PNG/SVG figures (requires annotation-qc-plots extra)",
    )


def register(subparsers):
    """Register summaries with the common CLI."""
    parser = subparsers.add_parser(
        "high-confidence-summary",
        help="Summarize evidence or compare prediction methods.",
    )
    add_arguments(parser)
    parser.set_defaults(func=run_cli)


def run_cli(args):
    """Parse sources, calculate metrics and pass them to report writers."""
    models = (
        parse_confidence_models(args.models).to_dict("records")
        if args.models
        else [
            {
                "model": args.label,
                "evidence": args.evidence,
                "manifest": args.manifest,
                "introns": args.introns,
            }
        ]
    )
    results = []
    for model in models:
        parsed = parse_confidence_run(
            model["evidence"], model.get("manifest"), model.get("introns")
        )
        result = summarize_confidence(
            parsed["evidence"],
            model["model"],
            parsed["manifest"],
            parsed["introns"],
            parsed["evidence_sha256"],
        )
        result["provenance"].update(parsed["source"])
        results.append(result)
    if args.models:
        write_comparison(results, args.output_dir)
    else:
        write_summary(results[0], args.output_dir)
    if args.plots:
        write_figures(
            prepare_confidence_plots(results), args.output_dir / "figures", args.species
        )
    for result in results:
        print(
            f"{result['model']}: {result['overview']['transcripts']} assessed genes/transcripts"
        )
        for warning in result["provenance"]["warnings"]:
            print("WARNING:", warning)
    return results


def main():
    """Support installed-module execution."""
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    run_cli(parser.parse_args())


if __name__ == "__main__":
    main()

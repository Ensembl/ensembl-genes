"""Run the current annotation feature and sequence-event QC metrics."""

from __future__ import annotations

import argparse
from pathlib import Path

from pyfaidx import Fasta

from ensembl.genes.annotation_qc.metrics.events import (
    compute_polyA_motif_metrics,
    compute_promoter_proximity_metrics,
    compute_splice_junction_metrics,
    compute_translation_metrics,
    derive_introns,
)
from ensembl.genes.annotation_qc.metrics.features import (
    compute_feature_metrics,
    summarize_feature_metrics,
)
from ensembl.genes.annotation_qc.parsers.annotation import parse_annotation


def _motifs(values: list[str] | None) -> dict[str, str]:
    """Parse repeated ``name=regex`` command-line motif arguments."""
    motifs = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"Motif must use NAME=REGEX syntax: {value}")
        name, pattern = value.split("=", 1)
        if not name or not pattern:
            raise ValueError(f"Motif must use NAME=REGEX syntax: {value}")
        motifs[name] = pattern
    return motifs


def run(
    annotation_path: str,
    genome_path: str,
    output_dir: str | Path,
    promoter_motifs: list[str] | None = None,
    promoter_upstream: int = 1000,
    polya_upstream: int = 50,
) -> dict[str, object]:
    """Calculate feature and genome-sequence event tables and write TSV files."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    annotation = parse_annotation(annotation_path)
    fasta = Fasta(genome_path)
    profile = {
        "promoter": {
            "motifs": _motifs(promoter_motifs),
            "upstream": promoter_upstream,
        },
        "polya": {
            "motifs": {"AATAAA": "AATAAA", "ATTAAA": "ATTAAA"},
            "upstream": polya_upstream,
        },
    }

    feature_metrics = compute_feature_metrics(annotation)
    feature_summary = summarize_feature_metrics(feature_metrics)
    introns = derive_introns(annotation)
    splice_metrics = compute_splice_junction_metrics(annotation, fasta)
    translation_metrics = compute_translation_metrics(
        annotation, fasta, profile=profile
    )
    promoter_metrics = compute_promoter_proximity_metrics(annotation, fasta, profile)
    polya_metrics = compute_polyA_motif_metrics(annotation, fasta, profile)

    outputs = {
        "feature_metrics": feature_metrics,
        "feature_summary": feature_summary,
        "intron_events": introns,
        "splice_events": splice_metrics,
        "translation_events": translation_metrics,
        "promoter_events": promoter_metrics,
        "polya_events": polya_metrics,
    }
    for name, table in outputs.items():
        table.to_csv(output / f"{name}.tsv", sep="\t", index=False)
    return outputs


def _add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--annotation", required=True, help="GFF3 or GTF annotation.")
    parser.add_argument("--genome", required=True, help="Indexed genome FASTA.")
    parser.add_argument("--outdir", default=".", help="Output directory.")
    parser.add_argument(
        "--promoter-motif",
        action="append",
        help="Promoter-like motif as NAME=REGEX; may be repeated.",
    )
    parser.add_argument(
        "--promoter-upstream",
        type=int,
        default=1000,
        help="Number of upstream bases to scan for promoter motifs.",
    )
    parser.add_argument(
        "--polya-upstream",
        type=int,
        default=50,
        help="Number of transcript-oriented upstream bases to scan for poly-A motifs.",
    )


def register(subparsers) -> None:
    """Register the feature/event QC subcommand."""
    parser = subparsers.add_parser(
        "feature-events",
        help="Calculate feature, splice, translation, promoter, and poly-A metrics.",
    )
    _add_arguments(parser)
    parser.set_defaults(func=_run)


def _run(args) -> dict[str, object]:
    return run(
        args.annotation,
        args.genome,
        args.outdir,
        args.promoter_motif,
        args.promoter_upstream,
        args.polya_upstream,
    )


def main() -> None:
    """Run the feature/event QC command standalone."""
    parser = argparse.ArgumentParser(description=__doc__)
    _add_arguments(parser)
    _run(parser.parse_args())


if __name__ == "__main__":
    main()

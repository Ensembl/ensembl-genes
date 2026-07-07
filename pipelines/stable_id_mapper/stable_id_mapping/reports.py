"""Report helpers for stable-ID mapping runs."""

from __future__ import annotations

from pathlib import Path

from .gff3 import parse_gff3


def missing_gene_ids_from_projection(
    ref_gff: str | Path,
    projected_gff: str | Path,
) -> set[str]:
    ref_features = parse_gff3(ref_gff)
    projected_features = parse_gff3(projected_gff)
    return set(ref_features["gene"]) - set(projected_features["gene"])


def write_missing_gene_report(
    ref_gff: str | Path,
    projected_gff: str | Path,
    output_report: str | Path,
) -> set[str]:
    ref_features = parse_gff3(ref_gff)
    projected_features = parse_gff3(projected_gff)
    missing_gene_ids = set(ref_features["gene"]) - set(projected_features["gene"])
    output_report = Path(output_report)
    output_report.parent.mkdir(parents=True, exist_ok=True)
    with output_report.open("w") as handle:
        handle.write("Stable ID mapper report\n")
        handle.write("=======================\n\n")
        handle.write(f"Genes in reference: {len(ref_features['gene'])}\n")
        handle.write(f"Projected by LiftOn: {len(projected_features['gene'])}\n")
        handle.write(f"Missing: {len(missing_gene_ids)}\n\n")
        handle.write("Missing gene IDs:\n")
        for stable_id in sorted(missing_gene_ids):
            handle.write(f"  gene:{stable_id}\n")
    return missing_gene_ids

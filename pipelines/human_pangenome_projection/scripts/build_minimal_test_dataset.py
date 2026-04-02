#!/usr/bin/env python3
"""Build a minimal local test dataset around a target gene from cluster test data.

Default data sources:
- cluster_test_data/Homo_sapiens.GRCh38.115.chr.gff3
- cluster_test_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
- cluster_test_data/CHM13.unmasked.fa

Given a gene symbol or Ensembl gene ID, this script selects +/- N neighboring genes
in GFF3 file order (default N=5), writes their full gene blocks to test.gff3,
and writes chromosome FASTA sequences for those genes to reference.fa and target.fa.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pysam

GENE_FEATURE_TYPES = {"gene", "pseudogene", "ncRNA_gene"}


@dataclass
class GeneMeta:
    """Metadata for one top-level gene feature in GFF3 order."""

    index: int
    seqid: str
    gene_id: Optional[str]
    gene_id_raw: Optional[str]
    gene_name: Optional[str]


def parse_attributes(attr_field: str) -> Dict[str, str]:
    """Parse a GFF3 attributes column into a dictionary."""
    attrs: Dict[str, str] = {}
    for token in attr_field.strip().split(";"):
        if not token or "=" not in token:
            continue
        key, value = token.split("=", 1)
        attrs[key] = value
    return attrs


def normalize_gene_id(value: Optional[str]) -> Optional[str]:
    """Normalize IDs like 'gene:ENSG...' -> 'ENSG...' for matching."""
    if not value:
        return None
    if value.startswith("gene:"):
        return value.split(":", 1)[1]
    return value


def detect_gene_line(parts: Sequence[str]) -> bool:
    """Return True when a parsed GFF3 row is a top-level gene line."""
    return len(parts) >= 9 and parts[2] in GENE_FEATURE_TYPES


def scan_gene_metadata(gff_path: Path) -> Tuple[List[GeneMeta], str]:
    """Scan GFF3 once to collect top-level gene metadata in file order."""
    genes: List[GeneMeta] = []
    gff_version = "##gff-version 3"

    with gff_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##gff-version"):
                gff_version = line
                continue
            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if not detect_gene_line(parts):
                continue

            attrs = parse_attributes(parts[8])
            gene_id_raw = attrs.get("ID") or attrs.get("gene_id")
            gene_id = normalize_gene_id(attrs.get("gene_id") or attrs.get("ID"))
            gene_name = attrs.get("gene_name") or attrs.get("Name")
            genes.append(
                GeneMeta(
                    index=len(genes),
                    seqid=parts[0],
                    gene_id=gene_id,
                    gene_id_raw=gene_id_raw,
                    gene_name=gene_name,
                )
            )

    return genes, gff_version


def pick_target_gene(genes: Sequence[GeneMeta], query: str) -> GeneMeta:
    """Find target gene by ID or name with deterministic fallback."""
    query_norm = normalize_gene_id(query)

    # Priority 1: exact normalized gene ID
    id_hits = [
        g for g in genes
        if g.gene_id == query_norm
        or g.gene_id_raw == query
    ]
    if len(id_hits) == 1:
        return id_hits[0]
    if len(id_hits) > 1:
        return id_hits[0]

    # Priority 2: exact symbol
    name_hits_exact = [g for g in genes if g.gene_name == query]
    if len(name_hits_exact) == 1:
        return name_hits_exact[0]

    # Priority 3: case-insensitive symbol
    query_lower = query.lower()
    name_hits_ci = [
        g for g in genes
        if g.gene_name and g.gene_name.lower() == query_lower
    ]
    if len(name_hits_ci) == 1:
        return name_hits_ci[0]
    if len(name_hits_ci) > 1:
        raise ValueError(
            f"Ambiguous gene name '{query}' ({len(name_hits_ci)} hits). "
            "Please pass a gene ID instead."
        )

    raise ValueError(f"Gene '{query}' not found in {len(genes)} GFF3 genes")


def write_selected_gene_blocks(
    gff_path: Path,
    out_gff: Path,
    selected_indices: Set[int],
    gff_version_line: str,
) -> int:
    """Write selected top-level gene blocks (gene + child lines) to output GFF3.

    Child lines are retained only when they explicitly carry a Parent attribute.
    This avoids pulling unrelated intergenic annotations that may appear between
    two adjacent gene blocks in Ensembl GFF3.
    """
    kept_lines = 0
    current_gene_index = -1
    keep_current_block = False

    with gff_path.open("r", encoding="utf-8") as in_handle, out_gff.open(
        "w", encoding="utf-8"
    ) as out_handle:
        out_handle.write(gff_version_line + "\n")

        for raw in in_handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            is_gene_line = detect_gene_line(parts)
            if is_gene_line:
                current_gene_index += 1
                keep_current_block = current_gene_index in selected_indices

            if not keep_current_block:
                continue

            attrs = parse_attributes(parts[8])
            if is_gene_line or attrs.get("Parent"):
                out_handle.write(line + "\n")
                kept_lines += 1

    return kept_lines


def resolve_target_seq_name(seq_name: str, references: Set[str]) -> Optional[str]:
    """Resolve chromosome naming differences (chr1 vs 1) if present."""
    if seq_name in references:
        return seq_name

    if seq_name.startswith("chr") and seq_name[3:] in references:
        return seq_name[3:]

    candidate = "chr" + seq_name
    if candidate in references:
        return candidate

    return None


def write_chrom_subset_fasta(
    source_fasta: Path,
    out_fasta: Path,
    chromosomes: Sequence[str],
) -> Dict[str, str]:
    """Write whole chromosome sequences for selected chromosomes.

    Returns:
        Mapping of requested chromosome -> source FASTA reference name used.
    """
    resolved: Dict[str, str] = {}

    with pysam.FastaFile(str(source_fasta)) as fasta, out_fasta.open(
        "w", encoding="utf-8"
    ) as out_handle:
        refs = set(fasta.references)
        for chrom in chromosomes:
            ref_name = resolve_target_seq_name(chrom, refs)
            if ref_name is None:
                raise ValueError(
                    f"Chromosome '{chrom}' not found in FASTA {source_fasta.name}"
                )

            seq = fasta.fetch(ref_name)
            resolved[chrom] = ref_name
            out_handle.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                out_handle.write(seq[i : i + 60] + "\n")

    return resolved


def default_outdir_for_gene(base_dir: Path, query: str) -> Path:
    """Build a deterministic output directory name."""
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in query)
    return base_dir / f"mini_{safe}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a minimal local test dataset around one gene",
    )
    parser.add_argument(
        "gene",
        help="Gene symbol or Ensembl gene ID (e.g. NOC2L or ENSG00000188976)",
    )
    parser.add_argument(
        "--flank-genes",
        type=int,
        default=5,
        help="Number of genes before/after target in GFF3 order (default: 5)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: cluster_test_data/Test/mini_<gene>)",
    )
    parser.add_argument(
        "--ref-gff",
        type=Path,
        default=Path("cluster_test_data/Homo_sapiens.GRCh38.115.chr.gff3"),
        help="Reference GFF3 path",
    )
    parser.add_argument(
        "--ref-fasta",
        type=Path,
        default=Path("cluster_test_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        help="Reference FASTA path",
    )
    parser.add_argument(
        "--target-fasta",
        type=Path,
        default=Path("cluster_test_data/CHM13.unmasked.fa"),
        help="Target FASTA path",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output directory if it already exists",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.flank_genes < 0:
        raise ValueError("--flank-genes must be >= 0")

    for path in (args.ref_gff, args.ref_fasta, args.target_fasta):
        if not path.exists():
            raise FileNotFoundError(f"Missing input: {path}")

    out_dir = args.output_dir
    if out_dir is None:
        out_dir = default_outdir_for_gene(Path("cluster_test_data/Test"), args.gene)

    if out_dir.exists():
        if not args.force:
            raise FileExistsError(
                f"Output directory exists: {out_dir} (use --force to overwrite)"
            )
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genes, gff_version = scan_gene_metadata(args.ref_gff)
    if not genes:
        raise RuntimeError("No genes found in GFF3")

    target = pick_target_gene(genes, args.gene)
    start_idx = max(0, target.index - args.flank_genes)
    end_idx = min(len(genes) - 1, target.index + args.flank_genes)
    selected = genes[start_idx : end_idx + 1]
    selected_indices = {g.index for g in selected}

    chromosomes_in_order: List[str] = []
    seen: Set[str] = set()
    for gene in selected:
        if gene.seqid not in seen:
            chromosomes_in_order.append(gene.seqid)
            seen.add(gene.seqid)

    out_gff = out_dir / "test.gff3"
    out_ref = out_dir / "reference.fa"
    out_target = out_dir / "target.fa"
    summary_json = out_dir / "dataset_summary.json"

    kept_lines = write_selected_gene_blocks(
        args.ref_gff,
        out_gff,
        selected_indices,
        gff_version_line=gff_version,
    )

    ref_resolved = write_chrom_subset_fasta(
        args.ref_fasta,
        out_ref,
        chromosomes_in_order,
    )
    target_resolved = write_chrom_subset_fasta(
        args.target_fasta,
        out_target,
        chromosomes_in_order,
    )

    summary = {
        "query": args.gene,
        "target_gene": {
            "index": target.index,
            "seqid": target.seqid,
            "gene_id": target.gene_id,
            "gene_name": target.gene_name,
        },
        "flank_genes": args.flank_genes,
        "selected_gene_count": len(selected),
        "selected_gene_index_range": [start_idx, end_idx],
        "selected_chromosomes": chromosomes_in_order,
        "gff_feature_lines_written": kept_lines,
        "outputs": {
            "gff3": str(out_gff),
            "reference_fasta": str(out_ref),
            "target_fasta": str(out_target),
        },
        "chromosome_name_resolution": {
            "reference": ref_resolved,
            "target": target_resolved,
        },
    }

    with summary_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Created minimal test dataset in: {out_dir}")
    print(f"  Target gene: {target.gene_name or target.gene_id or target.gene_id_raw}")
    print(f"  Gene window: {start_idx}..{end_idx} ({len(selected)} genes)")
    print(f"  Chromosomes: {', '.join(chromosomes_in_order)}")
    print(f"  GFF lines: {kept_lines}")
    print(f"  GFF3: {out_gff}")
    print(f"  REF : {out_ref}")
    print(f"  TGT : {out_target}")
    print(f"  Summary: {summary_json}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

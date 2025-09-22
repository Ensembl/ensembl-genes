#!/usr/bin/env python3
"""
Test script to convert a RefSeq GFF3 file to Ensembl-style naming using an assembly report.
"""

import os
from refseq2ensembl import convert_gff_to_ensembl

def main():
    # Inputs
    base_dir = "refseq_data/GCF/000/001/635/GCF_000001635.27/"
    gff_in = os.path.join(base_dir, "GCF_000001635.27_GRCm39_genomic.gff.gz")
    report_in = os.path.join(base_dir, "GCF_000001635.27_GRCm39_assembly_report.txt")
    gff_out = os.path.join(base_dir, "GCF_000001635.27_GRCm39_ensembl_style_all.gff3")

#    restrict_to = ["10"]
    restrict_to = None

    print("→ Converting GFF3 to Ensembl-style...")
    convert_gff_to_ensembl(
        gff_in,
        report_in,
        gff_out,
        restrict_to
    )
    print(f"[✓] Converted GFF3 written to: {gff_out}")

if __name__ == "__main__":
    main()

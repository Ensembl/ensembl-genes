#!/usr/bin/env python3

"""
Test script for loading a RefSeq-converted GFF3 into an Ensembl-style core DB.
Assumes `load_to_ensembl_core()` exists in refseq2ensembl module.
"""

import os
from refseq2ensembl import load_to_ensembl_core

def main():
    # Paths to required files
#    base_dir = "refseq_data/GCF/000/001/635/GCF_000001635.27/"
    base_dir = "/nfs/production/flicek/ensembl/genebuild/fergal/refseq_loading_test/"
#    gff_path = os.path.join(base_dir, "GCF_000001635.27_GRCm39_genomic.gff")
    gff_path = os.path.join(base_dir, "NC_small.gff3")
#    assembly_report = os.path.join(base_dir, "GCF_000001635.27_GRCm39_assembly_report.txt")
    assembly_report = "/nfs/production/flicek/ensembl/genebuild/fergal/refseq_loading_test/GCF_000001635.27_GRCm39_assembly_report.txt"
#    converted_fna_path = os.path.join(base_dir, "GCF_000001635.27_GRCm39_genomic_simple.fna")
    converted_fna_path = os.path.join(base_dir, "NC_small.fna")
    # Metadata
    species_name = "Mus musculus"
    assembly_accession = "GCF_000001635.27"

    # Database connection (adjust to your environment)
    db_host = "mysql-ens-genebuild-prod-1"
    db_user = "ensadmin"
    db_port = "4527"
    db_password = "ensembl"
    schema_sql_path="/nfs/production/flicek/ensembl/genebuild/fergal/refseq_loading_test/ensembl/sql/table.sql"

    # Summary
    print(f"→ Loading into core DB for {species_name} ({assembly_accession})...")
    print(f"   GFF3: {gff_path}")
    print(f"   Assembly report: {assembly_report}")

    load_to_ensembl_core(
        converted_gff_path=gff_path,
        converted_fna_path=converted_fna_path,
        assembly_report_path=assembly_report,
        species_name=species_name,
        assembly_accession=assembly_accession,
        db_host=db_host,
        db_user=db_user,
        db_password=db_password,
        db_port=db_port,
        schema_sql_path=schema_sql_path
    )

    print("✓ Done loading into core DB.")

if __name__ == "__main__":
    main()

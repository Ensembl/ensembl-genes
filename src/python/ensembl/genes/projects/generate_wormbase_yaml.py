#!/usr/bin/env python3
"""
Main entry point for generating WormBase ParaSite project YAML.

Usage:
  python -m ensembl.genes.projects.generate_wormbase_yaml \
    --release WBPS19 \
    --output wormbase_species.yaml \
    --audit-file wormbase_audit.tsv
"""

import argparse
import sys
import yaml

from ensembl.genes.projects.wormbase_ftp import WormBaseHTTP
from ensembl.genes.projects.wormbase_renderer import WormBaseRenderer


def main():
    parser = argparse.ArgumentParser(description="Generate WormBase ParaSite species.yaml")
    parser.add_argument("--release", required=True, help="WormBase ParaSite release (e.g., WBPS19)")
    parser.add_argument("--output", default="wormbase_species.yaml", help="Output YAML file")
    parser.add_argument("--audit-file", help="Optional output TSV file for audit logs")
    args = parser.parse_args()

    http_client = WormBaseHTTP(args.release)
    renderer = WormBaseRenderer(args.release)

    species_list = http_client.get_species()
    if not species_list:
        print(f"No species found for release {args.release}. Check if the release exists.", file=sys.stderr)
        sys.exit(1)

    yaml_docs = []
    audit_rows = []
    audit_header = "species_slug\tspecies\tbioproject\tftp_dumps\tincluded\tfiles_found\tmissing_expected_files\treason\n"

    for species_slug in species_list:
        bioprojects = http_client.get_bioprojects(species_slug)
        for bioproject in bioprojects:
            files = http_client.get_files(species_slug, bioproject)
            
            row, missing = renderer.render_row(species_slug, bioproject, files)
            
            is_valid = renderer.is_valid_row(row)
            
            files_found = len(files)
            missing_str = ",".join(missing) if missing else "none"
            
            if is_valid:
                yaml_docs.append(row)
                audit_rows.append(
                    f"{species_slug}\t{row['species']}\t{bioproject}\t{row['ftp_dumps']}\t"
                    f"included\t{files_found}\t{missing_str}\tValid core files found\n"
                )
            else:
                audit_rows.append(
                    f"{species_slug}\t{row['species']}\t{bioproject}\t{row['ftp_dumps']}\t"
                    f"excluded\t{files_found}\t{missing_str}\tNo core files found\n"
                )

    # Sort outputs
    yaml_docs.sort(key=lambda x: (x["species"].lower(), x["bioproject"]))

    with open(args.output, "w", encoding="utf-8") as f:
        for i, doc in enumerate(yaml_docs):
            if i > 0:
                f.write("\n")
            yaml.dump([doc], f, default_flow_style=False, sort_keys=False)

    if args.audit_file:
        with open(args.audit_file, "w", encoding="utf-8") as f:
            f.write(audit_header)
            for row in audit_rows:
                f.write(row)

    print(f"Successfully generated {args.output} with {len(yaml_docs)} genomes.")

if __name__ == "__main__":
    main()

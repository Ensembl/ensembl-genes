#!/usr/bin/env python3
"""
Main entry point for generating project-specific genome metadata YAML files.

Usage:
  python generate_project_yaml.py --project vgp input_dbs.txt
"""
import argparse
import sys
import yaml
import json
from pathlib import Path
from typing import List

from ensembl.genes.projects.config import get_project_config
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.yaml_renderer import YamlRenderer
from ensembl.genes.projects.registry.metadata_db import MetadataDbClient
from ensembl.genes.projects.registry.gb_tracker import GbTrackerClient
from ensembl.genes.projects.registry.ncbi_entrez import patch_ncbi_data

def _load_server_config() -> dict:
    config_path = Path(__file__).parent / "server_config.json"
    with open(config_path, "r") as f:
        return json.load(f)

def main():
    parser = argparse.ArgumentParser(description="Generate Ensembl project species.yaml")
    parser.add_argument("--project", required=True, help="Project name (e.g. vgp, hprc, mouse_genomes)")
    parser.add_argument("input_file", help="File containing list of DB names, GUUIDs, or Accessions")
    parser.add_argument("--output", default="species.yaml", help="Output YAML file")
    
    args = parser.parse_args()
    
    config = get_project_config(args.project)
    renderer = YamlRenderer(config)
    server_conf = _load_server_config()
    meta_conf = server_conf["meta_beta"]
    
    metadata_client = MetadataDbClient(
        host=meta_conf["db_host"],
        port=meta_conf["db_port"],
        user=meta_conf["db_user"],
        dbname=meta_conf["db_name"]
    )
    
    gb_client = GbTrackerClient(
        host=meta_conf["db_host"],
        port=meta_conf["db_port"],
        user=meta_conf["db_user"]
    )
    
    import logging
    logger = logging.getLogger(__name__)
    
    # Read inputs
    try:
        with open(args.input_file, 'r') as f:
            identifiers = [line.strip().split('\t')[0] for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)
        
    yaml_docs: List[dict] = []
    failed = []
    
    for identifier in identifiers:
        # Prefer metadata DB
        meta = metadata_client.fetch_by_identifier(identifier)
        
        # Fallback to gb_schema tracking for pre-release
        if not meta:
            meta = gb_client.fetch_by_identifier(identifier)
            if meta:
                logger.info(f"Routed '{identifier}' through gb_schema tracking as pre-release data.")
        
        # Explicit fallback logging for remaining failures requiring core DB
        if not meta:
            logger.warning(
                f"Validation failed for '{identifier}'. Core DB fallback required but deprecated. "
                f"Ensure this accession is indexed in ensembl_metadata_qrp or gb_schema."
            )
            failed.append(identifier)
            continue
            
        patch_ncbi_data(meta, config)
        doc = renderer.render(meta)
        yaml_docs.append(doc)
            
    # Write output
    with open(args.output, 'w') as f:
        yaml.dump(yaml_docs, f, default_flow_style=False, sort_keys=False)
        
    print(f"Successfully generated {args.output} for {len(yaml_docs)} genomes using {config.schema_type} schema format.")
    if failed:
        print(f"Warning: Failed to retrieve metadata for {len(failed)} identifiers: {failed}", file=sys.stderr)

if __name__ == "__main__":
    main()

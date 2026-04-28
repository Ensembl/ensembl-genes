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
from dataclasses import dataclass

@dataclass
class Candidate:
    identifier: str
    meta: GenomeMetadata
    doc: dict
    audit_decision: str
    audit_reason: str
    audit_resolved_date: str
    input_order: int

def _load_server_config() -> dict:
    config_path = Path(__file__).parent / "server_config.json"
    with open(config_path, "r") as f:
        return json.load(f)

def main():
    parser = argparse.ArgumentParser(description="Generate Ensembl project species.yaml")
    parser.add_argument("--project", required=True, help="Project name (e.g. vgp, hprc, mouse_genomes)")
    parser.add_argument("input_file", help="File containing list of DB names, GUUIDs, or Accessions")
    parser.add_argument("--output", default="species.yaml", help="Output YAML file")
    parser.add_argument("--audit-file", help="Optional output TSV file for audit logs")
    
    args = parser.parse_args()
    
    config = get_project_config(args.project)
    from ensembl.genes.projects.write_yaml import EnsemblFTP
    ftp_client = EnsemblFTP(timeout=30)
    renderer = YamlRenderer(config, ftp_client)
    
    server_conf = _load_server_config()
    meta_conf = server_conf["meta_beta"]
    
    metadata_client = MetadataDbClient(
        host=meta_conf["db_host"],
        port=meta_conf["db_port"],
        user=meta_conf["db_user"],
        dbname=meta_conf["db_name"]
    )
    
    gb_conf = server_conf["gb1"]
    gb_client = GbTrackerClient(
        host=gb_conf["db_host"],
        port=gb_conf["db_port"],
        user=gb_conf["db_user"]
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
    emitted_accessions = set()
    
    candidates: List[Candidate] = []
    
    for idx, identifier in enumerate(identifiers):
        if "-" in identifier and len(identifier) == 36:
            # Explicitly a Genome UUID -> Released genome track
            meta = metadata_client.fetch_by_identifier(identifier)
            if not meta:
                logger.warning(
                    f"Validation failed for '{identifier}'. Could not find UUID in metadata DB. "
                    f"Ensure this UUID is indexed in ensembl_metadata_qrp."
                )
                failed.append(identifier)
                continue
        else:
            # Treat as core DB name -> Pre-release genome track
            meta = gb_client.fetch_by_identifier(identifier)
            if meta:
                logger.info(f"Routed '{identifier}' through gb_schema tracking as pre-release data.")
            else:
                logger.warning(
                    f"Validation failed for '{identifier}'. Could not find DB in gb_schema. "
                    f"Core DB fallback is deprecated."
                )
                failed.append(identifier)
                continue
            
        patch_ncbi_data(meta, config)
        doc = renderer.render(meta)
        
        audit_decision = doc.pop("__audit_decision__", "excluded")
        audit_reason = doc.pop("__audit_reason__", "No document returned")
        audit_resolved_date = doc.pop("__audit_resolved_date__", "")
        
        candidates.append(Candidate(
            identifier=identifier,
            meta=meta,
            doc=doc,
            audit_decision=audit_decision,
            audit_reason=audit_reason,
            audit_resolved_date=audit_resolved_date,
            input_order=idx
        ))
        
        if doc and audit_decision != "excluded":
            emitted_accessions.add(meta.accession)
        else:
            # Exclude and warn
            logger.warning(f"Excluding {identifier}: {audit_reason}")
            failed.append(identifier)
            
    # GB Registry Discovery Pass
    gb_candidates = gb_client.fetch_project_pre_releases(config)
    discovered_count = 0
    # continue using input_order offset
    start_idx = len(identifiers)
    for idx, meta in enumerate(gb_candidates, start=start_idx):
        if meta.accession in emitted_accessions:
            continue
            
        patch_ncbi_data(meta, config)
        doc = renderer.render(meta)
        
        audit_decision = doc.pop("__audit_decision__", "excluded")
        audit_reason = doc.pop("__audit_reason__", "No document returned")
        audit_resolved_date = doc.pop("__audit_resolved_date__", "")
        
        identifier = f"discovered_gb_{meta.accession}"
        candidates.append(Candidate(
            identifier=identifier,
            meta=meta,
            doc=doc,
            audit_decision=audit_decision,
            audit_reason=audit_reason,
            audit_resolved_date=audit_resolved_date,
            input_order=idx
        ))
        
        if doc and audit_decision != "excluded":
            emitted_accessions.add(meta.accession)
            discovered_count += 1
            logger.info(f"Discovered and appended GB-only pre-release: {meta.accession}")
        else:
            logger.debug(f"Discovered GB-only pre-release {meta.accession} skipped: {audit_reason}")
            
    from collections import defaultdict
    
    grouped_candidates = defaultdict(list)
    for c in candidates:
        if c.audit_decision == "excluded":
            continue
        species_name = c.doc.get("species", c.doc.get("assembly", c.meta.species_name))
        source = (c.meta.annotation_source or "ensembl").lower().strip()
        key = (c.meta.accession, species_name, source)
        grouped_candidates[key].append(c)
        
    for key, group in grouped_candidates.items():
        if len(group) == 1:
            yaml_docs.append(group[0].doc)
            continue
            
        subgroups = defaultdict(list)
        for c in group:
            subgroups[(c.audit_resolved_date, c.doc.get("annotation_gtf", ""))].append(c)
            
        for subkey, subgroup in subgroups.items():
            if len(subgroup) == 1:
                yaml_docs.append(subgroup[0].doc)
                continue
                
            def score(c: Candidate):
                meta_date_normalized = c.meta.annotation_date.replace("-", "_") if c.meta.annotation_date else ""
                matches_date = (meta_date_normalized == c.audit_resolved_date)
                is_released = c.meta.is_released
                return (-int(matches_date), -int(is_released), c.input_order)
                
            subgroup.sort(key=score)
            kept = subgroup[0]
            kept.audit_decision = "kept_duplicate"
            kept.audit_reason = "Kept preferred duplicate"
            yaml_docs.append(kept.doc)
            
            for excluded in subgroup[1:]:
                excluded.audit_decision = "excluded_duplicate"
                excluded.audit_reason = f"same accession/species/source/date as {kept.identifier}"
                
    if args.audit_file:
        with open(args.audit_file, "a") as af:
            for c in candidates:
                af.write(f"{c.identifier}\t{c.meta.accession}\t{c.meta.species_name}\t{c.meta.annotation_date}\t{c.audit_resolved_date}\t{c.meta.annotation_source}\t{c.audit_decision}\t{c.audit_reason}\n")
            
    # Write output
    yaml_docs.sort(key=lambda x: x.get('species', x.get('assembly', '')))
    
    with open(args.output, 'w') as f:
        for i, doc in enumerate(yaml_docs):
            if i > 0:
                f.write("\n")
            yaml.dump([doc], f, default_flow_style=False, sort_keys=False)
        
    print(f"Successfully generated {args.output} for {len(yaml_docs)} genomes using {config.schema_type} schema format.")
    if failed:
        print(f"Warning: Failed to retrieve metadata for {len(failed)} identifiers: {failed}", file=sys.stderr)

if __name__ == "__main__":
    main()

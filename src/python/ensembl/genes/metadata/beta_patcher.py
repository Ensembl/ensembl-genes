#!/usr/bin/env python3
"""
Beta Metadata Patcher Script

Generates standardized SQL patch files for fixing metadata issues in:
1. Production metadata database (ensembl_genome_metadata)
2. Core database (meta table)

This ensures metadata consistency between both databases and provides
a standardized workflow for creating and applying patches.

Requires genome_uuid for all patches.

Usage:
    python beta_patcher.py patches.csv --jira-ticket EBD-1111 --output-dir ./patches/
    python beta_patcher.py patches.csv --jira-ticket EBD-1111 --core-suffix _core_115_1
"""

import argparse
import csv
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from ensembl.production.metadata.api.adaptors.genome import GenomeAdaptor
    METADATA_API_AVAILABLE = True
except ImportError:
    METADATA_API_AVAILABLE = False
    logging.warning(
        "ensembl-metadata-api not available. Install with:\n"
        "  git clone https://github.com/Ensembl/ensembl-metadata-api.git\n"
        "  cd ensembl-metadata-api\n"
        "  pip install -r requirements.txt && pip install -e ."
    )


def setup_logging(output_dir: Path, genome_uuid: str) -> logging.Logger:
    """Set up logging configuration."""
    log_file = output_dir / f"patch_{genome_uuid}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def get_genome_adaptor() -> Optional[GenomeAdaptor]:
    """
    Create a GenomeAdaptor instance using environment variables.

    Requires:
        METADATA_URI: Connection string for metadata database
        TAXONOMY_URI: Connection string for taxonomy database

    Returns:
        GenomeAdaptor instance or None if not available
    """
    if not METADATA_API_AVAILABLE:
        logging.error("ensembl-metadata-api package not available")
        return None

    metadata_uri = os.getenv('METADATA_URI')
    taxonomy_uri = os.getenv('TAXONOMY_URI')

    if not metadata_uri or not taxonomy_uri:
        logging.error(
            "METADATA_URI and TAXONOMY_URI environment variables required. "
            "Format: mysql+pymysql://user:pass@host:port/database"
        )
        return None

    try:
        return GenomeAdaptor(metadata_uri=metadata_uri, taxonomy_uri=taxonomy_uri)
    except Exception as e:
        logging.error(f"Failed to create GenomeAdaptor: {e}")
        return None


def get_genome_by_uuid(genome_uuid: str) -> Optional[Dict]:
    """
    Fetch genome details by UUID using ensembl-metadata-api.

    Args:
        genome_uuid: Genome UUID to query

    Returns:
        Dict containing genome information or None if not found
    """
    adaptor = get_genome_adaptor()
    if not adaptor:
        return None

    try:
        results = adaptor.fetch_genomes_by_genome_uuid(genome_uuid=genome_uuid)
        if results:
            genome, organism, assembly, release, site = results[0]
            return {
                "genome_uuid": genome.genome_uuid,
                "production_name": genome.production_name,
                "genebuild_version": genome.genebuild_version,
                "genebuild_date": genome.genebuild_date,
                "organism": {
                    "organism_uuid": organism.organism_uuid,
                    "taxonomy_id": organism.taxonomy_id,
                    "scientific_name": organism.scientific_name,
                    "strain": organism.strain,
                    "biosample_id": organism.biosample_id
                },
                "assembly": {
                    "assembly_uuid": assembly.assembly_uuid,
                    "accession": assembly.accession,
                    "name": assembly.name,
                    "level": assembly.level
                }
            }
        return None
    except Exception as e:
        logging.error(f"Failed to fetch genome {genome_uuid}: {e}")
        return None


def get_team_responsible_for_genome(genome_uuid: str, dataset_type: str = "genebuild") -> Optional[str]:
    """
    Fetch genebuild.team_responsible attribute for a genome.

    Args:
        genome_uuid: Genome UUID to query
        dataset_type: Dataset type (default: genebuild)

    Returns:
        Team responsible string or None if not found
    """
    adaptor = get_genome_adaptor()
    if not adaptor:
        return None

    try:
        # Use the metadata API to fetch the attribute value
        # This queries: dataset_attribute JOIN attribute JOIN dataset JOIN genome_dataset JOIN genome
        from sqlalchemy import text

        query = text("""
            SELECT dataset_attribute.value
            FROM dataset_attribute
            JOIN attribute ON dataset_attribute.attribute_id = attribute.attribute_id
            JOIN dataset ON dataset_attribute.dataset_id = dataset.dataset_id
            JOIN genome_dataset ON dataset.dataset_id = genome_dataset.dataset_id
            JOIN genome ON genome_dataset.genome_id = genome.genome_id
            WHERE genome.genome_uuid = :genome_uuid
            AND dataset.name = :dataset_type
            AND attribute.name = 'genebuild.team_responsible'
        """)

        with adaptor.metadata_db.connect() as conn:
            result = conn.execute(query, {"genome_uuid": genome_uuid, "dataset_type": dataset_type})
            row = result.fetchone()
            if row:
                return row[0]
        return None
    except Exception as e:
        logging.warning(f"Failed to fetch team_responsible for genome {genome_uuid}: {e}")
        return None


def get_affected_genomes_and_teams(genome_uuid: str, table_location: str, dataset_type: str = "genebuild") -> List[Dict]:
    """
    Get all genomes that share the same organism or assembly, and their team_responsible values.

    Args:
        genome_uuid: Source genome UUID
        table_location: Table being modified (organism or assembly)
        dataset_type: Dataset type (default: genebuild)

    Returns:
        List of dicts with genome_uuid, production_name, and team_responsible
    """
    if table_location not in ['organism', 'assembly']:
        return []

    adaptor = get_genome_adaptor()
    if not adaptor:
        return []

    try:
        from sqlalchemy import text

        if table_location == 'organism':
            # Find all genomes sharing the same organism
            query = text("""
                SELECT g2.genome_uuid, g2.production_name
                FROM genome g1
                JOIN genome g2 ON g1.organism_id = g2.organism_id
                WHERE g1.genome_uuid = :genome_uuid
            """)
        else:  # assembly
            # Find all genomes sharing the same assembly
            query = text("""
                SELECT g2.genome_uuid, g2.production_name
                FROM genome g1
                JOIN genome g2 ON g1.assembly_id = g2.assembly_id
                WHERE g1.genome_uuid = :genome_uuid
            """)

        affected_genomes = []
        with adaptor.metadata_db.connect() as conn:
            result = conn.execute(query, {"genome_uuid": genome_uuid})
            for row in result:
                uuid, prod_name = row
                team = get_team_responsible_for_genome(uuid, dataset_type)
                affected_genomes.append({
                    "genome_uuid": uuid,
                    "production_name": prod_name,
                    "team_responsible": team or "UNKNOWN"
                })

        return affected_genomes
    except Exception as e:
        logging.warning(f"Failed to fetch affected genomes for {genome_uuid}: {e}")
        return []



#Helper functions to ensure SELECTs actually validate UPDATEs
def _metadata_db_joins() -> str:
    """Generate standardized JOIN clauses for metadata DB queries."""
    return (
        "JOIN attribute ON dataset_attribute.attribute_id = attribute.attribute_id\n"
        "JOIN dataset ON dataset_attribute.dataset_id = dataset.dataset_id\n"
        "JOIN genome_dataset ON dataset.dataset_id = genome_dataset.dataset_id\n"
        "JOIN genome ON genome_dataset.genome_id = genome.genome_id"
    )


def _metadata_db_where(genome_uuid: str, dataset_type: str, attribute_name: str) -> str:
    """Generate standardized WHERE clause for metadata DB queries."""
    return (
        f"WHERE genome.genome_uuid = '{genome_uuid}'\n"
        f"AND dataset.name = '{dataset_type}'\n"
        f"AND attribute.name = '{attribute_name}'"
    )


def write_metadata_patch_for_genome(
    validate_file,
    patch_file,
    genome_uuid: str,
    patches: List[Tuple[str, str, str]],
    dataset_type: str = "genebuild",
    logger: logging.Logger = None,
    team_filter: str = None
):
    """
    Write validation and patch SQL for a single genome to open file handles.

    Args:
        validate_file: Open file handle for validation SQL
        patch_file: Open file handle for patch SQL
        genome_uuid: Genome UUID
        patches: List of (attribute_name, new_value, table_location) tuples
        dataset_type: Dataset type (genebuild, assembly, etc.)
        logger: Logger instance for warnings
        team_filter: Optional team name to filter by (e.g., 'genebuild')
    """
    # Valid table locations
    VALID_TABLES = {'dataset_attribute', 'genome', 'organism', 'assembly'}

    for attribute_name, new_value, table_location in patches:
        escaped_value = new_value.replace("'", "''")

        # Validate table_location
        if table_location not in VALID_TABLES:
            validate_file.write(f"-- ERROR: Invalid table_location '{table_location}' for {attribute_name}\n")
            validate_file.write(f"-- Valid tables: {', '.join(VALID_TABLES)}\n\n")
            patch_file.write(f"-- ERROR: Invalid table_location '{table_location}' for {attribute_name}\n")
            patch_file.write(f"-- Valid tables: {', '.join(VALID_TABLES)}\n\n")
            continue

        # Check for organism/assembly modifications and warn about affected genomes
        if table_location in ['organism', 'assembly']:
            affected_genomes = get_affected_genomes_and_teams(genome_uuid, table_location, dataset_type)
            if affected_genomes and len(affected_genomes) > 0:
                # Group by team_responsible
                teams_map = {}
                for genome_data in affected_genomes:
                    team = genome_data['team_responsible']
                    if team not in teams_map:
                        teams_map[team] = []
                    teams_map[team].append(f"{genome_data['production_name']} ({genome_data['genome_uuid']})")

                # Check team filter - skip if not all genomes match the filter
                if team_filter:
                    teams_set = set(teams_map.keys())
                    if teams_set != {team_filter}:
                        skip_msg = f"-- SKIPPED: {table_location} modification affects teams {teams_set}, filter requires only '{team_filter}'\n"
                        validate_file.write(skip_msg)
                        patch_file.write(skip_msg)
                        for team, genomes in teams_map.items():
                            team_line = f"--   Team: {team}\n"
                            validate_file.write(team_line)
                            patch_file.write(team_line)
                            for genome_str in genomes:
                                genome_line = f"--     - {genome_str}\n"
                                validate_file.write(genome_line)
                                patch_file.write(genome_line)
                        validate_file.write("\n")
                        patch_file.write("\n")
                        if logger:
                            logger.warning(f"SKIPPED: {attribute_name} for {genome_uuid} - affects teams {teams_set}, filter requires '{team_filter}'")
                        continue

                # Write warning to both files
                warning_header = f"-- WARNING: Modifying shared {table_location} affects {len(affected_genomes)} genome(s)\n"
                validate_file.write(warning_header)
                patch_file.write(warning_header)

                for team, genomes in teams_map.items():
                    team_line = f"--   Team: {team}\n"
                    validate_file.write(team_line)
                    patch_file.write(team_line)
                    for genome_str in genomes:
                        genome_line = f"--     - {genome_str}\n"
                        validate_file.write(genome_line)
                        patch_file.write(genome_line)

                if team_filter:
                    validate_file.write(f"-- Team filter '{team_filter}' PASSED - all affected genomes owned by this team\n")
                    patch_file.write(f"-- Team filter '{team_filter}' PASSED - all affected genomes owned by this team\n")
                else:
                    validate_file.write("-- ENSURE CONSULTATION WITH ALL AFFECTED TEAMS BEFORE APPLYING\n")
                    patch_file.write("-- ENSURE CONSULTATION WITH ALL AFFECTED TEAMS BEFORE APPLYING\n")

                # Log warning to CLI
                if logger:
                    logger.warning(f"Modifying {table_location} for {genome_uuid} affects {len(affected_genomes)} genome(s)")
                    for team, genomes in teams_map.items():
                        logger.warning(f"  Team {team}: {len(genomes)} genome(s)")
                    if team_filter:
                        logger.info(f"Team filter '{team_filter}' passed for {attribute_name}")

        # Write validation SQL
        validate_file.write(f"-- {genome_uuid} | {attribute_name} (table: {table_location})\n")
        if table_location == 'dataset_attribute':
            validate_file.write("SELECT genome.genome_uuid, dataset_attribute.value AS current_value, ")
            validate_file.write(f"'{escaped_value}' AS proposed_value\n")
            validate_file.write("FROM dataset_attribute\n")
            validate_file.write(_metadata_db_joins())
            validate_file.write("\n")
            validate_file.write(_metadata_db_where(genome_uuid, dataset_type, attribute_name))
            validate_file.write(";\n\n")
        elif table_location == 'genome':
            column_name = attribute_name.split('.')[-1]
            validate_file.write(f"SELECT genome.genome_uuid, genome.{column_name} AS current_value, ")
            validate_file.write(f"'{escaped_value}' AS proposed_value\n")
            validate_file.write("FROM genome\n")
            validate_file.write(f"WHERE genome.genome_uuid = '{genome_uuid}';\n\n")
        elif table_location == 'organism':
            column_name = attribute_name.split('.')[-1]
            validate_file.write(f"SELECT genome.genome_uuid, organism.{column_name} AS current_value, ")
            validate_file.write(f"'{escaped_value}' AS proposed_value\n")
            validate_file.write("FROM organism\n")
            validate_file.write("JOIN genome ON organism.organism_id = genome.organism_id\n")
            validate_file.write(f"WHERE genome.genome_uuid = '{genome_uuid}';\n\n")
        elif table_location == 'assembly':
            column_name = attribute_name.split('.')[-1]
            validate_file.write(f"SELECT genome.genome_uuid, assembly.{column_name} AS current_value, ")
            validate_file.write(f"'{escaped_value}' AS proposed_value\n")
            validate_file.write("FROM assembly\n")
            validate_file.write("JOIN genome ON assembly.assembly_id = genome.assembly_id\n")
            validate_file.write(f"WHERE genome.genome_uuid = '{genome_uuid}';\n\n")

        # Write patch SQL
        patch_file.write(f"-- {genome_uuid} | {attribute_name} (table: {table_location})\n")
        if table_location == 'dataset_attribute':
            # DELETE existing attribute value (if exists)
            patch_file.write("DELETE dataset_attribute FROM dataset_attribute\n")
            patch_file.write(_metadata_db_joins())
            patch_file.write("\n")
            patch_file.write(_metadata_db_where(genome_uuid, dataset_type, attribute_name))
            patch_file.write(";\n\n")

            # INSERT new attribute value
            patch_file.write("INSERT INTO dataset_attribute (dataset_id, attribute_id, value)\n")
            patch_file.write("SELECT dataset.dataset_id, attribute.attribute_id, ")
            patch_file.write(f"'{escaped_value}'\n")
            patch_file.write("FROM dataset\n")
            patch_file.write("JOIN genome_dataset ON dataset.dataset_id = genome_dataset.dataset_id\n")
            patch_file.write("JOIN genome ON genome_dataset.genome_id = genome.genome_id\n")
            patch_file.write("JOIN attribute ON attribute.name = ")
            patch_file.write(f"'{attribute_name}'\n")
            patch_file.write(f"WHERE genome.genome_uuid = '{genome_uuid}'\n")
            patch_file.write(f"AND dataset.name = '{dataset_type}'")
            patch_file.write(";\n\n")

        elif table_location == 'genome':
            column_name = attribute_name.split('.')[-1]
            patch_file.write("UPDATE genome\n")
            patch_file.write(f"SET {column_name} = '{escaped_value}'\n")
            patch_file.write(f"WHERE genome_uuid = '{genome_uuid}';\n\n")

        elif table_location == 'organism':
            column_name = attribute_name.split('.')[-1]
            patch_file.write("UPDATE organism\n")
            patch_file.write("JOIN genome ON organism.organism_id = genome.organism_id\n")
            patch_file.write(f"SET organism.{column_name} = '{escaped_value}'\n")
            patch_file.write(f"WHERE genome.genome_uuid = '{genome_uuid}';\n\n")

        elif table_location == 'assembly':
            column_name = attribute_name.split('.')[-1]
            patch_file.write("UPDATE assembly\n")
            patch_file.write("JOIN genome ON assembly.assembly_id = genome.assembly_id\n")
            patch_file.write(f"SET assembly.{column_name} = '{escaped_value}'\n")
            patch_file.write(f"WHERE genome.genome_uuid = '{genome_uuid}';\n\n")


def _core_db_where(meta_key: str, species_id: int) -> str:
    """Generate standardized WHERE clause for core DB queries."""
    return (
        f"WHERE meta_key = '{meta_key}'\n"
        f"AND species_id = {species_id}"
    )


def write_core_patch_for_genome(
    validate_file,
    patch_file,
    database: str,
    patches: List[Tuple[str, str, str]],
    species_id: int = 1
):
    """
    Write validation and patch SQL for a single core database to open file handles.

    Args:
        validate_file: Open file handle for validation SQL
        patch_file: Open file handle for patch SQL
        database: Core database name
        patches: List of (meta_key, new_value, table_location) tuples (table_location ignored for core DB)
        species_id: Species ID
    """
    for meta_key, new_value, _table_location in patches:
        escaped_value = new_value.replace("'", "''")

        # Write validation SQL
        validate_file.write(f"-- {database} | {meta_key}\n")
        validate_file.write(f"USE {database};\n")
        validate_file.write(f"SELECT '{database}' AS database_name, ")
        validate_file.write(f"'{meta_key}' AS meta_key, ")
        validate_file.write("meta_value AS current_value, ")
        validate_file.write(f"'{escaped_value}' AS proposed_value\n")
        validate_file.write("FROM meta\n")
        validate_file.write(_core_db_where(meta_key, species_id))
        validate_file.write(";\n\n")

        # Write patch SQL
        patch_file.write(f"-- {database} | {meta_key}\n")
        patch_file.write(f"USE {database};\n")
        patch_file.write("UPDATE meta\n")
        patch_file.write(f"SET meta_value = '{escaped_value}'\n")
        patch_file.write(_core_db_where(meta_key, species_id))
        patch_file.write(";\n\n")
        patch_file.write(f"INSERT IGNORE INTO meta (species_id, meta_key, meta_value)\n")
        patch_file.write(f"VALUES ({species_id}, '{meta_key}', '{escaped_value}');\n\n")


def read_csv_patches(csv_file: Path, core_suffix: str = "_core_114_1") -> List[Dict]:
    """
    Read patches from CSV file.

    CSV columns:
    - genome_uuid (required)
    - meta_key (required)
    - desired_meta_value (required)
    - dataset_type (optional, default: genebuild)
    - species_id (optional, default: 1)
    - table_location (optional, default: dataset_attribute)
      Valid values: dataset_attribute, genome, organism, assembly

    Args:
        csv_file: Path to CSV file
        core_suffix: Suffix to append to production_name for core DB name

    Returns:
        List of patch dictionaries

    Raises:
        ValueError if CSV format is invalid
    """
    patches = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)

        # Validate required columns
        required_cols = {'genome_uuid', 'meta_key', 'desired_meta_value'}

        if not required_cols.issubset(reader.fieldnames):
            raise ValueError(f"CSV must contain columns: {required_cols}")

        for row_num, row in enumerate(reader, start=2):  # Start at 2 for header + 1-indexed
            # Check required fields
            if not row.get('genome_uuid'):
                logging.error(f"Row {row_num}: Missing required field genome_uuid")
                continue

            if not row.get('meta_key') or not row.get('desired_meta_value'):
                logging.error(f"Row {row_num}: Missing required fields meta_key or desired_meta_value")
                continue

            patch = {
                'genome_uuid': row['genome_uuid'].strip(),
                'meta_key': row['meta_key'].strip(),
                'desired_meta_value': row['desired_meta_value'].strip(),
                'dataset_type': row.get('dataset_type', 'genebuild').strip(),
                'species_id': int(row.get('species_id', 1)),
                'table_location': row.get('table_location', 'dataset_attribute').strip(),
                'core_suffix': core_suffix,
                'row_num': row_num
            }

            patches.append(patch)

    return patches


def check_thoas_requirements(patches: List[Dict], logger: logging.Logger) -> bool:
    """
    Check if any patches require THOAS (taxonomic heritage) updates.

    Args:
        patches: List of patch dictionaries
        logger: Logger instance
    Returns:
        True if THOAS updates are required, False otherwise
    """
    thoas_loaded_fields = {
        'assembly.accession',
        'assembly.date',
        'assembly.name',
        'assembly.tol_id',
        'organism.taxonomy_id',
        'organism.species_taxonomy_id',
        'organism.common_name',
        'organism.scientific_name',
        'organism.scientific_parlance_name',
        'genome.genome_uuid',
    }

    for patch in patches:
        if patch['meta_key'] in thoas_loaded_fields:
            logger.info(f"Patch for genome_uuid {patch['genome_uuid']} requires THOAS update")
            logger.info(f"  Thoas-loaded field: {patch['meta_key']}")
            return True
    return False


def resolve_genome_info(patch: Dict, logger: logging.Logger) -> Optional[Dict]:
    """
    Resolve production name from genome UUID.

    Args:
        patch: Patch dictionary from CSV
        logger: Logger instance

    Returns:
        Dict with genome_uuid, production_name, core_db_name, or None if failed
    """
    row_num = patch['row_num']
    genome_uuid = patch['genome_uuid']

    logger.info(f"Row {row_num}: Fetching production_name for genome_uuid: {genome_uuid}")
    genome_info = get_genome_by_uuid(genome_uuid)

    if not genome_info:
        logger.error(f"Row {row_num}: Could not fetch genome information for UUID: {genome_uuid}")
        return None

    production_name = genome_info['production_name']
    logger.info(f"Row {row_num}: Found production_name: {production_name}")

    # Build core database name
    core_db_name = f"{production_name}{patch['core_suffix']}"

    return {
        'genome_uuid': genome_uuid,
        'production_name': production_name,
        'core_db_name': core_db_name
    }


def group_patches_by_genome(patches: List[Dict], logger: logging.Logger, jira_ticket: str = "") -> Dict[str, Dict]:
    """
    Group patches by genome UUID.

    Args:
        patches: List of patch dictionaries from CSV
        logger: Logger instance
        jira_ticket: Jira ticket reference

    Returns:
        Dict mapping genome_uuid to genome info and list of patches
    """
    grouped = {}

    for patch in patches:
        # Resolve genome information
        genome_info = resolve_genome_info(patch, logger)

        if not genome_info:
            logger.warning(f"Row {patch['row_num']}: Skipping due to resolution failure")
            continue

        genome_uuid = genome_info['genome_uuid']

        # Initialize genome entry if not exists
        if genome_uuid not in grouped:
            grouped[genome_uuid] = {
                'genome_uuid': genome_uuid,
                'production_name': genome_info['production_name'],
                'core_db_name': genome_info['core_db_name'],
                'dataset_type': patch['dataset_type'],
                'species_id': patch['species_id'],
                'jira_ticket': jira_ticket,
                'patches': [],
                'row_nums': []
            }

        # Add patch to genome's list (include table_location)
        grouped[genome_uuid]['patches'].append((patch['meta_key'], patch['desired_meta_value'], patch['table_location']))
        grouped[genome_uuid]['row_nums'].append(patch['row_num'])

    return grouped


def generate_all_patches(
    grouped_patches: Dict[str, Dict],
    output_dir: Path,
    jira_ticket: str,
    logger: logging.Logger,
    team_filter: str = None
) -> bool:
    """
    Generate consolidated patch files for all genomes.

    Args:
        grouped_patches: Dictionary of genome patches grouped by UUID
        output_dir: Output directory for patch files
        jira_ticket: Jira ticket reference
        logger: Logger instance
        team_filter: Optional team name to filter organism/assembly patches

    Returns:
        True if successful, False otherwise
    """

    # Metadata DB files
    metadata_patch_file = output_dir / f"patch_metadata_{jira_ticket}.sql"
    metadata_validate_file = output_dir / f"validate_metadata_{jira_ticket}.sql"

    # Core DB files
    core_patch_file = output_dir / f"patch_core_{jira_ticket}.sql"
    core_validate_file = output_dir / f"validate_core_{jira_ticket}.sql"

    try:
        # Write metadata patches
        with open(metadata_validate_file, 'w') as val_f, open(metadata_patch_file, 'w') as patch_f:
            val_f.write(f"-- Validation: Metadata DB | {jira_ticket} | {datetime.now().isoformat()}\n")
            val_f.write("USE ensembl_genome_metadata;\n\n")

            patch_f.write(f"-- Patch: Metadata DB | {jira_ticket} | {datetime.now().isoformat()}\n")
            patch_f.write(f"-- Validate first: {metadata_validate_file.name}\n")
            patch_f.write("USE ensembl_genome_metadata;\n\n")

            for genome_uuid, genome_data in grouped_patches.items():
                production_name = genome_data['production_name']
                patches = genome_data['patches']
                logger.info(f"  Metadata: {production_name} ({genome_uuid}): {len(patches)} patches")

                write_metadata_patch_for_genome(
                    val_f, patch_f, genome_uuid, patches, genome_data['dataset_type'], logger, team_filter
                )

        # Write core patches
        with open(core_validate_file, 'w') as val_f, open(core_patch_file, 'w') as patch_f:
            val_f.write(f"-- Validation: Core DBs | {jira_ticket} | {datetime.now().isoformat()}\n\n")

            patch_f.write(f"-- Patch: Core DBs | {jira_ticket} | {datetime.now().isoformat()}\n")
            patch_f.write(f"-- Validate first: {core_validate_file.name}\n\n")

            for genome_uuid, genome_data in grouped_patches.items():
                production_name = genome_data['production_name']
                core_db_name = genome_data['core_db_name']
                patches = genome_data['patches']
                logger.info(f"  Core: {core_db_name}: {len(patches)} patches")

                write_core_patch_for_genome(
                    val_f, patch_f, core_db_name, patches, genome_data['species_id']
                )

        logger.info(f"Generated files:")
        logger.info(f"  {metadata_validate_file}")
        logger.info(f"  {metadata_patch_file}")
        logger.info(f"  {core_validate_file}")
        logger.info(f"  {core_patch_file}")

        return True

    except Exception as e:
        logger.error(f"Failed to generate patches: {e}")
        return False


def main():
    """Main entry point for the beta patcher script."""
    parser = argparse.ArgumentParser(
        description="Generate standardized SQL patches for Ensembl beta metadata issues from CSV input",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process patches from CSV file
  python beta_patcher.py patches.csv --jira-ticket EBD-1111 --output-dir ./patches/

  # Process CSV with custom core suffix
  python beta_patcher.py patches.csv --jira-ticket EBD-1111 --core-suffix _core_115_1

  # Validation files are automatically generated
  # Run validate_*.sql files before applying patch_*.sql files

CSV Format:
  Required columns: genome_uuid, meta_key, desired_meta_value
  Optional columns: dataset_type, species_id, table_location

  table_location valid values: dataset_attribute (default), genome, organism, assembly

  Example:
    genome_uuid,meta_key,desired_meta_value,dataset_type,species_id,table_location
    a7335667-93e7-11ec-a39d-005056b38ce3,assembly.name,GRCh38.p14,genebuild,1,dataset_attribute
    a7335667-93e7-11ec-a39d-005056b38ce3,genebuild_version,2024-01,genebuild,1,genome
    a7335667-93e7-11ec-a39d-005056b38ce3,strain,reference,genebuild,1,organism
    a7335667-93e7-11ec-a39d-005056b38ce3,level,chromosome,genebuild,1,assembly

See patches_template.csv for a complete example.
        """
    )

    # Positional argument
    parser.add_argument(
        'csv_file',
        type=Path,
        help='CSV file with patches (see patches_template.csv for format)'
    )

    # Jira ticket (required)
    parser.add_argument(
        '--jira-ticket',
        type=str,
        required=True,
        help='Jira ticket reference (e.g., EBD-1111)'
    )

    # Core database suffix
    parser.add_argument(
        '--core-suffix',
        type=str,
        default='_core_114_1',
        help='Suffix to append to production_name for core DB name (default: _core_114_1)'
    )

    # Database connection (uses environment variables METADATA_URI and TAXONOMY_URI)
    parser.add_argument(
        '--metadata-uri',
        type=str,
        help='Metadata database URI (overrides METADATA_URI env var)'
    )

    parser.add_argument(
        '--taxonomy-uri',
        type=str,
        help='Taxonomy database URI (overrides TAXONOMY_URI env var)'
    )

    # Output options
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        default=Path.cwd(),
        help='Output directory for patch files (default: current directory)'
    )

    # Team filter option
    parser.add_argument(
        '--team-filter',
        type=str,
        help='Only apply organism/assembly patches if all affected genomes belong to this team (e.g., "genebuild")'
    )

    args = parser.parse_args()

    # Validate CSV file exists
    if not args.csv_file.exists():
        parser.error(f"CSV file not found: {args.csv_file}")

    # Override environment variables if provided
    if args.metadata_uri:
        os.environ['METADATA_URI'] = args.metadata_uri
    if args.taxonomy_uri:
        os.environ['TAXONOMY_URI'] = args.taxonomy_uri

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Set up logging
    log_file = args.output_dir / f"patch_csv_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    try:
        patches = read_csv_patches(args.csv_file, args.core_suffix)
        logger.info(f"Processing {len(patches)} patches from {args.csv_file}")
    except Exception as e:
        logger.error(f"Failed to read CSV: {e}")
        return 1
    
    if check_thoas_requirements(patches, logger):
        logger.error("One or more patches require THOAS updates. Please handle these separately.")
        return 1

    grouped_patches = group_patches_by_genome(patches, logger, args.jira_ticket)
    logger.info(f"Grouped into {len(grouped_patches)} genome(s)")

    if not grouped_patches:
        logger.error("No valid patches to process after grouping")
        return 1

    # Generate consolidated patch files
    team_filter = getattr(args, 'team_filter', None)
    if generate_all_patches(grouped_patches, args.output_dir, args.jira_ticket, logger, team_filter):
        logger.info(f"Success! Log: {log_file}")
        return 0
    else:
        logger.error(f"Failed. Log: {log_file}")
        return 1


if __name__ == '__main__':
    sys.exit(main())

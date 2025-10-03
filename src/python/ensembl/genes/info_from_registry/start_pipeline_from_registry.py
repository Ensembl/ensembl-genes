"""
This script connects to a genomic assembly registry, extracts assembly metadata,
and initializes annotation pipelines based on predefined configurations.

It supports both database-driven metadata retrieval and custom INI file loading,
handles directory setup, and creates necessary configuration and command files
for running genome annotation pipelines.

Typical usage:
    python start_pipeline_from_registry.py --gcas gca_list.txt ---settings_file settings.json
"""

import os
from pathlib import Path
import json
import sys
from typing import Optional, Dict, Any
import shutil
import pymysql # type: ignore
import argparse
import logging
from build_anno_commands import (build_annotation_commands)
from check_if_annotated import (check_if_annotated)
from mysql_helper import mysql_fetch_data
from taxonomy_helper import (assign_clade,assign_clade_info_custom_loading, get_parent_taxon)
from create_pipe_reg import create_registry_entry
from create_config import edit_config_anno, edit_config_main
from assign_species_prefix import get_species_prefix # type: ignore
from assign_stable_space import get_stable_space
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from metrics.busco_lineage_selector import get_dataset_match

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline_setup.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def load_settings(settings_file) -> dict:
    """
    Load JSON-formatted pipeline settings from a given file path.

    Args:
        settings_file (str): Path to the JSON settings file.

    Returns:
        dict: Parsed settings dictionary.

    Raises:
        FileNotFoundError: If the settings file does not exist.
        json.JSONDecodeError: If the file contents are not valid JSON.
    """

    logger.info(f"Loading settings from: {settings_file}")
    settings_path = Path(settings_file)
    if not settings_path.exists():
        raise FileNotFoundError(f"Settings file not found: {settings_path}")
    with open(settings_path, "r") as f:
        return json.load(f)


def load_anno_settings() -> dict:
    """
    Load the annotation-specific settings from a hardcoded path relative to
    the environment variable 'ENSCODE'.

    Returns:
        dict: Parsed annotation settings dictionary.

    Raises:
        FileNotFoundError: If the anno_settings.json file is not found.
        json.JSONDecodeError: If the file contents are not valid JSON.
    """
    logger.info("Loading anno settings json")
    settings = os.path.join(
         os.environ.get("ENSCODE"),
         "ensembl-genes",
         "src",
         "python",
         "ensembl",
         "genes",
         "info_from_registry",
         "anno_settings.json"
     )
    with open(settings, "r") as f:
        return json.load(f)

def load_main_settings() -> dict:
    """
    Load the annotation-specific settings from a hardcoded path relative to
    the environment variable 'ENSCODE'.

    Returns:
        dict: Parsed annotation settings dictionary.

    Raises:
        FileNotFoundError: If the anno_settings.json file is not found.
        json.JSONDecodeError: If the file contents are not valid JSON.
    """
    logger.info("Loading main settings json")
    settings = os.path.join(
         os.environ.get("ENSCODE"),
         "ensembl-genes",
         "src",
         "python",
         "ensembl",
         "genes",
         "info_from_registry",
         "main_settings.json"
     )
    with open(settings, "r") as f:
        return json.load(f)

def get_server_settings_anno(settings: dict) -> dict:
    """
    Determine and return server connection settings for pipeline and core databases.

    Uses either custom server settings from the configuration or falls back
    to environment-variable-based defaults depending on the 'server_set' parameter.

    Args:
        settings (dict): The pipeline settings dictionary.

    Returns:
        dict: Nested dictionary with keys 'pipeline_db' and 'core_db', each containing
              connection parameters like 'db_host', 'db_user', 'db_port', and 'db_password'.

    Raises:
        ValueError: If the 'server_set' value is unknown or unsupported.
    """

    logger.info("Getting server settings")
    custom = settings.get("custom_server", {})
    # Check if all custom_server values are non-empty
    if all(custom.get(k) for k in ["pipeline_db_host", "pipeline_db_port", "core_db_host", "core_db_port"]):
        logger.info("Custom server settings detected")
        return {
            "pipeline_db": {
                "db_host": custom["pipeline_db_host"],
                "db_user": settings["user"],
                "db_port": custom["pipeline_db_port"],
                "db_password": settings["password"],
            },
            "core_db": {
                "db_host": custom["core_db_host"],
                "db_user": settings["user"],
                "db_port": custom["core_db_port"],
                "db_password": settings["password"]}
        }
    else:
        # Fallback based on server_set
        server_set = str(settings.get("server_set", "1"))  # default to "1" if missing

        if server_set == "1":
            logger.info("Server set 1 detected")
            return {
                "pipeline_db": {
                "db_host": os.environ.get("GBS4"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP4")),
                "db_password": settings["password"],
            },
            "core_db": {
                "db_host": os.environ.get("GBS3"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP3")),
                "db_password": settings["password"]}
            }

        elif server_set == "2":
            logger.info("Server set 2 detected")
            return {
                "pipeline_db": {
                "db_host": os.environ.get("GBS7"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP7")),
                "db_password": settings["password"],
            },
            "core_db": {
                "db_host": os.environ.get("GBS6"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP6")),
                "db_password": settings["password"],}
            }
        else:
            raise ValueError(f"Unknown server_set value: {server_set}")

def get_server_settings_main(settings: dict) -> dict:
    """
    Determine and return server connection settings for pipeline and core databases.

    Uses either custom server settings from the configuration or falls back
    to environment-variable-based defaults depending on the 'server_set' parameter.

    Args:
        settings (dict): The pipeline settings dictionary.

    Returns:
        dict: Nested dictionary with keys 'pipeline_db' and 'core_db', each containing
              connection parameters like 'db_host', 'db_user', 'db_port', and 'db_password'.

    Raises:
        ValueError: If the 'server_set' value is unknown or unsupported.
    """

    logger.info("Getting server settings")
    custom = settings.get("custom_server", {})
    # Check if all custom_server values are non-empty
    if all(custom.get(k) for k in ["pipeline_db_host", "pipeline_db_port", "core_db_host", "core_db_port", "databases_host", "databases_port"]):
        logger.info("Custom server settings detected")
        return {
            "pipeline_db": {
                "db_host": custom["pipeline_db_host"],
                "db_user": settings["user"],
                "db_port": custom["pipeline_db_port"],
                "db_password": settings["password"],
            },
            "core_db": {
                "db_host": custom["core_db_host"],
                "db_user": settings["user"],
                "db_port": custom["core_db_port"],
                "db_password": settings["password"]},
            "databases": {
                "db_host": custom["databases_host"],
                "db_user": settings["user"],
                "db_port": custom["databases_port"],
                "db_password": settings["password"]}
        }
    else:
        # Fallback based on server_set
        server_set = str(settings.get("server_set", "1"))  # default to "1" if missing

        if server_set == "1":
            logger.info("Server set 1 detected")
            return {
                "pipeline_db": {
                "db_host": os.environ.get("GBS4"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP4")),
                "db_password": settings["password"],
            },
            "core_db": {
                "db_host": os.environ.get("GBS2"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP2")),
                "db_password": settings["password"]},
            "databases": {
                "db_host": os.environ.get("GBS3"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP3")),
                "db_password": settings["password"]
            }
            }


        elif server_set == "2":
            logger.info("Server set 2 detected")
            return {
                "pipeline_db": {
                "db_host": os.environ.get("GBS7"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP7")),
                "db_password": settings["password"],
            },
            "core_db": {
                "db_host": os.environ.get("GBS5"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP5")),
                "db_password": settings["password"]},
            "databases": {
                "db_host": os.environ.get("GBS6"),
                "db_user": settings["user"],
                "db_port": int(os.environ.get("GBP6")),
                "db_password": settings["password"]
            }
            }
        else:
            raise ValueError(f"Unknown server_set value: {server_set}")



def get_metadata_from_registry(server_info: dict, assembly_accession, settings: dict) -> Optional[Dict[str, Any]]:
    """
    Retrieve registry metadata for a given genome assembly accession or list of accessions.

    Supports two modes:
    - If 'init_file' is specified in settings, loads registry info via custom INI file parsing.
    - Otherwise, queries the MySQL database for metadata matching the accession(s).

    Args:
        server_info (dict): MySQL server connection info with 'registry' key.
        assembly_accession (str or list): Single GCA accession or list of GCAs.
        settings (dict): Additional parameters, may include 'init_file'.

    Returns:
        dict or None: Assembly metadata dictionary if found, None if no record matches,
                      or empty dict if an error occurs.

    Raises:
        FileNotFoundError: If 'init_file' is specified but not found.
    """

    # Check if using custom init file
    if "init_file" in settings and settings["init_file"]:
        logger.info("Initialization file detected")
        init_file = Path(settings["init_file"])
        if init_file.exists():
            logger.info("Loading info from init file")
            registry_info = custom_loading(settings)
            registry_info = assign_clade_info_custom_loading(registry_info)
            return registry_info
        raise FileNotFoundError(f"INI file not found: {init_file}")

    try:
        if isinstance(assembly_accession, str):
            # single accession → wrap in list
            assembly_accessions = [assembly_accession]
        else:
            assembly_accessions = assembly_accession

        # Build placeholders for SQL query
        placeholders = ",".join(["%s"] * len(assembly_accessions))

        registry_query = f"""
            SELECT 
                s.species_taxon_id, 
                a.lowest_taxon_id AS taxon_id,
                a.asm_name AS assembly_name,
                s.common_name, 
                a.refseq_accession AS assembly_refseq_accession,
                a.release_date AS assembly_date,
                s.scientific_name AS species_name,
                a.assembly_id, 
                mb.bioproject_name AS assembly_group
            FROM assembly a
            JOIN bioproject b ON a.assembly_id = b.assembly_id
            JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
            LEFT JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
            WHERE CONCAT(a.gca_chain, '.', a.gca_version) IN ({placeholders})
        """

        registry_info = mysql_fetch_data(
            registry_query,
            host=server_info["registry"]["db_host"],
            user=server_info["registry"]["db_user"],
            port=server_info["registry"]["db_port"],
            database=server_info["registry"]["db_name"],
            password= "",
            params=assembly_accessions,
        )
        if not registry_info:
            raise ValueError(f"No registry data found for accessions: {assembly_accessions}")

        logger.info(f"Registry query successful for {assembly_accessions}")
        return registry_info[0]

    except pymysql.Error as err:
        logger.error("MySQL error: %s", err)
        return {}

def add_generated_data(server_info: dict, assembly_accession: str, settings: dict) -> dict:
    """
    Enrich registry metadata with clade assignments and derived pipeline variables.

    Args:
        server_info (dict): Server connection details.
        assembly_accession (str): Genome assembly accession.
        settings (dict): Pipeline settings.

    Returns:
        dict: Enriched metadata dictionary with added fields for pipeline usage.
    """

    registry_info = get_metadata_from_registry(server_info, assembly_accession, settings)
    logger.info(f"Data collected from registry {assembly_accession}")
    clade, genus_id, clade_metadata = assign_clade(server_info,registry_info)
    registry_info["clade"] = clade
    registry_info["genus_taxon_id"] = genus_id

    if clade_metadata:
        registry_info.update(clade_metadata)

    info_dict = registry_info
    # Create variables for pipeline
    info_dict["strain_type"] = "strain"

    if "alternate_haplotype" in registry_info.get("assembly_name", ""):
        info_dict["common_name"] = "alternate haplotype"
        info_dict["species_strain"] = "alternate haplotype"

    if not registry_info.get("common_name"):
        info_dict["common_name"] = "NA"

    info_dict["species_display_name"] = (
        f"{registry_info['species_name']} ({registry_info['common_name']}) - {assembly_accession}"
    )
    info_dict["species_strain"] = "reference"

    raw_species = (
        registry_info.get("species_name", "").strip().lower().replace(" ", "_")
    )
    species_name = raw_species.rstrip("_")  # Remove trailing underscore

    # Extract binomial name
    parts = species_name.split("_")
    if len(parts) >= 2:
        p1, p2 = parts[:2]
        max_len = 15
        production_name = f"{p1[:max_len]}_{p2[:max_len]}"
    else:
        production_name = ""

    production_gca = assembly_accession.replace('.', 'v').replace('_', '').lower()
    production_name += f"_{production_gca}"

    # Update dictionary
    info_dict["species_name"] = species_name
    info_dict["production_name"] = production_name
    info_dict["species_strain_group"] = production_name
    info_dict["species_url"] = (
        f"{registry_info['species_name'].capitalize()}_{assembly_accession}"
    )
    info_dict["core_dbname"] = f"{settings['dbowner']}_{production_gca}_core_{settings['release_number']}_1"

    logger.info(f"Values formatted for {assembly_accession}")

    return info_dict

def get_rna_and_busco_check_threshold(anno_settings: dict) -> dict:
    """
    Extract thresholds for RNA-seq and BUSCO checks from the pipeline settings.

    Args:
        anno_settings (dict): Pipeline settings dictionary.

    Returns:
        dict: Dictionary with keys 'busco_threshold', 'busco_lower_threshold',
              'busco_difference_threshold', 'rnaseq_main_file_min_lines',
              and 'rnaseq_genus_file_min_lines'.
    """

    return {
    "busco_threshold": anno_settings["busco_threshold"],
    "busco_lower_threshold": anno_settings["busco_lower_threshold"],
    "busco_difference_threshold": anno_settings["busco_difference_threshold"],
    "rnaseq_main_file_min_lines": anno_settings["rnaseq_main_file_min_lines"],
    "rnaseq_genus_file_min_lines": anno_settings["rnaseq_genus_file_min_lines"],
}


def get_info_for_pipeline_anno(settings: dict, info_dict: dict, assembly_accession: str, anno_settings: dict) -> dict:
    """
    Prepare and add file paths and pipeline-specific parameters needed for running annotation.

    Args:
        settings (dict): Pipeline settings.
        info_dict (dict): Metadata dictionary for the assembly.
        assembly_accession (str): Assembly accession string.
        anno_settings (dict): Annotation-specific settings.

    Returns:
        dict: Updated info_dict with added file paths and pipeline parameters.
    """

    logger.info("Getting info for pipeline settings for %s", assembly_accession)
    output_path = Path(settings["base_output_dir"]) / assembly_accession
    genome_files_dir = output_path / "genome_files"
    toplevel_genome_file = (
        output_path / f"{info_dict['species_name']}_toplevel.fa"
    )

    short_read_dir = output_path / "short_read_fastq"
    if "use_existing_short_read_dir" in anno_settings and os.path.isdir(
        anno_settings["use_existing_short_read_dir"]
    ):
        short_read_dir = anno_settings["use_existing_short_read_dir"]

    long_read_dir = output_path / "long_read_fastq"
    gst_dir = output_path / "gst"
    rnaseq_summary_file = (
        short_read_dir / f"{info_dict['production_name']}.csv"
    )
    rnaseq_summary_file_genus = (
        short_read_dir / f"{info_dict['production_name']}_gen.csv"
    )
    long_read_summary_file = (
        short_read_dir / f"{info_dict['production_name']}_long_read.csv"
    )
    reheadered_toplevel_genome_file = (
            output_path
            / f"{info_dict['species_name']}_reheadered_toplevel.fa"
    )
    diamond_validation_db = Path(anno_settings["diamond_validation_db"])
    current_genebuild = settings["current_genebuild"]
    num_threads = anno_settings["num_threads"]
    ensembl_release = settings["release_number"]


    # Add values back to dictionary
    info_dict.update(
        {
            "output_path": str(output_path),
            "genome_files_dir": str(genome_files_dir),
            "toplevel_genome_file": str(toplevel_genome_file),
            "reheadered_toplevel_genome_file": str(reheadered_toplevel_genome_file),
            "short_read_dir": str(short_read_dir),
            "long_read_dir": str(long_read_dir),
            "gst_dir": str(gst_dir),
            "diamond_validation_db": Path(diamond_validation_db),
            "rnaseq_summary_file": str(rnaseq_summary_file),
            "rnaseq_summary_file_genus": str(rnaseq_summary_file_genus),
            "long_read_summary_file": str(long_read_summary_file),
            "current_genebuild": current_genebuild,
            "assembly_accession": str(assembly_accession),
            "num_threads":str(num_threads),
            "ensembl_release": ensembl_release
        }
    )

    return info_dict

def get_info_for_pipeline_main(settings: dict, info_dict: dict, assembly_accession: str) -> dict:
    """
    Prepare and add file paths and pipeline-specific parameters needed for running annotation.

    Args:
        settings (dict): Pipeline settings.
        info_dict (dict): Metadata dictionary for the assembly.
        assembly_accession (str): Assembly accession string.

    Returns:
        dict: Updated info_dict with added file paths and pipeline parameters.
    """

    logger.info("Getting info for pipeline settings for GCA %s", assembly_accession)
    output_path = Path(settings["base_output_dir"]) / info_dict["species_name"] / assembly_accession
    rnaseq_dir = output_path / "rnaseq"
    long_read_dir = output_path / "long_read"
    gst_dir = output_path / "gst"
    rnaseq_summary_file = rnaseq_dir / f"{info_dict['species_name']}.csv"
    rnaseq_summary_file_genus = rnaseq_dir / f"{info_dict['species_name']}_gen.csv"
    long_read_fastq_dir = long_read_dir / "input"
    current_genebuild = settings["current_genebuild"]
    registry_file = output_path / "Databases.pm"
    release_number = settings["release_number"]
    dbname_accession = assembly_accession.replace('.', 'v').replace('_', '').lower()
    email_address = settings["email"]
    dbowner = settings["dbowner"]
    user_r = settings["user_r"]
    user = settings["user"]
    password = settings["password"]
    server_set = settings["server_set"]
    pipeline_name = settings["pipeline_name"]
    long_read_summary_file = long_read_dir / f"{info_dict['species_name']}_long_read.csv"
    long_read_summary_file_genus = long_read_dir / f"{info_dict['species_name']}_long_read_gen.csv"
    pipe_db_name = f"{dbowner}_{dbname_accession}_pipe_{release_number}"



    # Add values back to dictionary
    info_dict.update(
        {
            "output_path": str(output_path),
            "rnaseq_dir": str(rnaseq_dir),
            "long_read_fastq_dir": str(long_read_fastq_dir),
            "registry_file": str(registry_file),
            "long_read_dir": str(long_read_dir),
            "gst_dir": str(gst_dir),
            "rnaseq_summary_file": str(rnaseq_summary_file),
            "rnaseq_summary_file_genus": str(rnaseq_summary_file_genus),
            "current_genebuild": current_genebuild,
            "assembly_accession": str(assembly_accession),
            "release_number": release_number,
            "dbname_accession": str(dbname_accession),
            "email_address":str(email_address),
            "dbowner": str(dbowner),
            "user_r": str(user_r),
            "user": str(user),
            "password": str(password),
            "server_set": server_set,
            "pipeline_name": str(pipeline_name),
            "long_read_summary_file": str(long_read_summary_file),
            "long_read_summary_file_genus": str(long_read_summary_file_genus),
            "pipe_db_name": str(pipe_db_name)
        }
    )

    return info_dict

def copy_general_module():
    logger.info("Copying general module")
    enscode = os.environ.get("ENSCODE")
    if not enscode:
        raise OSError("Environment variable 'ENSCODE' is not set")

    analysis_path = os.path.join(enscode, "ensembl-analysis")
    general_file = Path(analysis_path) / "modules" / "Bio" / "EnsEMBL" / "Analysis" / "Config" / "General.pm"
    example_file = general_file.with_suffix(".pm.example")

    # Copy if missing
    if not general_file.exists():
        if example_file.exists():
            shutil.copy(example_file, general_file)
            logger.info(f"Copied {example_file} → {general_file}")
        else:
            raise FileNotFoundError(f"Missing example file: {example_file}")
    else:
        logger.info(f"{general_file} already exists, nothing to do.")

def current_projection_source_db(projection_source_production_name, user):
    """
    Find the latest core database for the given projection_source_production_name on mysql-ens-mirror-1.
    If not found, fall back to homo_sapiens core.
    Args:
        projection_source_production_name (str): Name to look for in core name (comes from clade_settings).
        user (str): Read-only username.
    Returns:
        str: Projection source db name.
    Raises:
        RuntimeError: If database could not be found.

    """
    conn = pymysql.connect(
        host="mysql-ens-mirror-1",
        user=user,
        password="",
        port=4240
    )
    try:
        with conn.cursor() as cur:
            # First try the species-specific database
            cur.execute("SHOW DATABASES LIKE %s", (f"{projection_source_production_name}_core%",))
            out = [row[0] for row in cur.fetchall()]

            # If none, fall back to human core
            if not out:
                logger.info(f"No core DB for {projection_source_production_name}, falling back to Homo sapiens")
                cur.execute("SHOW DATABASES LIKE 'homo_sapiens_core%'")
                out = [row[0] for row in cur.fetchall()]

        if out:
            return out[-1]  # return the last (most recent) DB
        else:
            raise RuntimeError("No suitable projection database found.")
    finally:
        conn.close()


def custom_loading(settings):
    """
    Load custom key-value pairs from an INI-style initialization file specified in settings.

    Reads the file line-by-line, extracting key=value pairs into a dictionary.
    Lines that do not conform to key=value format or are blank are skipped with a message.

    Args:
        settings (dict): Pipeline settings dictionary that must include
                         the key 'init_file' pointing to the INI file path.

    Returns:
        dict: Dictionary containing key-value pairs loaded from the INI file.

    Raises:
        Exception: If the file cannot be opened or read.
    """
    logger.info("Custom loading started")
    init_file = Path(settings["init_file"])

    # Initialize variables
    custom_dict = {}
    custom_loading = False

    if init_file.exists():
        try:
            with init_file.open("r") as f:
                print("Using custom loading .ini file.")
                custom_loading = True

                for line in f:
                    line = line.strip()

                    if "=" in line:
                        key, value = map(str.strip, line.split("=", 1))
                        print(f"Found key/value pair: {key} => {value}")
                        custom_dict[key] = value
                    elif line == "":
                        continue
                    else:
                        print(f"Line format not recognised. Skipping line:\n{line}")

        except OSError as e:
            raise Exception(f"Could not open or read {init_file}") from e

    return custom_dict


def create_dir(path, mode=None):
    """
    Create a directory and optionally set its permissions.

    This function creates the directory at 'path' including any necessary
    parent directories. If 'mode' is provided, it changes the directory's permissions.

    Args:
        path (str or Path): The directory path to create.
        mode (int, optional): File system permissions mode (e.g., 0o775). Defaults to None.

    Raises:
        RuntimeError: If the directory creation or permission change fails.
    """

    try:
        os.makedirs(path, exist_ok=True)
        if mode is not None:
            os.chmod(path, mode)
    except Exception as e:
        raise RuntimeError(f"Failed to create dir: {path}") from e

def main(gcas, settings_file):
    """
    Create parameters for annotation pipeline.

    Args:
        gcas (str): Path to the list of GCAs
        settings_file : Pipeline settings JSON

    Returns:
        all_output_params (dict): Dictionary of all GCA output parameters
        saved_paths (dict): Paths to saved JSONs, with separate entries for 'anno' and 'main'
    """
    import json, os
    from pathlib import Path

    # Load settings
    settings = load_settings(settings_file)

    # Read in GCAs from file
    with open(gcas, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    gca_dict = {gca: {} for gca in lines}
    logger.info(f"Found {len(gca_dict)} GCAs")

    # Check if init_file exists
    has_init_file = settings.get("init_file") and os.path.isfile(settings["init_file"])

    all_output_params = {}

    # Loop through GCAs
    for gca in gca_dict:
        gca_dict[gca]["assembly_accession"] = gca

        # Registry info
        server_info = {
            "registry": {
                "db_host": os.environ.get("GBS1"),
                "db_user": settings["user_r"],
                "db_user_w": settings["user"],
                "db_port": int(os.environ.get("GBP1")),
                "db_name": "gb_metadata_start_test",
                "password": settings["password"]
            }
        }
        logger.info(server_info["registry"])

        # Check if it has annotation history
        if not has_init_file and int(settings.get("current_genebuild", 0)) == 0:
            check_if_annotated(gca, server_info)

        # Create output params
        info_dict = add_generated_data(server_info, gca, settings)

        if 'repbase_library' in info_dict and isinstance(info_dict['repbase_library'], str):
            # Vertebrate (main) pipeline
            pipeline_type = "main"

            server_settings = get_server_settings_main(settings)
            server_info.update(server_settings)
            main_settings = load_main_settings()

            info_dict["databases_host"] = server_info["databases"]["db_host"]
            info_dict["databases_port"] = server_info["databases"]["db_port"]
            info_dict["registry_host"] = server_info["registry"]["db_host"]
            info_dict["registry_port"] = server_info["registry"]["db_port"]
            info_dict["registry_db"] = server_info["registry"]["db_name"]
            info_dict["pipe_db_host"] = server_info["pipeline_db"]["db_host"]
            info_dict["pipe_db_port"] = server_info["pipeline_db"]["db_port"]
            info_dict["dna_db_host"] = server_info["core_db"]["db_host"]
            info_dict["dna_db_port"] = server_info["core_db"]["db_port"]


            if main_settings.get("replace_repbase_with_red_to_mask") == 1:
                info_dict["first_choice_repeat"] = "repeatdetector"

            # Protein BLAST DB based on clade
            clade = info_dict.get("clade", "").lower()
            if clade in ["mammalia", "rodentia", "primates", "marsupials"]:
                info_dict["protein_blast_db_file"] = "uniprot_mammalia_sp"
            elif clade in ["teleostei", "sharks"]:
                info_dict["protein_blast_db_file"] = "uniprot_vertebrataSP_plus_fishTR"
            else:
                info_dict["protein_blast_db_file"] = "uniprot_vertebrata_sp"

            # RepeatModeler library
            parent_name = get_parent_taxon(server_info, info_dict["species_taxon_id"])
            parent_name = parent_name.lower().replace(" ", "_")
            info_dict["repeatmodeler_library"] = os.path.join(
                os.environ["REPEATMODELER_DIR"], "species", parent_name, f"{parent_name}.repeatmodeler.fa"
            )

            output_params = get_info_for_pipeline_main(settings, info_dict, gca)
            create_dir(output_params["output_path"])
            edit_config_main(main_settings, output_params, pipeline_type)
            output_params["projection_source_db_name"] = current_projection_source_db(
                output_params["projection_source_production_name"], settings["user_r"]
            )
            output_params["stable_id_prefix"] = get_species_prefix(output_params["taxon_id"], server_info)
            output_params["stable_id_start"] = get_stable_space(
                output_params["taxon_id"], gca, output_params["assembly_id"], server_info
            )
            edit_config_main(main_settings, output_params, pipeline_type)
            copy_general_module()

        else:
            # Non-vertebrate (anno) pipeline
            pipeline_type = "anno"
            gca_dict[gca]["pipe_db_name"] = f"{settings['dbowner']}_{settings['pipeline_name']}_pipe_{settings['release_number']}"
            server_settings = get_server_settings_anno(settings)
            server_info.update(server_settings)

            server_info.setdefault("pipeline_db", {})["db_name"] = gca_dict[gca]["pipe_db_name"]
            server_info.setdefault("core_db", {})["db_name"] = info_dict["core_dbname"]
            info_dict["core_db"] = server_info["core_db"]
            info_dict["registry_db"] = server_info["registry"]
            
            # Assign BUSCO lineage
            busco_lineage_file = os.path.join(
            os.environ.get("ENSCODE"),
            "ensembl-genes",
            "src",
            "python",
            "ensembl",
            "genes",
            "metrics",
            "busco_lineage.json")

            with open(busco_lineage_file, "r") as f:
                dataset = json.load(f)

            ncbi_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{info_dict['taxon_id']}/dataset_report"

            busco_group_find = get_dataset_match(ncbi_url, dataset)

            if busco_group_find is not None:
                logger.info(f"Closest BUSCO lineage identified as {busco_group_find} for taxon ID {info_dict['taxon_id']}")
                info_dict["busco_group"] = busco_group_find
            else:
                logger.info(f"Falling back on BUSCO lineage from clade settings {info_dict['busco_group']}")

            # Load anno settings
            anno_settings = load_anno_settings()

            # DB adaptors
            core_adaptor = {
                "host": server_info["core_db"]["db_host"],
                "port": server_info["core_db"]["db_port"],
                "dbname": info_dict["core_dbname"],
                "user": settings["user_r"],
                "pass": settings["password"],
                "species": info_dict["production_name"],
                "group": "core",
            }
            registry_path = create_registry_entry(settings, server_info, core_adaptor)

            output_params = get_info_for_pipeline_anno(settings, info_dict, gca, anno_settings)
            output_params["registry_file"] = Path(registry_path)

            build_annotation_commands(core_adaptor, output_params, anno_settings, settings)

            rna_busco_settings = get_rna_and_busco_check_threshold(anno_settings)
            output_params.update(rna_busco_settings)
            output_params["stable_id_prefix"] = get_species_prefix(output_params["taxon_id"], server_info)
            output_params["stable_id_start"] = get_stable_space(
                output_params["taxon_id"], gca, output_params["assembly_id"], server_info
            )

        # Store output params
        output_params["pipeline"] = pipeline_type
        all_output_params[gca] = output_params
        logger.info(f"Finished with {gca}")

    # -------------------------
    # Save JSONs
    # -------------------------
    saved_paths = {}

    # Save anno GCAs as a single file
    anno_params = {g: p for g, p in all_output_params.items() if p["pipeline"] == "anno"}

    #Edit anno config
    if anno_params:
        # Call edit_config_anno ONCE here with only the generic settings
        first_gca_params = next(iter(anno_params.values()))
        edit_config_anno(
            anno_settings, 
            settings, 
            first_gca_params,   # just need a representative dict with registry_file, db names, etc.
            "anno", 
            server_settings
        )

        anno_json_path = Path(settings["base_output_dir"]) / "non_vert_pipeline_params.json"
        anno_json_path.parent.mkdir(parents=True, exist_ok=True)
        with anno_json_path.open("w") as f:
            json.dump(anno_params, f, indent=2, default=str)
        logger.info(f"Saved anno parameters to: {anno_json_path}")
        saved_paths["anno"] = anno_json_path

    # Save main GCAs as separate files
    main_params = {g: p for g, p in all_output_params.items() if p["pipeline"] == "main"}
    main_json_paths = {}
    for gca, params in main_params.items():
        path = Path(params["output_path"]) / "main_pipeline_params.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            json.dump({gca: params}, f, indent=2, default=str)
        logger.info(f"Saved main parameters for {gca} to: {path}")
        main_json_paths[gca] = path

    if main_json_paths:
        saved_paths["main"] = main_json_paths

    # Log if mixed pipelines are present
    pipelines_present = set(p["pipeline"] for p in all_output_params.values())
    if "anno" in pipelines_present and "main" in pipelines_present:
        msg = "ATTENTION: Both 'anno' and 'main' GCAs detected in input!"
        logger.info(msg)
        print(msg)

    logger.info("DONE")
    return server_info, all_output_params, saved_paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract metadata and initialize annotation pipelines."
    )
    parser.add_argument(
        "--gcas",
        type=str,
        required=True,
        help="Path to file containing GCA accessions (one per line)."
    )

    parser.add_argument(
        "--settings_file",
        type=str,
        required=True,
        help="Path to file containing edited settings."
    )

    args = parser.parse_args()

    server_info, all_output_params, saved_paths = main(args.gcas, args.settings_file)
    print("\n=== NEXT STEP INITIALISE PIPELINE ===")
    print("\n=== FINALLY RUN SEED IF NONVERT ===")
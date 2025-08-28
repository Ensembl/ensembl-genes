"""
This script connects to a genomic assembly registry, extracts assembly metadata,
and initializes annotation pipelines based on predefined configurations.

It supports both database-driven metadata retrieval and custom INI file loading,
handles directory setup, and creates necessary configuration and command files
for running genome annotation pipelines.

Typical usage:
    python start_pipeline_from_registry.py --gcas gca_list.txt --pipeline anno --settings_file settings.json
"""

import os
from pathlib import Path
import json
import pymysql # type: ignore
import argparse
import logging
from build_anno_commands import (build_annotation_commands)
from check_if_annotated import (check_if_annotated)
from mysql_helper import mysql_fetch_data
from assign_clade_based_on_tax import (assign_clade,assign_clade_info_custom_loading)
from create_pipe_reg import create_registry_entry
from create_config import edit_config
from assign_species_prefix import get_species_prefix # type: ignore
from assign_stable_space import get_stable_space

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
    settings = "src/python/ensembl/genes/info_from_registry/anno_settings.json"
    # settings = os.path.join(
    #     os.environ.get("ENSCODE"),
    #     "ensembl-genes",
    #     "src",
    #     "python",
    #     "ensembl",
    #     "genes",
    #     "info_from_registry",
    #     "anno_settings.json"
    # )
    with open(settings, "r") as f:
        return json.load(f)

def get_server_settings(settings: dict) -> dict:
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
    # Check if all custom_server values are non-empty (you can change logic if needed)
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



def get_metadata_from_registry(server_info: dict, assembly_accession, settings: dict) -> dict | None:
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
            JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
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
        if registry_info:
            return registry_info[0]

        return None

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

def get_rna_and_busco_check_threshold(settings: dict) -> dict:
    """
    Extract thresholds for RNA-seq and BUSCO checks from the pipeline settings.

    Args:
        settings (dict): Pipeline settings dictionary.

    Returns:
        dict: Dictionary with keys 'busco_threshold', 'busco_lower_threshold',
              'busco_difference_threshold', 'rnaseq_main_file_min_lines',
              and 'rnaseq_genus_file_min_lines'.
    """

    return {
    "busco_threshold": settings["busco_threshold"],
    "busco_lower_threshold": settings["busco_lower_threshold"],
    "busco_difference_threshold": settings["busco_difference_threshold"],
    "rnaseq_main_file_min_lines": settings["rnaseq_main_file_min_lines"],
    "rnaseq_genus_file_min_lines": settings["rnaseq_genus_file_min_lines"],
}


def get_info_for_pipeline(settings: dict, info_dict: dict, assembly_accession: str, anno_settings: dict) -> dict:
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

    logger.info("Getting info for pipeline settings for GCA %s", assembly_accession)
    output_path = Path(settings["base_output_dir"]) / assembly_accession
    genome_files_dir = output_path / "genome_files"
    toplevel_genome_file = (
        output_path / f"{info_dict['species_name']}_toplevel.fa"
    )

    short_read_dir = output_path / "short_read_fastq"
    if "use_existing_short_read_dir" in settings and os.path.isdir(
        settings["use_existing_short_read_dir"]
    ):
        short_read_dir = settings["use_existing_short_read_dir"]

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

def main(gcas, pipeline, settings_file):
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
    # Load settings
    settings = load_settings(settings_file)

    # Read in GCAs from file
    with open(gcas, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Create a dictionary with each GCA as a top-level key
    gca_dict = {gca: {} for gca in lines}
    logger.info(f"Found {len(gca_dict)} GCAs")

    # Check if init_file exists
    has_init_file = settings.get("init_file") and os.path.isfile(settings["init_file"])

    all_output_params = {}

    for gca in gca_dict:
        # Fill in basic info
        gca_dict[gca]["assembly_accession"] = gca
        gca_dict[gca]["pipe_db_name"] = f"{settings['dbowner']}_{settings['pipeline_name']}_pipe_{settings['release_number']}"

        # Get server settings
        server_info = get_server_settings(settings)

        # Add registry info
        server_info["registry"] = {
            "db_host": os.environ.get("GBS1"),
            "db_user": settings["user_r"],
            "db_user_w": settings["user"],
            "db_port": int(os.environ.get("GBP1")),
            "db_name": "gb_assembly_metadata_status_test",
            "password": settings["password"]
        }
        print(server_info["registry"])

        # Assign database names
        server_info.setdefault("pipeline_db", {})["db_name"] = gca_dict[gca]["pipe_db_name"]

        # Check if GCA is annotated
        logger.info(f"Checking {gca_dict[gca]['assembly_accession']} annotation status")
        if not has_init_file:
            if int(settings["current_genebuild"]) == 0:
                check_if_annotated(gca_dict[gca]["assembly_accession"], server_info)
            else:
                logger.info(f"Skipping annotation check {gca_dict[gca]['assembly_accession']}")

        # Create output params
        info_dict = add_generated_data(server_info, gca_dict[gca]["assembly_accession"], settings)
        server_info.setdefault("core_db", {})["db_name"] = info_dict["core_dbname"]
        info_dict["core_db"] = server_info["core_db"]
        info_dict["registry_db"] = server_info["registry"]

        if pipeline == "anno":
            logger.info("Anno setting detected")
            logger.info("Copy and modify pipeline config")
            #edit_config(settings)

            anno_settings = load_anno_settings()
            output_params = get_info_for_pipeline(settings, info_dict, gca_dict[gca]["assembly_accession"], anno_settings)

            # Create directories
            dirs_to_create = []
            # Create output dir with 775 permissions this was for BRAKER
            #create_dir(info_dict["output_path"], mode=0o775)

            # Add other dirs to the list
            dirs_to_create.extend(
                [
                    info_dict["genome_files_dir"],
                    info_dict["short_read_dir"],
                    info_dict["long_read_dir"],
                    info_dict["gst_dir"],
                ]
            )

            # Create the other directories without changing mode (default permissions)
            #for dir_path in dirs_to_create:
            #    create_dir(dir_path)
            logger.info("Created directories")

            # Add db adaptors to pipeline registry
            core_adaptor = {
                "host": server_info["core_db"]["db_host"],
                "port": server_info["core_db"]["db_port"],
                "dbname": info_dict["core_dbname"],
                "user": settings["user_r"],
                "pass": settings["password"],
                "species": info_dict["production_name"],
                "group": "core",
            }

            # Create a local copy of the registry and update the pipeline's resources with the new path
            #registry_path = create_registry_entry(settings, server_info, core_adaptor)
            #output_params["registry_file"] = Path(registry_path)
            # Build anno commands and add to dictionary
            build_annotation_commands(core_adaptor, output_params, anno_settings, settings)
            logger.info(f"Created anno commands for {gca_dict[gca]['assembly_accession']}")

            logger.info("Getting RNA and BUSCO thresholds")
            rna_busco_settings = get_rna_and_busco_check_threshold(settings)
            output_params.update(rna_busco_settings)

            #Assign species prefix
            output_params["stable_id_prefix"] = get_species_prefix(output_params["taxon_id"], server_info)
            logger.info(f"✅ Assigned prefix {output_params['stable_id_prefix']} for taxon {output_params.get('taxon_id')}")

            #Assign stable ID
            output_params["stable_id_start"] = get_stable_space(output_params["taxon_id"], gca_dict[gca]['assembly_accession'], output_params["assembly_id"], server_info)

            # Store the output_params for this GCA
            all_output_params[gca] = output_params
            logger.info(f"Finished with {gca_dict[gca]['assembly_accession']}")

    # Save all_output_params to output directory
    output_json_path = Path(settings["base_output_dir"]) / "non_vert_pipeline_params.json"

    try:
        with output_json_path.open("w") as f:
            json.dump(all_output_params, f, indent=2, default=str)
        logger.info("Saved all_output_params to: %s", output_json_path)
    except Exception as e:
        logger.error("Failed to save output dictionary: %s", e)

    logger.info("DONE")

    return all_output_params, output_json_path


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
        "--pipeline",
        type=str,
        choices=["anno", "main"],  # Add more valid pipeline types as needed
        required=True,
        help="Pipeline type to initialize (e.g., 'anno')."
    )

    parser.add_argument(
        "--settings_file",
        type=str,
        required=True,
        help="Path to file containing edited settings."
    )

    args = parser.parse_args()

    results, output_json_path = main(args.gcas, args.pipeline, args.settings_file)
    print("\n=== NEXT STEP INITIALISE PIPELINE ===")
    print("\n=== FINALLY RUN SEED NONVERT ===")
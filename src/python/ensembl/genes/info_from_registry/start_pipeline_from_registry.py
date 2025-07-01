"""
This script connects to the registry, extracts information,
and initializes annotation pipelines based on predefined configurations.
"""

import os
from pathlib import Path
import sys
import json
import pymysql

from src.python.ensembl.genes.info_from_registry.build_anno_commands import (
    build_annotation_commands,
)
from src.python.ensembl.genes.info_from_registry.check_if_annotated import (
    check_if_annotated,
)
from src.python.ensembl.genes.info_from_registry.create_pipe_reg import (
    create_registry_entry,
)
from src.python.ensembl.genes.info_from_registry.mysql_helper import mysql_fetch_data
from src.python.ensembl.genes.info_from_registry.assign_clade_based_on_tax import (
    assign_clade,
    assign_clade_info_custom_loading,
)


def get_metadata_from_registry(server_info, assembly_accession, param):
    """
    Retrieves registry metadata for a given genome assembly accession.

    This function supports two modes of operation:
    1. If a custom `init_file` is provided in `param`, it loads registry metadata using custom logic
       (e.g., from an INI file).
    2. Otherwise, it queries a MySQL registry database to fetch metadata related to the given
       `assembly_accession` or list of accessions.

    Args:
        server_info (dict): A dictionary containing MySQL server connection info under the
                            key 'registry', including 'db_host', 'db_user', and 'db_port'.
        assembly_accession (str or list): A single GCA accession (as a string) or a list of accessions
                                          to query from the registry.
        param (dict): A dictionary of additional parameters. If it includes a valid 'init_file' path,
                      custom loading is used instead of database queries.

    Returns:
        dict or None: A dictionary containing assembly metadata if found, or `None` if no matching
                      assembly is found. If an error occurs, an empty dictionary is returned.

    Raises:
        FileNotFoundError: If the provided 'init_file' path does not exist.
    """
    # Check if using custom init file
    if "init_file" in param and param["init_file"]:
        init_file = Path(param["init_file"])
        if init_file.exists():
            registry_info = custom_loading(param)
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
                mb.bioproject_name AS assembly_group,
                CONCAT(a.gca_chain, '.', a.gca_version) AS gca
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
            database="gb_assembly_metadata",
            params=assembly_accessions,
        )

        if registry_info:
            return registry_info[0]  # return single dict
        return None

    except pymysql.Error as err:
        print(f"MySQL error: {err}")
        return {}


def add_generated_data(server_info, assembly_accession, param):
    registry_info = get_metadata_from_registry(server_info, assembly_accession, param)

    clade, genus_id, clade_metadata = assign_clade(
        server_info,
        registry_info,
    )

    registry_info["clade"] = clade
    registry_info["genus_taxon_id"] = genus_id

    print("=== CLADE METADATA ===")
    print(json.dumps(clade_metadata, indent=2, default=str))
    print("======================")

    if clade_metadata:
        registry_info.update(clade_metadata)

    info_dict = registry_info
    print("=== INFO DICT ===")
    print(json.dumps(info_dict, indent=2, default=str))
    print("======================")
    # Create variables for pipeline
    info_dict["strain_type"] = "strain"

    if "alternate_haplotype" in registry_info.get("assembly_name", ""):
        info_dict["common_name"] = "alternate haplotype"
        info_dict["species_strain"] = "alternate haplotype"

    if not registry_info.get("common_name"):
        info_dict["common_name"] = "NA"

    info_dict["species_url"] = (
        f"{registry_info['species_name']}_{param['assembly_accession']}"
    )
    info_dict["species_display_name"] = (
        f"{registry_info['species_name']} ({registry_info['common_name']}) - {param['assembly_accession']}"
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
        binomial_species_name = f"{p1}_{p2}"
        max_len = 15
        production_name = f"{p1[:max_len]}_{p2[:max_len]}"
    else:
        binomial_species_name = ""
        production_name = ""

    production_gca = assembly_accession.replace('.', 'v').replace('_', '').lower()
    production_name += f"_{production_gca}"

    # Update dictionary
    info_dict["species_name"] = species_name
    info_dict["binomial_species_name"] = binomial_species_name
    info_dict["production_name"] = production_name
    info_dict["species_strain_group"] = production_name

    return info_dict


def get_info_from_params(param, info_dict):
    ensembl_release = param["ensembl_release"]
    output_path = Path(param["base_output_dir"]) / param["assembly_accession"]
    genome_files_dir = Path(param["output_path"]) / "genome_files"
    toplevel_genome_file = (
        Path(param["output_path"]) / f"{info_dict['species_name']}_toplevel.fa"
    )
    reheadered_toplevel_genome_file = (
        Path(param["output_path"])
        / f"{info_dict['species_name']}_reheadered_toplevel.fa"
    )
    short_read_dir = Path(param["output_path"]) / "short_read_fastq"
    if "use_existing_short_read_dir" in param and os.path.isdir(
        param["use_existing_short_read_dir"]
    ):
        short_read_dir = param["use_existing_short_read_dir"]
    long_read_dir = Path(param["output_path"]) / "long_read_fastq"
    gst_dir = Path(param["output_path"]) / "gst"
    rnaseq_summary_file = (
        Path(param["short_read_dir"]) / f"{info_dict['production_name']}.csv"
    )
    rnaseq_summary_file_genus = (
        Path(param["short_read_dir"]) / f"{info_dict['production_name']}_gen.csv"
    )
    long_read_summary_file = (
        Path(param["long_read_dir"]) / f"{info_dict['production_name']}_long_read.csv"
    )
    repeatmodeler_library = Path(param["repeatmodeler_library"])
    current_genebuild = param["current_genebuild"]
    assembly_accession = param["assembly_accession"]
    core_dbname = param["core_dbname"]
    num_threads = param["num_threads"]


    # Add values back to dictionary
    info_dict.update(
        {
            "ensembl_release": ensembl_release,
            "output_path": str(output_path),
            "genome_files_dir": str(genome_files_dir),
            "toplevel_genome_file": str(toplevel_genome_file),
            "reheadered_toplevel_genome_file": str(reheadered_toplevel_genome_file),
            "short_read_dir": str(short_read_dir),
            "long_read_dir": str(long_read_dir),
            "gst_dir": str(gst_dir),
            "rnaseq_summary_file": str(rnaseq_summary_file),
            "rnaseq_summary_file_genus": str(rnaseq_summary_file_genus),
            "long_read_summary_file": str(long_read_summary_file),
            "repeatmodeler_library": str(repeatmodeler_library),
            "current_genebuild": current_genebuild,
            "assembly_accession": str(assembly_accession),
            "core_dbname": str(core_dbname),
            "num_threads":str(num_threads),
        }
    )

    return info_dict


def custom_loading(param):

    init_file = Path(param["init_file"])

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
    try:
        os.makedirs(path, exist_ok=True)
        if mode is not None:
            os.chmod(path, mode)
    except Exception as e:
        raise RuntimeError(f"Failed to create dir: {path}") from e


def main():
    # Inherit params from perl
    param = json.load(sys.stdin)
    assembly_accession = param["assembly_accession"]
    # dna_db_name = f"{param['dbowner']}_{param['production_name']}{param['production_name_modifier']}_core_{param['release_number']}"
    param.update(
        {
            "pipe_db_name": f"{param['dbowner']}_{param['pipeline_name']}_pipe_{param['release_number']}",
            "core_dbname": f"{param['dbowner']}_{param['production_gca']}_core_{param['ensembl_release']}_1"        }
    )

    # Create server info dictionary
    server_info = {
        "registry": {
            "db_host": os.environ.get("GBS1"),
            "db_user": param["user_r"],
            "db_port": int(os.environ.get("GBP1")),
        },
        "pipeline_db": {
            "db_host": os.environ.get("GBS7"),
            "db_user": param["user"],
            "db_port": int(os.environ.get("GBP7")),
            "db_password": param["password"],
            "db_name": param["pipe_db_name"],
        },
        "core_db": {
            "db_host": os.environ.get("GBP6"),
            "db_user": param["user"],
            "db_port": int(os.environ.get("GBP6")),
            "db_password": param["password"],
            "db_name": param["core_dbname"],
        }
    }

    # Check if GCA annotation exists for GCA and stop the pipeline if it does (unless custom loading)
    if param.get("init_file") and os.path.isfile(param["init_file"]):
        check_if_annotated(assembly_accession, server_info)

    # Create output params
    info_dict = add_generated_data(server_info, assembly_accession, param)
    info_dict["core_db"] = server_info["core_db"]
    output_params = get_info_from_params(param, info_dict)


    # Create directories
    dirs_to_create = []
    # Create output dir with 775 permissions this was for BRAKER
    create_dir(info_dict["output_path"], mode=0o775)


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
    for dir_path in dirs_to_create:
        create_dir(dir_path)

    # Add db adaptors to pipeline registry
    adaptors = {
        "core_string": {
            "host": os.environ.get("GBS7"),
            "port": int(os.environ.get("GBP7")),
            "dbname": param["core_dbname"],
            "user": param["user_r"],
            "pass": param["password"],
            "species": info_dict["production_name"],
            "group": "core",
        },
    }

    # Create a local copy of the registry and update the pipeline's resources with the new path
    registry_path = create_registry_entry(param, server_info, adaptors)
    output_params["registry_file"] = str(registry_path)


    # Build anno commands and add to dictionary
    build_annotation_commands(adaptors, output_params)

    return output_params



if __name__ == "__main__":
    result = main()
    print(json.dumps(result))
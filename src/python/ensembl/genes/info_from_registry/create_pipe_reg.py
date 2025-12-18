"""
create_pipe_reg.py

This module updates the Ensembl pipeline registry configuration and inserts
the appropriate database adaptor entries for use in a genebuild pipeline.

It performs two key tasks:
1. Updates references to the shared registry file within the pipeline's resource
   description to use a copied, local version.
2. Appends new database connection entries into the copied registry file for use
   by Hive pipeline workers.

Functions:
----------
- update_registry_path_and_create_entry: Updates resource_description in pipeline DB.
- create_registry_entry: Adds new DBAdaptor entries into the local registry file.
"""

import shutil
from pathlib import Path
import os
import logging
from typing import Dict, Any


from ensembl.genes.info_from_registry.mysql_helper import mysql_update

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("pipeline_setup.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def update_registry_path_in_pipedb(
    parent_dir: str, server_info: Dict[str, Dict[str, Any]]
) -> None:
    """
    Updates the resource_description table in the pipeline database by replacing
    the shared registry path with a local copy. This ensures worker nodes access
    the correct local registry.

    Parameters:
    ----------
    settings : dict
        Must contain:
            - 'base_output_dir': Path to where the new registry will be copied.

    server_info : dict
        Connection info for the pipeline database under:
            - server_info["pipeline_db"]["db_host"]
            - server_info["pipeline_db"]["db_user"]
            - server_info["pipeline_db"]["db_port"]
            - server_info["pipeline_db"]["db_password"]
            - server_info["pipeline_db"]["db_name"]
    """
    logger.info("Updating registry path in pipeline DB")
    registry_file = os.path.join(
        str(os.environ.get("ENSCODE")),
        "ensembl-analysis",
        "scripts",
        "genebuild",
        "gbiab",
        "support_files",
        "Databases.pm",
    )
    logger.debug(  # pylint: disable=logging-fstring-interpolation
        f"parent_dir : {parent_dir}"
    )
    new_registry_file = Path(parent_dir).resolve() / "Databases.pm"

    if not new_registry_file.exists():
        raise OSError(f"Registry file doesn't exist: {new_registry_file}")

    logger.info(  # pylint: disable=logging-fstring-interpolation
        f"Registry file exist at {new_registry_file}"
    )

    update_resources_query = """
		UPDATE resource_description 
		SET worker_cmd_args = REPLACE(worker_cmd_args, %s, %s);
	"""

    success = mysql_update(
        query=update_resources_query,
        host=server_info["pipeline_db"]["db_host"],
        user=server_info["pipeline_db"]["db_user"],
        port=server_info["pipeline_db"]["db_port"],
        password=server_info["pipeline_db"]["db_password"],
        database=server_info["pipeline_db"]["db_name"],
        params=(registry_file, new_registry_file),
    )

    if success:
        logger.info("Registry path updated successfully.")
    else:
        raise RuntimeError("Failed to update registry path.")


def create_registry_entry(
    settings: Dict[str, Any],
    server_info: Dict[str, Dict[str, Any]],  # pylint: disable=unused-argument
    core_adaptor: Dict[str, str],
) -> Path:
    """
    Updates the local registry file with new DBAdaptor connection details for a
    core database, using the provided connection settings.

    Parameters:
    ----------
    settings : dict
        Must contain:
            - 'base_output_dir': str, where the local Databases.pm file is placed.

    server_info : dict
        Used to pass to `update_registry_path_and_create_entry`.

    core_adaptor : dict
        Must contain:
            - 'host', 'port', 'dbname', 'user', 'pass', 'species', 'group'

    Returns:
    -------
    Path to the modified local registry file (Databases.pm)

    Raises:
    ------
    FileNotFoundError:
        If the registry file does not exist after copying.

    RuntimeError:
        If there is an error overwriting the original registry file.
    """
    registry_path = Path(settings["base_output_dir"]) / "Databases.pm"
    registry_file = os.path.join(
        str(os.environ.get("ENSCODE")),
        "ensembl-analysis",
        "scripts",
        "genebuild",
        "gbiab",
        "support_files",
        "Databases.pm",
    )

    if not registry_path.exists():
        shutil.copy(registry_file, registry_path)

    if not registry_path.exists():
        raise FileNotFoundError(f"A registry file was not found at: {registry_path}")

    def db_string(db_details):
        return (
            f"Bio::EnsEMBL::DBSQL::DBAdaptor->new(\n"
            f"  -host => '{db_details['host']}',\n"
            f"  -port => '{db_details['port']}',\n"
            f"  -dbname => '{db_details['dbname']}',\n"
            f"  -user => '{db_details['user']}',\n"
            f"  -pass => '{db_details['pass']}',\n"
            f"  -species => '{db_details['species']}',\n"
            f"  -group => '{db_details['group']}',\n"
            f");\n"
        )

    # Generate registry entries
    core_string = db_string(core_adaptor)

    # Read the original registry content
    lines = registry_path.read_text().splitlines(keepends=True)

    # Insert database connection strings after first `{`
    new_lines = []
    inserted = False
    for line in lines:
        new_lines.append(line)
        if not inserted and "{" in line:
            new_lines.append(core_string)
            inserted = True

    # Write to a temporary file
    tmp_path = registry_path.with_suffix(".pm.tmp")
    tmp_path.write_text("".join(new_lines))

    # Overwrite the original registry file
    try:
        shutil.move(str(tmp_path), str(registry_path))
    except Exception as e:
        raise RuntimeError(  # pylint: disable=raise-missing-from
            f"Issue overwriting the old registry: {e}"
        )

    return registry_path

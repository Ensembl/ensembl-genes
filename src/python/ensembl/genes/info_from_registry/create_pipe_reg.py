import re
import shutil
from pathlib import Path
import os
import logging


from mysql_helper import mysql_update

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


def update_registry_path_and_create_entry(settings, server_info):
    new_registry_file = Path(settings["base_output_dir"]) / "Databases.pm"
    registry_file = os.path.join(
        os.environ.get("ENSCODE"),
        "ensembl-analysis",
        "scripts",
        "genebuild",
        "gbiab",
        "support_files",
        "Databases.pm"
    )

    if not new_registry_file.exists():
        shutil.copy(registry_file, new_registry_file)

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
        logger.error("Failed to update registry path.")


def create_registry_entry(settings, server_info, core_adaptor):
    update_registry_path_and_create_entry(settings, server_info)
    registry_path = Path(settings["base_output_dir"]) / "Databases.pm"

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
        raise RuntimeError(f"Issue overwriting the old registry: {e}")

    return registry_path
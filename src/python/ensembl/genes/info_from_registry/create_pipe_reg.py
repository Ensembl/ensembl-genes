import re
import shutil
from pathlib import Path

from src.python.ensembl.genes.info_from_registry.mysql_helper import mysql_update


def update_registry_path_and_create_entry(server_info, base_output_dir, registry_file):
    new_registry_file = base_output_dir / "Databases.pm"
    registry_file = re.sub(r"/+", "/", registry_file)

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
        print("Registry path updated successfully.")
    else:
        print("Failed to update registry path.")


def create_registry_entry(server_info, adaptors, base_output_dir, registry_file):
    update_registry_path_and_create_entry(server_info, base_output_dir, registry_file)
    registry_path = base_output_dir / "Databases.pm"

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
    core_string = db_string(adaptors["core_string"])

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
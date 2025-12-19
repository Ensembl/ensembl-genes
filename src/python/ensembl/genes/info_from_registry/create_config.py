#!/usr/bin/env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
create_config.py

This module handles copying and editing Ensembl Hive pipeline configuration
files based on runtime settings. It supports substituting key pipeline
parameters like user credentials, pipeline name, and release version.

Functions:
- copy_config: Copies the specified pipeline config file to the working directory.
- edit_config: Edits configuration values in the copied file using the provided settings.

Example usage:
    settings = {
        "config": "MyConfig.pm",
        "base_output_dir": "/my/output/path",
        "current_genebuild": "my_genebuild_104",
        "dbowner": "my_db_user",
        "pipeline_name": "my_pipeline",
        "password": "secure_pw",
        "user": "readonly_user",
        "user_r": "reader_user",
        "release_number": "104"
    }
    edit_config(settings)
"""

import os
import shutil
import re
import logging
from pathlib import Path
from typing import Dict, Any


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("pipeline_setup.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def copy_config(
    settings: Dict[str, Any], info_dict: Dict[str, Any], pipeline: str
) -> str:
    """
    Copies a configuration file from the ENSCODE environment into a local output directory.

    Parameters:
    ----------
    settings : dict
        Dictionary that must contain:
            - 'config': str, filename of the config (e.g., 'MyConfig.pm')
            - 'base_output_dir': str, path where the config will be copied

    Returns:
    -------
    str
        Path to the copied local config file.

    Raises:
    ------
    RuntimeError
        If the file cannot be copied for any reason (e.g., file not found, permissions).
    """

    original_config = os.path.join(
        str(os.environ.get("ENSCODE")),
        "ensembl-analysis",
        "modules",
        "Bio",
        "EnsEMBL",
        "Analysis",
        "Hive",
        "Config",
        settings["config"],
    )

    if pipeline == "anno":
        anno_parent = str(Path(info_dict["output_path"]).parent)
        local_config = os.path.join(anno_parent, settings["config"])
    elif pipeline == "main":
        local_config = os.path.join(info_dict["output_path"], settings["config"])
    else:
        raise ValueError(
            f"Unknown pipeline type: {pipeline}. Can't copy config. Please check taxon ID."
        )

    try:
        shutil.copy2(original_config, local_config)
        logging.info(  # pylint: disable=logging-fstring-interpolation
            f"Copied file from {original_config} to {local_config}"
        )
    except Exception as e:
        raise RuntimeError(  # pylint: disable=raise-missing-from
            f"Error copying file: {e}"
        )

    return local_config


def edit_config_anno(
    anno_settings: Dict[str, Any],
    settings: Dict[str, Any],
    info_dict: Dict[str, Any],
    pipeline: str,
    server_settings: Dict[str, Dict[str, Any]],
) -> None:
    """
    Edits specific parameter values in a copied Ensembl Hive config file.

    The function performs in-place substitution of key parameters such as:
    - current_genebuild
    - dbowner
    - pipeline_name
    - password
    - user
    - user_r
    - release_number

    Parameters:
    ----------
    anno_settings: dict
    settings : dict
        Dictionary with required keys:
            - 'config': str, name of the config file (e.g., 'MyConfig.pm')
            - 'base_output_dir': str, destination folder for the config
            - 'current_genebuild': str
            - 'dbowner': str
            - 'pipeline_name': str
            - 'password': str
            - 'user': str
            - 'user_r': str
            - 'release_number': str or int
    info_dict : dict
    pipeline : str

    Raises:
    ------
    FileNotFoundError
        If the copied file cannot be opened.
    KeyError
        If required settings keys are missing.
    """
    local_config = copy_config(anno_settings, info_dict, pipeline)

    with open(local_config, "r") as f:  # pylint: disable=unspecified-encoding
        content = f.read()

    # Replace the lines
    content = re.sub(
        r"'current_genebuild'\s*=>\s*[^,]+,",
        f"'current_genebuild'            => '{settings['current_genebuild']}',",
        content,
    )
    content = re.sub(
        r"'dbowner'\s*=>\s*[^,]+,",
        f"'dbowner'                      => '{settings['dbowner']}',",
        content,
    )
    content = re.sub(
        r"'pipeline_name'\s*=>\s*[^,]+,",
        f"'pipeline_name'                => '{settings['pipeline_name']}',",
        content,
    )
    content = re.sub(
        r"'password'\s*=>\s*[^,]+,",
        f"'password'                     => '{settings['password']}',",
        content,
    )
    content = re.sub(
        r"'user'\s*=>\s*[^,]+,",
        f"'user'                     => '{settings['user']}',",
        content,
    )
    content = re.sub(
        r"'user_r'\s*=>\s*[^,]+,",
        f"'user_r'                     => '{settings['user_r']}',",
        content,
    )

    content = re.sub(
        r"'release_number'\s*=>\s*[^,]+,",
        f"'release_number'                     => '{settings['release_number']}',",
        content,
    )

    content = re.sub(
        r"'pipe_db_server'\s*=>\s*[^,]+,",
        f"'pipe_db_server'                     => '{server_settings['pipeline_db']['db_host']}',",
        content,
    )

    content = re.sub(
        r"'pipe_db_port'\s*=>\s*[^,]+,",
        f"'pipe_db_port'                     => '{server_settings['pipeline_db']['db_port']}',",
        content,
    )
    content = re.sub(
        r"'dna_db_server'\s*=>\s*[^,]+,",
        f"'dna_db_server'                     => '{server_settings['core_db']['db_host']}',",
        content,
    )
    content = re.sub(
        r"'dna_db_port'\s*=>\s*[^,]+,",
        f"'dna_db_port'                     => '{server_settings['core_db']['db_port']}',",
        content,
    )

    content = re.sub(
        r"'registry_file'\s*=>\s*[^,]+,",
        f"'registry_file'                     => '{info_dict['registry_file']}',",
        content,
    )

    with open(local_config, "w") as f:  # pylint: disable=unspecified-encoding
        f.write(content)


def edit_config_main(  # pylint: disable=too-many-locals, too-many-branches
    settings: Dict[str, Any], info_dict: Dict[str, Any], pipeline: str
) -> str:
    """
    Edits specific parameter values in a copied Ensembl Hive config file safely,
    preserving comments and avoiding Perl syntax errors.

    Parameters
    ----------
    settings : dict
        Must contain 'config' and 'base_output_dir'.
    info_dict : dict
        Keys and values to replace in the config. Keys with value None are skipped.
    pipeline : str
        Type of pipeline, e.g., "anno" or "main".

    Raises
    ------
    FileNotFoundError
        If the copied config file cannot be found.
    ValueError
        If pipeline type is unknown.
    """

    local_config = copy_config(settings, info_dict, pipeline)

    stop_marker = "# No option below this mark should be modified"

    with open(local_config, "r") as f:  # pylint: disable=unspecified-encoding
        lines = f.readlines()

    new_lines = []
    stop_modifying = False

    for line in lines:
        # Replace the package declaration at the top
        if line.strip().startswith("package Bio::EnsEMBL::Analysis::Hive::Config::"):
            line = "package Genome_annotation_conf;\n"

        if stop_marker in line:
            stop_modifying = True

        if not stop_modifying:
            updated = False
            for key, value in info_dict.items():
                # Skip keys with None value
                if value is None:
                    continue

                # Pattern matches key => old_value, optionally with spaces and comments
                pattern = rf"^\s*{re.escape(key)}\s*=>\s*[^,]*,"

                if re.match(pattern, line):
                    # Preserve comment if present
                    if "#" in line:
                        parts = line.split("#", 1)
                        line_content = parts[0]
                        comment = "#" + parts[1].rstrip("\n")
                    else:
                        line_content = line
                        comment = ""

                    # Convert Python value to Perl representation
                    if isinstance(value, str):
                        perl_value = f"'{value}'"
                    elif isinstance(value, bool):
                        perl_value = "1" if value else "0"
                    else:
                        perl_value = str(value)

                    # Replace only the value before the comma
                    new_line = re.sub(r"=>\s*[^,]*,", f"=> {perl_value},", line_content)
                    if comment:
                        new_line = new_line.rstrip() + " " + comment + "\n"
                    else:
                        new_line = new_line.rstrip() + "\n"

                    new_lines.append(new_line)
                    updated = True
                    break

            if not updated:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Write back the updated config safely
    with open(local_config, "w", encoding="utf8") as f:
        f.writelines(new_lines)

    logging.info(  # pylint: disable=logging-fstring-interpolation
        f"Config edited successfully: {local_config}"
    )
    return local_config

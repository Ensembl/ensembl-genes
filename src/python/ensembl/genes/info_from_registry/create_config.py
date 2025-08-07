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


def copy_config(settings):
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
        os.environ.get("ENSCODE"),
        "ensembl-analysis",
        "modules",
        "Bio",
        "EnsEMBL",
        "Analysis",
        "Hive",
        "Config",
        settings["config"],
    )

    local_config = os.path.join(settings["base_output_dir"], settings["config"])

    try:
        shutil.copy2(original_config, local_config)
        logging.info(f"Copied file from {original_config} to {local_config}")
    except Exception as e:
        raise RuntimeError(f"Error copying file: {e}")

    return local_config


def edit_config(settings):
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

        Raises:
        ------
        FileNotFoundError
            If the copied file cannot be opened.
        KeyError
            If required settings keys are missing.
        """
    local_config = copy_config(settings)

    with open(local_config, "r") as f:
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

    with open(local_config, "w") as f:
        f.write(content)
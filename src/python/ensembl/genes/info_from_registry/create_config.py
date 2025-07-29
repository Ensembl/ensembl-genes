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
    Edits specified values in the copied file using a provided settings dictionary.

    Parameters:
        settings (dict): Dictionary with replacement values. Expected keys:
                         - "current_genebuild"
                         - "dbowner"
                         - "pipeline_name"
                         - "password"
                         - "user_r"
                         - "user"

    Raises:
        FileNotFoundError: If the file doesn't exist.
        KeyError: If required keys are missing in settings.
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
"""
Script to extract metadata, initialise an annotation pipeline using eHive,
and seed the pipeline with jobs based on a generated JSON file.

Modules used:
    - start_pipeline_from_registry.py (main)
    - seed_nonvert.py (seed_jobs_from_json)

Usage:
    python genebuild_init_pipeline.py --gcas path/to/gcas.txt --settings_file path/to/user_pipeline_settings.json
"""


import argparse
import logging
import re
from start_pipeline_from_registry import main as info
from seed_nonvert import seed_jobs_from_json
import subprocess
from pathlib import Path
from create_pipe_reg import update_registry_path_in_pipedb
import json



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

def init_pipeline_anno(config_file: str, hive_force_init: int = 1) -> str:
    """
    Initialize an eHive pipeline using a given config file.

    Args:
        config_file (str): Path to the .pm configuration file (e.g. EnsemblAnnoHelixer_conf.pm).
        hive_force_init (int): Whether to force initialization (default is 1).

    Returns:
        str: Extracted eHive database URL if successful.

    Raises:
        subprocess.CalledProcessError: If the `init_pipeline.pl` command fails.
        RuntimeError: If the expected MySQL URL is not found in the output.
    """
    logger.info("Innitialising anno pipeline")
    cmd = [
        "init_pipeline.pl",
        config_file,
        "--hive_force_init", str(hive_force_init)
    ]

    logger.info("Running: %s", " ".join(str(x) for x in cmd))
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        output = result.stdout
        for line in output.splitlines():
            match = re.search(r"(mysql://[^\s]+)", line)
            if match:
                return match.group(1)  # Just the URL part
        return 
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running init_pipeline.pl: {e.stdout}")
        raise

def init_pipeline_main(config_file: str, hive_force_init: int = 1) -> str:
    """
    Initialize an eHive pipeline using a given config file.

    Args:
        config_file (str): Path to the .pm configuration file.
        hive_force_init (int): Whether to force initialization (default is 1).

    Returns:
        str: Extracted eHive database URL if successful.

    Raises:
        subprocess.CalledProcessError: If the `init_pipeline.pl` command fails.
        RuntimeError: If the expected MySQL URL is not found in the output.
    """
    config_path = Path(config_file)
    config_dir = config_path.parent
    config_name = config_path.stem  # Gets filename without .pm extension
    
    cmd = [
        "init_pipeline.pl",
        config_name,  # Just the module name, not the full path
        "--hive_force_init", str(hive_force_init)
    ]

    logger.info("Running: %s", " ".join(str(x) for x in cmd))
    logger.info(f"Working directory: {config_dir}")
    
    try:
        result = subprocess.run(
            cmd, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True,
            cwd=str(config_dir)  # Run from the config directory
        )
        output = result.stdout
        for line in output.splitlines():
            match = re.search(r"(mysql://[^\s]+)", line)
            if match:
                return match.group(1)
        return None
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running init_pipeline.pl: {e.stdout}")
        raise
    

def main(gcas: str, settings_file: str, seed_url: str):
    """
    Main execution logic: extract metadata, initialize the pipeline, and seed jobs.

    Args:
        gcas (str): Path to file containing GCA accessions (one per line).
        settings_file (str): Path to settings file (JSON).
        seed_url (str): Existing URL to seed the pipeline. If provided initialisation will be skipped.

    Raises:
        FileNotFoundError: If no config file is found.
        RuntimeError: If multiple config files are found.
    """

    #Get info and create input JSON
    server_info, all_output_params, saved_paths = info(gcas, settings_file)


    #Save EHIVE URLs
    ehive_urls = {}

    #Anno
    if saved_paths.get("anno"):
        json_file = saved_paths["anno"]

        if seed_url:
            logger.info("Skipping pipeline initialisation.")
            logger.info("Using provided seed URL: %s", seed_url)
            ehive_urls["anno"] = seed_url
        else:
            logger.info("No seed URL provided. Initialising new pipeline.")
            first_key = next(iter(all_output_params))
            output_path = Path(all_output_params[first_key]["output_path"])
            parent_dir = output_path.parent
            conf_files = list(parent_dir.glob("*_conf.pm"))

            if not conf_files:
                raise FileNotFoundError(f"No .conf file found in {parent_dir}")
            if len(conf_files) > 1:
                raise RuntimeError(f"Multiple .conf files found in {parent_dir}: {conf_files}")

            conf_path = conf_files[0]
            logger.info(f"Config found at {conf_path}")

            ehive_url = init_pipeline_anno(conf_path)
            ehive_urls["anno"] = ehive_url
            logger.info(f"Extracted EHIVE_URL: {ehive_url}")

            update_registry_path_in_pipedb(parent_dir, server_info)
            seed_url = ehive_url  # unify below

        # Common seeding step
        logger.info(f"Seeding non-vertebrate pipeline from {json_file}")
        with open(json_file) as f:
            params = json.load(f)

        seed_jobs_from_json(
            json_file=params,
            analysis_id=1,
            ehive_url=seed_url)


    # --- Handle main (one JSON per GCA) ---
    if saved_paths.get("main"):
        ehive_urls["main"] = {}
        for gca, json_path in saved_paths["main"].items():
            logger.info(f"Initialising vertebrate pipeline for {gca} from {json_path}")

            output_path = Path(all_output_params[gca]["output_path"])
            conf_files = list(output_path.glob("*_conf.pm"))

            if not conf_files:
                raise FileNotFoundError(f"No .conf file found in {parent_dir}")
            elif len(conf_files) > 1:
                raise RuntimeError(f"Multiple .conf files found in {parent_dir}: {conf_files}")

            conf_path = conf_files[0]

            ehive_url = init_pipeline_main(str(conf_path))
            ehive_urls["main"][gca] = ehive_url
            logger.info(f"Extracted EHIVE_URL for {gca}: {ehive_url}")

    logger.info(f"Initalised pipelines. EHIVE_URL(s): {ehive_urls}")

    output_dir_current = Path.cwd()
    ehive_url_file = output_dir_current / "ehive_urls.json"

    with open(ehive_url_file, "w") as f:
        json.dump(ehive_urls, f, indent=2)

    logger.info(f"Saved EHIVE URLs to {ehive_url_file}")

    return ehive_urls

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract metadata, initialize and seed annotation pipelines."
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

    parser.add_argument(
        "--seed_url",
        type=str,
        required=False,
        help="Pipeline seed URL. If provided initialisation will be skipped."
    )

    args = parser.parse_args()

    ehive_urls = main(args.gcas, args.settings_file, args.seed_url)


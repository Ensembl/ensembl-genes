"""
Script to extract metadata, initialise an annotation pipeline using eHive,
and seed the pipeline with jobs based on a generated JSON file.

Modules used:
    - start_pipeline_from_registry.py (main)
    - seed_nonvert.py (seed_jobs_from_json)

Usage:
    python genebuild_nonvert.py --gcas path/to/gcas.txt --pipeline anno --settings_file non_vert_pipeline_params.json
"""


import argparse
import logging
import re
from start_pipeline_from_registry import main as info
from seed_nonvert import seed_jobs_from_json
import subprocess
from pathlib import Path



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

def init_pipeline(config_file: str, hive_force_init: int = 1) -> str:
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
    

def main(gcas: str, pipeline: str, settings_file: str):
    """
    Main execution logic: extract metadata, initialize the pipeline, and seed jobs.

    Args:
        gcas (str): Path to file containing GCA accessions (one per line).
        pipeline (str): Pipeline type to initialize (e.g., "anno").
        settings_file (str): Path to settings file (JSON).

    Raises:
        FileNotFoundError: If no config file is found.
        RuntimeError: If multiple config files are found.
    """

    #Get info and create input JSON
    all_output_params, output_json_path = info(gcas, pipeline, settings_file)

    #Initialise pipeline
    if pipeline == "anno":
        first_key = next(iter(all_output_params))
        output_path = Path(all_output_params[first_key]["output_path"])
        parent_dir = output_path.parent
        conf_files = list(parent_dir.glob("*_conf.pm"))

        if not conf_files:
            raise FileNotFoundError(f"No .conf file found in {parent_dir}")
        elif len(conf_files) > 1:
            raise RuntimeError(f"Multiple .conf files found in {parent_dir}: {conf_files}")
        
        conf_path = conf_files[0]

        ehive_url = init_pipeline(conf_path)
        logger.info(f"Extracted EHIVE_URL: {ehive_url}")
    
        #Seed pipeline
        seed_jobs_from_json(
            json_file=str(output_json_path),
            analysis_id=1,
            ehive_url=ehive_url)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract metadata, initialize annotation pipeline and seed pipeline."
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

    results = main(args.gcas, args.pipeline, args.settings_file)


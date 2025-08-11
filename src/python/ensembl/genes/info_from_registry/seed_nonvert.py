"""
seed_nonvert.py

This script allows you to seed jobs into an eHive pipeline using a JSON file or a dictionary of job parameters.
It converts Python dictionaries to Perl hash syntax and invokes the `seed_pipeline.pl` script via subprocess.

Typical usage example:
    python seed_jobs.py -j jobs.json -a 1 -u mysql://user:pass@host/db

Functions:
    - dict_to_perl_hash: Recursively converts a Python dictionary into a Perl-style hash string.
    - seed_jobs_from_json: Loads parameters and invokes eHive seeding commands.
    - main: Parses command-line arguments and calls the seeding function.
"""
import json
import logging
import subprocess
import argparse
from typing import Union

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

def dict_to_perl_hash(d):
    """
    Recursively convert a Python dictionary to a Perl hash string.

    Args:
        d (dict): Python dictionary to convert.

    Returns:
        str: A string representing the dictionary in Perl hash format.
    """
    items = []
    for k, v in d.items():
        key_str = f"'{k}'"
        if isinstance(v, dict):
            val_str = dict_to_perl_hash(v)
        elif isinstance(v, str):
            val_str = f"'{v}'"
        elif v is None:
            val_str = "undef"
        else:
            val_str = str(v)
        items.append(f"{key_str} => {val_str}")
    return "{{{}}}".format(", ".join(items))


def seed_jobs_from_json(
    json_file: Union[str, dict],
    analysis_id: int,
    ehive_url: str,
):
    """
    Seed jobs in an eHive pipeline using job parameters from a JSON file or dictionary.

    Args:
        json_file (Union[str, dict]): Path to a JSON file or a Python dictionary with job parameters.
        analysis_id (int): eHive analysis ID to assign jobs to.
        ehive_url (str): URL of the eHive database (e.g., mysql://user:pass@host/db).

    Raises:
        subprocess.CalledProcessError: If the seeding script fails.
    """

    if isinstance(json_file, str):
        with open(json_file) as f:
            params = json.load(f)
    else:
        params = json_file  # assume dict already

    for input_id, param_dict in params.items():
        perl_hash = dict_to_perl_hash(param_dict)
        cmd = [
            "seed_pipeline.pl",
            "-analysis_id", str(analysis_id),
            "-input_id", perl_hash,
            "-url", ehive_url
        ]
        subprocess.run(cmd, check=True)
        logging.info("Seeding complete")


def main():
    parser = argparse.ArgumentParser(description="Seed eHive pipeline jobs from a JSON file")
    parser.add_argument(
        "-j", "--json_file",
        required=True,
        help="Path to the JSON file with job parameters"
    )
    parser.add_argument(
        "-a", "--analysis_id",
        type=int,
        default=1,
        help=f"Analysis ID to seed (default: 1)"
    )
    parser.add_argument(
        "-u", "--url",
        required=True,
        help="EHIVE URL"
    )
    args = parser.parse_args()

    seed_jobs_from_json(
        json_file=args.json_file,
        analysis_id=args.analysis_id,
        ehive_url=args.url
    )

if __name__ == "__main__":
    main()
import os
import json
import subprocess
import argparse

ANALYSIS_ID = 1  # Change as appropriate

def dict_to_perl_hash(d):
    # Recursively convert a Python dictionary to a Perl hash string with single quotes
    items = []
    for k, v in d.items():
        if isinstance(v, dict):
            v_str = dict_to_perl_hash(v)
        elif isinstance(v, str):
            v_str = f"'{v}'"
        else:
            v_str = str(v)
        items.append(f"'{k}' => {v_str}")
    return "{" + ", ".join(items) + "}"

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
        default=ANALYSIS_ID,
        help=f"Analysis ID to seed (default: {ANALYSIS_ID})"
    )

    args = parser.parse_args()

    json_file = args.json_file
    analysis_id = args.analysis_id

    with open(json_file) as f:
        params = json.load(f)

    ehive_url = os.environ.get("EHIVE_URL")
    if not ehive_url:
        raise RuntimeError("EHIVE_URL environment variable not set")

    for input_id, param_dict in params.items():
        perl_hash = dict_to_perl_hash(param_dict)
        # Optionally include the input_id itself as a parameter:
        # perl_hash = dict_to_perl_hash({**param_dict, "input_id": input_id})
        cmd = [
            "seed_pipeline.pl",
            "-analysis_id", str(analysis_id),
            "-input_id", perl_hash,
            "-url", ehive_url
        ]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()

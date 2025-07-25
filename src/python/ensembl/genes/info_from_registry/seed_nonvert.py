import json
import subprocess
import argparse

def dict_to_perl_hash(d):
    """Recursively convert a Python dict to a Perl hash string."""
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
    json_file = args.json_file
    analysis_id = args.analysis_id
    ehive_url = args.url

    # Load job parameters from JSON
    with open(json_file) as f:
        params = json.load(f)

    for input_id, param_dict in params.items():
        perl_hash = dict_to_perl_hash(param_dict)
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
#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
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

import json
import argparse
import re
import sys


def extract_completeness(busco_str) -> float:
    """
    Extract the BUSCO completeness value from a string like:
    "C:99.2%[S:98.4%,D:0.8%],F:0.6%,M:0.2%,n:3640"
    and return it as a float.
    If the format is invalid, raise a ValueError.

    Args:
        busco_str (str): The BUSCO string to parse.
    Return: The completeness percentage as a float.
    """
    match = re.search(r"C:([\d.]+)%", busco_str)
    if match:
        return float(match.group(1))
    raise ValueError(f"Invalid BUSCO format: {busco_str}")


def evaluate_busco(
    genome_json_path: str,
    protein_json_path: str,
    min_range_protein_score: int,
    max_range_protein_score: int,
    diff_prot_gen_mode: int,
) -> None:
    """Evaluate BUSCO scores from genome and protein JSON files.
    This function loads the BUSCO JSON files for genome and protein,
    extracts the completeness scores, and compares them to determine
    if the protein BUSCO score is significantly higher than the genome BUSCO score.
    The criteria for significance are:
    - If the protein BUSCO score is >= 70%, it is considered significant.
    - If the protein BUSCO score is between 50% and 70%, it is considered significant
      if the difference between the protein and genome BUSCO scores is >= 10%.
    If the genome or protein BUSCO JSON files cannot be loaded, an error message is printed
    and the function returns False.

    Args:
        genome_json_path (str): _path to the genome BUSCO JSON file_
        protein_json_path (str): _path to the protein BUSCO JSON file_
        min_range_protein_score (int): Lowest threshold to analyse busco score in protein mode
        max_range_protein_score (int): Highest threshold to analyse busco score in protein mode
        diff_prot_gen_mode (int): Max difference between Busco in protein and in genome mode

    Raises:
        KeyError: If the expected '.busco' key is not found in either JSON file.

    Returns:
        exit 0 if the protein BUSCO score soddisfy criteria, otherwise exit 42
    """
    try:
        with open(genome_json_path) as g_file:  # pylint:disable=unspecified-encoding
            genome_data = json.load(g_file)
    except Exception as err_msg:  # pylint: disable=broad-except
        print(f"Failed to load genome BUSCO JSON: {err_msg}")
        sys.exit(42)

    try:
        with open(protein_json_path) as p_file:  # pylint:disable=unspecified-encoding
            protein_data = json.load(p_file)
    except Exception as err_msg:  # pylint: disable=broad-except
        print(f"Failed to load protein BUSCO JSON: {err_msg}")
        sys.exit(42)

    genome_busco_key = next((k for k in genome_data if k.endswith(".busco")), None)
    protein_busco_key = next((k for k in protein_data if k.endswith(".busco")), None)

    if not genome_busco_key or not protein_busco_key:
        raise KeyError("Missing '.busco' key in one of the JSON files")

    genome_busco_score = extract_completeness(genome_data[genome_busco_key])
    protein_busco_score = extract_completeness(protein_data[protein_busco_key])

    print(f"Genome BUSCO completeness: {genome_busco_score}%")
    print(f"Protein BUSCO completeness: {protein_busco_score}%")

    if protein_busco_score >= max_range_protein_score:
        print(
            f"Protein BUSCO completeness {protein_busco_score}% is above threshold (70%)"
        )
        sys.exit(0)
    if min_range_protein_score < protein_busco_score < max_range_protein_score:
        difference = protein_busco_score - genome_busco_score
        print(f"Difference (protein - genome): {difference:.2f}%")
        if difference >= diff_prot_gen_mode:
            sys.exit(0)
        else:
            sys.exit(42)
    print(f"Protein BUSCO completeness {protein_busco_score}% is too low")
    sys.exit(42)


def main():
    """Main function to parse command line arguments and evaluate BUSCO scores."""
    parser = argparse.ArgumentParser(
        description="Evaluate BUSCO scores from genome and protein JSON files."
    )
    parser.add_argument(
        "--genome", required=True, help="Path to genome BUSCO JSON file"
    )
    parser.add_argument(
        "--protein", required=True, help="Path to protein BUSCO JSON file"
    )
    parser.add_argument(
        "--min_range_protein_score",
        type=int,
        required=False,
        default=50,
        help="Lowest threshold to analyse busco score in protein mode",
    )
    parser.add_argument(
        "--max_range_protein_score",
        type=int,
        required=False,
        default=70,
        help="Highest threshold to analyse busco score in protein mode",
    )
    parser.add_argument(
        "--diff_prot_gen_mode",
        type=int,
        required=False,
        default=10,
        help="Max difference between Busco in protein and in genome mode",
    )
    args = parser.parse_args()

    try:  # pylint: disable=broad-except
        evaluate_busco(
            args.genome,
            args.protein,
            args.min_range_protein_score,
            args.max_range_protein_score,
            args.diff_prot_gen_mode,
        )

    except Exception as err_msg:  # pylint: disable=broad-except
        print(f"ERROR: {err_msg}")
        sys.exit(42)


if __name__ == "__main__":
    main()

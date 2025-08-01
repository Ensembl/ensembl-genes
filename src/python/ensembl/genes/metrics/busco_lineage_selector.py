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

import urllib.request
import argparse
from pathlib import Path
from typing import Any, List
#import requests
import json

def get_dataset_match(ncbi_url: str, dataset: list) -> List[Any]:
    """
    Get taxonomy tree from ncbi taxonomy datasets and find the closest match with the input list


    Args:
        ncbi_url (str): Ncbi dataset url
        dataset (list): list of data to match

    Returns:
        str: closest match in the dataset list

    Raises:
        requests.HTTPError: If an HTTP error occurs during the API request.
        Exception: If any other error occurs during the function's operation.

    """

    try:
        # Fetch data from the URL
        with urllib.request.urlopen(ncbi_url, timeout=10) as response:
            # Read the response and decode it
            data = response.read().decode('utf-8')
            # Parse the JSON data
            json_data = json.loads(data)

            # Extract classification names
            parents = json_data["reports"][0]["taxonomy"]["parents"]
            # Variable to store the matched result
            matched_value = None

            # Match parents against the dictionary
            for parent_id in reversed(parents):
                parent_id_str = str(parent_id)
                if parent_id_str in dataset:
                    matched_value = dataset[parent_id_str]
                    break       
    except urllib.error.URLError as url_err:
        print(f"URL error occurred: {url_err}")
    except json.JSONDecodeError as json_err:
        print(f"Error decoding JSON: {json_err}")
    #print (matched_value)    
    return matched_value


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Clade selector arguments")
    parser.add_argument(
        "-d",
        "--datasets",
        type=str,
        help="Path to file containing list of datasets (one per line)",
        required=True,
    )
    parser.add_argument("-t", "--taxon_id", type=str, help="Taxon id ", required=True)
    parser.add_argument("--output", type=str, help="Output file", default="stdout")
    parser.add_argument(
        "--ncbi_url",
        type=str,
        help="NCBI dataset url",
        default="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/",
    )
    return parser.parse_args()


def main():
    """Entry-point."""
    args = parse_args()

    ncbi_url = f"{args.ncbi_url}/{args.taxon_id}/dataset_report"

    with open(Path(args.datasets), "r") as file:
        #datasets = [line[: max(line.find(" "), 0) or None] for line in file]
        datasets = json.load(file)
    clade_match = get_dataset_match(ncbi_url, datasets)

    if not clade_match:
        raise ValueError("No match found")

    if args.output == "stdout":  # pylint:disable=no-else-return
        #print(clade_match[0].strip("\n"))
        print(clade_match)
    else:
        with open(args.output, "w+") as output:
            if clade_match[0] == args.species:
                output.write(clade_match[1])
            else:
                output.write(clade_match[0])

    return None


if __name__ == "__main__":
    main()
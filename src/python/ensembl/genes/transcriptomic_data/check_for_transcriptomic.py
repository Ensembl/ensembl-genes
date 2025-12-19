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
"""Check the availability for short and long read data from ENA website given a taxon id"""
import csv
from pathlib import Path
import argparse
import requests


def ena_rest_api(query: str) -> int:
    """Call to ENA API

    Args:
        query (str): query string to search ENA database

    Returns:
        int: number of runs found
    """
    search_url = f"https://www.ebi.ac.uk/ena/portal/api/search?display=report&query={query}&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform"  # pylint: disable=line-too-long
    search_result = requests.get(search_url, timeout=60)
    results = search_result.text.strip().split("\n")[1:]
    return len(results)


def check_data_from_ena(  # pylint: disable=too-many-locals
    taxon_id: int,
    tree: bool,
) -> dict:
    """
    Query ENA API to get short or long read data

    Args:
        taxon_id (int): NCBI taxon id
        tree (bool): whether to include subordinate taxa

    Returns:
        dict: number of runs found for each data type
    """

    TEXT_FORMAT: dict[str, str] = {  # pylint:disable=invalid-name
        "BOLD": "\033[1m",
        "UNDERLINE": "\033[4m",
        "END": "\033[0m",
    }
    if tree:
        query = f"tax_tree({taxon_id})"
    else:
        query = f"tax_eq({taxon_id})"

    query_short_paired = (
        query
        + " AND instrument_platform=ILLUMINA AND library_layout=PAIRED"
        + " AND library_source=TRANSCRIPTOMIC"
    )
    query_short_single = (
        query
        + " AND instrument_platform=ILLUMINA AND library_layout=SINGLE"
        + " AND library_source=TRANSCRIPTOMIC"
    )
    query_pacbio = (
        query + " AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC"
    )
    query_onp = (
        query
        + " AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC"
    )

    short_paired_runs = ena_rest_api(query_short_paired)
    short_single_runs = ena_rest_api(query_short_single)
    pacbio_read_runs = ena_rest_api(query_pacbio)
    onp_read_runs = ena_rest_api(query_onp)

    print(
        TEXT_FORMAT["BOLD"]
        + "Short-read paired-end illumina data available! "
        + TEXT_FORMAT["END"]
        + f"Found {short_paired_runs} runs."
    )
    print(
        TEXT_FORMAT["BOLD"]
        + "Short-read single-end illumina data available! "
        + TEXT_FORMAT["END"]
        + f"Found {short_single_runs} runs."
    )
    print(
        TEXT_FORMAT["BOLD"]
        + "Long-read PacBio data available! "
        + TEXT_FORMAT["END"]
        + f"Found {pacbio_read_runs} runs."
    )
    print(
        TEXT_FORMAT["BOLD"]
        + "Long_read ONP data available! "
        + TEXT_FORMAT["END"]
        + f"Found {onp_read_runs} runs."
    )

    return {
        "Short-read paired-end illumina": short_paired_runs,
        "Short-read single-end illumina": short_single_runs,
        "Long-read PacBio": pacbio_read_runs,
        "Long-read ONP": onp_read_runs,
    }


class InputSchema(argparse.ArgumentParser):
    """Input arguments"""

    def __init__(self):
        super().__init__()

        self.add_argument("-t", "--taxon_id", type=str, required=False, help="Taxon id")

        self.add_argument(
            "--tree",
            action="store_true",
            required=False,
            help="Turn on the 'Include subordinate taxa' option in your query to ENA",
        )
        self.add_argument(
            "--file",
            type=str,
            required=False,
            help="Path to the file containing a list of taxon ids",
        )
        self.add_argument("--output_dir", required=True, help="Output directory path")


def main() -> None:
    """Entrypoint"""
    parser = InputSchema()
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    taxon_ids = [args.taxon_id]  # Start with the single taxon_id

    if args.file:
        with open(args.file, "r", encoding="utf-8") as input_file:
            taxon_ids = input_file.read().splitlines()

    # Prepare CSV file
    csv_file = output_dir / "taxon_summary.csv"
    with open(csv_file, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)

        # Write the header row
        writer.writerow(
            [
                "Taxon ID",
                "Short-read paired-end illumina",
                "Short-read single-end illumina",
                "Long-read PacBio",
                "Long-read ONP",
            ]
        )

        for taxon_id in taxon_ids:
            results = check_data_from_ena(int(taxon_id), args.tree)
            writer.writerow(
                [
                    taxon_id,
                    results["Short-read paired-end illumina"],
                    results["Short-read single-end illumina"],
                    results["Long-read PacBio"],
                    results["Long-read ONP"],
                ]
            )


if __name__ == "__main__":
    main()

# See the NOTICE file distributed with this work for additional information #pylint: disable=missing-module-docstring
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
"""Download short and long read data from ENA website"""
import os.path
from pathlib import Path
from typing import List
import argparse
import requests
from requests.exceptions import RequestException, HTTPError, ConnectionError, Timeout
import re

def get_sample_info(accession: str) -> List:
    """Get info about sample name and description for the run accession"""
    biosample_url = f"https://www.ebi.ac.uk/biosamples/samples/{accession}"

    try:
        response = requests.get(biosample_url)
        response.raise_for_status()  # Raise an HTTPError if the request was not successful

        biosample_data = response.json()
        
        if "characteristics" in biosample_data and "tissue" in biosample_data["characteristics"]:
            sample = biosample_data["characteristics"]["tissue"][0]["text"].lower()
        elif "characteristics" in biosample_data and "organism_part" in biosample_data["characteristics"]:
            sample = biosample_data["characteristics"]["organism_part"][0]["text"].lower()
        else:
            sample = accession

        if "characteristics" in biosample_data and "description" in biosample_data["characteristics"]:
            description = biosample_data["characteristics"]["description"][0]["text"]
            if len(description) > 250:
                description = accession
        else:
            description = accession

        sample = re.sub(r'[ ;\(\)\/\\]', '_', sample)
        #remove punctuation
        sample = re.sub(r'[!\"#$%&()*\+,\-\'.\/:;<=>?@\[\]^`{|}~]', '', sample)
        multi_tissuesREGEX = ("([a-zA-Z]+\,)+")
        if re.search(multi_tissuesREGEX, sample):
            sample = "mixed_tissues"
        
        return (sample, description)

    except (RequestException, HTTPError, ConnectionError, Timeout) as e:
        print(f"An error occurred while fetching data from {biosample_url}: {str(e)}")
        # Handle the error here, you can log it or take other appropriate actions.
        return ("unknown", "unknown")

def get_data_from_ena(taxon_id: int, read_type: str, tree: bool) -> List[str]:
    """Query ENA API to get short or long read data"""
    csv_data = []

    if tree:
        query = f"tax_tree({taxon_id})"
    else:
        query = f"tax_eq({taxon_id})"
        
    if read_type == "short":
        query += " AND instrument_platform=ILLUMINA AND library_layout=PAIRED"
    else:
        query += " AND (instrument_platform=OXFORD_NANOPORE OR instrument_platform=PACBIO_SMRT)"

    query += " AND library_source=TRANSCRIPTOMIC"

    search_url = f"https://www.ebi.ac.uk/ena/portal/api/search?display=report&query={query}&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform,fastq_md5"

    search_result = requests.get(search_url)
    results = search_result.text.strip().split("\n")[1:]
    
    is_paired = "1"
    is_mate_1 = "-1"
    read_length = "1"
    is_plus_13 = "0"
    centre = "ENA"

    for row in results:
        row_data = row.split('\t')
        sample_accession = row_data[0]
        run_accession = row_data[1]

        multi_samples = sample_accession.split(';')
        if len(multi_samples) > 1:
            sample, description = "multiple", "multiple"
        else:
            sample, description = get_sample_info(sample_accession)
        
        try:
            read_count = int(row_data[3])
            instrument_platform = row_data[4]
        except ValueError:
            read_count = "0"
            instrument_platform = row_data[3]

        file_entries = row_data[2].split(";")
        md5_entries = row_data[5].split(";")

        if len(file_entries) == 2:
            pass  # Nothing special to do, continue as normal
        elif len(file_entries) == 3:
            # Identify the unwanted entry and remove it
            file_entries = [f for f in file_entries if "_1.fastq.gz" in f or "_2.fastq.gz" in f]
            # Assuming the order of md5 corresponds to files and the unwanted file is in the middle
            md5_entries = [md5_entries[i] for i, f in enumerate(row_data[2].split(";")) if "_1.fastq.gz" in f or "_2.fastq.gz" in f]
        else:
            print("Warning: Unexpected number of file entries, skipping "+run_accession)
            next  # Skip further processing for this row

        # Only proceed if we have exactly 2 entries after any necessary filtering
        if len(file_entries) == 2:
            if "ftp" in row_data[2] and ";" in row_data[5]:
                for file, md5_file in zip(file_entries, md5_entries):
                    file_path = os.path.basename(file)
                    md5_file_value = md5_file
                    csv_data.append(
                        (
                            sample,
                            run_accession,
                            is_paired,
                            file_path,
                            is_mate_1,
                            read_length,
                            is_plus_13,
                            centre,
                            instrument_platform,
                            description,
                            file,
                            md5_file_value,
                        )
                    )

    return csv_data

class InputSchema(argparse.ArgumentParser):
    """Input arguments"""

    def __init__(self):
        super().__init__()

        self.add_argument(
            "-t", "--taxon_id", type=int, required=True, help="Taxon id"
        )
        self.add_argument(
            "--tree",
            action='store_true',
            required=False,
            help="Turn on the 'Include subordinate taxa' option in your query to ENA",
        )
        self.add_argument(
            "-f", "--csv_file", type=str, required=True, help="Output file path (csv format)"
        )
        self.add_argument(
            "-r", "--read_type",
            choices=["short", "long"],
            required=True,
            help="Specify the type of transcriptomic data to download ['short', 'long']",
        )
        self.add_argument(
            "-l", "--limit",
            type=int,
            required=False,
            help="The number of runs to be included in your csv file - consider that 1 run = 2 files, so setting '-l 50' will result in 100 lines in your csv (WARNING: this limit does not consider quality of data, it is a simple subsampling)",
        )
        
def main() -> None:
    """Entrypoint"""
    parser = InputSchema()
    args = parser.parse_args()

    if os.path.isfile(args.csv_file) and os.stat(args.csv_file).st_size > 0:
        print("File " + args.csv_file + " exists, and is not empty, will not overwrite!")
        exit
    else:
        try:
            csv_data = get_data_from_ena(args.taxon_id, args.read_type, args.tree)

            if args.limit:
                csv_data = csv_data[:(args.limit * 2)]
            
            with open(Path(args.csv_file), "w", encoding="utf8") as csv_file:
                for row in csv_data:
                    csv_file.write("\t".join(row) + "\n")
        except (RequestException, HTTPError, ConnectionError, Timeout) as e:
            print(f"An error occurred during the data retrieval process: {str(e)}")
            # Handle the error here, you can log it or take other appropriate actions.

if __name__ == "__main__":
    main()

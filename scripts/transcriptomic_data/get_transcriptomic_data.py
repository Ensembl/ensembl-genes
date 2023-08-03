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


def get_sample_info(accession: str) -> List:
    """Get info about sample name and description for the run accession"""
    biosample_url = f"https://www.ebi.ac.uk/biosamples/samples/{accession}"
    biosample_data = requests.get(biosample_url).json()

    # sample name will be set as the tissue type or organism part fields from ENA BioSample
    # if neither exist, sample name is "unknown"
    # this requires a project for working out the best way to find sample names that make sense

    if "characteristics" in biosample_data and "tissue_type" in biosample_data["characteristics"]:
        sample = biosample_data["characteristics"]["tissue_type"][0]["text"]
    elif "characteristics" in biosample_data and "organism_part" in biosample_data["characteristics"]:
        sample = biosample_data["characteristics"]["organism_part"][0]["text"]
    else:
        sample = "unknown"
    # description is pulled from ENA BioSample
    # if the description does not exist, the sample accession is provided for tracibility

    # try:
    #    description = biosample_data['characteristics']['description'][0]['text']
    # except KeyError:
    #    description = accession

    if "characteristics" in biosample_data and "description" in biosample_data["characteristics"]:
        description = biosample_data["characteristics"]["description"][0]["text"]
    else:
        description = accession

    # we will likely come across more characters that should be removed from sample names
    # improvements on this should come with the sample names project
    replace_chars = '()/\\'
    for i in replace_chars:
         sample = sample.replace(i, "")
    sample = sample.replace(" ", "_")
    return (sample, description)


def get_data_from_ena(# pylint: disable=too-many-locals
    taxon_id: int, read_type: str
) -> List[str]:
    """Query ENA API to get short or long read data"""
    csv_data = []
    query = f"tax_eq({taxon_id})"
    if read_type == "short":
        query += " AND instrument_platform=ILLUMINA AND library_layout=PAIRED"
    else:
        query += " AND (instrument_platform=OXFORD_NANOPORE OR instrument_platform=PACBIO_SMRT)"

    query += " AND library_source=TRANSCRIPTOMIC"

    search_url = f"https://www.ebi.ac.uk/ena/portal/api/search?display=report&query={query}&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform"  # pylint: disable=line-too-long
    search_result = requests.get(search_url)
    results = search_result.text.strip().split("\n")[1:]

    # fields for the csv that we don't need any more but could break
    # downstream things so we'll keep them for now
    is_paired = "1"  # we always used paire-end data
    is_mate_1 = "-1"
    # read_length is always incorrect when calculated from ENA and is not needed with STAR anyway
    read_length = "1"
    is_plus_13 = "0"  # outdated sequencing info
    centre = "ENA"  # we get everything from ENA anyway
    
    for row in results:
        row_data = row.split()
        run_accession = row_data[0]
        sample_accession = row_data[1]
        sample, description = get_sample_info(sample_accession)
        try:
            read_count = int(row_data[3])
            instrument_platform = row_data[4]
        except ValueError:
            read_count = "0" #read count not always available in meta data on ena so will set to 0 in these cases - we need a better solution for this!
            instrument_platform = row_data[3]
        for file in row_data[2].split(";"):
            file_path=os.path.basename(file)
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
                )
            )
    return csv_data


class InputSchema(argparse.ArgumentParser):
    """Input arguments"""

    def __init__(self):
        super().__init__()

        self.add_argument(
            "-t", "--taxon_id", type=str, required=True, help="Taxon id"
        )
        self.add_argument(
            "-f", "--csv_file", type=str, required=True, help="Output file path (csv format)"
        )
        self.add_argument(
            "--read_type",
            choices=["short", "long"],
            required=True,
            help="Specify the type of transcriptomic data to download ['short', 'long']",
        )

def main() -> None:
    """Entrypoint"""
    parser=InputSchema()
    args = parser.parse_args()

    if os.path.isfile(args.csv_file):
        print ("File "+args.csv_file+" exists, will not overwrite!")
        exit
    else:
        csv_data = get_data_from_ena(args.taxon_id, args.read_type)

        with open(Path(args.csv_file), "w", encoding="utf8") as csv_file:
            for row in csv_data:
                csv_file.write(line + "\n")

if __name__ == "__main__":
    main()

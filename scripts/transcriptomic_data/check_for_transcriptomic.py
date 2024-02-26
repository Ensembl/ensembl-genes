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
"""Check the availability for short and long read data from ENA website given a taxon id"""
import os.path
from pathlib import Path
from typing import List
import argparse
import requests

def ena_rest_api(query:str)->int:
    search_url = f"https://www.ebi.ac.uk/ena/portal/api/search?display=report&query={query}&domain=read&result=read_run&fields=sample_accession,run_accession,fastq_ftp,read_count,instrument_platform"  # pylint: disable=line-too-long
    search_result = requests.get(search_url)
    results = search_result.text.strip().split("\n")[1:]
    return len(results)

def check_data_from_ena(# pylint: disable=too-many-locals
        taxon_id: int,
        tree: bool
) -> None:
    """Query ENA API to get short or long read data"""

    if tree:
        query = f"tax_tree({taxon_id})"
    else:
        query = f"tax_eq({taxon_id})"
    
    query_short_paired=query+f" AND instrument_platform=ILLUMINA AND library_layout=PAIRED AND library_source=TRANSCRIPTOMIC"
    query_short_single=query+f" AND instrument_platform=ILLUMINA AND library_layout=SINGLE AND library_source=TRANSCRIPTOMIC"
    query_pacbio=query+f" AND instrument_platform=PACBIO_SMRT AND library_source=TRANSCRIPTOMIC"
    query_onp=query+f" AND instrument_platform=OXFORD_NANOPORE AND library_source=TRANSCRIPTOMIC"
    
    short_paired_runs=ena_rest_api(query_short_paired)
    short_single_runs=ena_rest_api(query_short_single)
    pacbio_read_runs =ena_rest_api(query_pacbio)
    onp_read_runs=ena_rest_api(query_onp)


    print (text.BOLD+f"Short-read paired-end illumina data available! "+text.END+f"Found {short_paired_runs} runs.")
    print (text.BOLD+f"Short-read single-end illumina data available! "+text.END+f"Found {short_single_runs} runs.")
    print (text.BOLD+f"Long-read PacBio data available! "+text.END+f"Found {pacbio_read_runs} runs.")
    print (text.BOLD+f"Long_read ONP data available! "+text.END+f"Found {onp_read_runs} runs.")
        
class text:
    """formatting set"""

    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


class InputSchema(argparse.ArgumentParser):
    """Input arguments"""
    def __init__(self):
        super().__init__()

        self.add_argument(
            "-t", "--taxon_id", type=str, required=True, help="Taxon id"
        )

        self.add_argument(
            "--tree", action='store_true', required=False, help="Turn on the 'Include subordinate taxa' option in your query to ENA"
        )
                
def main() -> None:
    """Entrypoint"""
    parser=InputSchema()
    args = parser.parse_args()

    check_data_from_ena(args.taxon_id, args.tree)

if __name__ == "__main__":
    main()


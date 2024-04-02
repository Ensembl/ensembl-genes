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

import sys
import argparse
import requests
import pymysql
import json
from pathlib import Path
from collections import Counter
from typing import List, Tuple, Any, Dict

with open("./bioproject_tracking_config.json", "r") as f:
    config = json.load(f)

def mysql_fetch_data(query: str, database: str, host: str, port: int, user: str) -> List[Tuple[Any, ...]]:
    """
    Executes a given SQL query on a MySQL database and fetches the result.

    This function establishes a connection to the MySQL database using provided connection details.
    It then executes the given query using a cursor obtained from the connection. After executing the query,
    it fetches all the rows of the query result and returns them. The function handles any errors that might
    occur during the process and ensures that the database connection is closed before returning the result.

    Args:
        query (str): The SQL query to be executed.
        database (str): The name of the database to connect to.
        host (str): The host name or IP address of the MySQL server.
        port (int): The port number to use for the connection.
        user (str): The username to use for the database connection.

    Returns:
        tuple: A tuple of tuples containing the rows returned by the query execution.

    Note:
        This function does not handle database password authentication. Ensure that the provided user
        has the necessary permissions and that the database is configured to allow password-less connections
        from the given host.
    """
    try:
        conn = pymysql.connect(
            host=host, user=user, port=port, database=database.strip()
        )

        cursor = conn.cursor()
        cursor.execute(query)
        info = cursor.fetchall()

    except pymysql.Error as err:
        print(err)

    cursor.close()
    conn.close()
    try: 
        return info
    except UnboundLocalError:
        print(f"\nNothing returned for SQL query: {query}\n")
        sys.exit()
    

def get_assembly_accessions(bioproject_id: str, only_haploid: bool = False) -> List[str]:
    """
    Fetches assembly accessions from a given NCBI BioProject ID.

    Args:
    bioproject_id (str): The NCBI BioProject ID.
    only_haploid (bool): If True, fetch only haploid assemblies. Defaults to False.

    Returns:
    list: A list of assembly accessions.
    """
    base_url = config["urls"]["datasets"]["bioproject"]
    next_page_token = None
    assembly_accessions = {}
    
    while True:
        url = f"{base_url}/{bioproject_id}/dataset_report"
        if next_page_token:
            url += f"?page_token={next_page_token}"
            
        try:
            response = requests.get(url)
            response.raise_for_status()  # This will raise an exception for 4XX and 5XX errors
            data = response.json()
            
            assemblies = data.get('reports', [])
            for assembly in assemblies:
                assembly_info = assembly.get('assembly_info', {})
                if only_haploid and assembly_info.get('assembly_type') != 'haploid':
                    continue
                assembly_accession = assembly.get('accession')
                taxon_id = assembly.get('organism').get('tax_id')
                if assembly_accession:
                    assembly_accessions[assembly_accession] = {"taxon_id": taxon_id}
                    
                # Check for the next page token or equivalent
            next_page_token = data.get('next_page_token')
            if not next_page_token:
                break
            
        except requests.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
            break  # Or handle it in some other way, e.g., retry
        except Exception as err:
            print(f"An error occurred: {err}")
            break  # Or handle it differently, maybe a retry logic
        
    return assembly_accessions
    
def get_ensembl_live(bioproject_accessions_taxon: Dict[str, Dict[str, int]]) -> Dict[str, Dict[str, str]]:
    """
    Retrieves live Ensembl database names for a set of bioproject accessions.

    This function constructs and executes a SQL query to fetch the latest Ensembl database names (denoted as 'dbname')
    associated with a list of assembly accessions. The query selects databases where the assembly accession matches
    those provided and filters for 'core' databases in the latest data release. The information is then used to
    update and return a dictionary mapping each accession to its corresponding taxonomy information and live Ensembl 
    'dbname'.

    Args:
        bioproject_accessions_taxon (dict): A dictionary where keys are assembly accessions and values are dictionaries 
                                            containing taxonomy information for those accessions.

    Returns:
        dict: A dictionary where each key is an assembly accession, and each value is a dictionary with the original 
              taxonomy information plus a 'dbname' key containing the name of the live Ensembl database for that accession.

    Requires:
        - A working connection to the Ensembl metadata database specified by the configuration in `config`.
        - The `mysql_fetch_data` function to execute the query and fetch data.

    Example:
        >>> bioproject_accessions_taxon = {'GCA_123456.1': {'taxon_id': 9606}}
        >>> live_annotations = get_ensembl_live(bioproject_accessions_taxon)
        >>> print(live_annotations)
        {'GCA_123456.1': {'taxon_id': 9606, 'dbname': 'homo_sapiens_core_104_38'}}

    Note:
        The configuration for the database connection (`config`) must be defined externally with keys for 
        'db_host', 'db_port', and 'db_user' within a nested 'server_details'->'meta' structure.
    """
    accessions = list(bioproject_accessions_taxon.keys())
    formatted_accessions = ",".join(f'"{item}"' for item in accessions)
    data_query = (
        "SELECT assembly.assembly_accession, dbname FROM assembly JOIN genome USING (assembly_id) JOIN genome_database USING (genome_id) WHERE genome.data_release_id=(SELECT MAX(data_release_id) FROM genome) AND assembly.assembly_accession in (" + formatted_accessions + ") AND dbname like '%core%';"
    )
    data_fetch = mysql_fetch_data(
        data_query,
        config["server_details"]["meta"]["rapid"]["db_name"],
        config["server_details"]["meta"]["rapid"]["db_host"],
        config["server_details"]["meta"]["rapid"]["db_port"],
        config["server_details"]["meta"]["rapid"]["db_user"],
    )
    live_annotations = {}
    for tuple in data_fetch:
        accession = tuple[0]
        live_info = {"dbname" : tuple[1]}
        live_annotations[accession] = bioproject_accessions_taxon[accession]
        live_annotations[accession].update(live_info)
    return(live_annotations)

def get_taxonomy_info(
            live_annotations: Dict[str, Dict[str, str]],
            bioproject_accessions_taxon: Dict[str, Dict[str, int]],
            rank: str
        ) -> Dict[str, Dict[str, str]]:
    """
    Updates the live_annotations dictionary with taxonomy information for a specified rank.
    
    This function iterates over each accession in live_annotations, retrieves the corresponding taxon ID,
    and makes a request to the NCBI Datasets API to fetch taxonomy information. It then updates the 
    live_annotations dictionary with the name of the specified rank (e.g., 'order', 'family') for each accession.
    
    Args:
        live_annotations (dict): A dictionary where each key is an accession number and its value is another 
                                  dictionary with various annotation details. This dictionary is updated in-place.
        bioproject_accessions_taxon (dict): A dictionary mapping accession numbers to their respective taxon 
                                            information, including taxon IDs.
        rank (str): The taxonomic rank for which the name should be retrieved and added to live_annotations 
                    (e.g., 'order', 'family').

    Returns:
        dict: The updated live_annotations dictionary with the added taxonomy information for the specified rank.

    Raises:
        requests.HTTPError: If an HTTP error occurs during the API request.
        Exception: If any other error occurs during the function's operation.

    Note:
        The function updates live_annotations in-place and also returns it for convenience. Each accession in 
        live_annotations is updated with a new key-value pair, where the key is the specified rank and the value 
        is the name of that rank from the taxonomy data.

    Example:
        >>> live_annotations = {'accession1': {'some_annotation': 'value'}}
        >>> bioproject_accessions_taxon = {'accession1': {'taxon_id': '12345'}}
        >>> rank = 'order'
        >>> updated_annotations = get_taxonomy_info(live_annotations, bioproject_accessions_taxon, rank)
        >>> print(updated_annotations)
        {'accession1': {'some_annotation': 'value', 'order': 'SomeOrderName'}}
    """
    for accession in live_annotations:
        taxon_id = bioproject_accessions_taxon[accession]["taxon_id"]
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{taxon_id}/dataset_report"
        
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
            data = response.json()
            taxonomies = data.get('reports', [])
            for taxonomy in taxonomies:
                rank_name = taxonomy.get('taxonomy', {}).get('classification', {}).get(rank, {}).get('name')
                taxonomy_info = {rank : rank_name}
                live_annotations[accession].update(taxonomy_info)
                
        except requests.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
        except Exception as err:
            print(f"An error occurred: {err}")
           
    return(live_annotations)
        
def write_report(live_annotations, report_file, rank):
    with open(Path(report_file), "w", encoding="utf8") as report_file:
        for accession in live_annotations:
            try:
                report_file.write(accession + '\t' + live_annotations[accession]["dbname"] + '\t' + live_annotations[accession][rank] + '\n')
            except TypeError:
                report_file.write(accession + '\t' + live_annotations[accession]["dbname"] + '\t unknown \n')
                print(live_annotations[accession])
def main():
    """
    Main function to handle command-line arguments and output the result.
    """
    parser = argparse.ArgumentParser(description='Fetch assembly accessions from NCBI BioProject and report the number of corresponding annotations in rapid.ensembl.org.')
    parser.add_argument('bioproject_id', type=str, help='NCBI BioProject ID')
    parser.add_argument('--haploid', action='store_true', help='Fetch only haploid assemblies')
    parser.add_argument('--report_file', type=str, help='Where to write report to', default='./report_file.csv')
    parser.add_argument('--rank', type=str, help='Taxonomic rank to classify', default='order')
    
    args = parser.parse_args()

    bioproject_accessions_taxon = get_assembly_accessions(args.bioproject_id, args.haploid)
    live_annotations = get_ensembl_live(bioproject_accessions_taxon)
    live_annotations_classified = get_taxonomy_info(live_annotations, bioproject_accessions_taxon, args.rank)
    unique_taxon_ids = {details['taxon_id'] for details in live_annotations.values()}
    print(f"Found {len(bioproject_accessions_taxon)} assemblies under BioProject ID {args.bioproject_id}")
    print(f"Found {len(live_annotations)} annotations in rapid.ensembl.org for {len(unique_taxon_ids)} unique species")

    rank_values = [details[args.rank] for details in live_annotations.values()]
    rank_counts = Counter(rank_values)
    print("\nBreakdown:")
    print(rank_counts)
    write_report(live_annotations, args.report_file, args.rank)
    
if __name__ == "__main__":
    main()

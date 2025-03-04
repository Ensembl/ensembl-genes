import argparse
import requests
import pymysql
import json
import logging
from pathlib import Path
from collections import Counter
from typing import List, Tuple, Any, Dict

# Enable logging for debugging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Load configuration
with open("./bioproject_tracking_config.json", "r") as f:
    config = json.load(f)


def mysql_fetch_data(query: str, params: Tuple = ()) -> List[Tuple[Any, ...]]:
    """
    Executes a SQL query with optional parameters to fetch results from the MySQL database.

    This function connects to a MySQL database using the credentials specified in the
    global `config` variable. It then executes the provided SQL query (along with any
    parameters) and returns all of the fetched results as a list of tuples. If an error
    occurs during the database operation, the error is logged and an empty list is returned.

    Args:
        query (str): The SQL query to be executed.
        params (Tuple, optional): A tuple of parameters to be used with the query.
            Defaults to an empty tuple ().

    Returns:
        List[Tuple[Any, ...]]: A list of tuples containing the rows returned by the query.
        If an error occurs or no data is found, an empty list is returned.
    """
    try:
        conn = pymysql.connect(
            host=config["server_details"]["meta"]["beta"]["db_host"],
            user=config["server_details"]["meta"]["beta"]["db_user"],
            port=config["server_details"]["meta"]["beta"]["db_port"],
            database=config["server_details"]["meta"]["beta"]["db_name"],
        )
        with conn.cursor() as cursor:
            cursor.execute(query, params)
            result = cursor.fetchall()
        conn.close()
        return result
    except pymysql.Error as err:
        logging.error(f"MySQL Error: {err}")
        return []


def get_assembly_accessions(query_id: str, query_type: str, only_haploid: bool = False) -> Dict[str, Dict[str, int]]:
    """
    Fetches assembly accessions from NCBI API.

    This function queries the NCBI Datasets API to retrieve assembly accessions for the
    provided `query_id`. The `query_type` indicates whether the ID belongs to a
    "bioproject" or a "taxon". The function supports pagination by automatically requesting
    subsequent pages using the `next_page_token`. If `only_haploid` is True, only assemblies
    of type "haploid" will be included.

    Args:
        query_id (str): The BioProject ID or Taxonomy ID to query.
        query_type (str): The type of the query ("bioproject" or "taxon").
        only_haploid (bool, optional): If True, only haploid assemblies are fetched.
            Defaults to False.

    Raises:
        ValueError: If an invalid `query_type` is passed.

    Returns:
        Dict[str, Dict[str, int]]: A dictionary where keys are assembly accession strings,
        and values are dictionaries with relevant information (currently just "taxon_id").
        Example:
            {
                "GCA_000001405.39": {"taxon_id": 9606},
                ...
            }
    """
    base_url = config["urls"]["datasets"].get(query_type)
    if not base_url:
        raise ValueError("Invalid query_type. Must be 'bioproject' or 'taxon'.")

    assembly_accessions = {}
    page_size = 5000
    next_page_token = None

    while True:
        url = f"{base_url}/{query_id}/dataset_report?page_size={page_size}"
        if next_page_token:
            url += f"&page_token={next_page_token}"

        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()

            for assembly in data.get("reports", []):
                assembly_info = assembly.get("assembly_info", {})
                if only_haploid and assembly_info.get("assembly_type") != "haploid":
                    continue

                accession = assembly.get("accession")
                taxon_id = assembly.get("organism", {}).get("tax_id")

                if accession:
                    assembly_accessions[assembly["accession"]] = {"taxon_id": taxon_id}

            next_page_token = data.get("next_page_token")
            if not next_page_token:
                break
        except requests.RequestException as e:
            logging.error(f"API error: {e}")
            break

    return assembly_accessions


def get_ensembl_live(accessions_taxon: Dict[str, Dict[str, int]]) -> Dict[str, Dict[str, str]]:
    """
    Fetches Ensembl live database names for a list of assembly accessions.

    This function queries a local MySQL database (using `mysql_fetch_data`) to find
    Ensembl release information that corresponds to each given assembly accession.
    It then merges the returned data (genome UUID and database name) with the existing
    taxon information, storing them under the keys "guuid" and "dbname".

    Args:
        accessions_taxon (Dict[str, Dict[str, int]]): A dictionary where keys are
            assembly accessions, and values contain taxonomic information such
            as "taxon_id". Example:
            {
                "GCA_000001405.39": {"taxon_id": 9606},
                ...
            }

    Returns:
        Dict[str, Dict[str, str]]: A dictionary with assembly accessions as keys and
        updated values that contain "guuid", "dbname", and existing taxon info. Example:
            {
                "GCA_000001405.39": {
                    "taxon_id": 9606,
                    "guuid": "some-uuid",
                    "dbname": "ensembl_core"
                },
                ...
            }
        If no accessions are found or an error occurs, an empty dictionary is returned.
    """
    accessions = list(accessions_taxon.keys())
    if not accessions:
        return {}

    query = """
        SELECT assembly.accession, genome.genome_uuid, dataset_source.name AS database_name
        FROM genome
        JOIN assembly ON genome.assembly_id = assembly.assembly_id
        JOIN genome_dataset ON genome.genome_id = genome_dataset.genome_id
        JOIN dataset ON genome_dataset.dataset_id = dataset.dataset_id
        JOIN dataset_source ON dataset.dataset_source_id = dataset_source.dataset_source_id
        JOIN genome_release ON genome.genome_id = genome_release.genome_id
        JOIN ensembl_release ON genome_release.release_id = ensembl_release.release_id
        WHERE dataset.name = 'assembly'
        AND assembly.accession IN ({})
    """.format(", ".join(["%s"] * len(accessions)))

    data_fetch = mysql_fetch_data(query, tuple(accessions))

    live_annotations = {}
    for accession, guuid, dbname in data_fetch:
        live_annotations[accession] = accessions_taxon[accession]
        live_annotations[accession].update({"guuid": guuid, "dbname": dbname})

    return live_annotations


def get_taxonomy_info(
    live_annotations: Dict[str, Dict[str, str]],
    accessions_taxon: Dict[str, Dict[str, int]],
    rank: str
) -> Dict[str, Dict[str, str]]:
    """
    Fetches taxonomy information for a given rank from the NCBI Datasets API.

    For each accession in `live_annotations`, this function retrieves its `taxon_id`
    and queries the NCBI taxonomy endpoint to get the classification at the specified rank.
    The retrieved rank name is then added to the corresponding annotation data.

    Args:
        live_annotations (Dict[str, Dict[str, str]]): A dictionary where each key is
            an assembly accession and its value contains annotation details (including
            at least "taxon_id").
        accessions_taxon (Dict[str, Dict[str, int]]): A dictionary mapping accessions
            to basic taxon information, such as {"GCA_000001405.39": {"taxon_id": 9606}}.
        rank (str): The taxonomic rank to retrieve (e.g., "order", "class", "phylum").

    Returns:
        Dict[str, Dict[str, str]]: The updated `live_annotations` dictionary with an
        additional key for the requested rank. For example:
            {
                "GCA_000001405.39": {
                    "taxon_id": 9606,
                    "guuid": "some-uuid",
                    "dbname": "ensembl_core",
                    "order": "Primates"
                },
                ...
            }
        If the API query fails or the rank is not found, the value for that rank is not added
        (or set to None/unknown).
    """
    for accession in live_annotations:
        taxon_id = accessions_taxon[accession]["taxon_id"]
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{taxon_id}/dataset_report"

        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()

            taxonomies = data.get('reports', [])
            for taxonomy in taxonomies:
                rank_name = taxonomy.get('taxonomy', {}).get('classification', {}).get(rank, {}).get('name')
                taxonomy_info = {rank: rank_name}
                live_annotations[accession].update(taxonomy_info)

        except requests.HTTPError as http_err:
            logging.error(f"HTTP error occurred: {http_err}")
        except Exception as err:
            logging.error(f"An error occurred: {err}")

    return live_annotations


def add_ftp(live_annotations: Dict[str, Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    """
    Adds FTP links to the live annotations dictionary.

    This function queries the local MySQL database (using `mysql_fetch_data`) to retrieve
    `genebuild_date` and `scientific_name` for each genome UUID. Using this information, it
    constructs an FTP link for the relevant Ensembl data and adds it to the `live_annotations`
    under the key "ftp".

    Args:
        live_annotations (Dict[str, Dict[str, str]]): A dictionary where each key is
            an assembly accession, and its value contains annotation details (including
            at least "guuid").

    Returns:
        Dict[str, Dict[str, str]]: The updated `live_annotations` with the "ftp" key.
        Example:
            {
                "GCA_000001405.39": {
                    "taxon_id": 9606,
                    "guuid": "some-uuid",
                    "dbname": "ensembl_core",
                    "ftp": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/..."
                },
                ...
            }
        If no matching information is found for a UUID, the "ftp" key is not added.
    """
    for accession in live_annotations:
        guuid = live_annotations[accession]["guuid"]

        query = """
            SELECT genome.genebuild_date, organism.scientific_name
            FROM genome
            JOIN organism USING(organism_id)
            WHERE genome.genome_uuid = %s
        """
        data_fetch = mysql_fetch_data(query, (guuid,))

        if data_fetch:
            date, scientific_name = data_fetch[0]
            date = date.replace("-", "_")
            scientific_name = scientific_name.replace(" ", "_")

            # Construct FTP link
            ftp_link = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{scientific_name}/{accession}/ensembl/geneset/{date}/"
            live_annotations[accession]["ftp"] = ftp_link

    return live_annotations


def write_report(
    live_annotations: Dict[str, Dict[str, str]],
    report_file: str,
    rank: str,
    include_ftp: bool = False
):
    """
    Writes the report file with annotations.

    This function takes the final dictionary of `live_annotations`, iterates through each
    accession, and writes a row of annotation data to the specified `report_file`.
    By default, the columns include:
      - Accession
      - Genome UUID (guuid)
      - Database name (dbname)
      - Taxonomic rank (e.g., order)
    If `include_ftp` is True, it also appends the constructed FTP link.

    Args:
        live_annotations (Dict[str, Dict[str, str]]): A dictionary containing annotation
            information (keys are accessions, values include guuid, dbname, and rank info).
        report_file (str): Path to the output file where the report should be written.
        rank (str): The taxonomic rank to include in the report (e.g., "order").
        include_ftp (bool, optional): If True, the FTP link is appended to each row.
            Defaults to False.
    """
    report_path = Path(report_file)
    with open(report_path, "w", encoding="utf-8") as file:
        for accession, details in live_annotations.items():
            # Convert each value to string and replace None with a fallback value
            row = [
                str(accession),
                str(details.get("guuid") or "unknown"),
                str(details.get("dbname") or "unknown"),
                str(details.get(rank) or "unknown"),
            ]
            if include_ftp:
                row.append(str(details.get("ftp") or "N/A"))
                
            file.write("\t".join(row) + "\n")

    logging.info(f"Report written to {report_path.resolve()}.")
                
                
def main():
    """
    Main entry point of the script.

    Parses command-line arguments to determine if a BioProject ID or Taxon ID is provided,
    whether to filter only haploid assemblies, the desired output file path, the taxonomic
    rank to retrieve, and whether to include FTP links. It then:

      1. Fetches assembly accessions from the NCBI API based on the provided ID.
      2. Fetches corresponding Ensembl live database information from the local MySQL database.
      3. Retrieves taxonomy information for the specified rank from the NCBI API.
      4. Optionally adds FTP links if requested.
      5. Writes all collected data into a tab-separated report file.

    Usage Example:
        python script_name.py --bioproject_id PRJNA12345 --haploid --rank order --ftp

    Returns:
        None
    """
    parser = argparse.ArgumentParser(description="Fetch assembly data and report annotations.")
    parser.add_argument("--bioproject_id", type=str, help="NCBI BioProject ID")
    parser.add_argument("--taxon_id", type=str, help="Taxonomy ID")
    parser.add_argument("--haploid", action="store_true", help="Fetch only haploid assemblies")
    parser.add_argument("--report_file", type=str, default="./report_file.csv")
    parser.add_argument("--rank", type=str, default="order")
    parser.add_argument("--ftp", action="store_true", help="Include FTP links in the report")

    args = parser.parse_args()

    # Determine which ID is provided and which query type is relevant
    accessions_taxon = get_assembly_accessions(
        args.bioproject_id or args.taxon_id,
        "bioproject" if args.bioproject_id else "taxon",
        args.haploid
    )

    # Fetch additional annotations
    live_annotations = get_ensembl_live(accessions_taxon)
    live_annotations_classified = get_taxonomy_info(live_annotations, accessions_taxon, args.rank)

    # Optionally add FTP links
    if args.ftp:
        live_annotations = add_ftp(live_annotations)

    unique_taxon_ids = {details['taxon_id'] for details in live_annotations.values()}
    if (args.bioproject_id):
        print(f"Found {len(accessions_taxon)} assemblies under BioProject ID {args.bioproject_id}")
        
    elif (args.taxon_id):
        print(f"Found {len(accessions_taxon)} assemblies for taxon ID {args.taxon_id}")
    print(f"Found {len(live_annotations)} annotations in beta.ensembl.org for {len(unique_taxon_ids)} unique species")
        
    rank_values = [details[args.rank] for details in live_annotations.values()]
    rank_counts = Counter(rank_values)
    print("\nBreakdown:")
    print(rank_counts)
            
    # Write final report
    write_report(live_annotations, args.report_file, args.rank, include_ftp=args.ftp)


if __name__ == "__main__":
    main()

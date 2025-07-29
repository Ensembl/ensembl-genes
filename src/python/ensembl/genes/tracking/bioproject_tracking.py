import argparse
import requests
import pymysql
import json
import logging
import re
from pathlib import Path
from collections import Counter
from typing import List, Tuple, Any, Dict, Optional

# Enable logging for debugging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Load configuration
with open("./bioproject_tracking_config.json", "r") as f:
    config = json.load(f)


def mysql_fetch_data(
            query: str,
            params: Tuple = (),
            server_group: str = "meta",
            server_name: str = "beta",
            db_name: Optional[str] = None
        ) -> List[Tuple[Any, ...]]:
    """
    Executes a SQL query with optional parameters to fetch results from a specified MySQL server.

    The function supports dynamic selection of MySQL servers based on the `config` structure,
    which organizes server credentials under groups like "meta", "data", or "pre-release".

    Args:
        query (str): The SQL query to execute.
        params (Tuple, optional): Parameters to substitute into the SQL query. Default is ().
        server_group (str): The group of servers to use (e.g., "meta", "data", "pre-release").
        server_name (str): The specific server key within the group (e.g., "beta", "gb1").
        db_name (Optional[str]): Optional explicit database name. If not provided, uses the one
                                 from the config (if available).

    Returns:
        List[Tuple[Any, ...]]: Fetched query results. Returns an empty list on error or no data.
    """
    try:
        server_config = config["server_details"][server_group][server_name]
        
        connection = pymysql.connect(
            host=server_config["db_host"],
            user=server_config["db_user"],
            port=server_config["db_port"],
            database=db_name or server_config.get("db_name", ""),  # allow empty db if not needed
        )
        with connection.cursor() as cursor:
            cursor.execute(query, params)
            result = cursor.fetchall()
            connection.close()
            return result
        
    except KeyError as key_err:
        logging.error(f"Invalid server group or name in config: {key_err}")
    except pymysql.Error as sql_err:
        logging.error(f"MySQL Error: {sql_err}")

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

    # Compute accessions that were not annotated
    missing_annotations = [acc for acc in accessions if acc not in live_annotations]
    return live_annotations, missing_annotations


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


def add_ftp(
        annotations: Dict[str, Dict[str, str]],
        release_type: str = "live"
) -> Dict[str, Dict[str, str]]:
    for accession, annotation in annotations.items():
        if release_type == "pre":
            dbname = annotation.get("dbname")
            if not dbname:
                continue  # Can't query without db name
            
            # Query the pre-release DB itself for scientific name
            query = """
            SELECT meta_value 
            FROM meta 
            WHERE meta_key = 'organism.scientific_name'
            """
            result = mysql_fetch_data(
                query,
                (),
                server_group="pre-release",
                server_name="gb1",
                db_name=dbname
            )
                
            if not result:
                continue
            
            scientific_name = result[0][0].replace(" ", "_")
            ftp_link = f"https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/{scientific_name}/{accession}"
            annotation["ftp"] = ftp_link
                
        elif release_type == "live":
            guuid = annotation.get("guuid")
            if not guuid or guuid == "unknown":
                continue
            
            query = """
            SELECT genome.genebuild_date, organism.scientific_name
            FROM genome
            JOIN organism USING(organism_id)
            WHERE genome.genome_uuid = %s
            """
            data_fetch = mysql_fetch_data(query, (guuid,))
            
            if not data_fetch:
                continue
            
            date, scientific_name = data_fetch[0]
            date = date.replace("-", "_")
            scientific_name = scientific_name.replace(" ", "_")
            
            ftp_link = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{scientific_name}/{accession}/ensembl/geneset/{date}/"
            annotation["ftp"] = ftp_link
                
    return annotations

def get_pre_release(missing_annotations: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Checks for the existence of pre-release Ensembl database schemas for a list of missing accessions.

    This function takes a list of GenBank/RefSeq assembly accessions that were not found in the
    Ensembl live database lookup (e.g., from `get_ensembl_live`). It reformats each accession
    to match the expected pre-release Ensembl database schema naming convention (e.g.,
    'GCA_000001405.39' → 'gca000001405v39') and queries the MySQL server's 
    `information_schema.schemata` table to check if a database with that name exists.

    If such a schema is found, the function records the corresponding accession and
    database name in the returned dictionary.

    Args:
        missing_annotations (List[str]): A list of assembly accessions (e.g., 'GCA_000001405.39')
            for which no live annotation database was found.

    Returns:
        Dict[str, Dict[str, str]]: A dictionary mapping each accession with a found pre-release
        database to a nested dictionary containing:
            - "dbname": The name of the pre-release schema (e.g., 'gca000001405v39').

        Example:
            {
                "GCA_000001405.39": {
                    "dbname": "gca000001405v39"
                },
                ...
            }

    Notes:
        - Only accessions matching the pattern 'GCA/GCF_XXXXXXXXX.XX' are processed.
        - If no matching schema is found for an accession, it is omitted from the result.
        - The function uses `mysql_fetch_data` to query the MySQL server's schema metadata.
    """
    pre_release_annotations = {}

    for accession in missing_annotations:
        match = re.match(r"(GCA|GCF)_(\d+)\.(\d+)", accession)
        if not match:
            continue  # skip malformed accession
        
        prefix, number, version = match.groups()
        lookup_string = f"%{prefix.lower()}{number}v{version}%"
        
        query = """
        SELECT SCHEMA_NAME
        FROM information_schema.schemata
        WHERE SCHEMA_NAME like %s
        """
        result = mysql_fetch_data(
            query,
            (lookup_string,),
            server_group="pre-release",
            server_name="gb1"
        )
            
        if result:
            pre_release_annotations[accession] = {
                "guiid": 'unknown',
                "dbname": result[0][0]  # extract the SCHEMA_NAME from the tuple
            }
    return pre_release_annotations
            
def write_report(
    live_annotations: Dict[str, Dict[str, str]],
    report_file: str,
    rank: str,
    include_ftp: bool = False
):
    """
    Writes the report file with annotations.

    This function takes the final dictionary of live_annotations and writes each 
    annotation to a tab-separated report file specified by the report_file path.
    Each row in the file contains the following columns:
      - Accession
      - Genome UUID (guuid)
      - Database name (dbname)
      - The specified taxonomic rank (e.g., "order", "class", etc.)
    Optionally, if include_ftp is True, an additional column is added for the FTP link.

    Each value is converted to a string. If a value is missing (i.e., is None or falsy),
    it will be replaced by a default value:
      - For guuid, dbname, and the rank column, the default is "unknown".
      - For the FTP link, the default is "N/A".

    Parameters
    ----------
    live_annotations : Dict[str, Dict[str, str]]
        A dictionary where each key is an assembly accession and its value is a dictionary
        containing annotation details (including at least "guuid", "dbname", and the taxonomic rank).
    report_file : str
        The file path where the report should be written.
    rank : str
        The taxonomic rank to be included in the report (e.g., "order", "class").
    include_ftp : bool, optional
        If True, an additional column with FTP links is added to each row (default is False).

    Returns
    -------
    None
        The function writes the report directly to the specified file and logs the output.
    """
    report_path = Path(report_file)
    lines_written = 0

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
            lines_written += 1

    logging.info(f"Report written to {report_path.resolve()} with {lines_written} lines.")



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
    parser.add_argument("--pre_release", action="store_true", help="Include list of pre-release databases in the report")

    args = parser.parse_args()

    # Determine which ID is provided and which query type is relevant
    accessions_taxon = get_assembly_accessions(
        args.bioproject_id or args.taxon_id,
        "bioproject" if args.bioproject_id else "taxon",
        args.haploid
    )

    # Fetch annotations
    live_annotations, missing_annotations = get_ensembl_live(accessions_taxon)
    # Ensure taxonomic rank is added before any logic that depends on it
    get_taxonomy_info(live_annotations, accessions_taxon, args.rank)
    
    # Optionally add FTP links
    if args.ftp:
        live_annotations = add_ftp(live_annotations, 'live')

    unique_taxon_ids = {details['taxon_id'] for details in live_annotations.values()}

    # Optionally add Pre-release data
    if args.pre_release:
        pre_release_annotations = get_pre_release(missing_annotations)

        for accession in pre_release_annotations:
            if accession in accessions_taxon:
                pre_release_annotations[accession]["taxon_id"] = accessions_taxon[accession]["taxon_id"]
        get_taxonomy_info(pre_release_annotations, accessions_taxon, args.rank)

        if args.ftp:
            pre_release_annotations = add_ftp(pre_release_annotations, 'pre')

        all_annotations = {**live_annotations, **pre_release_annotations}

    else:
        all_annotations = live_annotations
        
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
    write_report(all_annotations, args.report_file, args.rank, include_ftp=args.ftp)
    


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# pylint:disable=line-too-long, logging-fstring-interpolation
"""Script to track assemblies from NCBI BioProject/Taxon IDs and their presence in Ensembl.
Generates a report of assembly accessions, Ensembl database info, taxonomy classification,
and FTP links if available.
"""

import argparse


import json
import logging
import re
from pathlib import Path
from collections import Counter
from typing import List, Sequence, Tuple, Any, Dict, Optional
import subprocess
import shutil
import requests
import pymysql

# -----------------------------------
# Logging
# -----------------------------------
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# -----------------------------------
# Config
# -----------------------------------
with open(  # pylint:disable=unspecified-encoding
    "./bioproject_tracking_config.json", "r"
) as f:
    config = json.load(f)


# -----------------------------------
# DB helper
# -----------------------------------
def mysql_fetch_data(
    query: str,
    params: Tuple = (),
    server_group: str = "meta",
    server_name: str = "beta",
    db_name: Optional[str] = None,
) -> Sequence[Tuple[Any, ...]]:
    """
    Executes a SQL query with optional parameters to fetch results from a specified MySQL server.

    Args:
        query (str): The SQL query to execute.
        params (Tuple, optional): Parameters to pass to the SQL query. Defaults to ().
        server_group (str, optional): The server group as defined in the config. Defaults to "meta".
        server_name (str, optional): The server name within the group as defined in the config. Defaults to "beta".
        db_name (Optional[str], optional): The database name to connect to. If None, uses the default from config. Defaults to None.
    
    Returns:
        Sequence[Tuple[Any, ...]]: The fetched rows from the query result.
    
    Raises:
        KeyError: If the server group or name is not found in the config.
        pymysql.Error: If there is an error connecting to the database or executing the query.
    """
    try:
        server_config = config["server_details"][server_group][server_name]
        connection = pymysql.connect(
            host=server_config["db_host"],
            user=server_config["db_user"],
            port=server_config["db_port"],
            database=db_name or server_config.get("db_name", ""),
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


# -----------------------------------
# NCBI accessions helper
# -----------------------------------
def get_assembly_accessions(  # pylint:disable=too-many-branches, too-many-statements, too-many-locals
    query_id: str, query_type: str, only_haploid: bool = False
) -> Dict[str, Dict[str, int]]:
    """
    Fetches assembly accessions from NCBI Datasets API based on BioProject or Taxon ID.
    
    Args:
        query_id (str): The BioProject or Taxon ID to query.
        query_type (str): The type of query, either "bioproject" or "taxon".
        only_haploid (bool): If True, only include haploid assemblies. Defaults to False.
        
    Returns:
        Dict[str, Dict[str, int]]: A dictionary mapping assembly accessions to their taxon IDs.
    """

    def _parse_reports(lines: List[str]) -> Dict[str, Dict[str, int]]:
        accs: Dict[str, Dict[str, int]] = {}
        for line in lines:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            report = obj.get("report", obj)
            assembly_info = report.get("assembly_info") or {}
            if only_haploid and assembly_info.get("assembly_type") != "haploid":
                continue
            accession = report.get("accession")
            taxon_id = (report.get("organism") or {}).get("tax_id")
            if not isinstance(taxon_id, int):
                continue
            if accession:
                accs[accession] = {"taxon_id": taxon_id}
        return accs

    # --- 1) Try CLI (gets ALL versions) ---
    datasets_cli = shutil.which("datasets")
    if datasets_cli:
        try:
            if query_type == "bioproject":
                cmd = [
                    datasets_cli,
                    "summary",
                    "genome",
                    "accession",
                    str(query_id),
                    "--assembly-version",
                    "all",
                    "--as-json-lines",
                ]
            elif query_type == "taxon":
                cmd = [
                    datasets_cli,
                    "summary",
                    "genome",
                    "taxon",
                    str(query_id),
                    "--assembly-version",
                    "all",
                    "--as-json-lines",
                ]
            else:
                logging.error(f"Invalid query_type '{query_type}'")
                return {}

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            parsed = _parse_reports(result.stdout.splitlines())
            if parsed:  # pylint:disable=no-else-return
                logging.info(
                    f"[Datasets CLI] Retrieved {len(parsed)} assemblies (all versions)."
                )
                return parsed
            else:
                logging.warning(
                    "[Datasets CLI] Parsed zero rows from stdout; falling back to API."
                )
        except subprocess.CalledProcessError as e:
            logging.warning(
                f"[Datasets CLI] Exit {e.returncode}. STDERR:\n{e.stderr}\nFalling back to API."
            )
        except Exception as e:  # pylint:disable=broad-exception-caught
            logging.warning(
                f"[Datasets CLI] Unexpected error: {e}. Falling back to API."
            )

    # --- 2) Fallback: API (latest only) ---
    try:
        base_url = config["urls"]["datasets"].get(query_type)
    except Exception:  # pylint:disable=broad-exception-caught
        base_url = None

    if not base_url:
        logging.error(
            "[Datasets API] Base URL not found in config; returning empty set."
        )
        return {}

    assembly_accessions: Dict[str, Dict[str, int]] = {}
    page_size = 5000
    next_page_token = None

    while True:
        url = f"{base_url}/{query_id}/dataset_report?page_size={page_size}"
        if next_page_token:
            url += f"&page_token={next_page_token}"
        try:
            response = requests.get(url)  # pylint:disable=missing-timeout
            response.raise_for_status()
            data = response.json()
            for assembly in data.get("reports", []):
                assembly_info = assembly.get("assembly_info", {}) or {}
                if only_haploid and assembly_info.get("assembly_type") != "haploid":
                    continue
                accession = assembly.get("accession")
                taxon_id = (assembly.get("organism") or {}).get("tax_id")
                if not isinstance(taxon_id, int):
                    continue
                if accession:
                    assembly_accessions[accession] = {"taxon_id": taxon_id}
            next_page_token = data.get("next_page_token")
            if not next_page_token:
                break
        except requests.RequestException as e:
            logging.error(f"[Datasets API] Error: {e}")
            break

    logging.info(
        f"[Datasets API] Retrieved {len(assembly_accessions)} assemblies (latest only)."
    )
    return assembly_accessions


# -----------------------------------
# Ensembl live lookup (collect all + best)
# -----------------------------------
_CORE_VER_RE = re.compile(r"_core_(\d+)_", re.IGNORECASE)


def _score_dbname(dbname: str) -> Tuple[int, int, str]:
    """
    Higher tuple is better.
    1) Higher core version
    2) Prefer non-'cm'
    3) Lexicographic tiebreak
    
    Args:
        dbname (str): The database name to score.
        
    Returns:
        Tuple[int, int, str]: A tuple representing the score.
    """
    m = _CORE_VER_RE.search(dbname or "")
    core_ver = int(m.group(1)) if m else -1
    non_cm = 1 if "_cm_" not in (dbname or "").lower() else 0
    return (core_ver, non_cm, dbname or "")


def get_ensembl_live(
    accessions_taxon: Dict[str, Dict[str, int]],
) -> Tuple[Dict[str, Dict[str, Any]], List[str]]:
    """
    Fetches Ensembl live database annotations for a list of assembly accessions.
    
    Args:
        accessions_taxon (Dict[str, Dict[str, int]]): A dictionary mapping assembly accessions
            to basic taxon information, such as {"GCA_000001405.39": {"taxon_id": 9606}}.
            
    Returns:
        Tuple[Dict[str, Dict[str, Any]], List[str]]: A tuple containing:
        1) A dictionary where each key is an assembly accession and its value is another dictionary with:
           - taxon_id (from input)
           - matches: List[{"guuid","dbname"}]
           - guuid, dbname: best match (latest core; prefer non-cm)
        2) A list of accessions that were not found in the Ensembl live database.
    """
    accessions = list(accessions_taxon.keys())
    if not accessions:
        return {}, []
    placeholders = ", ".join(["%s"] * len(accessions))
    query = f"""
    SELECT DISTINCT
        assembly.accession,
        genome.genome_uuid,
        dataset_source.name AS database_name
    FROM genome
    JOIN assembly ON genome.assembly_id = assembly.assembly_id
    JOIN genome_dataset ON genome.genome_id = genome_dataset.genome_id
    JOIN dataset ON genome_dataset.dataset_id = dataset.dataset_id
    JOIN dataset_source ON dataset.dataset_source_id = dataset_source.dataset_source_id
    JOIN genome_release ON genome.genome_id = genome_release.genome_id
    JOIN ensembl_release ON genome_release.release_id = ensembl_release.release_id
    WHERE dataset.name = 'genebuild'
      AND assembly.accession IN ({placeholders})
    """

    rows = mysql_fetch_data(query, tuple(accessions))

    # Seed with taxon info so missing ones can still be tracked
    live_annotations: Dict[str, Dict[str, Any]] = {}

    for accession, guuid, dbname in rows:
        entry = live_annotations.setdefault(
            accession, dict(accessions_taxon.get(accession, {}))
        )
        entry.setdefault("matches", [])
        # Avoid duplicates across joins if any
        if not any(
            m["guuid"] == guuid and m["dbname"] == dbname for m in entry["matches"]
        ):
            entry["matches"].append({"guuid": guuid, "dbname": dbname})

    # Choose best match for backward compatibility fields
    for accession, entry in live_annotations.items():
        matches = entry.get("matches", [])
        if matches:
            best = max(matches, key=lambda m: _score_dbname(m["dbname"]))
            entry["guuid"] = best["guuid"]
            entry["dbname"] = best["dbname"]

    # Compute missing (accessions with no rows at all)
    matched_accessions = set(live_annotations.keys())
    missing_annotations = [acc for acc in accessions if acc not in matched_accessions]

    return live_annotations, missing_annotations


# -----------------------------------
# Taxonomy info
# -----------------------------------
def get_taxonomy_info(
    live_annotations: Dict[str, Dict[str, Any]],
    accessions_taxon: Dict[str, Dict[str, int]],
    rank: str,
) -> Dict[str, Dict[str, Any]]:
    """
    Fetches taxonomy information for a given rank from the NCBI Datasets API.
    For each accession in `live_annotations`, this function retrieves its `taxon_id`\
    and queries the NCBI taxonomy endpoint to get the classification at the specified rank.
    The retrieved rank name is then added to the corresponding annotation data.
    
    Args:
        live_annotations (Dict[str, Dict[str, str]]): A dictionary where each key is \
            an assembly accession and its value contains annotation details (including \
            at least "taxon_id").
        accessions_taxon (Dict[str, Dict[str, int]]): A dictionary mapping accessions \
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
            response = requests.get(url)  # pylint:disable=missing-timeout
            response.raise_for_status()
            data = response.json()

            taxonomies = data.get("reports", [])
            for taxonomy in taxonomies:
                rank_name = (
                    taxonomy.get("taxonomy", {})
                    .get("classification", {})
                    .get(rank, {})
                    .get("name")
                )
                taxonomy_info = {rank: rank_name}
                live_annotations[accession].update(taxonomy_info)

        except requests.HTTPError as http_err:
            logging.error(f"HTTP error occurred: {http_err}")
        except Exception as err:  # pylint:disable=broad-exception-caught
            logging.error(f"An error occurred: {err}")

    return live_annotations


# -----------------------------------
# Add FTP (per match + best)
# -----------------------------------
def add_ftp(  # pylint:disable=too-many-locals, too-many-branches
    annotations: Dict[str, Dict[str, Any]], release_type: str = "live"
) -> Dict[str, Dict[str, Any]]:
    """
    Adds FTP links to the annotations dictionary based on the release type.
    
    Args:
        annotations (Dict[str, Dict[str, Any]]): A dictionary where each key is \
            an assembly accession and its value contains annotation details, \
            including matches with database names or genome UUIDs.
        release_type (str): The type of release to determine FTP link construction.\
            Can be either "live" or "pre".
            
    Returns:
        Dict[str, Dict[str, Any]]: The updated annotations dictionary with FTP links \
        added to each match and the top-level annotation where applicable.
    """
    for accession, annotation in annotations.items():
        matches = annotation.get("matches") or []

        if release_type == "pre":
            dbnames = set(m.get("dbname") for m in matches if m.get("dbname"))
            if annotation.get("dbname"):
                dbnames.add(annotation["dbname"])

            for dbname in dbnames:
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
                    db_name=dbname,
                )

                if not result:
                    continue  # skip if no scientific name found

                scientific_name = result[0][0].replace(" ", "_").replace(".", "")
                ftp_link = f"https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/{scientific_name}/{accession}"

                # attach to matching entries
                for m in matches:
                    if m.get("dbname") == dbname:
                        m["ftp"] = ftp_link
                if annotation.get("dbname") == dbname:
                    annotation["ftp"] = ftp_link  # top-level

        elif release_type == "live":
            # Build map guuid->match for easy updates
            by_guuid = {m.get("guuid"): m for m in matches if m.get("guuid")}
            guuids = set(by_guuid.keys())
            if annotation.get("guuid"):
                guuids.add(annotation["guuid"])

            for guuid in guuids:
                if not guuid or guuid == "unknown":
                    continue
                query = """
                SELECT genome.genebuild_date, organism.scientific_name, dataset_attribute.value
                FROM genome
                JOIN organism USING(organism_id)
                JOIN genome_dataset USING(genome_id) 
                JOIN dataset_attribute USING(dataset_id)
                WHERE genome.genome_uuid = %s
                AND attribute_id=169
                """
                data_fetch = mysql_fetch_data(query, (guuid,))
                if not data_fetch:
                    continue

                date, scientific_name, source = data_fetch[0]
                date = str(date).replace("-", "_")
                scientific_name = scientific_name.replace(" ", "_").replace(".", "")
                source = source.lower()
                ftp_link = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{scientific_name}/{accession}/{source}/geneset/{date}/"

                if guuid in by_guuid:
                    by_guuid[guuid]["ftp"] = ftp_link
                    by_guuid[guuid]["date"] = date
                if annotation.get("guuid") == guuid:
                    annotation["ftp"] = ftp_link
                    annotation["date"] = date

    return annotations


def get_pre_release(missing_annotations: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Checks for the existence of pre-release Ensembl database schemas for a \
        list of missing accessions.
    This function takes a list of GenBank/RefSeq assembly accessions that were not found in the \
    Ensembl live database lookup (e.g., from `get_ensembl_live`). It reformats each accession \
    to match the expected pre-release Ensembl database schema naming convention (e.g.,\
    'GCA_000001405.39' → 'gca000001405v39') and queries the MySQL server's \
    `information_schema.schemata` table to check if a database with that name exists.
    If such a schema is found, the function records the corresponding accession and \
    database name in the returned dictionary.
    
    Args:
        missing_annotations (List[str]): A list of assembly accessions (e.g., 'GCA_000001405.39') \
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
    pre_release_annotations: Dict[str, Dict[str, str]] = {}

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
            query, (lookup_string,), server_group="pre-release", server_name="gb1"
        )

        if result:
            pre_release_annotations[accession] = {
                "guuid": "unknown",  # typo fixed from guiid -> guuid
                "dbname": result[0][0],
            }
    return pre_release_annotations


# -----------------------------------
# Report writer (explode when multiple matches)
# -----------------------------------
def write_report(
    live_annotations: Dict[str, Dict[str, Any]],
    report_file: str,
    rank: str,
    include_class: bool = False,
    include_ftp: bool = False,
):
    """
    Writes a tab-separated report. If an accession has multiple matches,
    produces one line per match; otherwise one line.
    
    Args:
        live_annotations (Dict[str, Dict[str, Any]]): A dictionary where each key is \
            an assembly accession and its value contains annotation details, \
            including matches with database names or genome UUIDs.
        report_file (str): The path to the output report file.
        rank (str): The taxonomic rank to include if `include_class` is True.
        include_class (bool): If True, includes the taxonomic classification at the specified rank.
        include_ftp (bool): If True, includes FTP links in the report.
        
    Returns:
        None
    """
    report_path = Path(report_file)
    lines_written = 0

    with open(report_path, "w", encoding="utf-8") as file:
        for accession, details in live_annotations.items():
            matches = details.get("matches")

            # If we have multiple, write one row per match
            rows = (
                matches
                if matches
                else [
                    {
                        "guuid": details.get("guuid"),
                        "dbname": details.get("dbname"),
                        "ftp": details.get("ftp"),
                        "date": details.get("date"),
                    }
                ]
            )

            for m in rows:
                row = [
                    str(accession),
                    str(m.get("guuid") or "unknown"),
                    str(m.get("dbname") or "unknown"),
                ]
                if include_class:
                    row.append(str(details.get(rank) or "unknown"))
                if include_ftp:
                    row.append(str(m.get("date") or details.get("date") or ""))
                    row.append(str(m.get("ftp") or details.get("ftp") or "N/A"))

                file.write("\t".join(row) + "\n")
                lines_written += 1

    logging.info(
        f"Report written to {report_path.resolve()} with {lines_written} lines."
    )


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
    parser = argparse.ArgumentParser(
        description="Fetch assembly data and report annotations."
    )
    parser.add_argument("--bioproject_id", type=str, help="NCBI BioProject ID")
    parser.add_argument("--taxon_id", type=str, help="Taxonomy ID")
    parser.add_argument(
        "--haploid", action="store_true", help="Fetch only haploid assemblies"
    )
    parser.add_argument("--report_file", type=str, default="./report_file.tsv")
    parser.add_argument(
        "--classification",
        action="store_true",
        help="Provide breakdown of taxonomic classification",
    )
    parser.add_argument("--rank", type=str, default="order")
    parser.add_argument(
        "--ftp", action="store_true", help="Include FTP links in the report"
    )
    parser.add_argument(
        "--pre_release",
        action="store_true",
        help="Include list of pre-release databases in the report",
    )

    args = parser.parse_args()

    # Determine which ID is provided and which query type is relevant
    if not args.bioproject_id and not args.taxon_id:
        parser.error("You must provide either --bioproject_id or --taxon_id")

    accessions_taxon = get_assembly_accessions(
        args.bioproject_id or args.taxon_id,
        "bioproject" if args.bioproject_id else "taxon",
        args.haploid,
    )

    # Fetch annotations (live)
    live_annotations, missing_annotations = get_ensembl_live(accessions_taxon)

    # Taxonomy classification (top-level only; rows inherit from details)
    if args.classification:
        get_taxonomy_info(live_annotations, accessions_taxon, args.rank)

    # Optionally add FTP links
    if args.ftp:
        live_annotations = add_ftp(live_annotations, "live")

    # Unique species count from live annotations
    unique_taxon_ids = {
        details["taxon_id"]
        for details in live_annotations.values()
        if "taxon_id" in details
    }

    # Optionally add Pre-release data
    if args.pre_release and missing_annotations:
        pre_release_annotations = get_pre_release(missing_annotations)

        for (
            accession
        ) in pre_release_annotations:  # pylint:disable=consider-using-dict-items
            if accession in accessions_taxon:
                pre_release_annotations[accession]["taxon_id"] = accessions_taxon[
                    accession
                ]["taxon_id"]

        if args.classification and pre_release_annotations:
            get_taxonomy_info(pre_release_annotations, accessions_taxon, args.rank)

        if args.ftp and pre_release_annotations:
            pre_release_annotations = add_ftp(pre_release_annotations, "pre")

        all_annotations: Dict[str, Dict[str, Any]] = {
            **live_annotations,
            **pre_release_annotations,
        }
    else:
        all_annotations = live_annotations

    # Console summary
    if args.bioproject_id:
        print(
            f"Found {len(accessions_taxon)} assemblies under BioProject ID {args.bioproject_id}"
        )
    elif args.taxon_id:
        print(f"Found {len(accessions_taxon)} assemblies for taxon ID {args.taxon_id}")

    print(
        f"Found {len(live_annotations)} annotations in beta.ensembl.org \
            for {len(unique_taxon_ids)} unique species"
    )

    if args.classification:
        rank_values = [
            details.get(args.rank, "unknown") for details in live_annotations.values()
        ]
        rank_counts = Counter(rank_values)
        print("\nBreakdown:")
        print(rank_counts)

    # Write final report
    write_report(
        all_annotations,
        args.report_file,
        args.rank,
        include_class=args.classification,
        include_ftp=args.ftp,
    )


if __name__ == "__main__":
    main()

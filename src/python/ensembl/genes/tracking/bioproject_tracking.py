#!/usr/bin/env python3
import argparse
import requests
import pymysql
import json
import logging
import re
from pathlib import Path
from collections import Counter
from typing import List, Tuple, Any, Dict, Optional
import subprocess
import shutil

# -----------------------------------
# Logging
# -----------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# -----------------------------------
# Config
# -----------------------------------
with open("./bioproject_tracking_config.json", "r") as f:
    config = json.load(f)

# -----------------------------------
# DB helper
# -----------------------------------
def mysql_fetch_data(
    query: str,
    params: Tuple = (),
    server_group: str = "meta",
    server_name: str = "beta",
    db_name: Optional[str] = None
) -> List[Tuple[Any, ...]]:
    """
    Executes a SQL query with optional parameters to fetch results from a specified MySQL server.

    Returns empty list on error or no data.
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
def get_assembly_accessions(
    query_id: str,
    query_type: str,
    only_haploid: bool = False
) -> Dict[str, Dict[str, int]]:
    """
    Prefer NCBI Datasets CLI (all versions), fall back to Datasets API (latest only).
    Supports query_type in {"bioproject", "taxon"}.
    Always returns a dict (possibly empty): {accession: {"taxon_id": int}}
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
            if accession:
                accs[accession] = {"taxon_id": taxon_id}
        return accs

    # --- 1) Try CLI (gets ALL versions) ---
    datasets_cli = shutil.which("datasets")
    if datasets_cli:
        try:
            if query_type == "bioproject":
                cmd = [
                    datasets_cli, "summary", "genome", "accession", str(query_id),
                    "--assembly-version", "all",
                    "--as-json-lines",
                ]
            elif query_type == "taxon":
                cmd = [
                    datasets_cli, "summary", "genome", "taxon", str(query_id),
                    "--assembly-version", "all",
                    "--as-json-lines",
                ]
            else:
                logging.error(f"Invalid query_type '{query_type}'")
                return {}

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            parsed = _parse_reports(result.stdout.splitlines())
            if parsed:
                logging.info(f"[Datasets CLI] Retrieved {len(parsed)} assemblies (all versions).")
                return parsed
            else:
                logging.warning("[Datasets CLI] Parsed zero rows from stdout; falling back to API.")
        except subprocess.CalledProcessError as e:
            logging.warning(f"[Datasets CLI] Exit {e.returncode}. STDERR:\n{e.stderr}\nFalling back to API.")
        except Exception as e:
            logging.warning(f"[Datasets CLI] Unexpected error: {e}. Falling back to API.")

    # --- 2) Fallback: API (latest only) ---
    try:
        base_url = config["urls"]["datasets"].get(query_type)
    except Exception:
        base_url = None

    if not base_url:
        logging.error("[Datasets API] Base URL not found in config; returning empty set.")
        return {}

    assembly_accessions: Dict[str, Dict[str, int]] = {}
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
                assembly_info = assembly.get("assembly_info", {}) or {}
                if only_haploid and assembly_info.get("assembly_type") != "haploid":
                    continue
                accession = assembly.get("accession")
                taxon_id = (assembly.get("organism") or {}).get("tax_id")
                if accession:
                    assembly_accessions[accession] = {"taxon_id": taxon_id}
            next_page_token = data.get("next_page_token")
            if not next_page_token:
                break
        except requests.RequestException as e:
            logging.error(f"[Datasets API] Error: {e}")
            break

    logging.info(f"[Datasets API] Retrieved {len(assembly_accessions)} assemblies (latest only).")
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
    """
    m = _CORE_VER_RE.search(dbname or "")
    core_ver = int(m.group(1)) if m else -1
    non_cm = 1 if "_cm_" not in (dbname or "").lower() else 0
    return (core_ver, non_cm, dbname or "")

def get_ensembl_live(accessions_taxon: Dict[str, Dict[str, int]]):
    """
    Returns (annotated_dict, missing_list).
    annotated_dict[accession] contains:
      - taxon_id (from input)
      - matches: List[{"guuid","dbname", optional "ftp","date"}]
      - guuid, dbname: best match (latest core; prefer non-cm)
    """
    accessions = list(accessions_taxon.keys())
    if not accessions:
        return {}, []

    query = """
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
      AND assembly.accession IN ({})
    """.format(", ".join(["%s"] * len(accessions)))

    rows = mysql_fetch_data(query, tuple(accessions))

    # Seed with taxon info so missing ones can still be tracked
    live_annotations: Dict[str, Dict[str, Any]] = {}

    for accession, guuid, dbname in rows:
        entry = live_annotations.setdefault(accession, dict(accessions_taxon.get(accession, {})))
        entry.setdefault("matches", [])
        # Avoid duplicates across joins if any
        if not any(m["guuid"] == guuid and m["dbname"] == dbname for m in entry["matches"]):
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
    rank: str
) -> Dict[str, Dict[str, Any]]:
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

# -----------------------------------
# Add FTP (per match + best)
# -----------------------------------
def add_ftp(
    annotations: Dict[str, Dict[str, Any]],
    release_type: str = "live"
) -> Dict[str, Dict[str, Any]]:
    for accession, annotation in annotations.items():
        matches = annotation.get("matches") or []

        if release_type == "pre":
            # Use dbname(s) to fetch scientific name
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
                    db_name=dbname
                )
                if not result:
                    continue
                scientific_name = result[0][0].replace(" ", "_")
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
                scientific_name = scientific_name.replace(" ", "_")
                source = source.lower()
                ftp_link = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{scientific_name}/{accession}/{source}/geneset/{date}/"

                if guuid in by_guuid:
                    by_guuid[guuid]["ftp"] = ftp_link
                    by_guuid[guuid]["date"] = date
                if annotation.get("guuid") == guuid:
                    annotation["ftp"] = ftp_link
                    annotation["date"] = date

    return annotations

# -----------------------------------
# Pre-release lookup
# -----------------------------------
def get_pre_release(missing_annotations: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Checks information_schema for pre-release databases matching the accession pattern.
    Returns dict: {accession: {"guuid": "unknown", "dbname": "<schema>"}}
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
            query,
            (lookup_string,),
            server_group="pre-release",
            server_name="gb1"
        )

        if result:
            pre_release_annotations[accession] = {
                "guuid": 'unknown',  # typo fixed from guiid -> guuid
                "dbname": result[0][0]
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
    include_ftp: bool = False
):
    """
    Writes a tab-separated report. If an accession has multiple matches,
    produces one line per match; otherwise one line.
    """
    report_path = Path(report_file)
    lines_written = 0

    with open(report_path, "w", encoding="utf-8") as file:
        for accession, details in live_annotations.items():
            matches = details.get("matches")

            # If we have multiple, write one row per match
            rows = matches if matches else [{
                "guuid": details.get("guuid"),
                "dbname": details.get("dbname"),
                "ftp": details.get("ftp"),
                "date": details.get("date")
            }]

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

    logging.info(f"Report written to {report_path.resolve()} with {lines_written} lines.")

# -----------------------------------
# Main
# -----------------------------------
def main():
    """
    CLI entrypoint.
    """
    parser = argparse.ArgumentParser(description="Fetch assembly data and report annotations.")
    parser.add_argument("--bioproject_id", type=str, help="NCBI BioProject ID")
    parser.add_argument("--taxon_id", type=str, help="Taxonomy ID")
    parser.add_argument("--haploid", action="store_true", help="Fetch only haploid assemblies")
    parser.add_argument("--report_file", type=str, default="./report_file.tsv")
    parser.add_argument("--classification", action="store_true", help="Provide breakdown of taxonomic classification")
    parser.add_argument("--rank", type=str, default="order")
    parser.add_argument("--ftp", action="store_true", help="Include FTP links in the report")
    parser.add_argument("--pre_release", action="store_true", help="Include list of pre-release databases in the report")

    args = parser.parse_args()

    # Determine which ID is provided and which query type is relevant
    if not args.bioproject_id and not args.taxon_id:
        parser.error("You must provide either --bioproject_id or --taxon_id")

    accessions_taxon = get_assembly_accessions(
        args.bioproject_id or args.taxon_id,
        "bioproject" if args.bioproject_id else "taxon",
        args.haploid
    )

    # Fetch annotations (live)
    live_annotations, missing_annotations = get_ensembl_live(accessions_taxon)

    # Taxonomy classification (top-level only; rows inherit from details)
    if args.classification:
        get_taxonomy_info(live_annotations, accessions_taxon, args.rank)

    # Optionally add FTP links
    if args.ftp:
        live_annotations = add_ftp(live_annotations, 'live')

    # Unique species count from live annotations
    unique_taxon_ids = {details['taxon_id'] for details in live_annotations.values() if 'taxon_id' in details}

    # Optionally add Pre-release data
    if args.pre_release and missing_annotations:
        pre_release_annotations = get_pre_release(missing_annotations)

        for accession in pre_release_annotations:
            if accession in accessions_taxon:
                pre_release_annotations[accession]["taxon_id"] = accessions_taxon[accession]["taxon_id"]

        if args.classification and pre_release_annotations:
            get_taxonomy_info(pre_release_annotations, accessions_taxon, args.rank)

        if args.ftp and pre_release_annotations:
            pre_release_annotations = add_ftp(pre_release_annotations, 'pre')

        all_annotations: Dict[str, Dict[str, Any]] = {**live_annotations, **pre_release_annotations}
    else:
        all_annotations = live_annotations

    # Console summary
    if args.bioproject_id:
        print(f"Found {len(accessions_taxon)} assemblies under BioProject ID {args.bioproject_id}")
    elif args.taxon_id:
        print(f"Found {len(accessions_taxon)} assemblies for taxon ID {args.taxon_id}")

    print(f"Found {len(live_annotations)} annotations in beta.ensembl.org for {len(unique_taxon_ids)} unique species")

    if args.classification:
        rank_values = [details.get(args.rank, "unknown") for details in live_annotations.values()]
        rank_counts = Counter(rank_values)
        print("\nBreakdown:")
        print(rank_counts)

    # Write final report
    write_report(all_annotations, args.report_file, args.rank, include_class=args.classification, include_ftp=args.ftp)

if __name__ == "__main__":
    main()

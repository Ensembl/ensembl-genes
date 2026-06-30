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
# pylint: disable=missing-module-docstring, consider-using-with, logging-not-lazy, logging-fstring-interpolation, consider-using-dict-items, unspecified-encoding,invalid-name , broad-exception-caught, line-too-long, redefined-outer-name, missing-timeout, unused-argument, missing-function-docstring
import argparse
import json
import logging
import logging.config
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pymysql
import requests
import xmltodict

try:
    from ensembl.genes.metadata.bioproject_from_registry import get_bioproject_names
except ImportError:
    from bioproject_from_registry import get_bioproject_names

# Module logger (configured in __main__ via logging.config)
logger: logging.Logger = logging.getLogger(__name__)


def mysql_fetch_data(
    query: str, database: str, host: str, port: int, user: str
) -> List[Tuple[Any, ...]]:
    """
    Run a simple SELECT query and return fetched rows.
    Returns an empty list on error.

    Args:
        query (str): SQL SELECT query to execute.
        database (str): Name of the database.
        host (str): Host address of the MySQL server.
        port (int): Port number of the MySQL server.
        user (str): Username to connect to the database.

    Returns:
        List of tuples representing query results.
    """
    conn = None
    cursor = None
    info: List[Tuple[Any, ...]] = []
    try:
        conn = pymysql.connect(
            host=host, user=user, port=port, database=database.strip()
        )
        cursor = conn.cursor()
        cursor.execute(query)
        info = list(cursor.fetchall())
    except pymysql.Error:
        logger.exception("MySQL error while executing query: %s", query)
    finally:
        if cursor is not None:
            try:
                cursor.close()
            except Exception:
                logger.exception("Error closing cursor")
        if conn is not None:
            try:
                conn.close()
            except Exception:
                logger.exception("Error closing connection")
    return info


def get_ena_metadata(accession: str, truth_dict: Dict[str, Any]) -> Dict[str, str]:
    """Fetch assembly metadata from ENA XML API.

    Args:
        accession (str): Assembly accession (e.g., GCA_000001405.28)
        truth_dict (Dict[str, Any]): Dictionary containing truth values for metadata.

    Returns:
        Dict[str, str]: Dictionary containing metadata from ENA.
    """
    return_dict: Dict[str, str] = {}
    assembly_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    assembly_xml = requests.get(assembly_url)
    assembly_dict = xmltodict.parse(assembly_xml.text)

    assembly_attribs = assembly_dict["ASSEMBLY_SET"]["ASSEMBLY"]["ASSEMBLY_ATTRIBUTES"][
        "ASSEMBLY_ATTRIBUTE"
    ]
    # fix if single or multiple items
    if isinstance(assembly_attribs, dict):
        assembly_attribs = [assembly_attribs]

    # assembly metadata
    return_dict["assembly.name"] = (
        (
            assembly_dict.get("ASSEMBLY_SET", {}).get("ASSEMBLY", {}).get("NAME", "")
        ).replace(" ", "_")
        if assembly_dict.get("ASSEMBLY_SET", {}).get("ASSEMBLY", {}).get("NAME")
        else ""
    )

    return_dict["assembly.level"] = (
        assembly_dict.get("ASSEMBLY_SET", {})
        .get("ASSEMBLY", {})
        .get("ASSEMBLY_LEVEL", "")
    )

    for attrib_set in assembly_attribs:
        # assembly date
        if attrib_set.get("TAG") == "ENA-LAST-UPDATED":
            raw_date = attrib_set.get("VALUE", "")
            # normalise from YYYY-MM-DD to YYYY-MM (or leave shorter strings as-is)
            if len(raw_date) >= 7:
                return_dict["assembly.date"] = raw_date[:7]
            else:
                return_dict["assembly.date"] = raw_date

    # organism metadata: sample id
    if "organism.biosample_id" not in truth_dict:
        biosample_id = (
            assembly_dict.get("ASSEMBLY_SET", {})
            .get("ASSEMBLY", {})
            .get("SAMPLE_REF", {})
            .get("IDENTIFIERS", {})
            .get("PRIMARY_ID", "")
        )
        if biosample_id:
            return_dict["organism.biosample_id"] = biosample_id
        else:
            logger.critical(
                " | BIOSAMPLE_ID | organism.biosample_id could not be found in the ENA metadata, this is a required key!"
            )

    # taxonomy id (workaround for bad records)
    if accession in ("GCA_944452655.1", "GCA_944452715.1"):
        return_dict["organism.taxonomy_id"] = "1539398"
    else:
        tax_id = (
            assembly_dict.get("ASSEMBLY_SET", {})
            .get("ASSEMBLY", {})
            .get("TAXON", {})
            .get("TAXON_ID", "")
        )
        if tax_id:
            return_dict["organism.taxonomy_id"] = tax_id
        else:
            logger.warning(
                " | TAXONOMY_ID | organism.taxonomy_id could not be found in the ENA metadata"
            )

    return return_dict


def get_ncbi_metadata(
    accession: str, assembly_name: str, scientific_name: str, search: str
) -> Dict[str, str]:
    """Fetch assembly metadata from NCBI assembly report.

    Args:
        accession (str): Assembly accession (e.g., GCA_000001405.28)
        assembly_name (str): Assembly name (e.g., GRCh38)
        scientific_name (str): Scientific name of the organism
        search (str): Type of metadata to search for ("ucsc" or "biosample")

    Returns:
        Dict[str, str]: Dictionary containing metadata from NCBI.
    """
    organism = f"{accession}_{assembly_name}"
    ncbi_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[0:3]}/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/{organism}/{organism}_assembly_report.txt"
    ncbi_return = requests.get(ncbi_url).text.splitlines()

    return_dict: Dict[str, str] = {
        "assembly.ucsc_alias": "",
        "organism.strain": "",
        "organism.strain_type": "",
    }

    if search == "ucsc":
        for line in ncbi_return:
            ucsc_match = re.search(r"# Synonyms:\s*([A-Za-z0-9]+)", line)
            if ucsc_match:
                return_dict["assembly.ucsc_alias"] = ucsc_match.group(1)

    elif search == "biosample":
        for line in ncbi_return:
            strain_match = re.search(
                r"# Infraspecific name:\s*([A-Za-z]+)=([A-Za-z0-9 \-\./]+)", line
            )
            if strain_match:
                return_dict["organism.strain_type"] = strain_match.group(1)
                return_dict["organism.strain"] = strain_match.group(2)

    return return_dict


def get_biosample_metadata(
    biosample_id: str, assembly_accession: str, assembly_name: str, scientific_name: str
) -> Dict[str, str]:
    """Fetch assembly metadata from EBI BioSample API.

    Args:
        biosample_id (str): BioSample accession (e.g., SAMN00000001)
        assembly_accession (str): Assembly accession (e.g., GCA_000001405.28)
        assembly_name (str): Assembly name (e.g., GRCh38)
        scientific_name (str): Scientific name of the organism

    Returns:
        Dict[str, str]: Dictionary containing metadata from BioSample.
    """

    biosample_url = f"https://www.ebi.ac.uk/biosamples/samples/{biosample_id}"
    biosample_return = requests.get(biosample_url).text
    return_dict: Dict[str, str] = {"organism.strain": "", "organism.strain_type": ""}

    try:
        biosample_data = json.loads(biosample_return)

        strain_types = {"population", "race", "ecotype", "breed", "strain", "cultivar"}

        for strain_type in strain_types:
            try:
                return_dict["organism.strain"] = biosample_data["characteristics"][
                    strain_type
                ][0]["text"]
                return_dict["organism.strain_type"] = strain_type
                break
            except KeyError:
                continue

        if return_dict["organism.strain"] == "Caucasian":
            return_dict["organism.strain"] = "European"

        try:
            return_dict["assembly.tol_id"] = biosample_data["characteristics"]["tolid"][
                0
            ]["text"]
        except KeyError:
            return_dict["assembly.tol_id"] = ""

    except json.decoder.JSONDecodeError:
        # fallback to NCBI
        biosample_dict = get_ncbi_metadata(
            assembly_accession, assembly_name, scientific_name, "biosample"
        )
        return_dict["organism.strain"] = biosample_dict.get("organism.strain", "")
        return_dict["organism.strain_type"] = biosample_dict.get(
            "organism.strain_type", ""
        )
        return_dict["assembly.tol_id"] = ""

    return return_dict


def normalise_text(value: Any) -> str:
    """Convert generated metadata values to the string representation used by meta."""
    if value is None:
        return ""
    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d")
    return str(value)


def serialise_json(value: Any) -> Any:
    """Make generated metadata safe for JSON output."""
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: serialise_json(val) for key, val in value.items()}
    if isinstance(value, list):
        return [serialise_json(item) for item in value]
    if isinstance(value, tuple):
        return [serialise_json(item) for item in value]
    return value


def sql_escape(value: Any) -> str:
    """Escape a value for the generated SQL patch."""
    return normalise_text(value).replace("'", "''")


def write_metadata_json(  # pylint: disable=too-many-arguments
    json_path: Path,
    db_name: str,
    species_id: int,
    core_dict: Dict[str, Any],
    truth_dict: Dict[str, Any],
    sql_actions: List[Dict[str, str]],
) -> None:
    """Write a JSON artifact with all metadata generated by the legacy script."""
    payload = {
        "database": db_name,
        "species_id": species_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "core_metadata": {
            key: normalise_text(value) for key, value in core_dict.items()
        },
        "metadata": {key: normalise_text(value) for key, value in truth_dict.items()},
        "sql_actions": sql_actions,
    }

    json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(json_path, "w", encoding="utf-8") as json_out:
        json.dump(serialise_json(payload), json_out, indent=2, sort_keys=True)
        json_out.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare SQL updates for core dbs")
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default=".",
        help="Path where the output and temp files will write to. \
        Uses current dir by default",
    )
    parser.add_argument(
        "-d",
        "--db_name",
        help="Database name",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--host",
        help="Host server",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--port",
        help="Host server port",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--production_name",
        help="species.production_name",
        required=False,
    )
    parser.add_argument(
        "-t",
        "--team",
        required=True,
        type=lambda x: x.capitalize(),
        help="Team responsible for the database",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output (check that all required keys are not NULL/ empty)",
    )
    parser.add_argument(
        "--json_output",
        help="Path for JSON metadata output. Defaults to <output_name>.json in output_dir.",
    )

    args = parser.parse_args()

    server_info = {
        "staging": {
            "db_host": args.host,
            "db_port": int(args.port),
            "db_user": "ensro",
        },
        "meta": {
            "db_host": "mysql-ens-meta-prod-1.ebi.ac.uk",
            "db_port": 4483,
            "db_user": "ensro",
            "db_pass": "",
        },
    }

    metadata_dir = Path(__file__).parent
    provider_static_file = metadata_dir / "provider_static.txt"
    ref_static_file = metadata_dir / "ref_static.txt"
    snp_static_file = metadata_dir / "snp_static.txt"
    url_static_file = metadata_dir / "url_static.txt"

    db = args.db_name
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_name = args.production_name or db
    sql_out = open(output_dir / f"{output_name}.sql", "w")
    json_path = (
        Path(args.json_output)
        if args.json_output
        else output_dir / f"{output_name}.json"
    )

    print(f"Working on database: {db}")
    print(f"USE {db};", file=sql_out)

    # set up logger
    if args.production_name:
        log_file_path = output_dir / f"{args.production_name}_metadata.log"
    else:
        log_file_path = output_dir / f"{db}_metadata.log"
    log_ini_path = metadata_dir / "logging.conf"
    logging.config.fileConfig(
        log_ini_path,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=True,
    )
    logger = logging.getLogger()
    logger.propagate = False

    core_dict = {}
    if args.production_name:
        core_dict["species.production_name"] = args.production_name
        # getting species_id via the production_name
        species_query = f"SELECT species_id FROM meta WHERE meta_key='species.production_name' AND meta_value='{core_dict['species.production_name']}';"
        species_meta = mysql_fetch_data(
            species_query,
            host=server_info["staging"]["db_host"],
            user=server_info["staging"]["db_user"],
            port=server_info["staging"]["db_port"],
            database=db,
        )
        species_id = species_meta[0][0]
        print("species ID: " + str(species_id))
    else:
        species_id = 1
        print(
            "WARNING: no production_name provided, using default species ID ="
            + str(species_id)
        )

    # get all existing assembly, species and genebuild metadata from the core db
    core_query = f"SELECT meta_key,meta_value FROM meta WHERE species_id = {species_id} AND (meta_key LIKE 'assembly%' OR meta_key LIKE 'species%' OR meta_key LIKE 'genebuild%' OR meta_key LIKE 'genome%' OR meta_key LIKE 'organism%' OR meta_key LIKE 'sample%' OR meta_key LIKE 'annotation%' OR meta_key LIKE 'gencode%');"
    print(core_query)
    core_meta = mysql_fetch_data(
        core_query,
        host=server_info["staging"]["db_host"],
        user=server_info["staging"]["db_user"],
        port=server_info["staging"]["db_port"],
        database=db,
    )
    core_genome_groups = []
    for meta_pair in core_meta:
        core_dict[meta_pair[0]] = meta_pair[1]
        if meta_pair[0] == "genome.genome_group":
            core_genome_groups.append(meta_pair[1])

    # get the assembly metadata from the sources of truth (sources of truth in parentheses)
    # expected assembly meta_keys: assembly.accession (from core), assembly.date (from ena), assembly.is_reference (static), assembly.name (ena), assembly.provider_name (core or default), assembly.provider_url (core or default), assembly.level (ena), assembly.tolid (biosample), assembly.ucsc_alias (ncbi), assembly.long_name, assembly.url_name (static)
    # expected genome meta_keys: genome.genome_group (registry)
    # expected organism meta_keys: organism.taxonomy_id (ena), organism.species_taxonomy_id (taxonomy db), organism.common_name (taxonomy db), organism.strain (biosample), organism.scientific_name (taxonomy db), organism.scientific_parlance_name (static), organism.strain_type (biosample), organism.sample_accession (ena)
    # expected genebuild meta_keys: genebuild.initial_release_date, genebuild.last_geneset_update, genebuild.level, genebuild.method, genebuild.method_display, genebuild.start_date, genebuild.version (create and check and required), genebuild.sample_gene (core), genebuild.sample_location (core), genebuild.id, genebuild.projection_source_db, genebuild.havana_datafreeze_date, genebuild.provider_name (static or core or default), genebuild.provider_url (static or core or default), genebuild.annotation_source (core or default)
    truth_dict = {}

    # now some assembly.accession values will be GCFs - that breaks finding things based on a GCA
    gca_accession = core_dict["assembly.accession"]
    if core_dict["assembly.accession"].startswith("GCF"):
        gca_accession = core_dict["assembly.alt_accession"]
        truth_dict["assembly.alt_accession"] = core_dict["assembly.alt_accession"]
        truth_dict["assembly.accession_refseq"] = core_dict["assembly.accession"]
        truth_dict["assembly.accession_body"] = "RefSeq"
    elif "assembly.accession_refseq" in core_dict:
        # if the core comes from the main site, we need to put the GCF back as the main accession
        truth_dict["assembly.alt_accession"] = core_dict["assembly.alt_accession"]
        core_dict["assembly.accession"] = core_dict["assembly.accession_refseq"]
        truth_dict["assembly.accession_body"] = "RefSeq"

    # Some goddamn ENA assembly records do not have the BioSample ID, so I'm hardcoding them, I swear to jaysus, I'm so done with metadata!!!
    if db == "caenorhabditis_elegans_core_57_110_282":
        truth_dict["organism.biosample_id"] = "SAMN04256190"
    elif db == "ciona_intestinalis_core_110_3":
        truth_dict["organism.biosample_id"] = "SAMD00414333"
    elif db == "homo_sapiens_37_core_110_37":
        truth_dict["organism.biosample_id"] = "SAMN12121739"
    elif db == "homo_sapiens_core_110_38":
        truth_dict["organism.biosample_id"] = "SAMN12121739"
    elif db == "saccharomyces_cerevisiae_core_57_110_4":
        truth_dict["organism.biosample_id"] = "SAMEA3184125"
    elif db == "mus_musculus_core_110_39":
        truth_dict["organism.biosample_id"] = "SAMN26853311"
    elif db == "bos_taurus_core_110_1":
        truth_dict["organism.biosample_id"] = "SAMN03145444"

    bioproject_names = get_bioproject_names(
        gca_accession, user=server_info["meta"]["db_user"]
    )

    # get metadata from ENA records
    try:
        truth_dict.update(get_ena_metadata(gca_accession, truth_dict))
    except KeyError:
        logger.critical("No assembly accession found, cannot process any further")

    # get common and scientific names from NCBI taxonomy (seems inefficient to query the db twice, but I don't think it returns ordered results and I don't want to risk mixing them up)
    try:
        name_query = f"SELECT name FROM ncbi_taxa_name WHERE taxon_id={truth_dict['organism.taxonomy_id']} AND name_class='genbank common name';"
        name_info = mysql_fetch_data(
            name_query,
            host=server_info["meta"]["db_host"],
            user=server_info["meta"]["db_user"],
            port=server_info["meta"]["db_port"],
            database="ncbi_taxonomy",
        )
        truth_dict["organism.common_name"] = (name_info[0][0]).capitalize()
    except IndexError:  # not everything has a genbank common name
        truth_dict["organism.common_name"] = ""
    s_name_query = f"SELECT name FROM ncbi_taxa_name WHERE taxon_id={truth_dict['organism.taxonomy_id']} AND name_class='scientific name';"
    s_name_info = mysql_fetch_data(
        s_name_query,
        host=server_info["meta"]["db_host"],
        user=server_info["meta"]["db_user"],
        port=server_info["meta"]["db_port"],
        database="ncbi_taxonomy",
    )
    truth_dict["organism.scientific_name"] = (s_name_info[0][0]).capitalize()

    # get metadata from NCBI taxonomy
    taxonomy_query = f"SELECT rank, parent_id FROM ncbi_taxa_node WHERE taxon_id={truth_dict['organism.taxonomy_id']};"
    rank_info = mysql_fetch_data(
        taxonomy_query,
        host=server_info["meta"]["db_host"],
        user=server_info["meta"]["db_user"],
        port=server_info["meta"]["db_port"],
        database="ncbi_taxonomy",
    )
    rank_dict = {}
    for pair in rank_info:
        rank_dict["rank"] = pair[0]
        rank_dict["parent_id"] = str(pair[1])

    if rank_dict["rank"] == "species":
        truth_dict["organism.species_taxonomy_id"] = truth_dict["organism.taxonomy_id"]
    else:
        while rank_dict["rank"] != "species":
            current_id = str(rank_dict["parent_id"])
            parent_taxonomy_query = f"SELECT rank, parent_id FROM ncbi_taxa_node WHERE taxon_id={rank_dict['parent_id']};"
            parent_rank_info = mysql_fetch_data(
                parent_taxonomy_query,
                host=server_info["meta"]["db_host"],
                user=server_info["meta"]["db_user"],
                port=server_info["meta"]["db_port"],
                database="ncbi_taxonomy",
            )
            for pair in parent_rank_info:
                rank_dict["rank"] = pair[0]
                rank_dict["parent_id"] = str(pair[1])

            truth_dict["organism.species_taxonomy_id"] = current_id

    # get metadata from NCBI records
    truth_dict.update(
        get_ncbi_metadata(
            core_dict["assembly.accession"],
            truth_dict["assembly.name"],
            truth_dict["organism.scientific_name"],
            "ucsc",
        )
    )

    # get metadata from BioSample records
    if truth_dict["organism.biosample_id"] != "":
        truth_dict.update(
            get_biosample_metadata(
                truth_dict["organism.biosample_id"],
                gca_accession,
                truth_dict["assembly.name"],
                truth_dict["organism.scientific_name"],
            )
        )
    else:
        # there's probably a better source for ToLIDs - DToL portal?
        truth_dict["assembly.tol_id"] = ""

    # now to get some metadata that can only come from static files
    # scientific_parlance_name
    snp_list = open(snp_static_file).readlines()
    for line in snp_list:
        if truth_dict["organism.scientific_name"] in line:
            snp = line.split("\t")[1].strip()
            truth_dict["organism.scientific_parlance_name"] = snp

    # assembly.url_name
    url_list = open(url_static_file).readlines()
    for line in url_list:
        if gca_accession in line:
            url = line.split("\t")[1].strip()
            truth_dict["assembly.url_name"] = url

    # assembly.is_reference
    ref_list = open(ref_static_file).readlines()
    for line in ref_list:
        if truth_dict["organism.scientific_name"] in line:
            ref_accession = line.split("\t")[1].strip()
            if gca_accession == ref_accession:
                truth_dict["assembly.is_reference"] = 1

    # now to create some values
    # assembly provider and url - if not in core already, set to default provider,"ENA"
    if "assembly.provider_name" not in core_dict:
        truth_dict["assembly.provider_name"] = "ENA"
        truth_dict["assembly.provider_url"] = "https://www.ebi.ac.uk/ena/browser/home"
        logger.warning(
            " | ASSEMBLY_PROVIDER | No assembly provider information found, using default, ENA."
        )

    provider_dict = {}
    with open(provider_static_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > 2:
                provider_dict[parts[0]] = {"name": parts[1], "url": parts[2]}

    # genebuild.provider_name and _url keys: check if annotation.provider_name and _url keys exist, if not check spreadsheet, else set to default "Ensembl", "Ensembl url" (full_genebuild, anno, braker, hprc)
    # genebuild.version key is now being used for metadata loading, we are setting it to the first genebuild.version for everything (unless there is already a value set for this key), data teams need to be aware when handing over an updated annotation, i.e. assembly is same as an existing genome, but the gene set has been updated (new data, or a fix)
    # maybe add a check on the metadata database here - does the label GCA_XXXX_ENSXX match an existing dataset, then warn! (could be part of the check for required keys!)
    if "gencode.version" in core_dict:
        truth_dict["genebuild.version"] = core_dict["gencode.version"].replace(" ", "")
        truth_dict["genebuild.method_display"] = "Manual annotation"
    elif "genebuild.method" in core_dict:
        # a quick check where sample genes have BRAKER stable id prefixes, because some braker annotations had incorrect value for genebuild.method (I'm doing this with the sample gene text as it will contain the BRAKER stable id prefix)
        if "BRAKER" in core_dict["sample.gene_text"]:
            truth_dict["genebuild.version"] = "BRK01"
            truth_dict["genebuild.method"] = "braker"
            truth_dict["genebuild.method_display"] = "BRAKER2"
            truth_dict["genebuild.annotation_source"] = "braker"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = (
                "https://beta.ensembl.org/help/articles/braker-2-genome-annotation"
            )
        elif "HELIXER" in core_dict["sample.gene_text"]:
            truth_dict["genebuild.version"] = "HLX01"
            truth_dict["genebuild.method"] = "helixer"
            truth_dict["genebuild.method_display"] = "Helixer"
            truth_dict["genebuild.annotation_source"] = "helixer"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = (
                "https://beta.ensembl.org/help/articles/helixer-genome-annotation"
            )
        # otherwise check the genebuild.method for annotation key updates/additions
        elif core_dict["genebuild.method"] == "full_genebuild":
            truth_dict["genebuild.version"] = "ENS01"
            truth_dict["genebuild.method_display"] = "Ensembl Genebuild"
            truth_dict["genebuild.annotation_source"] = "ensembl"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = (
                "https://beta.ensembl.org/help/articles/vertebrate-genome-annotation"
            )
        elif core_dict["genebuild.method"] == "anno":
            truth_dict["genebuild.version"] = "ENS01"
            truth_dict["genebuild.method_display"] = "Ensembl Genebuild"
            truth_dict["genebuild.annotation_source"] = "ensembl"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            truth_dict["genebuild.provider_url"] = (
                "https://beta.ensembl.org/help/articles/non-vertebrate-genome-annotation"
            )
        elif core_dict["genebuild.method"] == "projection_build":
            truth_dict["genebuild.version"] = "ENS01"
            truth_dict["genebuild.annotation_source"] = "ensembl"
            truth_dict["genebuild.provider_name"] = "Ensembl"
            if core_dict["species.scientific_name"] == "homo sapiens":
                truth_dict["genebuild.provider_url"] = (
                    "https://beta.ensembl.org/help/articles/human-genome-automated-annotation"
                )
                truth_dict["genebuild.method_display"] = "Mapping from GRCh38"
            else:
                truth_dict["genebuild.provider_url"] = (
                    "https://beta.ensembl.org/help/articles/vertebrate-genome-annotation"
                )
                truth_dict["genebuild.method_display"] = "Mapping from reference"
        # for non-genebuilds, I don't have a set of rules (as above), so I rely on the core meta keys and the static
        elif core_dict["genebuild.method"] == "import":
            if "genebuild.method_display" not in core_dict:
                truth_dict["genebuild.method_display"] = "Import"
            if "genebuild.version" not in core_dict:
                truth_dict["genebuild.version"] = "EXT01"
            # try to get the annotation source
            if (
                "species.annotation_source" in core_dict
                and "genebuild.annotation_source" not in core_dict
            ):
                truth_dict["genebuild.annotation_source"] = core_dict[
                    "species.annotation_source"
                ]
            elif "genebuild.annotation_source" not in core_dict:
                logger.critical(
                    " | GENEBUILD_PROVIDER_URL | No genebuild.annotation_source could be found in the core db, this is a required key!"
                )
            # try to get the annotation provider name from core
            if (
                "annotation.provider_name" in core_dict
                and "genebuild.provider_name" not in core_dict
            ):
                truth_dict["genebuild.provider_name"] = core_dict[
                    "annotation.provider_name"
                ]
            # otherwise get it from provider_static.txt
            elif (
                core_dict["species.production_name"] in provider_dict
                and "genebuild.provider_name" not in core_dict
            ):
                truth_dict["genebuild.provider_name"] = provider_dict[
                    core_dict["species.production_name"]["name"]
                ]
            elif "genebuild.provider_name" not in core_dict:
                logger.critical(
                    " | GENEBUILD_PROVIDER_NAME | No genebuild.provider_name could be found either in the core db or in provider_static.txt, this is a required key!"
                )
            # try to get the annotation provider url from core
            if (
                "annotation.provider_url" in core_dict
                and "genebuild.provider_url" not in core_dict
            ):
                truth_dict["genebuild.provider_url"] = core_dict[
                    "annotation.provider_url"
                ]
            # otherwise get it from provider_static.txt
            elif (
                core_dict["species.production_name"] in provider_dict
                and "genebuild.provider_url" not in core_dict
            ):
                truth_dict["genebuild.provider_url"] = provider_dict[
                    core_dict["species.production_name"]["url"]
                ]
            elif "genebuild.provider_url" not in core_dict:
                logger.critical(
                    " | GENEBUILD_PROVIDER_URL | No genebuild.provider_url could be found either in the core db or in provider_static.txt, this is a required key!"
                )

    else:
        logger.warning(
            " | GENEBUILD_METHOD | No genebuild.method, assuming this is an Ensembl annotation"
        )
        truth_dict["genebuild.version"] = "ENS01"
        truth_dict["genebuild.annotation_source"] = "ensembl"
        truth_dict["genebuild.provider_name"] = "Ensembl"
        truth_dict["genebuild.provider_url"] = (
            "https://beta.ensembl.org/help/articles/vertebrate-genome-annotation"
        )

    # if the genebuild.version already exists in the core db, I'll just leave that value
    if "genebuild.version" in core_dict and core_dict["genebuild.version"] != "":
        truth_dict["genebuild.version"] = core_dict["genebuild.version"]

    # if the sample gene info is in the core db, update the meta key names
    try:
        truth_dict["genebuild.sample_gene"] = core_dict["sample.gene_param"]
    except KeyError:
        logger.critical(
            " | SAMPLE.GENE_PARAM | No sample.gene_param could be found in the core db, this is a required key!"
        )
    try:
        truth_dict["genebuild.sample_location"] = core_dict["sample.location_param"]
    except KeyError:
        logger.critical(
            " | SAMPLE.LOCATION_PARAM | No sample.location_param could be found in the core db, this is a required key!"
        )

    # let's do a check for the genebuild.last_geneset_update key because it's required by web
    try:
        truth_dict["genebuild.last_geneset_update"] = core_dict[
            "genebuild.last_geneset_update"
        ]
    except KeyError:
        logger.warning(
            " | GENEBUILD.LAST_GENESET_UPDATE | No genebuild.last_geneset_update could be found in the core db, this is a required key, I'm setting it to today."
        )
        truth_dict["genebuild.last_geneset_update"] = datetime.now().strftime("%Y-%m")

    # update the remaining species keys -> organisms keys
    try:
        truth_dict["organism.production_name"] = core_dict["species.production_name"]
    except KeyError:
        logger.critical(
            " | SPECIES.PRODUCTION_NAME | No species.production_name could be found in the core db, this is a required key!"
        )

    # Set the team responsible for this genome
    truth_dict["genebuild.team_responsible"] = args.team

    # genome.genome_group can have multiple values, so handle it outside the
    # single-value truth_dict update path.
    for bioproject_name in bioproject_names:
        if bioproject_name not in core_genome_groups:
            meta_value = bioproject_name.replace("'", "''")
            print(
                f"INSERT IGNORE INTO meta (species_id, meta_key, meta_value) "
                f"VALUES({species_id}, 'genome.genome_group', '{meta_value}');",
                file=sql_out,
            )
    sql_actions = []

    # if the expected meta_key does not exist in the core metadata but there is a truth value, INSERT
    # if the expected meta_key exists in the core metadata but the meta_value does not match the truth value and truth value is not NULL, UPDATE
    # if the expected meta_key exists in the core metadata and matches the truth value, NO ACTION
    for meta_key, raw_value in truth_dict.items():
        # Escape single quotes in truth_dict[meta_key]
        meta_value = normalise_text(raw_value)
        escaped_value = sql_escape(meta_value)

        if meta_key not in core_dict:
            if meta_value:
                sql = (
                    "INSERT IGNORE INTO meta (species_id, meta_key, meta_value) "
                    f"VALUES({species_id}, '{meta_key}', '{escaped_value}');"
                )
                print(sql, file=sql_out)
                sql_actions.append(
                    {
                        "action": "insert",
                        "meta_key": meta_key,
                        "meta_value": meta_value,
                        "sql": sql,
                    }
                )
        elif meta_value != normalise_text(core_dict[meta_key]) and meta_value != "":
            sql = (
                f"UPDATE meta SET meta_value='{escaped_value}' "
                f"WHERE meta_key='{meta_key}';"
            )
            print(sql, file=sql_out)
            sql_actions.append(
                {
                    "action": "update",
                    "meta_key": meta_key,
                    "meta_value": meta_value,
                    "sql": sql,
                }
            )

    # do a check on required keys
    required_meta_keys = [
        "organism.production_name",
        "organism.taxonomy_id",
        "organism.scientific_name",
        "assembly.accession",
        "assembly.name",
        "genebuild.version",
        "genebuild.method",
        "genebuild.method_display",
        "genebuild.start_date",
        "genebuild.annotation_source",
        "genebuild.provider_name",
        "genebuild.provider_url",
        "genebuild.sample_gene",
        "genebuild.sample_location",
        "genebuild.last_geneset_update",
        "organism.biosample_id",
    ]
    for required_key in required_meta_keys:
        if required_key not in core_dict:
            if required_key not in truth_dict:
                logger.critical(f"You are missing required meta key: {required_key}")

    # report if there has been a change in common name
    if "species.common_name" in core_dict:
        if (
            core_dict["species.common_name"].lower()
            != truth_dict["organism.common_name"].lower()
        ):
            logger.warning(
                ' | COMMON NAME | The value for species.common_name in your meta table: "'
                + core_dict["species.common_name"]
                + '" does not match the value that I am assigning to organism.common_name: "'
                + truth_dict["organism.common_name"]
                + '"'
            )

    meta_keys_to_remove = [
        "species.strain",
        "strain.type",
    ]
    for remove_key in meta_keys_to_remove:
        sql = f"DELETE from meta WHERE meta_key='{remove_key}';"
        print(sql, file=sql_out)
        sql_actions.append(
            {
                "action": "delete",
                "meta_key": remove_key,
                "meta_value": "",
                "sql": sql,
            }
        )

    sql_out.close()
    write_metadata_json(json_path, db, species_id, core_dict, truth_dict, sql_actions)

    if args.verbose:
        # Print with the 'required' column if -v flag is present
        print(
            "\n".join(
                f"{key:<28}   {val:<20}   {'required' if key in required_meta_keys else ''}"
                for key, val in truth_dict.items()
            )
        )

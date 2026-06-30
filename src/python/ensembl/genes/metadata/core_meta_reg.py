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
# pylint: disable=missing-module-docstring, logging-not-lazy, logging-fstring-interpolation, line-too-long, too-many-lines, too-many-locals, too-many-branches, too-many-statements, broad-exception-caught

import argparse
import json
import logging
import logging.config
import os
from datetime import date, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Tuple, TypeAlias

import pymysql
from ensembl.genes.info_from_registry.registry_helper import (
    fetch_current_genebuild_record,
)

if TYPE_CHECKING:
    RegistryConnection: TypeAlias = pymysql.connections.Connection[
        pymysql.cursors.DictCursor
    ]
else:
    RegistryConnection = pymysql.connections.Connection

logger: logging.Logger = logging.getLogger(__name__)

DEFAULT_REGISTRY_DB = "gb_assembly_metadata"
DEFAULT_ASSEMBLY_PROVIDER_NAME = "ENA"
DEFAULT_ASSEMBLY_PROVIDER_URL = "https://www.ebi.ac.uk/ena/browser/home"
DEFAULT_ENSEMBL_PROVIDER_NAME = "Ensembl"

DEFAULT_VERTEBRATE_GENEBUILD_URL = (
    "https://beta.ensembl.org/help/articles/vertebrate-genome-annotation"
)
DEFAULT_NONVERTEBRATE_GENEBUILD_URL = (
    "https://beta.ensembl.org/help/articles/non-vertebrate-genome-annotation"
)
DEFAULT_BRAKER_URL = "https://beta.ensembl.org/help/articles/braker-2-genome-annotation"
DEFAULT_HELIXER_URL = "https://beta.ensembl.org/help/articles/helixer-genome-annotation"
DEFAULT_HUMAN_MAPPING_URL = (
    "https://beta.ensembl.org/help/articles/human-genome-automated-annotation"
)

BIOSAMPLE_OVERRIDES = {
    "caenorhabditis_elegans_core_57_110_282": "SAMN04256190",
    "ciona_intestinalis_core_110_3": "SAMD00414333",
    "homo_sapiens_37_core_110_37": "SAMN12121739",
    "homo_sapiens_core_110_38": "SAMN12121739",
    "saccharomyces_cerevisiae_core_57_110_4": "SAMEA3184125",
    "mus_musculus_core_110_39": "SAMN26853311",
    "bos_taurus_core_110_1": "SAMN03145444",
}

REQUIRED_META_KEYS = [
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

META_KEYS_TO_REMOVE = [
    "species.strain",
    "strain.type",
]


def mysql_fetch_data(  # pylint: disable=too-many-arguments
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str = "",
    params: Optional[Tuple[Any, ...] | List[Any]] = None,
) -> List[Dict[str, Any]]:
    """
    Run a SELECT query and return fetched rows as dictionaries.
    Returns an empty list on error.
    """
    conn = None
    cursor = None
    info: List[Dict[str, Any]] = []
    try:
        conn = pymysql.connect(
            host=host,
            user=user,
            port=int(port),
            password=password,
            database=database.strip(),
            cursorclass=pymysql.cursors.DictCursor,
        )
        cursor = conn.cursor()
        cursor.execute(query, params or ())
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


def normalise_text(value: Any) -> str:
    """Convert registry/core values to the string representation used by meta."""
    if value is None:
        return ""
    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d")
    if isinstance(value, date):
        return value.strftime("%Y-%m-%d")
    return str(value)


def normalise_month(value: Any) -> str:
    """Normalise a date-like value to YYYY-MM, matching the old script output."""
    text = normalise_text(value)
    if len(text) >= 7:
        return text[:7]
    return text


def capitalise_first(value: Any) -> str:
    """Match the old metadata script's capitalization behaviour."""
    text = normalise_text(value)
    return text.capitalize() if text else ""


def serialise_json(value: Any) -> Any:
    """Make database values safe for JSON output."""
    if isinstance(value, (datetime, date)):
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


def read_two_column_static_file(path: Path) -> Dict[str, str]:
    """Load static metadata files keyed by their first tab-delimited column."""
    result: Dict[str, str] = {}
    if not path.exists():
        logger.warning("Static metadata file does not exist: %s", path)
        return result

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [part.strip() for part in line.split("\t")]
            if len(parts) >= 2 and parts[0]:
                result[parts[0]] = parts[1]
    return result


def read_provider_static_file(path: Path) -> Dict[str, Dict[str, str]]:
    """Load provider fallback values, handling mixed tab/space separated rows."""
    result: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        logger.warning("Provider static file does not exist: %s", path)
        return result

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            if "\t" in line:
                parts = [part.strip() for part in line.split("\t") if part.strip()]
                if len(parts) >= 3:
                    result[parts[0]] = {"name": parts[-2], "url": parts[-1]}
                continue

            parts = line.split()
            url_index = next(
                (index for index, part in enumerate(parts) if part.startswith("http")),
                None,
            )
            if url_index and url_index >= 2:
                key = " ".join(parts[: url_index - 1])
                result[key] = {"name": parts[url_index - 1], "url": parts[url_index]}

    return result


def setup_logging(output_dir: Path, output_name: str) -> None:
    """Configure the script logger using the local logging.conf file."""
    global logger  # pylint:disable=global-statement,invalid-name

    metadata_dir = Path(__file__).parent
    log_file_path = output_dir / f"{output_name}_metadata.log"
    log_ini_path = metadata_dir / "logging.conf"
    logging.config.fileConfig(
        log_ini_path,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    logger = logging.getLogger(__name__)
    logger.disabled = False


def get_core_metadata(  # pylint: disable=too-many-arguments
    db_name: str,
    host: str,
    port: int,
    user: str,
    password: str,
    production_name: Optional[str],
) -> Tuple[int, Dict[str, str]]:
    """Fetch existing core meta values used to decide INSERT/UPDATE output."""
    core_dict: Dict[str, str] = {}

    if production_name:
        core_dict["species.production_name"] = production_name
        species_query = """
            SELECT species_id
            FROM meta
            WHERE meta_key = 'species.production_name'
              AND meta_value = %s
        """
        species_meta = mysql_fetch_data(
            species_query,
            host=host,
            user=user,
            port=port,
            password=password,
            database=db_name,
            params=(production_name,),
        )
        if not species_meta:
            raise RuntimeError(
                f"Could not find species_id for production name {production_name}"
            )
        species_id = int(species_meta[0]["species_id"])
        print("species ID: " + str(species_id))
    else:
        species_id = 1
        print(
            "WARNING: no production_name provided, using default species ID ="
            + str(species_id)
        )

    core_query = """
        SELECT meta_key, meta_value
        FROM meta
        WHERE species_id = %s
          AND (
            meta_key LIKE 'assembly%%'
            OR meta_key LIKE 'species%%'
            OR meta_key LIKE 'genebuild%%'
            OR meta_key LIKE 'organism%%'
            OR meta_key LIKE 'sample%%'
            OR meta_key LIKE 'annotation%%'
            OR meta_key LIKE 'gencode%%'
          )
    """
    print(core_query)
    core_meta = mysql_fetch_data(
        core_query,
        host=host,
        user=user,
        port=port,
        password=password,
        database=db_name,
        params=(species_id,),
    )
    for meta_pair in core_meta:
        core_dict[normalise_text(meta_pair["meta_key"])] = normalise_text(
            meta_pair["meta_value"]
        )

    return species_id, core_dict


def registry_accession_candidates(
    core_dict: Dict[str, str], assembly_accession: Optional[str]
) -> List[str]:
    """Return accessions to try against the registry, preserving order."""
    candidates = [
        assembly_accession,
        core_dict.get("assembly.accession"),
        core_dict.get("assembly.alt_accession"),
        core_dict.get("assembly.accession_refseq"),
    ]
    seen = set()
    ordered_candidates = []
    for candidate in candidates:
        if candidate and candidate not in seen:
            ordered_candidates.append(candidate)
            seen.add(candidate)
    return ordered_candidates


def fetch_registry_assembly(
    connection: RegistryConnection,
    accession_candidates: Iterable[str],
) -> Tuple[str, Dict[str, Any]]:
    """Fetch assembly, organism, and species registry data for an accession."""
    query = """
        SELECT
            a.assembly_id,
            CONCAT(a.gca_chain, '.', a.gca_version) AS assembly_accession,
            a.lowest_taxon_id AS taxonomy_id,
            a.asm_name AS assembly_name,
            a.asm_level AS assembly_level,
            a.release_date AS assembly_release_date,
            a.refseq_accession AS assembly_refseq_accession,
            a.submitter AS assembly_submitter,
            a.is_current,
            a.asm_type,
            s.species_taxon_id,
            s.scientific_name,
            s.common_name,
            s.parlance_name,
            s.species_prefix,
            s.clade,
            o.organism_id,
            o.biosample_id,
            o.bioproject_id AS organism_bioproject_id,
            o.tol_prefix,
            o.infra_type,
            o.infra_name
        FROM assembly a
        JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
        LEFT JOIN organism o ON a.assembly_id = o.assembly_id
        WHERE CONCAT(a.gca_chain, '.', a.gca_version) = %s
           OR a.refseq_accession = %s
        LIMIT 1
    """

    with connection.cursor(pymysql.cursors.DictCursor) as cursor:
        for accession in accession_candidates:
            cursor.execute(query, (accession, accession))
            registry_row = cursor.fetchone()
            if registry_row:
                return accession, registry_row

    raise RuntimeError(
        "No registry assembly record found for accessions: "
        + ", ".join(accession_candidates)
    )


def fetch_registry_bioprojects(
    connection: RegistryConnection,
    assembly_id: int,
) -> List[Dict[str, Any]]:
    """Fetch all bioproject rows for a registry assembly."""
    query = """
        SELECT
            b.bioproject_id,
            mb.bioproject_name
        FROM bioproject b
        LEFT JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
        WHERE b.assembly_id = %s
        ORDER BY b.bioproject_id
    """
    with connection.cursor(pymysql.cursors.DictCursor) as cursor:
        cursor.execute(query, (assembly_id,))
        return list(cursor.fetchall())


def fetch_registry_genebuilds(
    connection: RegistryConnection,
    registry_gca: str,
    assembly_id: int,
    genebuilder: Optional[str],
) -> Tuple[Optional[Dict[str, Any]], List[Dict[str, Any]]]:
    """Fetch the selected current genebuild and all current genebuild attempts."""
    all_query = """
        SELECT *
        FROM genebuild_status
        WHERE assembly_id = %s
          AND last_attempt = 1
        ORDER BY COALESCE(date_status_update, date_started) DESC,
                 genebuild_status_id DESC
    """
    with connection.cursor(pymysql.cursors.DictCursor) as cursor:
        cursor.execute(all_query, (assembly_id,))
        all_current = list(cursor.fetchall())

    selected_record = None
    if genebuilder:
        helper_record = fetch_current_genebuild_record(
            connection, registry_gca, genebuilder
        )
        if helper_record:
            selected_id = helper_record["genebuild_status_id"]
            selected_record = next(
                (
                    row
                    for row in all_current
                    if row["genebuild_status_id"] == selected_id
                ),
                None,
            )
        if not selected_record:
            logger.warning(
                "No current genebuild_status record found for %s and genebuilder %s",
                registry_gca,
                genebuilder,
            )
    elif all_current:
        selected_record = all_current[0]
        if len(all_current) > 1:
            logger.warning(
                "Multiple current genebuild_status records found for %s; using genebuild_status_id %s. Use --genebuilder to choose explicitly.",
                registry_gca,
                selected_record["genebuild_status_id"],
            )

    return selected_record, all_current


def fetch_registry_metadata(
    registry_config: Dict[str, Any],
    accession_candidates: Iterable[str],
    genebuilder: Optional[str],
) -> Dict[str, Any]:
    """Fetch all registry metadata needed by this script."""
    connection = pymysql.connect(
        database=registry_config["db_name"],
        host=registry_config["db_host"],
        port=int(registry_config["db_port"]),
        user=registry_config["db_user"],
        password=registry_config.get("db_password", ""),
        cursorclass=pymysql.cursors.DictCursor,
    )

    try:
        queried_accession, assembly = fetch_registry_assembly(
            connection, accession_candidates
        )
        bioprojects = fetch_registry_bioprojects(
            connection, int(assembly["assembly_id"])
        )
        selected_genebuild, current_genebuilds = fetch_registry_genebuilds(
            connection,
            normalise_text(assembly["assembly_accession"]),
            int(assembly["assembly_id"]),
            genebuilder,
        )
    finally:
        connection.close()

    return {
        "queried_accession": queried_accession,
        "assembly": assembly,
        "bioprojects": bioprojects,
        "selected_genebuild": selected_genebuild,
        "current_genebuilds": current_genebuilds,
    }


def apply_accession_logic(
    core_dict: Dict[str, str],
    truth_dict: Dict[str, Any],
    registry_assembly: Dict[str, Any],
) -> None:
    """Preserve the old GCA/GCF accession handling."""
    core_accession = core_dict.get("assembly.accession", "")
    registry_gca = normalise_text(registry_assembly.get("assembly_accession"))
    registry_refseq = normalise_text(registry_assembly.get("assembly_refseq_accession"))

    if core_accession.startswith("GCF"):
        truth_dict["assembly.alt_accession"] = (
            core_dict.get("assembly.alt_accession") or registry_gca
        )
        truth_dict["assembly.accession_refseq"] = core_accession
        truth_dict["assembly.accession_body"] = "RefSeq"
    elif "assembly.accession_refseq" in core_dict:
        truth_dict["assembly.alt_accession"] = (
            core_dict.get("assembly.alt_accession") or registry_gca
        )
        core_dict["assembly.accession"] = core_dict["assembly.accession_refseq"]
        truth_dict["assembly.accession_body"] = "RefSeq"
    elif registry_refseq and registry_refseq == core_accession:
        truth_dict["assembly.alt_accession"] = registry_gca
        truth_dict["assembly.accession_refseq"] = registry_refseq
        truth_dict["assembly.accession_body"] = "RefSeq"


def add_registry_assembly_metadata(
    db_name: str,
    registry_assembly: Dict[str, Any],
    truth_dict: Dict[str, Any],
) -> None:
    """Map assembly/organism/species registry columns to core meta keys."""
    truth_dict["assembly.name"] = normalise_text(
        registry_assembly.get("assembly_name")
    ).replace(" ", "_")
    truth_dict["assembly.level"] = normalise_text(
        registry_assembly.get("assembly_level")
    ).lower()
    truth_dict["assembly.date"] = normalise_month(
        registry_assembly.get("assembly_release_date")
    )

    if db_name in BIOSAMPLE_OVERRIDES:
        truth_dict["organism.biosample_id"] = BIOSAMPLE_OVERRIDES[db_name]
    else:
        truth_dict["organism.biosample_id"] = normalise_text(
            registry_assembly.get("biosample_id")
        )

    truth_dict["organism.taxonomy_id"] = normalise_text(
        registry_assembly.get("taxonomy_id")
    )
    truth_dict["organism.species_taxonomy_id"] = normalise_text(
        registry_assembly.get("species_taxon_id")
    )
    truth_dict["organism.common_name"] = capitalise_first(
        registry_assembly.get("common_name")
    )
    truth_dict["organism.scientific_name"] = capitalise_first(
        registry_assembly.get("scientific_name")
    )

    strain = normalise_text(registry_assembly.get("infra_name"))
    if strain == "Caucasian":
        strain = "European"
    truth_dict["organism.strain"] = strain
    truth_dict["organism.strain_type"] = normalise_text(
        registry_assembly.get("infra_type")
    )
    truth_dict["assembly.tol_id"] = normalise_text(registry_assembly.get("tol_prefix"))

    # The current registry has no UCSC alias column. Keep the old "no action when
    # absent" behavior by carrying an empty retrieved value.
    truth_dict["assembly.ucsc_alias"] = ""


def add_static_metadata(
    metadata_dir: Path,
    registry_assembly: Dict[str, Any],
    truth_dict: Dict[str, Any],
) -> Dict[str, Dict[str, str]]:
    """Add static-file values that are still not represented in the registry."""
    provider_static_file = metadata_dir / "provider_static.txt"
    ref_static_file = metadata_dir / "ref_static.txt"
    snp_static_file = metadata_dir / "snp_static.txt"
    url_static_file = metadata_dir / "url_static.txt"

    snp_values = read_two_column_static_file(snp_static_file)
    url_values = read_two_column_static_file(url_static_file)
    ref_values = read_two_column_static_file(ref_static_file)
    provider_values = read_provider_static_file(provider_static_file)

    scientific_name = truth_dict.get("organism.scientific_name", "")
    registry_parlance = normalise_text(registry_assembly.get("parlance_name"))
    if registry_parlance:
        truth_dict["organism.scientific_parlance_name"] = registry_parlance
    elif scientific_name in snp_values:
        truth_dict["organism.scientific_parlance_name"] = snp_values[scientific_name]

    registry_gca = normalise_text(registry_assembly.get("assembly_accession"))
    if registry_gca in url_values:
        truth_dict["assembly.url_name"] = url_values[registry_gca]

    if scientific_name in ref_values and registry_gca == ref_values[scientific_name]:
        truth_dict["assembly.is_reference"] = "1"

    return provider_values


def provider_from_static(
    provider_values: Dict[str, Dict[str, str]], production_name: str
) -> Tuple[str, str]:
    """Return provider fallback values for a production name."""
    provider = provider_values.get(production_name)
    if provider:
        return provider.get("name", ""), provider.get("url", "")
    return "", ""


def add_assembly_provider_metadata(
    core_dict: Dict[str, str], truth_dict: Dict[str, Any]
) -> None:
    """Add assembly provider defaults when missing from the core DB."""
    if "assembly.provider_name" not in core_dict:
        truth_dict["assembly.provider_name"] = DEFAULT_ASSEMBLY_PROVIDER_NAME
        truth_dict["assembly.provider_url"] = DEFAULT_ASSEMBLY_PROVIDER_URL
        logger.warning(
            " | ASSEMBLY_PROVIDER | No assembly provider information found, using default, ENA."
        )


def add_import_provider_metadata(
    core_dict: Dict[str, str],
    truth_dict: Dict[str, Any],
    provider_values: Dict[str, Dict[str, str]],
) -> None:
    """Apply old import-provider fallback rules."""
    production_name = core_dict.get("species.production_name", "")

    if (
        "annotation.provider_name" in core_dict
        and "genebuild.provider_name" not in core_dict
    ):
        truth_dict["genebuild.provider_name"] = core_dict["annotation.provider_name"]
    elif "genebuild.provider_name" not in core_dict:
        provider_name, _provider_url = provider_from_static(
            provider_values, production_name
        )
        if provider_name:
            truth_dict["genebuild.provider_name"] = provider_name
        else:
            logger.critical(
                " | GENEBUILD_PROVIDER_NAME | No genebuild.provider_name could be found either in the core db or in provider_static.txt, this is a required key!"
            )

    if (
        "annotation.provider_url" in core_dict
        and "genebuild.provider_url" not in core_dict
    ):
        truth_dict["genebuild.provider_url"] = core_dict["annotation.provider_url"]
    elif "genebuild.provider_url" not in core_dict:
        _provider_name, provider_url = provider_from_static(
            provider_values, production_name
        )
        if provider_url:
            truth_dict["genebuild.provider_url"] = provider_url
        else:
            logger.critical(
                " | GENEBUILD_PROVIDER_URL | No genebuild.provider_url could be found either in the core db or in provider_static.txt, this is a required key!"
            )


def add_genebuild_method_defaults(
    method: str,
    registry_source: str,
    core_dict: Dict[str, str],
    truth_dict: Dict[str, Any],
    provider_values: Dict[str, Dict[str, str]],
) -> None:
    """Map registry genebuild method/source to old web-facing meta values."""
    scientific_name = truth_dict.get("organism.scientific_name", "").lower()

    if method == "manual_annotation":
        truth_dict.setdefault("genebuild.method_display", "Manual annotation")
        if registry_source:
            truth_dict.setdefault("genebuild.annotation_source", registry_source)

    elif method == "braker":
        truth_dict.setdefault("genebuild.version", "BRK01")
        truth_dict["genebuild.method_display"] = "BRAKER2"
        truth_dict["genebuild.annotation_source"] = "braker"
        truth_dict["genebuild.provider_name"] = DEFAULT_ENSEMBL_PROVIDER_NAME
        truth_dict["genebuild.provider_url"] = DEFAULT_BRAKER_URL

    elif method == "helixer":
        truth_dict.setdefault("genebuild.version", "HLX01")
        truth_dict["genebuild.method_display"] = "Helixer"
        truth_dict["genebuild.annotation_source"] = "helixer"
        truth_dict["genebuild.provider_name"] = DEFAULT_ENSEMBL_PROVIDER_NAME
        truth_dict["genebuild.provider_url"] = DEFAULT_HELIXER_URL

    elif method == "full_genebuild":
        truth_dict.setdefault("genebuild.version", "ENS01")
        truth_dict["genebuild.method_display"] = "Ensembl Genebuild"
        truth_dict["genebuild.annotation_source"] = registry_source or "ensembl"
        truth_dict["genebuild.provider_name"] = DEFAULT_ENSEMBL_PROVIDER_NAME
        truth_dict["genebuild.provider_url"] = DEFAULT_VERTEBRATE_GENEBUILD_URL

    elif method == "anno":
        truth_dict.setdefault("genebuild.version", "ENS01")
        truth_dict["genebuild.method_display"] = "Ensembl Genebuild"
        truth_dict["genebuild.annotation_source"] = registry_source or "ensembl"
        truth_dict["genebuild.provider_name"] = DEFAULT_ENSEMBL_PROVIDER_NAME
        truth_dict["genebuild.provider_url"] = DEFAULT_NONVERTEBRATE_GENEBUILD_URL

    elif method == "projection_build":
        truth_dict.setdefault("genebuild.version", "ENS01")
        truth_dict["genebuild.annotation_source"] = registry_source or "ensembl"
        truth_dict["genebuild.provider_name"] = DEFAULT_ENSEMBL_PROVIDER_NAME
        if scientific_name == "homo sapiens":
            truth_dict["genebuild.provider_url"] = DEFAULT_HUMAN_MAPPING_URL
            truth_dict["genebuild.method_display"] = "Mapping from GRCh38"
        else:
            truth_dict["genebuild.provider_url"] = DEFAULT_VERTEBRATE_GENEBUILD_URL
            truth_dict["genebuild.method_display"] = "Mapping from reference"

    elif method in {"import", "external_annotation_import"}:
        truth_dict.setdefault("genebuild.method_display", "Import")
        truth_dict.setdefault("genebuild.version", "EXT01")
        if registry_source:
            truth_dict["genebuild.annotation_source"] = registry_source
        elif (
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
        add_import_provider_metadata(core_dict, truth_dict, provider_values)

    elif method == "mixed_strategy_build":
        truth_dict.setdefault("genebuild.version", "ENS01")
        truth_dict["genebuild.method_display"] = "Mixed strategy build"
        truth_dict["genebuild.annotation_source"] = registry_source or "ensembl"
        truth_dict["genebuild.provider_name"] = DEFAULT_ENSEMBL_PROVIDER_NAME
        truth_dict["genebuild.provider_url"] = DEFAULT_VERTEBRATE_GENEBUILD_URL

    elif method == "pending":
        logger.warning(
            " | GENEBUILD_METHOD | Registry genebuild method is pending; provider/display metadata may need manual review."
        )
        if registry_source:
            truth_dict["genebuild.annotation_source"] = registry_source


def add_genebuild_metadata(
    selected_genebuild: Optional[Dict[str, Any]],
    core_dict: Dict[str, str],
    truth_dict: Dict[str, Any],
    provider_values: Dict[str, Dict[str, str]],
) -> None:
    """Add genebuild metadata from registry, with old core/default fallbacks."""
    method = ""
    registry_source = ""
    if selected_genebuild:
        method = normalise_text(selected_genebuild.get("annotation_method"))
        registry_source = normalise_text(selected_genebuild.get("annotation_source"))
        registry_version = normalise_text(selected_genebuild.get("genebuild_version"))

        if method:
            truth_dict["genebuild.method"] = method
        if registry_version:
            truth_dict["genebuild.version"] = registry_version
        if selected_genebuild.get("date_started"):
            truth_dict["genebuild.start_date"] = normalise_month(
                selected_genebuild["date_started"]
            )
        if selected_genebuild.get("last_genebuild_update"):
            truth_dict["genebuild.last_geneset_update"] = normalise_month(
                selected_genebuild["last_genebuild_update"]
            )
        if selected_genebuild.get("release_date"):
            truth_dict["genebuild.initial_release_date"] = normalise_month(
                selected_genebuild["release_date"]
            )

        add_genebuild_method_defaults(
            method, registry_source, core_dict, truth_dict, provider_values
        )

    if "genebuild.method_display" not in truth_dict:
        if "gencode.version" in core_dict:
            truth_dict["genebuild.version"] = core_dict["gencode.version"].replace(
                " ", ""
            )
            truth_dict.setdefault("genebuild.method", "manual_annotation")
            truth_dict["genebuild.method_display"] = "Manual annotation"
        elif "genebuild.method" in core_dict:
            method = core_dict["genebuild.method"]
            truth_dict["genebuild.method"] = method
            add_genebuild_method_defaults(
                method, registry_source, core_dict, truth_dict, provider_values
            )
        else:
            logger.warning(
                " | GENEBUILD_METHOD | No registry/core genebuild.method, assuming this is an Ensembl annotation"
            )
            truth_dict["genebuild.version"] = "ENS01"
            truth_dict["genebuild.method"] = "full_genebuild"
            truth_dict["genebuild.annotation_source"] = "ensembl"
            truth_dict["genebuild.provider_name"] = DEFAULT_ENSEMBL_PROVIDER_NAME
            truth_dict["genebuild.provider_url"] = DEFAULT_VERTEBRATE_GENEBUILD_URL
            truth_dict["genebuild.method_display"] = "Ensembl Genebuild"

    if "genebuild.version" not in truth_dict and core_dict.get("genebuild.version"):
        truth_dict["genebuild.version"] = core_dict["genebuild.version"]
    if "genebuild.start_date" not in truth_dict and core_dict.get(
        "genebuild.start_date"
    ):
        truth_dict["genebuild.start_date"] = core_dict["genebuild.start_date"]
    if "genebuild.method" not in truth_dict and core_dict.get("genebuild.method"):
        truth_dict["genebuild.method"] = core_dict["genebuild.method"]


def add_core_fallback_metadata(
    core_dict: Dict[str, str], truth_dict: Dict[str, Any], team: str
) -> None:
    """Add metadata that is still sourced from core values or script arguments."""
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

    if "genebuild.last_geneset_update" not in truth_dict:
        if core_dict.get("genebuild.last_geneset_update"):
            truth_dict["genebuild.last_geneset_update"] = core_dict[
                "genebuild.last_geneset_update"
            ]
        else:
            logger.warning(
                " | GENEBUILD.LAST_GENESET_UPDATE | No genebuild.last_geneset_update could be found in the registry or core db, this is a required key, I'm setting it to today."
            )
            truth_dict["genebuild.last_geneset_update"] = datetime.now().strftime(
                "%Y-%m"
            )

    try:
        truth_dict["organism.production_name"] = core_dict["species.production_name"]
    except KeyError:
        logger.critical(
            " | SPECIES.PRODUCTION_NAME | No species.production_name could be found in the core db, this is a required key!"
        )

    truth_dict["genebuild.team_responsible"] = team


def build_truth_metadata(
    db_name: str,
    core_dict: Dict[str, str],
    registry_payload: Dict[str, Any],
    metadata_dir: Path,
    team: str,
) -> Dict[str, Any]:
    """Build the full metadata dictionary used for SQL and JSON output."""
    truth_dict: Dict[str, Any] = {}
    registry_assembly = registry_payload["assembly"]

    apply_accession_logic(core_dict, truth_dict, registry_assembly)
    add_registry_assembly_metadata(db_name, registry_assembly, truth_dict)
    provider_values = add_static_metadata(metadata_dir, registry_assembly, truth_dict)
    add_assembly_provider_metadata(core_dict, truth_dict)
    add_genebuild_metadata(
        registry_payload.get("selected_genebuild"),
        core_dict,
        truth_dict,
        provider_values,
    )
    add_core_fallback_metadata(core_dict, truth_dict, team)

    return truth_dict


def build_sql_actions(
    truth_dict: Dict[str, Any], core_dict: Dict[str, str], species_id: int
) -> List[Dict[str, str]]:
    """Create SQL patch statements with the old insert/update/delete behavior."""
    actions: List[Dict[str, str]] = []

    for meta_key, raw_value in truth_dict.items():
        meta_value = normalise_text(raw_value)
        escaped_value = sql_escape(meta_value)

        if meta_key not in core_dict:
            if meta_value:
                sql = (
                    "INSERT IGNORE INTO meta (species_id, meta_key, meta_value) "
                    f"VALUES({species_id}, '{meta_key}', '{escaped_value}');"
                )
                actions.append(
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
            actions.append(
                {
                    "action": "update",
                    "meta_key": meta_key,
                    "meta_value": meta_value,
                    "sql": sql,
                }
            )

    for remove_key in META_KEYS_TO_REMOVE:
        actions.append(
            {
                "action": "delete",
                "meta_key": remove_key,
                "meta_value": "",
                "sql": f"DELETE from meta WHERE meta_key='{remove_key}';",
            }
        )

    return actions


def check_required_metadata(
    core_dict: Dict[str, str], truth_dict: Dict[str, Any]
) -> None:
    """Log missing required metadata keys."""
    for required_key in REQUIRED_META_KEYS:
        if required_key not in core_dict and not normalise_text(
            truth_dict.get(required_key)
        ):
            logger.critical(f"You are missing required meta key: {required_key}")

    if "species.common_name" in core_dict:
        organism_common_name = normalise_text(truth_dict.get("organism.common_name"))
        if core_dict["species.common_name"].lower() != organism_common_name.lower():
            logger.warning(
                ' | COMMON NAME | The value for species.common_name in your meta table: "'
                + core_dict["species.common_name"]
                + '" does not match the value that I am assigning to organism.common_name: "'
                + organism_common_name
                + '"'
            )


def write_sql(sql_path: Path, db_name: str, sql_actions: List[Dict[str, str]]) -> None:
    """Write the SQL patch file."""
    with open(sql_path, "w", encoding="utf-8") as sql_out:
        print(f"USE {db_name};", file=sql_out)
        for action in sql_actions:
            print(action["sql"], file=sql_out)


def write_metadata_json(  # pylint: disable=too-many-arguments
    json_path: Path,
    db_name: str,
    species_id: int,
    registry_config: Dict[str, Any],
    registry_payload: Dict[str, Any],
    core_dict: Dict[str, str],
    truth_dict: Dict[str, Any],
    sql_actions: List[Dict[str, str]],
) -> None:
    """Write a JSON artifact with all metadata retrieved and generated."""
    safe_registry_config = {
        "db_host": registry_config.get("db_host"),
        "db_port": registry_config.get("db_port"),
        "db_name": registry_config.get("db_name"),
        "db_user": registry_config.get("db_user"),
    }
    payload = {
        "database": db_name,
        "species_id": species_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "registry_connection": safe_registry_config,
        "registry_metadata": registry_payload,
        "core_metadata": core_dict,
        "metadata": {key: normalise_text(value) for key, value in truth_dict.items()},
        "sql_actions": sql_actions,
    }

    with open(json_path, "w", encoding="utf-8") as json_out:
        json.dump(serialise_json(payload), json_out, indent=2, sort_keys=True)
        json_out.write("\n")


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the command line parser."""
    parser = argparse.ArgumentParser(
        description="Prepare SQL updates for core dbs using the assembly registry"
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default=".",
        help="Path where the output and temp files will write to. Uses current dir by default",
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
        help="Core database host server",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--port",
        help="Core database host server port",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--db_user",
        default="ensro",
        help="Core database user",
    )
    parser.add_argument(
        "--db_password",
        default="",
        help="Core database password",
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
        "-a",
        "--assembly_accession",
        help="Assembly accession to query in the registry. Defaults to assembly.accession from the core DB.",
    )
    parser.add_argument(
        "--registry_host",
        default=os.environ.get("GBS1"),
        help="Registry database host. Defaults to GBS1.",
    )
    parser.add_argument(
        "--registry_port",
        default=os.environ.get("GBP1"),
        type=int,
        help="Registry database port. Defaults to GBP1.",
    )
    parser.add_argument(
        "--registry_user",
        default="ensro",
        help="Registry database user",
    )
    parser.add_argument(
        "--registry_password",
        default="",
        help="Registry database password",
    )
    parser.add_argument(
        "--registry_db",
        default=DEFAULT_REGISTRY_DB,
        help=f"Registry database name. Defaults to {DEFAULT_REGISTRY_DB}.",
    )
    parser.add_argument(
        "--genebuilder",
        help="Optional registry genebuilder to use when multiple current genebuild records exist.",
    )
    parser.add_argument(
        "--json_output",
        help="Path for JSON metadata output. Defaults to <output_name>.json in output_dir.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output (check that all required keys are not NULL/ empty)",
    )
    return parser


def main() -> None:
    """Script entry point."""
    parser = build_arg_parser()
    args = parser.parse_args()

    if not args.registry_host or not args.registry_port:
        raise RuntimeError(
            "Registry connection details are required. Provide --registry_host/--registry_port or set GBS1/GBP1."
        )

    db_name = args.db_name
    output_name = args.production_name or db_name
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir, output_name)

    sql_path = output_dir / f"{output_name}_registry.sql"
    json_path = (
        Path(args.json_output)
        if args.json_output
        else output_dir / f"{output_name}_registry.json"
    )

    print(f"Working on database: {db_name}")

    species_id, core_dict = get_core_metadata(
        db_name=db_name,
        host=args.host,
        port=args.port,
        user=args.db_user,
        password=args.db_password,
        production_name=args.production_name,
    )

    registry_config = {
        "db_host": args.registry_host,
        "db_port": int(args.registry_port),
        "db_user": args.registry_user,
        "db_password": args.registry_password,
        "db_name": args.registry_db,
    }
    accession_candidates = registry_accession_candidates(
        core_dict, args.assembly_accession
    )
    if not accession_candidates:
        raise RuntimeError(
            "No assembly accession found. Provide --assembly_accession or ensure the core DB has assembly.accession."
        )

    registry_payload = fetch_registry_metadata(
        registry_config, accession_candidates, args.genebuilder
    )
    metadata_dir = Path(__file__).parent
    truth_dict = build_truth_metadata(
        db_name, core_dict, registry_payload, metadata_dir, args.team
    )
    check_required_metadata(core_dict, truth_dict)

    sql_actions = build_sql_actions(truth_dict, core_dict, species_id)
    write_sql(sql_path, db_name, sql_actions)
    write_metadata_json(
        json_path,
        db_name,
        species_id,
        registry_config,
        registry_payload,
        core_dict,
        truth_dict,
        sql_actions,
    )

    if args.verbose:
        print(
            "\n".join(
                f"{key:<35}   {normalise_text(val):<25}   {'required' if key in REQUIRED_META_KEYS else ''}"
                for key, val in truth_dict.items()
            )
        )

    print(f"Wrote SQL patch: {sql_path}")
    print(f"Wrote metadata JSON: {json_path}")


if __name__ == "__main__":
    main()

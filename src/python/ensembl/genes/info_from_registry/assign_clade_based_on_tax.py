import os

import pymysql
import json
from mysql_helper import mysql_fetch_data
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline_setup.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def create_tax_dictionary_from_registry(server_info, registry_info):
    """Get taxonomy info for the single taxon_id from the registry_info dict."""

    taxon_id = registry_info["taxon_id"]

    taxonomy_query = """
        SELECT lowest_taxon_id, taxon_class_id, taxon_class
        FROM taxonomy
        WHERE lowest_taxon_id = %s
        ORDER BY FIELD(taxon_class, 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom');
    """

    try:
        taxonomy_info = mysql_fetch_data(
            taxonomy_query,
            host=server_info["registry"]["db_host"],
            user=server_info["registry"]["db_user"],
            port=server_info["registry"]["db_port"],
            database=server_info["registry"]["registry_db"],
            password= "",
            params=taxon_id,
        )

        taxonomy_dict = {}
        for row in taxonomy_info:
            lowest_taxon_id = row["lowest_taxon_id"]
            if lowest_taxon_id not in taxonomy_dict:
                taxonomy_dict[lowest_taxon_id] = []
            taxonomy_dict[lowest_taxon_id].append(
                {
                    "taxon_class_id": row["taxon_class_id"],
                    "taxon_class": row["taxon_class"],
                }
            )

        logger.info(f"Found taxonomy for {taxon_id}")
        return taxonomy_dict

    except pymysql.Error as err:
        logger.error("Error while fetching taxonomy info: %s", err)
        return {}


def load_clade_data():
    """Hardcoded path for clade settings."""
    json_file = os.path.join(
        os.environ.get("ENSCODE"),
        "ensembl-genes",
        "src",
        "python",
        "ensembl",
        "genes",
        "info_from_registry",
        "clade_settings.json"
    )

    with open(json_file, "r") as f:
        logging.info("Loading clade settings json file.")
        return json.load(f)


def assign_clade(server_info, registry_info):
    """
    Given a lowest taxon ID, assign the clade and return clade details and genus_taxon_id.
    """
    clade_data = load_clade_data()
    taxonomy_dict = create_tax_dictionary_from_registry(server_info, registry_info)

    # Accept both int and str keys
    lowest_taxon_id = registry_info["taxon_id"]
    taxonomy_hierarchy = taxonomy_dict.get(str(lowest_taxon_id)) or taxonomy_dict.get(lowest_taxon_id, [])

    if not taxonomy_hierarchy:
        logging.warning(f"Taxonomy hierarchy not found for taxon ID {registry_info['taxon_id']}")
        return "Unassigned", None, None

    internal_clade = "Unassigned"
    genus_taxon_id = None
    clade_details = None

    # Get genus taxon_id
    for taxon in taxonomy_hierarchy:
        if taxon["taxon_class"] == "genus":
            genus_taxon_id = taxon["taxon_class_id"]

    # Step 1: try exact match on the lowest_taxon_id
    for clade_name, details in clade_data.items():
        clade_taxon_id = int(details.get("taxon_id", -1))
        if clade_taxon_id == lowest_taxon_id:
            internal_clade = clade_name
            clade_details = {k: v for k, v in details.items() if k != "taxon_id"}
            logging.info(f"Exact match: Assigned clade '{internal_clade}' for taxon {registry_info['taxon_id']}")
            return internal_clade, genus_taxon_id, clade_details

    # Step 2: Walk up the taxonomy hierarchy
    taxon_classes_order = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

    for taxon_class in taxon_classes_order:
        matching = next((t for t in taxonomy_hierarchy if t["taxon_class"] == taxon_class), None)
        if not matching:
            continue
        current_taxon_id = int(matching["taxon_class_id"])

        for clade_name, details in clade_data.items():
            clade_taxon_id = int(details.get("taxon_id", -1))
            if clade_taxon_id == current_taxon_id:
                internal_clade = clade_name
                clade_details = {k: v for k, v in details.items() if k != "taxon_id"}
                logging.info(f"Hierarchy match: Assigned clade '{internal_clade}' via {taxon_class} taxon_id {current_taxon_id}")
                return internal_clade, genus_taxon_id, clade_details

    # No match found
    logging.error(f"No clade found for taxon {registry_info['taxon_id']} in full hierarchy.")
    return "Unassigned", genus_taxon_id, None


def assign_clade_info_custom_loading(registry_info):
    """
    Look for a specific clade in the JSON data and return all values except for taxon_id.

    Args:
            registry_info (dict): Dictionary containing clade information with 'clade' key

    Returns:
            dict or None: Dictionary containing all clade details except taxon_id,
                                     or None if clade not found
    """
    clade_data = load_clade_data()

    # Extract the clade name from the dictionary
    clade_name = registry_info.get("clade")

    if not clade_name:
        logging.error("No 'clade' key found in provided dictionary")
        return None

    # Look for the specific clade in the dictionary
    if clade_name in clade_data:
        details = clade_data[clade_name]
        # Return all details except taxon_id
        clade_details = {k: v for k, v in details.items() if k != "taxon_id"}
        return clade_details
    else:
        logging.warning(f"Clade '{clade_name}' not found in clade data")
        return None

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
# pylint: disable=logging-fstring-interpolation
"""This module handles the assignment of unique species prefixes
based on taxon IDs by interacting with the assembly registry and
metadata databases.
"""
import json
import logging
import os
import random
import string
from typing import Optional
import pymysql  # type: ignore

from ensembl.genes.info_from_registry.mysql_helper import mysql_fetch_data


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("pipeline_setup.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def get_special_cases() -> dict[str, str]:
    """Retrieve special cases for species prefixes from a JSON file.
    Returns:
        dict: A dictionary mapping taxon IDs to their special species prefixes.
    """
    enscode = os.environ.get("ENSCODE")
    if not enscode:
        raise EnvironmentError("Environment variable ENSCODE is not set")
    json_path = os.path.join(
        enscode,
        "ensembl-genes",
        "src",
        "python",
        "ensembl",
        "genes",
        "info_from_registry",
        "anno_settings.json",
    )
    with open(json_path, "r") as file:  # pylint: disable=unspecified-encoding
        special_cases = json.load(file)
    return special_cases


def exiting_prefix(server_info: dict) -> list[str]:
    """Get a list of existing species prefixes from the gb assembly registry and metadata databases.

    Args:
        server_info (dict): Information about the database server.

    Returns:
        list[str]: A list of existing species prefixes.
    """
    # Getting existing prefix from registry db. To be removed when the registry is updated.
    prefix_registry_query = f"SELECT DISTINCT species_prefix FROM assembly ;"  # pylint: disable=f-string-without-interpolation
    output_registry = mysql_fetch_data(
        prefix_registry_query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database=server_info["registry"]["db_name"],
        password="",
    )
    # Getting existing prefix from metadata db
    prefix_metadata_query = f"SELECT DISTINCT prefix FROM species_prefix ;"  # pylint: disable=f-string-without-interpolation
    output_metadata = mysql_fetch_data(
        prefix_metadata_query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database=server_info["registry"]["db_name"],
        password="",
    )
    prefix_list = [
        list(item.values())[0] for item in list(output_registry) + list(output_metadata)
    ]
    existing_prefix = list(set(prefix_list))
    logger.debug(f"Num Existing prefix: {len(existing_prefix)}")
    return existing_prefix


def generate_random_prefix(existing_prefix: list[str]) -> str:
    """Generate a random species prefix that does not already exist in the provided list.

    Args:
        existing_prefix (list[str]): A list of existing species prefixes.

    Returns:
        str: A random three/four letter prefix.
    """
    letters = string.ascii_uppercase
    if len(existing_prefix) >= 26**3:
        logger.info("Creating prefix: four random letters")
        length = 4
    else:
        logger.info("Creating prefix: three random letters")
        length = 3
    while True:
        candidate = "ENS" + "".join(random.choices(letters, k=length))
        if candidate not in existing_prefix:
            return candidate


def insert_prefix_into_db(
    prefix: str,
    taxon_id: int,
    conn: pymysql.connections.Connection,
    store_new_registry: bool = False,
) -> bool:
    """Insert a new species prefix into the database.
    Args:
        prefix (str): The species prefix to insert.
        taxon_id (int): The lowest taxon ID associated with the prefix.
        conn (pymysql.connections.Connection): The database connection object.
        store_new_registry (bool): Whether to allow duplicated exception if
        the prefix already exists in the database.
    Returns:
        bool: True if the prefix was successfully inserted, False if it already exists.
    Raises:
        pymysql.err.IntegrityError: If there is a duplicate entry for the prefix.
    """
    logger.info(f"Inserting prefix {prefix} for taxon id {taxon_id}")
    try:
        with conn.cursor() as cursor:
            query = (
                "INSERT INTO species_prefix (lowest_taxon_id, prefix) "
                "VALUES (%s, %s);"
            )
            cursor.execute(query, (taxon_id, prefix))
        return True
    except pymysql.err.IntegrityError as err:
        if err.args[0] == 1062:
            if store_new_registry:
                logger.info(
                    f"Existing prefix {prefix} already stored in the new registry. Allow exception"
                )
                return True
            return False
        raise


def create_prefix(existing_prefix: list[str], taxon_id: int, server_info: dict) -> str:
    """Create a new species prefix for a given taxon ID and insert it into the database.
    It will retry up to 10,000 times to ensure uniqueness.

    Args:
        existing_prefix (list[str]): A list of existing species prefixes.
        taxon_id (int): The taxon ID for the new prefix.
        server_info (dict): Information about the database server.

    Raises:
        RuntimeError: If a unique prefix cannot be generated after many attempts.

    Returns:
        str: The newly created species prefix.
    """
    logger.info(f"Creating new prefix for taxon ID: {taxon_id}")
    conn = pymysql.connect(
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user_w"],
        port=server_info["registry"]["db_port"],
        password=server_info["registry"]["password"],
        database=server_info["registry"]["db_name"],
    )

    with conn:
        for _ in range(10):  # max attempts
            prefix = generate_random_prefix(existing_prefix)
            if insert_prefix_into_db(prefix, taxon_id, conn):
                logger.info(f"Successfully inserted: {prefix}")
                return prefix
            existing_prefix.append(prefix)  # optimize by adding to "existing prefixes"
    raise RuntimeError("Failed to generate unique prefix after many attempts.")


def get_species_prefix(taxon_id: int, server_info: dict) -> Optional[str]:
    """
    This function retrieves the species prefix from the assembly registry
    and metadata databases.
    If the prefix is not found, it creates a new one. There are special cases
    where the prefix is predefined.
    - Canis lupus (wolf) -> ENSCAF
    - Canis lupus familiaris (Domestic dog) -> ENSCAF
    - Heterocephalus glaber (naked mole rat) -> ENSHGL

    Args:
        taxon_id (int): lowest taxon id

    Raises:
        ValueError: If the taxon ID is not found.

    Returns:
        str: unique species prefix that already exist or a new one when no prefix is found.
    """

    # Special cases
    special_cases = get_special_cases()

    if str(taxon_id) in special_cases:
        logger.info(f"The prefix is a special case for taxon ID: {taxon_id}")
        species_prefix = special_cases.get(str(taxon_id))
    else:

        logger.info(f"Searching for prefix for taxon ID: {taxon_id}")

        prefix_registry_query = (
            f"SELECT DISTINCT species_prefix FROM assembly WHERE taxonomy = {taxon_id}"
        )
        output_registry = mysql_fetch_data(
            prefix_registry_query,
            host=server_info["registry"]["db_host"],
            user=server_info["registry"]["db_user"],
            port=server_info["registry"]["db_port"],
            database="gb_assembly_registry",
            password="",
        )
        if output_registry:
            logger.info(f"Prefix found in old registry: {output_registry}")

        prefix_metadata_query = f"SELECT DISTINCT prefix FROM species_prefix \
            WHERE lowest_taxon_id = {taxon_id}"
        output_metadata = mysql_fetch_data(
            prefix_metadata_query,
            host=server_info["registry"]["db_host"],
            user=server_info["registry"]["db_user"],
            port=server_info["registry"]["db_port"],
            database=server_info["registry"]["db_name"],
            password="",
        )
        if output_metadata:
            logger.info(f"Prefix found in new metadata registry: {output_metadata}")

        # Combine output and get list of unique values
        output = [
            list(item.values())[0]
            for item in list(output_registry) + list(output_metadata)
        ]
        prefix_list = list(set(output))

        # no prefix, create new prefix
        if len(prefix_list) == 0:
            logger.info(f"Creating a new prefix for taxon ID: {taxon_id}")
            existing_prefix = exiting_prefix(server_info)
            species_prefix = create_prefix(existing_prefix, taxon_id, server_info)

        # unique prefix detected
        elif len(prefix_list) == 1:
            logger.info(f"Unique prefix detected for taxon ID: {taxon_id}")
            species_prefix = str(prefix_list[0])
            logger.info(f"Saving exting prefix {species_prefix} in new registry")
            # Open a new connection
            conn = pymysql.connect(
                host=server_info["registry"]["db_host"],
                user=server_info["registry"]["db_user_w"],
                port=server_info["registry"]["db_port"],
                password=server_info["registry"]["password"],
                database=server_info["registry"]["db_name"],
            )

            if insert_prefix_into_db(
                species_prefix, taxon_id, conn, store_new_registry=True
            ):
                logger.info(f"Successfully inserted: {species_prefix}")
                conn.close()

        else:
            raise ValueError(
                f"The taxon {taxon_id} is already registered and multiple prefix \
                    were detected: {prefix_list}"
            )

    return species_prefix

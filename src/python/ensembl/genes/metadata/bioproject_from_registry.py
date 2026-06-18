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

import logging
import os
from typing import List, Optional

import pymysql

logger = logging.getLogger(__name__)

DEFAULT_REGISTRY_DB = "gb_assembly_metadata"


def normalise_genome_group_name(genome_group_name: str) -> str:
    """Return the meta-safe genome group name."""
    return genome_group_name.strip().lower().replace("/", "-")


def get_bioproject_names(
    assembly_accession: str,
    user: str,
    host: Optional[str] = None,
    port: Optional[int] = None,
    database: str = DEFAULT_REGISTRY_DB,
) -> List[str]:
    """Return main BioProject genome group names for an assembly accession.

    Args:
        assembly_accession: GCA assembly accession including version.
        user: Registry database user.
        host: Registry host. Defaults to the GBS1 environment variable.
        port: Registry port. Defaults to the GBP1 environment variable.
        database: Registry database name.

    Returns:
        Main BioProject genome group names, or an empty list if they are unavailable.
    """
    registry_host = host or os.environ.get("GBS1")
    registry_port = port or os.environ.get("GBP1")

    if not user:
        logger.warning(" | GENOME.GENOME_GROUP | Registry user not set")
        return []

    if not registry_host or not registry_port:
        logger.warning(
            " | GENOME.GENOME_GROUP | Registry host/port not set; expected GBS1 and GBP1"
        )
        return []

    try:
        registry_port = int(registry_port)
    except (TypeError, ValueError):
        logger.warning(
            " | GENOME.GENOME_GROUP | Invalid registry port for BioProject lookup: %s",
            registry_port,
        )
        return []

    query = """
        SELECT DISTINCT mb.bioproject_name AS genome_group
        FROM assembly a
        JOIN bioproject b ON a.assembly_id = b.assembly_id
        JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
        WHERE CONCAT(a.gca_chain, '.', a.gca_version) = %s
          AND mb.bioproject_name IS NOT NULL
        ORDER BY genome_group
    """

    conn = None
    try:
        conn = pymysql.connect(
            host=registry_host,
            user=user,
            port=registry_port,
            database=database,
            cursorclass=pymysql.cursors.DictCursor,
        )
        with conn.cursor() as cursor:
            cursor.execute(query, (assembly_accession,))
            results = cursor.fetchall()
    except pymysql.Error:
        logger.exception(
            " | GENOME.GENOME_GROUP | Could not fetch genome group names from registry for %s",
            assembly_accession,
        )
        return []
    finally:
        if conn is not None:
            try:
                conn.close()
            except pymysql.Error:
                logger.exception("Error closing registry connection")

    registry_genome_group_names = [
        normalise_genome_group_name(str(result["genome_group"]))
        for result in results
        if result.get("genome_group")
    ]

    if not registry_genome_group_names:
        logger.warning(
            " | GENOME.GENOME_GROUP | No genome group name found in registry for %s",
            assembly_accession,
        )
        return []

    genome_group_names = [
        genome_group_name
        for genome_group_name in registry_genome_group_names
        if genome_group_name != "hprc"
    ]

    return genome_group_names


def get_bioproject_name(
    assembly_accession: str,
    user: str,
    host: Optional[str] = None,
    port: Optional[int] = None,
    database: str = DEFAULT_REGISTRY_DB,
) -> str:
    """Return the first main BioProject genome group name for an assembly accession."""
    bioproject_names = get_bioproject_names(
        assembly_accession,
        user=user,
        host=host,
        port=port,
        database=database,
    )
    return bioproject_names[0] if bioproject_names else ""

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


def get_bioproject_names(
    assembly_accession: str,
    user: str,
    host: Optional[str] = None,
    port: Optional[int] = None,
    database: str = DEFAULT_REGISTRY_DB,
) -> List[str]:
    """Return the main BioProject names for an assembly accession.

    Args:
        assembly_accession: GCA assembly accession including version.
        user: Registry database user.
        host: Registry host. Defaults to the GBS1 environment variable.
        port: Registry port. Defaults to the GBP1 environment variable.
        database: Registry database name.

    Returns:
        BioProject names, or an empty list if they are unavailable.
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
        SELECT DISTINCT mb.bioproject_name
        FROM assembly a
        JOIN bioproject b ON a.assembly_id = b.assembly_id
        JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
        WHERE CONCAT(a.gca_chain, '.', a.gca_version) = %s
          AND mb.bioproject_name IS NOT NULL
        ORDER BY mb.bioproject_name
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
            " | GENOME.GENOME_GROUP | Could not fetch BioProject name from registry for %s",
            assembly_accession,
        )
        return []
    finally:
        if conn is not None:
            try:
                conn.close()
            except pymysql.Error:
                logger.exception("Error closing registry connection")

    bioproject_names = [
        str(result["bioproject_name"])
        for result in results
        if result.get("bioproject_name")
    ]

    if not bioproject_names:
        logger.warning(
            " | GENOME.GENOME_GROUP | No BioProject name found in registry for %s",
            assembly_accession,
        )
        return []

    return bioproject_names


def get_bioproject_name(
    assembly_accession: str,
    user: str,
    host: Optional[str] = None,
    port: Optional[int] = None,
    database: str = DEFAULT_REGISTRY_DB,
) -> str:
    """Return the first main BioProject name for an assembly accession."""
    bioproject_names = get_bioproject_names(
        assembly_accession,
        user=user,
        host=host,
        port=port,
        database=database,
    )
    return bioproject_names[0] if bioproject_names else ""

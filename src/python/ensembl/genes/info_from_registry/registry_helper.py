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
"""
Registry Helper Module

This module provides utility functions for common registry database queries
that are shared across multiple scripts.
"""

import re
from typing import cast, Optional, Any
import pymysql


def fetch_assembly_id(
    connection: pymysql.connections.Connection, assembly: str
) -> Optional[int]:
    """
    Fetch the assembly ID for a given assembly accession.

    Args:
        connection: MySQL connection object
        assembly (str): Assembly Accession (GCA format like GCA_123456789.1)
    Returns:
        int: Assembly ID if found, else None
    """
    query = """
    SELECT assembly_id FROM assembly
    WHERE CONCAT(gca_chain, '.', gca_version) = %s
    """

    with connection.cursor(pymysql.cursors.DictCursor) as cursor:
        cursor.execute(query, (assembly,))
        result = cursor.fetchone()

    if result is None:
        return None
    result_dict = cast(dict[str, int], result)
    return result_dict["assembly_id"]


def fetch_current_genebuild_record(
    connection: pymysql.connections.Connection,
    assembly: str,
    genebuilder: Optional[str] = None,
) -> Optional[dict[str, Any]]:
    """
    Fetch the full current genebuild record for a given assembly.

    Args:
        connection: MySQL connection object
        assembly (str): Assembly Accession (GCA format)
        genebuilder (str, optional): Genebuilder name to filter by.\
            If None, fetches any active record.
    Returns:
        dict: Record with genebuild_status_id, gb_status, genebuilder,\
            annotation_method, genebuild_version, and other fields if found, else None
    """
    params: tuple[Any, ...] = ()
    if genebuilder:
        query = """
        SELECT genebuild_status_id, gb_status, genebuilder, annotation_method, genebuild_version
        FROM genebuild_status
        WHERE gca_accession = %s AND genebuilder = %s AND last_attempt = 1
        """
        params = (assembly, genebuilder)
    else:
        query = """
        SELECT genebuild_status_id, gb_status, genebuilder, annotation_method, genebuild_version
        FROM genebuild_status
        WHERE gca_accession = %s AND last_attempt = 1
        """
        params = (assembly,)

    with connection.cursor(pymysql.cursors.DictCursor) as cursor:
        cursor.execute(query, params)
        result = cursor.fetchone()

    if result is None:
        return None
    return cast(dict[str, Any], result)


def fetch_genebuild_status_id(
    connection: pymysql.connections.Connection, assembly: str
) -> Optional[int]:
    """
    Fetch the current genebuild status ID for a given assembly.

    Args:
        connection: MySQL connection object
        assembly (str): Assembly Accession (GCA format)
    Returns:
        int: genebuild_status_id if found (where last_attempt=1), else None
    """
    record = fetch_current_genebuild_record(connection, assembly)
    return record["genebuild_status_id"] if record else None


def fetch_highest_genebuild_version(
    connection: pymysql.connections.Connection, assembly: str
) -> Optional[str]:
    """
    Fetch the highest genebuild version for a given assembly across all genebuilders.
    This ensures version uniqueness per assembly, regardless of who owns the genebuild.
    Args:
        connection: MySQL connection object
        assembly (str): Assembly Accession (GCA format)
    Returns:
        str: Highest genebuild version (e.g., "ENS02") if any exists, else None
    """
    query = """
    SELECT genebuild_version
    FROM genebuild_status
    WHERE gca_accession = %s
    ORDER BY genebuild_version DESC
    LIMIT 1
    """

    with connection.cursor(pymysql.cursors.DictCursor) as cursor:
        cursor.execute(query, (assembly,))
        result = cursor.fetchone()

    return result["genebuild_version"] if result else None


def increment_genebuild_version(version: str) -> str:
    """
    Increment a genebuild version string.
    Examples:
        ENS01 -> ENS02
        ENS09 -> ENS10
        BRK01 -> BRK02
        HLX99 -> HLX100
    Args:
        version (str): Current version (e.g., "ENS01")
    Returns:
        str: Incremented version (e.g., "ENS02")
    """

    # Extract prefix and numeric suffix
    match = re.match(r"^([A-Z]+)(\d+)$", version)
    if not match:
        raise ValueError(f"Invalid genebuild version format: {version}")

    prefix = match.group(1)
    number = int(match.group(2))

    # Increment and preserve zero-padding
    new_number = number + 1
    # Preserve the original zero-padding width
    width = len(match.group(2))

    return f"{prefix}{new_number:0{width}d}"


def fetch_registry_ids(
    connection: pymysql.connections.Connection,
    assembly: str,
    genebuilder: Optional[str] = None,
) -> tuple[int, Optional[int]]:
    """
    Fetch both assembly_id and genebuild_status_id for a given assembly.

    Args:
        connection: MySQL connection object
        assembly (str): Assembly Accession (GCA format)
        genebuilder (str, optional): Genebuilder name to filter by. If None and multiple
                                     active records exist, returns the first one found.
    Returns:
        tuple: (assembly_id, genebuild_status_id) where genebuild_status_id may be None
    Raises:
        ValueError: If assembly is not found in registry
    """
    assembly_id = fetch_assembly_id(connection, assembly)
    if not assembly_id:
        raise ValueError(f"Assembly not found in registry: {assembly}")

    record = fetch_current_genebuild_record(connection, assembly, genebuilder)
    genebuild_status_id = record["genebuild_status_id"] if record else None

    return assembly_id, genebuild_status_id

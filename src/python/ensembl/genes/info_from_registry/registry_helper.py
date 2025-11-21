"""
Registry Helper Module

This module provides utility functions for common registry database queries
that are shared across multiple scripts.
"""

from typing import Optional
import pymysql


def fetch_assembly_id(connection: pymysql.connections.Connection, assembly: str) -> Optional[int]:
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

    with connection.cursor() as cursor:
        cursor.execute(query, (assembly,))
        result = cursor.fetchone()

    return result["assembly_id"] if result else None


def fetch_current_genebuild_record(
    connection: pymysql.connections.Connection, assembly: str, genebuilder: Optional[str] = None
) -> Optional[dict[str, any]]:
    """
    Fetch the full current genebuild record for a given assembly.

    Args:
        connection: MySQL connection object
        assembly (str): Assembly Accession (GCA format)
        genebuilder (str, optional): Genebuilder name to filter by. If None, fetches any active record.
    Returns:
        dict: Record with genebuild_status_id, gb_status, genebuilder, annotation_method, genebuild_version, and other fields if found, else None
    """
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

    with connection.cursor() as cursor:
        cursor.execute(query, params)
        result = cursor.fetchone()

    return result


def fetch_genebuild_status_id(connection: pymysql.connections.Connection, assembly: str) -> Optional[int]:
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


def fetch_registry_ids(
    connection: pymysql.connections.Connection, assembly: str, genebuilder: Optional[str] = None
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

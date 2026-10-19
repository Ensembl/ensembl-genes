"""Reusable queries and utilities for the Ensembl registry database."""

import re
from typing import Any, TypeAlias, cast

import pymysql

from ensembl.genes.mysql_helper import MySQLConnection

RegistryConnection: TypeAlias = MySQLConnection


def fetch_assembly_id(connection: MySQLConnection, assembly: str) -> int | None:
    """Return the registry assembly ID for an accession."""

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
    connection: RegistryConnection,
    assembly: str,
    genebuilder: str | None = None,
) -> dict[str, Any] | None:
    """Return the active genebuild record for an assembly."""

    if genebuilder:
        query = """
        SELECT genebuild_status_id, gb_status, genebuilder, annotation_method, genebuild_version
        FROM genebuild_status
        WHERE gca_accession = %s AND genebuilder = %s AND last_attempt = 1
        """
        params: tuple[Any, ...] = (assembly, genebuilder)
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

    return cast(dict[str, Any], result) if result is not None else None


def fetch_genebuild_status_id(
    connection: RegistryConnection, assembly: str
) -> int | None:
    """Return the active genebuild status ID for an assembly."""

    record = fetch_current_genebuild_record(connection, assembly)
    return record["genebuild_status_id"] if record else None


def fetch_highest_genebuild_version(
    connection: RegistryConnection, assembly: str
) -> str | None:
    """Return the highest genebuild version recorded for an assembly."""

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
    """Increment a version such as ``ENS01`` to ``ENS02``."""

    match = re.match(r"^([A-Z]+)(\d+)$", version)
    if not match:
        raise ValueError(f"Invalid genebuild version format: {version}")

    prefix = match.group(1)
    number = int(match.group(2))
    width = len(match.group(2))
    return f"{prefix}{number + 1:0{width}d}"


def fetch_registry_ids(
    connection: RegistryConnection,
    assembly: str,
    genebuilder: str | None = None,
) -> tuple[int, int | None]:
    """Return the assembly and active genebuild status IDs."""

    assembly_id = fetch_assembly_id(connection, assembly)
    if not assembly_id:
        raise ValueError(f"Assembly not found in registry: {assembly}")

    record = fetch_current_genebuild_record(connection, assembly, genebuilder)
    genebuild_status_id = record["genebuild_status_id"] if record else None
    return assembly_id, genebuild_status_id

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
Write metrics from core database to registry database.

Reads metrics from core DB meta table and writes them to registry tables:
 - assembly.* keys -> assembly_metrics
 - genebuild.* keys -> annotation_metrics

Uses DELETE then INSERT pattern to avoid stale/duplicated metrics.
"""

# pylint:disable=f-string-without-interpolation, broad-exception-raised
import argparse
import sys
from typing import Optional
import pymysql
from ensembl.genes.info_from_registry.mysql_helper import mysql_get_connection
from ensembl.genes.info_from_registry.registry_helper import fetch_registry_ids


def fetch_core_metrics(
    core_connection: pymysql.connections.Connection, species_id: int
) -> list[dict[str, str]]:
    """
    Fetch metrics from core database meta table.

    Args:
        core_connection: MySQL connection object for core DB
        species_id (int): Species ID in core meta table
    Returns:
        list: List of dictionaries with meta_key and meta_value
    """
    query = """
    SELECT meta_key, meta_value
    FROM meta
    WHERE species_id=%s AND (
    meta_key LIKE 'assembly.busco%%' OR
    meta_key LIKE 'assembly.stats.%%' OR
    meta_key LIKE 'genebuild.busco%%' OR
    meta_key LIKE 'genebuild.stats.%%' OR
        meta_key = 'genebuild.last_geneset_update'
    )
    """
    with core_connection.cursor() as cursor:
        cursor.execute(query, (species_id,))
        return cursor.fetchall()  # type: ignore[return-value]


def partition_metrics(
    rows: list[dict[str, str]],
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """
    Partition metrics into assembly and genebuild categories.

    Args:
        rows (list): List of dicts with meta_key and meta_value
    Returns:
        tuple: (assembly_rows, genebuild_rows) as lists of (key, value) tuples
    """
    assembly_rows = []
    genebuild_rows = []

    for row in rows:
        key = row["meta_key"]
        value = row["meta_value"]

        if key.startswith("assembly.busco") or key.startswith("assembly.stats."):
            assembly_rows.append((key, value))
        elif (
            key.startswith("genebuild.busco")
            or key.startswith("genebuild.stats.")
            or key == "genebuild.last_geneset_update"
        ):
            genebuild_rows.append((key, value))

    return assembly_rows, genebuild_rows


def write_assembly_metrics(
    registry_connection: pymysql.connections.Connection,
    assembly_id: int,
    rows: list[tuple[str, str]],
    dev: bool,
) -> None:
    """
    Write assembly metrics to registry using DELETE then INSERT pattern.

    Args:
        registry_connection: MySQL connection object for registry DB
        assembly_id (int): Assembly ID in registry
        rows (list): List of (metric_name, metric_value) tuples
        dev (bool): If True, only print SQL statements without executing them

    Returns: None
    """
    if not rows:
        print("No assembly metrics to write")
        return

    # Extract metric names for deletion
    metric_names = [name for name, value in rows]

    # Single DELETE for all metrics at once
    placeholders = ",".join(["%s"] * len(metric_names))
    delete_query = f"""
    DELETE FROM assembly_metrics
    WHERE assembly_id=%s AND metrics_name IN ({placeholders})
    """

    insert_query = """
    INSERT INTO assembly_metrics (assembly_id, metrics_name, metrics_value)
    VALUES (%s, %s, %s)
    """

    if dev:
        print("Would execute:")
        print(delete_query, (assembly_id, *metric_names))
        for name, value in rows:
            print(insert_query, (assembly_id, name, value))
    else:
        with registry_connection.cursor() as cursor:
            print(delete_query, (assembly_id, *metric_names))

            # Single DELETE for all metrics
            cursor.execute(delete_query, (assembly_id, *metric_names))
            deleted = cursor.rowcount
            print(
                f"Deleted {deleted} existing assembly metrics for assembly_id {assembly_id}"
            )

            # Bulk INSERT
            cursor.executemany(
                insert_query, [(assembly_id, name, value) for name, value in rows]
            )
        print(f"Wrote {len(rows)} assembly metrics for assembly_id {assembly_id}")


def write_genebuild_metrics(
    registry_connection: pymysql.connections.Connection,
    genebuild_status_id: Optional[int],
    assembly_id: int,
    rows: list[tuple[str, str]],
    dev: bool,
) -> None:
    """
    Write genebuild metrics to registry using DELETE then INSERT pattern.

    Args:
        registry_connection: MySQL connection object for registry DB
        genebuild_status_id (int): Genebuild status ID in registry
        assembly_id (int): Assembly ID in registry
        rows (list): List of (metric_name, metric_value) tuples
        dev (bool): If True, only print SQL statements without executing them

    Returns: None
    """
    if not rows:
        print("No genebuild metrics to write")
        return

    if not genebuild_status_id:
        print("No genebuild_status_id available; skipping genebuild metrics")
        return

    metric_names = [name for name, value in rows]

    placeholders = ",".join(["%s"] * len(metric_names))
    delete_query = f"""
    DELETE FROM annotation_metrics
    WHERE genebuild_status_id=%s AND metrics_name IN ({placeholders})
    """

    insert_query = """
    INSERT INTO annotation_metrics (genebuild_status_id, assembly_id, metrics_name, metrics_value)
    VALUES (%s, %s, %s, %s)
    """

    if dev:
        print("Would execute:")
        print(delete_query, (genebuild_status_id, *metric_names))
        for name, value in rows:
            print(insert_query, (genebuild_status_id, assembly_id, name, value))
    else:
        with registry_connection.cursor() as cursor:
            cursor.execute(delete_query, (genebuild_status_id, *metric_names))
            deleted = cursor.rowcount
            print(
                f"Deleted {deleted} existing genebuild metrics for \
                    genebuild_status_id {genebuild_status_id}"
            )

            cursor.executemany(
                insert_query,
                [
                    (genebuild_status_id, str(assembly_id), name, value)
                    for name, value in rows
                ],
            )
        print(
            f"Wrote {len(rows)} genebuild metrics for genebuild_status_id {genebuild_status_id}"
        )


def main(  # pylint:disable=too-many-arguments, too-many-locals
    registry_host: str,
    registry_port: int,
    registry_user: str,
    registry_password: str,
    registry_db: str,
    core_host: str,
    core_port: int,
    core_user: str,
    core_password: Optional[str],
    core_db: str,
    assembly: str,
    species_id: int,
    genebuilder: str,
    dev: bool,
) -> None:
    """
    Main function to write metrics from core DB to registry DB.

    Args:
        registry_host, registry_port, registry_user, registry_password,\
            registry_db: Registry DB connection
        core_host, core_port, core_user, core_password, core_db: Core DB connection
        assembly (str): Assembly Accession (GCA format)
        species_id (int): Species ID in core meta table
        genebuilder (str): Genebuilder name to identify which genebuild record to update
        dev (bool): If True, only print SQL statements without executing
    """
    registry_connection = None
    core_connection = None

    try:
        registry_connection = mysql_get_connection(
            database=registry_db,
            host=registry_host,
            port=registry_port,
            user=registry_user,
            password=registry_password,
        )

        if not registry_connection:
            raise Exception("Failed to connect to the registry database")

        core_connection = mysql_get_connection(
            database=core_db,
            host=core_host,
            port=core_port,
            user=core_user,
            password=core_password,
        )

        if not core_connection:
            raise Exception("Failed to connect to the core database")

        # Start transaction on registry
        registry_connection.begin()

        print(
            f"Fetching registry IDs for assembly {assembly} and genebuilder {genebuilder}..."
        )
        assembly_id, genebuild_status_id = fetch_registry_ids(
            registry_connection, assembly, genebuilder
        )
        print(
            f"Found assembly_id: {assembly_id}, genebuild_status_id: {genebuild_status_id}"
        )

        print(f"Fetching metrics from core database...")
        meta_rows = fetch_core_metrics(core_connection, species_id)
        if not meta_rows:
            print("No metrics found in core database")
            sys.exit(0)

        # Partition metrics
        assembly_rows, genebuild_rows = partition_metrics(meta_rows)
        print(
            f"Found {len(assembly_rows)} assembly metrics and \
                {len(genebuild_rows)} genebuild metrics"
        )

        # Write metrics to registry
        write_assembly_metrics(registry_connection, assembly_id, assembly_rows, dev)
        write_genebuild_metrics(
            registry_connection, genebuild_status_id, assembly_id, genebuild_rows, dev
        )

        if not dev:
            # Commit transaction
            registry_connection.commit()
            print("Metrics written successfully")
        else:
            print("DEV MODE: No changes were made to the registry database")

    except Exception as e:  # pylint:disable=broad-exception-caught
        if registry_connection:
            registry_connection.rollback()
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)

    finally:
        if registry_connection:
            registry_connection.close()
        if core_connection:
            core_connection.close()

    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Write metrics from core database to registry database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Fetches metrics from core meta table and writes to registry tables.

            Examples:
            %(prog)s --registry_host localhost --registry_user myuser --registry_password mypass \\
                     --registry_db registry --core_host corehost --core_user myuser \\
                     --core_password mypass --core_db my_core_db --assembly GCA_123456789.1 \\
                     --genebuilder john_doe

            %(prog)s --registry_host localhost --registry_user myuser --registry_password mypass \\
                     --registry_db registry --core_host corehost --core_user myuser \\
                     --core_password mypass --core_db my_core_db --assembly GCA_123456789.1 \\
                     --genebuilder john_doe --dev
        """,
    )

    # Registry connection
    parser.add_argument("--registry_host", required=True, help="Registry MySQL host")
    parser.add_argument(
        "--registry_port", type=int, default=3306, help="Registry MySQL port"
    )
    parser.add_argument("--registry_user", required=True, help="Registry MySQL user")
    parser.add_argument(
        "--registry_password", required=True, help="Registry MySQL password"
    )
    parser.add_argument(
        "--registry_db", required=True, help="Registry MySQL database name"
    )

    # Core connection
    parser.add_argument("--core_host", required=True, help="Core MySQL host")
    parser.add_argument("--core_port", type=int, default=3306, help="Core MySQL port")
    parser.add_argument("--core_user", required=True, help="Core MySQL user")
    parser.add_argument("--core_password", default=None, help="Core MySQL password")
    parser.add_argument("--core_db", required=True, help="Core MySQL database name")

    # Context
    parser.add_argument(
        "--assembly",
        required=True,
        help="Assembly Accession (GCA format, e.g., GCA_123456789.1)",
    )
    parser.add_argument(
        "--genebuilder",
        required=True,
        help="Genebuilder name (person who owns the genebuild record)",
    )
    parser.add_argument(
        "--species_id",
        type=int,
        default=1,
        help="Species ID in core meta table (default: 1)",
    )

    # Options
    parser.add_argument(
        "--dev",
        action="store_true",
        help="Development mode: print SQL statements without executing them",
    )

    args = parser.parse_args()

    main(
        registry_host=args.registry_host,
        registry_port=args.registry_port,
        registry_user=args.registry_user,
        registry_password=args.registry_password,
        registry_db=args.registry_db,
        core_host=args.core_host,
        core_port=args.core_port,
        core_user=args.core_user,
        core_password=args.core_password,
        core_db=args.core_db,
        assembly=args.assembly,
        genebuilder=args.genebuilder,
        species_id=args.species_id,
        dev=args.dev,
    )

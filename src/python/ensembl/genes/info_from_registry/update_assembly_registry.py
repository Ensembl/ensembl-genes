"""
Utility to update the status of a genebuild in the registry db
"""

# pylint: disable=f-string-without-interpolation, broad-exception-raised, broad-exception-caught
import argparse
from datetime import datetime
import sys
from typing import Optional
import pymysql
from ensembl.genes.info_from_registry.mysql_helper import mysql_get_connection
from ensembl.genes.info_from_registry.registry_helper import (
    fetch_assembly_id,
    fetch_current_genebuild_record,
)


def ensure_genebuilder_exists(
    connection: pymysql.connections.Connection, genebuilder: str
) -> None:
    """
    Ensure genebuilder exists in the genebuilder table.

    Args:
        connection: MySQL connection object
        genebuilder (str): Genebuilder name
    """
    check_query = "SELECT genebuilder FROM genebuilder WHERE genebuilder = %s"

    with connection.cursor() as cursor:
        cursor.execute(check_query, (genebuilder,))
        if not cursor.fetchone():
            raise ValueError(
                f"Genebuilder '{genebuilder}' does not exist in the\
                    genebuilder table. Please add it first."
            )


# fetch_current_record is now fetch_current_genebuild_record imported from registry_helper
def insert_new_record(  # pylint:disable=too-many-arguments, too-many-positional-arguments
    connection: pymysql.connections.Connection,
    assembly_id: int,
    assembly: str,
    genebuilder: str,
    status: str,
    annotation_source: str,
    annotation_method: str,
    current_date: str,
    release_type: str,
    genebuild_version: str,
    dev: bool,
) -> None:
    """
    Insert a new genebuild status record.

    Args:
        connection: MySQL connection object
        assembly_id (int): Assembly ID
        assembly (str): Assembly accession
        genebuilder (str): Genebuilder name
        status (str): Status to set
        annotation_source (str): Annotation source
        annotation_method (str): Annotation method
        current_date (str): Current date
        release_type (str): Release type (default: "not_available")
        dev (bool): If True, only print SQL without executing
    """
    # date completed is being used to track the last time this record was updated
    # so we set it to the current date
    date_status_update = current_date

    # Add not avaialable for release type
    query = """
    INSERT INTO genebuild_status (
        assembly_id, gca_accession, gb_status, last_attempt, 
        genebuilder, annotation_source, annotation_method,
        date_started, date_status_update, release_type,
        genebuild_version
    )
    VALUES (%s, %s, %s, 1, %s, %s, %s, %s, %s, %s, %s)
    """

    params = (
        assembly_id,
        assembly,
        status,
        genebuilder,
        annotation_source,
        annotation_method,
        current_date,
        date_status_update,
        release_type,
        genebuild_version,
    )

    if dev:
        print("Would execute:")
        print(query % params)
    else:
        with connection.cursor() as cursor:
            cursor.execute(query, params)
        print(f"Inserted new record for GCA {assembly} with status '{status}'")


def update_existing_record(  # pylint:disable=too-many-arguments, too-many-positional-arguments
    connection: pymysql.connections.Connection,
    record_id: int,
    status: str,
    current_date: str,
    dev: bool,
    annotation_method: Optional[str] = None,
    annotation_source: Optional[str] = None,
    genebuild_version: Optional[str] = None,
) -> None:
    """
    Update an existing genebuild status record.

    Args:
        connection: MySQL connection object
        record_id (int): genebuild_status_id to update
        status (str): Status to set
        current_date (str): Current date
        dev (bool): If True, only print SQL without executing
        annotation_method (str, optional): New annotation method
        annotation_source (str, optional): New annotation source
        genebuild_version (str, optional): New genebuild version
    """
    query_parts = ["gb_status = %s", "date_status_update = %s"]
    params = [status, current_date]
    if annotation_method:
        query_parts.append("annotation_method = %s")
        params.append(annotation_method)
    if annotation_source:
        query_parts.append("annotation_source = %s")
        params.append(annotation_source)
    if genebuild_version:
        query_parts.append("genebuild_version = %s")
        params.append(genebuild_version)
    query = f"""
UPDATE genebuild_status
SET {', '.join(query_parts)}
WHERE genebuild_status_id = %s
"""
    params.append(str(record_id))

    if dev:
        print("Would execute:")
        print(query, params)
    else:
        with connection.cursor() as cursor:
            cursor.execute(query, params)
        print(f"Updated record {record_id} to status '{status}'")


def set_old_record_historical(
    connection: pymysql.connections.Connection, record_id: int, dev: bool
) -> None:
    """
    Set an existing record to historical (last_attempt = 0).

    Args:
        connection: MySQL connection object
        record_id (int): genebuild_status_id to update
        dev (bool): If True, only print SQL without executing
    """
    query = (
        "UPDATE genebuild_status SET last_attempt = 0 WHERE genebuild_status_id = %s"
    )

    if dev:
        print("Would execute:")
        print(query % (record_id,))
    else:
        with connection.cursor() as cursor:
            cursor.execute(query, (record_id,))
        print(f"Set record {record_id} to historical")


def main(  # pylint:disable=too-many-arguments, too-many-statements, too-many-branches, too-many-locals, too-many-positional-arguments
    host: str,
    port: int,
    user: str,
    password: str,
    database: str,
    assembly: str,
    status: str,
    genebuilder: str,
    annotation_source: Optional[str],
    annotation_method: Optional[str],
    release_type: str,
    genebuild_version: Optional[str],
    dev: bool = False,
) -> None:
    """
    Main function to update the genebuild status in the registry database.

    Args:
        host (str): MySQL host
        port (int): MySQL port
        user (str): MySQL user
        password (str): MySQL password
        database (str): MySQL database name
        assembly (str): Assembly Accession (GCA format)
        status (str): Status to update to
        genebuilder (str): Genebuilder name
        annotation_source (str): Annotation source
        annotation_method (str): Annotation method
        release_type (str): Release type (default: "not_available")
        genebuild_version (str): Genebuild version
        dev (bool): If True, only print SQL statements without executing
    """
    connection = None

    try:
        # Connect to database
        connection = mysql_get_connection(
            database=database, host=host, port=port, user=user, password=password
        )

        if not connection:
            raise Exception("Failed to connect to the database")

        # Start transaction
        connection.begin()

        ensure_genebuilder_exists(connection, genebuilder)

        assembly_id = fetch_assembly_id(connection, assembly)
        if not assembly_id:
            raise Exception(f"Assembly not found for GCA: {assembly}")

        current_date = datetime.now().strftime("%Y-%m-%d")
        existing_record = fetch_current_genebuild_record(
            connection, assembly, genebuilder
        )

        if not existing_record:
            # No existing record for this genebuilder - INSERT new
            print(f"No existing record found for {assembly} by {genebuilder}")
            print(f"Creating new record with status: {status}")

            # Use defaults for new records if not specified
            method_to_insert = annotation_method if annotation_method else "pending"
            source_to_insert = annotation_source if annotation_source else "ensembl"
            version_to_insert = genebuild_version if genebuild_version else "ENS01"

            insert_new_record(
                connection,
                assembly_id,
                assembly,
                genebuilder,
                status,
                source_to_insert,
                method_to_insert,
                current_date,
                release_type,
                version_to_insert,
                dev,
            )

        else:
            # Existing record found for this genebuilder
            current_status = existing_record["gb_status"]
            current_method = existing_record.get("annotation_method")
            current_source = existing_record.get("annotation_source")
            current_version = existing_record.get("genebuild_version")
            record_id = existing_record["genebuild_status_id"]

            print(f"Found existing record for {assembly} by {genebuilder}")
            print(f"Current status: '{current_status}', requested status: '{status}'")

            # Status categories
            terminal_statuses = ["live", "pre_released", "handed_over", "archive"]
            active_statuses = ["insufficient_data", "in_progress", "check_busco"]
            completed_status = "completed"  # terminal-like status

            # Determine if method/source/version should be updated
            # (preserve existing if not specified)
            method_to_update = (
                annotation_method if annotation_method is not None else None
            )
            source_to_update = (
                annotation_source if annotation_source is not None else None
            )
            version_to_update = (
                genebuild_version
                if (
                    genebuild_version is not None
                    and genebuild_version != current_version
                )
                else None
            )

            # Same status - check for method/source/version changes
            if current_status == status:
                # Check if annotation_method, annotation_source, or
                # genebuild_version is provided and different
                method_changed = (
                    annotation_method and annotation_method != current_method
                )
                source_changed = (
                    annotation_source and annotation_source != current_source
                )
                version_changed = (
                    genebuild_version is not None
                    and genebuild_version != current_version
                )

                if method_changed or source_changed or version_changed:
                    changes = []
                    if method_changed:
                        changes.append(
                            f"method from '{current_method}' to '{annotation_method}'"
                        )
                    if source_changed:
                        changes.append(
                            f"source from '{current_source}' to '{annotation_source}'"
                        )
                    if version_changed:
                        changes.append(
                            f"version from '{current_version}' to '{genebuild_version}'"
                        )
                    print(
                        f"Status is already '{status}', but updating {' and '.join(changes)}"
                    )
                    update_existing_record(
                        connection,
                        record_id,
                        status,
                        current_date,
                        dev,
                        annotation_method if method_changed else None,
                        annotation_source if source_changed else None,
                        genebuild_version if version_changed else None,
                    )
                else:
                    print(f"Status is already '{status}'. No changes needed.")
                    sys.exit(0)

            # Block backwards transitions from completed
            elif current_status == completed_status and status in active_statuses:
                print(f"ERROR: Cannot move backwards from 'completed' to '{status}'")
                print(
                    f"Completed is a terminal-like status. To restart work, \
                        use a new genebuild_version."
                )
                sys.exit(1)

            # Block backwards transitions from terminal statuses
            elif current_status in terminal_statuses and status in (
                active_statuses + [completed_status]
            ):
                # Check if version changed - if so, allow new attempt
                # If genebuild_version not provided, use current version (no version change)
                effective_version = (
                    genebuild_version if genebuild_version else current_version
                )

                if effective_version != current_version:
                    print(
                        f"Moving from terminal status '{current_status}' to \
                            '{status}' with new version {effective_version}"
                    )
                    print(f"Creating new attempt.")

                    method_to_insert = (
                        annotation_method if annotation_method else "pending"
                    )
                    source_to_insert = (
                        annotation_source if annotation_source else current_source
                    )

                    set_old_record_historical(connection, record_id, dev)
                    insert_new_record(
                        connection,
                        assembly_id,
                        assembly,
                        genebuilder,
                        status,
                        source_to_insert,
                        method_to_insert,
                        current_date,
                        release_type,
                        effective_version,
                        dev,
                    )
                else:
                    print(
                        f"ERROR: Cannot move from terminal status '{current_status}' \
                            to '{status}' with same genebuild_version '{current_version}'"
                    )
                    print(
                        f"To restart work, you must provide a new \
                            --genebuild_version (e.g., ENS02, ENS03, etc.)"
                    )
                    sys.exit(1)

            # Terminal to terminal transition - UPDATE same record
            elif current_status in terminal_statuses and status in terminal_statuses:
                print(
                    f"Moving from terminal status '{current_status}' to terminal status '{status}'"
                )
                print(f"Updating existing record.")
                update_existing_record(
                    connection,
                    record_id,
                    status,
                    current_date,
                    dev,
                    method_to_update,
                    source_to_update,
                    version_to_update,
                )

            # Completed to terminal transition - UPDATE same record
            elif current_status == completed_status and status in terminal_statuses:
                print(f"Moving from 'completed' to terminal status '{status}'")
                print(f"Updating existing record.")
                update_existing_record(
                    connection,
                    record_id,
                    status,
                    current_date,
                    dev,
                    method_to_update,
                    source_to_update,
                    version_to_update,
                )

            # Active to active or active to completed - UPDATE same record
            elif current_status in active_statuses and (
                status in active_statuses or status == completed_status
            ):
                print(
                    f"Current status '{current_status}' is active. Updating to '{status}'."
                )
                update_existing_record(
                    connection,
                    record_id,
                    status,
                    current_date,
                    dev,
                    method_to_update,
                    source_to_update,
                    version_to_update,
                )

            # Active to terminal - UPDATE same record
            elif current_status in active_statuses and status in terminal_statuses:
                print(
                    f"Moving from active status '{current_status}' to terminal status '{status}'"
                )
                print(f"Updating existing record.")
                update_existing_record(
                    connection,
                    record_id,
                    status,
                    current_date,
                    dev,
                    method_to_update,
                    source_to_update,
                    version_to_update,
                )

            else:
                raise ValueError(
                    f"Unexpected status transition: '{current_status}' to '{status}'"
                )

        if not dev:
            # Commit transaction
            connection.commit()
            print("Transaction completed successfully")
        else:
            print("DEV MODE: No changes were made to the database")

    except Exception as e:
        if connection:
            connection.rollback()
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)

    finally:
        if connection:
            connection.close()

    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Update the status of a genebuild in the registry database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Status values should match the gb_status enum:
            in_progress, insufficient_data, check_busco, completed, 
            pre_released, handed_over, archive

            Examples:
            %(prog)s --host localhost --user myuser --password mypass --database registry 
                    --assembly GCA_123456789.1 --status completed --genebuilder john_doe

            %(prog)s --host localhost --user myuser --password mypass --database registry 
                    --assembly GCA_123456789.1 --status in_progress --genebuilder jane_smith --dev
        """,
    )

    parser.add_argument("--host", required=True, help="MySQL host")
    parser.add_argument("--port", type=int, default=3306, help="MySQL port")
    parser.add_argument("--user", required=True, help="MySQL user")
    parser.add_argument("--password", required=True, help="MySQL password")
    parser.add_argument("--database", required=True, help="MySQL database name")
    parser.add_argument(
        "--assembly",
        required=True,
        help="Assembly Accession (GCA format, e.g., GCA_123456789.1)",
    )
    parser.add_argument(
        "--status",
        required=True,
        choices=[
            "in_progress",
            "insufficient_data",
            "check_busco",
            "completed",
            "pre_released",
            "handed_over",
            "live",
            "archive",
        ],
        help="Status to update to",
    )
    parser.add_argument(
        "--genebuilder",
        required=True,
        help="Genebuilder name (person doing the annotation)",
    )
    parser.add_argument(
        "--dev",
        action="store_true",
        help="Development mode: print SQL statements without executing them",
    )

    parser.add_argument(
        "--annotation_source",
        default=None,
        choices=[
            "ensembl",
            "helixer",
            "external",
            "import_refseq",
            "import_community",
            "import_wormbase",
            "import_flybase",
            "import_genbank",
            "import_noninsdc",
        ],
        help="Annotation source (default: preserve existing value, \
            or 'ensembl' for new records)",
    )

    parser.add_argument(
        "--annotation_method",
        default=None,
        choices=[
            "pending",
            "full_genebuild",
            "anno",
            "braker",
            "helixer",
            "projection_build",
            "mixed_strategy_build",
            "import",
            "external_annotation_import",
        ],
        help="Annotation method (default: preserve existing value)",
    )

    parser.add_argument(
        "--release_type",
        default="not_available",
        choices=["main", "beta", "not_available"],
        help="Release type (default: not_available)",
    )

    parser.add_argument(
        "--genebuild_version",
        default=None,
        help="Genebuild version (default: preserve existing value, \
            or 'ENS01' for new records)",
    )

    args = parser.parse_args()

    main(
        host=args.host,
        port=args.port,
        user=args.user,
        password=args.password,
        database=args.database,
        assembly=args.assembly,
        status=args.status,
        genebuilder=args.genebuilder,
        annotation_source=args.annotation_source,
        annotation_method=args.annotation_method,
        release_type=args.release_type,
        genebuild_version=args.genebuild_version,
        dev=args.dev,
    )

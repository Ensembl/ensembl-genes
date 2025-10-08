"""
Utility to update the status of a genebuild in the registry db 
"""

from mysql_helper import mysql_get_connection
import argparse
from datetime import datetime


def ensure_genebuilder_exists(connection, genebuilder):
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
                f"Genebuilder '{genebuilder}' does not exist in the genebuilder table. Please add it first."
            )


def fetch_current_record(connection, assembly):
    """
    Fetch the current genebuild record for a given assembly.

    Args:
        connection: MySQL connection object
        assembly (str): Assembly Accession (GCA)
    Returns:
        dict: Record with genebuild_status_id and gb_status if found, else None
    """
    query = """
    SELECT genebuild_status_id, gb_status 
    FROM genebuild_status 
    WHERE gca_accession = %s AND last_attempt = 1
    """

    with connection.cursor() as cursor:
        cursor.execute(query, (assembly,))
        result = cursor.fetchone()

    return result


def fetch_assembly_id(connection, assembly):
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


def insert_new_record(
    connection,
    assembly_id,
    assembly,
    genebuilder,
    status,
    annotation_source,
    annotation_method,
    current_date,
    release_type,
    genebuild_version,
    dev,
):
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


def update_existing_record(connection, record_id, status, current_date, dev, annotation_method=None):
    query_parts = ["gb_status = %s", "date_status_update = %s"]
    params = [status, current_date]
    if annotation_method:
        query_parts.append("annotation_method = %s")
        params.append(annotation_method)
    query = f"""
UPDATE genebuild_status
SET {', '.join(query_parts)}
WHERE genebuild_status_id = %s
"""
    params.append(record_id)

    if dev:
        print("Would execute:")
        print(query, params)
    else:
        with connection.cursor() as cursor:
            cursor.execute(query, params)
        print(f"Updated record {record_id} to status '{status}'")


def set_old_record_historical(connection, record_id, dev):
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


def main(
    host,
    port,
    user,
    password,
    database,
    assembly,
    status,
    genebuilder,
    annotation_source,
    annotation_method,
    release_type,
    genebuild_version,
    dev=False,
):
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

        # Ensure genebuilder exists
        ensure_genebuilder_exists(connection, genebuilder)

        # Get assembly_id
        assembly_id = fetch_assembly_id(connection, assembly)
        if not assembly_id:
            raise Exception(f"Assembly not found for GCA: {assembly}")

        # Get current date
        current_date = datetime.now().strftime("%Y-%m-%d")

        # Check for existing current record
        existing_record = fetch_current_record(connection, assembly)

        if not existing_record:
            # No existing record - INSERT new
            print(f"No existing record found for {assembly}")
            print(f"Creating new record with status: {status}")

            insert_new_record(
                connection,
                assembly_id,
                assembly,
                genebuilder,
                status,
                annotation_source,
                annotation_method,
                current_date,
                release_type,
                genebuild_version,
                dev,
            )

        else:
            # Existing record found
            current_status = existing_record["gb_status"]
            record_id = existing_record["genebuild_status_id"]

            print(f"Found existing record for {assembly} by {genebuilder}")
            print(f"Current status: '{current_status}', requested status: '{status}'")

            # Map status names to match schema enum values
            terminal_statuses = ["completed", "pre_released", "handed_over", "archive"]
            active_statuses = ["insufficient_data", "in_progress", "check_busco"]

            if current_status == status:
                if current_status in active_statuses:
                    print("Status unchanged and active; updating method/timestamp if provided.")
                    update_existing_record(
                        connection, record_id, status, current_date, dev, annotation_method
                    )
                else:
                    print(f"Status is already '{status}'. No changes needed.")
                return

            if current_status in terminal_statuses:
                # Terminal status - create new attempt
                print(
                    f"Current status '{current_status}' is terminal. Creating new attempt."
                )

                set_old_record_historical(connection, record_id, dev)
                insert_new_record(
                    connection,
                    assembly_id,
                    assembly,
                    genebuilder,
                    status,
                    annotation_source,
                    annotation_method,
                    current_date,
                    release_type,
                    genebuild_version,
                    dev,
                )

            elif current_status in active_statuses:
                # Active status - update existing record
                print(
                    f"Current status '{current_status}' is active. Updating existing record."
                )

                update_existing_record(connection, record_id, status, current_date, dev, annotation_method)

            else:
                print(
                    f"Unexpected current status: '{current_status}'. No action taken."
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
        print(f"Error occurred: {str(e)}")
        raise

    finally:
        if connection:
            connection.close()


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
        default="ensembl",
        choices=[
            "ensembl",
            "external",
            "import_refseq",
            "import_community",
            "import_wormbase",
            "import_flybase",
            "import_genbank",
            "import_noninsdc",
        ],
        help="Annotation source (default: ensembl)",
    )

    parser.add_argument(
        "--annotation_method",
        default="pending",
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
        help="Annotation method (default: pending)",
    )

    parser.add_argument(
        "--release_type",
        default="not_available",
        choices=['main','beta','not_available'],
        help="Release type (default: not_available)",
    )

    parser.add_argument(
        "--genebuild_version",
        default="ENS01",
        help="Genebuild version (default: ENS01)",
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

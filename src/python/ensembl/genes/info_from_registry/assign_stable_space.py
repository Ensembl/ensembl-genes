import logging
from typing import Any
import pymysql  # type: ignore
from mysql_helper import mysql_fetch_data
from check_stable_space_old_registry import get_old_stable_space_info

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("pipeline_setup.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def insert_to_db(
    insert_query, conn: pymysql.connections.Connection, store_new_registry: bool
) -> bool:
    logger.info(f"Inserting data: {insert_query}")
    try:
        with conn.cursor() as cursor:
            cursor.execute(insert_query)
        return True
    except pymysql.err.IntegrityError as e:
        if e.args[0] == 1062:
            if store_new_registry:
                logger.info(
                    f"Duplicate entry found, but store_new_registry is True. Report as success."
                )
                return True
            else:
                return False
        raise


def stable_space_per_taxon(taxon_id: int, server_info: dict) -> int:
    """Get the next stable space ID for a given taxon.

    Args:
        taxon_id (int): The taxon identifier for which to find the stable space ID.
        server_info (dict): The server information for database connection.

    Returns:
        int: The next available stable space ID.
    """
    logger.info(f"Get the next stable space ID for taxon ID {taxon_id}")
    space_query = f"SELECT MAX(stable_space_id) as max_stable_id FROM stable_space_species_log WHERE lowest_taxon_id = {taxon_id};"
    output_query = mysql_fetch_data(
        space_query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database=server_info["registry"]["db_name"],
        password="",
    )

    stable_space_tmp = output_query[0].get("max_stable_id", None)

    ####################
    # remove when all pipelines have moved to the new registry:
    old_stable_space_info = get_old_stable_space_info(taxon_id, server_info)
    ####################

    if old_stable_space_info is None:
        if stable_space_tmp is None:
            stable_space_id = 1
            logger.info(
                f"No stable space found for taxon ID {taxon_id}. "
                f"Assigning new stable space ID {stable_space_id}."
            )
        else:
            stable_space_id = stable_space_tmp + 1
            logger.info(
                f"No old stable space for taxon ID {taxon_id}. "
                f"Using tmp value {stable_space_tmp}, assigning new ID {stable_space_id}."
            )

    else:  # old_stable_space_info is not None
        if stable_space_tmp is None:
            stable_space_id = old_stable_space_info + 1
            logger.info(
                f"Found old stable space {old_stable_space_info} for taxon ID {taxon_id}. "
                f"Assigning new ID {stable_space_id}."
            )
        elif stable_space_tmp == old_stable_space_info:
            stable_space_id = stable_space_tmp + 1
            logger.info(
                f"Stable space {stable_space_tmp} matches old record for taxon ID {taxon_id}. "
                f"Assigning new ID {stable_space_id}."
            )
        else:
            msg = (
                f"Conflicting stable space for taxon ID {taxon_id}: "
                f"{old_stable_space_info} != {stable_space_tmp}"
            )
            logger.error(msg)
            raise Exception(msg)

    return stable_space_id


def stable_space_range(stable_space_id: int, server_info: dict) -> bool:
    """Check if a stable space range exists for the given stable space ID.

    Args:
        stable_space_id (int): The stable space ID to check.
        server_info (dict): The server information for database connection.

    Returns:
        bool: True if the stable space range exists, False otherwise.
    """
    logger.info(f"Get stable space renage for {stable_space_id}")
    query = f"SELECT * FROM stable_space WHERE stable_space_id = '{stable_space_id}';"
    output_query = mysql_fetch_data(
        query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database=server_info["registry"]["db_name"],
        password="",
    )

    if output_query:
        logger.info(
            f"Stable space ID {stable_space_id} already exists with range: {output_query[0]}."
        )
        return True

    else:
        logger.info(
            f"No existing stable space range found for ID {stable_space_id}. Creating a new range space."
        )
        previous_space_id = stable_space_id - 1
        query = f"""SELECT * FROM stable_space WHERE stable_space_id = '{previous_space_id}';"""
        output_query = mysql_fetch_data(
            query,
            host=server_info["registry"]["db_host"],
            user=server_info["registry"]["db_user"],
            port=server_info["registry"]["db_port"],
            database=server_info["registry"]["db_name"],
            password="",
        )
        logger.info(
            f"Output query for previous stable space ID {previous_space_id}: {output_query}"
        )

        if output_query:
            new_start = (
                output_query[0].get("stable_space_end") + 1
            )  # Increment the end of the previous stable space by 1
            new_end = new_start + 4999999
            logger.info(
                f"New stable space range for ID {stable_space_id} will be from {new_start} to {new_end}."
            )

            insert_query = f"INSERT INTO stable_space (stable_space_id, stable_space_start, stable_space_end) VALUES ({stable_space_id}, {new_start}, {new_end});"
            conn = pymysql.connect(
                host=server_info["registry"]["db_host"],
                user=server_info["registry"]["db_user_w"],
                port=server_info["registry"]["db_port"],
                password=server_info["registry"]["password"],
                database=server_info["registry"]["db_name"],
            )

            if insert_to_db(insert_query, conn, store_new_registry=True):
                logger.info(
                    f"Successfully created stable space range for ID {stable_space_id}."
                )
                conn.close()
                return new_start
            else:
                logger.error(
                    f"Failed to insert stable space range for ID {stable_space_id}."
                )
                return False
        else:
            logger.error(
                f"Failed to create stable space range for ID {stable_space_id}. Previous space ID {previous_space_id} not found."
            )
            return False


def assign_stable_id(
    taxon_id: int, gca_accession: str, assembly_id: int, server_info: dict
) -> tuple[bool, int, int] | tuple[bool, None, None]:
    """Assign a stable space ID for a given taxon and GCA accession.

    Args:
        taxon_id (int): The taxon identifier for which to assign a stable space ID.
        gca_accession (str): The GCA accession number for the assembly.
        assembly_id (int): The assembly identifier.
        server_info (dict): The server information for database connection.

    Returns:
        tuple[bool, int|None, int|None]: A tuple containing a success flag and the assigned stable space ID or None if failed.
    """

    stable_space_id = stable_space_per_taxon(taxon_id, server_info)
    # Check if stable space range exists
    logger.info(f"Check if stable space range exists {stable_space_id}")
    stable_space_start = stable_space_range(stable_space_id, server_info)

    if stable_space_start is not False:
        logger.info(
            f"Assigned stable space ID {stable_space_id} for {gca_accession} and taxon ID {taxon_id}."
        )

        insert_query = f"""INSERT INTO stable_space_species_log (stable_space_id, lowest_taxon_id, gca_accession, assembly_id) 
            VALUES ({stable_space_id}, {taxon_id}, '{gca_accession}', {assembly_id});
            """
        conn = pymysql.connect(
            host=server_info["registry"]["db_host"],
            user=server_info["registry"]["db_user_w"],
            port=server_info["registry"]["db_port"],
            password=server_info["registry"]["password"],
            database=server_info["registry"]["db_name"],
        )

        if insert_to_db(insert_query, conn, store_new_registry=True):
            logger.info(
                f"Successfully inserted stable space ID {stable_space_id} for GCA {gca_accession}."
            )
            conn.close()
            return True, stable_space_id, stable_space_start
        else:
            logger.error(
                f"Failed to insert stable space ID {stable_space_id} for GCA {gca_accession}."
            )
            return False, None, None
    else:
        logger.error(
            f"Failed to assign stable space ID {stable_space_id} for GCA {gca_accession}."
        )
        return False, None, None


def get_stable_space(
    taxon_id: int, gca_accession: str, assembly_id: int, server_info: dict
) -> int:
    """Get the stable space ID for a given taxon and GCA accession.

    Args:
        taxon_id (int): The taxon identifier for which to find the stable space ID.
        gca_accession (str): The GCA accession number for the assembly.
        assembly_id (int): The assembly identifier.
        server_info (dict): The server information for database connection.

    Returns:
        int: assigned stable space ID or existing one if already assigned.
    """

    # Check if GCA already has assigned a stable space
    logger.info(f"Check if {gca_accession} already has assigned a stable space")
    space_gca_query = f"SELECT stable_space_species_log.stable_space_id, stable_space.stable_space_start FROM stable_space_species_log JOIN stable_space ON stable_space.stable_space_id = stable_space_species_log.stable_space_id  WHERE stable_space_species_log.gca_accession = '{gca_accession}';"
    output_query = mysql_fetch_data(
        space_gca_query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database=server_info["registry"]["db_name"],
        password="",
    )

    if output_query:
        stable_space_id = output_query[0].get("stable_space_id")
        stable_space_start = output_query[0].get("stable_space_start")
        logger.info(
            f"Stable space {stable_space_id} already assigned for GCA {gca_accession}."
        )
        return stable_space_start

    else:
        logger.info(
            f"No stable space found for assembly {gca_accession}. Checking taxon ID {taxon_id} for existing stable space."
        )

        success = False
        while not success:
            logger.info(
                f"Assigning stable space for taxon ID {taxon_id} and {gca_accession}."
            )
            success, stable_space_id, stable_space_start = assign_stable_id(
                taxon_id, gca_accession, assembly_id, server_info
            )

        return int(stable_space_start)

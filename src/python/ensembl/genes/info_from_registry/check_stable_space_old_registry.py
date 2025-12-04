"""
check_stable_space_from_old_registry.py

This module gets legacy stable id from old registry. Module will be deprecated when
all pipelines have moved to the new registry. It's main functions is to facilitate
using the new registry to start anno runs.
"""
import logging
from mysql_helper import mysql_fetch_data
from typing import Dict, Any, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler("pipeline_setup.log"), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def get_old_stable_space_info(
    taxon_id: int, server_info: Dict[str, Dict[str, Any]]
) -> Optional[int]:
    """Get the stable space that has been assigned to old registry.

    Args:
        taxon_id (int): The taxon identifier for which to find the stable space ID.
        server_info (dict): The server information for database connection.

    Returns:
        int: Assigned stable space.
    """
    logger.info(f"Check if taxon ID {taxon_id} has stable space in old registry")
    space_query = (
        f"SELECT MAX(stable_id_space_id) as old_stable_id "
        f"FROM assembly "
        f"WHERE taxonomy = {taxon_id};"
    )

    output_query = mysql_fetch_data(
        space_query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database="gb_assembly_registry",
        password="",
    )

    old_stable_id = output_query[0].get("max_stable_id", None)

    if old_stable_id is None:
        logger.info(f"Stable space not found in old registry for taxon ID {taxon_id}.")

    logger.info(
        f"Max stable space in the old registry for taxon ID {taxon_id} is {old_stable_id}."
    )

    return old_stable_id

"""
check_if_annotated.py

This module defines a function to verify if a given genome assembly accession
is already annotated in the genebuild registry. It logs annotation status and raises
an error if any assembly has already been annotated.
"""

from mysql_helper import mysql_fetch_data
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline_setup.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def check_if_annotated(assembly_accession, server_info):
    """
     Check if a genome assembly is already annotated in the genebuild registry.

     Parameters:
     ----------
     assembly_accession : str
         The GCA accession of the assembly to check (e.g., "GCA_000001405.28").

     server_info : dict
         Dictionary containing server connection details, with a 'registry' key
         pointing to another dict with the following keys:
             - db_host : str
             - db_user : str
             - db_port : int
             - db_name : str

     Raises:
     ------
     RuntimeError
         If the given assembly accession is already present in the registry.

     Returns:
     -------
     None
     """

    registry_query = f"""
        SELECT 
            gca_accession, 
            gb_status, 
            genebuilder
        FROM genebuild 
        WHERE gca_accession =  %s
    """

    registry_rows = mysql_fetch_data(
        registry_query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database=server_info["registry"]["db_name"],
        password="",
        params=(assembly_accession,),
    )

    if registry_rows:
        for row in registry_rows:
            gca = row.get("gca_accession")
            gb_status = row.get("gb_status")
            genebuilder = row.get("genebuilder")
            logger.error(f"\n Annotation already exists for GCA: {gca}")
            logger.error(f"   Status     : {gb_status}")
            logger.error(f"   Genebuilder: {genebuilder}")
        raise RuntimeError("Terminating: One or more assemblies are already annotated.")

    logger.info("Start check complete: %s", assembly_accession)
    return None

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
    registry_query = f"""
        SELECT 
            assembly_accession, 
            progress_status, 
            genebuilder
        FROM genebuild_status 
        WHERE assembly_accession =  %s
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

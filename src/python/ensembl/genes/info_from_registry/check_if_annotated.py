from src.python.ensembl.genes.info_from_registry.mysql_helper import mysql_fetch_data


def check_if_annotated(assembly_accession, server_info):
    if isinstance(assembly_accession, str):
        assembly_accessions = [assembly_accession]
    else:
        assembly_accessions = assembly_accession

    placeholders = ",".join(["%s"] * len(assembly_accessions))

    registry_query = f"""
        SELECT 
            gca_accession, 
            gb_status, 
            genebuilder
        FROM genebuild 
        WHERE gca_accession IN ({placeholders})
    """

    registry_rows = mysql_fetch_data(
        registry_query,
        host=server_info["registry"]["db_host"],
        user=server_info["registry"]["db_user"],
        port=server_info["registry"]["db_port"],
        database="gb_assembly_metadata",
        params=assembly_accessions,
    )

    if registry_rows:
        for row in registry_rows:
            gca = row.get("gca_accession")
            gb_status = row.get("gb_status")
            genebuilder = row.get("genebuilder")
            print(f"\n Annotation already exists for GCA: {gca}")
            print(f"   Status     : {gb_status}")
            print(f"   Genebuilder: {genebuilder}")
        raise RuntimeError("Terminating: One or more assemblies are already annotated.")

    print("No existing annotation found for provided GCA(s).")
    return None

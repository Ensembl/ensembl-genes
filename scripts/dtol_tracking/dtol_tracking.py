import os.path, sys, getopt
import argparse
import requests
import pymysql


class text:
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


def mysql_fetch_data(query, database, host, port, user, password, multi):
    try:
        conn = pymysql.connect(
            host=host, user=user, passwd=password, port=port, database=database.strip()
        )
        cursor = conn.cursor()
        info = []
        if multi:
            for stmt in query.split(";"):
                if stmt.strip():
                    cursor.execute(stmt)
                    info.append(cursor.fetchall())
        else:
            cursor.execute(query)
            info = cursor.fetchall()

    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist on the server")
        else:
            print(err)

    cursor.close()
    conn.close()
    return info


def check_site(gca_list, site):
    on_list = []
    not_on_list = []
    if site == "rapid":
        species_list_url = "https://rapid.ensembl.org/info/about/species.html"
        species_list_response = requests.get(species_list_url)
    elif site == "projects":
        species_list_url = "https://projects.ensembl.org/darwin-tree-of-life/"
        species_list_response = requests.get(species_list_url)
    for gca in gca_list:
        if gca in species_list_response.text:
            on_list.append(gca)
        else:
            not_on_list.append(gca)
    return (on_list, not_on_list)


if __name__ == "__main__":

    # registry db details
    db_name = "gb_assembly_registry"
    db_host = "mysql-ens-genebuild-prod-1.ebi.ac.uk"
    db_port = 4527
    db_user = "ensro"
    db_pass = ""

    # rapid meta db details
    meta_db_name = "ensembl_metadata_qrp"
    meta_db_host = "mysql-ens-meta-prod-1.ebi.ac.uk"
    meta_db_port = 4483
    meta_db_user = "ensro"
    meta_db_pass = ""

    parser = argparse.ArgumentParser(
        description="Track status of DToL assemblies in Ensembl."
    )
    parser.add_argument(
        "-summary",
        help="Provide a summary of the status of DToL assemblies in Ensembl",
        required=False,
    )
    parser.add_argument(
        "-in_progress",
        help="Provide a summary of the status of DToL assemblies in Ensembl",
        required=False,
    )
    parser.add_argument(
        "-unannotated",
        help="Provide a list of GCAs for the DToL assemblies that have yet to be annotated.",
        required=False,
    )
    parser.add_argument(
        "-projects_missing",
        help="Provide a list of the core dbs for the Ensembl annotated DToL assemblies that have been release on rapid.ensembl.org",
        required=False,
    )

    args = parser.parse_args()
    summary = args.summary
    in_progress = args.in_progress
    unannotated = args.unannotated
    projects_missing = args.projects_missing

    dtol_count_query = 'SELECT count(*) FROM assembly JOIN meta USING(assembly_id) WHERE meta.assembly_group="dtol";SELECT count(*) FROM assembly JOIN meta USING(assembly_id) WHERE meta.assembly_group="dtol" AND meta.assembly_name LIKE "%alternate%";SELECT count(*) FROM assembly JOIN meta USING(assembly_id) WHERE meta.assembly_group="dtol" AND assembly.annotated_status="ensembl";SELECT count(*) FROM assembly JOIN meta USING(assembly_id) WHERE meta.assembly_group="dtol" AND assembly.annotated_status="ensembl" AND meta.assembly_name LIKE "%alternate%";'
    dtol_count_fetch = mysql_fetch_data(
        dtol_count_query,
        db_name,
        db_host,
        db_port,
        db_user,
        db_pass,
        1,
    )
    dtol_count = dtol_count_fetch[0][0][0]
    dtol_alternates = dtol_count_fetch[1][0][0]
    dtol_annotated = dtol_count_fetch[2][0][0]
    dtol_annotated_alternates = dtol_count_fetch[3][0][0]

    dtol_annotated_gcas_query = 'SELECT assembly.chain, assembly.version FROM assembly JOIN meta USING (assembly_id) WHERE meta.assembly_group="dtol" AND assembly.annotated_status="ensembl";'
    dtol_annotated_gcas_return = mysql_fetch_data(
        dtol_annotated_gcas_query,
        db_name,
        db_host,
        db_port,
        db_user,
        db_pass,
        0,
    )
    dtol_annotated_gcas = []
    for tuple in dtol_annotated_gcas_return:
        dtol_annotated_gcas.append(str(tuple[0]) + "." + str(tuple[1]))

    if summary:
        print(
            text.BOLD
            + "\nNumber of DToL assemblies in public archives: "
            + str(dtol_count)
            + text.END
            + " (of which are alternate assemblies: "
            + str(dtol_alternates)
            + ")\n"
            + text.BOLD
            + "Number of Ensembl annotated DToL assemblies: "
            + str(dtol_annotated)
            + ""
            + text.END
            + " (of which are alternate assemblies: "
            + str(dtol_annotated_alternates)
            + ")\n"
        )

    if in_progress:
        dtol_in_progress_gcas_query = 'SELECT assembly.chain, assembly.version FROM assembly JOIN meta USING (assembly_id) JOIN genebuild_status USING (assembly_id) WHERE meta.assembly_group="dtol" AND genebuild_status.progress_status="in_progress";'
        dtol_in_progress_gcas_return = mysql_fetch_data(
            dtol_in_progress_gcas_query,
            db_name,
            db_host,
            db_port,
            db_user,
            db_pass,
            0,
        )
        dtol_in_progress_gcas = []
        for tuple in dtol_in_progress_gcas_return:
            dtol_in_progress_gcas.append(str(tuple[0]) + "." + str(tuple[1]))

        print(
            text.BOLD
            + "Annotations for "
            + str(len(dtol_in_progress_gcas))
            + " DToL assemblies are currently in progress.\n"
            + text.END
            + "Assembly accessions:\n"
        )
        for gca in dtol_in_progress_gcas:
            print(gca)

    if unannotated:
        dtol_unannotated_gcas_query = 'SELECT assembly.chain, assembly.version FROM assembly JOIN meta USING (assembly_id) WHERE meta.assembly_group="dtol" AND assembly.annotated_status="unannotated";'
        dtol_unannotated_gcas_return = mysql_fetch_data(
            dtol_unannotated_gcas_query,
            db_name,
            db_host,
            db_port,
            db_user,
            db_pass,
            0,
        )
        dtol_unannotated_gcas = []
        for tuple in dtol_unannotated_gcas_return:
            dtol_unannotated_gcas.append(str(tuple[0]) + "." + str(tuple[1]))
        print(
            text.BOLD
            + str(len(dtol_unannotated_gcas))
            + " DToL assemblies are yet to be annotated.\n"
            + text.END
            + "Assembly accessions:\n"
        )
        for gca in dtol_unannotated_gcas:
            print(gca)

    if projects_missing:
        projects_released, projects_unreleased = check_site(
            dtol_annotated_gcas, "projects"
        )
        unreleased_dbs = []
        for gca in projects_unreleased:
            unreleased_gca_query = (
                'SELECT dbname FROM assembly JOIN genome USING (assembly_id) JOIN genome_database USING (genome_id) WHERE genome.data_release_id=(SELECT MAX(data_release_id) FROM genome) AND assembly.assembly_accession="'
                + gca
                + '";'
            )
            unreleased_gca_return = mysql_fetch_data(
                unreleased_gca_query,
                meta_db_name,
                meta_db_host,
                meta_db_port,
                meta_db_user,
                meta_db_pass,
                0,
            )

            unreleased_dbs.append(unreleased_gca_return[0][0])

        print(
            text.BOLD
            + str(len(unreleased_dbs))
            + " DToL assemblies have been annotated but do not appear on https://projects.ensembl.org/darwin-tree-of-life/.\n"
            + text.END
            + "Core databases:\n"
        )
        for db in unreleased_dbs:
            print(db)

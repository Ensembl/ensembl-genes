import os.path, sys, getopt
import pymysql
import requests
import argparse
from ftplib import FTP
from ftplib import error_temp, error_perm
import re
from collections import OrderedDict

class EnsemblFTP:
    def __init__(self):
        # Initialize both FTP connections
        self.ensembl_ftp = FTP("ftp.ensembl.org")
        self.ensembl_ftp.login()
        self.ebi_ftp = FTP("ftp.ebi.ac.uk")
        self.ebi_ftp.login()

        # Define paths for both FTPs
        self.ensembl_ftp_path = "https://ftp.ensembl.org/"
        self.ebi_ftp_path = "https://ftp.ebi.ac.uk/"

    def return_to_root(self, ftp_connection):
        """Return to the root directory to reset the FTP session state."""
        ftp_connection.cwd("/")  # Reset the FTP directory to root

    def check_for_file(self, species_name, prod_name, accession, source, file_type):
        # Choose the correct FTP connection based on file_type
        if file_type == "repeatmodeler":
            ftp_connection = self.ebi_ftp
            ftp_path = self.ebi_ftp_path
            path = (
                "pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/"
                + species_name
                + "/"
            )
            file_name = accession + ".repeatmodeler.fa"
        elif file_type == "busco":
            ftp_connection = self.ensembl_ftp
            ftp_path = self.ensembl_ftp_path
            path = (
                "pub/rapid-release/species/"
                + species_name
                + "/"
                + accession
                + "/"
                + source
                + "/statistics/"
            )
            file_name = prod_name + "_busco_short_summary.txt"
        
        print(f"Checking path: {path}")  # Debugging to check the path

        try:
            # Ensure the FTP connection is in the root directory before changing to the desired directory
            self.return_to_root(ftp_connection)
            ftp_connection.cwd(path)  # Change to the desired directory
            if file_name in ftp_connection.nlst():
                return ftp_path + path + file_name
            else:
                return 0
        except error_perm as e:
            # Ignore the "550 Failed to change directory." error, silently return 0
            if "550" in str(e) and "Failed to change directory" in str(e):
                return 0
            else:
                print(f"FTP permission error: {e}")
                return 0
        except error_temp as e:
            if "421" in str(e):
                print("Too many users connected. Retrying...")
                time.sleep(5)  # Wait for 5 seconds before retrying
                return self.check_for_file(species_name, prod_name, accession, source, file_type)
            else:
                print(f"FTP temporary error: {e}")
                return 0
        except Exception as e:
            print(f"Error: {e}")
            return 0

    def close_connections(self):
        """Close both FTP connections."""
        if self.ensembl_ftp:
            self.ensembl_ftp.quit()
        if self.ebi_ftp:
            self.ebi_ftp.quit()

def mysql_fetch_data(query, database, host, port, user, password):
    try:
        conn = pymysql.connect(
            host=host, user=user, passwd=password, port=port, database=database.strip()
        )

        cursor = conn.cursor()
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


def write_yaml(info_dict, icon, yaml_out, project, use_server, alternate, ftp_client):
    # there are species on main for which the upper case production name is used in the url instead of the upper case species name
    prod_url_list = ["bos_taurus_hybrid", "bos_indicus_hybrid"]

    assembly_name = info_dict["assembly.name"].replace(" ", "_")
    date = info_dict["genebuild.last_geneset_update"].replace("-", "_")
    species_name = info_dict["species.scientific_name"].replace(" ", "_")
    rm_species_name = (
        info_dict["species.scientific_name"].replace(" ", "_")
    ).lower()  # set this before adding strain/breed for the repeatmodeler ftp path
    if "species.strain" in info_dict:
        if info_dict["species.strain"] != "reference":
            info_dict["species.scientific_name"] = (
                info_dict["species.scientific_name"]
                + " ("
                + info_dict["species.strain"]
                + ")"
            )
    lc_species_name = (info_dict["species.scientific_name"].replace(" ", "_")).lower()
    uc_prod_name = (info_dict["species.production_name"]).capitalize()

    assembly_report_url = (
        "https://www.ncbi.nlm.nih.gov/assembly/" + info_dict["assembly.accession"]
    )
    assembly_report_response = requests.get(assembly_report_url)
    submitter_match = re.search(
        "Submitter: </dt><dd>([^<]*)</dd><dt>", assembly_report_response.text
    )
    if submitter_match:
        submitter = submitter_match.group(1)
    else:
        submitter = "unknown"
    # for the chicken reference the submitter info is missing, set it manually
    if info_dict["assembly.accession"] == "GCA_000002315.5":
        submitter = "Genome Reference Consortium"

    # 01-12-22: Add annotation source for new rapid FTP structur
    if "species.annotation_source" in info_dict:
        source = info_dict["species.annotation_source"]
    else:
        source = 'ensembl'

    if use_server == "rapid":
        ftp_base = "https://ftp.ensembl.org/pub/rapid-release/species"

        yaml = "- species: " + info_dict["species.scientific_name"] + "\n"

        if project in ("vgp", "dtol", "erga", "cbp", "bge", "asg"):
            yaml += "  image: " + icon + "\n"
        else:
            yaml += "  submitted_by: " + submitter + "\n"

        yaml += "  accession: " + info_dict["assembly.accession"] + "\n"

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            # 07-06-22: Add column for "Annotation method" so user can clearly see whether Ensembl genebuild or BRAKER2 annotation
            if "genebuild.method_display" in info_dict:
                yaml += "  annotation_method: " + info_dict["genebuild.method_display"] + "\n"
            else:
                yaml += "  annotation_method: BRAKER2\n"
                
        yaml += (
            "  annotation_gtf: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/"
            + source
            + "/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-genes.gtf.gz\n"
        )
        yaml += (
            "  annotation_gff3: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/"
            + source
            + "/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-genes.gff3.gz\n"
        )
        yaml += (
            "  proteins: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/"
            + source
            + "/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-pep.fa.gz\n"
        )
        yaml += (
            "  transcripts: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/"
            + source
            + "/geneset/"
            + date
            + "/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-"
            + date
            + "-cdna.fa.gz\n"
        )
        yaml += (
            "  softmasked_genome: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/"
            + source
            + "/genome/"
            + species_name
            + "-"
            + info_dict["assembly.accession"]
            + "-softmasked.fa.gz\n"
        )
        rm_file = ftp_client.check_for_file(
            rm_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler",
        )
        if rm_file:
            yaml += "  repeat_library: " + rm_file + "\n"

        yaml += (
            "  ftp_dumps: "
            + ftp_base
            + "/"
            + species_name
            + "/"
            + info_dict["assembly.accession"]
            + "/"
            + source
            + "\n"
        )

        main_species_url = "http://www.ensembl.org/info/about/species.html"
        main_species_response = requests.get(main_species_url)
        if info_dict["assembly.accession"] in main_species_response.text:
            yaml += (
                "  ensembl_link: https://www.ensembl.org/"
                + species_name
                + "/Info/Index\n"
            )
        else:
            yaml += (
                "  rapid_link: https://rapid.ensembl.org/"
                + uc_prod_name
                + "/Info/Index\n"
            )

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            # 10-05-22: Add column for busco score files for DToL only (the ftp will soon be moved from temp DToL FTP to RR FTP)
            busco_file = ftp_client.check_for_file(
                species_name,
                info_dict["species.production_name"],
                info_dict["assembly.accession"],
                source,
                "busco",
            )
            if busco_file:
                yaml += "  busco_score: " + busco_file + "\n"
            # 10-06-22: Add column for alternates (where available) for DToL only
            if alternate:
                alternate_url = "https://rapid.ensembl.org/" + alternate + "/Info/Index"
                yaml += "  alternate: " + alternate_url + "\n"

    elif use_server == "main":
        # 12-03-21: GENE-SWitCH project is frozen on Ensembl version 102 for the time being!
        if project == "geneswitch":
            release = "release-102"
            release_number = "102"
        else:
            release = "release-" + info_dict["schema_version"]
            release_number = info_dict["schema_version"]
        ftp_base = "https://ftp.ensembl.org/pub/" + release

        yaml = "- species: " + info_dict["species.scientific_name"] + "\n"
        if project in ("vgp", "dtol", "erga", "cbp", "bge", "asg"):
            yaml += "  image: " + icon + "\n"
        else:
            yaml += "  submitted_by: " + submitter + "\n"

        yaml += "  accession: " + info_dict["assembly.accession"] + "\n"

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            # 07-06-22: Add column for "Annotation method" so user can clearly see whether Ensembl genebuild or BRAKER2 annotation
            yaml += "  annotation_method: Ensembl genebuild\n"

        # 01-12-22: Add annotation source for new rapid FTP structure (not needed for main, but is used in check_for_file function for retrieving BUSCO file - should clean this up!)
        source = 'ensembl'
            
        yaml += (
            "  annotation_gtf: "
            + ftp_base
            + "/gtf/"
            + info_dict["species.production_name"]
            + "/"
            + uc_prod_name
            + "."
            + assembly_name
            + "."
            + release_number
            + ".gtf.gz\n"
        )
        yaml += (
            "  annotation_gff3: "
            + ftp_base
            + "/gff3/"
            + info_dict["species.production_name"]
            + "/"
            + uc_prod_name
            + "."
            + assembly_name
            + "."
            + release_number
            + ".gff3.gz\n"
        )
        yaml += (
            "  proteins: "
            + ftp_base
            + "/fasta/"
            + info_dict["species.production_name"]
            + "/pep/"
            + uc_prod_name
            + "."
            + assembly_name
            + ".pep.all.fa.gz\n"
        )
        yaml += (
            "  transcripts: "
            + ftp_base
            + "/fasta/"
            + info_dict["species.production_name"]
            + "/cdna/"
            + uc_prod_name
            + "."
            + assembly_name
            + ".cdna.all.fa.gz\n"
        )
        yaml += (
            "  softmasked_genome: "
            + ftp_base
            + "/fasta/"
            + info_dict["species.production_name"]
            + "/dna/"
            + uc_prod_name
            + "."
            + assembly_name
            + ".dna_sm.toplevel.fa.gz\n"
        )
        rm_file = ftp_client.check_for_file(
            rm_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler",
        )
        if rm_file:
            yaml += "  repeat_library: " + rm_file + "\n"

        yaml += "  ftp_dumps: " + ftp_base + "\n"
        if info_dict["species.production_name"] in prod_url_list:
            yaml += (
                "  ensembl_link: https://www.ensembl.org/"
                + uc_prod_name
                + "/Info/Index\n"
            )
        elif project == "geneswitch":
            # 12-03-21: GENE-SWitCH project is frozen on Ensembl version 102 for the time being!
            yaml += (
                "  ensembl_link: https://e102.ensembl.org/"
                + species_name
                + "/Info/Index\n"
            )
        else:
            yaml += (
                "  ensembl_link: https://www.ensembl.org/"
                + species_name
                + "/Info/Index\n"
            )

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            # 10-05-22: Add column for busco score files for DToL only (the ftp will soon be moved from temp DToL FTP to RR FTP)
            busco_file = ftp_client.check_for_file(
                lc_species_name,
                info_dict["species.production_name"],
                info_dict["assembly.accession"],
                source,
                "busco",
            )
            if busco_file:
                yaml += "  busco_score: " + busco_file + "\n"
            # 10-06-22: Add column for alternates (where available) for DToL only
            if alternate:
                alternate_url = "https://rapid.ensembl.org/" + alternate + "/Info/Index"
                yaml += "  alternate: " + alternate_url + "\n"
                
    print(yaml, file=yaml_out)


if __name__ == "__main__":

    server_dict = {
        "rapid": {
            "db_host": "mysql-ens-sta-5.ebi.ac.uk",
            "db_port": 4684,
            "db_user": "ensro",
            "db_pass": "",
        },
        "main": {
            "db_host": "mysql-ens-mirror-1.ebi.ac.uk",
            "db_port": 4240,
            "db_user": "ensro",
            "db_pass": "",
        },
        "meta": {
            "db_host": "mysql-ens-meta-prod-1.ebi.ac.uk",
            "db_port": 4483,
            "db_user": "ensro",
            "db_pass": "",
        },
    }

    icon_list = open("icons.txt").readlines()
    icon_dict = {}
    for line in icon_list:
        icon_dict[line.split()[0]] = line.split()[1]

    parser = argparse.ArgumentParser(
        description="Create species.yaml file for a given project page."
    )
    parser.add_argument(
        "-f",
        "--db_file",
        help="Name for file containing list of VGP databases on the Rapid Release or Main server",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--project",
        choices=["aquafaang", "asg", "bge", "bovreg", "cbp", "dtol", "erga", "geneswitch", "vgp"],
        help="Name of the project this set of database belongs to",
        required=True,
    )
    args = parser.parse_args()

    db_list = open(args.db_file).readlines()
    project = args.project
    yaml_out = open(project + "_species.yaml", "w")

    sorted_db_list = sorted(db_list)

    if project == "aquafaang":
        # move danio rerio reference dbs to end of the list
        for db in sorted_db_list:
            if "danio_rerio_core" in db:
                sorted_db_list.append(sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "bovreg":
        # move bos taurus reference dbs to top of the list
        for db in sorted_db_list:
            if "bos_taurus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "geneswitch":
        # move sus scrofa or gallus gallus reference dbs to top of the list
        for db in sorted_db_list:
            if "sus_scrofa_core" in db or "gallus_gallus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

    #open the FTP connection
    ftp_client = EnsemblFTP()

    for db in sorted_db_list:
        db = db.strip()
        # check if the db is on rapid
        conn = pymysql.connect(
            host=server_dict["rapid"]["db_host"],
            user=server_dict["rapid"]["db_user"],
            passwd=server_dict["rapid"]["db_pass"],
            port=server_dict["rapid"]["db_port"],
        )

        cur = conn.cursor()
        exists_query = "SHOW DATABASES LIKE '" + db.strip() + "'"
        exists_rapid = cur.execute(exists_query)

        if exists_rapid:
            use_server = "rapid"
        # if it's not on the rapid server check main
        else:
            conn = pymysql.connect(
                host=server_dict["main"]["db_host"],
                user=server_dict["main"]["db_user"],
                passwd=server_dict["main"]["db_pass"],
                port=server_dict["main"]["db_port"],
            )

            cur = conn.cursor()
            exists_query = "SHOW DATABASES LIKE '" + db.strip() + "'"
            exists_main = cur.execute(exists_query)
            if exists_main:
                use_server = "main"
            else:
                raise Exception(
                    "Unable to find database "
                    + db.strip()
                    + " on rapid or main servers!"
                )

        if use_server:
            # retrieve the species name, assembly accession and assembly name from the database
            info_query = "SELECT meta_key,meta_value FROM meta WHERE meta_key in ('species.scientific_name','assembly.accession','assembly.name','species.production_name','species.strain','schema_version','genebuild.last_geneset_update','species.annotation_source') OR meta_key like 'genebuild.method%'"
            info = mysql_fetch_data(
                info_query,
                db,
                server_dict[use_server]["db_host"],
                server_dict[use_server]["db_port"],
                server_dict[use_server]["db_user"],
                server_dict[use_server]["db_pass"],
            )

            info_dict = {}
            for tuple in info:
                info_dict[tuple[0]] = tuple[1]

            # check if the assembly has an alternate via the metadata db
            alternate_assembly_name = info_dict["assembly.name"]+"_alternate_haplotype"
            alternate_query = (
                "SELECT organism.url_name from assembly join genome using (assembly_id) join organism using (organism_id) where assembly.assembly_name='" + alternate_assembly_name + "' limit 1;"
            )
            alternate_fetch = mysql_fetch_data(
                alternate_query,
                "ensembl_metadata_qrp",
                server_dict["meta"]["db_host"],
                server_dict["meta"]["db_port"],
                server_dict["meta"]["db_user"],
                server_dict["meta"]["db_pass"],
            )
            try:
                alternate = alternate_fetch[0][0]
            except:
                alternate = ""
                
            # retrieve the species classification info from the database. in order to assign an icon
            class_query = (
                "SELECT meta_value FROM meta WHERE meta_key='species.classification'"
            )
            classifications = mysql_fetch_data(
                class_query,
                db,
                server_dict[use_server]["db_host"],
                server_dict[use_server]["db_port"],
                server_dict[use_server]["db_user"],
                server_dict[use_server]["db_pass"],
            )

            class_list = []
            for tuple in classifications:
                class_list.append(tuple[0])

            icon = "Metazoa.png"
            chordate = 0
            if "Chordata" in class_list:
                chordate = 1
            for classification in class_list:
                try:
                    icon = icon_dict[classification]
                    break
                except KeyError:
                    continue
            if chordate and icon == "Metazoa.png":
                icon = "Chordates.png"

            write_yaml(info_dict, icon, yaml_out, project, use_server, alternate, ftp_client)
        else:
            print(
                "Could not find database "
                + db.strip()
                + " on mirror or rapid release servers!\n"
            )
    ftp_client.close_connections()

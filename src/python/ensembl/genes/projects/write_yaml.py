#!/usr/bin/env python3

import argparse
import json
import os
import re
import sys
import time
from ftplib import FTP, error_perm, error_temp
from typing import Dict, List, Tuple, Optional

import pymysql
import requests


class EnsemblFTP:
    """
    A class to handle FTP connections and file lookups on Ensembl FTP servers.
    """

    def __init__(self) -> None:
        """
        Initialize FTP connections to the Ensembl and EBI servers.
        """
        self.ensembl_ftp = FTP("ftp.ensembl.org")
        self.ensembl_ftp.login()
        self.ebi_ftp = FTP("ftp.ebi.ac.uk")
        self.ebi_ftp.login()

        self.ensembl_ftp_path = "https://ftp.ensembl.org/"
        self.ebi_ftp_path = "https://ftp.ebi.ac.uk/"

    def return_to_root(self, ftp_connection: FTP) -> None:
        """
        Return to the root directory of the FTP connection to reset its state.
        
        Parameters
        ----------
        ftp_connection : FTP
            The FTP connection object.
        """
        ftp_connection.cwd("/")

    def check_for_file(
            self,
            species_name: str,
            prod_name: str,
            accession: str,
            source: str,
            file_type: str
    ) -> str:
        """
        Check if a specific file (repeatmodeler or busco) is present on the FTP site.
        
        Parameters
        ----------
        species_name : str
        The species scientific name (used in paths).
        prod_name : str
        The production name of the species (used in file naming for busco).
        accession : str
        The assembly accession ID.
        source : str
        Annotation source ('ensembl' or other).
        file_type : str
        The type of file to check: 'repeatmodeler' or 'busco'.
        
        Returns
        -------
        str
        The full FTP URL if the file is found, otherwise an empty string.
        """
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
            # Two possible file name formats:
            file_name_protein = prod_name + "_protein_busco_short_summary.txt"
            file_name_alternative = prod_name + "_busco_short_summary.txt"
        else:
            return ""
    
        # For logging/debugging, print the expected primary file name URL
        #if file_type == "busco":
        #    print(ftp_path + path + file_name_protein)
        #else:
        #    print(ftp_path + path + file_name)
        
        try:
            # Reset FTP directory state
            self.return_to_root(ftp_connection)
            ftp_connection.cwd(path)
            files_list = ftp_connection.nlst()
            
            if file_type == "busco":
                if file_name_protein in files_list:
                    return ftp_path + path + file_name_protein
                elif file_name_alternative in files_list:
                    return ftp_path + path + file_name_alternative
                else:
                    return ""
            else:
                if file_name in files_list:
                    return ftp_path + path + file_name
                else:
                    return ""
        except error_perm as e:
            # If directory doesn't exist, just return empty string
            if "550" in str(e) and "Failed to change directory" in str(e):
                return ""
            else:
                print(f"FTP permission error: {e}", file=sys.stderr)
                return ""
        except error_temp as e:
            if "421" in str(e):
                print("Too many users connected. Retrying...", file=sys.stderr)
                time.sleep(5)
                return self.check_for_file(species_name, prod_name, accession, source, file_type)
            else:
                print(f"FTP temporary error: {e}", file=sys.stderr)
                return ""
        except Exception as e:
            print(f"Error while checking FTP file: {e}", file=sys.stderr)
            return ""


    def check_pre_release_file(
        self,
        species_name: str,
        accession: str,
        extension: str
    ) -> str:
        """
        Look in the 'pre-release' FTP folder for first file ending in `extension`
        (e.g. '.gtf.gz', '.gff3.gz', '.dna.softmasked.fa.gz').
        Return its HTTPS URL or empty string.
        """
        ftp = self.ebi_ftp
        base = self.ebi_ftp_path
        path = f"pub/databases/ensembl/pre-release/{species_name}/{accession}/"
        try:
            self.return_to_root(ftp)
            ftp.cwd(path)
            for fname in ftp.nlst():
                if fname.lower().endswith(extension):
                    return base + path + fname
        except error_perm:
            return ""
        return ""

    def close_connections(self) -> None:
        """
        Close both FTP connections.
        """
        if self.ensembl_ftp:
            self.ensembl_ftp.quit()
        if self.ebi_ftp:
            self.ebi_ftp.quit()



            
def mysql_fetch_data(
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str
) -> Tuple:
    """
    Fetch data from a MySQL database using a given query.

    Parameters
    ----------
    query : str
        The SQL query to execute.
    database : str
        The database name.
    host : str
        The database host.
    port : int
        The database port.
    user : str
        The username for the database.
    password : str
        The password for the database.

    Returns
    -------
    Tuple
        A tuple of results fetched from the database query.
    """
    conn = pymysql.connect(
        host=host, user=user, passwd=password, port=port, database=database.strip()
    )
    cursor = conn.cursor()
    cursor.execute(query)
    info = cursor.fetchall()
    cursor.close()
    conn.close()
    return info

def check_url_status(url):
    """
    Checks if a given URL is reachable.
    
    Args:
        url (str): The URL to check.
    
    Returns:
        bool: True if the URL returns a 200 status code, False otherwise.
    """
    try:
        response = requests.head(url, allow_redirects=True, timeout=5)  # Use HEAD for efficiency
        return response.status_code == 200
    except requests.RequestException as e:
        print(f"Error checking URL {url}: {e}")
        return False

def write_yaml(
        info_dict: Dict[str, str],
    icon: str,
    yaml_out,
    project: str,
    use_server: str,
    alternate: str,
        guuid: str,
    ftp_client: EnsemblFTP
) -> None:
    """
    Write YAML content for a species entry based on the provided metadata and project type.

    Parameters
    ----------
    info_dict : Dict[str, str]
        Dictionary containing species metadata.
    icon : str
        Icon filename for the species.
    yaml_out : file-like object
        The output file handle for writing the YAML.
    project : str
        The project name (e.g., 'aquafaang', 'asg', 'bge', ...).
    use_server : str
        Which server the data was fetched from ('st5', 'st6' or 'main').
    alternate : str
        Alternate assembly name, if available.
    ftp_client : EnsemblFTP
        An initialized FTP client instance for checking extra files.
    """
    prod_url_list = ["bos_taurus_hybrid", "bos_indicus_hybrid"]

    assembly_name = info_dict["assembly.name"].replace(" ", "_")
    date = info_dict["genebuild.last_geneset_update"].replace("-", "_")
    species_name = info_dict["species.scientific_name"].replace(" ", "_")
    rm_species_name = species_name.lower()

    if "species.strain" in info_dict and info_dict["species.strain"] != "reference":
        info_dict["species.scientific_name"] += f" ({info_dict['species.strain']})"

    lc_species_name = info_dict["species.scientific_name"].replace(" ", "_").lower()
    uc_prod_name = info_dict["species.production_name"].capitalize()

    # Get submitter from assembly report
    assembly_report_url = f"https://www.ncbi.nlm.nih.gov/assembly/{info_dict['assembly.accession']}"
    assembly_report_response = requests.get(assembly_report_url)
    submitter_match = re.search(
        r"Submitter: </dt><dd>([^<]*)</dd><dt>", assembly_report_response.text
    )
    submitter = submitter_match.group(1) if submitter_match else "unknown"
    if info_dict["assembly.accession"] == "GCA_000002315.5":
        submitter = "Genome Reference Consortium"

    source = info_dict.get("species.annotation_source", "ensembl")

    # Start building YAML content
    yaml = f"- species: {info_dict['species.scientific_name']}\n"

    if use_server == "st5" or use_server == "st6":
        ftp_base = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"

        if project in ("vgp", "dtol", "erga", "cbp", "bge", "asg"):
            yaml += f"  image: {icon}\n"
        else:
            yaml += f"  submitted_by: {submitter}\n"

        yaml += f"  accession: {info_dict['assembly.accession']}\n"

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            # Add annotation method if available, else default to BRAKER2
            method = info_dict.get("genebuild.method_display", "BRAKER2")
            yaml += f"  annotation_method: {method}\n"

        # ==== annotation_gtf ====
        gtf = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/" \
              f"{source}/geneset/{date}/genes.gtf.gz"
        if check_url_status(gtf):
            yaml += f"  annotation_gtf: {gtf}\n"
        else:
            fb = ftp_client.check_pre_release_file(species_name,
                                                  info_dict["assembly.accession"],
                                                  ".gtf.gz")
            if fb:
                yaml += f"  annotation_gtf: {fb}\n"
        # ==== annotation_gff3: try gzipped, then uncompressed, then pre‑release ====
        gff_gz = (
            f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/"
            f"{source}/geneset/{date}/genes.gff3.gz"
        )
        gff    = (
            f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/"
            f"{source}/geneset/{date}/genes.gff3"
        )
        if check_url_status(gff_gz):
            yaml += f"  annotation_gff3: {gff_gz}\n"
        elif check_url_status(gff):
            yaml += f"  annotation_gff3: {gff}\n"
        else:
            # fallback: look in pre‑release for gzipped first, then uncompressed
            fb = ftp_client.check_pre_release_file(
                species_name,
                info_dict["assembly.accession"],
                ".gff3.gz"
            )
            if not fb:
                fb = ftp_client.check_pre_release_file(
                    species_name,
                    info_dict["assembly.accession"],
                    ".gff3"
                )
            if fb:
                yaml += f"  annotation_gff3: {fb}\n"
        # ==== proteins ====
        pep = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/" \
              f"{source}/geneset/{date}/pep.fa.gz"
        if check_url_status(pep):
            yaml += f"  proteins: {pep}\n"
        # ==== transcripts ====
        cdna = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/" \
               f"{source}/geneset/{date}/cdna.fa.gz"
        if check_url_status(cdna):
            yaml += f"  transcripts: {cdna}\n"
        # ==== softmasked genome ====
        soft = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/genome/" \
               f"softmasked.fa.gz"
        if check_url_status(soft):
            yaml += f"  softmasked_genome: {soft}\n"
        else:
            fb = ftp_client.check_pre_release_file(species_name,
                                                  info_dict["assembly.accession"],
                                                  ".dna.softmasked.fa.gz")
            if fb:
                yaml += f"  softmasked_genome: {fb}\n"
        rm_file = ftp_client.check_for_file(
            rm_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler"
        )
        if rm_file:
            yaml += f"  repeat_library: {rm_file}\n"

        yaml += f"  ftp_dumps: {ftp_base}/{species_name}/{info_dict['assembly.accession']}/\n"

        main_species_url = "http://www.ensembl.org/info/about/species.html"
        main_species_response = requests.get(main_species_url)
        if info_dict["assembly.accession"] in main_species_response.text:
            yaml += f"  ensembl_link: https://www.ensembl.org/{species_name}/Info/Index\n"
        else:
            beta_link = f"https://beta.ensembl.org/species/{guuid}"
            if check_url_status(beta_link):
                yaml += f"  beta_link: https://beta.ensembl.org/species/{guuid}\n"
            else:
                yaml += f"  beta_link: Coming soon!\n"

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            try:
                yaml += f"  busco_score: {info_dict['genebuild.busco']}\n"
            except KeyError:
                busco_file = ftp_client.check_for_file(
                    species_name,
                    info_dict["species.production_name"],
                    info_dict["assembly.accession"],
                    source,
                    "busco"
                )
                if busco_file:
                    yaml += f"  busco_score: {busco_file}\n"
            if alternate:
                alternate_url = f"https://rapid.ensembl.org/{alternate}/Info/Index"
                yaml += f"  alternate: {alternate_url}\n"

    elif use_server == "main":
        # Handle "geneswitch" project special case with release 102
        if project == "geneswitch":
            release = "release-102"
            release_number = "102"
        else:
            release = "release-" + info_dict["schema_version"]
            release_number = info_dict["schema_version"]

        ftp_base = f"https://ftp.ensembl.org/pub/{release}"

        if project in ("vgp", "dtol", "erga", "cbp", "bge", "asg"):
            yaml += f"  image: {icon}\n"
        else:
            yaml += f"  submitted_by: {submitter}\n"

        yaml += f"  accession: {info_dict['assembly.accession']}\n"

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            # Annotation method is always Ensembl genebuild on main
            yaml += "  annotation_method: Ensembl genebuild\n"

        source = 'ensembl'  # for consistency with busco retrieval logic

        yaml += (
            f"  annotation_gtf: {ftp_base}/gtf/{info_dict['species.production_name']}/"
            f"{uc_prod_name}.{assembly_name}.{release_number}.gtf.gz\n"
        )
        yaml += (
            f"  annotation_gff3: {ftp_base}/gff3/{info_dict['species.production_name']}/"
            f"{uc_prod_name}.{assembly_name}.{release_number}.gff3.gz\n"
        )
        yaml += (
            f"  proteins: {ftp_base}/fasta/{info_dict['species.production_name']}/pep/"
            f"{uc_prod_name}.{assembly_name}.pep.all.fa.gz\n"
        )
        yaml += (
            f"  transcripts: {ftp_base}/fasta/{info_dict['species.production_name']}/cdna/"
            f"{uc_prod_name}.{assembly_name}.cdna.all.fa.gz\n"
        )
        yaml += (
            f"  softmasked_genome: {ftp_base}/fasta/{info_dict['species.production_name']}/dna/"
            f"{uc_prod_name}.{assembly_name}.dna_sm.toplevel.fa.gz\n"
        )
        rm_file = ftp_client.check_for_file(
            rm_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler"
        )
        if rm_file:
            yaml += f"  repeat_library: {rm_file}\n"

        yaml += f"  ftp_dumps: {ftp_base}\n"
        if info_dict["species.production_name"] in prod_url_list:
            yaml += f"  ensembl_link: https://www.ensembl.org/{uc_prod_name}/Info/Index\n"
        elif project == "geneswitch":
            # geneswitch project frozen at e102
            yaml += f"  ensembl_link: https://e102.ensembl.org/{species_name}/Info/Index\n"
        else:
            yaml += f"  ensembl_link: https://www.ensembl.org/{species_name}/Info/Index\n"

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            try:
                yaml += f"  busco_score: {info_dict['genebuild.busco']}\n"
            except KeyError:
                busco_file = ftp_client.check_for_file(
                    lc_species_name,
                    info_dict["species.production_name"],
                    info_dict["assembly.accession"],
                    source,
                    "busco"
                )
                if busco_file:
                    yaml += f"  busco_score: {busco_file}\n"
            if alternate:
                alternate_url = f"https://rapid.ensembl.org/{alternate}/Info/Index"
                yaml += f"  alternate: {alternate_url}\n"

    print(yaml, file=yaml_out)

import pymysql
import json

def check_database_on_server(db, server_key, server_dict):
    """
    Checks if a database exists on a given server.

    Args:
        db (str): The name of the database to check.
        server_key (str): The key of the server in the config.
        server_dict (dict): Dictionary containing server connection details.

    Returns:
        bool: True if the database exists, False otherwise.
    """
    try:
        conn = pymysql.connect(
            host=server_dict[server_key]["db_host"],
            user=server_dict[server_key]["db_user"],
            passwd=server_dict[server_key]["db_pass"],
            port=server_dict[server_key]["db_port"],
        )
        with conn.cursor() as cur:
            #print(f"Checking for DB '{db}' on server '{server_key}'", file=sys.stderr)
            cur.execute("SELECT SCHEMA_NAME FROM information_schema.SCHEMATA WHERE SCHEMA_NAME = %s", (db,))
            result = cur.fetchone()
            #print(f"Query result on {server_key}: {result}", file=sys.stderr)
            return result is not None
        
    except pymysql.MySQLError as e:
        print(f"Error connecting to {server_key}: {e}")
        return False
    finally:
        if 'conn' in locals() and conn:
            conn.close()
            
            
def find_database_server(db, server_dict):
    """
    Determines which server the given database exists on.
    
    Args:
        db (str): The name of the database to check.
        server_dict (dict): Dictionary containing server connection details from the config file.

    Returns:
        str: The name of the server where the database is found.

    Raises:
        Exception: If the database is not found on any of the servers.
    """
    # Check servers in order: st5 -> st6 -> main
    for server_key in ["st5", "st6", "main"]:
        if check_database_on_server(db, server_key, server_dict):
            return server_key
        
    # If no server found, raise an error
    raise Exception(f"Unable to find database {db} on any configured servers!")
        
    
def main() -> None:
    """
    Main function that parses arguments, reads server config from JSON,
    queries MySQL databases, and writes YAML content.
    """
    parser = argparse.ArgumentParser(
        description="Create a species.yaml file for a given project page."
    )
    parser.add_argument(
        "-f",
        "--db_file",
        help="Name of the file containing list of databases on Beta or Main servers and their Genome UUIDs",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--project",
        choices=["aquafaang", "asg", "bge", "bovreg", "cbp", "dtol", "erga", "geneswitch", "vgp"],
        help="Name of the project this set of databases belongs to",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--config_file",
        help="Path to the JSON configuration file containing server connection details",
        required=True,
    )
    args = parser.parse_args()
    project = args.project

    # Load server connection details from JSON
    with open(args.config_file) as config_f:
        server_dict = json.load(config_f)

    # Read icons
    icon_dict = {}
    with open("icons.txt") as icon_file:
        for line in icon_file:
            fields = line.split()
            if len(fields) >= 2:
                icon_dict[fields[0]] = fields[1]

    # Read database list
    with open(args.db_file) as db_file:
        db_list = db_file.read().strip().split("\n")

    sorted_db_list = sorted(db_list)

    # Special rearrangements for certain projects
    if project == "aquafaang":
        # Move danio rerio reference DBs to the end
        for db in sorted_db_list:
            if "danio_rerio_core" in db:
                sorted_db_list.append(sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "bovreg":
        # Move bos taurus reference DBs to the top
        for db in sorted_db_list:
            if "bos_taurus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "geneswitch":
        # Move sus scrofa or gallus gallus reference DBs to top
        for db in sorted_db_list:
            if "sus_scrofa_core" in db or "gallus_gallus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

    # Initialize FTP connection
    ftp_client = EnsemblFTP()

    # Open the output YAML file
    with open(f"{project}_species.yaml", "w") as yaml_out:
        #will replace with guiid eventually
        for line in sorted_db_list:
            db = line.split('\t')[0].strip()
            guuid = line.split('\t')[1].strip()
            use_server = find_database_server(db, server_dict)

            # Retrieve metadata from the chosen server
            # update this with query from the metadata db
            info_query = (
                "SELECT meta_key, meta_value FROM meta "
                "WHERE meta_key IN "
                "('species.scientific_name','assembly.accession','assembly.name',"
                "'species.production_name','species.strain','schema_version',"
                "'genebuild.last_geneset_update','species.annotation_source','genebuild.busco') "
                "OR meta_key LIKE 'genebuild.method%'"
            )
            info = mysql_fetch_data(
                info_query,
                db,
                server_dict[use_server]["db_host"],
                server_dict[use_server]["db_port"],
                server_dict[use_server]["db_user"],
                server_dict[use_server]["db_pass"],
            )

            info_dict = {row[0]: row[1] for row in info}

            # Check if assembly has an alternate
            # update this with query from the metadata db
            alternate_assembly_name = info_dict["assembly.name"] + "_alternate_haplotype"
            alternate_query = (
                "SELECT organism.url_name "
                "FROM assembly "
                "JOIN genome USING (assembly_id) "
                "JOIN organism USING (organism_id) "
                f"WHERE assembly.assembly_name='{alternate_assembly_name}' LIMIT 1;"
            )
            alternate_fetch = mysql_fetch_data(
                alternate_query,
                "ensembl_metadata_qrp",
                server_dict["meta"]["db_host"],
                server_dict["meta"]["db_port"],
                server_dict["meta"]["db_user"],
                server_dict["meta"]["db_pass"],
            )
            alternate = alternate_fetch[0][0] if alternate_fetch else ""

            # Retrieve classification info
            # replace this with call to datasets
            class_query = "SELECT meta_value FROM meta WHERE meta_key='species.classification'"
            classifications = mysql_fetch_data(
                class_query,
                db,
                server_dict[use_server]["db_host"],
                server_dict[use_server]["db_port"],
                server_dict[use_server]["db_user"],
                server_dict[use_server]["db_pass"],
            )

            class_list = [c[0] for c in classifications]

            icon = "Metazoa.png"
            chordate = "Chordata" in class_list
            for classification in class_list:
                if classification in icon_dict:
                    icon = icon_dict[classification]
                    break

            if chordate and icon == "Metazoa.png":
                icon = "Chordates.png"

            write_yaml(info_dict, icon, yaml_out, project, use_server, alternate, guuid, ftp_client)

    ftp_client.close_connections()


if __name__ == "__main__":
    main()

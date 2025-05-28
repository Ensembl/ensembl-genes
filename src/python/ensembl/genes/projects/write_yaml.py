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
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Create a global session with a retry strategy
session = requests.Session()
retry_strategy = Retry(
    total=3,                # Number of retries
    backoff_factor=1,       # Delay between retries (1, 2, 4 seconds, etc.)
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["HEAD", "GET", "OPTIONS"]
)
adapter = HTTPAdapter(max_retries=retry_strategy)
session.mount("https://", adapter)
session.mount("http://", adapter)


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


def check_url_status(url, timeout=10):
    """
    Checks if a given URL is reachable with a longer timeout and retries.

    Args:
        url (str): The URL to check.
        timeout (int): Timeout in seconds for the request.

    Returns:
        bool: True if the URL returns a 200 status code, False otherwise.
    """
    try:
        response = session.head(url, allow_redirects=True, timeout=timeout)
        return response.status_code == 200
    except requests.RequestException as e:
        print(f"Error checking URL {url}: {e}")
        return False


def validate_ftp_link(
    url: str,
    link_type: str,
    scientific_name: Optional[str] = None,
    assembly_accession: Optional[str] = None
) -> str:
    """
    Validates the given FTP URL by checking if it exists.

    Returns:
      - the URL if it exists;
      - if the URL is missing and link_type == "ftp_dumps":
          URL for pre-release dumps (requires scientific_name & assembly_accession)
    """
    # if the URL is already valid, just return it
    if check_url_status(url):
        return url

    # otherwise, only handle the ftp_dumps case specially
    if link_type == "ftp_dumps":
        if not scientific_name or not assembly_accession:
            raise ValueError(
                "scientific_name and assembly_accession must be provided when link_type='ftp_dumps'"
            )
        ftp_folder = "https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release"
        return f"{ftp_folder}/{scientific_name}/{assembly_accession}"

    # for any other link_type
    return "Coming soon!"



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
            method = info_dict.get("genebuild.method_display", "BRAKER2")
            yaml += f"  annotation_method: {method}\n"

        annotation_gtf = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/{source}/geneset/{date}/genes.gtf.gz"
        yaml += f"  annotation_gtf: {validate_ftp_link(annotation_gtf, 'annotation_gtf')}\n"

        annotation_gff3 = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/{source}/geneset/{date}/genes.gff3"
        yaml += f"  annotation_gff3: {validate_ftp_link(annotation_gff3, 'annotation_gff3')}\n"

        proteins = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/{source}/geneset/{date}/pep.fa.gz"
        yaml += f"  proteins: {validate_ftp_link(proteins, 'proteins')}\n"

        transcripts = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/{source}/geneset/{date}/cdna.fa.gz"
        yaml += f"  transcripts: {validate_ftp_link(transcripts, 'transcripts')}\n"

        softmasked_genome = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/genome/softmasked.fa.gz"
        yaml += f"  softmasked_genome: {validate_ftp_link(softmasked_genome, 'softmasked_genome')}\n"

        rm_file = ftp_client.check_for_file(
            rm_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler"
        )
        if rm_file:
            yaml += f"  repeat_library: {validate_ftp_link(rm_file, 'repeat_library')}\n"

        ftp_dumps = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/"
        yaml += f"  ftp_dumps: {validate_ftp_link(ftp_dumps, 'ftp_dumps', species_name, info_dict['assembly.accession'])}\n"

        main_species_url = "http://www.ensembl.org/info/about/species.html"
        main_species_response = requests.get(main_species_url)
        if info_dict["assembly.accession"] in main_species_response.text:
            ensembl_link = f"https://www.ensembl.org/{species_name}/Info/Index"
            yaml += f"  ensembl_link: {validate_ftp_link(ensembl_link,'ensembl_link')}\n"
        else:
            beta_link = f"https://beta.ensembl.org/species/{guuid}"
            yaml += f"  beta_link: {validate_ftp_link(beta_link, 'beta_link')}\n"

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
                    yaml += f"  busco_score: {validate_ftp_link(busco_file, 'busco')}\n"
            if alternate:
                alternate_url = f"https://rapid.ensembl.org/{alternate}/Info/Index"
                yaml += f"  alternate: {validate_ftp_link(alternate_url, 'alternate')}\n"

    elif use_server == "main":
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
            yaml += "  annotation_method: Ensembl genebuild\n"

        source = 'ensembl'

        annotation_gtf = f"{ftp_base}/gtf/{info_dict['species.production_name']}/{uc_prod_name}.{assembly_name}.{release_number}.gtf.gz"
        yaml += f"  annotation_gtf: {validate_ftp_link(annotation_gtf, 'annotation_gtf')}\n"

        annotation_gff3 = f"{ftp_base}/gff3/{info_dict['species.production_name']}/{uc_prod_name}.{assembly_name}.{release_number}.gff3.gz"
        yaml += f"  annotation_gff3: {validate_ftp_link(annotation_gff3, 'annotation_gff3')}\n"

        proteins = f"{ftp_base}/fasta/{info_dict['species.production_name']}/pep/{uc_prod_name}.{assembly_name}.pep.all.fa.gz"
        yaml += f"  proteins: {validate_ftp_link(proteins, 'proteins')}\n"

        transcripts = f"{ftp_base}/fasta/{info_dict['species.production_name']}/cdna/{uc_prod_name}.{assembly_name}.cdna.all.fa.gz"
        yaml += f"  transcripts: {validate_ftp_link(transcripts, 'transcripts')}\n"

        softmasked_genome = f"{ftp_base}/fasta/{info_dict['species.production_name']}/dna/{uc_prod_name}.{assembly_name}.dna_sm.toplevel.fa.gz"
        yaml += f"  softmasked_genome: {validate_ftp_link(softmasked_genome, 'softmasked_genome')}\n"

        rm_file = ftp_client.check_for_file(
            rm_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler"
        )
        if rm_file:
            yaml += f"  repeat_library: {validate_ftp_link(rm_file, 'repeat_library')}\n"

        yaml += f"  ftp_dumps: {validate_ftp_link(ftp_base, 'ftp_dumps')}\n"
        if info_dict["species.production_name"] in prod_url_list:
            ensembl_link = f"https://www.ensembl.org/{uc_prod_name}/Info/Index"
            yaml += f"  ensembl_link: {validate_ftp_link(ensembl_link, 'ensembl_link')}\n"
        elif project == "geneswitch":
            ensembl_link = f"https://e102.ensembl.org/{species_name}/Info/Index"
            yaml += f"  ensembl_link: {validate_ftp_link(ensembl_link, 'ensembl_link')}\n"
        else:
            ensembl_link = f"https://www.ensembl.org/{species_name}/Info/Index"
            yaml += f"  ensembl_link: {validate_ftp_link(ensembl_link, 'ensembl_link')}\n"

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
                    yaml += f"  busco_score: {validate_ftp_link(busco_file, 'busco')}\n"
            if alternate:
                alternate_url = f"https://rapid.ensembl.org/{alternate}/Info/Index"
                yaml += f"  alternate: {validate_ftp_link(alternate_url, 'alternate')}\n"

    print(yaml, file=yaml_out)


def check_database_on_server(db, server_key, server_dict):
    """
    Checks if a database exists on a given server.
    """
    try:
        conn = pymysql.connect(
            host=server_dict[server_key]["db_host"],
            user=server_dict[server_key]["db_user"],
            passwd=server_dict[server_key]["db_pass"],
            port=server_dict[server_key]["db_port"],
        )
        with conn.cursor() as cur:
            cur.execute("SHOW DATABASES LIKE %s", (db,))
            return cur.fetchone() is not None
    except pymysql.MySQLError as e:
        print(f"Error connecting to {server_key}: {e}")
        return False
    finally:
        if 'conn' in locals() and conn:
            conn.close()


def find_database_server(db, server_dict):
    """
    Determines which server the given database exists on.
    """
    for server_key in ["st5", "st6", "main"]:
        if check_database_on_server(db, server_key, server_dict):
            return server_key
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

    with open(args.config_file) as config_f:
        server_dict = json.load(config_f)

    icon_dict = {}
    with open("icons.txt") as icon_file:
        for line in icon_file:
            fields = line.split()
            if len(fields) >= 2:
                icon_dict[fields[0]] = fields[1]

    with open(args.db_file) as db_file:
        db_list = db_file.read().strip().split("\n")

    sorted_db_list = sorted(db_list)

    if project == "aquafaang":
        for db in sorted_db_list:
            if "danio_rerio_core" in db:
                sorted_db_list.append(sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "bovreg":
        for db in sorted_db_list:
            if "bos_taurus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

    if project == "geneswitch":
        for db in sorted_db_list:
            if "sus_scrofa_core" in db or "gallus_gallus_core" in db:
                sorted_db_list.insert(0, sorted_db_list.pop(sorted_db_list.index(db)))

    ftp_client = EnsemblFTP()

    with open(f"{project}_species.yaml", "w") as yaml_out:
        for line in sorted_db_list:
            db = line.split('\t')[0].strip()
            guuid = line.split('\t')[1].strip()
            use_server = find_database_server(db, server_dict)

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

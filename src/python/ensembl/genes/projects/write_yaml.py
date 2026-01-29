#!/usr/bin/env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# pylint: disable=missing-module-docstring,consider-using-in,  too-many-locals, too-many-return-statements, missing-timeout, no-else-return, too-many-statements, too-many-arguments, too-many-branches, too-many-nested-blocks, logging-not-lazy, logging-fstring-interpolation, consider-using-dict-items, broad-exception-caught, line-too-long, unused-argument, missing-function-docstring, unspecified-encoding, broad-exception-raised
import argparse
import json
import re
import sys
import time
from ftplib import FTP, error_perm, error_temp
from typing import Any, Dict, Optional, Tuple
import socket
import pymysql
import requests


class EnsemblFTP:
    """
    Robust FTP client for Ensembl/EBI with reconnect + tolerant shutdown.
    """

    def __init__(
        self, timeout: int = 30, max_retries: int = 2, retry_sleep: float = 3.0
    ) -> None:
        """
        Initialize EnsemblFTP client and connect to Ensembl and EBI FTP servers.

        Args:
            timeout (int, optional): Timeout for FTP connections. Defaults to 30.
            max_retries (int, optional): Maximum number of retries for failed operations. Defaults to 2.
            retry_sleep (float, optional): Sleep duration between retries. Defaults to 3.0.
        """
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_sleep = retry_sleep

        self.ensembl_ftp: Optional[FTP] = None
        self.ebi_ftp: Optional[FTP] = None

        self._connect_ensembl()
        self._connect_ebi()

        self.ensembl_ftp_path = "https://ftp.ensembl.org/"
        self.ebi_ftp_path = "https://ftp.ebi.ac.uk/"

    # ---- connections ----
    def _connect_ensembl(self):
        """Connect to Ensembl FTP server."""
        try:
            self.ensembl_ftp = FTP("ftp.ensembl.org", timeout=self.timeout)
            self.ensembl_ftp.set_pasv(True)
            self.ensembl_ftp.login()
        except Exception:
            self.ensembl_ftp = None
            raise

    def _connect_ebi(self):
        """Connect to EBI FTP server."""
        try:
            self.ebi_ftp = FTP("ftp.ebi.ac.uk", timeout=self.timeout)
            self.ebi_ftp.set_pasv(True)
            self.ebi_ftp.login()
        except Exception:
            self.ebi_ftp = None
            raise

    def _retry(self, fn, which: str, *args, **kwargs) -> Any:
        """
        Retry FTP operation; on failure, reconnect that endpoint and try again.
        """
        last_exc: Optional[BaseException] = None
        for _ in range(self.max_retries):
            try:
                return fn(*args, **kwargs)
            except (
                ConnectionResetError,
                EOFError,
                OSError,
                error_temp,
                socket.timeout,
            ) as e:
                last_exc = e
                try:
                    if which == "ensembl":
                        self._connect_ensembl()
                    else:
                        self._connect_ebi()
                except Exception:
                    pass
                time.sleep(self.retry_sleep)
        assert last_exc is not None
        raise last_exc

    # ---- utils ----
    def _require_ftp(self, ftp: Optional[FTP]) -> FTP:
        if ftp is None:
            raise RuntimeError("FTP connection is not available")
        return ftp

    def return_to_root(self, ftp_connection: FTP) -> None:
        """Return to root

        Args:
            ftp_connection (FTP): Ftp path
        """
        which = "ensembl" if ftp_connection is self.ensembl_ftp else "ebi"
        self._retry(ftp_connection.cwd, which, "/")

    # ---- lookups ----
    def check_for_file(
        self,
        species_name: str,
        prod_name: str,
        accession: str,
        source: str,
        file_type: str,
    ) -> str:
        """
        Check for existence of specific file types on FTP servers.
        Args:
            species_name (str): _species_name_
            prod_name (str): _prod_name_
            accession (str): GCA accession
            source (str): _source_
            file_type (str): type of file to check for: "repeatmodeler" or "busco"

        Returns:
            str: URL to the file if found, else empty string
        """
        if file_type == "repeatmodeler":
            ftp_connection = self.ebi_ftp
            which = "ebi"
            ftp_path = self.ebi_ftp_path
            path = (
                "pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/"
                + species_name
                + "/"
            )
            file_name = accession + ".repeatmodeler.fa"
        elif file_type == "busco":
            ftp_connection = self.ensembl_ftp
            which = "ensembl"
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
            file_name_protein = prod_name + "_protein_busco_short_summary.txt"
            file_name_alternative = prod_name + "_busco_short_summary.txt"
        else:
            return ""

        try:
            ftp_connection = self._require_ftp(ftp_connection)
            self.return_to_root(ftp_connection)
            self._retry(ftp_connection.cwd, which, path)
            files_list = self._retry(ftp_connection.nlst, which)

            if file_type == "busco":
                if file_name_protein in files_list:
                    return ftp_path + path + file_name_protein
                if file_name_alternative in files_list:
                    return ftp_path + path + file_name_alternative
                return ""
            else:
                return ftp_path + path + file_name if file_name in files_list else ""
        except error_perm as e:
            if "550" in str(e):
                return ""
            print(f"FTP permission error: {e}", file=sys.stderr)
            return ""
        except error_temp as e:
            print(f"FTP temporary error: {e}", file=sys.stderr)
            return ""
        except Exception as e:
            print(f"Error while checking FTP file: {e}", file=sys.stderr)
            return ""

    def check_pre_release_file(
        self, species_name: str, accession: str, extension: str
    ) -> str:
        """
        Check for pre-release files on EBI FTP server.

        Args:
            species_name (str): _species_name_
            accession (str): GCA accession
            extension (str): File extension to look for (e.g., ".gtf.gz")

        Returns:
            str: URL to the pre-release file if found, else empty string
        """
        ftp = self.ebi_ftp
        which = "ebi"
        base = self.ebi_ftp_path
        path = f"pub/databases/ensembl/pre-release/{species_name}/{accession}/"
        try:

            ftp = self._require_ftp(ftp)
            self.return_to_root(ftp)
            self._retry(ftp.cwd, which, path)
            for fname in self._retry(ftp.nlst, which):
                if fname.lower().endswith(extension):
                    return base + path + fname
        except error_perm:
            return ""
        except Exception as e:
            print(f"Error while checking pre-release file: {e}", file=sys.stderr)
        return ""

    def close_connections(self) -> None:
        """
        Best-effort shutdown: if QUIT fails because the server dropped us,
        fall back to close() and never raise.
        """
        for conn in (self.ensembl_ftp, self.ebi_ftp):
            if not conn:
                continue
            try:
                conn.quit()
            except Exception:
                try:
                    conn.close()
                except Exception:
                    pass


def mysql_fetch_data(
    query: str, database: str, host: str, port: int, user: str, password: str
) -> Tuple:
    """
    Fetch data from a MySQL database using a given query.

    Args:
        query (str): SQL query to execute.
        database (str): Name of the database.
        host (str): Database host.
        port (int): Database port.
        user (str): Database user.
        password (str): Database password.

    Returns:
    Tuple containing the fetched rows from the database.
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
        response = requests.head(
            url, allow_redirects=True, timeout=5
        )  # Use HEAD for efficiency
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
    ftp_client: EnsemblFTP,
) -> None:
    """
    Write YAML content for a species entry based on the provided metadata and project type.

    Args:
        info_dict (Dict[str, str]): Metadata dictionary for the species.
        icon (str): Icon URL or identifier for the species.
        yaml_out: Output stream to write the YAML content.
        project (str): Project type (e.g., "vgp", "dtol", etc.).
        use_server (str): Server to use ("st5", "st6", "main").
        alternate (str): Alternate assembly name if applicable.
        guuid (str): Genome unique identifier.
        ftp_client (EnsemblFTP): FTP client instance for file checks.

    Returns:
        None
    """
    prod_url_list = ["bos_taurus_hybrid", "bos_indicus_hybrid"]

    assembly_name = info_dict["assembly.name"].replace(" ", "_")
    try:
        date = info_dict["genebuild.last_geneset_update"].replace("-", "_")
    except KeyError:
        date = ""
    species_name = info_dict["species.scientific_name"].replace(" ", "_")
    species_name = species_name.replace(".", "")
    lc_species_name = species_name.lower()

    if "species.strain" in info_dict and info_dict["species.strain"] != "reference":
        info_dict["species.scientific_name"] += f" ({info_dict['species.strain']})"

    uc_prod_name = info_dict["species.production_name"].capitalize()

    # Get submitter from assembly report
    assembly_report_url = (
        f"https://www.ncbi.nlm.nih.gov/assembly/{info_dict['assembly.accession']}"
    )
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

    if use_server == "st5" or use_server == "st6" or use_server == "gb1":
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
        gtf = (
            f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/"
            f"{source}/geneset/{date}/genes.gtf.gz"
        )
        if check_url_status(gtf):
            yaml += f"  annotation_gtf: {gtf}\n"
        else:
            fb = ftp_client.check_pre_release_file(
                species_name, info_dict["assembly.accession"], ".gtf.gz"
            )
            if fb:
                yaml += f"  annotation_gtf: {fb}\n"
        # ==== annotation_gff3: try gzipped, then uncompressed, then pre‑release ====
        gff_gz = (
            f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/"
            f"{source}/geneset/{date}/genes.gff3.gz"
        )
        gff = (
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
                species_name, info_dict["assembly.accession"], ".gff3.gz"
            )
            if not fb:
                fb = ftp_client.check_pre_release_file(
                    species_name, info_dict["assembly.accession"], ".gff3"
                )
            if fb:
                yaml += f"  annotation_gff3: {fb}\n"
        # ==== proteins ====
        pep = (
            f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/"
            f"{source}/geneset/{date}/pep.fa.gz"
        )
        if check_url_status(pep):
            yaml += f"  proteins: {pep}\n"
        # ==== transcripts ====
        cdna = (
            f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/"
            f"{source}/geneset/{date}/cdna.fa.gz"
        )
        if check_url_status(cdna):
            yaml += f"  transcripts: {cdna}\n"
        # ==== softmasked genome ====
        soft = (
            f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}/genome/"
            f"softmasked.fa.gz"
        )
        if check_url_status(soft):
            yaml += f"  softmasked_genome: {soft}\n"
        else:
            fb = ftp_client.check_pre_release_file(
                species_name, info_dict["assembly.accession"], ".dna.softmasked.fa.gz"
            )
            if fb:
                yaml += f"  softmasked_genome: {fb}\n"
        rm_file = ftp_client.check_for_file(
            lc_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler",
        )
        if rm_file:
            yaml += f"  repeat_library: {rm_file}\n"
        # ==== ftp dumps ====
        dumps = f"{ftp_base}/{species_name}/{info_dict['assembly.accession']}"
        if check_url_status(dumps):
            yaml += f"  ftp_dumps: {ftp_base}/{species_name}/{info_dict['assembly.accession']}/\n"
        else:
            yaml += f"  ftp_dumps: https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/{species_name}/{info_dict['assembly.accession']}/\n"

        main_species_url = "http://www.ensembl.org/info/about/species.html"
        main_species_response = requests.get(main_species_url)
        if info_dict["assembly.accession"] in main_species_response.text:
            yaml += (
                f"  ensembl_link: https://www.ensembl.org/{species_name}/Info/Index\n"
            )
        else:
            beta_link = f"https://beta.ensembl.org/species/{guuid}"
            if check_url_status(beta_link):
                yaml += f"  beta_link: https://beta.ensembl.org/species/{guuid}\n"
            else:
                yaml += "  beta_link: Coming soon!\n"

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            try:
                yaml += f"  busco_score: {info_dict['genebuild.busco']}\n"
            except KeyError:
                busco_file = ftp_client.check_for_file(
                    species_name,
                    info_dict["species.production_name"],
                    info_dict["assembly.accession"],
                    source,
                    "busco",
                )
                if busco_file:
                    yaml += f"  busco_score: {busco_file}\n"
            try:
                yaml += f"  busco_lineage: {info_dict['genebuild.busco_dataset']}\n"
            except KeyError:
                if busco_file:
                    yaml += f"  busco_lineage: {busco_file}\n"
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

        source = "ensembl"  # for consistency with busco retrieval logic

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
            lc_species_name,
            info_dict["species.production_name"],
            info_dict["assembly.accession"],
            source,
            "repeatmodeler",
        )
        if rm_file:
            yaml += f"  repeat_library: {rm_file}\n"

        yaml += f"  ftp_dumps: {ftp_base}\n"
        if info_dict["species.production_name"] in prod_url_list:
            yaml += (
                f"  ensembl_link: https://www.ensembl.org/{uc_prod_name}/Info/Index\n"
            )
        elif project == "geneswitch":
            # geneswitch project frozen at e102
            yaml += (
                f"  ensembl_link: https://e102.ensembl.org/{species_name}/Info/Index\n"
            )
        else:
            yaml += (
                f"  ensembl_link: https://www.ensembl.org/{species_name}/Info/Index\n"
            )

        if project in ("dtol", "erga", "cbp", "bge", "asg"):
            try:
                yaml += f"  busco_score: {info_dict['genebuild.busco']}\n"
            except KeyError:
                busco_file = ftp_client.check_for_file(
                    lc_species_name,
                    info_dict["species.production_name"],
                    info_dict["assembly.accession"],
                    source,
                    "busco",
                )
                if busco_file:
                    yaml += f"  busco_score: {busco_file}\n"
            try:
                yaml += f"  busco_lineage: {info_dict['genebuild.busco_dataset']}\n"
            except KeyError:
                if busco_file:
                    yaml += f"  busco_lineage: {busco_file}\n"
            if alternate:
                alternate_url = f"https://rapid.ensembl.org/{alternate}/Info/Index"
                yaml += f"  alternate: {alternate_url}\n"

    print(yaml, file=yaml_out)


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
            # print(f"Checking for DB '{db}' on server '{server_key}'", file=sys.stderr)
            cur.execute(
                "SELECT SCHEMA_NAME FROM information_schema.SCHEMATA WHERE SCHEMA_NAME = %s",
                (db,),
            )
            result = cur.fetchone()
            # print(f"Query result on {server_key}: {result}", file=sys.stderr)
            return result is not None

    except pymysql.MySQLError as e:
        print(f"Error connecting to {server_key}: {e}")
        return False
    finally:
        if "conn" in locals() and conn:
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
        choices=[
            "aquafaang",
            "asg",
            "bge",
            "bovreg",
            "cbp",
            "dtol",
            "erga",
            "geneswitch",
            "vgp",
        ],
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
    ftp_client = EnsemblFTP(timeout=30, max_retries=2, retry_sleep=3.0)

    # Open the output YAML file
    with open(f"{project}_species.yaml", "w") as yaml_out:
        # will replace with guiid eventually
        for line in sorted_db_list:
            db = line.split("\t")[0].strip()
            guuid = line.split("\t")[1].strip()
            if guuid == "unknown":
                use_server = "gb1"
            else:
                use_server = find_database_server(db, server_dict)

            # Retrieve metadata from the chosen server
            # update this with query from the metadata db
            info_query = (
                "SELECT meta_key, meta_value FROM meta "
                "WHERE meta_key IN "
                "('species.scientific_name','assembly.accession','assembly.name',"
                "'species.production_name','species.strain','schema_version',"
                "'genebuild.last_geneset_update','species.annotation_source','genebuild.busco','genebuild.busco_dataset') "
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
            alternate_assembly_name = (
                info_dict["assembly.name"] + "_alternate_haplotype"
            )
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

            class_list = [c[0] for c in classifications]

            icon = "Metazoa.png"
            chordate = "Chordata" in class_list
            for classification in class_list:
                if classification in icon_dict:
                    icon = icon_dict[classification]
                    break

            if chordate and icon == "Metazoa.png":
                icon = "Chordates.png"

            write_yaml(
                info_dict,
                icon,
                yaml_out,
                project,
                use_server,
                alternate,
                guuid,
                ftp_client,
            )

    ftp_client.close_connections()


if __name__ == "__main__":
    main()

"""FTP/HTTP helpers for checking Ensembl/EBI file and beta-page availability."""

import socket
import sys
import time
from ftplib import FTP, error_perm, error_temp
from typing import Any, Optional

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
            max_retries (int, optional): Maximum number of retries for failed
                operations. Defaults to 2.
            retry_sleep (float, optional): Sleep duration between retries.
                Defaults to 3.0.
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
                except Exception:  # pylint: disable=broad-exception-caught
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
    def check_for_file(  # pylint: disable=too-many-locals, too-many-arguments, too-many-return-statements
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
        # Initialize all file names upfront to avoid possibly-used-before-assignment
        file_name_protein: Optional[str] = None
        file_name_alternative: Optional[str] = None
        file_name: Optional[str] = None

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
                if file_name_protein and file_name_protein in files_list:
                    return ftp_path + path + file_name_protein
                if file_name_alternative and file_name_alternative in files_list:
                    return ftp_path + path + file_name_alternative
                return ""
            return (
                ftp_path + path + file_name
                if file_name and file_name in files_list
                else ""
            )
        except error_perm as e:
            if "550" in str(e):
                return ""
            print(f"FTP permission error: {e}", file=sys.stderr)
            return ""
        except error_temp as e:
            print(f"FTP temporary error: {e}", file=sys.stderr)
            return ""
        except Exception as e:  # pylint: disable=broad-exception-caught
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
        except Exception as e:  # pylint: disable=broad-exception-caught
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
            except Exception:  # pylint: disable=broad-exception-caught
                try:
                    conn.close()
                except Exception:  # pylint: disable=broad-exception-caught
                    pass


def check_url_status(url: str) -> bool:
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
        if response.status_code == 200:
            return True
        if response.status_code in [404, 403]:
            return False
        # Fallback to GET for other HTTP errors (e.g. 405 Method Not Allowed)
        response = requests.get(url, allow_redirects=True, stream=True, timeout=5)
        response.close()
        return response.status_code == 200
    except requests.RequestException as e:
        print(f"Error checking URL {url}: {e}")
        return False


# Markers in the beta.ensembl.org page body that indicate the species page
# is not actually available (the error page still returns HTTP 200).
_BETA_UNAVAILABLE_MARKERS = (
    "We do not recognise the species identified",
    "Find available species in the Species selector",
)


def check_beta_species_status(genome_uuid: str) -> str:
    """
    Check whether a beta.ensembl.org species page is actually usable.

    The beta site returns HTTP 200 even for unrecognised species (showing an
    error page), so a status-code check alone is insufficient — the response
    body must be inspected for known "not available" markers.

    Args:
        genome_uuid (str): The genome UUID to check.

    Returns:
        str: One of:
            - "available"   : a normal, usable beta species page
            - "unavailable" : page loads (or non-200) but species not recognised
            - "error"       : network/HTTP error, or missing/unknown UUID
    """
    if not genome_uuid or genome_uuid == "unknown":
        return "error"
    url = f"https://beta.ensembl.org/species/{genome_uuid}"
    try:
        response = requests.get(url, allow_redirects=True, timeout=10)
        # Non-200 responses are treated as unavailable (page not served).
        if response.status_code != 200:
            return "unavailable"
        body = response.text or ""
        if any(marker in body for marker in _BETA_UNAVAILABLE_MARKERS):
            return "unavailable"
        return "available"
    except requests.RequestException as e:
        print(f"Error checking beta species {genome_uuid}: {e}", file=sys.stderr)
        return "error"

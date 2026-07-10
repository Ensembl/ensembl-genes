"""HTTP client for discovering WormBase ParaSite files from EBI FTP."""

import re
import sys
from typing import List

import requests


class WormBaseHTTP:
    """HTTP client to discover WormBase ParaSite data over EBI's HTTP-exposed FTP."""

    def __init__(self, release: str):
        self.release = release
        self.base_url = f"https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/{release}/species/"

    def _get_links(self, url: str) -> List[str]:
        """Fetch a directory listing and extract href links."""
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Match hrefs, ignoring parent directory links
            links = re.findall(r'href="([^"/]+)/?"', response.text)
            return [link for link in set(links) if link and not link.startswith("?") and link != "Parent Directory"]
        except requests.RequestException as e:
            print(f"Error fetching directory {url}: {e}", file=sys.stderr)
            return []

    def get_species(self) -> List[str]:
        """List all species slugs for the release."""
        return self._get_links(self.base_url)

    def get_bioprojects(self, species_slug: str) -> List[str]:
        """List all BioProjects for a given species slug."""
        url = f"{self.base_url}{species_slug}/"
        return [link for link in self._get_links(url) if link.startswith("PRJ")]

    def get_files(self, species_slug: str, bioproject: str) -> List[str]:
        """List all files for a given species and BioProject."""
        url = f"{self.base_url}{species_slug}/{bioproject}/"
        return self._get_links(url)

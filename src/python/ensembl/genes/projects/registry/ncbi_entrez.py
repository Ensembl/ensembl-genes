"""
Data fetcher for NCBI via Entrez HTTP endpoints.
Performs web scraping for submitters and HPRC population data.
"""

import logging
import re
from typing import Optional

import requests

from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.config import ProjectConfig

logger = logging.getLogger(__name__)


def _fetch_assembly_report_type(accession: str, assembly_name: str) -> Optional[str]:
    """Downloads the authoritative assembly report from NCBI genomes FTP to parse maternal/paternal."""
    if not accession or not assembly_name:
        return None
    parts = accession.split("_")
    if len(parts) < 2:
        return None
    num = parts[1].split(".")[0]
    if len(num) < 9:
        return None

    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{parts[0]}/{num[0:3]}/{num[3:6]}/{num[6:9]}/{accession}_{assembly_name.replace(' ', '_')}/{accession}_{assembly_name.replace(' ', '_')}_assembly_report.txt"
    try:
        res = requests.get(url, timeout=10)
        if res.status_code == 200:
            for line in res.text.splitlines():
                if line.startswith("# Assembly type:"):
                    lower = line.lower()
                    if "maternal" in lower:
                        return "maternal"
                    if "paternal" in lower:
                        return "paternal"
    except Exception as e:
        logger.warning(f"Failed to fetch assembly report for {accession}: {e}")
    return None


def patch_ncbi_data(meta: GenomeMetadata, config: ProjectConfig) -> None:
    """Modifies the GenomeMetadata inline with NCBI scraped data."""
    if not (config.scrape_ncbi_submitter or config.scrape_ncbi_population):
        return

    # First, always get the assembly page to find submitters and the biosample ID
    assembly_url = f"https://www.ncbi.nlm.nih.gov/assembly/{meta.accession}"
    try:
        assembly_html = requests.get(assembly_url, timeout=10).text.splitlines()
    except requests.RequestException as e:
        logger.warning(f"Failed to fetch assembly page for {meta.accession}: {e}")
        return

    biosample_id = None

    for line in assembly_html:
        # Find Submitter for standard projects
        if config.scrape_ncbi_submitter and not meta.assembly_submitter:
            # Parse the submitter from the NCBI assembly page HTML
            sub_regex = re.search(r"Submitter<\/dt><dd>([^<]+)<\/dd>", line)
            if sub_regex:
                meta.assembly_submitter = sub_regex.group(1).title()

        # Find BioSample link
        if config.scrape_ncbi_population and not biosample_id:
            bio_regex = re.search(r"\/biosample\/([A-Z0-9]+)\/\"", line)
            if bio_regex:
                biosample_id = bio_regex.group(1)

    # Find parent_of_origin from authoritative FTP Assembly Report semantics
    if not meta.parent_of_origin:
        meta.parent_of_origin = _fetch_assembly_report_type(
            meta.accession, meta.assembly_name
        )

    # Fallback heuristic for parent_of_origin if authoritative metadata is missing
    if not meta.parent_of_origin and meta.assembly_name:
        asm_lower = meta.assembly_name.lower()
        if "_mat" in asm_lower or "maternal" in asm_lower:
            meta.parent_of_origin = "maternal"
            logger.info(
                f"Fallback heuristic used for {meta.accession}: inferred maternal parent_of_origin from assembly name."
            )
        elif "_pat" in asm_lower or "paternal" in asm_lower:
            meta.parent_of_origin = "paternal"
            logger.info(
                f"Fallback heuristic used for {meta.accession}: inferred paternal parent_of_origin from assembly name."
            )

    # Hardcoded override for HPRC project
    if config.schema_type == "hprc":
        if meta.accession == "GCA_009914755.4":
            meta.assembly_submitter = "T2T Consortium"
        else:
            meta.assembly_submitter = "UCSC Genomics Institute"

    # Scrape BioSample page if needed
    if config.scrape_ncbi_population and biosample_id:
        _scrape_biosample_population(meta, biosample_id)


def _scrape_biosample_population(meta: GenomeMetadata, biosample_id: str) -> None:
    """Replicates the HPRC BioSample population scraping logic."""
    biosample_url = f"https://www.ncbi.nlm.nih.gov/biosample/{biosample_id}"
    try:
        biosample_html = requests.get(biosample_url, timeout=10).text.splitlines()
    except requests.RequestException as e:
        logger.warning(f"Failed to fetch biosample page {biosample_id}: {e}")
        return

    pop_string = ""
    for line in biosample_html:
        pop_desc_regex = re.search(
            r"Population Description<\/th><td>([A-Za-z0-9 ]+)", line
        )
        pop_regex = re.search(r"population=([A-Za-z0-9 ,]+)", line)
        race_regex = re.search(r"race<\/th><td>([A-Za-z0-9 ]+)", line)

        if pop_desc_regex:
            pop_string += pop_desc_regex.group(1)
        elif pop_regex:
            pop_string += pop_regex.group(1)
        elif race_regex:
            pop_string += race_regex.group(1)

    if pop_string:
        pop_string = pop_string.title().replace("Usa", "USA")
        meta.population = pop_string

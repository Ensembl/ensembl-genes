"""
Project configuration mappings determining data fetch behavior and YAML schema rules.
"""

from dataclasses import dataclass
from typing import List, Optional


@dataclass
class ProjectConfig:  # pylint: disable=too-many-instance-attributes
    """
    Directs the pipeline on what data sources are allowed and how to map
    GenomeMetadata to the final YAML schema format.
    """

    project_name: str

    # Data Fetching Directives
    allow_core_db_fallback: bool = True
    allow_beta_urls: bool = True
    scrape_ncbi_submitter: bool = False
    scrape_ncbi_population: bool = False
    check_ftp_repeats: bool = False
    check_ftp_variants: bool = False

    # Formatting Directives
    schema_type: str = "standard"  # "standard", "hprc", "mouse"

    # Pre-release Registry Scoping
    bioproject_scoping: Optional[List[str]] = None
    custom_group_scoping: Optional[List[str]] = None

    # Base URL prefixes (if customized for project)
    ftp_base_url: str = "https://ftp.ensembl.org/pub"
    rapid_base_url: str = "https://rapid.ensembl.org"


def get_project_config(project_name: str) -> ProjectConfig:
    """Factory to retrieve configuration for a specific project."""
    name_lower = project_name.lower()

    if name_lower == "hprc":
        return ProjectConfig(
            project_name="hprc",
            allow_core_db_fallback=False,
            allow_beta_urls=False,
            scrape_ncbi_submitter=True,
            scrape_ncbi_population=True,
            check_ftp_variants=True,
            schema_type="hprc",
            bioproject_scoping=["HPRC"],
        )
    if name_lower in ("mouse_genomes", "mouse"):
        return ProjectConfig(
            project_name="mouse_genomes",
            schema_type="mouse",
            allow_beta_urls=True,
            # mouse genomes generally don't have bioproject scoping in gb_schema
        )
    if name_lower == "aquafaang":
        return ProjectConfig(
            project_name=name_lower,
            schema_type="standard",
            custom_group_scoping=["AQUA-FAANG"],
        )
    if name_lower in ["vgp", "darwin_tree_of_life", "erga", "cbp", "bge", "asg"]:
        bp_scope = (
            ["DToL"] if name_lower == "darwin_tree_of_life" else [name_lower.upper()]
        )
        return ProjectConfig(
            project_name=name_lower,
            schema_type="standard",
            check_ftp_repeats=(name_lower in ["vgp", "darwin_tree_of_life", "erga"]),
            scrape_ncbi_submitter=(
                name_lower not in ["vgp", "darwin_tree_of_life", "erga"]
            ),
            bioproject_scoping=bp_scope,
        )
    # Default behavior for unknown standard projects
    return ProjectConfig(
        project_name=name_lower,
        schema_type="standard",
        scrape_ncbi_submitter=True,
        bioproject_scoping=[name_lower.upper()],
    )

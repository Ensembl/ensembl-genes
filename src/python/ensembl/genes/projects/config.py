"""
Project configuration mappings determining data fetch behavior and YAML schema rules.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Callable
from ensembl.genes.projects.models import GenomeMetadata

@dataclass
class ProjectConfig:
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
    
    # Base URL prefixes (if customized for project)
    ftp_base_url: str = "https://ftp.ensembl.org/pub"
    rapid_base_url: str = "https://rapid.ensembl.org"
    
def get_project_config(project_name: str) -> ProjectConfig:
    """Factory to retrieve configuration for a specific project."""
    name_lower = project_name.lower()
    
    if name_lower == "hprc":
        return ProjectConfig(
            project_name="hprc",
            allow_core_db_fallback=False,  # Enforce metadata DB usage where possible
            allow_beta_urls=False,
            scrape_ncbi_submitter=True,   # HPRC currently needs NCBI biosample scraping (Phase 6 target to fix)
            scrape_ncbi_population=True,
            check_ftp_variants=True,
            schema_type="hprc"
        )
    elif name_lower == "mouse_genomes" or name_lower == "mouse":
        return ProjectConfig(
            project_name="mouse_genomes",
            schema_type="mouse",
            allow_beta_urls=True
        )
    elif name_lower in ["vgp", "darwin_tree_of_life", "erga"]:
        # Standard projects with optional icons
        return ProjectConfig(
            project_name=name_lower,
            schema_type="standard",
            check_ftp_repeats=True
        )
    else:
        # Default behavior for unknown standard projects
        return ProjectConfig(
            project_name=name_lower,
            schema_type="standard",
            scrape_ncbi_submitter=True
        )

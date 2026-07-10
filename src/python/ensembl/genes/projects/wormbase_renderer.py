"""Renderer for WormBase ParaSite project YAML."""

from typing import Dict, List, Tuple

from ensembl.genes.projects.ftp_client import check_url_status


class WormBaseRenderer:
    """Processes discovered WormBase files and generates project YAML rows."""

    def __init__(self, release: str):
        self.release = release
        self.base_url = f"https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/{release}/species/"

    def _derive_species_name(self, species_slug: str) -> str:
        """Convert 'acanthocheilonema_viteae' to 'Acanthocheilonema viteae'."""
        parts = species_slug.split("_")
        if not parts:
            return species_slug
        parts[0] = parts[0].capitalize()
        return " ".join(parts)

    def render_row(
        self, species_slug: str, bioproject: str, files: List[str]
    ) -> Tuple[Dict, List[str]]:
        """
        Render a YAML row for a species and bioproject.
        Returns the row dictionary and a list of missing expected files.
        """
        prefix = f"{species_slug}.{bioproject}.{self.release}"
        dir_url = f"{self.base_url}{species_slug}/{bioproject}/"

        file_mappings = {
            "genome": f"{prefix}.genomic.fa.gz",
            "genome_masked": f"{prefix}.genomic_masked.fa.gz",
            "genome_softmasked": f"{prefix}.genomic_softmasked.fa.gz",
            "annotation_gff3": f"{prefix}.annotations.gff3.gz",
            "annotation_gtf": f"{prefix}.canonical_geneset.gtf.gz",
            "proteins": f"{prefix}.protein.fa.gz",
            "transcripts_mrna": f"{prefix}.mRNA_transcripts.fa.gz",
            "transcripts_cds": f"{prefix}.CDS_transcripts.fa.gz",
            "orthologues": f"{prefix}.orthologs.tsv.gz",
            "paralogues": f"{prefix}.paralogs.tsv.gz",
        }

        row = {
            "species": self._derive_species_name(species_slug),
            "bioproject": bioproject,
            "bioproject_link": f"https://www.ncbi.nlm.nih.gov/bioproject/{bioproject}",
            "release": self.release,
        }

        missing = []
        for key, expected_filename in file_mappings.items():
            if expected_filename in files:
                row[key] = f"{dir_url}{expected_filename}"
            else:
                missing.append(expected_filename)

        # Check for species page link
        species_page_url = f"https://parasite.wormbase.org/{species_slug}_{bioproject.lower()}"
        if check_url_status(species_page_url):
            row["species_page"] = species_page_url

        row["ftp_dumps"] = dir_url

        return row, missing

    def is_valid_row(self, row: Dict) -> bool:
        """A row is valid if it has at least one core file."""
        core_fields = [
            "genome",
            "genome_masked",
            "genome_softmasked",
            "annotation_gff3",
            "annotation_gtf",
            "proteins",
            "transcripts_mrna",
            "transcripts_cds"
        ]
        return any(field in row for field in core_fields)

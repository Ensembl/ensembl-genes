"""
Converts internal GenomeMetadata objects into specific project YAML schemas.
"""
import copy
from typing import Dict, Any, Optional
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.config import ProjectConfig

class YamlRenderer:
    """Renders GenomeMetadata into validated dictionary structures for YAML output."""
    
    def __init__(self, config: ProjectConfig):
        self.config = config
        self.icons = self._load_icons()

    def _load_icons(self) -> Dict[str, str]:
        import os
        icon_dict = {}
        icon_path = os.path.join(os.path.dirname(__file__), "icons.txt")
        if os.path.exists(icon_path):
            with open(icon_path) as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        icon_dict[parts[0]] = parts[1]
        return icon_dict

    def render(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Dispatches to the correct schema renderer based on project config."""
        if self.config.schema_type == "hprc":
            return self._render_hprc(meta)
        elif self.config.schema_type == "mouse":
            return self._render_mouse(meta)
        else:
            return self._render_standard(meta)

    def _render_standard(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema A: Standard Projects (VGP, DToL, ERGA)"""
        doc: Dict[str, Any] = {}
        
        # Core fields
        species_entry = meta.species_name
        if meta.strain:
            species_entry += f" {meta.strain}"
        doc["species"] = species_entry
        
        # Icon mapping (default Metazoa)
        if self.config.project_name in ["vgp", "dtol", "erga", "darwin_tree_of_life"]:
            doc["image"] = "Metazoa.png" # Need taxonomy lookup for deeper mapping
        elif self.config.scrape_ncbi_submitter and meta.assembly_submitter:
            doc["submitted_by"] = meta.assembly_submitter
            
        doc["accession"] = meta.accession
        
        # Annotations
        doc["annotation_method"] = meta.annotation_method or "BRAKER2"
            
        doc["annotation_gtf"] = self._build_ftp_url(meta, "geneset", "genes.gtf.gz")
        doc["annotation_gff3"] = self._build_ftp_url(meta, "geneset", "genes.gff3.gz")
        doc["proteins"] = self._build_ftp_url(meta, "geneset", "pep.fa.gz")
        doc["transcripts"] = self._build_ftp_url(meta, "geneset", "cdna.fa.gz")
        doc["softmasked_genome"] = self._build_ftp_url(meta, "genome", "softmasked.fa.gz")
        doc["ftp_dumps"] = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{meta.species_name}/{meta.accession}/"
        
        if self.config.check_ftp_repeats:
            pass # In reality, we'd check FTP existence
            
        # Linkages
        if meta.is_on_rapid:
            doc["ensembl_link"] = f"https://rapid.ensembl.org/{meta.species_name}"
        elif self.config.allow_beta_urls:
            doc["beta_link"] = f"https://beta.ensembl.org/species/{meta.genome_uuid}"
            
        # Quality
        if meta.busco_score:
            doc["busco_score"] = meta.busco_score
        if meta.busco_lineage:
            doc["busco_lineage"] = meta.busco_lineage
            
        # Relationships
        if meta.alternate_of:
            doc["alternate"] = meta.alternate_of

        return {k: v for k, v in doc.items() if v is not None}

    def _render_hprc(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema B: HPRC"""
        doc: Dict[str, Any] = {}
        
        doc["assembly"] = meta.assembly_name
        if meta.population:
            doc["population"] = meta.population
            
        if meta.parent_of_origin:
            doc["parent_of_origin"] = meta.parent_of_origin
            
        doc["assembly_accession"] = meta.accession
        doc["assembly_link"] = f"https://www.ebi.ac.uk/ena/browser/view/{meta.accession}"
        
        if meta.assembly_submitter:
            doc["assembly_submitter"] = meta.assembly_submitter
            
        doc["annotation_gtf"] = self._build_rapid_ftp_url(meta, "gtf")
        doc["annotation_gff3"] = self._build_rapid_ftp_url(meta, "gff3")
        doc["proteins"] = self._build_rapid_ftp_url(meta, "pep")
        doc["transcripts"] = self._build_rapid_ftp_url(meta, "cdna")
        doc["ftp_dumps"] = self._build_rapid_ftp_url(meta, "base")
        
        doc["rapid_link"] = f"{self.config.rapid_base_url}/{meta.species_name}"
        
        if self.config.check_ftp_variants:
            if meta.has_variants_clinvar:
                doc["variants_clinvar"] = self._build_rapid_ftp_url(meta, "clinvar")
            if meta.has_variants_gnomad:
                doc["variants_gnomad"] = self._build_rapid_ftp_url(meta, "gnomad")
            if meta.has_variants_vep:
                doc["variants_vep"] = self._build_rapid_ftp_url(meta, "vep")

        return {k: v for k, v in doc.items() if v is not None}

    def _render_mouse(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema C: Mouse Genomes"""
        doc: Dict[str, Any] = {}
        
        doc["species"] = meta.species_name
        if meta.strain:
            doc["strain"] = meta.strain
            
        doc["accession"] = meta.accession
        
        doc["annotation_gtf"] = self._build_ftp_url(meta, "geneset", "genes.gtf.gz")
        doc["annotation_gff3"] = self._build_ftp_url(meta, "geneset", "genes.gff3.gz")
        doc["proteins"] = self._build_ftp_url(meta, "geneset", "pep.fa.gz")
        doc["transcripts"] = self._build_ftp_url(meta, "geneset", "cdna.fa.gz")
        doc["softmasked_genome"] = self._build_ftp_url(meta, "genome", "softmasked.fa.gz")
        doc["ftp_dumps"] = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{meta.species_name}/{meta.accession}/"
        
        if self.config.allow_beta_urls:
            doc["beta_link"] = f"https://beta.ensembl.org/species/{meta.genome_uuid}"
            
        return {k: v for k, v in doc.items() if v is not None}

    def _build_ftp_url(self, meta: GenomeMetadata, category: str, file_suffix: str) -> str:
        """Helper to build standard FTP URLs for standard/mouse schemas."""
        date = meta.annotation_date or "2024_01"
        source = "ensembl"
        return f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{meta.species_name}/{meta.accession}/{source}/{category}/{date}/{file_suffix}"
        
    def _build_rapid_ftp_url(self, meta: GenomeMetadata, resource_type: str) -> str:
        """Helper to build Rapid Release FTP URLs"""
        base = f"https://ftp.ensembl.org/pub/rapid-release/species/{meta.species_name}"
        return f"{base}/{resource_type}"

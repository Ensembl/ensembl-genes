"""
Converts internal GenomeMetadata objects into specific project YAML schemas.
"""
import copy
from typing import Dict, Any, Optional
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.config import ProjectConfig
from ensembl.genes.projects.write_yaml import check_url_status
import requests
import xmltodict
import logging

logger = logging.getLogger(__name__)

class YamlRenderer:
    """Renders GenomeMetadata into validated dictionary structures for YAML output."""
    
    def __init__(self, config: ProjectConfig, ftp_client=None):
        self.config = config
        self.ftp_client = ftp_client
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

    def _fetch_taxonomy_classes(self, taxon_id: int) -> list:
        """Fetches lineage from NCBI to replicate core DB species.classification."""
        if not taxon_id:
            return []
        try:
            res = requests.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxon_id}&retmode=xml", timeout=5)
            if res.status_code == 200:
                data = xmltodict.parse(res.text)
                taxa_set = data.get('TaxaSet', {}).get('Taxon', {})
                if isinstance(taxa_set, list): taxa_set = taxa_set[0]
                lineage = taxa_set.get('LineageEx', {}).get('Taxon', [])
                if isinstance(lineage, dict): lineage = [lineage]
                return [t.get('ScientificName') for t in lineage if 'ScientificName' in t]
        except Exception as e:
            logger.warning(f"Failed to fetch taxonomy for {taxon_id}: {e}")
        return []

    def _render_standard(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema A: Standard Projects (VGP, DToL, ERGA)"""
        doc: Dict[str, Any] = {}
        
        # Core fields
        species_entry = meta.species_name
        if meta.strain:
            species_entry += f" {meta.strain}"
        doc["species"] = species_entry
        
        # Icon mapping
        if self.config.project_name in ["vgp", "dtol", "erga", "darwin_tree_of_life", "cbp", "bge", "asg"]:
            icon = "Metazoa.png"
            class_list = self._fetch_taxonomy_classes(meta.taxon_id)
            for classification in class_list:
                if classification in self.icons:
                    icon = self.icons[classification]
            doc["image"] = icon
        elif self.config.scrape_ncbi_submitter and meta.assembly_submitter:
            doc["submitted_by"] = meta.assembly_submitter
            
        doc["accession"] = meta.accession
        doc["annotation_method"] = meta.annotation_method or "BRAKER2"
        
        ftp_species_name = meta.species_name.capitalize().replace(" ", "_")
        
        # Determine tracking paths and test fallbacks (Fix 1 and 2)
        target_Released = meta.is_released
        
        # If released, we build normal FTPs. If they 404, we degrade to pre_release if possible
        if target_Released:
            released_gtf = self._build_ftp_url(meta, "geneset", "genes.gtf.gz", ftp_species_name)
            if not check_url_status(released_gtf):
                # Degrade to Pre-release
                target_Released = False
                
        # If still released (or not degraded)
        if target_Released:
            doc["annotation_gtf"] = self._build_ftp_url(meta, "geneset", "genes.gtf.gz", ftp_species_name)
            doc["annotation_gff3"] = self._build_ftp_url(meta, "geneset", "genes.gff3.gz", ftp_species_name)
            doc["proteins"] = self._build_ftp_url(meta, "geneset", "pep.fa.gz", ftp_species_name)
            doc["transcripts"] = self._build_ftp_url(meta, "geneset", "cdna.fa.gz", ftp_species_name)
            doc["softmasked_genome"] = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/genome/softmasked.fa.gz"
            
            # Additional pre-release specific check on GFF uncompressed legacy rule
            if not check_url_status(doc["annotation_gff3"]):
                uncompressed_gff = self._build_ftp_url(meta, "geneset", "genes.gff3", ftp_species_name)
                if check_url_status(uncompressed_gff):
                    doc["annotation_gff3"] = uncompressed_gff
        else:
            # Pre-release logic
            fb_gtf = self.ftp_client.check_pre_release_file(ftp_species_name, meta.accession, ".gtf.gz") if self.ftp_client else ""
            if not fb_gtf:
                return {} # Suppress entirely if neither exists!
                
            doc["annotation_gtf"] = fb_gtf
            
            fb_gff = self.ftp_client.check_pre_release_file(ftp_species_name, meta.accession, ".gff3.gz") if self.ftp_client else ""
            if not fb_gff:
                fb_gff = self.ftp_client.check_pre_release_file(ftp_species_name, meta.accession, ".gff3") if self.ftp_client else ""
            if fb_gff:
                doc["annotation_gff3"] = fb_gff
                
            doc["proteins"] = f"https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/{ftp_species_name}/{meta.accession}/{meta.accession}.pep.fa.gz"
            doc["transcripts"] = f"https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/{ftp_species_name}/{meta.accession}/{meta.accession}.cdna.fa.gz"
            
            fb_soft = self.ftp_client.check_pre_release_file(ftp_species_name, meta.accession, ".dna.softmasked.fa.gz") if self.ftp_client else ""
            if fb_soft:
                doc["softmasked_genome"] = fb_soft
        
        if target_Released:
            doc["ftp_dumps"] = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/"
        else:
            doc["ftp_dumps"] = f"https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/{ftp_species_name}/{meta.accession}/"
            
        # Linkages
        if meta.is_on_rapid:
            doc["ensembl_link"] = f"https://rapid.ensembl.org/{meta.species_name}"
        elif self.config.allow_beta_urls:
            if target_Released:
                doc["beta_link"] = f"https://beta.ensembl.org/species/{meta.genome_uuid}"
            else:
                doc["beta_link"] = "Coming soon!"
            
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
            
        doc["assembly_accession"] = meta.accession
        
        if meta.assembly_submitter:
            doc["assembly_submitter"] = meta.assembly_submitter
            
        doc["annotation_gtf"] = self._build_rapid_ftp_url(meta, "gtf")
        doc["annotation_gff3"] = self._build_rapid_ftp_url(meta, "gff3")
        doc["proteins"] = self._build_rapid_ftp_url(meta, "pep")
        doc["transcripts"] = self._build_rapid_ftp_url(meta, "cdna")
        doc["ftp_dumps"] = self._build_rapid_ftp_url(meta, "base")
        
        doc["rapid_link"] = f"https://beta.ensembl.org/species/{meta.genome_uuid}"
        
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

    def _build_ftp_url(self, meta: GenomeMetadata, category: str, file_suffix: str, ftp_species_name: str) -> str:
        """Helper to build standard FTP URLs for standard/mouse schemas."""
        date = meta.annotation_date
        if date:
            date = date.replace("-", "_")
        if not date:
            date = "unknown_date"
        source = "ensembl"
        return f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/{source}/{category}/{date}/{file_suffix}"
        
    def _build_rapid_ftp_url(self, meta: GenomeMetadata, resource_type: str) -> str:
        """Helper to build Rapid Release FTP URLs"""
        base = f"https://ftp.ensembl.org/pub/rapid-release/species/{meta.species_name}"
        return f"{base}/{resource_type}"

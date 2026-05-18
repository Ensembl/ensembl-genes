"""
Converts internal GenomeMetadata objects into specific project YAML schemas.
"""
import copy
from typing import Dict, Any, Optional, List
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.config import ProjectConfig
from ensembl.genes.projects.ftp_client import check_url_status
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
                # NCBI LineageEx is ordered root→leaf; reverse to leaf→root so that
                # the first icon match is always the most-specific classification,
                # exactly replicating the legacy core-DB species.classification ordering.
                lineage = list(reversed(lineage))
                return [t.get('ScientificName') for t in lineage if 'ScientificName' in t]
        except Exception as e:
            logger.warning(f"Failed to fetch taxonomy for {taxon_id}: {e}")
        return []

    def _normalise_species_for_ftp(self, species_name: str) -> List[str]:
        import re
        variants = []
        
        # Ensure first letter is capitalized without lowercasing the rest
        if not species_name:
            return variants
        capitalized_species = species_name[0].upper() + species_name[1:]

        # A. Canonical Ensembl-style: replace spaces, hyphens, dots with underscores
        canonical = re.sub(r'[ \-\.]', '_', capitalized_species)
        # C. Collapse multiple underscores
        canonical = re.sub(r'_+', '_', canonical)
        variants.append(canonical)

        # B. Minimal normalisation (current behavior fallback): replace spaces with underscores only
        minimal = capitalized_species.replace(" ", "_")
        if minimal not in variants:
            variants.append(minimal)

        # C. Collapse multiple underscores for minimal
        minimal_collapsed = re.sub(r'_+', '_', minimal)
        if minimal_collapsed not in variants:
            variants.append(minimal_collapsed)

        # D. Lowercase variants (for pre-release/repeat paths if needed)
        lowercase_variants = [v.lower() for v in variants]
        for lv in lowercase_variants:
            if lv not in variants:
                variants.append(lv)

        # Preserve legacy fallbacks just in case
        legacy_v2 = minimal.replace(".", "")
        if legacy_v2 not in variants:
            variants.append(legacy_v2)

        legacy_v3 = minimal.replace(".", "_")
        if legacy_v3 not in variants:
            variants.append(legacy_v3)

        legacy_v4 = re.sub(r'_+', '_', legacy_v3)
        if legacy_v4 not in variants:
            variants.append(legacy_v4)
            
        # Reorder if we have a cached success for this species
        if hasattr(self, '_ftp_species_cache') and species_name in self._ftp_species_cache:
            known_good = self._ftp_species_cache[species_name]
            if known_good in variants:
                variants.remove(known_good)
            variants.insert(0, known_good)

        return variants

    def _resolve_ftp_assets(self, meta: GenomeMetadata) -> Dict[str, Any]:
        import re
        if not hasattr(self, '_ftp_species_cache'):
            self._ftp_species_cache = {}
            
        # Original old base logic for checking if we used a fallback
        ftp_species_name_base = meta.species_name.capitalize().replace(" ", "_")
        variants = self._normalise_species_for_ftp(meta.species_name)
        
        metadata_date = meta.annotation_date.replace("-", "_") if meta.annotation_date else "unknown_date"
        source = (meta.annotation_source or "ensembl").lower().strip()
        
        target_Released = meta.is_released
        resolved_rel_variant = None
        resolved_date = metadata_date
        audit_decision = "excluded"
        audit_reason = "No released or pre-release FTP assets found."
        
        # 1. Try Released logic
        if target_Released:
            for variant in variants:
                base_url = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{variant}/{meta.accession}/{source}/geneset/"
                
                # 1a. Try metadata date first
                gtf_test = f"{base_url}{metadata_date}/genes.gtf.gz"
                gff3_test = f"{base_url}{metadata_date}/genes.gff3.gz"
                if check_url_status(gtf_test) and check_url_status(gff3_test):
                    resolved_rel_variant = variant
                    resolved_date = metadata_date
                    break
                    
                # 1b. If metadata date fails, list directories
                try:
                    response = requests.get(base_url, timeout=10)
                    if response.status_code == 200:
                        matches = re.findall(r'href="(\d{4}_\d{2})/?', response.text)
                        dates = sorted(list(set(matches)), reverse=True)
                        for d in dates:
                            gtf_d = f"{base_url}{d}/genes.gtf.gz"
                            gff3_d = f"{base_url}{d}/genes.gff3.gz"
                            if check_url_status(gtf_d) and check_url_status(gff3_d):
                                resolved_rel_variant = variant
                                resolved_date = d
                                logger.info(f"Resolved FTP date {d} differs from metadata {metadata_date}")
                                break
                except Exception as e:
                    logger.debug(f"Error fetching directory {base_url}: {e}")
                    
                if resolved_rel_variant:
                    break
                    
            if resolved_rel_variant:
                self._ftp_species_cache[meta.species_name] = resolved_rel_variant
                if resolved_rel_variant != ftp_species_name_base:
                    logger.info(f"Resolved FTP species name:\n      input=\"{meta.species_name}\"\n      used=\"{resolved_rel_variant}\"")
                
                audit_decision = "included_released"
                audit_reason = "Found released FTP assets."
                return {
                    "is_released": True,
                    "ftp_species_name": resolved_rel_variant,
                    "resolved_date": resolved_date,
                    "audit_decision": audit_decision,
                    "audit_reason": audit_reason
                }
        
        # 2. Try Pre-release Fallback
        resolved_pre_variant = None
        pre_release_urls = {}
        if self.ftp_client:
            for variant in variants:
                fb_gtf = self.ftp_client.check_pre_release_file(variant, meta.accession, ".gtf.gz")
                if fb_gtf:
                    resolved_pre_variant = variant
                    pre_release_urls["annotation_gtf"] = fb_gtf
                    
                    fb_gff = self.ftp_client.check_pre_release_file(variant, meta.accession, ".gff3.gz")
                    if not fb_gff:
                        fb_gff = self.ftp_client.check_pre_release_file(variant, meta.accession, ".gff3")
                    if fb_gff: pre_release_urls["annotation_gff3"] = fb_gff
                    
                    fb_pep = self.ftp_client.check_pre_release_file(variant, meta.accession, ".pep.fa.gz")
                    if fb_pep: pre_release_urls["proteins"] = fb_pep
                    
                    fb_cdna = self.ftp_client.check_pre_release_file(variant, meta.accession, ".cdna.fa.gz")
                    if fb_cdna: pre_release_urls["transcripts"] = fb_cdna
                    
                    fb_soft = self.ftp_client.check_pre_release_file(variant, meta.accession, ".dna.softmasked.fa.gz")
                    if fb_soft: pre_release_urls["softmasked_genome"] = fb_soft
                    break

        if resolved_pre_variant:
            self._ftp_species_cache[meta.species_name] = resolved_pre_variant
            if resolved_pre_variant != ftp_species_name_base:
                logger.info(f"Resolved FTP species name:\n      input=\"{meta.species_name}\"\n      used=\"{resolved_pre_variant}\"")
                
            audit_decision = "included_prerelease"
            audit_reason = "Found pre-release FTP assets."
            return {
                "is_released": False,
                "ftp_species_name": resolved_pre_variant,
                "resolved_date": metadata_date,
                "audit_decision": audit_decision,
                "audit_reason": audit_reason,
                "pre_release_urls": pre_release_urls
            }
            
        return {
            "is_released": False,
            "ftp_species_name": ftp_species_name_base,
            "resolved_date": metadata_date,
            "audit_decision": "excluded",
            "audit_reason": audit_reason
        }

    def _render_standard(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema A: Standard Projects (VGP, DToL, ERGA)"""
        doc: Dict[str, Any] = {}
        
        # species display name must only come from the scientific name — never from strain,
        # sample description, habitat text, or any other free-text metadata field.
        # (strain is displayed separately in _render_mouse; it must never appear here.)
        doc["species"] = meta.species_name
        
        # Icon mapping — mirrors legacy write_yaml.py priority logic exactly:
        # class_list is already leaf→root (reversed in _fetch_taxonomy_classes), so
        # first-match-wins picks the most-specific mapped classification.
        # Chordata fallback: if nothing more specific matched, chordates get Chordates.png
        # rather than the generic Metazoa.png.
        if self.config.project_name in ["vgp", "dtol", "erga", "darwin_tree_of_life", "cbp", "bge", "asg"]:
            icon = "Metazoa.png"
            class_list = self._fetch_taxonomy_classes(meta.taxon_id)
            for classification in class_list:
                if classification in self.icons:
                    icon = self.icons[classification]
                    break  # first match wins (list is leaf→root, so most-specific class wins)
            if icon == "Metazoa.png" and "Chordata" in class_list:
                icon = "Chordates.png"
            doc["image"] = icon
        elif self.config.scrape_ncbi_submitter and meta.assembly_submitter:
            doc["submitted_by"] = meta.assembly_submitter
            
        doc["accession"] = meta.accession
        doc["annotation_method"] = meta.annotation_method or "BRAKER2"
        
        ftp_resolution = self._resolve_ftp_assets(meta)
        target_Released = ftp_resolution["is_released"]
        ftp_species_name = ftp_resolution["ftp_species_name"]
        
        doc["__audit_decision__"] = ftp_resolution["audit_decision"]
        doc["__audit_reason__"] = ftp_resolution["audit_reason"]
        doc["__audit_resolved_date__"] = ftp_resolution["resolved_date"]
        
        if ftp_resolution["audit_decision"] == "excluded":
            return doc # Returns only audit keys

        if target_Released:
            meta.annotation_date = ftp_resolution["resolved_date"].replace("_", "-")
            doc["annotation_gtf"] = self._build_ftp_url(meta, "geneset", "genes.gtf.gz", ftp_species_name)
            doc["annotation_gff3"] = self._build_ftp_url(meta, "geneset", "genes.gff3.gz", ftp_species_name)
            doc["proteins"] = self._build_ftp_url(meta, "geneset", "pep.fa.gz", ftp_species_name)
            doc["transcripts"] = self._build_ftp_url(meta, "geneset", "cdna.fa.gz", ftp_species_name)
            doc["softmasked_genome"] = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/genome/softmasked.fa.gz"
            
            # repeat_library — checked before emitting; omitted if file does not exist
            repeat_species = ftp_species_name.lower()
            repeat_url = f"https://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/{repeat_species}/{meta.accession}.repeatmodeler.fa"
            if check_url_status(repeat_url):
                doc["repeat_library"] = repeat_url
            
            if not check_url_status(doc["annotation_gff3"]):
                uncompressed_gff = self._build_ftp_url(meta, "geneset", "genes.gff3", ftp_species_name)
                if check_url_status(uncompressed_gff):
                    doc["annotation_gff3"] = uncompressed_gff
        else:
            pre_urls = ftp_resolution.get("pre_release_urls", {})
            for k, v in pre_urls.items():
                doc[k] = v
        
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
        if meta.parent_of_origin:
            doc["parent_of_origin"] = meta.parent_of_origin
            
        doc["assembly_accession"] = meta.accession
        doc["assembly_link"] = f"https://www.ebi.ac.uk/ena/browser/view/{meta.accession}"
        
        if meta.assembly_submitter:
            doc["assembly_submitter"] = meta.assembly_submitter
            
        ftp_resolution = self._resolve_ftp_assets(meta)
        target_Released = ftp_resolution["is_released"]
        ftp_species_name = ftp_resolution["ftp_species_name"]
        
        doc["__audit_decision__"] = ftp_resolution["audit_decision"]
        doc["__audit_reason__"] = ftp_resolution["audit_reason"]
        doc["__audit_resolved_date__"] = ftp_resolution["resolved_date"]
        
        if ftp_resolution["audit_decision"] == "excluded":
            return doc # Returns only audit keys

        if target_Released:
            meta.annotation_date = ftp_resolution["resolved_date"].replace("_", "-")
            doc["annotation_gtf"] = self._build_ftp_url(meta, "geneset", "genes.gtf.gz", ftp_species_name)
            doc["annotation_gff3"] = self._build_ftp_url(meta, "geneset", "genes.gff3.gz", ftp_species_name)
            doc["proteins"] = self._build_ftp_url(meta, "geneset", "pep.fa.gz", ftp_species_name)
            doc["transcripts"] = self._build_ftp_url(meta, "geneset", "cdna.fa.gz", ftp_species_name)
        else:
            pre_urls = ftp_resolution.get("pre_release_urls", {})
            for k, v in pre_urls.items():
                doc[k] = v
        
        vep_url = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/vep/ensembl/geneset/"
        if check_url_status(vep_url):
            doc["variants_vep"] = vep_url
            
        doc["ftp_dumps"] = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/"
        doc["beta_link"] = f"https://beta.ensembl.org/species/{meta.genome_uuid}"

        return {k: v for k, v in doc.items() if v is not None}

    def _render_mouse(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema C: Mouse Genomes"""
        doc: Dict[str, Any] = {}
        
        doc["species"] = meta.species_name
        if meta.strain:
            doc["strain"] = meta.strain
            
        doc["accession"] = meta.accession
        
        ftp_resolution = self._resolve_ftp_assets(meta)
        target_Released = ftp_resolution["is_released"]
        ftp_species_name = ftp_resolution["ftp_species_name"]
        
        doc["__audit_decision__"] = ftp_resolution["audit_decision"]
        doc["__audit_reason__"] = ftp_resolution["audit_reason"]
        doc["__audit_resolved_date__"] = ftp_resolution["resolved_date"]
        
        if ftp_resolution["audit_decision"] == "excluded":
            return doc # Returns only audit keys

        if target_Released:
            meta.annotation_date = ftp_resolution["resolved_date"].replace("_", "-")
            doc["annotation_gtf"] = self._build_ftp_url(meta, "geneset", "genes.gtf.gz", ftp_species_name)
            doc["annotation_gff3"] = self._build_ftp_url(meta, "geneset", "genes.gff3.gz", ftp_species_name)
            doc["proteins"] = self._build_ftp_url(meta, "geneset", "pep.fa.gz", ftp_species_name)
            doc["transcripts"] = self._build_ftp_url(meta, "geneset", "cdna.fa.gz", ftp_species_name)
            doc["softmasked_genome"] = self._build_ftp_url(meta, "genome", "softmasked.fa.gz", ftp_species_name)
        else:
            pre_urls = ftp_resolution.get("pre_release_urls", {})
            for k, v in pre_urls.items():
                doc[k] = v
                
        if target_Released:
            doc["ftp_dumps"] = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/"
        else:
            doc["ftp_dumps"] = f"https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/{ftp_species_name}/{meta.accession}/"
        
        if self.config.allow_beta_urls:
            doc["beta_link"] = f"https://beta.ensembl.org/species/{meta.genome_uuid}"
            
        if meta.alternate_of:
            doc["alternate"] = meta.alternate_of
            
        return {k: v for k, v in doc.items() if v is not None}

    def _build_ftp_url(self, meta: GenomeMetadata, category: str, file_suffix: str, ftp_species_name: str) -> str:
        """Helper to build standard FTP URLs for standard/mouse schemas."""
        date = meta.annotation_date
        if date:
            date = date.replace("-", "_")
        if not date:
            date = "unknown_date"
        source = (meta.annotation_source or "ensembl").lower().strip()
        return f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{ftp_species_name}/{meta.accession}/{source}/{category}/{date}/{file_suffix}"
        
    def _build_rapid_ftp_url(self, meta: GenomeMetadata, resource_type: str) -> str:
        """Helper to build Rapid Release FTP URLs"""
        base = f"https://ftp.ensembl.org/pub/rapid-release/species/{meta.species_name}"
        return f"{base}/{resource_type}"

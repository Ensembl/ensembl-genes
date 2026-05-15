import re
import requests
from ensembl.genes.projects.write_yaml import check_url_status
import logging

logger = logging.getLogger(__name__)

def _resolve_ftp_assets(self, meta):
    """
    Shared FTP resolver that:
    - tries metadata date first
    - lists available YYYY_MM geneset dirs if metadata date fails
    - respects annotation_source such as ensembl / braker
    - tries pre-release fallback
    - returns resolved URLs plus an audit decision
    """
    ftp_species_name_base = meta.species_name.capitalize().replace(" ", "_")
    variants = self._get_ftp_species_variants(meta.species_name)
    
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
        audit_decision = "included_prerelease"
        audit_reason = "Found pre-release FTP assets."
        return {
            "is_released": False,
            "ftp_species_name": resolved_pre_variant,
            "resolved_date": metadata_date, # Date is typically unknown/irrelevant for pre-release but keep metadata
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


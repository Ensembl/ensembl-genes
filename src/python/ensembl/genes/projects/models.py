"""
Unified state and domain models for the genome tracking and YAML generation pipeline.
"""
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any


@dataclass
class GenomeMetadata:
    """
    Superset of all possible metadata attributes needed to generate project YAMLs.
    
    This acts as the standardized internal payload, abstracting away the differences
    between HPRC, Mouse Pangenome, and standard projects like VGP/DToL.
    """
    # Core indentifiers
    genome_uuid: str
    dbname: str
    accession: str
    
    # Required Names
    species_name: str
    assembly_name: str
    
    # Optional Taxonomy
    common_name: Optional[str] = None
    strain: Optional[str] = None
    taxon_id: Optional[int] = None
    
    # Optional relationships
    assembly_submitter: Optional[str] = None
    alternate_of: Optional[str] = None  # URL/name string for alternate haplotype
    parent_of_origin: Optional[str] = None  # maternal or paternal
    population: Optional[str] = None    # used by HPRC
    
    # Annotation properties
    annotation_method: Optional[str] = None
    annotation_date: Optional[str] = None
    
    # Quality metrics
    busco_score: Optional[str] = None
    busco_lineage: Optional[str] = None
    
    # External Server Links (Calculated in the renderer, or boolean presence here)
    is_on_rapid: bool = False
    is_on_beta: bool = False
    is_on_main: bool = False
    
    # FTP Resources Validation (Checked existence)
    has_repeat_library: bool = False
    has_variants_clinvar: bool = False
    has_variants_gnomad: bool = False
    has_variants_vep: bool = False
    
    # Extensible payload for unforeseen additions
    extra: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        # Infer parent of origin from assembly name if not already set
        if not self.parent_of_origin and self.assembly_name:
            asm_lower = self.assembly_name.lower()
            if "_mat" in asm_lower or "maternal" in asm_lower:
                self.parent_of_origin = "maternal"
            elif "_pat" in asm_lower or "paternal" in asm_lower:
                self.parent_of_origin = "paternal"

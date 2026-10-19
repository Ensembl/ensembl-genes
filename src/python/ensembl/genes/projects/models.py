"""
Unified state and domain models for the genome tracking and YAML generation pipeline.
"""

from dataclasses import dataclass, field
from typing import Any


@dataclass
class GenomeMetadata:  # pylint: disable=too-many-instance-attributes
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
    common_name: str | None = None
    strain: str | None = None
    taxon_id: int | None = None
    taxonomy_lineage: list[str] | None = None  # leaf→root classification names

    # Optional relationships
    assembly_submitter: str | None = None
    alternate_of: str | None = None  # URL/name string for alternate haplotype
    parent_of_origin: str | None = None  # maternal or paternal
    population: str | None = None  # used by HPRC

    # Annotation properties
    annotation_source: str | None = None
    annotation_method: str | None = None
    annotation_date: str | None = None

    # Quality metrics
    busco_score: str | None = None
    busco_lineage: str | None = None

    # External Server Links (Calculated in the renderer, or boolean presence here)
    is_on_rapid: bool = False
    is_on_beta: bool = False
    is_on_main: bool = False
    is_released: bool = False

    # FTP Resources Validation (Checked existence)
    has_repeat_library: bool = False
    has_variants_clinvar: bool = False
    has_variants_gnomad: bool = False
    has_variants_vep: bool = False

    # Extensible payload for unforeseen additions
    extra: dict[str, Any] = field(default_factory=dict)

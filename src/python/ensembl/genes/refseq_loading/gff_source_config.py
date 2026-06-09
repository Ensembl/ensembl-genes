"""Source-specific configuration for loading GFF3 into Ensembl core databases."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping


@dataclass(frozen=True)
class GffSourceConfig:
    """Configuration describing how one GFF3 source should be interpreted."""

    name: str
    source_label: str
    analysis_logic_name: str
    analysis_program: str
    parsed_gene_feature_types: frozenset[str]
    parsed_transcript_feature_types: frozenset[str]
    biotype_transcript_feature_types: frozenset[str]
    transcript_feature_biotype_map: Mapping[str, str] = field(default_factory=dict)
    exon_gbkey_biotype_map: Mapping[str, str] = field(default_factory=dict)
    gene_biotype_overrides: Mapping[str, str] = field(default_factory=dict)
    transcript_biotype_overrides: Mapping[str, str] = field(default_factory=dict)
    id_prefixes_to_strip: tuple[str, ...] = ("gene-", "rna-", "cds-", "exon-")
    gbkey_attribute: str = "gbkey"
    pseudo_attribute: str = "pseudo"
    transcript_biotype_attribute: str = "transcript_biotype"
    gene_biotype_attribute: str = "gene_biotype"
    gene_name_attributes: tuple[str, ...] = ("Name", "gene")
    transcript_stable_id_attributes: tuple[str, ...] = ("Name",)
    exon_stable_id_attributes: tuple[str, ...] = ()
    translation_stable_id_attributes: tuple[str, ...] = ()
    gene_xref_prefix: str = "GeneID:"
    transcribed_pseudogene_gbkey_token: str = "Transcribed_Pseudogene"
    segment_gbkey_suffix: str = "_segment"
    segment_biotype_prefix: str = "IG"
    default_biotype: str = "protein_coding"
    toplevel_attrib_type_id: int = 6


REFSEQ_CONFIG = GffSourceConfig(
    name="refseq",
    source_label="refseq",
    analysis_logic_name="refseq_import",
    analysis_program="NCBI_RefSeq",
    parsed_gene_feature_types=frozenset({"gene", "pseudogene"}),
    parsed_transcript_feature_types=frozenset(
        {
            "mRNA",
            "transcript",
            "lnc_RNA",
            "snRNA",
            "rRNA",
            "snoRNA",
            "ncRNA",
            "antisense_RNA",
            "scRNA",
            "telomerase_RNA",
            "RNase_P_RNA",
            "SRP_RNA",
            "RNase_MRP_RNA",
            "piRNA",
            "siRNA",
            "tRNA",
            "pseudogenic_tRNA",
            "D_gene_segment",
            "V_gene_segment",
            "J_gene_segment",
            "C_gene_segment",
            "C_region",
        }
    ),
    biotype_transcript_feature_types=frozenset(
        {
            "mRNA",
            "transcript",
            "lnc_RNA",
            "snRNA",
            "rRNA",
            "snoRNA",
            "ncRNA",
            "antisense_RNA",
            "scRNA",
            "telomerase_RNA",
            "RNase_P_RNA",
            "SRP_RNA",
            "RNase_MRP_RNA",
            "piRNA",
            "siRNA",
            "tRNA",
            "pseudogenic_tRNA",
            "C_region",
            "precursor_RNA",
        }
    ),
    transcript_feature_biotype_map={
        "mRNA": "protein_coding",
        "C_region": "IG_C_gene",
    },
    exon_gbkey_biotype_map={
        "ncRNA": "ncRNA",
        "precursor_RNA": "precursor_RNA",
        "C_region": "IG_C_gene",
    },
    gene_biotype_overrides={
        "V_segment": "IG_V_gene",
        "D_segment": "IG_D_gene",
        "J_segment": "IG_J_gene",
        "C_segment": "IG_C_gene",
        "C_region": "IG_C_gene",
    },
    transcript_biotype_overrides={
        "lncRNA": "lncRNA",
        "antisense_RNA": "antisense_RNA",
        "pseudogene": "pseudogene",
        "transcribed_pseudogene": "transcribed_pseudogene",
        "rRNA": "rRNA",
        "snRNA": "snRNA",
        "snoRNA": "snoRNA",
        "tRNA": "tRNA",
        "miRNA": "miRNA",
        "ncRNA": "ncRNA",
        "misc_RNA": "misc_RNA",
        "telomerase_RNA": "telomerase_RNA",
        "RNase_P_RNA": "RNase_P_RNA",
        "SRP_RNA": "SRP_RNA",
        "RNase_MRP_RNA": "RNase_MRP_RNA",
        "IG_V_gene": "IG_V_gene",
        "IG_D_gene": "IG_D_gene",
        "IG_J_gene": "IG_J_gene",
        "IG_C_gene": "IG_C_gene",
    },
)

GENERIC_GFF_CONFIG = GffSourceConfig(
    name="generic",
    source_label="gff",
    analysis_logic_name="gff_import",
    analysis_program="GFF3",
    parsed_gene_feature_types=frozenset({"gene"}),
    parsed_transcript_feature_types=frozenset(
        {
            "mRNA",
            "transcript",
            "lnc_RNA",
            "snRNA",
            "rRNA",
            "snoRNA",
            "ncRNA",
            "antisense_RNA",
            "scRNA",
            "piRNA",
            "siRNA",
            "tRNA",
        }
    ),
    biotype_transcript_feature_types=frozenset(
        {
            "mRNA",
            "transcript",
            "lnc_RNA",
            "snRNA",
            "rRNA",
            "snoRNA",
            "ncRNA",
            "antisense_RNA",
            "scRNA",
            "piRNA",
            "siRNA",
            "tRNA",
        }
    ),
    transcript_feature_biotype_map={"mRNA": "protein_coding"},
    id_prefixes_to_strip=(),
)

ENSEMBL_GFF_CONFIG = GffSourceConfig(
    name="ensembl",
    source_label="ensembl",
    analysis_logic_name="ensembl_gff_import",
    analysis_program="Ensembl_GFF3",
    parsed_gene_feature_types=frozenset(
        {
            "gene",
            "ncRNA_gene",
            "pseudogene",
        }
    ),
    parsed_transcript_feature_types=frozenset(
        {
            "mRNA",
            "transcript",
            "pseudogenic_transcript",
            "lnc_RNA",
            "snRNA",
            "rRNA",
            "snoRNA",
            "ncRNA",
            "antisense_RNA",
            "scRNA",
            "piRNA",
            "siRNA",
            "miRNA",
            "tRNA",
            "vault_RNA",
            "Y_RNA",
            "RNase_MRP_RNA",
            "RNase_P_RNA",
            "telomerase_RNA",
            "D_gene_segment",
            "V_gene_segment",
            "J_gene_segment",
            "C_gene_segment",
            "C_region",
        }
    ),
    biotype_transcript_feature_types=frozenset(
        {
            "mRNA",
            "transcript",
            "pseudogenic_transcript",
            "lnc_RNA",
            "snRNA",
            "rRNA",
            "snoRNA",
            "ncRNA",
            "antisense_RNA",
            "scRNA",
            "piRNA",
            "siRNA",
            "miRNA",
            "tRNA",
            "vault_RNA",
            "Y_RNA",
            "RNase_MRP_RNA",
            "RNase_P_RNA",
            "telomerase_RNA",
            "D_gene_segment",
            "V_gene_segment",
            "J_gene_segment",
            "C_gene_segment",
            "C_region",
        }
    ),
    transcript_feature_biotype_map={
        "mRNA": "protein_coding",
        "pseudogenic_transcript": "pseudogene",
        "C_region": "IG_C_gene",
    },
    gene_biotype_overrides={
        "V_segment": "IG_V_gene",
        "D_segment": "IG_D_gene",
        "J_segment": "IG_J_gene",
        "C_segment": "IG_C_gene",
        "C_region": "IG_C_gene",
    },
    transcript_biotype_overrides={
        "lncRNA": "lncRNA",
        "antisense_RNA": "antisense_RNA",
        "pseudogene": "pseudogene",
        "transcribed_pseudogene": "transcribed_pseudogene",
        "rRNA": "rRNA",
        "snRNA": "snRNA",
        "snoRNA": "snoRNA",
        "tRNA": "tRNA",
        "miRNA": "miRNA",
        "ncRNA": "ncRNA",
        "misc_RNA": "misc_RNA",
        "IG_V_gene": "IG_V_gene",
        "IG_D_gene": "IG_D_gene",
        "IG_J_gene": "IG_J_gene",
        "IG_C_gene": "IG_C_gene",
    },
    id_prefixes_to_strip=(
        "gene:",
        "transcript:",
        "exon:",
        "CDS:",
        "protein:",
        "chromosome:",
    ),
    transcript_biotype_attribute="biotype",
    gene_biotype_attribute="biotype",
    gene_name_attributes=("Name", "gene_id"),
    transcript_stable_id_attributes=("transcript_id", "Name"),
    exon_stable_id_attributes=("exon_id", "Name", "ID"),
    translation_stable_id_attributes=("protein_id", "Name", "ID"),
)

SOURCE_CONFIGS: dict[str, GffSourceConfig] = {
    ENSEMBL_GFF_CONFIG.name: ENSEMBL_GFF_CONFIG,
    GENERIC_GFF_CONFIG.name: GENERIC_GFF_CONFIG,
    REFSEQ_CONFIG.name: REFSEQ_CONFIG,
}


def register_source_config(source_config: GffSourceConfig) -> None:
    """Register a source configuration for later lookup by name."""

    SOURCE_CONFIGS[source_config.name.lower()] = source_config


def available_source_configs() -> tuple[str, ...]:
    """Return the names of registered GFF source configurations."""

    return tuple(sorted(SOURCE_CONFIGS))


def get_source_config(source_name: str) -> GffSourceConfig:
    """Return a registered source configuration by name."""

    try:
        return SOURCE_CONFIGS[source_name.lower()]
    except KeyError as exc:
        available = ", ".join(available_source_configs())
        raise ValueError(
            f"Unknown GFF source config '{source_name}'. Available configs: {available}"
        ) from exc

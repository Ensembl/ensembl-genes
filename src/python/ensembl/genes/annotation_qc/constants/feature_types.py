"""Shared feature-type groupings for standardized annotation QC tables."""

GENE_FEATURE_TYPES = frozenset(
    {"gene", "ncRNA_gene", "pseudogene", "transposable_element_gene"}
)

TRANSCRIPT_FEATURE_TYPES = frozenset(
    {
        "mRNA",
        "transcript",
        "snoRNA",
        "snRNA",
        "miRNA",
        "rRNA",
        "tRNA",
        "lnc_RNA",
        "lncRNA",
        "ncRNA",
        "antisense_RNA",
        "sRNA",
        "scaRNA",
        "pseudogenic_transcript",
        "scRNA",
        "Y_RNA",
    }
)

FIVE_PRIME_UTR_FEATURE_TYPES = frozenset({"five_prime_UTR"})
THREE_PRIME_UTR_FEATURE_TYPES = frozenset({"three_prime_UTR"})
CDS_FEATURE_TYPES = frozenset({"CDS"})
EXON_FEATURE_TYPES = frozenset({"exon"})

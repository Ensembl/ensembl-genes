"""Configuration and constants for the pangenome mapping system."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Set


# Canonical splice sites
CANONICAL_DONOR_SITES = {"GT", "GC"}
CANONICAL_ACCEPTOR_SITES = {"AG"}

# Codon definitions
START_CODONS = {"ATG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Default thresholds
DEFAULT_MIN_IDENTITY = 0.95
DEFAULT_MIN_BLOCK_LENGTH = 10000  # 10kb minimum syntenic block
DEFAULT_MIN_FEATURE_IDENTITY = 0.90
DEFAULT_MIN_MAPQ = 10
DEFAULT_CNV_MAX_TOTAL_COPIES = 3
DEFAULT_CNV_MIN_EXPANSION_COVERAGE = 0.60
DEFAULT_CNV_MIN_SCORE_RATIO = 0.85
DEFAULT_CNV_MAX_LOCUS_OVERLAP = 0.30
DEFAULT_CNV_GROUP_MIN_RECIPROCAL_OVERLAP = 0.50
DEFAULT_CNV_AMBIGUITY_DISTANCE_BP = 10000
DEFAULT_CNV_AMBIGUITY_SCORE_DELTA = 0.01
DEFAULT_ENABLE_PARALOG_REASSIGNMENT = True
DEFAULT_PARALOG_OVERLAP_THRESHOLD = 0.50
DEFAULT_PARALOG_MAX_LOCUS_CANDIDATES = 8
DEFAULT_PARALOG_MAX_ROUNDS = 2
DEFAULT_ENABLE_MISSED_LOCUS_RECOVERY = True
DEFAULT_MISSED_LOCUS_MAX_GENES = 1000
DEFAULT_MISSED_LOCUS_SEARCH_PADDING_BP = 200000
DEFAULT_MISSED_LOCUS_MAX_WINDOW_BP = 800000
DEFAULT_MISSED_LOCUS_MIN_IDENTITY = 0.85
DEFAULT_MISSED_LOCUS_MIN_COVERAGE = 0.60
DEFAULT_ENABLE_BOUNDARY_REFINEMENT = True
DEFAULT_BOUNDARY_REFINE_MIN_COVERAGE = 0.98
DEFAULT_BOUNDARY_REFINE_WINDOW_BP = 500
DEFAULT_BOUNDARY_REFINE_ANCHOR_BP = 80
DEFAULT_BOUNDARY_REFINE_MAX_UNALIGNED_BP = 8
DEFAULT_ENABLE_VALIDATION_RESCUE = True
DEFAULT_VALIDATION_RESCUE_WINDOW_BP = 1200
DEFAULT_VALIDATION_RESCUE_ANCHOR_BP = 120
DEFAULT_VALIDATION_RESCUE_MAX_GENES = 2000
DEFAULT_ENABLE_INTERNAL_REALIGN_RESCUE = True
DEFAULT_INTERNAL_REALIGN_MAX_GENES = 3000
DEFAULT_INTERNAL_REALIGN_REF_FLANK_BP = 3000
DEFAULT_INTERNAL_REALIGN_TARGET_FLANK_BP = 10000
DEFAULT_INTERNAL_REALIGN_MIN_PATH_COVERAGE = 0.98
DEFAULT_INTERNAL_REALIGN_MAX_INTERNAL_GAP_BP = 6
DEFAULT_INTERNAL_REALIGN_TRIGGER_INDEL_JUMP_BP = 6
DEFAULT_INTERNAL_REALIGN_FALLBACK_BACKEND = "edlib"
DEFAULT_INTERNAL_REALIGN_ACCEPT_ONLY_IF_IMPROVED = True
DEFAULT_ENABLE_SPLICE_SHIFT_RESCUE = True
DEFAULT_SPLICE_SHIFT_MAX_BP = 6
DEFAULT_ENABLE_CODON_FRAME_RESCUE = True
DEFAULT_CODON_RESCUE_MAX_BP = 6
DEFAULT_ENABLE_PROTEIN_QC = True
DEFAULT_PROTEIN_QC_MIN_IDENTITY = 0.90
DEFAULT_PROTEIN_QC_MIN_COVERAGE = 0.90


@dataclass
class MappingConfig:
    """Configuration for the mapping pipeline."""
    
    # Input files
    ref_fasta: Path
    ref_gff: Path
    target_fasta: Path
    
    # Output files
    output_gff: Path
    output_stats: Optional[Path] = None
    output_report: Optional[Path] = None
    
    # Filtering
    chromosomes: Optional[Set[str]] = None  # None = all chromosomes
    biotypes: Optional[Set[str]] = None  # None = all biotypes
    
    # Alignment parameters
    min_identity: float = DEFAULT_MIN_IDENTITY
    min_block_length: int = DEFAULT_MIN_BLOCK_LENGTH
    min_mapq: int = DEFAULT_MIN_MAPQ
    primary_only: bool = True
    
    # Feature mapping parameters
    min_feature_identity: float = DEFAULT_MIN_FEATURE_IDENTITY
    validate_splice_sites: bool = True
    validate_codons: bool = True
    
    # CNV handling parameters
    cnv_max_total_copies: int = DEFAULT_CNV_MAX_TOTAL_COPIES
    cnv_min_expansion_coverage: float = DEFAULT_CNV_MIN_EXPANSION_COVERAGE
    cnv_min_score_ratio: float = DEFAULT_CNV_MIN_SCORE_RATIO
    cnv_max_locus_overlap: float = DEFAULT_CNV_MAX_LOCUS_OVERLAP
    cnv_group_min_reciprocal_overlap: float = DEFAULT_CNV_GROUP_MIN_RECIPROCAL_OVERLAP
    cnv_ambiguity_distance_bp: int = DEFAULT_CNV_AMBIGUITY_DISTANCE_BP
    cnv_ambiguity_score_delta: float = DEFAULT_CNV_AMBIGUITY_SCORE_DELTA

    # Synteny-aware paralog reassignment
    enable_paralog_reassignment: bool = DEFAULT_ENABLE_PARALOG_REASSIGNMENT
    paralog_overlap_threshold: float = DEFAULT_PARALOG_OVERLAP_THRESHOLD
    paralog_max_locus_candidates: int = DEFAULT_PARALOG_MAX_LOCUS_CANDIDATES
    paralog_max_rounds: int = DEFAULT_PARALOG_MAX_ROUNDS

    # Second-pass local recovery for missed loci
    enable_missed_locus_recovery: bool = DEFAULT_ENABLE_MISSED_LOCUS_RECOVERY
    missed_locus_max_genes: int = DEFAULT_MISSED_LOCUS_MAX_GENES
    missed_locus_search_padding_bp: int = DEFAULT_MISSED_LOCUS_SEARCH_PADDING_BP
    missed_locus_max_window_bp: int = DEFAULT_MISSED_LOCUS_MAX_WINDOW_BP
    missed_locus_min_identity: float = DEFAULT_MISSED_LOCUS_MIN_IDENTITY
    missed_locus_min_coverage: float = DEFAULT_MISSED_LOCUS_MIN_COVERAGE
    
    # Boundary refinement for uncertain projected edges
    enable_boundary_refinement: bool = DEFAULT_ENABLE_BOUNDARY_REFINEMENT
    boundary_refine_min_coverage: float = DEFAULT_BOUNDARY_REFINE_MIN_COVERAGE
    boundary_refine_window_bp: int = DEFAULT_BOUNDARY_REFINE_WINDOW_BP
    boundary_refine_anchor_bp: int = DEFAULT_BOUNDARY_REFINE_ANCHOR_BP
    boundary_refine_max_unaligned_bp: int = DEFAULT_BOUNDARY_REFINE_MAX_UNALIGNED_BP
    
    # Validation-guided rescue for likely method-induced boundary errors
    enable_validation_rescue: bool = DEFAULT_ENABLE_VALIDATION_RESCUE
    validation_rescue_window_bp: int = DEFAULT_VALIDATION_RESCUE_WINDOW_BP
    validation_rescue_anchor_bp: int = DEFAULT_VALIDATION_RESCUE_ANCHOR_BP
    validation_rescue_max_genes: int = DEFAULT_VALIDATION_RESCUE_MAX_GENES
    enable_internal_realign_rescue: bool = DEFAULT_ENABLE_INTERNAL_REALIGN_RESCUE
    internal_realign_max_genes: int = DEFAULT_INTERNAL_REALIGN_MAX_GENES
    internal_realign_ref_flank_bp: int = DEFAULT_INTERNAL_REALIGN_REF_FLANK_BP
    internal_realign_target_flank_bp: int = DEFAULT_INTERNAL_REALIGN_TARGET_FLANK_BP
    internal_realign_min_path_coverage: float = DEFAULT_INTERNAL_REALIGN_MIN_PATH_COVERAGE
    internal_realign_max_internal_gap_bp: int = DEFAULT_INTERNAL_REALIGN_MAX_INTERNAL_GAP_BP
    internal_realign_trigger_indel_jump_bp: int = DEFAULT_INTERNAL_REALIGN_TRIGGER_INDEL_JUMP_BP
    internal_realign_fallback_backend: str = DEFAULT_INTERNAL_REALIGN_FALLBACK_BACKEND
    internal_realign_accept_only_if_improved: bool = DEFAULT_INTERNAL_REALIGN_ACCEPT_ONLY_IF_IMPROVED
    enable_splice_shift_rescue: bool = DEFAULT_ENABLE_SPLICE_SHIFT_RESCUE
    splice_shift_max_bp: int = DEFAULT_SPLICE_SHIFT_MAX_BP
    enable_codon_frame_rescue: bool = DEFAULT_ENABLE_CODON_FRAME_RESCUE
    codon_rescue_max_bp: int = DEFAULT_CODON_RESCUE_MAX_BP
    enable_protein_qc: bool = DEFAULT_ENABLE_PROTEIN_QC
    protein_qc_min_identity: float = DEFAULT_PROTEIN_QC_MIN_IDENTITY
    protein_qc_min_coverage: float = DEFAULT_PROTEIN_QC_MIN_COVERAGE
    
    # Performance
    threads: int = 4
    
    # Temporary files
    temp_dir: Optional[Path] = None
    keep_temp: bool = False
    
    def __post_init__(self):
        """Convert string paths to Path objects."""
        if isinstance(self.ref_fasta, str):
            self.ref_fasta = Path(self.ref_fasta)
        if isinstance(self.ref_gff, str):
            self.ref_gff = Path(self.ref_gff)
        if isinstance(self.target_fasta, str):
            self.target_fasta = Path(self.target_fasta)
        if isinstance(self.output_gff, str):
            self.output_gff = Path(self.output_gff)
        if self.output_stats and isinstance(self.output_stats, str):
            self.output_stats = Path(self.output_stats)
        if self.output_report and isinstance(self.output_report, str):
            self.output_report = Path(self.output_report)
        if self.temp_dir and isinstance(self.temp_dir, str):
            self.temp_dir = Path(self.temp_dir)
        
        if self.cnv_max_total_copies < 1:
            raise ValueError("cnv_max_total_copies must be >= 1")
        if self.cnv_min_expansion_coverage < 0:
            raise ValueError("cnv_min_expansion_coverage must be >= 0")
        if self.cnv_min_score_ratio < 0:
            raise ValueError("cnv_min_score_ratio must be >= 0")
        if self.cnv_max_locus_overlap < 0:
            raise ValueError("cnv_max_locus_overlap must be >= 0")
        if self.cnv_group_min_reciprocal_overlap < 0:
            raise ValueError("cnv_group_min_reciprocal_overlap must be >= 0")
        if self.cnv_ambiguity_distance_bp < 0:
            raise ValueError("cnv_ambiguity_distance_bp must be >= 0")
        if self.cnv_ambiguity_score_delta < 0:
            raise ValueError("cnv_ambiguity_score_delta must be >= 0")
        if self.paralog_overlap_threshold < 0:
            raise ValueError("paralog_overlap_threshold must be >= 0")
        if self.paralog_max_locus_candidates < 1:
            raise ValueError("paralog_max_locus_candidates must be >= 1")
        if self.paralog_max_rounds < 1:
            raise ValueError("paralog_max_rounds must be >= 1")
        if self.missed_locus_max_genes < 0:
            raise ValueError("missed_locus_max_genes must be >= 0")
        if self.missed_locus_search_padding_bp < 0:
            raise ValueError("missed_locus_search_padding_bp must be >= 0")
        if self.missed_locus_max_window_bp < 1000:
            raise ValueError("missed_locus_max_window_bp must be >= 1000")
        if self.missed_locus_min_identity < 0:
            raise ValueError("missed_locus_min_identity must be >= 0")
        if self.missed_locus_min_coverage < 0:
            raise ValueError("missed_locus_min_coverage must be >= 0")
        if self.boundary_refine_min_coverage < 0:
            raise ValueError("boundary_refine_min_coverage must be >= 0")
        if self.boundary_refine_window_bp < 0:
            raise ValueError("boundary_refine_window_bp must be >= 0")
        if self.boundary_refine_anchor_bp < 8:
            raise ValueError("boundary_refine_anchor_bp must be >= 8")
        if self.boundary_refine_max_unaligned_bp < 0:
            raise ValueError("boundary_refine_max_unaligned_bp must be >= 0")
        if self.validation_rescue_window_bp < 0:
            raise ValueError("validation_rescue_window_bp must be >= 0")
        if self.validation_rescue_anchor_bp < 8:
            raise ValueError("validation_rescue_anchor_bp must be >= 8")
        if self.validation_rescue_max_genes < 0:
            raise ValueError("validation_rescue_max_genes must be >= 0")
        if self.internal_realign_max_genes < 0:
            raise ValueError("internal_realign_max_genes must be >= 0")
        if self.internal_realign_ref_flank_bp < 0:
            raise ValueError("internal_realign_ref_flank_bp must be >= 0")
        if self.internal_realign_target_flank_bp < 0:
            raise ValueError("internal_realign_target_flank_bp must be >= 0")
        if not (0 <= self.internal_realign_min_path_coverage <= 1):
            raise ValueError("internal_realign_min_path_coverage must be between 0 and 1")
        if self.internal_realign_max_internal_gap_bp < 0:
            raise ValueError("internal_realign_max_internal_gap_bp must be >= 0")
        if self.internal_realign_trigger_indel_jump_bp < 0:
            raise ValueError("internal_realign_trigger_indel_jump_bp must be >= 0")
        if self.internal_realign_fallback_backend not in {"edlib"}:
            raise ValueError("internal_realign_fallback_backend must be 'edlib'")
        if self.splice_shift_max_bp < 0:
            raise ValueError("splice_shift_max_bp must be >= 0")
        if self.codon_rescue_max_bp < 0:
            raise ValueError("codon_rescue_max_bp must be >= 0")
        if self.protein_qc_min_identity < 0:
            raise ValueError("protein_qc_min_identity must be >= 0")
        if self.protein_qc_min_coverage < 0:
            raise ValueError("protein_qc_min_coverage must be >= 0")


# GFF3 feature types hierarchy
GFF3_GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene"}
GFF3_TRANSCRIPT_TYPES = {
    "mRNA", "transcript", "lnc_RNA", "ncRNA",
    "rRNA", "tRNA", "snRNA", "snoRNA", "miRNA",
    "pseudogenic_transcript", "unconfirmed_transcript",
    "V_gene_segment", "D_gene_segment", "J_gene_segment", "C_gene_segment"
}
GFF3_EXON_TYPE = "exon"
GFF3_CDS_TYPE = "CDS"
GFF3_UTR_TYPES = {"five_prime_UTR", "three_prime_UTR"}

# Mapping status codes
class MappingStatus:
    """Status codes for feature mapping results."""
    MAPPED = "mapped"
    PARTIAL = "partial"
    UNMAPPED = "unmapped"
    CONFLICT = "conflict"
    NOVEL = "novel"

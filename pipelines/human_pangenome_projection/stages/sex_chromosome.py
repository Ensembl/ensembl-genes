"""Sex chromosome detection and mapping constraints.

Detects X/Y chromosomes in target genome and prevents cross-chromosome
mismapping, especially in PAR (pseudoautosomal) regions.
"""

import gzip
import logging
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

logger = logging.getLogger(__name__)

# PAR regions on GRCh38 (1-based, inclusive)
# These regions are identical on X and Y
PAR1_START = 10001
PAR1_END = 2781479
PAR2_START = 155701383
PAR2_END = 156030895

# Sex-specific marker genes (sequences for alignment-based detection)
# Format: (gene_name, chromosome, sequence_snippet)
# These are short unique sequences from sex-specific genes
MARKER_GENES = {
    # Y-specific markers
    "SRY": ("Y", "ATGCAATCATATGCTTCTGCTATGTTAAGCGTATTCAACAGCGATGATTACAG"),
    "AMELY": ("Y", "CCCTGGGCTCTGTAAAGAATAGTGTGTTGATTCTTTATCCCAGATGTTTCTC"),
    "DDX3Y": ("Y", "ATGGCAGTTATGGACTGCAATGGAGTAGATGGAGATTATGGTGAAGATGGA"),
    
    # X-specific markers
    "XIST": ("X", "GGTTTGTTTTAGAGACAGGGTCTTGCTCTGTCACCCAGGCTGGAGTGCAGTG"),
    "AMELX": ("X", "CCCTGGGCTCTGTAAAGAATAGTGTGTTGATTCTTTATCTCAGATGTTTCTC"),
}

# Common sex chromosome naming patterns
X_PATTERNS = [
    r'^chrX$', r'^X$', r'^chrx$', r'^x$',
    r'^NC_000023\.\d+$',  # RefSeq X
]
Y_PATTERNS = [
    r'^chrY$', r'^Y$', r'^chry$', r'^y$', 
    r'^NC_000024\.\d+$',  # RefSeq Y
]


@dataclass
class SexChromosomeMap:
    """Classification of target sequences as X, Y, or autosomal."""
    
    x_sequences: Set[str] = field(default_factory=set)  # Sequence names classified as X
    y_sequences: Set[str] = field(default_factory=set)  # Sequence names classified as Y
    par_regions: List[Tuple[str, int, int]] = field(default_factory=list)  # (chrom, start, end)
    detection_method: str = "unknown"  # "named", "marker", "none"
    has_x: bool = False
    has_y: bool = False
    
    def get_sex_class(self, seq_name: str) -> Optional[str]:
        """Get sex chromosome classification for a sequence."""
        if seq_name in self.x_sequences:
            return "X"
        if seq_name in self.y_sequences:
            return "Y"
        return None
    
    def is_par_region(self, seq_name: str, start: int, end: int) -> bool:
        """Check if a region is in a PAR."""
        sex_class = self.get_sex_class(seq_name)
        if sex_class not in ("X", "Y"):
            return False
        
        # Check against known PAR coordinates
        if (start <= PAR1_END and end >= PAR1_START):
            return True
        if (start <= PAR2_END and end >= PAR2_START):
            return True
        return False
    
    def should_allow_mapping(self, ref_chrom: str, target_chrom: str, 
                             ref_start: int = 0, ref_end: int = 0) -> bool:
        """Check if mapping from ref_chrom to target_chrom is allowed."""
        ref_sex = self._classify_ref_chrom(ref_chrom)
        target_sex = self.get_sex_class(target_chrom)
        
        # Non-sex chromosomes: allow any mapping
        if ref_sex is None:
            return True
        
        # Target is unknown/autosomal: allow if ref is sex chrom but target unknown
        if target_sex is None:
            return True
        
        # PAR regions can map to either X or Y
        if self.is_par_region(ref_chrom, ref_start, ref_end):
            return target_sex in ("X", "Y")
        
        # Strict X->X, Y->Y for non-PAR
        return ref_sex == target_sex
    
    def _classify_ref_chrom(self, chrom: str) -> Optional[str]:
        """Classify reference chromosome as X, Y, or None."""
        for pattern in X_PATTERNS:
            if re.match(pattern, chrom, re.IGNORECASE):
                return "X"
        for pattern in Y_PATTERNS:
            if re.match(pattern, chrom, re.IGNORECASE):
                return "Y"
        return None


def detect_named_sex_chromosomes(fasta_path: Union[str, Path]) -> SexChromosomeMap:
    """
    Tier 1: Detect sex chromosomes by sequence name.
    
    Args:
        fasta_path: Path to target FASTA file
        
    Returns:
        SexChromosomeMap with classified sequences
    """
    result = SexChromosomeMap()
    fasta_path = Path(fasta_path)
    
    # Read sequence names from FASTA headers (handles both plain and gzipped files)
    seq_names = []
    try:
        opener = gzip.open if str(fasta_path).endswith('.gz') else open
        with opener(fasta_path, 'rt') as f:
            for line in f:
                if line.startswith('>'):
                    # Extract sequence name (first word after >)
                    name = line[1:].split()[0].strip()
                    seq_names.append(name)
    except Exception as e:
        logger.warning(f"Error reading FASTA headers: {e}")
        return result
    
    # Check for X patterns
    for name in seq_names:
        for pattern in X_PATTERNS:
            if re.match(pattern, name, re.IGNORECASE):
                result.x_sequences.add(name)
                result.has_x = True
                break
    
    # Check for Y patterns
    for name in seq_names:
        for pattern in Y_PATTERNS:
            if re.match(pattern, name, re.IGNORECASE):
                result.y_sequences.add(name)
                result.has_y = True
                break
    
    if result.has_x or result.has_y:
        result.detection_method = "named"
        logger.info(f"Sex chromosomes detected by name: X={result.has_x}, Y={result.has_y}")
        if result.x_sequences:
            logger.info(f"  X sequences: {result.x_sequences}")
        if result.y_sequences:
            logger.info(f"  Y sequences: {result.y_sequences}")
    
    return result


def detect_marker_sex_chromosomes(
    fasta_path: Union[str, Path],
    threads: int = 4
) -> SexChromosomeMap:
    """
    Tier 2: Detect sex chromosomes by aligning marker gene sequences.
    
    Args:
        fasta_path: Path to target FASTA file
        threads: Number of threads for minimap2
        
    Returns:
        SexChromosomeMap with classified sequences
    """
    result = SexChromosomeMap()
    fasta_path = Path(fasta_path)
    
    # Create temp file with marker sequences
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        markers_fa = Path(f.name)
        for gene_name, (chrom, seq) in MARKER_GENES.items():
            f.write(f">{gene_name}_{chrom}\n{seq}\n")
    
    try:
        # Align markers to target
        with tempfile.NamedTemporaryFile(suffix='.paf', delete=False) as paf_f:
            paf_path = Path(paf_f.name)
        
        cmd = [
            "minimap2",
            "-x", "sr",  # Short read preset for short sequences
            "-t", str(threads),
            "-c",  # Output CIGAR
            str(fasta_path),
            str(markers_fa),
            "-o", str(paf_path)
        ]
        
        subprocess.run(cmd, capture_output=True, check=True)
        
        # Parse alignments
        x_hits = {}  # scaffold -> count
        y_hits = {}
        
        with open(paf_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                
                query = parts[0]  # e.g., "SRY_Y"
                target_name = parts[5]
                mapq = int(parts[11])
                
                if mapq < 20:  # Low quality
                    continue
                
                marker_chrom = query.split('_')[-1]  # "X" or "Y"
                
                if marker_chrom == "X":
                    x_hits[target_name] = x_hits.get(target_name, 0) + 1
                elif marker_chrom == "Y":
                    y_hits[target_name] = y_hits.get(target_name, 0) + 1
        
        # Classify scaffolds by hits
        for scaffold in x_hits:
            if scaffold not in y_hits or x_hits[scaffold] > y_hits.get(scaffold, 0):
                result.x_sequences.add(scaffold)
                result.has_x = True
        
        for scaffold in y_hits:
            if scaffold not in x_hits or y_hits[scaffold] > x_hits.get(scaffold, 0):
                result.y_sequences.add(scaffold)
                result.has_y = True
        
        if result.has_x or result.has_y:
            result.detection_method = "marker"
            logger.info(f"Sex chromosomes detected by markers: X={result.has_x}, Y={result.has_y}")
        
    except subprocess.CalledProcessError as e:
        logger.warning(f"Marker alignment failed: {e}")
    except Exception as e:
        logger.warning(f"Marker detection error: {e}")
    finally:
        # Cleanup temp files
        markers_fa.unlink(missing_ok=True)
        if 'paf_path' in locals():
            paf_path.unlink(missing_ok=True)
    
    return result


def detect_sex_chromosomes(
    fasta_path: Union[str, Path],
    threads: int = 4
) -> SexChromosomeMap:
    """
    Detect sex chromosomes using tiered approach.
    
    1. First try named detection (chrX, chrY)
    2. Fall back to marker-based detection if no named chromosomes found
    
    Args:
        fasta_path: Path to target FASTA file
        threads: Number of threads for minimap2
        
    Returns:
        SexChromosomeMap with classification
    """
    logger.info("Detecting sex chromosomes in target genome...")
    
    # Tier 1: Named detection
    result = detect_named_sex_chromosomes(fasta_path)
    
    if result.has_x or result.has_y:
        return result
    
    # Tier 2: Marker-based detection
    logger.info("No named sex chromosomes found, trying marker-based detection...")
    result = detect_marker_sex_chromosomes(fasta_path, threads)
    
    if not result.has_x and not result.has_y:
        logger.warning("Could not detect sex chromosomes in target genome")
        result.detection_method = "none"
    
    return result

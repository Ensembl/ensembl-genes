"""Stage 1: Synteny detection using whole-genome alignment.

Uses minimap2 to align reference and target genomes and identify syntenic blocks.
"""

import logging
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple, Union

from src.models import GenomicInterval, Strand, SyntenicBlock, SyntenicMap
from src.config import DEFAULT_MIN_BLOCK_LENGTH, DEFAULT_MIN_IDENTITY

logger = logging.getLogger(__name__)


@dataclass
class PAFRecord:
    """A parsed PAF (Pairwise Alignment Format) record."""
    
    # Query (target genome in our case)
    query_name: str
    query_length: int
    query_start: int  # 0-based
    query_end: int    # 0-based, exclusive
    
    # Target (reference genome in our case)
    strand: str
    target_name: str
    target_length: int
    target_start: int  # 0-based
    target_end: int    # 0-based, exclusive
    
    # Alignment stats
    matches: int
    alignment_length: int
    mapping_quality: int
    
    # Optional tags
    cs_tag: Optional[str] = None
    cigar: Optional[str] = None
    tp_tag: Optional[str] = None  # P(primary), S(secondary), I(inversion), ...
    chain_score: Optional[int] = None  # s1:i from minimap2
    
    @property
    def query_start_1based(self) -> int:
        return self.query_start + 1
    
    @property
    def query_end_1based(self) -> int:
        return self.query_end
    
    @property
    def target_start_1based(self) -> int:
        return self.target_start + 1
    
    @property
    def target_end_1based(self) -> int:
        return self.target_end
    
    @property
    def identity(self) -> float:
        """Calculate alignment identity."""
        if self.alignment_length == 0:
            return 0.0
        return self.matches / self.alignment_length
    
    @property
    def is_primary(self) -> bool:
        """Return True if alignment is primary (or tp tag missing)."""
        if self.tp_tag is None:
            return True
        return self.tp_tag == "P"
    
    @classmethod
    def from_line(cls, line: str) -> Optional["PAFRecord"]:
        """Parse a PAF line into a record.
        
        Args:
            line: A single PAF line
            
        Returns:
            PAFRecord or None if parsing fails
        """
        parts = line.rstrip().split("\t")
        if len(parts) < 12:
            return None
        
        try:
            record = cls(
                query_name=parts[0],
                query_length=int(parts[1]),
                query_start=int(parts[2]),
                query_end=int(parts[3]),
                strand=parts[4],
                target_name=parts[5],
                target_length=int(parts[6]),
                target_start=int(parts[7]),
                target_end=int(parts[8]),
                matches=int(parts[9]),
                alignment_length=int(parts[10]),
                mapping_quality=int(parts[11])
            )
            
            # Parse optional tags
            for tag in parts[12:]:
                if tag.startswith("cs:Z:"):
                    record.cs_tag = tag[5:]
                elif tag.startswith("cg:Z:"):
                    record.cigar = tag[5:]
                elif tag.startswith("tp:A:"):
                    record.tp_tag = tag[5:]
                elif tag.startswith("s1:i:"):
                    try:
                        record.chain_score = int(tag[5:])
                    except ValueError:
                        pass
            
            return record
        except (ValueError, IndexError) as e:
            logger.warning(f"Failed to parse PAF line: {e}")
            return None


def run_minimap2(
    ref_fasta: Union[str, Path],
    target_fasta: Union[str, Path],
    output_paf: Union[str, Path],
    threads: int = 4,
    preset: str = "asm5"
) -> bool:
    """Run minimap2 for whole-genome alignment.
    
    Args:
        ref_fasta: Path to reference FASTA
        target_fasta: Path to target FASTA
        output_paf: Path for output PAF file
        threads: Number of threads
        preset: minimap2 preset (asm5 for ~0.1% divergence, asm10 for ~1%)
        
    Returns:
        True if successful
    """
    cmd = [
        "minimap2",
        "-cx", preset,     # Preset for assembly-to-assembly alignment
        "--cs=long",       # Long-form cs tag for exact coordinate projection
        "-t", str(threads),
        str(ref_fasta),
        str(target_fasta)
    ]
    
    logger.info(f"Running minimap2: {' '.join(cmd)}")
    
    try:
        with open(output_paf, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        logger.info(f"minimap2 completed, output: {output_paf}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"minimap2 failed: {e.stderr.decode()}")
        return False
    except FileNotFoundError:
        logger.error("minimap2 not found. Please install minimap2.")
        return False


def parse_paf(paf_path: Union[str, Path]) -> Iterator[PAFRecord]:
    """Parse a PAF file.
    
    Args:
        paf_path: Path to PAF file
        
    Yields:
        PAFRecord objects
    """
    with open(paf_path, "r") as f:
        for line in f:
            if line.strip():
                record = PAFRecord.from_line(line)
                if record:
                    yield record


def _paf_has_cs_tag(
    paf_path: Union[str, Path],
    max_records: int = 1000,
) -> bool:
    """Quickly check whether a PAF contains cs tags."""
    checked = 0
    try:
        with open(paf_path, "r") as handle:
            for line in handle:
                if not line.strip():
                    continue
                checked += 1
                if "cs:Z:" in line:
                    return True
                if checked >= max_records:
                    break
    except OSError:
        return False
    return False


def build_syntenic_blocks(
    paf_path: Union[str, Path],
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_block_length: int = DEFAULT_MIN_BLOCK_LENGTH,
    min_mapq: int = 10,
    primary_only: bool = True,
    min_block_quality: Optional[float] = None,
    chromosomes: Optional[Set[str]] = None
) -> List[SyntenicBlock]:
    """Build syntenic blocks from PAF alignment.
    
    Args:
        paf_path: Path to PAF file from minimap2
        min_identity: Minimum alignment identity (0-1)
        min_block_length: Minimum block length in bases
        min_mapq: Minimum mapping quality
        primary_only: Keep only primary alignments from minimap2 tp tag
        min_block_quality: Minimum composite quality score based on
            (identity^2 * alignment_length). Defaults to a modestly stricter
            threshold than raw min_identity/min_block_length cutoffs.
        chromosomes: Optional set of chromosomes to include
        
    Returns:
        List of SyntenicBlock objects
    """
    blocks = []
    block_id = 0
    total_records = 0
    filtered_secondary = 0
    filtered_mapq = 0
    filtered_identity = 0
    filtered_length = 0
    filtered_quality = 0
    quality_threshold = (
        min_block_quality
        if min_block_quality is not None
        else ((min_identity ** 2) * min_block_length * 1.05)
    )
    
    for record in parse_paf(paf_path):
        total_records += 1
        
        if primary_only and not record.is_primary:
            filtered_secondary += 1
            continue
        
        if record.mapping_quality < min_mapq:
            filtered_mapq += 1
            continue
        
        # Apply filters
        # Note: large blocks can have slightly lower average identity because
        # minimap2 includes indels in alignment_length.  Skip the raw identity
        # cutoff when the composite quality score is ≥ 10× the threshold —
        # these blocks are reliable alignments whose cs tags are essential for
        # accurate coordinate projection.
        quality_score = (record.identity ** 2) * record.alignment_length
        high_quality_exempt = quality_score >= quality_threshold * 10

        if record.identity < min_identity and not high_quality_exempt:
            filtered_identity += 1
            continue
        
        if record.alignment_length < min_block_length:
            filtered_length += 1
            continue
        
        if quality_score < quality_threshold:
            filtered_quality += 1
            continue
        
        if chromosomes:
            # Note: in PAF, "target" is reference and "query" is the assembly being mapped
            if record.target_name not in chromosomes:
                continue
        
        # Determine strand
        ref_strand = Strand.PLUS
        target_strand = Strand.from_string(record.strand)
        
        block = SyntenicBlock(
            ref_interval=GenomicInterval(
                seq_region=record.target_name,  # Reference in minimap2 terms
                start=record.target_start_1based,
                end=record.target_end_1based,
                strand=ref_strand
            ),
            target_interval=GenomicInterval(
                seq_region=record.query_name,  # Query/target assembly
                start=record.query_start_1based,
                end=record.query_end_1based,
                strand=target_strand
            ),
            identity=record.identity,
            alignment_length=record.alignment_length,
            matches=record.matches,
            mismatches=record.alignment_length - record.matches,
            cs_tag=record.cs_tag,
            cigar=record.cigar,
            block_id=f"block_{block_id:06d}"
        )
        
        blocks.append(block)
        block_id += 1
    
    logger.info(
        "Built %s syntenic blocks from %s PAF records "
        "(filtered: secondary=%s, low_mapq=%s, low_identity=%s, short=%s, low_quality=%s)",
        len(blocks),
        total_records,
        filtered_secondary,
        filtered_mapq,
        filtered_identity,
        filtered_length,
        filtered_quality
    )
    return blocks


def merge_adjacent_blocks(
    blocks: List[SyntenicBlock],
    max_gap: int = 10000,
    max_identity_diff: float = 0.02
) -> List[SyntenicBlock]:
    """Merge adjacent syntenic blocks that are consistent.
    
    Args:
        blocks: List of syntenic blocks
        max_gap: Maximum gap between blocks to merge
        max_identity_diff: Maximum identity difference for merging
        
    Returns:
        Merged list of blocks
    """
    if not blocks:
        return []
    
    # Sort by reference position
    sorted_blocks = sorted(
        blocks,
        key=lambda b: (b.ref_interval.seq_region, b.ref_interval.start)
    )
    
    merged = []
    current = sorted_blocks[0]
    
    for next_block in sorted_blocks[1:]:
        # Check if blocks can be merged
        can_merge = (
            # Same chromosome
            current.ref_interval.seq_region == next_block.ref_interval.seq_region and
            current.target_interval.seq_region == next_block.target_interval.seq_region and
            # Same orientation
            current.is_inverted == next_block.is_inverted and
            # Close enough in reference
            next_block.ref_interval.start - current.ref_interval.end <= max_gap and
            # Similar identity
            abs(current.identity - next_block.identity) <= max_identity_diff
        )
        
        if can_merge:
            # Target positions should also be adjacent
            if current.is_inverted:
                target_gap = current.target_interval.start - next_block.target_interval.end
            else:
                target_gap = next_block.target_interval.start - current.target_interval.end
            
            if abs(target_gap) <= max_gap:
                # Merge blocks
                current = SyntenicBlock(
                    ref_interval=GenomicInterval(
                        seq_region=current.ref_interval.seq_region,
                        start=current.ref_interval.start,
                        end=next_block.ref_interval.end,
                        strand=current.ref_interval.strand
                    ),
                    target_interval=GenomicInterval(
                        seq_region=current.target_interval.seq_region,
                        start=min(current.target_interval.start, next_block.target_interval.start),
                        end=max(current.target_interval.end, next_block.target_interval.end),
                        strand=current.target_interval.strand
                    ),
                    identity=(current.identity * current.alignment_length + 
                              next_block.identity * next_block.alignment_length) / 
                             (current.alignment_length + next_block.alignment_length),
                    alignment_length=current.alignment_length + next_block.alignment_length,
                    matches=current.matches + next_block.matches,
                    mismatches=current.mismatches + next_block.mismatches,
                    cs_tag=None,  # Can't merge cs tags
                    cigar=None,
                    block_id=current.block_id
                )
                continue
        
        # Can't merge, save current and start new
        merged.append(current)
        current = next_block
    
    merged.append(current)
    
    logger.info(f"Merged {len(blocks)} blocks into {len(merged)} blocks")
    return merged


def create_syntenic_map(
    ref_fasta: Union[str, Path],
    target_fasta: Union[str, Path],
    threads: int = 4,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_block_length: int = DEFAULT_MIN_BLOCK_LENGTH,
    min_mapq: int = 10,
    primary_only: bool = True,
    min_block_quality: Optional[float] = None,
    chromosomes: Optional[Set[str]] = None,
    temp_dir: Optional[Path] = None,
    keep_paf: bool = False,
    paf_output: Optional[Path] = None,
    merge_blocks: bool = False  # Default False to preserve cs tags for exact projection
) -> SyntenicMap:
    """Create a syntenic map between reference and target genomes.
    
    Args:
        ref_fasta: Path to reference FASTA
        target_fasta: Path to target FASTA
        threads: Number of threads for minimap2
        min_identity: Minimum alignment identity
        min_block_length: Minimum block length
        min_mapq: Minimum mapping quality for candidate synteny blocks
        primary_only: Keep only primary minimap2 alignments
        min_block_quality: Minimum composite quality score for a block
        chromosomes: Optional set of chromosomes to include
        temp_dir: Temporary directory for intermediate files
        keep_paf: Whether to keep PAF file after processing
        paf_output: Specific path for PAF output (if keep_paf)
        merge_blocks: Whether to merge adjacent blocks (destroys cs tags)
        
    Returns:
        SyntenicMap with indexed blocks
    """
    # Create PAF file path
    if paf_output:
        paf_path = paf_output
    else:
        tmp_dir = temp_dir or Path(tempfile.gettempdir())
        paf_path = tmp_dir / f"alignment_{Path(ref_fasta).stem}_{Path(target_fasta).stem}.paf"
    
    try:
        # Check if PAF already exists and is valid (non-empty)
        if paf_path.exists() and paf_path.stat().st_size > 0:
            if _paf_has_cs_tag(paf_path):
                logger.info(f"Using existing PAF file: {paf_path}")
            else:
                logger.warning(
                    "Existing PAF lacks cs tags required for exact projection; regenerating: %s",
                    paf_path,
                )
                try:
                    paf_path.unlink()
                except OSError:
                    pass
                success = run_minimap2(
                    ref_fasta, target_fasta, paf_path, threads=threads
                )
                if not success:
                    raise RuntimeError("minimap2 alignment failed")
        else:
            # Run minimap2
            success = run_minimap2(
                ref_fasta, target_fasta, paf_path, threads=threads
            )
            
            if not success:
                raise RuntimeError("minimap2 alignment failed")
        
        # Build blocks from PAF
        blocks = build_syntenic_blocks(
            paf_path,
            min_identity=min_identity,
            min_block_length=min_block_length,
            min_mapq=min_mapq,
            primary_only=primary_only,
            min_block_quality=min_block_quality,
            chromosomes=chromosomes
        )
        
        # Optionally merge adjacent blocks (loses cs tags)
        if merge_blocks:
            blocks = merge_adjacent_blocks(blocks)
        else:
            logger.info(f"Keeping {len(blocks)} unmerged blocks with cs tags")
        
        # Create syntenic map
        syntenic_map = SyntenicMap(blocks=blocks)
        syntenic_map.build_index()
        
        return syntenic_map
    
    finally:
        # Clean up if not keeping
        if not keep_paf and paf_path.exists() and not paf_output:
            paf_path.unlink()


def get_synteny_statistics(syntenic_map: SyntenicMap) -> Dict:
    """Calculate statistics about the syntenic map.
    
    Args:
        syntenic_map: The syntenic map to analyze
        
    Returns:
        Dictionary of statistics
    """
    blocks = syntenic_map.blocks
    
    if not blocks:
        return {"total_blocks": 0}
    
    total_ref_coverage = sum(b.ref_interval.length for b in blocks)
    total_target_coverage = sum(b.target_interval.length for b in blocks)
    
    # Count by chromosome
    ref_by_chr: Dict[str, int] = {}
    target_by_chr: Dict[str, int] = {}
    
    for block in blocks:
        ref_chr = block.ref_interval.seq_region
        ref_by_chr[ref_chr] = ref_by_chr.get(ref_chr, 0) + block.ref_interval.length
        
        target_chr = block.target_interval.seq_region
        target_by_chr[target_chr] = target_by_chr.get(target_chr, 0) + block.target_interval.length
    
    # Count inversions
    inversions = sum(1 for b in blocks if b.is_inverted)
    
    # Identity stats
    identities = [b.identity for b in blocks]
    
    return {
        "total_blocks": len(blocks),
        "total_ref_coverage_bp": total_ref_coverage,
        "total_target_coverage_bp": total_target_coverage,
        "ref_coverage_by_chr": ref_by_chr,
        "target_coverage_by_chr": target_by_chr,
        "inversions": inversions,
        "inversion_rate": inversions / len(blocks) if blocks else 0,
        "mean_identity": sum(identities) / len(identities),
        "min_identity": min(identities),
        "max_identity": max(identities),
        "mean_block_length": total_ref_coverage / len(blocks)
    }

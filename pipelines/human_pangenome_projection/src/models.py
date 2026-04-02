"""Core data models for genomic features and syntenic blocks."""

from bisect import bisect_right
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple, Any


class Strand(Enum):
    """Strand orientation."""
    PLUS = "+"
    MINUS = "-"
    UNSTRANDED = "."
    
    @classmethod
    def from_string(cls, s: str) -> "Strand":
        if s == "+":
            return cls.PLUS
        elif s == "-":
            return cls.MINUS
        return cls.UNSTRANDED
    
    def __str__(self) -> str:
        return self.value


@dataclass
class GenomicInterval:
    """A genomic interval with coordinates."""
    seq_region: str
    start: int  # 1-based, inclusive
    end: int    # 1-based, inclusive
    strand: Strand = Strand.UNSTRANDED
    
    @property
    def length(self) -> int:
        return self.end - self.start + 1
    
    def overlaps(self, other: "GenomicInterval") -> bool:
        """Check if this interval overlaps another on the same strand."""
        if self.seq_region != other.seq_region:
            return False
        return self.start <= other.end and other.start <= self.end
    
    def contains(self, pos: int) -> bool:
        """Check if a position falls within this interval."""
        return self.start <= pos <= self.end
    
    def __lt__(self, other: "GenomicInterval") -> bool:
        if self.seq_region != other.seq_region:
            return self.seq_region < other.seq_region
        return (self.start, self.end) < (other.start, other.end)


@dataclass
class GenomicFeature:
    """Base class for all genomic features."""
    feature_id: str
    feature_type: str
    interval: GenomicInterval
    source: str = "."
    score: Optional[float] = None
    phase: Optional[int] = None  # 0, 1, 2 for CDS
    attributes: Dict[str, Any] = field(default_factory=dict)
    
    # Mapping metadata
    mapping_status: Optional[str] = None
    mapping_identity: Optional[float] = None
    mapped_from: Optional[str] = None  # Original feature ID if mapped
    
    @property
    def seq_region(self) -> str:
        return self.interval.seq_region
    
    @property
    def start(self) -> int:
        return self.interval.start
    
    @property
    def end(self) -> int:
        return self.interval.end
    
    @property
    def strand(self) -> Strand:
        return self.interval.strand
    
    @property
    def length(self) -> int:
        return self.interval.length


@dataclass
class Exon(GenomicFeature):
    """An exon within a transcript."""
    exon_number: Optional[int] = None
    
    def __post_init__(self):
        if self.feature_type != "exon":
            self.feature_type = "exon"


@dataclass
class CDS(GenomicFeature):
    """A coding sequence segment."""
    
    def __post_init__(self):
        if self.feature_type != "CDS":
            self.feature_type = "CDS"


@dataclass
class UTR(GenomicFeature):
    """A UTR (5' or 3')."""
    pass


@dataclass
class Transcript(GenomicFeature):
    """A transcript with its exons and CDS."""
    gene_id: Optional[str] = None
    transcript_name: Optional[str] = None
    biotype: Optional[str] = None
    
    exons: List[Exon] = field(default_factory=list)
    cds_list: List[CDS] = field(default_factory=list)
    utrs: List[UTR] = field(default_factory=list)
    
    # Validation results
    splice_sites_valid: Optional[bool] = None
    start_codon_valid: Optional[bool] = None
    stop_codon_valid: Optional[bool] = None
    frame_valid: Optional[bool] = None
    
    def __post_init__(self):
        if self.feature_type not in {"mRNA", "transcript", "lnc_RNA", "ncRNA", 
                                      "rRNA", "tRNA", "snRNA", "snoRNA", "miRNA",
                                      "pseudogenic_transcript"}:
            self.feature_type = "transcript"
    
    @property
    def is_coding(self) -> bool:
        """Check if transcript is protein-coding."""
        return len(self.cds_list) > 0
    
    @property
    def sorted_exons(self) -> List[Exon]:
        """Return exons sorted by genomic position, accounting for strand."""
        exons = sorted(self.exons, key=lambda e: e.start)
        if self.strand == Strand.MINUS:
            exons = list(reversed(exons))
        return exons
    
    @property
    def sorted_cds(self) -> List[CDS]:
        """Return CDS sorted by genomic position, accounting for strand."""
        cds = sorted(self.cds_list, key=lambda c: c.start)
        if self.strand == Strand.MINUS:
            cds = list(reversed(cds))
        return cds
    
    @property
    def cds_length(self) -> int:
        """Total CDS length in bases."""
        return sum(c.length for c in self.cds_list)
    
    def get_splice_sites(self) -> List[Tuple[int, int]]:
        """Get donor-acceptor pairs for each intron."""
        sorted_exons = sorted(self.exons, key=lambda e: e.start)
        sites = []
        for i in range(len(sorted_exons) - 1):
            donor = sorted_exons[i].end
            acceptor = sorted_exons[i + 1].start
            sites.append((donor, acceptor))
        return sites


@dataclass
class Gene(GenomicFeature):
    """A gene with its transcripts."""
    gene_name: Optional[str] = None
    biotype: Optional[str] = None
    description: Optional[str] = None
    
    transcripts: List[Transcript] = field(default_factory=list)
    
    def __post_init__(self):
        if self.feature_type not in {"gene", "pseudogene", "ncRNA_gene"}:
            self.feature_type = "gene"
    
    @property
    def is_coding(self) -> bool:
        """Check if gene has any coding transcripts."""
        return any(t.is_coding for t in self.transcripts)
    
    def get_transcript(self, transcript_id: str) -> Optional[Transcript]:
        """Get a transcript by ID."""
        for t in self.transcripts:
            if t.feature_id == transcript_id:
                return t
        return None


@dataclass
class SyntenicBlock:
    """A syntenic block between reference and target genomes."""
    ref_interval: GenomicInterval
    target_interval: GenomicInterval
    
    # Alignment information
    identity: float  # 0.0 to 1.0
    alignment_length: int
    matches: int
    mismatches: int
    
    # The original alignment info
    cs_tag: Optional[str] = None  # cigar-style difference string
    cigar: Optional[str] = None
    
    # Block ID for reference
    block_id: Optional[str] = None
    
    # Cached parsed cs tag offsets: list of (ref_pos, cumulative_offset)
    # Populated when build_index() is called on SyntenicMap
    cached_offset_map: Optional[List[Tuple[int, int]]] = None
    
    @property
    def is_inverted(self) -> bool:
        """Check if block is inverted between ref and target."""
        return self.ref_interval.strand != self.target_interval.strand
    
    @property
    def ref_length(self) -> int:
        return self.ref_interval.length
    
    @property
    def target_length(self) -> int:
        return self.target_interval.length
    
    def contains_ref_position(self, pos: int) -> bool:
        """Check if reference position is within this block."""
        return self.ref_interval.contains(pos)


@dataclass
class SyntenicMap:
    """Collection of syntenic blocks for coordinate mapping."""
    blocks: List[SyntenicBlock] = field(default_factory=list)
    
    # Index for fast lookup: ref_seq_region -> sorted blocks
    _ref_index: Dict[str, List[SyntenicBlock]] = field(default_factory=dict)
    _target_index: Dict[str, List[SyntenicBlock]] = field(default_factory=dict)
    _ref_starts: Dict[str, List[int]] = field(default_factory=dict)
    _target_starts: Dict[str, List[int]] = field(default_factory=dict)
    
    def build_index(self):
        """Build spatial index for fast lookups and pre-parse cs tags."""
        self._ref_index.clear()
        self._target_index.clear()
        self._ref_starts.clear()
        self._target_starts.clear()
        
        for block in self.blocks:
            # Pre-parse cs tag if present
            if block.cs_tag and block.cached_offset_map is None:
                block.cached_offset_map = self._parse_cs_to_offsets(block)
            
            ref_chr = block.ref_interval.seq_region
            if ref_chr not in self._ref_index:
                self._ref_index[ref_chr] = []
            self._ref_index[ref_chr].append(block)
            
            target_chr = block.target_interval.seq_region
            if target_chr not in self._target_index:
                self._target_index[target_chr] = []
            self._target_index[target_chr].append(block)
        
        # Sort by start position
        for chrom, chr_blocks in self._ref_index.items():
            chr_blocks.sort(key=lambda b: b.ref_interval.start)
            self._ref_starts[chrom] = [b.ref_interval.start for b in chr_blocks]
        for chrom, chr_blocks in self._target_index.items():
            chr_blocks.sort(key=lambda b: b.target_interval.start)
            self._target_starts[chrom] = [b.target_interval.start for b in chr_blocks]
    
    @staticmethod
    def _parse_cs_to_offsets(block: "SyntenicBlock") -> List[Tuple[int, int]]:
        """Parse cs tag into sparse offset map.
        
        Returns list of (ref_pos, cumulative_offset) tuples.
        """
        offset_changes = []
        cs_tag = block.cs_tag
        if not cs_tag:
            return offset_changes
        
        ref_pos = block.ref_interval.start
        cumulative_offset = 0
        offset_changes.append((ref_pos, cumulative_offset))
        
        pos = 0
        while pos < len(cs_tag):
            ch = cs_tag[pos]
            
            if ch == ':':
                # Short form: :N (N identical bases)
                end = pos + 1
                while end < len(cs_tag) and cs_tag[end].isdigit():
                    end += 1
                length = int(cs_tag[pos+1:end])
                ref_pos += length
                pos = end
                
            elif ch == '=':
                # Long form: =SEQ
                end = pos + 1
                while end < len(cs_tag) and cs_tag[end] in 'ACGTNacgtn':
                    end += 1
                ref_pos += end - pos - 1
                pos = end
                
            elif ch == '*':
                # Substitution: *XY
                ref_pos += 1
                pos += 3
                
            elif ch == '+':
                # Insertion: +SEQ
                end = pos + 1
                while end < len(cs_tag) and cs_tag[end] in 'ACGTNacgtn':
                    end += 1
                ins_len = end - pos - 1
                cumulative_offset += ins_len
                offset_changes.append((ref_pos, cumulative_offset))
                pos = end
                
            elif ch == '-':
                # Deletion: -SEQ
                end = pos + 1
                while end < len(cs_tag) and cs_tag[end] in 'ACGTNacgtn':
                    end += 1
                del_len = end - pos - 1
                ref_pos += del_len
                cumulative_offset -= del_len
                offset_changes.append((ref_pos, cumulative_offset))
                pos = end
                
            elif ch == '~':
                # Intron: ~nnXnn
                end = pos + 1
                while end < len(cs_tag) and (cs_tag[end].isdigit() or cs_tag[end] in 'gtag'):
                    end += 1
                pos = end
                
            else:
                pos += 1
        
        return offset_changes
    
    def find_blocks_for_ref_region(
        self, 
        seq_region: str, 
        start: int, 
        end: int
    ) -> List[SyntenicBlock]:
        """Find all syntenic blocks overlapping a reference region."""
        if seq_region not in self._ref_index:
            return []

        chr_blocks = self._ref_index[seq_region]
        chr_starts = self._ref_starts.get(seq_region, [])
        if not chr_starts:
            return []

        # Only blocks with block.start <= query_end can overlap
        stop_idx = bisect_right(chr_starts, end)
        overlapping = []
        for block in chr_blocks[:stop_idx]:
            if start <= block.ref_interval.end:
                overlapping.append(block)
        return overlapping
    
    def find_block_for_ref_position(
        self, 
        seq_region: str, 
        pos: int
    ) -> Optional[SyntenicBlock]:
        """Find the syntenic block containing a reference position."""
        blocks = self.find_blocks_for_ref_region(seq_region, pos, pos)
        # Return highest identity block if multiple
        if blocks:
            return max(blocks, key=lambda b: b.identity)
        return None


@dataclass
class MappingResult:
    """Result of mapping a feature from reference to target."""
    original_feature: GenomicFeature
    mapped_feature: Optional[GenomicFeature] = None
    
    status: str = "unmapped"  # mapped, partial, unmapped, conflict
    identity: Optional[float] = None
    
    # Validation details
    validation_passed: bool = False
    validation_errors: List[str] = field(default_factory=list)
    
    # For conflict resolution
    conflicting_features: List[GenomicFeature] = field(default_factory=list)
    
    # Source blocks used
    source_blocks: List[SyntenicBlock] = field(default_factory=list)

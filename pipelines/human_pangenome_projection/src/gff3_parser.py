"""GFF3 parser for Ensembl/GENCODE format.

Parses GFF3 files and builds a hierarchical structure of Gene -> Transcript -> Exon/CDS.
"""

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, TextIO, Tuple, Union

from .models import (
    CDS, Exon, Gene, GenomicFeature, GenomicInterval, Strand, Transcript, UTR
)
from .config import (
    GFF3_GENE_TYPES, GFF3_TRANSCRIPT_TYPES, GFF3_EXON_TYPE, 
    GFF3_CDS_TYPE, GFF3_UTR_TYPES
)


def parse_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF3 attribute string into a dictionary.
    
    Args:
        attr_string: Semicolon-separated key=value pairs
        
    Returns:
        Dictionary of attribute key-value pairs
    """
    attributes = {}
    if not attr_string or attr_string == ".":
        return attributes
    
    for item in attr_string.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            # URL decode the value
            value = value.replace("%3B", ";").replace("%3D", "=")
            value = value.replace("%26", "&").replace("%2C", ",")
            attributes[key] = value
    
    return attributes


def format_attributes(attributes: Dict[str, str]) -> str:
    """Format attribute dictionary as GFF3 attribute string.
    
    Args:
        attributes: Dictionary of attribute key-value pairs
        
    Returns:
        Semicolon-separated key=value string
    """
    if not attributes:
        return "."
    
    items = []
    # Order: ID, Parent, Name, then alphabetically
    priority_keys = ["ID", "Parent", "Name", "gene_id", "transcript_id"]
    
    for key in priority_keys:
        if key in attributes:
            value = str(attributes[key])
            # URL encode special characters
            value = value.replace(";", "%3B").replace("=", "%3D")
            value = value.replace("&", "%26").replace(",", "%2C")
            items.append(f"{key}={value}")
    
    for key in sorted(attributes.keys()):
        if key not in priority_keys:
            value = str(attributes[key])
            value = value.replace(";", "%3B").replace("=", "%3D")
            value = value.replace("&", "%26").replace(",", "%2C")
            items.append(f"{key}={value}")
    
    return ";".join(items)


class GFF3Record:
    """A single parsed GFF3 record."""
    
    __slots__ = ['seq_region', 'source', 'feature_type', 'start', 'end', 
                 'score', 'strand', 'phase', 'attributes']
    
    def __init__(
        self,
        seq_region: str,
        source: str,
        feature_type: str,
        start: int,
        end: int,
        score: Optional[float],
        strand: str,
        phase: Optional[int],
        attributes: Dict[str, str]
    ):
        self.seq_region = seq_region
        self.source = source
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
    
    @property
    def id(self) -> Optional[str]:
        return self.attributes.get("ID")
    
    @property
    def parent(self) -> Optional[str]:
        return self.attributes.get("Parent")
    
    @classmethod
    def from_line(cls, line: str) -> Optional["GFF3Record"]:
        """Parse a GFF3 line into a record.
        
        Args:
            line: A single GFF3 line (non-comment, non-empty)
            
        Returns:
            GFF3Record or None if parsing fails
        """
        parts = line.rstrip().split("\t")
        if len(parts) != 9:
            return None
        
        try:
            seq_region = parts[0]
            source = parts[1]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            score = float(parts[5]) if parts[5] != "." else None
            strand = parts[6]
            phase = int(parts[7]) if parts[7] != "." else None
            attributes = parse_attributes(parts[8])
            
            return cls(
                seq_region=seq_region,
                source=source,
                feature_type=feature_type,
                start=start,
                end=end,
                score=score,
                strand=strand,
                phase=phase,
                attributes=attributes
            )
        except (ValueError, IndexError):
            return None


class GFF3Parser:
    """Parser for GFF3 files with hierarchical feature building."""
    
    def __init__(
        self,
        chromosomes: Optional[Set[str]] = None,
        biotypes: Optional[Set[str]] = None,
        include_pseudogenes: bool = True
    ):
        """Initialize parser with optional filters.
        
        Args:
            chromosomes: If set, only parse features on these chromosomes
            biotypes: If set, only parse genes with these biotypes
            include_pseudogenes: Whether to include pseudogenes
        """
        self.chromosomes = chromosomes
        self.biotypes = biotypes
        self.include_pseudogenes = include_pseudogenes
        
        # Storage during parsing
        self._genes: Dict[str, Gene] = {}
        self._transcripts: Dict[str, Transcript] = {}
        self._orphan_features: List[GFF3Record] = []
    
    def _open_file(self, filepath: Union[str, Path]) -> TextIO:
        """Open a file, handling gzip compression."""
        filepath = Path(filepath)
        if filepath.suffix == ".gz":
            return gzip.open(filepath, "rt", encoding="utf-8")
        return open(filepath, "r", encoding="utf-8")
    
    def parse_lines(self, line_iterator: Iterator[str]) -> Dict[str, Gene]:
        """Parse GFF3 lines from an iterator.
        
        Args:
            line_iterator: Iterator yielding GFF3 lines
            
        Returns:
            Dictionary mapping gene_id to Gene objects
        """
        self._genes.clear()
        self._transcripts.clear()
        self._orphan_features.clear()
        
        # First pass: collect all records
        records_by_type: Dict[str, List[GFF3Record]] = defaultdict(list)
        
        for line in line_iterator:
            if line.startswith("#") or not line.strip():
                continue
            
            record = GFF3Record.from_line(line)
            if record is None:
                continue
            
            # Apply chromosome filter
            if self.chromosomes and record.seq_region not in self.chromosomes:
                continue
            
            records_by_type[record.feature_type].append(record)
        
        # Build genes
        for gene_type in GFF3_GENE_TYPES:
            for record in records_by_type.get(gene_type, []):
                self._build_gene(record)
        
        # Build transcripts and link to genes
        for tx_type in GFF3_TRANSCRIPT_TYPES:
            for record in records_by_type.get(tx_type, []):
                self._build_transcript(record)
        
        # Build exons and link to transcripts
        for record in records_by_type.get(GFF3_EXON_TYPE, []):
            self._build_exon(record)
        
        # Build CDS and link to transcripts
        for record in records_by_type.get(GFF3_CDS_TYPE, []):
            self._build_cds(record)
        
        # Build UTRs and link to transcripts
        for utr_type in GFF3_UTR_TYPES:
            for record in records_by_type.get(utr_type, []):
                self._build_utr(record)
        
        # Apply biotype filter if specified
        if self.biotypes:
            self._genes = {
                gid: gene for gid, gene in self._genes.items()
                if gene.biotype in self.biotypes
            }
        
        # Filter out pseudogenes if not wanted
        if not self.include_pseudogenes:
            self._genes = {
                gid: gene for gid, gene in self._genes.items()
                if gene.feature_type != "pseudogene" and gene.biotype != "pseudogene"
            }
        
        return self._genes

    def parse(self, filepath: Union[str, Path]) -> Dict[str, Gene]:
        """Parse a GFF3 file and return genes indexed by gene_id.
        
        Args:
            filepath: Path to GFF3 file (can be gzipped)
            
        Returns:
            Dictionary mapping gene_id to Gene objects
        """
        try:
            with self._open_file(filepath) as f:
                return self.parse_lines(f)
        except OSError as e:
            # Handle potential file errors gracefully if needed, 
            # though raising is appropriate here
            raise e
    
    def _build_gene(self, record: GFF3Record) -> Gene:
        """Build a Gene object from a GFF3 record."""
        gene_id = record.id or record.attributes.get("gene_id", f"gene_{len(self._genes)}")
        
        gene = Gene(
            feature_id=gene_id,
            feature_type=record.feature_type,
            interval=GenomicInterval(
                seq_region=record.seq_region,
                start=record.start,
                end=record.end,
                strand=Strand.from_string(record.strand)
            ),
            source=record.source,
            score=record.score,
            attributes=record.attributes,
            gene_name=record.attributes.get("Name") or record.attributes.get("gene_name"),
            biotype=record.attributes.get("gene_biotype") or record.attributes.get("biotype") or record.attributes.get("gene_type"),
            description=record.attributes.get("description")
        )
        
        self._genes[gene_id] = gene
        return gene
    
    def _build_transcript(self, record: GFF3Record) -> Optional[Transcript]:
        """Build a Transcript object from a GFF3 record."""
        tx_id = record.id or record.attributes.get("transcript_id", f"tx_{len(self._transcripts)}")
        parent_id = record.parent or record.attributes.get("gene_id")
        
        # Find parent gene
        parent_gene = self._genes.get(parent_id) if parent_id else None
        
        transcript = Transcript(
            feature_id=tx_id,
            feature_type=record.feature_type,
            interval=GenomicInterval(
                seq_region=record.seq_region,
                start=record.start,
                end=record.end,
                strand=Strand.from_string(record.strand)
            ),
            source=record.source,
            score=record.score,
            attributes=record.attributes,
            gene_id=parent_id,
            transcript_name=record.attributes.get("Name") or record.attributes.get("transcript_name"),
            biotype=record.attributes.get("transcript_biotype") or record.attributes.get("biotype")
        )
        
        self._transcripts[tx_id] = transcript
        
        if parent_gene:
            parent_gene.transcripts.append(transcript)
        
        return transcript
    
    def _build_exon(self, record: GFF3Record) -> Optional[Exon]:
        """Build an Exon object from a GFF3 record."""
        exon_id = record.id or record.attributes.get("exon_id", f"exon_{id(record)}")
        parent_id = record.parent or record.attributes.get("transcript_id")
        
        # Parse exon number if available
        exon_number = None
        if "exon_number" in record.attributes:
            try:
                exon_number = int(record.attributes["exon_number"])
            except ValueError:
                pass
        
        exon = Exon(
            feature_id=exon_id,
            feature_type="exon",
            interval=GenomicInterval(
                seq_region=record.seq_region,
                start=record.start,
                end=record.end,
                strand=Strand.from_string(record.strand)
            ),
            source=record.source,
            score=record.score,
            phase=record.phase,
            attributes=record.attributes,
            exon_number=exon_number
        )
        
        # Link to parent transcript(s)
        if parent_id:
            for pid in parent_id.split(","):
                pid = pid.strip()
                if pid in self._transcripts:
                    self._transcripts[pid].exons.append(exon)
        
        return exon
    
    def _build_cds(self, record: GFF3Record) -> Optional[CDS]:
        """Build a CDS object from a GFF3 record."""
        cds_id = record.id or record.attributes.get("ID", f"cds_{id(record)}")
        parent_id = record.parent or record.attributes.get("transcript_id")
        
        cds = CDS(
            feature_id=cds_id,
            feature_type="CDS",
            interval=GenomicInterval(
                seq_region=record.seq_region,
                start=record.start,
                end=record.end,
                strand=Strand.from_string(record.strand)
            ),
            source=record.source,
            score=record.score,
            phase=record.phase,
            attributes=record.attributes
        )
        
        # Link to parent transcript(s)
        if parent_id:
            for pid in parent_id.split(","):
                pid = pid.strip()
                if pid in self._transcripts:
                    self._transcripts[pid].cds_list.append(cds)
        
        return cds
    
    def _build_utr(self, record: GFF3Record) -> Optional[UTR]:
        """Build a UTR object from a GFF3 record."""
        utr_id = record.id or f"utr_{id(record)}"
        parent_id = record.parent or record.attributes.get("transcript_id")
        
        utr = UTR(
            feature_id=utr_id,
            feature_type=record.feature_type,
            interval=GenomicInterval(
                seq_region=record.seq_region,
                start=record.start,
                end=record.end,
                strand=Strand.from_string(record.strand)
            ),
            source=record.source,
            score=record.score,
            phase=record.phase,
            attributes=record.attributes
        )
        
        # Link to parent transcript(s)
        if parent_id:
            for pid in parent_id.split(","):
                pid = pid.strip()
                if pid in self._transcripts:
                    self._transcripts[pid].utrs.append(utr)
        
        return utr
    
    def iterate_records(self, filepath: Union[str, Path]) -> Iterator[GFF3Record]:
        """Iterate over raw GFF3 records without building hierarchy.
        
        Args:
            filepath: Path to GFF3 file
            
        Yields:
            GFF3Record objects
        """
        with self._open_file(filepath) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                
                record = GFF3Record.from_line(line)
                if record is None:
                    continue
                
                if self.chromosomes and record.seq_region not in self.chromosomes:
                    continue
                
                yield record


def parse_gff3(
    filepath: Union[str, Path],
    chromosomes: Optional[Set[str]] = None,
    biotypes: Optional[Set[str]] = None
) -> Dict[str, Gene]:
    """Convenience function to parse a GFF3 file.
    
    Args:
        filepath: Path to GFF3 file
        chromosomes: Optional set of chromosomes to include
        biotypes: Optional set of gene biotypes to include
        
    Returns:
        Dictionary mapping gene_id to Gene objects
    """
    parser = GFF3Parser(chromosomes=chromosomes, biotypes=biotypes)
    return parser.parse(filepath)

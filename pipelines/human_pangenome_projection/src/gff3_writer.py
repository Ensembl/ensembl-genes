"""GFF3 writer for Ensembl/GENCODE format output."""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Union

from .models import CDS, Exon, Gene, GenomicFeature, Strand, Transcript, UTR
from .gff3_parser import format_attributes


class GFF3Writer:
    """Writer for GFF3 format output."""
    
    def __init__(
        self,
        output_path: Union[str, Path],
        source: str = "pangenome_mapper",
        include_provenance: bool = True
    ):
        """Initialize GFF3 writer.
        
        Args:
            output_path: Path for output GFF3 file (will use gzip if .gz extension)
            source: Source field for GFF3 records
            include_provenance: Whether to add mapping provenance attributes
        """
        self.output_path = Path(output_path)
        self.source = source
        self.include_provenance = include_provenance
        self._handle: Optional[TextIO] = None
    
    def _open(self) -> TextIO:
        """Open output file for writing."""
        if self.output_path.suffix == ".gz":
            return gzip.open(self.output_path, "wt", encoding="utf-8")
        return open(self.output_path, "w", encoding="utf-8")
    
    def __enter__(self) -> "GFF3Writer":
        self._handle = self._open()
        self._write_header()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._handle:
            self._handle.close()
            self._handle = None
    
    def _write_header(self):
        """Write GFF3 header."""
        self._handle.write("##gff-version 3\n")

    @staticmethod
    def _externalize_id(raw_id: Optional[str]) -> Optional[str]:
        """Convert internal mapped IDs to output-facing IDs."""
        if raw_id is None:
            return None
        external = raw_id
        if external.startswith("mapped_"):
            external = external[len("mapped_"):]
        external = external.replace("__copy", "_copy")
        return external
    
    def _format_feature(
        self,
        feature: GenomicFeature,
        parent_id: Optional[str] = None,
        additional_attrs: Optional[Dict[str, str]] = None
    ) -> str:
        """Format a feature as a GFF3 line.
        
        Args:
            feature: The genomic feature to format
            parent_id: Parent ID to add to attributes
            additional_attrs: Extra attributes to include
            
        Returns:
            Formatted GFF3 line
        """
        # Build attributes
        attrs = dict(feature.attributes)
        attrs["ID"] = self._externalize_id(feature.feature_id)
        
        if parent_id:
            attrs["Parent"] = self._externalize_id(parent_id)
        
        # Add provenance if enabled
        if self.include_provenance:
            if feature.mapped_from:
                attrs["mapped_from"] = feature.mapped_from
            if feature.mapping_identity is not None:
                attrs["mapping_identity"] = f"{feature.mapping_identity:.4f}"
            if feature.mapping_status:
                attrs["mapping_status"] = feature.mapping_status
        
        if additional_attrs:
            attrs.update(additional_attrs)
        
        # Format score
        score = "." if feature.score is None else str(feature.score)
        
        # Format phase
        phase = "." if feature.phase is None else str(feature.phase)
        
        # Build line
        parts = [
            feature.seq_region,
            feature.source or self.source,
            feature.feature_type,
            str(feature.start),
            str(feature.end),
            score,
            str(feature.strand),
            phase,
            format_attributes(attrs)
        ]
        
        return "\t".join(parts)
    
    def write_gene(self, gene: Gene):
        """Write a gene and all its children to the GFF3 file.
        
        Args:
            gene: Gene object with transcripts
        """
        if self._handle is None:
            raise RuntimeError("Writer not opened. Use with statement or call __enter__")
        
        # Write gene
        self._handle.write(self._format_feature(gene) + "\n")
        
        # Write transcripts
        for transcript in gene.transcripts:
            self.write_transcript(transcript, parent_id=gene.feature_id)
    
    def write_transcript(self, transcript: Transcript, parent_id: Optional[str] = None):
        """Write a transcript and all its children.
        
        Args:
            transcript: Transcript object
            parent_id: Parent gene ID
        """
        if self._handle is None:
            raise RuntimeError("Writer not opened")
        
        # Write transcript
        tx_parent = parent_id or transcript.gene_id
        self._handle.write(self._format_feature(transcript, parent_id=tx_parent) + "\n")
        
        # Write exons (sorted by position)
        for exon in sorted(transcript.exons, key=lambda e: e.start):
            self._handle.write(
                self._format_feature(exon, parent_id=transcript.feature_id) + "\n"
            )
        
        # Write CDS (sorted by position)
        for cds in sorted(transcript.cds_list, key=lambda c: c.start):
            self._handle.write(
                self._format_feature(cds, parent_id=transcript.feature_id) + "\n"
            )
        
        # Write UTRs (sorted by position)
        for utr in sorted(transcript.utrs, key=lambda u: u.start):
            self._handle.write(
                self._format_feature(utr, parent_id=transcript.feature_id) + "\n"
            )
    
    def write_genes(self, genes: List[Gene]):
        """Write multiple genes.
        
        Args:
            genes: List of Gene objects
        """
        # Sort genes by chromosome and position
        sorted_genes = sorted(genes, key=lambda g: (g.seq_region, g.start))
        
        for gene in sorted_genes:
            self.write_gene(gene)
    
    def write_directive(self, directive: str):
        """Write a GFF3 directive line.
        
        Args:
            directive: Directive string (without ##)
        """
        if self._handle is None:
            raise RuntimeError("Writer not opened")
        self._handle.write(f"##{directive}\n")
    
    def write_sequence_region(self, seq_region: str, start: int, end: int):
        """Write a sequence-region directive.
        
        Args:
            seq_region: Chromosome/scaffold name
            start: Start coordinate (1-based)
            end: End coordinate
        """
        self.write_directive(f"sequence-region {seq_region} {start} {end}")


def write_genes_to_gff3(
    genes: List[Gene],
    output_path: Union[str, Path],
    source: str = "pangenome_mapper",
    include_provenance: bool = True
):
    """Convenience function to write genes to a GFF3 file.
    
    Args:
        genes: List of Gene objects
        output_path: Output file path
        source: Source field for GFF3
        include_provenance: Whether to include mapping metadata
    """
    with GFF3Writer(output_path, source=source, include_provenance=include_provenance) as writer:
        writer.write_genes(genes)

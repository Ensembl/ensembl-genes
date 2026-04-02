"""Structural validation of mapped features.

Validates splice sites, start/stop codons, and feature integrity.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from src.models import Gene, Strand, Transcript
from src.fasta_handler import FastaHandler
from src.config import (
    CANONICAL_DONOR_SITES, CANONICAL_ACCEPTOR_SITES,
    START_CODONS, STOP_CODONS
)

logger = logging.getLogger(__name__)


@dataclass
class SpliceSiteValidation:
    """Validation result for a splice site."""
    intron_index: int
    donor_pos: int
    acceptor_pos: int
    donor_dinuc: str
    acceptor_dinuc: str
    is_valid: bool
    error: Optional[str] = None


@dataclass
class CodonValidation:
    """Validation result for a codon."""
    codon_type: str  # "start" or "stop"
    position: int
    codon: str
    is_valid: bool
    error: Optional[str] = None


@dataclass
class TranscriptValidation:
    """Complete validation result for a transcript."""
    transcript_id: str
    
    # Splice sites
    splice_sites_checked: int = 0
    splice_sites_valid: int = 0
    splice_site_results: List[SpliceSiteValidation] = field(default_factory=list)
    
    # Codons
    start_codon_valid: bool = False
    stop_codon_valid: bool = False
    start_codon_result: Optional[CodonValidation] = None
    stop_codon_result: Optional[CodonValidation] = None
    
    # Frame
    cds_length_valid: bool = False
    cds_length: int = 0
    
    # Overall
    is_valid: bool = False
    errors: List[str] = field(default_factory=list)


@dataclass
class GeneValidation:
    """Complete validation result for a gene."""
    gene_id: str
    transcripts_checked: int = 0
    transcripts_valid: int = 0
    transcript_results: List[TranscriptValidation] = field(default_factory=list)
    is_valid: bool = False


class StructuralValidator:
    """Validates structural integrity of mapped features."""
    
    def __init__(self, fasta_handler: FastaHandler):
        """Initialize validator.
        
        Args:
            fasta_handler: Handler for target genome FASTA
        """
        self.fasta = fasta_handler
    
    def validate_splice_sites(
        self,
        transcript: Transcript
    ) -> List[SpliceSiteValidation]:
        """Validate splice sites for a transcript.
        
        Args:
            transcript: The transcript to validate
            
        Returns:
            List of SpliceSiteValidation results
        """
        results = []
        
        if len(transcript.exons) < 2:
            return results
        
        # Sort exons by genomic position
        sorted_exons = sorted(transcript.exons, key=lambda e: e.start)
        
        for i in range(len(sorted_exons) - 1):
            exon1 = sorted_exons[i]
            exon2 = sorted_exons[i + 1]
            
            # Get intron boundaries
            if transcript.strand == Strand.PLUS:
                donor_pos = exon1.end
                acceptor_pos = exon2.start
                
                # Get dinucleotides (first 2 bases of intron, last 2 bases of intron)
                donor_dinuc = self.fasta.fetch(
                    exon1.seq_region, donor_pos + 1, donor_pos + 2, Strand.PLUS
                )
                acceptor_dinuc = self.fasta.fetch(
                    exon2.seq_region, acceptor_pos - 2, acceptor_pos - 1, Strand.PLUS
                )
            else:
                # Minus strand: donor at higher coord end of exon
                donor_pos = exon2.start
                acceptor_pos = exon1.end
                
                donor_dinuc = self.fasta.fetch(
                    exon2.seq_region, donor_pos - 2, donor_pos - 1, Strand.MINUS
                )
                acceptor_dinuc = self.fasta.fetch(
                    exon1.seq_region, acceptor_pos + 1, acceptor_pos + 2, Strand.MINUS
                )
            
            is_valid = (
                donor_dinuc in CANONICAL_DONOR_SITES and
                acceptor_dinuc in CANONICAL_ACCEPTOR_SITES
            )
            
            error = None
            if not is_valid:
                error_parts = []
                if donor_dinuc not in CANONICAL_DONOR_SITES:
                    error_parts.append(f"non-canonical donor {donor_dinuc}")
                if acceptor_dinuc not in CANONICAL_ACCEPTOR_SITES:
                    error_parts.append(f"non-canonical acceptor {acceptor_dinuc}")
                error = "; ".join(error_parts)
            
            results.append(SpliceSiteValidation(
                intron_index=i,
                donor_pos=donor_pos,
                acceptor_pos=acceptor_pos,
                donor_dinuc=donor_dinuc,
                acceptor_dinuc=acceptor_dinuc,
                is_valid=is_valid,
                error=error
            ))
        
        return results
    
    def validate_start_codon(
        self,
        transcript: Transcript
    ) -> Optional[CodonValidation]:
        """Validate start codon for a coding transcript.
        
        Args:
            transcript: The transcript to validate
            
        Returns:
            CodonValidation or None if not coding
        """
        if not transcript.cds_list:
            return None
        
        # Get first CDS (in transcript direction)
        sorted_cds = sorted(transcript.cds_list, key=lambda c: c.start)
        
        if transcript.strand == Strand.PLUS:
            first_cds = sorted_cds[0]
            codon_pos = first_cds.start
        else:
            first_cds = sorted_cds[-1]
            codon_pos = first_cds.end
        
        codon = self.fasta.get_codon(
            first_cds.seq_region, codon_pos, transcript.strand
        )
        
        is_valid = codon in START_CODONS
        
        return CodonValidation(
            codon_type="start",
            position=codon_pos,
            codon=codon,
            is_valid=is_valid,
            error=None if is_valid else f"non-ATG start codon: {codon}"
        )
    
    def validate_stop_codon(
        self,
        transcript: Transcript
    ) -> Optional[CodonValidation]:
        """Validate stop codon for a coding transcript.
        
        Args:
            transcript: The transcript to validate
            
        Returns:
            CodonValidation or None if not coding
        """
        if not transcript.cds_list:
            return None
        
        # Get last CDS (in transcript direction)
        sorted_cds = sorted(transcript.cds_list, key=lambda c: c.start)
        
        if transcript.strand == Strand.PLUS:
            last_cds = sorted_cds[-1]
            # Stop codon is at or after CDS end
            codon_pos = last_cds.end - 2  # Last 3 bases
        else:
            last_cds = sorted_cds[0]
            codon_pos = last_cds.start + 2
        
        codon = self.fasta.get_codon(
            last_cds.seq_region, codon_pos, transcript.strand
        )
        
        is_valid = codon in STOP_CODONS
        
        return CodonValidation(
            codon_type="stop",
            position=codon_pos,
            codon=codon,
            is_valid=is_valid,
            error=None if is_valid else f"non-stop codon: {codon}"
        )
    
    def validate_cds_frame(
        self,
        transcript: Transcript
    ) -> Tuple[bool, int]:
        """Check if total CDS length is divisible by 3.
        
        Args:
            transcript: The transcript to validate
            
        Returns:
            Tuple of (is_valid, cds_length)
        """
        if not transcript.cds_list:
            return True, 0
        
        total_length = sum(cds.length for cds in transcript.cds_list)
        
        return total_length % 3 == 0, total_length
    
    def validate_transcript(
        self,
        transcript: Transcript
    ) -> TranscriptValidation:
        """Validate a transcript completely.
        
        Args:
            transcript: The transcript to validate
            
        Returns:
            TranscriptValidation result
        """
        result = TranscriptValidation(transcript_id=transcript.feature_id)
        
        # Validate splice sites
        splice_results = self.validate_splice_sites(transcript)
        result.splice_site_results = splice_results
        result.splice_sites_checked = len(splice_results)
        result.splice_sites_valid = sum(1 for s in splice_results if s.is_valid)
        
        # Validate codons if coding
        if transcript.is_coding:
            start_result = self.validate_start_codon(transcript)
            stop_result = self.validate_stop_codon(transcript)
            
            result.start_codon_result = start_result
            result.stop_codon_result = stop_result
            result.start_codon_valid = start_result.is_valid if start_result else False
            result.stop_codon_valid = stop_result.is_valid if stop_result else False
            
            # Validate frame
            frame_valid, cds_length = self.validate_cds_frame(transcript)
            result.cds_length_valid = frame_valid
            result.cds_length = cds_length
            
            # Collect errors
            if start_result and not start_result.is_valid:
                result.errors.append(start_result.error)
            if stop_result and not stop_result.is_valid:
                result.errors.append(stop_result.error)
            if not frame_valid:
                result.errors.append(f"CDS length {cds_length} not divisible by 3")
        
        for splice in splice_results:
            if not splice.is_valid:
                result.errors.append(f"Intron {splice.intron_index}: {splice.error}")
        
        # Overall validity
        splice_ok = result.splice_sites_valid == result.splice_sites_checked
        
        if transcript.is_coding:
            result.is_valid = (
                splice_ok and
                result.start_codon_valid and
                result.stop_codon_valid and
                result.cds_length_valid
            )
        else:
            result.is_valid = splice_ok
        
        return result
    
    def validate_gene(self, gene: Gene) -> GeneValidation:
        """Validate all transcripts in a gene.
        
        Args:
            gene: The gene to validate
            
        Returns:
            GeneValidation result
        """
        result = GeneValidation(gene_id=gene.feature_id)
        
        for transcript in gene.transcripts:
            tx_result = self.validate_transcript(transcript)
            result.transcript_results.append(tx_result)
            result.transcripts_checked += 1
            if tx_result.is_valid:
                result.transcripts_valid += 1
        
        result.is_valid = result.transcripts_valid == result.transcripts_checked
        
        return result


def validate_mapped_genes(
    genes: Dict[str, Gene],
    fasta_handler: FastaHandler,
    progress_callback=None
) -> Dict[str, GeneValidation]:
    """Validate all mapped genes.
    
    Args:
        genes: Dictionary of gene_id -> Gene
        fasta_handler: Handler for target FASTA
        progress_callback: Optional progress callback
        
    Returns:
        Dictionary of gene_id -> GeneValidation
    """
    validator = StructuralValidator(fasta_handler)
    results = {}
    
    for i, (gene_id, gene) in enumerate(genes.items()):
        results[gene_id] = validator.validate_gene(gene)
        
        if progress_callback and (i + 1) % 100 == 0:
            progress_callback(i + 1, len(genes))
    
    valid_count = sum(1 for r in results.values() if r.is_valid)
    logger.info(f"Validated {len(genes)} genes: {valid_count} fully valid")
    
    return results

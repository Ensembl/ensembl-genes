"""FASTA file handler with indexing support."""

from collections import OrderedDict
from pathlib import Path
from typing import Dict, Optional, Union

import pysam
from Bio.Seq import Seq

from .models import Strand


class FastaHandler:
    """Handler for indexed FASTA files using pysam."""
    
    def __init__(self, fasta_path: Union[str, Path], cache_size: int = 4096):
        """Initialize FASTA handler.
        
        Args:
            fasta_path: Path to FASTA file (will create .fai index if needed)
        """
        self.fasta_path = Path(fasta_path)
        self._fasta: Optional[pysam.FastaFile] = None
        self._sequence_lengths: Dict[str, int] = {}
        self.cache_size = max(0, cache_size)
        self._fetch_cache: "OrderedDict[tuple, str]" = OrderedDict()
    
    def open(self):
        """Open the FASTA file. Handles both plain and bgzipped (.fa.gz) FASTA."""
        if self._fasta is not None:
            return

        # Create index if needed
        # For bgzipped files, pysam.faidx() creates both .fai and .gzi indices
        index_path = Path(str(self.fasta_path) + ".fai")
        if not index_path.exists():
            pysam.faidx(str(self.fasta_path))

        self._fasta = pysam.FastaFile(str(self.fasta_path))
        
        # Cache sequence lengths
        for name, length in zip(self._fasta.references, self._fasta.lengths):
            self._sequence_lengths[name] = length
    
    def close(self):
        """Close the FASTA file."""
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None
        self._fetch_cache.clear()
    
    def __enter__(self) -> "FastaHandler":
        self.open()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    @property
    def sequences(self) -> list:
        """Get list of sequence names."""
        if self._fasta is None:
            self.open()
        return list(self._fasta.references)
    
    def get_length(self, seq_name: str) -> int:
        """Get length of a sequence.
        
        Args:
            seq_name: Sequence/chromosome name
            
        Returns:
            Length in bases
        """
        if self._fasta is None:
            self.open()
        return self._sequence_lengths.get(seq_name, 0)
    
    def fetch(
        self,
        seq_name: str,
        start: int,
        end: int,
        strand: Strand = Strand.PLUS
    ) -> str:
        """Fetch a sequence region.
        
        Args:
            seq_name: Sequence/chromosome name
            start: Start position (1-based, inclusive)
            end: End position (1-based, inclusive)
            strand: Strand (will reverse complement if MINUS)
            
        Returns:
            Sequence string
        """
        if self._fasta is None:
            self.open()

        cache_key = (seq_name, start, end, str(strand))
        if self.cache_size > 0 and cache_key in self._fetch_cache:
            cached = self._fetch_cache.pop(cache_key)
            self._fetch_cache[cache_key] = cached
            return cached
        
        # pysam uses 0-based, half-open coordinates
        seq = self._fasta.fetch(seq_name, start - 1, end)
        seq = seq.upper()
        
        if strand == Strand.MINUS:
            seq = str(Seq(seq).reverse_complement())

        if self.cache_size > 0:
            self._fetch_cache[cache_key] = seq
            if len(self._fetch_cache) > self.cache_size:
                self._fetch_cache.popitem(last=False)
        
        return seq
    
    def get_codon(
        self,
        seq_name: str,
        pos: int,
        strand: Strand = Strand.PLUS
    ) -> str:
        """Get a codon at a specific position.
        
        Args:
            seq_name: Sequence/chromosome name
            pos: Position of first base of codon (1-based)
            strand: Strand
            
        Returns:
            3-base codon sequence
        """
        if strand == Strand.MINUS:
            return self.fetch(seq_name, pos - 2, pos, strand)
        else:
            return self.fetch(seq_name, pos, pos + 2, strand)
    
    def get_splice_site_dinucleotide(
        self,
        seq_name: str,
        pos: int,
        is_donor: bool,
        strand: Strand = Strand.PLUS
    ) -> str:
        """Get splice site dinucleotide.
        
        Args:
            seq_name: Sequence/chromosome name
            pos: Position at intron/exon boundary (1-based)
            is_donor: True for donor (5'), False for acceptor (3')
            strand: Strand
            
        Returns:
            2-base splice site sequence (e.g., "GT", "AG")
        """
        if strand == Strand.PLUS:
            if is_donor:
                # Donor: first 2 bases of intron (after exon end)
                return self.fetch(seq_name, pos + 1, pos + 2, strand)
            else:
                # Acceptor: last 2 bases of intron (before exon start)
                return self.fetch(seq_name, pos - 2, pos - 1, strand)
        else:
            # Minus strand: positions are mirrored
            if is_donor:
                return self.fetch(seq_name, pos - 2, pos - 1, strand)
            else:
                return self.fetch(seq_name, pos + 1, pos + 2, strand)
    
    def translate_cds(self, cds_sequences: list, phase: int = 0) -> str:
        """Translate concatenated CDS sequences.
        
        Args:
            cds_sequences: List of CDS sequences in transcript order
            phase: Starting phase (0, 1, or 2)
            
        Returns:
            Protein sequence (amino acids)
        """
        full_seq = "".join(cds_sequences)
        
        # Adjust for phase
        if phase > 0:
            full_seq = full_seq[phase:]
        
        # Translate
        try:
            protein = str(Seq(full_seq).translate(to_stop=False))
            return protein
        except Exception:
            return ""

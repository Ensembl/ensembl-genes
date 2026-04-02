#!/usr/bin/env python3
"""
Comprehensive validation script for Pangenome Mapping.

Evaluates mapping quality by subsampling genes and performing detailed 
pairwise analysis of transcripts and proteins.

Features:
- Parallel processing
- Fast test mode (GFF subsetting)
- Detailed structural conservation metrics
- Formatted console output
"""

import argparse
import json
import logging
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from functools import partial
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Set

import Bio.Seq
from Bio import Align

# Add parent directory to path to import src modules
sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.models import Gene, Transcript, Strand, CDS
from src.gff3_parser import GFF3Parser
from src.fasta_handler import FastaHandler
from src.config import GFF3_TRANSCRIPT_TYPES

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Maximum sequence length for full alignment (longer sequences use k-mer fallback)
# Biopython PairwiseAligner is O(n*m), so 50kb*50kb = 2.5 billion operations
MAX_ALIGNABLE_LENGTH = 50000  # 50kb
COPY_SUFFIX_RE = re.compile(r"^(.*?)(?:_copy(\d+))$")


def split_copy_suffix(value: Optional[str]) -> Tuple[str, Optional[int]]:
    """Split terminal _copyN suffix from an identifier."""
    if not value:
        return "", None
    match = COPY_SUFFIX_RE.match(value)
    if not match:
        return value, None
    base = match.group(1)
    copy_idx = int(match.group(2)) if match.group(2) is not None else None
    return base, copy_idx


def strip_copy_suffix(value: Optional[str]) -> str:
    """Drop terminal _copyN suffix if present."""
    base, _ = split_copy_suffix(value)
    return base


def parse_copy_index(value: Optional[str]) -> Optional[int]:
    """Extract terminal copy index from _copyN suffix."""
    _, copy_idx = split_copy_suffix(value)
    return copy_idx


def normalize_gene_id(gid: Optional[str]) -> str:
    """Normalize gene-style identifiers for cross-file matching."""
    if not gid:
        return ""
    value = gid.strip()
    if value.startswith("Name="):
        value = value[5:]
    if value.startswith("mapped_"):
        value = value[len("mapped_"):]
    if value.startswith("mapped_gene:"):
        value = value[len("mapped_gene:"):]
    if value.startswith("gene:"):
        value = value[len("gene:"):]
    return value


def normalize_transcript_id(txid: Optional[str]) -> str:
    """Normalize transcript-style IDs for robust matching."""
    if not txid:
        return ""
    value = txid.strip()
    if value.startswith("mapped_"):
        value = value[len("mapped_"):]
    if value.startswith("mapped_transcript:"):
        value = value[len("mapped_transcript:"):]
    if value.startswith("transcript:"):
        value = value[len("transcript:"):]
    return value


def extract_gff_attr(line: str, key: str) -> Optional[str]:
    """Extract one attribute value from a raw GFF3 line."""
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 9:
        return None
    attrs = parts[8]
    token = f"{key}="
    if token not in attrs:
        return None
    return attrs.split(token, 1)[1].split(";", 1)[0]


def find_gene_feature_line(lines: List[str]) -> Optional[str]:
    """Return the top-level gene feature line from grouped GFF lines."""
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        if parts[2] in ("gene", "pseudogene", "ncRNA_gene"):
            return line
    return None


def fast_kmer_identity(seq1: str, seq2: str, k: int = 21) -> Tuple[float, float]:
    """
    Fast approximate identity using k-mer overlap.
    Returns (identity%, coverage%) - much faster than alignment for long sequences.
    """
    if not seq1 or not seq2:
        return 0.0, 0.0
    
    # Build k-mer sets
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1)) if len(seq1) >= k else {seq1}
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1)) if len(seq2) >= k else {seq2}
    
    if not kmers1:
        return 0.0, 0.0
    
    # Jaccard-like similarity
    # Note: k-mer method cannot distinguish identity from coverage without full alignment
    # Both metrics represent the fraction of reference k-mers found in target
    intersection = len(kmers1 & kmers2)
    identity = (intersection / len(kmers1)) * 100  # Approx: percent of ref k-mers found in target
    coverage = identity  # Same for k-mer approximation
    
    return identity, coverage


# --- Helper Functions (Picklable for Multiprocessing) ---

def get_sequence(transcript: Transcript, handler: FastaHandler) -> str:
    """Construct spliced transcript sequence."""
    sorted_exons = sorted(transcript.exons, key=lambda e: e.start)
    seq_parts = []
    for exon in sorted_exons:
        # fetch returns + strand
        seq = handler.fetch(transcript.seq_region, exon.start, exon.end)
        seq_parts.append(seq)
    full_seq = "".join(seq_parts)
    if transcript.strand == Strand.MINUS:
        return str(Bio.Seq.Seq(full_seq).reverse_complement())
    return full_seq

def get_cds_sequence(transcript: Transcript, handler: FastaHandler) -> str:
    sorted_cds = sorted(transcript.cds_list, key=lambda c: c.start)
    seq_parts = []
    for cds in sorted_cds:
        seq = handler.fetch(transcript.seq_region, cds.start, cds.end)
        seq_parts.append(seq)
    full_seq = "".join(seq_parts)
    if transcript.strand == Strand.MINUS:
        return str(Bio.Seq.Seq(full_seq).reverse_complement())
    return full_seq

def get_splice_windows(tx: Transcript, handler: FastaHandler, window: int = 10) -> List[Tuple[str, str]]:
    """Get list of (donor_window, acceptor_window) for all introns.
    Window is +/- 'window' bp centered on the splice junction.
    (10bp exon + 10bp intron).
    """
    windows = []
    sorted_exons = sorted(tx.exons, key=lambda e: e.start)
    if len(sorted_exons) < 2:
        return []
    
    for i in range(len(sorted_exons) - 1):
        e1 = sorted_exons[i]
        e2 = sorted_exons[i+1]
        intron_start = e1.end + 1
        intron_end = e2.start - 1
        if intron_end < intron_start: continue 
        
        # Determine Genomic Coordinates for Donor/Acceptor Junctions
        # Note: pysam fetch(start, end) is 0-based half-open [start, end)
        
        if tx.strand == Strand.MINUS:
            # Donor is at Intron Genomic End (Transcript 5')
            # Junction: Intron End | Exon Start (Genomically) -> but Transcript wise Exon-Intron.
            # Intron End `E`. Exon starts at `E+1`.
            # We want 10bp Exon (E+1...E+10), 10bp Intron (E-9...E).
            # RC of this gives Transcript 5'->3' order (Exon...Intron).
            
            # Fetch genomic range [E - window + 1 , E + window + 1] ?
            # Let's say window=10.
            # Intron 3' (Genomic End): index `intron_end - 1`.
            # We want indices `intron_end - 10` to `intron_end + 9` (20 bases).
            # `intron_end` (1-based) E.
            # `fetch(E - 10, E + 10)`. 
            #   -> [E-10, E-1] (10bp Intron end) + [E, E+9] (10bp Exon start).
            # Wait, Donor (Transcript 5') on Minus strand should be Exon->Intron.
            # Genomic: Exon [High] -> Intron [Low].
            # So Genomic Sequence [Low..High] is Intron...Exon.
            # RC(Intron...Exon) = RC(Exon)...RC(Intron). Correct.
            # So we fetch centered on Intron End.
            
            donor_genomic = handler.fetch(tx.seq_region, intron_end - window, intron_end + window)
            donor_seq = str(Bio.Seq.Seq(donor_genomic).reverse_complement())
            
            # Acceptor is at Intron Genomic Start (Transcript 3')
            # Genomic: Intron [Low..] <- Exon [..Low-1]
            # Junction at `intron_start`.
            # Genomic Window centered on `intron_start`.
            # [intron_start-10 ... intron_start+9]
            # Genomic: Exon...Intron.
            # RC: RC(Intron)...RC(Exon). (Intron->Exon). Correct.
            
            # intron_start S. 0-based S-1.
            # fetch(S-1 - 10, S-1 + 10).
            
            acceptor_genomic = handler.fetch(tx.seq_region, intron_start - 1 - window, intron_start - 1 + window)
            acceptor_seq = str(Bio.Seq.Seq(acceptor_genomic).reverse_complement())
            
        else:
            # Plus Strand
            # Donor: Intron Start. Exon->Intron.
            # Junction at `intron_start`.
            # Exon [..S-1] -> Intron [S..].
            # Window centered on S.
            # [S-11 ... S-2] (Exon) + [S-1 ... S+8] (Intron).
            # fetch(S-1 - 10, S-1 + 10).
            donor_genomic = handler.fetch(tx.seq_region, intron_start - 1 - window, intron_start - 1 + window)
            donor_seq = donor_genomic
            
            # Acceptor: Intron End. Intron->Exon.
            # Junction at `intron_end` (last base) | Exon start.
            # Or `intron_end + 1` is Exon start.
            # We want Intron...Exon.
            # Window centered on Exon Start (`intron_end + 1`).
            # Or `intron_end`?
            # E is last intron base. E+1 is first exon base.
            # fetch(E - 10, E + 10).
            # [E-10 ... E-1] (Intron). [E ... E+9] (Exon). Correct.
            
            acceptor_genomic = handler.fetch(tx.seq_region, intron_end - window, intron_end + window)
            acceptor_seq = acceptor_genomic
            
        windows.append((donor_seq, acceptor_seq))
    return windows

def check_seq_identity(s1: str, s2: str) -> float:
    """Calculate identity between two strings."""
    if not s1 or not s2: return 0.0
    # Lengths should ideally match for window comparison
    l = min(len(s1), len(s2))
    if l == 0: return 0.0
    matches = sum(1 for i in range(l) if s1[i] == s2[i])
    return matches / l

def worker_evaluate_pair(ref_gene, target_gene, ref_handler, target_handler, splice_window, splice_identity_threshold):
    """
    Worker function to evaluate a single gene pair.
    Args:
        ref_gene, target_gene: Gene objects
        ref_handler, target_handler: Pre-initialized FastaHandler instances
        splice_window, splice_identity_threshold: Parameters
    Returns:
        (metrics_dict, list_of_alignment_strings)
    """
    
    # Select representative transcript
    def get_representative(gene):
        if not gene.transcripts: return None
        return sorted(gene.transcripts, 
                     key=lambda t: (t.is_coding, t.cds_length, t.length), 
                     reverse=True)[0]
    
    ref_tx = get_representative(ref_gene)
    normalized_ref_gene_id = normalize_gene_id(ref_gene.feature_id)
    if not ref_tx:
        return ({
            "gene_id": normalized_ref_gene_id,
            "status": "unmapped_no_transcripts"
        }, [])

    # Find matching target transcript using stable identifiers.
    ref_tx_id = ref_tx.feature_id
    ref_tx_norm = normalize_transcript_id(ref_tx_id)
    ref_tx_base = strip_copy_suffix(ref_tx_norm)
    ref_tx_candidates = {
        ref_tx_id,
        ref_tx_norm,
        ref_tx_base,
        f"transcript:{ref_tx_norm}",
        f"mapped_{ref_tx_id}",
        f"mapped_transcript:{ref_tx_norm}",
    }

    target_tx = None
    target_match_strategy = "id_or_mapped_from"
    best_score = None
    for tx in target_gene.transcripts:
        tx_id = tx.feature_id
        tx_norm = normalize_transcript_id(tx_id)
        tx_base = strip_copy_suffix(tx_norm)
        mapped_from = (tx.attributes or {}).get("mapped_from", "")
        mapped_from_norm = normalize_transcript_id(mapped_from)
        mapped_from_base = strip_copy_suffix(mapped_from_norm)

        rank = None
        # Highest confidence: explicit mapped_from provenance.
        if mapped_from and (
            mapped_from == ref_tx_id
            or mapped_from_norm == ref_tx_norm
            or mapped_from_base == ref_tx_base
        ):
            rank = 0
        # Next best: exact/normalized transcript ID match.
        elif tx_id in ref_tx_candidates or tx_norm == ref_tx_norm or tx_base == ref_tx_base:
            rank = 1
        # Fallback: mapped_ forms after normalization.
        elif f"mapped_transcript:{tx_norm}" in ref_tx_candidates:
            rank = 2

        if rank is None:
            continue

        score = (
            rank,
            0 if tx.is_coding == ref_tx.is_coding else 1,
            abs(len(tx.exons) - len(ref_tx.exons)),
            abs(tx.length - ref_tx.length),
            tx_id,
        )
        if best_score is None or score < best_score:
            best_score = score
            target_tx = tx
            if rank == 0:
                break
    
    # Debug: Print first few mismatches
    # if not target_tx:
    #     print(f"TX MISMATCH: ref={ref_tx_id} tried={candidate_ids} target_tx_ids={[t.feature_id for t in target_gene.transcripts[:3]]}")
    
    if not target_tx:
        fallback_tx = get_representative(target_gene)
        if fallback_tx is not None:
            target_tx = fallback_tx
            target_match_strategy = "representative_fallback"
        else:
            return ({
                "gene_id": normalized_ref_gene_id,
                "transcript_id": ref_tx.feature_id,
                "biotype": ref_gene.biotype,
                "status": "unmapped_transcript_missing"
            }, [])

    # Setup aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    metrics = {
        "gene_id": normalized_ref_gene_id,
        "transcript_id": ref_tx.feature_id,
        "biotype": ref_gene.biotype,
        "status": "mapped",
        "target_transcript_match": target_match_strategy,
    }
    alignments_out = []

    # 1. Transcript Alignment
    ref_seq = get_sequence(ref_tx, ref_handler)
    tgt_seq = get_sequence(target_tx, target_handler)
    
    metrics["ref_length"] = len(ref_seq)
    metrics["target_length"] = len(tgt_seq)
    metrics["length_diff_pct"] = (len(tgt_seq) - len(ref_seq)) / len(ref_seq) * 100 if len(ref_seq) > 0 else 0

    # Alignment
    if not ref_seq or not tgt_seq:
        metrics["transcript_identity"] = 0.0
        metrics["transcript_coverage"] = 0.0
        metrics["alignment_method"] = "none"
        status_detail = "empty_seq_ref" if not ref_seq else "empty_seq_tgt"
        metrics["status"] = f"unmapped_{status_detail}"
        alignments_out.append(f"TRANSCRIPT {ref_tx.feature_id}: EMPTY SEQUENCE ({status_detail})")
    elif len(ref_seq) > MAX_ALIGNABLE_LENGTH or len(tgt_seq) > MAX_ALIGNABLE_LENGTH:
        # Use fast k-mer fallback for very long sequences
        identity, coverage = fast_kmer_identity(ref_seq, tgt_seq)
        metrics["transcript_identity"] = identity
        metrics["transcript_coverage"] = coverage
        metrics["alignment_method"] = "kmer"
        alignments_out.append(f"TRANSCRIPT {ref_tx.feature_id}: KMER FALLBACK (seq too long: {len(ref_seq)}/{len(tgt_seq)}bp)")
    else:
        alignment = aligner.align(ref_seq, tgt_seq)[0]
        
        # Identity: percentage of aligned positions that are identical
        identity = (alignment.counts().identities / alignment.length) * 100 if alignment.length > 0 else 0
        
        # Coverage: percentage of reference bases that are aligned (not gapped)
        # alignment.aligned[0] contains the ref coordinate ranges that are aligned
        ref_aligned_bases = sum(end - start for start, end in alignment.aligned[0])
        coverage = (ref_aligned_bases / len(ref_seq)) * 100 if len(ref_seq) > 0 else 0
        
        metrics["transcript_identity"] = identity
        metrics["transcript_coverage"] = coverage
        metrics["alignment_method"] = "full"
        alignments_out.append(format_alignment("TRANSCRIPT", ref_tx.feature_id, alignment, identity, coverage))

    # 2. Structural Conservation metrics
    
    # Exon Count
    metrics["ref_exon_count"] = len(ref_tx.exons)
    metrics["target_exon_count"] = len(target_tx.exons)
    metrics["exon_count_conserved"] = (len(ref_tx.exons) == len(target_tx.exons))

    # 3. Protein Alignment (if coding)
    if ref_tx.is_coding and target_tx.is_coding: # Only if both are coding
        ref_cds = get_cds_sequence(ref_tx, ref_handler)
        tgt_cds = get_cds_sequence(target_tx, target_handler)
        
        metrics["ref_cds_length"] = len(ref_cds)
        metrics["target_cds_length"] = len(tgt_cds)
        metrics["cds_length_conserved"] = (len(ref_cds) == len(tgt_cds))
        
        # Translation
        ref_prot = str(Bio.Seq.Seq(ref_cds).translate()).strip("*")
        tgt_prot = str(Bio.Seq.Seq(tgt_cds).translate()).strip("*")
        
        p_alignment = aligner.align(ref_prot, tgt_prot)[0]
        p_identity = (p_alignment.counts().identities / p_alignment.length) * 100 if p_alignment.length > 0 else 0
        p_coverage = (p_alignment.counts().identities / len(ref_prot)) * 100 if len(ref_prot) > 0 else 0
        
        metrics["protein_identity"] = p_identity
        metrics["protein_coverage"] = p_coverage
        alignments_out.append(format_alignment("PROTEIN", ref_tx.feature_id, p_alignment, p_identity, p_coverage))
        
        metrics["valid_start_codon_ref"] = ref_cds.startswith("ATG")
        metrics["valid_start_codon_target"] = tgt_cds.startswith("ATG")
        metrics["valid_stop_codon_ref"] = ref_cds.endswith(("TAA", "TAG", "TGA"))
        metrics["valid_stop_codon_target"] = tgt_cds.endswith(("TAA", "TAG", "TGA"))

    # 4. Splice Sites (Alignment Based)
    # Strategy: Map REF exon junctions to TARGET coordinates via alignment
    
    # 1. Get Junction Coords (Transcript Space)
    def get_junctions(tx):
        # Sort exons: Ascending for PLUS, Descending for MINUS (Transcript 5'->3')
        if tx.strand == Strand.MINUS:
            sorted_exons = sorted(tx.exons, key=lambda e: e.start, reverse=True)
        else:
            sorted_exons = sorted(tx.exons, key=lambda e: e.start)
            
        junctions = set()
        cum_len = 0
        for i, exon in enumerate(sorted_exons):
            cum_len += (exon.end - exon.start + 1)
            # Add junction if not last exon
            if i < len(sorted_exons) - 1:
                junctions.add(cum_len)
        return junctions
        
    ref_junctions = get_junctions(ref_tx)
    tgt_junctions = get_junctions(target_tx)
    
    metrics["ref_intron_count"] = len(ref_junctions)
    metrics["target_intron_count"] = len(tgt_junctions)
    
    # Only do alignment-based splice site mapping if we have a full alignment
    if metrics.get("alignment_method") == "full":
        # 2. Map Ref Junctions to Target via Alignment
        aligned_ref = alignment[0] # Sequence string with gaps
        aligned_tgt = alignment[1]
        
        conserved_count = 0
        
        # Tracking
        curr_ref_idx = 0
        curr_tgt_idx = 0
        
        # We need to map Ref Coordinate `R` to Target Coordinate `T`.
        # A Ref Junction `R` means "After base R".
        # Find the alignment column corresponding to Ref base `R`.
        # Map to Target base `T` at that column.
        # Check if `T` is in `tgt_junctions`.
        # Handling GAPs:
        # If Ref base `R` aligns to Gap in Target?
        #   Ref: ...A (R)
        #   Tgt: ...-
        #   The Ref base exists but has no Target counterpart. 
        #   Biologically, this segment is deleted in Target.
        #   The "Junction" at R is likely lost or shifted.
        #   Tolerant Check: If the Target has a junction at the *last non-gap base*?
        #   Tgt: ...G - (Gap at R)
        #   If Target Junction is after G?
        #   Then "After G" aligns to "After A".
        #   Wait, if Ref has A, Target has Gap.
        #   Ref: ...prev A | next...
        #   Tgt: ...prev - | next...
        #   Splice is between A and next.
        #   Splice in Target is between prev and next.
        #   If Ref(A) aligns to Gap.
        #   Ref(next) aligns to Tgt(next).
        #   If Tgt(prev) is Tgt's junction?
        #   Then yes, splice is "conserved".
        
        # Implementation:
        # Scan alignment column by column.
        
        # We want to check `curr_ref_idx` against `ref_junctions`.
        # Check triggers when `curr_ref_idx` EQUALS a junction coordinate.
        # But `curr_ref_idx` increments.
        # If R=10.
        # We reach col `c` where ref base is at position 10.
        # `curr_ref_idx` AFTER increment is 10.
        # So check happens just after we reach base 10.
        
        for i in range(len(aligned_ref)):
            r_char = aligned_ref[i]
            t_char = aligned_tgt[i]
            
            is_ref_base = (r_char != '-')
            if is_ref_base:
                curr_ref_idx += 1
            
            if t_char != '-':
                curr_tgt_idx += 1
                
            if is_ref_base and (curr_ref_idx in ref_junctions):
                # Case 1: Target has base here. `t_char != '-'`.
                #   We are at `curr_tgt_idx`.
                #   Is `curr_tgt_idx` a Target Junction?
                #   If yes -> Conserved.
                
                # Case 2: Target has Gap here. `t_char == '-'`.
                #   Ref base `R` maps to Gap.
                #   Ref: ... (R) ...
                #   Tgt: ... (-) ...
                #   We are at `curr_tgt_idx` (last non-gap target base).
                #   Is `curr_tgt_idx` a Target Junction?
                #   If yes -> Conserved (Deletion at end of exon).
                
                if curr_tgt_idx in tgt_junctions:
                    conserved_count += 1
                    
            # Handle GAP in Ref (Target Insertion)
            # Ref: ... (prev) - (next) ...
            # Tgt: ... (prev) I (next) ...
            # If Ref Junction was at `prev`.
            # We checked it when processing `prev`.
            # At that point, `curr_tgt_idx` was `prev_tgt`.
            # If Target has Insertion `I` *before* the junction?
            # i.e. Exon extension.
            # Ref: ...A | T...
            # Tgt: ...A G | T... (Insertion G).
            # Alignment:
            # Ref: A - T
            # Tgt: A G T
            # Col 1: Ref(A). Tgt(A). Ref=R. Ref is Junction.
            #   Check: Tgt=T. Is T a Junction? 
            #   No, Target junction is T+1 (after G).
            #   So strict logic fails. 
            #   BUT user wants "conserved around positions".
            #   If gap follows?
            #   Logic: "If we are at Ref Junction, check if Target Junction is `curr_tgt_idx` OR reachable via local gaps?"
            #   This requires lookahead.
            
            # Simpler: Allow a tolerance window?
            # User asked for "Check pairwise alignment... bases aligned...".
            # The strict check handles Deletions (Gap in Target) correctly IF the alignment places the gap *aligned* to the Ref base.
            # For Insertions (Gap in Ref), the check happens *before* the gap.
            # If we see Ref Junction, and Target matches `curr_tgt_idx` but Target Junction is `curr_tgt_idx + N`.
            # And the NEXT N columns are Gaps in Ref (Target Insertions).
            # Then Target Junction comes N bases later (in Target).
            # The Ref Junction aligns to a position *before* the Target Junction.
            # This is a shift, not a loss.
            
            # Tolerant: Lookahead.
            # If Ref(Junction at R). Check tgt_junctions for curr_tgt_idx.
            # If NOT found but if `aligned_ref[i+1]` is Gap?
            #   Scan ahead while ref is Gap.
            #   Check each `curr_tgt_idx + delta` for junction membership.
            pass 
        
        # Re-run loop with lookahead for insertions
        conserved_count = 0
        curr_ref_idx = 0
        curr_tgt_idx = 0
        
        for i in range(len(aligned_ref)):
            r_char = aligned_ref[i]
            t_char = aligned_tgt[i]
            
            is_ref_base = (r_char != '-')
            if is_ref_base:
                curr_ref_idx += 1
            
            if t_char != '-':
                curr_tgt_idx += 1
                
            if is_ref_base and (curr_ref_idx in ref_junctions):
                # Check 1: Exact match or Deletion-handling (Target Junction is HERE)
                if curr_tgt_idx in tgt_junctions:
                    conserved_count += 1
                    continue
                    
                # Check 2: Insertion-handling (Target Junction is AHEAD in gap block)
                # Scan ahead in alignment
                # While Ref is Gap, if we hit Target Junction?
                temp_tgt = curr_tgt_idx
                found_ahead = False
                for k in range(i + 1, len(aligned_ref)):
                    r_next = aligned_ref[k]
                    t_next = aligned_tgt[k]
                    
                    if r_next != '-':
                        # End of Ref Gap block (or no gap). Ref advances.
                        # Junction check window closed.
                        break
                    
                    if t_next != '-':
                        temp_tgt += 1
                        if temp_tgt in tgt_junctions:
                            found_ahead = True
                            break
                
                if found_ahead:
                    conserved_count += 1
                
        metrics["introns_conserved_count"] = conserved_count
    else:
        # K-mer fallback: Use simpler exon count comparison
        # If exon counts match, assume splice sites are conserved
        if len(ref_junctions) == len(tgt_junctions):
            metrics["introns_conserved_count"] = len(ref_junctions)
        else:
            # Pro-rate based on minimum count (pessimistic estimate)
            metrics["introns_conserved_count"] = min(len(ref_junctions), len(tgt_junctions))
    
    return metrics, alignments_out


def worker_batch_process(args_pack):
    """
    Process a batch of genes.
    args_pack: (batch_list, ref_fasta, tgt_fasta, split_win, splice_id, batch_idx)
    batch_list: list of tuples (ref_id, ref_lines, tgt_id, tgt_lines)
    """
    batch_list, ref_fasta, target_fasta, splice_w, splice_id, batch_idx = args_pack
    
    # Init FASTA handlers ONCE per batch for performance
    ref_handler = FastaHandler(ref_fasta)
    target_handler = FastaHandler(target_fasta)
    
    # Init parsers
    ref_parser = GFF3Parser()
    tgt_parser = GFF3Parser()
    
    results = []
    alignments = []
    
    for i, (ref_id, ref_lines, tgt_id, tgt_lines) in enumerate(batch_list):
        # Parse on demand
        try:
            ref_genes = ref_parser.parse_lines(iter(ref_lines))
            target_genes = tgt_parser.parse_lines(iter(tgt_lines))
        except Exception as e:
            print(f"Error parsing {ref_id}: {e}")
            continue
        
        if not ref_genes or not target_genes:
            continue
            
        # Get the genes
        ref_gene = next(iter(ref_genes.values()))
        target_gene = next(iter(target_genes.values()))
        
        if not ref_gene.transcripts:
            print(f"WARNING: Ref Gene {ref_id} has NO transcripts. Lines: {len(ref_lines)}")
            
        if not target_gene.transcripts:
             print(f"WARNING: Target Gene {tgt_id} has NO transcripts. Lines: {len(tgt_lines)}")
        
        # Call evaluator with cached handlers
        metric, aln = worker_evaluate_pair(
            ref_gene, target_gene, ref_handler, target_handler, splice_w, splice_id
        )
        results.append(metric)
        alignments.extend(aln)
    
    # Progress indicator
    print(f"Batch {batch_idx} complete: {len(results)} genes processed", flush=True)
        
    return results, alignments


# --- Synteny/Gene Order Conservation ---

class GeneNeighborIndex:
    """Pre-compute sorted gene lists per chromosome for efficient neighbor lookup."""
    
    def __init__(self, gff_lines_map: Dict[str, List[str]], biotype_fn):
        """
        Build index from GFF lines map.
        Args:
            gff_lines_map: {gene_id: [gff_lines]}
            biotype_fn: Function to determine biotype group from biotype string
        """
        self.genes = {}  # gene_id -> {chrom, start, end, strand, biotype, biotype_group}
        self.by_chrom = defaultdict(list)  # chrom -> sorted list of (start, end, gene_id)
        self.biotype_fn = biotype_fn
        
        for gene_id, lines in gff_lines_map.items():
            gene_info = self._parse_gene_info(gene_id, lines)
            if gene_info:
                self.genes[gene_id] = gene_info
                self.by_chrom[gene_info["chrom"]].append(
                    (gene_info["start"], gene_info["end"], gene_id)
                )
        
        # Sort each chromosome's genes by start position
        for chrom in self.by_chrom:
            self.by_chrom[chrom].sort(key=lambda x: x[0])
    
    def _parse_gene_info(self, gene_id: str, lines: List[str]) -> Optional[Dict]:
        """Parse minimal gene info from GFF lines."""
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            if feature_type in ("gene", "pseudogene", "ncRNA_gene"):
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                
                # Extract biotype
                attrs = parts[8]
                biotype = "unknown"
                for key in ["gene_biotype=", "biotype=", "gene_type="]:
                    if key in attrs:
                        biotype = attrs.split(key)[1].split(";")[0]
                        break
                
                return {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "biotype": biotype,
                    "biotype_group": self.biotype_fn(biotype)
                }
        return None
    
    def get_neighbors(self, gene_id: str, n: int = 1, protein_coding_only: bool = False) -> Tuple[List[str], List[str]]:
        """
        Get n nearest non-overlapping neighbors on each side.
        Returns: (list of 5' neighbor gene_ids, list of 3' neighbor gene_ids)
        5' = lower coordinates, 3' = higher coordinates (ignoring strand)
        """
        if gene_id not in self.genes:
            return [], []
        
        gene = self.genes[gene_id]
        chrom = gene["chrom"]
        gene_start = gene["start"]
        gene_end = gene["end"]
        
        chrom_genes = self.by_chrom.get(chrom, [])
        if not chrom_genes:
            return [], []
        
        # Find position of this gene in sorted list
        idx = None
        for i, (start, end, gid) in enumerate(chrom_genes):
            if gid == gene_id:
                idx = i
                break
        
        if idx is None:
            return [], []
        
        # Collect 5' neighbors (lower coordinates, going backwards)
        neighbors_5p = []
        for i in range(idx - 1, -1, -1):
            start, end, gid = chrom_genes[i]
            # Check non-overlapping
            if end < gene_start:
                if protein_coding_only:
                    if self.genes[gid]["biotype_group"] == "Protein Coding":
                        neighbors_5p.append(gid)
                else:
                    neighbors_5p.append(gid)
                if len(neighbors_5p) >= n:
                    break
        
        # Collect 3' neighbors (higher coordinates, going forwards)
        neighbors_3p = []
        for i in range(idx + 1, len(chrom_genes)):
            start, end, gid = chrom_genes[i]
            # Check non-overlapping
            if start > gene_end:
                if protein_coding_only:
                    if self.genes[gid]["biotype_group"] == "Protein Coding":
                        neighbors_3p.append(gid)
                else:
                    neighbors_3p.append(gid)
                if len(neighbors_3p) >= n:
                    break
        
        return neighbors_5p, neighbors_3p


def compute_synteny_metrics(
    ref_index: GeneNeighborIndex,
    target_index: GeneNeighborIndex,
    gene_ids: List[str],
    id_mapping: Dict[str, str]  # ref_gene_id -> target_gene_id
) -> List[Dict]:
    """
    Compute synteny metrics for a list of genes.
    Returns list of dicts with synteny scores per gene.
    """
    results = []
    
    for ref_gid in gene_ids:
        target_gid = id_mapping.get(ref_gid)
        
        if not target_gid or target_gid not in target_index.genes:
            # Gene not mapped - can't compute synteny
            results.append({
                "gene_id": ref_gid,
                "biotype_group": ref_index.genes.get(ref_gid, {}).get("biotype_group", "Other"),
                "synteny_status": "unmapped"
            })
            continue
        
        ref_info = ref_index.genes[ref_gid]
        metrics = {
            "gene_id": ref_gid,
            "biotype_group": ref_info["biotype_group"],
            "synteny_status": "mapped"
        }
        
        # Method 1: Strict Protein-Coding Neighbors
        ref_5p_pc, ref_3p_pc = ref_index.get_neighbors(ref_gid, n=1, protein_coding_only=True)
        tgt_5p_pc, tgt_3p_pc = target_index.get_neighbors(target_gid, n=1, protein_coding_only=True)
        
        # Map ref neighbor IDs to target
        ref_5p_pc_mapped = [id_mapping.get(g) for g in ref_5p_pc if id_mapping.get(g)]
        ref_3p_pc_mapped = [id_mapping.get(g) for g in ref_3p_pc if id_mapping.get(g)]
        
        pc_5p_match = (ref_5p_pc_mapped == tgt_5p_pc) if ref_5p_pc_mapped or tgt_5p_pc else True
        pc_3p_match = (ref_3p_pc_mapped == tgt_3p_pc) if ref_3p_pc_mapped or tgt_3p_pc else True
        
        # Handle boundaries
        if not ref_5p_pc and not tgt_5p_pc:
            pc_conserved = pc_3p_match
            metrics["pc_neighbors_valid"] = 1 if ref_3p_pc or tgt_3p_pc else 0
        elif not ref_3p_pc and not tgt_3p_pc:
            pc_conserved = pc_5p_match
            metrics["pc_neighbors_valid"] = 1 if ref_5p_pc or tgt_5p_pc else 0
        else:
            pc_conserved = pc_5p_match and pc_3p_match
            metrics["pc_neighbors_valid"] = 1
        
        metrics["pc_neighbors_conserved"] = 1 if pc_conserved else 0
        
        # Method 2: Any-Biotype Neighbors
        ref_5p_any, ref_3p_any = ref_index.get_neighbors(ref_gid, n=1, protein_coding_only=False)
        tgt_5p_any, tgt_3p_any = target_index.get_neighbors(target_gid, n=1, protein_coding_only=False)
        
        ref_5p_any_mapped = [id_mapping.get(g) for g in ref_5p_any if id_mapping.get(g)]
        ref_3p_any_mapped = [id_mapping.get(g) for g in ref_3p_any if id_mapping.get(g)]
        
        any_5p_match = (ref_5p_any_mapped == tgt_5p_any) if ref_5p_any_mapped or tgt_5p_any else True
        any_3p_match = (ref_3p_any_mapped == tgt_3p_any) if ref_3p_any_mapped or tgt_3p_any else True
        
        if not ref_5p_any and not tgt_5p_any:
            any_conserved = any_3p_match
            metrics["any_neighbors_valid"] = 1 if ref_3p_any or tgt_3p_any else 0
        elif not ref_3p_any and not tgt_3p_any:
            any_conserved = any_5p_match
            metrics["any_neighbors_valid"] = 1 if ref_5p_any or tgt_5p_any else 0
        else:
            any_conserved = any_5p_match and any_3p_match
            metrics["any_neighbors_valid"] = 1
        
        metrics["any_neighbors_conserved"] = 1 if any_conserved else 0
        
        # Method 3: Neighborhood Score (5 genes per side)
        ref_5p_n5, ref_3p_n5 = ref_index.get_neighbors(ref_gid, n=5, protein_coding_only=False)
        tgt_5p_n5, tgt_3p_n5 = target_index.get_neighbors(target_gid, n=5, protein_coding_only=False)
        
        def score_side(ref_neighbors, tgt_neighbors, id_map):
            """Score one side based on neighbor conservation."""
            if not ref_neighbors:
                return None  # Boundary - no neighbors expected
            
            ref_mapped = [id_map.get(g) for g in ref_neighbors if id_map.get(g)]
            if not ref_mapped:
                return 0.0  # All expected neighbors unmapped
            
            # Check how many are present and in order
            present_count = sum(1 for g in ref_mapped if g in tgt_neighbors)
            
            if present_count == 0:
                return 0.0  # None of the expected genes present
            
            # Check order preservation
            ref_order = [g for g in ref_mapped if g in tgt_neighbors]
            tgt_order = [g for g in tgt_neighbors if g in ref_mapped]
            
            if ref_order == tgt_order:
                # In order - score based on completeness
                if present_count == len(ref_mapped):
                    return 1.0  # All present and in order
                elif present_count >= len(ref_mapped) - 1:
                    return 0.4  # Missing one but others in order
                else:
                    return 0.1  # Only some present
            else:
                return 0.0  # Present but out of order
        
        score_5p = score_side(ref_5p_n5, tgt_5p_n5, id_mapping)
        score_3p = score_side(ref_3p_n5, tgt_3p_n5, id_mapping)
        
        # Average, handling boundaries
        if score_5p is None and score_3p is None:
            metrics["neighborhood_score"] = None  # Isolated gene
        elif score_5p is None:
            metrics["neighborhood_score"] = score_3p
        elif score_3p is None:
            metrics["neighborhood_score"] = score_5p
        else:
            metrics["neighborhood_score"] = (score_5p + score_3p) / 2
        
        results.append(metrics)
    
    return results


class ValidationRunner:
    def __init__(self, args):
        self.args = args

    def determine_biotype_group(self, biotype: str) -> str:
        if not biotype: return "Other"
        b = biotype.lower()
        if "protein_coding" in b: return "Protein Coding"
        if "lnc" in b or "linc" in b: return "lncRNA"
        if "pseudogene" in b: return "Pseudogene"
        if "ig_" in b or "tr_" in b or "v_gene" in b or "j_gene" in b: return "IG/TR"
        if any(x in b for x in ["mirna", "snrna", "snorna", "rrna", "trna"]): return "Small RNA"
        return "Other"

    # ... fast_load_genes and load_subset_gff omitted for brevity if unchanged ...

    # We need to include them to keep the file valid, or use multi_replace.
    # But since this is ReplaceFileContent (single block), I must match exact target. 
    # The target block is huge. I will try to use the boundaries defined in the request.
    
    # Wait, the instruction said "Replace check_splice_motifs...".
    # I should target just the function if possible, but the user requested changes in worker_evaluate_pair too.
    # And print_summary_table.
    # This is a multi-part change. 
    # I will use multi_replace_file_content instead?
    # No, the previous tool call was replace_file_content.
    # I will construct the block to cover from `check_splice_motifs` definition down to `print_summary_table` start?
    # That's too big.
    
    # I'll use multi_replace_file_content to be safe and precise.
    def _process_orphan_group(self, parent_id, target_gene_id, orphan_lines, gene_lines, transcript_to_gene):
        """Recursively process orphans when a parent is resolved."""
        if parent_id in orphan_lines:
            items = orphan_lines.pop(parent_id)
            for item in items:
                # item is (line, my_id)
                line, my_id = item
                gene_lines[target_gene_id].append(line)
                
                if my_id:
                    transcript_to_gene[my_id] = target_gene_id
                    # Recursively resolve children of this item
                    self._process_orphan_group(my_id, target_gene_id, orphan_lines, gene_lines, transcript_to_gene)

    def load_gff_lines(self, filepath):
        """Load GFF lines grouped by ID, handling out-of-order features.
        Returns: {gene_id: [lines]}
        """
        logger.info(f"Indexing GFF: {filepath}")
        gene_lines = defaultdict(list)
        orphan_lines = defaultdict(list) # ParentID -> [(line, my_id)]
        transcript_to_gene = {} 
        
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split("\t")
                if len(parts) < 9: continue
                
                ftype = parts[2]
                attrs = parts[8]
                
                curr_id = None
                parent = None
                
                if "ID=" in attrs:
                    curr_id = attrs.split("ID=")[1].split(";")[0]

                if "Parent=" in attrs:
                    parent = attrs.split("Parent=")[1].split(";")[0]

                # Identify Genes efficiently
                # Some GFFs use "gene", "pseudogene", etc.
                is_gene = ftype in ["gene", "pseudogene", "ncRNA_gene"]
                
                if is_gene and curr_id:
                    store_id = curr_id
                    if store_id.startswith("gene:"): store_id = store_id[5:]
                    
                    gene_lines[store_id].append(line)
                    
                    # Resolve orphans waiting for this gene
                    self._process_orphan_group(curr_id, store_id, orphan_lines, gene_lines, transcript_to_gene)
                    
                elif parent:
                    # Try to find GeneID
                    target_gene_id = None
                    
                    # Check maps (GeneID match)
                    store_parent = parent
                    if store_parent.startswith("gene:"): store_parent = store_parent[5:]
                    
                    if store_parent in gene_lines: target_gene_id = store_parent
                    elif parent in gene_lines: target_gene_id = parent
                    elif parent in transcript_to_gene: target_gene_id = transcript_to_gene[parent]
                    
                    if target_gene_id:
                        gene_lines[target_gene_id].append(line)
                        if curr_id:
                            transcript_to_gene[curr_id] = target_gene_id
                            # Resolve children
                            self._process_orphan_group(curr_id, target_gene_id, orphan_lines, gene_lines, transcript_to_gene)
                    else:
                        # Store as orphan
                        orphan_lines[parent].append((line, curr_id))
                        
        return gene_lines

    def _extract_biotype_from_lines(self, lines: List[str]) -> str:
        """Extract gene biotype from grouped GFF lines."""
        gline = find_gene_feature_line(lines)
        if gline is None and lines:
            gline = lines[0]
        if not gline:
            return "unknown"
        for key in ["gene_biotype=", "biotype=", "gene_type="]:
            if key in gline:
                return gline.split(key)[1].split(";", 1)[0]
        return "unknown"

    def _candidate_rank(self, ref_norm: str, candidate: Dict[str, Any]) -> Tuple[Any, ...]:
        """Deterministic ranking for choosing one primary target locus per reference gene."""
        role = candidate.get("copy_role", "")
        if role in {"primary", "reference"}:
            role_rank = 0
        elif role in {"", "unknown"}:
            role_rank = 1
        else:
            role_rank = 2

        mapped_from_rank = 0 if candidate.get("has_mapped_from") else 1
        confidence = candidate.get("cnv_confidence", "")
        confidence_rank = 0 if confidence in {"", "high"} else 1
        same_id_rank = 0 if candidate.get("target_base") == ref_norm else 1
        copy_idx = candidate.get("copy_index")
        copy_rank = copy_idx if copy_idx is not None else 10**9

        return (
            role_rank,
            mapped_from_rank,
            same_id_rank,
            copy_rank,
            confidence_rank,
            candidate.get("target_id", ""),
        )

    def build_target_gene_links(
        self,
        ref_lines_map: Dict[str, List[str]],
        target_lines_map: Dict[str, List[str]],
    ) -> Tuple[Dict[str, str], Dict[str, List[Dict[str, Any]]], Dict[str, Dict[str, Any]]]:
        """Build robust reference->target linkages using IDs and mapped_from provenance."""
        ref_norm_to_id: Dict[str, str] = {}
        for ref_id in ref_lines_map.keys():
            ref_norm = normalize_gene_id(ref_id)
            if ref_norm and ref_norm not in ref_norm_to_id:
                ref_norm_to_id[ref_norm] = ref_id

        ref_to_candidates: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        target_meta: Dict[str, Dict[str, Any]] = {}
        unresolved_target_genes = 0

        for target_id, lines in target_lines_map.items():
            gene_line = find_gene_feature_line(lines)

            mapped_from_raw = extract_gff_attr(gene_line, "mapped_from") if gene_line else None
            mapped_from_norm = normalize_gene_id(mapped_from_raw)
            mapped_from_base = strip_copy_suffix(mapped_from_norm)

            target_norm = normalize_gene_id(target_id)
            target_base = strip_copy_suffix(target_norm)

            copy_index_raw = extract_gff_attr(gene_line, "cnv_copy_index") if gene_line else None
            copy_index = None
            if copy_index_raw is not None:
                try:
                    copy_index = int(copy_index_raw)
                except ValueError:
                    copy_index = None
            if copy_index is None:
                copy_index = parse_copy_index(target_norm)

            copy_role = ((extract_gff_attr(gene_line, "cnv_copy_role") if gene_line else None) or "").strip().lower()
            cnv_confidence = ((extract_gff_attr(gene_line, "cnv_assignment_confidence") if gene_line else None) or "").strip().lower()

            candidate = {
                "target_id": target_id,
                "target_norm": target_norm,
                "target_base": target_base,
                "mapped_from_norm": mapped_from_norm,
                "mapped_from_base": mapped_from_base,
                "has_mapped_from": bool(mapped_from_base),
                "copy_index": copy_index,
                "copy_role": copy_role,
                "cnv_confidence": cnv_confidence,
            }
            target_meta[target_id] = candidate

            candidate_refs: List[str] = []
            if mapped_from_base:
                candidate_refs.append(mapped_from_base)
            elif target_base:
                candidate_refs.append(target_base)

            added = False
            for ref_norm in dict.fromkeys(candidate_refs):
                ref_id = ref_norm_to_id.get(ref_norm)
                if ref_id is None:
                    continue
                ref_to_candidates[ref_id].append(candidate)
                added = True

            if not added:
                unresolved_target_genes += 1

        ref_to_primary: Dict[str, str] = {}
        for ref_id in ref_lines_map.keys():
            candidates = ref_to_candidates.get(ref_id, [])
            if not candidates:
                continue
            ref_norm = normalize_gene_id(ref_id)
            best = min(candidates, key=lambda c: self._candidate_rank(ref_norm, c))
            ref_to_primary[ref_id] = best["target_id"]

        logger.info(
            "Resolved reference->target links for %d/%d genes (%d target genes unresolved)",
            len(ref_to_primary),
            len(ref_lines_map),
            unresolved_target_genes,
        )

        return ref_to_primary, ref_to_candidates, target_meta

    def compute_cnv_summary(
        self,
        ref_lines_map: Dict[str, List[str]],
        ref_to_candidates: Dict[str, List[Dict[str, Any]]],
    ) -> Dict[str, Any]:
        """Compute copy-number expansion/contraction summary against reference=1 copy."""
        total_ref_genes = len(ref_lines_map)
        mapped_ref_genes = 0
        unchanged_genes = 0
        expansion_genes = 0
        contraction_genes = 0
        total_extra_copies = 0
        max_target_copy_number = 0
        ambiguous_assignment_genes = 0

        top_changes: List[Dict[str, Any]] = []
        role_counter: Counter = Counter()

        for ref_id in ref_lines_map.keys():
            ref_norm = normalize_gene_id(ref_id)
            unique_candidates = {
                c["target_id"]: c for c in ref_to_candidates.get(ref_id, [])
            }
            candidates = list(unique_candidates.values())
            target_copy_number = len(candidates)
            max_target_copy_number = max(max_target_copy_number, target_copy_number)

            if target_copy_number > 0:
                mapped_ref_genes += 1
            if target_copy_number == 0:
                contraction_genes += 1
            elif target_copy_number == 1:
                unchanged_genes += 1
            else:
                expansion_genes += 1
                total_extra_copies += (target_copy_number - 1)

            if candidates and any(
                (c.get("cnv_confidence", "") not in {"", "high"}) for c in candidates
            ):
                ambiguous_assignment_genes += 1

            for c in candidates:
                role = c.get("copy_role")
                if role:
                    role_counter[role] += 1

            delta = target_copy_number - 1
            if delta != 0:
                top_changes.append({
                    "gene_id": ref_norm,
                    "target_copies": target_copy_number,
                    "delta": delta,
                })

        top_changes.sort(
            key=lambda x: (abs(x["delta"]), x["target_copies"], x["gene_id"]),
            reverse=True,
        )

        return {
            "total_ref_genes": total_ref_genes,
            "mapped_ref_genes": mapped_ref_genes,
            "unchanged_genes": unchanged_genes,
            "expansion_genes": expansion_genes,
            "contraction_genes": contraction_genes,
            "total_extra_copies": total_extra_copies,
            "max_target_copy_number": max_target_copy_number,
            "ambiguous_assignment_genes": ambiguous_assignment_genes,
            "role_counts": dict(role_counter),
            "top_copy_number_changes": top_changes[:20],
        }


    def fast_load_genes(self, gff_path: str, count_per_type: int = 20) -> Set[str]:
        """Fast scan of GFF to find first N genes of each major biotype."""
        target_biotypes = {
            "protein_coding": [],
            "lncRNA": [],
            "pseudogene": [],
            "miRNA": [],
            "IG_V_gene": []
        }
        
        # Mapping simple names to regex/check logic if needed, but Biotype attribute is usually clear
        found_ids = set()
        
        # Scan file
        cmd = ["grep", "biotype=", str(gff_path)] # Rough filter to speed up simple line read?
        # Actually just read python line by line is okay if we stop early?
        # But grep is faster.
        # Let's simple read lines. output of grep | head limits is messy because they are interleaved.
        
        # We'll read the file line by line but stop when all buckets are full.
        
        full_buckets = 0
        with open(gff_path, 'r') as f:
            for line in f:
                if "\tgene\t" not in line: continue
                
                # Extract biotype (support biotype= or gene_type=)
                try:
                    parts = line.split("\t")
                    attrs = parts[8]
                    
                    biotype = None
                    if "biotype=" in attrs:
                        biotype = attrs.split("biotype=")[1].split(";")[0]
                    elif "gene_type=" in attrs:
                        biotype = attrs.split("gene_type=")[1].split(";")[0]
                        
                    if biotype:
                        # Extract ID
                        if "ID=" in attrs:
                            gid = attrs.split("ID=")[1].split(";")[0]
                            if gid.startswith("gene:"): gid = gid[5:]
                            
                            bucket = None
                            if "protein_coding" in biotype: bucket = "protein_coding"
                            elif "lnc" in biotype: bucket = "lncRNA"
                            elif "pseudogene" in biotype: bucket = "pseudogene"
                            elif "miRNA" in biotype: bucket = "miRNA"
                            elif "IG_V_gene" in biotype: bucket = "IG_V_gene"
                            
                            if bucket and len(target_biotypes[bucket]) < count_per_type:
                                target_biotypes[bucket].append(gid)
                                if len(target_biotypes[bucket]) == count_per_type:
                                    full_buckets += 1
                except:
                    continue
                    
                if full_buckets >= 5:
                    break
                    
        # Flatten
        for ids in target_biotypes.values():
            found_ids.update(ids)
            
        return found_ids

    def load_subset_gff(self, gff_path: str, gene_ids: Set[str]) -> Dict[str, Gene]:
        """Extract lines for specific genes and parse them."""
        # Create temp file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp_path = tmp.name
            
        try:
            # Create pattern file for grep -F
            # We need to grep the Gene ID, but separate Transcript/Exon IDs won't match.
            # However, usually Parent=GeneID works for transcripts.
            # And Parent=TranscriptID works for exons.
            # This recursive grep is hard with one pass.
            
            # Alternative: If parsing the whole file is too slow, 
            # maybe for test mode we accept missing children if grep is simple?
            # Or we iterate the whole file in Python once and write matching lines?
            # Iterating once in Python is usually < 10s for 1GB file. The Parser overhead is object creation.
            # So let's iterate and filter lines, then feed to parser.
            
            with open(gff_path, 'r') as infile, open(tmp_path, 'w') as outfile:
                # print headers
                outfile.write("##gff-version 3\n")
                
                # Capture logic: we need children.
                # Simplest robust way: 1 pass, store lines? No memory.
                # Actually, GFF3Parser logs show "Loaded X items".
                # Let's assume for Test Mode we just let GFF3Parser run on a filtered file
                # constructed by grep or similar.
                
                # Grep strategy:
                # 1. Grep genes.
                # 2. Grep Parent=GeneID (Transcripts).
                # 3. Grep Parent=TranscriptID (Exons).
                
                # This requires multiple passes or knowing transcript IDs.
                # Let's try to just use GFF3Parser with a filter?
                # GFF3Parser doesn't have an ID filter.
                
                # Fast Python Filter:
                # Read file. If line matches gene_id -> keep. If Parent matches kept_id -> keep.
                
                kept_ids = set(gene_ids)
                
                for line in infile:
                    if line.startswith("#"): continue
                    parts = line.split("\t")
                    if len(parts) < 9: continue
                    
                    attrs = parts[8]
                    
                    # Check ID
                    curr_id = None
                    if "ID=" in attrs:
                        curr_id = attrs.split("ID=")[1].split(";")[0]
                        if curr_id.startswith("gene:") or curr_id.startswith("transcript:"):
                            clean_id = curr_id.split(":")[1] if ":" in curr_id else curr_id
                        else:
                            clean_id = curr_id
                    
                    keep = False
                    
                    # If it's one of our target genes
                    if curr_id and (clean_id in gene_ids or (curr_id in gene_ids) or (curr_id.replace("gene:", "") in gene_ids)):
                        keep = True
                        kept_ids.add(curr_id) # Add full ID to kept for children
                    
                    # If Parent is in kept_ids
                    if "Parent=" in attrs:
                        parent = attrs.split("Parent=")[1].split(";")[0]
                        if parent in kept_ids:
                            keep = True
                            if curr_id: kept_ids.add(curr_id) # Add transcript ID to kept for exons
                            
                    if keep:
                        outfile.write(line)

            # Now parse the temp file
            parser = GFF3Parser()
            return parser.parse(tmp_path)
            
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def run(self):
        logger.info(f"Starting validation (Threads: {self.args.threads})")
        
        # 1. Load GFF Lines (Memory Efficient)
        ref_lines_map = self.load_gff_lines(self.args.ref_gff)
        logger.info(f"Loaded {len(ref_lines_map)} reference gene groups.")
        
        target_lines_map = self.load_gff_lines(self.args.target_gff)
        logger.info(f"Loaded {len(target_lines_map)} target gene groups.")
        ref_to_primary, ref_to_candidates, _target_gene_meta = self.build_target_gene_links(
            ref_lines_map, target_lines_map
        )
        
        # 2. Select Genes
        all_ids = list(ref_lines_map.keys())
        selected_ids = []
        
        if self.args.test_mode:
            logger.info("TEST MODE: Scanning for subset of genes (prioritizing mapped)...")
            
            counts = defaultdict(int)
            limit = self.args.test_mode_count
            
            # Prioritize finding genes that ARE mapped
            mapped_candidates = []
            unmapped_candidates = []
            
            for gid in all_ids:
                # Check if mapped
                is_mapped = gid in ref_to_primary
                
                if is_mapped:
                    mapped_candidates.append(gid)
                else:
                    unmapped_candidates.append(gid)
                    
            # Select mix (mostly mapped for verification)
            # Shuffle slightly to avoid bias if sorted by chrom
            import random
            random.seed(42)
            random.shuffle(mapped_candidates)
            
            for gid in mapped_candidates:
                 # Group bio
                lines = ref_lines_map[gid]
                biotype = self._extract_biotype_from_lines(lines)
                group = self.determine_biotype_group(biotype)
                
                if counts[group] < limit:
                    selected_ids.append(gid)
                    counts[group] += 1
                    
                if sum(counts.values()) >= (limit * 6): 
                    break
        else:
            # Subsample
            import random
            random.seed(self.args.seed)
            if self.args.percent >= 100:
                selected_ids = all_ids
            else:
                sample_size = int(len(all_ids) * (self.args.percent / 100))
                selected_ids = random.sample(all_ids, sample_size)
            
        logger.info(f"Identified {len(selected_ids)} genes for verification.")
        
        # 3. Create Task Batches
        batch_size = 50 
        task_batches = []
        current_batch = []
        
        for gid in selected_ids:
            target_id = ref_to_primary.get(gid)
            
            if target_id and target_id in target_lines_map:
                current_batch.append((gid, ref_lines_map[gid], target_id, target_lines_map[target_id]))
            
            if len(current_batch) >= batch_size:
                task_batches.append(current_batch)
                current_batch = []
        
        if current_batch:
            task_batches.append(current_batch)
            
        logger.info(f"Prepared {len(task_batches)} batches for processing.")
        
        # 4. Run Parallel
        # Use 'spawn' context for cross-platform compatibility (avoids fork issues with pysam)
        all_results = []
        all_alignments = []
        
        try:
            ctx = multiprocessing.get_context('spawn')
        except ValueError:
            # Fallback for older Python
            ctx = multiprocessing
        
        with ctx.Pool(self.args.threads) as pool:
            worker_args = []
            for batch_idx, batch in enumerate(task_batches):
                worker_args.append((
                    batch, 
                    self.args.ref_fasta, 
                    self.args.target_fasta, 
                    self.args.splice_window, 
                    self.args.splice_identity,
                    batch_idx
                ))
            
            for res_list, aln_list in pool.imap_unordered(worker_batch_process, worker_args):
                # Annotate group
                for r in res_list:
                    r["biotype_group"] = self.determine_biotype_group(r.get("biotype"))
                all_results.extend(res_list)
                all_alignments.extend(aln_list)
        
        # Unmapped
        mapped_set = {normalize_gene_id(r["gene_id"]) for r in all_results}
        
        for gid in selected_ids:
            norm_gid = normalize_gene_id(gid)
            if norm_gid not in mapped_set:
                # Extract biotype from ref_lines_map
                biotype = self._extract_biotype_from_lines(ref_lines_map.get(gid, []))
                        
                all_results.append({
                    "gene_id": norm_gid,
                    "biotype": biotype,
                    "biotype_group": self.determine_biotype_group(biotype),
                    "status": "unmapped"
                })
        
        # Mapping Metrics Output
        self.print_summary_table(all_results)
        
        # Synteny Analysis
        logger.info("Computing synteny conservation metrics...")
        
        # Build neighbor indices
        ref_index = GeneNeighborIndex(ref_lines_map, self.determine_biotype_group)
        target_index = GeneNeighborIndex(target_lines_map, self.determine_biotype_group)
        
        # Build COMPLETE ID mapping: ref_gene_id -> selected primary target_gene_id
        # (Not just selected genes - neighbors need full mapping)
        id_mapping = {
            gid: target_id
            for gid, target_id in ref_to_primary.items()
            if target_id in target_lines_map
        }
        
        synteny_results = compute_synteny_metrics(ref_index, target_index, selected_ids, id_mapping)
        self.print_synteny_table(synteny_results)

        # CNV summary (expansions / contractions / top changes)
        cnv_summary = self.compute_cnv_summary(ref_lines_map, ref_to_candidates)
        self.print_cnv_table(cnv_summary)
        
        # Write outputs
        self.write_outputs(all_results, all_alignments, synteny_results, cnv_summary)


    def print_summary_table(self, results):
        """Print formatted table to STDOUT."""
        # Aggregate by group
        groups = defaultdict(lambda: {
            "total": 0, "mapped": 0, 
            "ident_sum": 0, "cov_sum": 0,
            "exon_cons": 0, "splice_valid": 0, 
            "splice_candidate": 0, # Denominator for splice stats (ref_introns > 0)
            "cds_cons": 0,
            "kmer_aligned": 0  # Count of sequences using k-mer fallback
        })
        
        # Add "Overall" group explicitly or compute it?
        # Let's compute key 'Overall' alongside per-group
        
        for r in results:
            g = r.get("biotype_group", "Other")
            
            for key in [g, "Overall"]:
                groups[key]["total"] += 1
                if r["status"] == "mapped":
                    groups[key]["mapped"] += 1
                    groups[key]["ident_sum"] += r.get("transcript_identity", 0)
                    groups[key]["cov_sum"] += r.get("transcript_coverage", 0)
                    if r.get("exon_count_conserved"): groups[key]["exon_cons"] += 1
                    if r.get("cds_length_conserved"): groups[key]["cds_cons"] += 1
                    if r.get("alignment_method") == "kmer": groups[key]["kmer_aligned"] += 1
                
                # Splice Check (Intron Level):
                # Accumulate counts
                ref_introns = r.get("ref_intron_count", 0)
                if ref_introns > 0:
                    groups[key]["splice_candidate"] += ref_introns # Denominator is total REF introns
                    groups[key]["splice_valid"] += r.get("introns_conserved_count", 0) # Numerator is conserved

        print("\n" + "="*125)
        print(f"{'Biotype Group':<20} | {'Genes':<8} | {'Mapped%':<8} | {'Ident%':<8} | {'Cov%':<8} | {'ExonCons%':<10} | {'SpliceOK%':<25}")
        print("-" * 125)
        
        # Sort keys but put Overall last
        sorted_keys = sorted([k for k in groups.keys() if k != "Overall"])
        if "Overall" in groups:
            sorted_keys.append("Overall")
            
        for g in sorted_keys:
            stats = groups[g]
            mapped_pct = (stats["mapped"] / stats["total"] * 100) if stats["total"] else 0
            ident_avg = (stats["ident_sum"] / stats["mapped"]) if stats["mapped"] else 0
            cov_avg = (stats["cov_sum"] / stats["mapped"]) if stats["mapped"] else 0
            exon_cons_pct = (stats["exon_cons"] / stats["mapped"] * 100) if stats["mapped"] else 0
            
            # Splice OK % ignores single-exon genes (where splice_candidate is 0)
            if stats["splice_candidate"] > 0:
                splice_ok_pct = (stats["splice_valid"] / stats["splice_candidate"] * 100)
                splice_str = f"{splice_ok_pct:.1f}% ({stats['splice_valid']}/{stats['splice_candidate']})"
                # Pad splice_str to fit alignment
                splice_str = f"{splice_str:<25}" 
            else:
                splice_str = f"{'N/A':<25}"
            
            if g == "Overall":
                 print("-" * 125)
            
            print(f"{g:<20} | {stats['total']:<8} | {mapped_pct:<8.1f} | {ident_avg:<8.1f} | {cov_avg:<8.1f} | {exon_cons_pct:<10.1f} | {splice_str}")
        print("="*125)
        
        # Show k-mer fallback count if any
        overall_kmer = groups.get("Overall", {}).get("kmer_aligned", 0)
        if overall_kmer > 0:
            print(f"Note: {overall_kmer} genes used k-mer approximation (sequences > {MAX_ALIGNABLE_LENGTH/1000:.0f}kb)")
        print()

    def print_synteny_table(self, synteny_results):
        """Print formatted synteny conservation table to STDOUT."""
        # Aggregate by biotype group
        groups = defaultdict(lambda: {
            "total": 0,
            "pc_valid": 0, "pc_conserved": 0,
            "any_valid": 0, "any_conserved": 0,
            "neighborhood_count": 0, "neighborhood_sum": 0.0
        })
        
        for r in synteny_results:
            if r.get("synteny_status") != "mapped":
                continue  # Skip unmapped genes
                
            g = r.get("biotype_group", "Other")
            
            for key in [g, "Overall"]:
                groups[key]["total"] += 1
                
                # Method 1: Protein-coding neighbors
                if r.get("pc_neighbors_valid"):
                    groups[key]["pc_valid"] += 1
                    groups[key]["pc_conserved"] += r.get("pc_neighbors_conserved", 0)
                
                # Method 2: Any-biotype neighbors
                if r.get("any_neighbors_valid"):
                    groups[key]["any_valid"] += 1
                    groups[key]["any_conserved"] += r.get("any_neighbors_conserved", 0)
                
                # Method 3: Neighborhood score
                if r.get("neighborhood_score") is not None:
                    groups[key]["neighborhood_count"] += 1
                    groups[key]["neighborhood_sum"] += r.get("neighborhood_score", 0.0)
        
        print("\nSYNTENY CONSERVATION (Gene Order)")
        print("="*110)
        print(f"{'Biotype Group':<20} | {'Genes':<8} | {'PC Neighbors%':<14} | {'Any Neighbors%':<15} | {'Neighborhood Score':<20}")
        print("-" * 110)
        
        # Sort keys but put Overall last
        sorted_keys = sorted([k for k in groups.keys() if k != "Overall"])
        if "Overall" in groups:
            sorted_keys.append("Overall")
        
        for g in sorted_keys:
            stats = groups[g]
            
            # Protein-coding neighbors
            if stats["pc_valid"] > 0:
                pc_pct = (stats["pc_conserved"] / stats["pc_valid"]) * 100
                pc_str = f"{pc_pct:.1f}% ({stats['pc_conserved']}/{stats['pc_valid']})"
            else:
                pc_str = "N/A"
            
            # Any-biotype neighbors
            if stats["any_valid"] > 0:
                any_pct = (stats["any_conserved"] / stats["any_valid"]) * 100
                any_str = f"{any_pct:.1f}% ({stats['any_conserved']}/{stats['any_valid']})"
            else:
                any_str = "N/A"
            
            # Neighborhood score
            if stats["neighborhood_count"] > 0:
                neighborhood_avg = stats["neighborhood_sum"] / stats["neighborhood_count"]
                neighborhood_str = f"{neighborhood_avg:.3f}"
            else:
                neighborhood_str = "N/A"
            
            if g == "Overall":
                print("-" * 110)
            
            print(f"{g:<20} | {stats['total']:<8} | {pc_str:<14} | {any_str:<15} | {neighborhood_str:<20}")
        
        print("="*110)
        print()

    def print_cnv_table(self, cnv_summary: Dict[str, Any]):
        """Print CNV expansion/contraction summary to STDOUT."""
        total = cnv_summary.get("total_ref_genes", 0)
        mapped = cnv_summary.get("mapped_ref_genes", 0)
        unchanged = cnv_summary.get("unchanged_genes", 0)
        expansion = cnv_summary.get("expansion_genes", 0)
        contraction = cnv_summary.get("contraction_genes", 0)
        extra = cnv_summary.get("total_extra_copies", 0)
        max_copies = cnv_summary.get("max_target_copy_number", 0)
        ambiguous = cnv_summary.get("ambiguous_assignment_genes", 0)
        role_counts = cnv_summary.get("role_counts", {})

        print("\nCNV COPY-NUMBER SUMMARY")
        print("=" * 110)
        print(f"Reference genes evaluated: {total}")
        print(f"Mapped genes: {mapped}")
        print(f"Unchanged (1 copy): {unchanged}")
        print(f"Expansions (>1 copies): {expansion}")
        print(f"Contractions (0 copies): {contraction}")
        print(f"Total extra copies above reference: {extra}")
        print(f"Maximum observed target copy number: {max_copies}")
        print(f"Genes with non-high CNV confidence assignments: {ambiguous}")
        if role_counts:
            formatted_roles = ", ".join(f"{k}={v}" for k, v in sorted(role_counts.items()))
            print(f"CNV copy roles: {formatted_roles}")

        top_changes = cnv_summary.get("top_copy_number_changes", [])
        if top_changes:
            print("-" * 110)
            print("Top copy-number changes (target copies vs reference=1):")
            print(f"{'Gene ID':<24} | {'Target Copies':<13} | {'Delta':<6}")
            print("-" * 110)
            for row in top_changes[:10]:
                print(f"{row['gene_id']:<24} | {row['target_copies']:<13} | {row['delta']:<6}")
        print("=" * 110)
        print()

    def write_outputs(self, results, alignment_dumps, synteny_results=None, cnv_summary=None):
        metrics_file = "validation_metrics.json"
        align_file = self.args.dump_alignments
        
        final_output = {
            "summary": "See STDOUT or analyze details.",
            "details": results,
            "synteny": synteny_results or [],
            "cnv_summary": cnv_summary or {},
        }
        
        with open(metrics_file, "w") as f:
            json.dump(final_output, f, indent=2)
        logger.info(f"Wrote metrics to {metrics_file}")
        
        if align_file:
            with open(align_file, "w") as f:
                f.write("\n".join(alignment_dumps))
            logger.info(f"Wrote alignments to {align_file}")


def format_alignment(type_label, id_label, alignment, identity, coverage):
    header = f"\n> {type_label} ALIGNMENT: {id_label}\n"
    header += f"> Identity: {identity:.2f}% | Coverage: {coverage:.2f}%\n"
    return header + str(alignment)

# ... (worker_evaluate_pair is already updated)

def main():
    parser = argparse.ArgumentParser(description="Pangenome Mapping Validation")
    parser.add_argument("--ref-fasta", required=True, help="Reference FASTA")
    parser.add_argument("--ref-gff", required=True, help="Reference GFF3")
    parser.add_argument("--target-fasta", required=True, help="Target FASTA")
    parser.add_argument("--target-gff", required=True, help="Target GFF3")
    
    parser.add_argument("--percent", type=float, default=10.0, help="Subsample percentage")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads")
    parser.add_argument("--dump-alignments", help="Output alignment file")
    parser.add_argument("--chromosomes", nargs="+", help="Filter by chromosomes")
    
    parser.add_argument("--test-mode", action="store_true", help="Fast test mode: analyse first N genes per biotype")
    parser.add_argument("--test-mode-count", type=int, default=20, help="Number of genes to analyse in test mode per biotype")
    
    # New Args
    parser.add_argument("--splice-window", type=int, default=10, help="Flanking window size for splice check (default 10 -> 20bp window)")
    parser.add_argument("--splice-identity", type=float, default=0.9, help="Minimum identity threshold for splice window conservation (default 0.9)")
    
    args = parser.parse_args()
    
    runner = ValidationRunner(args)
    runner.run()

if __name__ == "__main__":
    main()

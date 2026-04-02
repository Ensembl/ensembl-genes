#!/usr/bin/env python3
"""
Investigate problematic mapping cases.

Analyzes transcripts with poor mapping quality, uses BLAST for independent
verification, and classifies cases as recoverable vs genuinely divergent.
"""

import argparse
import json
import logging
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
import pysam

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

# BLAST path (homebrew install location)
BLASTN = "/opt/homebrew/bin/blastn"
MAKEBLASTDB = "/opt/homebrew/bin/makeblastdb"


@dataclass
class ExonAnalysis:
    """Analysis result for a single exon."""
    exon_num: int
    ref_start: int
    ref_end: int
    ref_length: int
    
    # From mapping
    mapped_start: Optional[int] = None
    mapped_end: Optional[int] = None
    mapped_length: Optional[int] = None
    
    # Alignment quality (from validation alignment)
    alignment_identity: Optional[float] = None
    alignment_coverage: Optional[float] = None
    
    # BLAST verification
    blast_found: bool = False
    blast_identity: Optional[float] = None
    blast_coverage: Optional[float] = None
    blast_location: Optional[str] = None  # chrom:start-end
    
    # Classification
    status: str = "unknown"  # well_mapped, poorly_mapped, missing, blast_recoverable


@dataclass 
class TranscriptAnalysis:
    """Complete analysis for one transcript."""
    transcript_id: str
    gene_id: str
    biotype: str
    biotype_group: str
    
    # Validation metrics
    val_identity: float
    val_coverage: float
    
    # Exon analysis
    exons: List[ExonAnalysis] = field(default_factory=list)
    
    # Summary
    exons_well_mapped: int = 0
    exons_poorly_mapped: int = 0
    exons_missing: int = 0
    exons_blast_recoverable: int = 0
    
    # Final classification
    classification: str = "unknown"  # mapping_error, genuinely_divergent, partially_recoverable, fully_recoverable


def load_validation_metrics(json_path: str) -> List[Dict]:
    """Load validation metrics JSON."""
    with open(json_path) as f:
        data = json.load(f)
    return data.get("details", [])


def select_problem_transcripts(
    metrics: List[Dict],
    max_identity: float = 90.0,
    max_coverage: float = 90.0,
    samples_per_biotype: int = 5
) -> Dict[str, List[Dict]]:
    """Select transcripts with poor mapping for investigation."""
    
    # Filter to problem cases
    problems = [
        m for m in metrics
        if m.get("transcript_identity") and (
            m["transcript_identity"] < max_identity or 
            m.get("transcript_coverage", 100) < max_coverage
        )
    ]
    
    # Group by biotype
    by_biotype = defaultdict(list)
    for m in problems:
        by_biotype[m.get("biotype_group", "Other")].append(m)
    
    # Sample from each biotype
    selected = {}
    for biotype, items in by_biotype.items():
        # Sort by identity (worst first)
        items.sort(key=lambda x: x.get("transcript_identity", 100))
        selected[biotype] = items[:samples_per_biotype]
    
    return selected


def parse_gff_for_exons(gff_path: str, transcript_ids: Set[str]) -> Dict[str, List[Dict]]:
    """Extract exon coordinates for specified transcripts from GFF3."""
    exons_by_tx = defaultdict(list)
    
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            if feature_type != "exon":
                continue
            
            attrs = parts[8]
            parent_match = re.search(r'Parent=([^;]+)', attrs)
            if not parent_match:
                continue
            
            parent_id = parent_match.group(1)
            # Normalize ID
            norm_id = parent_id
            if norm_id.startswith("transcript:"):
                norm_id = norm_id[11:]
            if norm_id.startswith("mapped_transcript:"):
                norm_id = norm_id[18:]
            
            # Check if this matches any of our target transcripts
            for tx_id in transcript_ids:
                tx_norm = tx_id
                if tx_norm.startswith("transcript:"):
                    tx_norm = tx_norm[11:]
                if tx_norm.startswith("mapped_transcript:"):
                    tx_norm = tx_norm[18:]
                
                if norm_id == tx_norm:
                    exons_by_tx[tx_id].append({
                        "chrom": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6]
                    })
                    break
    
    # Sort exons by position
    for tx_id in exons_by_tx:
        exons_by_tx[tx_id].sort(key=lambda x: x["start"])
    
    return dict(exons_by_tx)


def extract_exon_sequences(
    fasta_path: str,
    exons: List[Dict],
    transcript_id: str
) -> List[Tuple[str, str]]:
    """Extract exon sequences from FASTA. Returns [(exon_id, sequence), ...]"""
    sequences = []
    
    handler = pysam.FastaFile(fasta_path)
    
    for i, exon in enumerate(exons):
        try:
            seq = handler.fetch(exon["chrom"], exon["start"] - 1, exon["end"])
            if exon.get("strand") == "-":
                seq = str(Seq(seq).reverse_complement())
            
            exon_id = f"{transcript_id}_exon{i+1}"
            sequences.append((exon_id, seq))
        except Exception as e:
            logger.warning(f"Failed to extract exon {i+1} for {transcript_id}: {e}")
    
    handler.close()
    return sequences


def run_blast_verification(
    exon_sequences: List[Tuple[str, str]],
    target_fasta: str,
    db_path: Optional[str] = None
) -> Dict[str, Dict]:
    """Run BLAST to verify exons can be found in target genome."""
    
    if not exon_sequences:
        return {}
    
    # Check/create BLAST database
    if db_path is None:
        db_path = target_fasta + ".blastdb"
    
    if not os.path.exists(db_path + ".nin") and not os.path.exists(db_path + ".nsq"):
        logger.info(f"Creating BLAST database for {target_fasta}...")
        cmd = [MAKEBLASTDB, "-in", target_fasta, "-dbtype", "nucl", "-out", db_path]
        subprocess.run(cmd, capture_output=True, check=True)
    
    # Write exon sequences to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        query_file = f.name
        for exon_id, seq in exon_sequences:
            if len(seq) >= 20:  # Minimum length for BLAST
                f.write(f">{exon_id}\n{seq}\n")
    
    results = {}
    
    try:
        # Run BLAST
        cmd = [
            BLASTN,
            "-query", query_file,
            "-db", db_path,
            "-outfmt", "6 qseqid sseqid pident length qlen sstart send evalue",
            "-evalue", "1e-5",
            "-max_target_seqs", "5",
            "-num_threads", "4"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            align_len = int(parts[3])
            qlen = int(parts[4])
            sstart = int(parts[5])
            send = int(parts[6])
            
            coverage = (align_len / qlen) * 100 if qlen > 0 else 0
            
            # Keep best hit per exon
            if qseqid not in results or pident > results[qseqid].get("identity", 0):
                results[qseqid] = {
                    "identity": pident,
                    "coverage": coverage,
                    "target_chrom": sseqid,
                    "target_start": min(sstart, send),
                    "target_end": max(sstart, send)
                }
    
    finally:
        os.unlink(query_file)
    
    return results


def analyze_transcript(
    metrics: Dict,
    ref_gff: str,
    target_gff: str,
    ref_fasta: str,
    target_fasta: str,
    blast_db: Optional[str] = None
) -> TranscriptAnalysis:
    """Perform full analysis of a transcript."""
    
    tx_id = metrics["transcript_id"]
    
    analysis = TranscriptAnalysis(
        transcript_id=tx_id,
        gene_id=metrics.get("gene_id", ""),
        biotype=metrics.get("biotype", "unknown"),
        biotype_group=metrics.get("biotype_group", "Other"),
        val_identity=metrics.get("transcript_identity", 0),
        val_coverage=metrics.get("transcript_coverage", 0)
    )
    
    # Get ref exons
    ref_exons = parse_gff_for_exons(ref_gff, {tx_id})
    ref_exon_list = ref_exons.get(tx_id, [])
    
    if not ref_exon_list:
        logger.warning(f"No ref exons found for {tx_id}")
        analysis.classification = "no_ref_data"
        return analysis
    
    # Get target exons (mapped)
    target_tx_id = f"mapped_transcript:{tx_id.replace('transcript:', '')}"
    target_exons = parse_gff_for_exons(target_gff, {target_tx_id, tx_id})
    target_exon_list = target_exons.get(target_tx_id, target_exons.get(tx_id, []))
    
    # Extract ref exon sequences for BLAST
    ref_exon_seqs = extract_exon_sequences(ref_fasta, ref_exon_list, tx_id)
    
    # Run BLAST
    blast_results = run_blast_verification(ref_exon_seqs, target_fasta, blast_db)
    
    # Analyze each exon
    for i, ref_exon in enumerate(ref_exon_list):
        exon_analysis = ExonAnalysis(
            exon_num=i + 1,
            ref_start=ref_exon["start"],
            ref_end=ref_exon["end"],
            ref_length=ref_exon["end"] - ref_exon["start"] + 1
        )
        
        # Check if mapped
        if i < len(target_exon_list):
            tgt = target_exon_list[i]
            exon_analysis.mapped_start = tgt["start"]
            exon_analysis.mapped_end = tgt["end"]
            exon_analysis.mapped_length = tgt["end"] - tgt["start"] + 1
        
        # Check BLAST
        exon_key = f"{tx_id}_exon{i+1}"
        if exon_key in blast_results:
            br = blast_results[exon_key]
            exon_analysis.blast_found = True
            exon_analysis.blast_identity = br["identity"]
            exon_analysis.blast_coverage = br["coverage"]
            exon_analysis.blast_location = f"{br['target_chrom']}:{br['target_start']}-{br['target_end']}"
        
        # Classify exon
        if exon_analysis.mapped_length:
            len_ratio = exon_analysis.mapped_length / exon_analysis.ref_length
            if 0.95 <= len_ratio <= 1.05 and exon_analysis.blast_identity and exon_analysis.blast_identity >= 95:
                exon_analysis.status = "well_mapped"
                analysis.exons_well_mapped += 1
            elif exon_analysis.blast_found and exon_analysis.blast_identity and exon_analysis.blast_identity >= 90:
                exon_analysis.status = "poorly_mapped"
                analysis.exons_poorly_mapped += 1
                if exon_analysis.blast_coverage and exon_analysis.blast_coverage >= 90:
                    exon_analysis.status = "blast_recoverable"
                    analysis.exons_blast_recoverable += 1
            else:
                exon_analysis.status = "poorly_mapped"
                analysis.exons_poorly_mapped += 1
        else:
            if exon_analysis.blast_found and exon_analysis.blast_identity and exon_analysis.blast_identity >= 90:
                exon_analysis.status = "blast_recoverable"
                analysis.exons_blast_recoverable += 1
            else:
                exon_analysis.status = "missing"
                analysis.exons_missing += 1
        
        analysis.exons.append(exon_analysis)
    
    # Overall classification
    total_exons = len(analysis.exons)
    if total_exons == 0:
        analysis.classification = "no_exons"
    elif analysis.exons_well_mapped == total_exons:
        analysis.classification = "well_mapped"
    elif analysis.exons_blast_recoverable > 0:
        if analysis.exons_well_mapped + analysis.exons_blast_recoverable == total_exons:
            analysis.classification = "fully_recoverable"
        else:
            analysis.classification = "partially_recoverable"
    elif analysis.exons_poorly_mapped > 0 or analysis.exons_missing > 0:
        analysis.classification = "genuinely_divergent"
    else:
        analysis.classification = "mapping_error"
    
    return analysis


def print_summary(analyses: List[TranscriptAnalysis]):
    """Print summary statistics."""
    
    print("\n" + "="*80)
    print("MAPPING ISSUE INVESTIGATION SUMMARY")
    print("="*80)
    
    # By classification
    by_class = defaultdict(list)
    for a in analyses:
        by_class[a.classification].append(a)
    
    print("\nClassification Summary:")
    print("-" * 50)
    for cls in ["fully_recoverable", "partially_recoverable", "mapping_error", "genuinely_divergent", "well_mapped", "no_exons"]:
        if cls in by_class:
            print(f"  {cls:25s}: {len(by_class[cls]):4d}")
    
    # By biotype
    print("\nBy Biotype Group:")
    print("-" * 50)
    by_biotype = defaultdict(lambda: defaultdict(int))
    for a in analyses:
        by_biotype[a.biotype_group][a.classification] += 1
    
    for biotype in sorted(by_biotype.keys()):
        counts = by_biotype[biotype]
        total = sum(counts.values())
        recoverable = counts.get("fully_recoverable", 0) + counts.get("partially_recoverable", 0)
        print(f"  {biotype:20s}: {total:4d} total, {recoverable:4d} recoverable ({100*recoverable/total:.1f}%)")
    
    # Exon-level stats
    total_exons = sum(len(a.exons) for a in analyses)
    well_mapped = sum(a.exons_well_mapped for a in analyses)
    blast_rec = sum(a.exons_blast_recoverable for a in analyses)
    poor = sum(a.exons_poorly_mapped for a in analyses)
    missing = sum(a.exons_missing for a in analyses)
    
    print(f"\nExon-Level Statistics:")
    print("-" * 50)
    print(f"  Total exons analyzed:    {total_exons}")
    print(f"  Well mapped:             {well_mapped} ({100*well_mapped/total_exons:.1f}%)")
    print(f"  BLAST recoverable:       {blast_rec} ({100*blast_rec/total_exons:.1f}%)")
    print(f"  Poorly mapped:           {poor} ({100*poor/total_exons:.1f}%)")
    print(f"  Missing:                 {missing} ({100*missing/total_exons:.1f}%)")
    
    print("="*80 + "\n")


def write_detailed_report(analyses: List[TranscriptAnalysis], output_path: str):
    """Write detailed TSV report."""
    
    with open(output_path, 'w') as f:
        # Header
        f.write("transcript_id\tgene_id\tbiotype\tbiotype_group\t")
        f.write("val_identity\tval_coverage\t")
        f.write("total_exons\twell_mapped\tpoorly_mapped\tmissing\tblast_recoverable\t")
        f.write("classification\n")
        
        for a in analyses:
            f.write(f"{a.transcript_id}\t{a.gene_id}\t{a.biotype}\t{a.biotype_group}\t")
            f.write(f"{a.val_identity:.1f}\t{a.val_coverage:.1f}\t")
            f.write(f"{len(a.exons)}\t{a.exons_well_mapped}\t{a.exons_poorly_mapped}\t")
            f.write(f"{a.exons_missing}\t{a.exons_blast_recoverable}\t")
            f.write(f"{a.classification}\n")
    
    logger.info(f"Wrote detailed report to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Investigate problematic mapping cases")
    parser.add_argument("--validation-json", required=True, help="Path to validation_metrics.json")
    parser.add_argument("--ref-gff", required=True, help="Reference GFF3")
    parser.add_argument("--target-gff", required=True, help="Target (mapped) GFF3")
    parser.add_argument("--ref-fasta", required=True, help="Reference FASTA")
    parser.add_argument("--target-fasta", required=True, help="Target FASTA")
    parser.add_argument("--max-identity", type=float, default=90.0, help="Max identity threshold for problems")
    parser.add_argument("--samples-per-biotype", type=int, default=10, help="Samples per biotype group")
    parser.add_argument("--output", default="investigation_report.tsv", help="Output TSV path")
    
    args = parser.parse_args()
    
    # Load metrics
    logger.info(f"Loading validation metrics from {args.validation_json}")
    metrics = load_validation_metrics(args.validation_json)
    logger.info(f"Loaded {len(metrics)} transcript records")
    
    # Select problem cases
    selected = select_problem_transcripts(
        metrics, 
        max_identity=args.max_identity,
        samples_per_biotype=args.samples_per_biotype
    )
    
    total_selected = sum(len(v) for v in selected.values())
    logger.info(f"Selected {total_selected} transcripts for investigation:")
    for biotype, items in selected.items():
        logger.info(f"  {biotype}: {len(items)}")
    
    # Create BLAST database (once)
    blast_db = args.target_fasta + ".blastdb"
    
    # Analyze each transcript
    all_analyses = []
    for biotype, items in selected.items():
        logger.info(f"Analyzing {biotype} transcripts...")
        for m in items:
            try:
                analysis = analyze_transcript(
                    m,
                    args.ref_gff,
                    args.target_gff,
                    args.ref_fasta,
                    args.target_fasta,
                    blast_db
                )
                all_analyses.append(analysis)
            except Exception as e:
                logger.warning(f"Failed to analyze {m.get('transcript_id')}: {e}")
    
    # Output
    print_summary(all_analyses)
    write_detailed_report(all_analyses, args.output)


if __name__ == "__main__":
    main()

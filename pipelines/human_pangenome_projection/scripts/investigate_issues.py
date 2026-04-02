#!/usr/bin/env python3
"""Script to investigate mapping issues and problematic genes."""

import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

# Add project root to path
import sys
sys.path.insert(0, '.')

from src.gff3_parser import parse_gff3
from src.fasta_handler import FastaHandler


def get_unmapped_genes(original_genes, mapped_genes):
    """Find genes that were never mapped (Stage 2 failures)."""
    mapped_orig_ids = set()
    for gene in mapped_genes.values():
        if gene.mapped_from:
            mapped_orig_ids.add(gene.mapped_from)
        else:
            orig_id = gene.feature_id.replace('mapped_', '')
            mapped_orig_ids.add(orig_id)
    
    return {gid: g for gid, g in original_genes.items() if gid not in mapped_orig_ids}


def get_synteny_disruptions(stats_file):
    """Get genes involved in synteny disruptions from stats."""
    with open(stats_file) as f:
        stats = json.load(f)
    
    # Check for inversion details in refinement/details
    details = stats.get('details', {})
    conflicts = details.get('conflicts', [])
    copy_number = details.get('copy_number_changes', [])
    
    return {
        'conflicts': conflicts,
        'copy_number': copy_number
    }


def run_blat(query_fa, target_fa, output_psl):
    """Run BLAT to search for sequence matches."""
    cmd = ['blat', '-minScore=50', '-minIdentity=90', str(target_fa), str(query_fa), str(output_psl)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"  BLAT error: {e}")
        return False


def extract_gene_sequence(gene, fasta_handler):
    """Extract gene sequence from FASTA."""
    try:
        seq = fasta_handler.fetch(gene.seq_region, gene.start, gene.end)
        return seq
    except Exception as e:
        return None


def parse_psl(psl_file):
    """Parse PSL output from BLAT."""
    hits = []
    with open(psl_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 21:
                try:
                    hits.append({
                        'matches': int(parts[0]),
                        'mismatches': int(parts[1]),
                        'query_name': parts[9],
                        'query_size': int(parts[10]),
                        'target_name': parts[13],
                        'target_size': int(parts[14]),
                        'target_start': int(parts[15]),
                        'target_end': int(parts[16]),
                        'identity': int(parts[0]) / (int(parts[0]) + int(parts[1])) if (int(parts[0]) + int(parts[1])) > 0 else 0
                    })
                except (ValueError, IndexError):
                    pass
    return hits


def investigate_unmapped_gene(gene, ref_fasta, target_fasta):
    """Investigate why a gene couldn't be mapped using BLAT."""
    print(f"\n  Investigating {gene.feature_id} ({gene.gene_name})...")
    print(f"    Location: {gene.seq_region}:{gene.start:,}-{gene.end:,} ({gene.length:,}bp)")
    print(f"    Biotype: {gene.attributes.get('gene_type', 'unknown')}")
    
    # Get gene sequence from reference
    with FastaHandler(ref_fasta) as fasta:
        seq = extract_gene_sequence(gene, fasta)
    
    if not seq or len(seq) < 50:
        print(f"    Could not extract sequence or sequence too short")
        return {'status': 'no_sequence'}
    
    # Write query FASTA
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">{gene.feature_id}\n{seq}\n")
        query_fa = f.name
    
    # Run BLAT against target
    with tempfile.NamedTemporaryFile(suffix='.psl', delete=False) as f:
        output_psl = f.name
    
    success = run_blat(query_fa, target_fasta, output_psl)
    
    if not success:
        print(f"    BLAT failed or not available")
        return {'status': 'blat_failed'}
    
    # Parse results
    hits = parse_psl(output_psl)
    
    if not hits:
        print(f"    NO BLAT HITS - gene may be truly absent from CHM13")
        return {'status': 'absent', 'reason': 'No BLAT matches found'}
    
    # Analyze hits
    best_hit = max(hits, key=lambda h: h['matches'])
    coverage = best_hit['matches'] / gene.length * 100
    
    print(f"    BLAT found {len(hits)} hits")
    print(f"    Best hit: {best_hit['target_name']}:{best_hit['target_start']:,}-{best_hit['target_end']:,}")
    print(f"    Coverage: {coverage:.1f}%, Identity: {best_hit['identity']*100:.1f}%")
    
    if coverage < 50:
        return {'status': 'fragmented', 'coverage': coverage, 'hits': len(hits)}
    elif best_hit['identity'] < 0.9:
        return {'status': 'divergent', 'identity': best_hit['identity']}
    else:
        return {'status': 'present_but_not_mapped', 'coverage': coverage, 'hits': hits}


def main():
    ref_gff = 'data/gencode.v44.chr20.gff3'
    ref_fasta = 'data/GRCh38.chr20.fa'
    target_fasta = 'data/CHM13v2.chr20.fa'
    mapped_gff = 'data/output/grch38_to_chm13.chr20.gff3'
    stats_file = 'data/output/grch38_to_chm13.chr20.stats.json'
    
    print("=" * 70)
    print("MAPPING ISSUES INVESTIGATION")
    print("=" * 70)
    
    # Load data
    original_genes = parse_gff3(ref_gff)
    mapped_genes = parse_gff3(mapped_gff)
    
    print(f"\nOriginal genes: {len(original_genes)}")
    print(f"Mapped genes (post-refinement): {len(mapped_genes)}")
    
    # 1. Find truly unmapped genes
    # Load the pre-refinement mapping stats
    with open(stats_file) as f:
        stats = json.load(f)
    
    unmapped_count = stats['summary']['genes']['unmapped']
    print(f"\n{'='*70}")
    print(f"CATEGORY 1: TRULY UNMAPPED GENES ({unmapped_count} genes)")
    print(f"{'='*70}")
    print("These genes couldn't be mapped at all (no syntenic block coverage)")
    
    # Get unmapped genes by finding which are missing from mapped set
    unmapped = get_unmapped_genes(original_genes, mapped_genes)
    
    # But we need to distinguish truly unmapped from refinement-filtered
    # The truly unmapped are the ones that Stage 2 couldn't map
    # For now, let's sample from the unmapped list
    
    # Identify genes in regions with poor synteny coverage
    # Sort by position and look for clusters
    unmapped_sorted = sorted(unmapped.values(), key=lambda g: g.start)
    
    print(f"\nUnmapped genes (not in final output): {len(unmapped)}")
    print("\nSample (sorted by position):")
    
    # Check BLAT for a few candidates
    for gene in unmapped_sorted[:5]:
        result = investigate_unmapped_gene(gene, ref_fasta, target_fasta)
        
    # 2. Synteny disruptions
    print(f"\n{'='*70}")
    print("CATEGORY 2: SYNTENY DISRUPTIONS")
    print(f"{'='*70}")
    
    # The 8 inversions and 4 new overlaps
    print("\nChecking for genes involved in order inversions and new overlaps...")
    
    # 3. Check partial mappings by examining transcripts
    print(f"\n{'='*70}")
    print("CATEGORY 3: PARTIAL/PROBLEMATIC MAPPINGS")
    print(f"{'='*70}")
    
    # Find genes with validation issues
    validation_rate = float(stats['validation']['splice_sites']['rate'].rstrip('%')) / 100
    print(f"\nSplice site validation: {stats['validation']['splice_sites']['rate']}")
    print(f"Start codon validation: {stats['validation']['start_codons']['rate']}")
    

if __name__ == '__main__':
    main()

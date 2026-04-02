"""Secondary mapping for genes not covered by synteny blocks.

This module provides gap-filling functionality for genes that fall in
regions without synteny coverage (e.g., segmental duplications).

Memory-efficient approach: batches all genes into a single minimap2 call
to avoid loading multiple copies of the target genome.
"""

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from src.models import Gene, GenomicInterval, Strand, SyntenicBlock, SyntenicMap

logger = logging.getLogger(__name__)


def find_uncovered_genes(
    genes: Dict[str, Gene],
    syntenic_map: SyntenicMap
) -> List[Gene]:
    """Find genes not covered by any synteny block.
    
    Args:
        genes: Reference genes
        syntenic_map: Syntenic map with coverage blocks
        
    Returns:
        List of genes without synteny coverage
    """
    uncovered = []
    
    for gene_id, gene in genes.items():
        blocks = syntenic_map.find_blocks_for_ref_region(
            gene.seq_region, gene.start, gene.end
        )
        
        if not blocks:
            uncovered.append(gene)
        else:
            # Compute union coverage to avoid double-counting overlapping blocks
            overlap_intervals = []
            for block in blocks:
                overlap_start = max(gene.start, block.ref_interval.start)
                overlap_end = min(gene.end, block.ref_interval.end)
                if overlap_end >= overlap_start:
                    overlap_intervals.append((overlap_start, overlap_end))

            overlap_intervals.sort()
            merged = []
            for start, end in overlap_intervals:
                if not merged or start > merged[-1][1] + 1:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)

            covered_bp = sum(end - start + 1 for start, end in merged)
            coverage = covered_bp / gene.length
            if coverage < 0.5:  # Less than 50% covered
                uncovered.append(gene)
    
    logger.info(f"Found {len(uncovered)} genes without sufficient synteny coverage")
    return uncovered


def _batch_align_genes(
    genes: List[Gene],
    ref_fasta: Path,
    target_fasta: Path,
    threads: int = 4,
    min_identity: float = 0.7,  # Relaxed from 0.9 for difficult regions
    min_coverage: float = 0.5,  # Relaxed from 0.8 for difficult regions
    timeout_per_gene: int = 10,
    temp_dir: Optional[Path] = None
) -> Dict[str, dict]:
    """Align all genes in a single minimap2 call.
    
    Uses minimap2's internal threading (-t) to parallelize alignment
    while only loading the target genome once.
    
    Args:
        genes: List of genes to align
        ref_fasta: Reference FASTA file
        target_fasta: Target FASTA file
        threads: Number of threads for minimap2
        min_identity: Minimum alignment identity
        min_coverage: Minimum query coverage
        timeout_per_gene: Timeout per gene (multiplied by gene count)
        
    Returns:
        Dictionary mapping gene_id to hit data
    """
    from src.fasta_handler import FastaHandler
    
    if not genes:
        return {}
    
    # Extract all gene sequences
    gene_data = {}
    with FastaHandler(ref_fasta) as fasta:
        for gene in genes:
            try:
                seq = fasta.fetch(gene.seq_region, gene.start, gene.end)
                # Keep short genes too - they might be important (e.g. minimal exons)
                if len(seq) >= 30:  # Relaxed from 50
                    gene_data[gene.feature_id] = {
                        'seq': seq,
                        'gene': gene
                    }
            except Exception as e:
                logger.debug(f"Could not extract sequence for {gene.feature_id}: {e}")
    
    if not gene_data:
        return {}
    
    logger.info(f"Extracted sequences for {len(gene_data)} genes")
    
    # Write all gene sequences to a single FASTA
    query_fa = None
    try:
        # Use temp_dir if provided, otherwise system temp
        temp_path = temp_dir or Path(tempfile.gettempdir())
        query_fa = temp_path / f"gap_fill_queries_{os.getpid()}.fa"
        with open(query_fa, 'w') as f:
            for gene_id, data in gene_data.items():
                f.write(f">{gene_id}\n{data['seq']}\n")
        
        # Run single minimap2 with internal threading
        # This loads the target genome only once
        cmd = [
            'minimap2',
            '-c',                    # Output PAF with CIGAR
            '-x', 'asm20',          # Use asm20 (up to 20% divergence) for difficult regions
            '--cs=long',            # Long cs tag for coordinate projection
            '-N', '20',             # Report up to 20 secondary alignments (increased from 5 for repetitive families)
            '-t', str(threads),     # Use multiple threads internally
            str(target_fasta),
            query_fa
        ]
        
        # Calculate timeout based on gene count - but cap it reasonably
        timeout = min(max(300, len(genes) * timeout_per_gene), 7200)  # 5 min to 2 hours
        
        logger.info(f"Running minimap2 with {threads} threads for {len(gene_data)} genes...")
        
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                timeout=timeout
            )
        except subprocess.TimeoutExpired:
            logger.warning(f"Batch alignment timed out after {timeout}s")
            return {}
        
        if result.returncode != 0 and result.stderr:
            logger.warning(f"minimap2 stderr: {result.stderr[:500]}")
        
        # Parse results - find best hit for each gene
        gene_hits = {}
        
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 12:
                continue
            
            try:
                gene_id = parts[0]
                if gene_id not in gene_data:
                    continue
                
                query_len = int(parts[1])
                q_start = int(parts[2])
                q_end = int(parts[3])
                strand = parts[4]
                target_chr = parts[5]
                t_start = int(parts[7])
                t_end = int(parts[8])
                matches = int(parts[9])
                align_len = int(parts[10])
                
                coverage = (q_end - q_start) / query_len
                identity = matches / align_len if align_len > 0 else 0
                
                if coverage >= min_coverage and identity >= min_identity:
                    score = matches * identity
                    
                    cs_tag = None
                    for tag in parts[12:]:
                        if tag.startswith('cs:'):
                            cs_tag = tag[5:]  # Remove "cs:Z:"
                            break
                    
                    gene = gene_data[gene_id]['gene']
                    hit = {
                        'gene_id': gene_id,
                        'target_chr': target_chr,
                        't_start': t_start + 1,  # Convert to 1-based
                        't_end': t_end,
                        'strand': strand,
                        'identity': identity,
                        'coverage': coverage,
                        'matches': matches,
                        'align_len': align_len,
                        'cs_tag': cs_tag,
                        'seq_region': gene.seq_region,
                        'start': gene.start,
                        'end': gene.end,
                        'gene_strand': '+' if gene.strand == Strand.PLUS else '-',
                        'score': score
                    }
                    
                    # Collect ALL hits per gene for later locus assignment
                    if gene_id not in gene_hits:
                        gene_hits[gene_id] = []
                    gene_hits[gene_id].append(hit)
                    
            except (ValueError, IndexError) as e:
                continue
        
        logger.info(f"Found alignments for {len(gene_hits)}/{len(gene_data)} genes")
        return gene_hits
        
    finally:
        if query_fa and os.path.exists(query_fa):
            try:
                os.unlink(query_fa)
            except:
                pass


def fill_synteny_gaps(
    genes: Dict[str, Gene],
    syntenic_map: SyntenicMap,
    ref_fasta: Path,
    target_fasta: Path,
    threads: int = 4,
    temp_dir: Optional[Path] = None
) -> Tuple[SyntenicMap, int]:
    """Fill synteny gaps by creating local blocks for uncovered genes.
    
    Uses a single minimap2 call with all genes batched together,
    with internal threading for parallelization. This is memory-efficient
    as the target genome is only loaded once.
    
    Args:
        genes: Reference genes
        syntenic_map: Existing syntenic map
        ref_fasta: Reference FASTA
        target_fasta: Target FASTA
        threads: Number of threads for minimap2
        temp_dir: Directory for temp files (defaults to output dir)
        
    Returns:
        Tuple of (updated SyntenicMap, number of blocks added)
    """
    uncovered = find_uncovered_genes(genes, syntenic_map)
    
    if not uncovered:
        logger.info("No genes need gap filling")
        return syntenic_map, 0
    
    logger.info(f"Attempting to fill gaps for {len(uncovered)} genes using {threads} threads...")
    
    # Sort genes by ID for deterministic batching
    uncovered.sort(key=lambda g: g.feature_id)
    
    # Process in chunks to avoid timeouts with large gene sets
    chunk_size = 100
    all_gene_hits = {}
    
    total_chunks = (len(uncovered) + chunk_size - 1) // chunk_size
    
    for i in range(0, len(uncovered), chunk_size):
        chunk = uncovered[i:i + chunk_size]
        logger.info(f"Processing gap filling batch {i//chunk_size + 1}/{total_chunks} ({len(chunk)} genes)...")
        
        batch_hits = _batch_align_genes(
            chunk,
            ref_fasta,
            target_fasta,
            threads=threads,
            temp_dir=temp_dir
        )
        all_gene_hits.update(batch_hits)
    
    # Assign genes to unique loci to prevent paralog collapse
    assigned_hits = _assign_genes_to_loci(all_gene_hits, uncovered)
    
    # Create synteny blocks from assigned hits
    blocks_added = 0
    for gene_id, hit in assigned_hits.items():
        block = SyntenicBlock(
            ref_interval=GenomicInterval(
                seq_region=hit['seq_region'],
                start=hit['start'],
                end=hit['end'],
                strand=Strand.PLUS if hit['gene_strand'] == '+' else Strand.MINUS
            ),
            target_interval=GenomicInterval(
                seq_region=hit['target_chr'],
                start=hit['t_start'],
                end=hit['t_end'],
                strand=Strand.PLUS if hit['strand'] == '+' else Strand.MINUS
            ),
            identity=hit['identity'],
            alignment_length=hit['align_len'],
            matches=hit['matches'],
            mismatches=hit['align_len'] - hit['matches'],
            cs_tag=hit['cs_tag'],
            block_id=f"local_{gene_id}"
        )
        syntenic_map.blocks.append(block)
        blocks_added += 1
    
    # Rebuild index if blocks were added
    if blocks_added > 0:
        syntenic_map.build_index()
        logger.info(f"Added {blocks_added} local synteny blocks for gap regions")
    
    failed = len(uncovered) - blocks_added
    if failed > 0:
        logger.info(f"{failed} genes could not be mapped (low identity/coverage or unmappable)")
    
    return syntenic_map, blocks_added


def _assign_genes_to_loci(
    all_gene_hits: Dict[str, List[dict]],
    genes: List[Gene]
) -> Dict[str, dict]:
    """Assign genes to unique target loci using greedy matching.
    
    For paralogous genes that could map to multiple loci, this ensures
    each gene gets assigned to a different target locus where possible,
    sorted by reference position to maintain synteny.
    
    Args:
        all_gene_hits: Dict mapping gene_id to list of candidate hits
        genes: List of genes (for reference position sorting)
        
    Returns:
        Dict mapping gene_id to the assigned hit
    """
    # Build gene order by reference position
    gene_order = {g.feature_id: (g.seq_region, g.start) for g in genes}
    
    # Sort genes by reference position
    genes_with_hits = [
        g for g in genes 
        if g.feature_id in all_gene_hits
    ]
    genes_sorted = sorted(
        genes_with_hits, 
        key=lambda g: gene_order.get(g.feature_id, ('', 0))
    )
    
    # Track which target loci (intervals) are already used
    # List of (chr, start, end) tuples
    used_intervals: List[Tuple[str, int, int]] = []
    
    assigned = {}
    
    for gene in genes_sorted:
        gene_id = gene.feature_id
        hits = all_gene_hits.get(gene_id, [])
        
        if not hits:
            continue
        
        # Sort hits by score (best first)
        hits_sorted = sorted(hits, key=lambda h: h['score'], reverse=True)
        
        # Find best available hit (not overlapping with used intervals)
        for hit in hits_sorted:
            hit_chr = hit['target_chr']
            hit_start = hit['t_start']
            hit_end = hit['t_end']
            
            # Check if this locus overlaps with any used interval
            locus_used = False
            for used_chr, used_start, used_end in used_intervals:
                if used_chr == hit_chr:
                    # Check for actual coordinate overlap
                    if hit_start <= used_end and hit_end >= used_start:
                        locus_used = True
                        break
            
            if not locus_used:
                # Assign gene to this locus
                assigned[gene_id] = hit
                used_intervals.append((hit_chr, hit_start, hit_end))
                break
        else:
            # No available locus - use best hit anyway (will conflict)
            assigned[gene_id] = hits_sorted[0]
    
    return assigned

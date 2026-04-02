"""Stage 3: Refinement and conflict resolution.

Handles:
- Detecting unmapped regions for novel features
- Resolving overlapping/conflicting mappings
- Copy number analysis
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from src.models import Gene, GenomicInterval, Strand, SyntenicMap, Transcript
from src.config import MappingStatus
from validation.structural import GeneValidation

logger = logging.getLogger(__name__)


@dataclass
class ConflictRegion:
    """A region where multiple features mapped but shouldn't overlap."""
    interval: GenomicInterval
    conflicting_genes: List[Gene]
    selected_gene: Optional[Gene] = None
    selection_reason: str = ""


@dataclass
class CopyNumberChange:
    """A gene with copy number difference between reference and target."""
    gene_id: str
    gene_name: Optional[str]
    ref_copy_count: int
    target_copy_count: int
    change_type: str  # "expansion", "contraction", "normal"
    target_gene_ids: List[str] = field(default_factory=list)
    inference: str = "direct"
    supporting_genes: List[str] = field(default_factory=list)
    labeling_strategy: str = "synteny_first_parsimony"
    cluster_id: str = ""
    ambiguity_tag: str = ""
    contraction_label: str = ""


@dataclass
class UnmappedRegion:
    """A region in target that has no mapped features from reference."""
    interval: GenomicInterval
    nearest_upstream_gene: Optional[str] = None
    nearest_downstream_gene: Optional[str] = None
    potential_novel: bool = False


class ConflictResolver:
    """Resolves conflicts between overlapping mapped features."""
    
    def __init__(
        self,
        mapped_genes: Dict[str, Gene],
        original_genes: Optional[Dict[str, Gene]] = None,
        validation_results: Optional[Dict[str, GeneValidation]] = None,
        overlap_threshold: float = 0.5
    ):
        """Initialize resolver.
        
        Args:
            mapped_genes: Dictionary of mapped genes
            original_genes: Original reference genes (for checking reference overlaps)
            overlap_threshold: Minimum reciprocal overlap to consider conflict
        """
        self.mapped_genes = mapped_genes
        self.original_genes = original_genes or {}
        self.validation_results = validation_results or {}
        self.overlap_threshold = overlap_threshold
        self.conflicts: List[ConflictRegion] = []
    
    def find_overlapping_genes(self) -> Dict[str, List[Gene]]:
        """Find genes that overlap on the target genome.
        
        Returns:
            Dictionary of chromosome -> list of overlapping gene groups
        """
        # Group genes by chromosome
        genes_by_chr: Dict[str, List[Gene]] = defaultdict(list)
        for gene in self.mapped_genes.values():
            genes_by_chr[gene.seq_region].append(gene)
        
        overlapping = {}
        
        for chrom, genes in genes_by_chr.items():
            # Sort by start position
            sorted_genes = sorted(genes, key=lambda g: g.start)
            
            # Find overlapping groups
            groups = []
            current_group = []
            current_end = 0
            
            for gene in sorted_genes:
                if gene.start > current_end and current_group:
                    if len(current_group) > 1:
                        groups.append(current_group)
                    current_group = [gene]
                    current_end = gene.end
                else:
                    current_group.append(gene)
                    current_end = max(current_end, gene.end)
            
            if len(current_group) > 1:
                groups.append(current_group)
            
            if groups:
                overlapping[chrom] = groups
        
        return overlapping
    
    def _calculate_overlap(self, gene1: Gene, gene2: Gene) -> float:
        """Calculate reciprocal overlap between two genes."""
        if gene1.seq_region != gene2.seq_region:
            return 0.0
        
        if gene1.strand != gene2.strand:
            return 0.0
        
        overlap_start = max(gene1.start, gene2.start)
        overlap_end = min(gene1.end, gene2.end)
        
        if overlap_start > overlap_end:
            return 0.0
        
        overlap_len = overlap_end - overlap_start + 1
        min_len = min(gene1.length, gene2.length)
        
        return overlap_len / min_len if min_len > 0 else 0.0
    
    def _did_genes_overlap_in_reference(self, gene1: Gene, gene2: Gene) -> bool:
        """Check if two genes overlapped in the reference genome.
        
        Uses the mapped_from attribute to trace back to original coordinates.
        Returns True if genes overlapped in reference (nested genes are valid).
        """
        orig1_id = gene1.mapped_from
        orig2_id = gene2.mapped_from
        
        if not orig1_id or not orig2_id:
            return False
        
        # Look up original genes
        orig1 = self.original_genes.get(orig1_id)
        orig2 = self.original_genes.get(orig2_id)
        
        if not orig1 or not orig2:
            return False
        
        # Check if on same chromosome
        if orig1.seq_region != orig2.seq_region:
            return False
        
        # Check if intervals overlap
        return orig1.start <= orig2.end and orig2.start <= orig1.end
    
    def _gene_quality_score(
        self,
        gene: Gene,
        ref_idx: int,
        target_idx: int
    ) -> Tuple[float, int, float, int, int]:
        """Composite score for one conflict candidate (higher is better)."""
        gene_val = self.validation_results.get(gene.feature_id)
        if gene_val:
            tx_checked = max(1, gene_val.transcripts_checked)
            structure_score = gene_val.transcripts_valid / tx_checked
            if gene_val.is_valid:
                structure_score += 1.0
        else:
            structure_score = 0.0

        cnv_role = gene.attributes.get("cnv_copy_role")
        if cnv_role == "primary":
            cnv_score = 2
        elif cnv_role == "expansion":
            cnv_score = 0
        else:
            cnv_score = 1

        identity = gene.mapping_identity or 0.0
        order_delta = abs(ref_idx - target_idx)
        confidence = gene.attributes.get("projection_confidence_class", "")
        if confidence == "high":
            confidence_score = 2
        elif confidence in {"medium", "biological_divergence_suspected"}:
            confidence_score = 1
        else:
            confidence_score = 0

        return (
            structure_score,
            cnv_score,
            identity,
            confidence_score,
            -order_delta,
        )
    
    def _get_reference_order(self, genes: List[Gene]) -> Dict[str, int]:
        """Get the reference genome order for a list of genes.
        
        Returns dict mapping gene_id to reference position index.
        """
        ref_positions = {}
        for gene in genes:
            orig_id = gene.mapped_from or gene.feature_id
            orig = self.original_genes.get(orig_id)
            if orig:
                ref_positions[gene.feature_id] = (orig.seq_region, orig.start)
            else:
                ref_positions[gene.feature_id] = (gene.seq_region, gene.start)
        
        # Sort by reference position and assign indices
        sorted_genes = sorted(genes, key=lambda g: ref_positions[g.feature_id])
        return {g.feature_id: i for i, g in enumerate(sorted_genes)}
    
    def _resolve_overlap_group_synteny_aware(
        self,
        group: List[Gene],
        chrom: str
    ) -> Tuple[List[Gene], List[Tuple[Gene, str]]]:
        """Resolve an overlapping group preserving synteny.
        
        Only removes genes that DIRECTLY overlap with kept genes.
        Non-overlapping genes are always kept regardless of order.
        
        Args:
            group: List of genes in overlapping region
            chrom: Chromosome name
            
        Returns:
            Tuple of (genes to keep, list of (removed gene, removal reason))
        """
        if len(group) <= 1:
            return group, []
        
        # Get reference order for these genes
        ref_order = self._get_reference_order(group)
        
        # Sort group by target position
        sorted_by_target = sorted(group, key=lambda g: g.start)
        
        target_rank = {g.feature_id: i for i, g in enumerate(sorted_by_target)}

        # Process each gene; when overlaps happen, pick one-best using quality+synteny+CNV.
        kept_genes = []
        removed_genes = []
        
        for gene in sorted_by_target:
            ref_idx = ref_order[gene.feature_id]
            
            # Check if this gene overlaps with any already-kept gene
            overlapping_kept = None
            for kept in kept_genes:
                overlap = self._calculate_overlap(gene, kept)
                if overlap >= self.overlap_threshold:
                    # Check if they overlapped in reference (valid nested genes)
                    if not self._did_genes_overlap_in_reference(gene, kept):
                        overlapping_kept = kept
                        break
            
            if overlapping_kept:
                candidate_score = self._gene_quality_score(
                    gene,
                    ref_idx=ref_idx,
                    target_idx=target_rank[gene.feature_id],
                )
                kept_score = self._gene_quality_score(
                    overlapping_kept,
                    ref_idx=ref_order[overlapping_kept.feature_id],
                    target_idx=target_rank[overlapping_kept.feature_id],
                )

                if candidate_score > kept_score:
                    kept_genes = [
                        g for g in kept_genes
                        if g.feature_id != overlapping_kept.feature_id
                    ]
                    removed_genes.append((
                        overlapping_kept,
                        (
                            "overlap_conflict: replaced by "
                            f"{gene.feature_id} "
                            f"(score={candidate_score} > {kept_score})"
                        ),
                    ))
                    kept_genes.append(gene)
                else:
                    removed_genes.append((
                        gene,
                        (
                            "overlap_conflict: suppressed by "
                            f"{overlapping_kept.feature_id} "
                            f"(score={candidate_score} <= {kept_score})"
                        ),
                    ))
            else:
                # No conflict - ALWAYS keep this gene
                kept_genes.append(gene)
        
        return kept_genes, removed_genes
    
    def resolve_conflicts(self) -> Tuple[Dict[str, Gene], List[ConflictRegion]]:
        """Resolve conflicts using synteny-aware assignment.
        
        Preserves reference gene order where possible. Handles CNV regions
        by keeping the assignment that best maintains synteny.
        
        Returns:
            Tuple of (filtered_genes, conflict_regions)
        """
        overlapping = self.find_overlapping_genes()
        
        genes_to_remove: Set[str] = set()
        removal_reasons: Dict[str, str] = {}
        conflicts = []
        
        for chrom, groups in overlapping.items():
            for group in groups:
                kept, removed = self._resolve_overlap_group_synteny_aware(group, chrom)
                
                for gene, reason in removed:
                    genes_to_remove.add(gene.feature_id)
                    removal_reasons[gene.feature_id] = reason
                
                if removed:
                    conflict = ConflictRegion(
                        interval=GenomicInterval(
                            seq_region=chrom,
                            start=min(g.start for g in group),
                            end=max(g.end for g in group),
                            strand=group[0].strand
                        ),
                        conflicting_genes=group,
                        selected_gene=kept[0] if kept else None,
                        selection_reason=f"kept {len(kept)}, removed {len(removed)} (synteny-aware)"
                    )
                    conflicts.append(conflict)
        
        # Create filtered gene set
        filtered_genes = {
            gid: gene for gid, gene in self.mapped_genes.items()
            if gid not in genes_to_remove
        }
        
        logger.info(
            f"Resolved {len(conflicts)} conflict regions, removed {len(genes_to_remove)} genes"
        )
        
        # Store removal reasons for reporting
        self.removal_reasons = removal_reasons
        self.conflicts = conflicts
        return filtered_genes, conflicts


class CopyNumberAnalyzer:
    """Analyzes copy number changes between reference and target."""
    
    def __init__(
        self,
        original_genes: Dict[str, Gene],
        mapped_genes: Dict[str, Gene]
    ):
        """Initialize analyzer.
        
        Args:
            original_genes: Original reference genes
            mapped_genes: Mapped target genes
        """
        self.original_genes = original_genes
        self.mapped_genes = mapped_genes
    
    def analyze(self) -> List[CopyNumberChange]:
        """Analyze copy number changes.
        
        Returns:
            List of CopyNumberChange records
        """
        # Track how many times each reference gene was mapped
        ref_to_target: Dict[str, List[str]] = defaultdict(list)
        
        for target_id, target_gene in self.mapped_genes.items():
            ref_id = target_gene.mapped_from
            if ref_id:
                ref_to_target[ref_id].append(target_id)
        
        # Reference neighborhood index for contraction inference.
        ref_by_chr: Dict[str, List[str]] = defaultdict(list)
        for ref_id, ref_gene in self.original_genes.items():
            ref_by_chr[ref_gene.seq_region].append(ref_id)
        for chrom in ref_by_chr:
            ref_by_chr[chrom].sort(key=lambda rid: self.original_genes[rid].start)
        
        expanded_ref_ids = {
            ref_id for ref_id, tgt in ref_to_target.items()
            if len(tgt) > 1
        }
        
        changes = []
        
        for ref_id, ref_gene in self.original_genes.items():
            target_copies = ref_to_target.get(ref_id, [])
            target_count = len(target_copies)
            
            if target_count == 0:
                change_type = "contraction"
            elif target_count == 1:
                change_type = "normal"
            else:
                change_type = "expansion"
            
            if change_type != "normal":
                inference = "direct"
                supporting = []
                cnv_cluster_id = ref_id
                ambiguity_tag = ""
                contraction_label = ""
                
                if change_type == "expansion":
                    copy_confidences = []
                    for tgt_id in target_copies:
                        tgt_gene = self.mapped_genes.get(tgt_id)
                        if tgt_gene:
                            copy_confidences.append(
                                tgt_gene.attributes.get("cnv_assignment_confidence", "high")
                            )
                            if tgt_gene.attributes.get("cnv_ambiguity_tag"):
                                ambiguity_tag = tgt_gene.attributes["cnv_ambiguity_tag"]
                    if any(c == "ambiguous" for c in copy_confidences):
                        inference = "ambiguous_expansion_assignment"
                    else:
                        inference = "synteny_supported_expansion"
                    cluster_sources = set([ref_id])
                    for tgt_id in target_copies:
                        tgt_gene = self.mapped_genes.get(tgt_id)
                        if not tgt_gene:
                            continue
                        cluster_id = tgt_gene.attributes.get("cnv_cluster_id", "")
                        if cluster_id:
                            for src in cluster_id.split(","):
                                if src:
                                    cluster_sources.add(src)
                    cnv_cluster_id = ",".join(sorted(cluster_sources))
                
                if change_type == "contraction":
                    chrom_ids = ref_by_chr.get(ref_gene.seq_region, [])
                    contraction_label = "single_gene_loss"
                    if ref_id in chrom_ids:
                        idx = chrom_ids.index(ref_id)
                        neighbors = []
                        for j in range(max(0, idx - 2), min(len(chrom_ids), idx + 3)):
                            if j == idx:
                                continue
                            nbr = chrom_ids[j]
                            if nbr in expanded_ref_ids:
                                neighbors.append(nbr)
                        if neighbors:
                            inference = "ambiguous_cluster_contraction"
                            supporting = neighbors
                            contraction_label = "cluster_ambiguous_loss"
                            ambiguity_tag = "cluster_context"
                            cnv_cluster_id = ",".join(sorted(set([ref_id] + neighbors)))
                        else:
                            inference = "high_confidence_contraction"
                
                changes.append(CopyNumberChange(
                    gene_id=ref_id,
                    gene_name=ref_gene.gene_name,
                    ref_copy_count=1,
                    target_copy_count=target_count,
                    change_type=change_type,
                    target_gene_ids=target_copies,
                    inference=inference,
                    supporting_genes=supporting,
                    labeling_strategy="deterministic_cluster_order",
                    cluster_id=cnv_cluster_id,
                    ambiguity_tag=ambiguity_tag,
                    contraction_label=contraction_label
                ))
        
        logger.info(
            f"Copy number analysis: {sum(1 for c in changes if c.change_type == 'expansion')} expansions, "
            f"{sum(1 for c in changes if c.change_type == 'contraction')} contractions"
        )
        
        return changes


class UnmappedRegionFinder:
    """Finds regions in target without mapped features."""
    
    def __init__(
        self,
        mapped_genes: Dict[str, Gene],
        target_chromosome_lengths: Dict[str, int],
        min_gap_size: int = 50000
    ):
        """Initialize finder.
        
        Args:
            mapped_genes: Mapped genes on target
            target_chromosome_lengths: Length of each target chromosome
            min_gap_size: Minimum gap size to report
        """
        self.mapped_genes = mapped_genes
        self.target_chromosome_lengths = target_chromosome_lengths
        self.min_gap_size = min_gap_size
    
    def find_gaps(self) -> List[UnmappedRegion]:
        """Find large gaps between mapped genes.
        
        Returns:
            List of UnmappedRegion
        """
        gaps = []
        
        # Group genes by chromosome
        genes_by_chr: Dict[str, List[Gene]] = defaultdict(list)
        for gene in self.mapped_genes.values():
            genes_by_chr[gene.seq_region].append(gene)
        
        for chrom, chrom_length in self.target_chromosome_lengths.items():
            genes = sorted(genes_by_chr.get(chrom, []), key=lambda g: g.start)
            
            # Check gap at start
            if genes and genes[0].start > self.min_gap_size:
                gaps.append(UnmappedRegion(
                    interval=GenomicInterval(
                        seq_region=chrom,
                        start=1,
                        end=genes[0].start - 1,
                        strand=Strand.UNSTRANDED
                    ),
                    nearest_downstream_gene=genes[0].feature_id
                ))
            
            # Check gaps between genes
            for i in range(len(genes) - 1):
                gap_start = genes[i].end + 1
                gap_end = genes[i+1].start - 1
                gap_size = gap_end - gap_start + 1
                
                if gap_size >= self.min_gap_size:
                    gaps.append(UnmappedRegion(
                        interval=GenomicInterval(
                            seq_region=chrom,
                            start=gap_start,
                            end=gap_end,
                            strand=Strand.UNSTRANDED
                        ),
                        nearest_upstream_gene=genes[i].feature_id,
                        nearest_downstream_gene=genes[i+1].feature_id
                    ))
            
            # Check gap at end
            if genes and chrom_length - genes[-1].end > self.min_gap_size:
                gaps.append(UnmappedRegion(
                    interval=GenomicInterval(
                        seq_region=chrom,
                        start=genes[-1].end + 1,
                        end=chrom_length,
                        strand=Strand.UNSTRANDED
                    ),
                    nearest_upstream_gene=genes[-1].feature_id
                ))
        
        logger.info(f"Found {len(gaps)} unmapped regions >= {self.min_gap_size}bp")
        
        return gaps


def refine_mappings(
    original_genes: Dict[str, Gene],
    mapped_genes: Dict[str, Gene],
    target_chromosome_lengths: Dict[str, int],
    validation_results: Optional[Dict[str, GeneValidation]] = None,
    resolve_conflicts: bool = True,
    analyze_copy_number: bool = True,
    find_unmapped_regions: bool = True
) -> Tuple[Dict[str, Gene], Dict]:
    """Run Stage 3 refinement on mapped genes.
    
    Args:
        original_genes: Original reference genes
        mapped_genes: Mapped target genes
        target_chromosome_lengths: Target chromosome lengths
        validation_results: Optional structural validation for conflict scoring
        resolve_conflicts: Whether to resolve conflicts
        analyze_copy_number: Whether to analyze copy number
        find_unmapped_regions: Whether to find unmapped regions
        
    Returns:
        Tuple of (refined_genes, refinement_report)
    """
    report = {
        "conflicts": [],
        "copy_number_changes": [],
        "unmapped_regions": [],
        "conflict_filtered_genes": []  # Genes that mapped but were filtered due to overlap
    }
    
    refined_genes = dict(mapped_genes)
    conflict_filtered_ids: Set[str] = set()  # Track genes removed by conflict resolution
    
    # Resolve conflicts
    if resolve_conflicts:
        resolver = ConflictResolver(
            refined_genes,
            original_genes,
            validation_results=validation_results,
        )
        refined_genes, conflicts = resolver.resolve_conflicts()
        
        # Track which original genes were filtered due to conflicts
        # Use detailed removal reasons from synteny-aware resolver
        removal_reasons = getattr(resolver, 'removal_reasons', {})
        
        for gene_id, reason in removal_reasons.items():
            gene = mapped_genes.get(gene_id)
            if gene:
                conflict_filtered_ids.add(gene.mapped_from or gene_id)
                report["conflict_filtered_genes"].append({
                    "gene_id": gene.mapped_from or gene_id,
                    "gene_name": gene.gene_name,
                    "reason": reason,
                    "mapped_position": f"{gene.seq_region}:{gene.start}-{gene.end}"
                })
        
        report["conflicts"] = [
            {
                "region": f"{c.interval.seq_region}:{c.interval.start}-{c.interval.end}",
                "conflicting_count": len(c.conflicting_genes),
                "selected": c.selected_gene.feature_id if c.selected_gene else None,
                "reason": c.selection_reason
            }
            for c in conflicts
        ]
    
    # Analyze copy number - use PRE-CONFLICT mapped genes to avoid counting
    # conflict-filtered genes as contractions
    if analyze_copy_number:
        # Use original mapped_genes, not refined_genes
        analyzer = CopyNumberAnalyzer(original_genes, mapped_genes)
        cn_changes = analyzer.analyze()
        
        # Filter out genes that were filtered due to conflicts - these are NOT true contractions
        true_contractions = [
            c for c in cn_changes 
            if c.change_type == "contraction" and c.gene_id not in conflict_filtered_ids
        ]
        other_changes = [c for c in cn_changes if c.change_type != "contraction"]
        
        report["copy_number_changes"] = [
            {
                "gene_id": c.gene_id,
                "gene_name": c.gene_name,
                "ref_copies": c.ref_copy_count,
                "target_copies": c.target_copy_count,
                "change_type": c.change_type,
                "target_gene_ids": c.target_gene_ids,
                "inference": c.inference,
                "supporting_genes": c.supporting_genes,
                "labeling_strategy": c.labeling_strategy,
                "cluster_id": c.cluster_id,
                "ambiguity_tag": c.ambiguity_tag,
                "contraction_label": c.contraction_label
            }
            for c in true_contractions + other_changes
        ]
        
        # Log accurate numbers
        logger.info(
            f"Copy number: {len(true_contractions)} true contractions, "
            f"{len(conflict_filtered_ids)} filtered by conflicts"
        )
    
    # Find unmapped regions
    if find_unmapped_regions and target_chromosome_lengths:
        finder = UnmappedRegionFinder(refined_genes, target_chromosome_lengths)
        unmapped = finder.find_gaps()
        report["unmapped_regions"] = [
            {
                "region": f"{u.interval.seq_region}:{u.interval.start}-{u.interval.end}",
                "size": u.interval.length,
                "upstream": u.nearest_upstream_gene,
                "downstream": u.nearest_downstream_gene
            }
            for u in unmapped
        ]
    
    return refined_genes, report

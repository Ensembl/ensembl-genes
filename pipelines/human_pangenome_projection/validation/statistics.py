"""Statistics and reporting for mapping results."""

import json
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from src.models import Gene, Transcript
from src.config import MappingStatus
from stages.mapping import GeneMappingResult
from validation.structural import GeneValidation

logger = logging.getLogger(__name__)


# Ensembl major biotype groups
# See: https://www.ensembl.org/info/genome/genebuild/biotypes.html
BIOTYPE_GROUPS = {
    'protein_coding': 'coding',
    'IG_C_gene': 'coding',
    'IG_D_gene': 'coding',
    'IG_J_gene': 'coding',
    'IG_V_gene': 'coding',
    'TR_C_gene': 'coding',
    'TR_D_gene': 'coding',
    'TR_J_gene': 'coding',
    'TR_V_gene': 'coding',
    
    'lncRNA': 'long_noncoding',
    'lincRNA': 'long_noncoding',
    'antisense': 'long_noncoding',
    'sense_intronic': 'long_noncoding',
    'sense_overlapping': 'long_noncoding',
    'processed_transcript': 'long_noncoding',
    'macro_lncRNA': 'long_noncoding',
    'bidirectional_promoter_lncRNA': 'long_noncoding',
    
    'miRNA': 'small_noncoding',
    'rRNA': 'small_noncoding',
    'snRNA': 'small_noncoding',
    'snoRNA': 'small_noncoding',
    'tRNA': 'small_noncoding',
    'scRNA': 'small_noncoding',
    'scaRNA': 'small_noncoding',
    'vaultRNA': 'small_noncoding',
    'sRNA': 'small_noncoding',
    'misc_RNA': 'small_noncoding',
    'ribozyme': 'small_noncoding',
    'Mt_tRNA': 'small_noncoding',
    'Mt_rRNA': 'small_noncoding',
    
    'pseudogene': 'pseudogene',
    'processed_pseudogene': 'pseudogene',
    'unprocessed_pseudogene': 'pseudogene',
    'transcribed_processed_pseudogene': 'pseudogene',
    'transcribed_unprocessed_pseudogene': 'pseudogene',
    'transcribed_unitary_pseudogene': 'pseudogene',
    'unitary_pseudogene': 'pseudogene',
    'polymorphic_pseudogene': 'pseudogene',
    'IG_C_pseudogene': 'pseudogene',
    'IG_J_pseudogene': 'pseudogene',
    'IG_V_pseudogene': 'pseudogene',
    'IG_pseudogene': 'pseudogene',
    'TR_J_pseudogene': 'pseudogene',
    'TR_V_pseudogene': 'pseudogene',
    'rRNA_pseudogene': 'pseudogene',
    'translated_processed_pseudogene': 'pseudogene',
    'translated_unprocessed_pseudogene': 'pseudogene',
}

# Human-readable names for major biotype groups
BIOTYPE_GROUP_NAMES = {
    'coding': 'Protein Coding',
    'long_noncoding': 'Long Non-coding RNA',
    'small_noncoding': 'Small Non-coding RNA', 
    'pseudogene': 'Pseudogene',
    'other': 'Other'
}


def get_biotype_group(biotype: Optional[str]) -> str:
    """Get the major biotype group for a given biotype."""
    if biotype is None:
        return 'other'
    return BIOTYPE_GROUPS.get(biotype, 'other')


@dataclass
class BiotypeStat:
    """Statistics for a single biotype or biotype group."""
    name: str
    total_genes: int = 0
    mapped_genes: int = 0
    total_transcripts: int = 0
    mapped_transcripts: int = 0
    
    @property
    def gene_mapping_rate(self) -> float:
        if self.total_genes == 0:
            return 1.0
        return self.mapped_genes / self.total_genes
    
    @property
    def transcript_mapping_rate(self) -> float:
        if self.total_transcripts == 0:
            return 1.0
        return self.mapped_transcripts / self.total_transcripts
    
    def to_dict(self) -> Dict:
        return {
            'name': self.name,
            'genes': {
                'total': self.total_genes,
                'mapped': self.mapped_genes,
                'rate': self.gene_mapping_rate
            },
            'transcripts': {
                'total': self.total_transcripts,
                'mapped': self.mapped_transcripts,
                'rate': self.transcript_mapping_rate
            }
        }


@dataclass
class GeneOrderAnalysis:
    """Analysis of gene order conservation between assemblies."""
    
    # Total gene pairs analyzed
    total_pairs: int = 0
    
    # Pairs where relative order is preserved
    order_preserved: int = 0
    
    # Pairs where relative order is inverted
    order_inverted: int = 0
    
    # Pairs where genes are on different chromosomes in target
    different_chromosomes: int = 0
    
    # Pairs that couldn't be analyzed (unmapped genes)
    unmappable_pairs: int = 0
    
    # Detailed inversions: list of (gene1_id, gene2_id, ref_order, target_order)
    inversion_details: List[Dict] = field(default_factory=list)
    
    @property
    def order_conservation_rate(self) -> float:
        analyzable = self.total_pairs - self.unmappable_pairs
        if analyzable == 0:
            return 1.0
        return self.order_preserved / analyzable


@dataclass
class OverlapConflict:
    """A conflict where non-overlapping reference genes overlap in target."""
    gene1_id: str
    gene1_name: Optional[str]
    gene2_id: str
    gene2_name: Optional[str]
    
    # Reference positions (non-overlapping)
    ref_gene1_start: int
    ref_gene1_end: int
    ref_gene2_start: int
    ref_gene2_end: int
    
    # Target positions (overlapping)
    target_gene1_start: int
    target_gene1_end: int
    target_gene2_start: int
    target_gene2_end: int
    
    overlap_bp: int  # Number of overlapping bases
    
    def to_dict(self) -> Dict:
        return {
            "gene1": {"id": self.gene1_id, "name": self.gene1_name},
            "gene2": {"id": self.gene2_id, "name": self.gene2_name},
            "reference": {
                "gene1": f"{self.ref_gene1_start}-{self.ref_gene1_end}",
                "gene2": f"{self.ref_gene2_start}-{self.ref_gene2_end}",
                "overlap": False
            },
            "target": {
                "gene1": f"{self.target_gene1_start}-{self.target_gene1_end}",
                "gene2": f"{self.target_gene2_start}-{self.target_gene2_end}",
                "overlap_bp": self.overlap_bp
            }
        }


@dataclass
class SyntenyAnalysis:
    """Combined synteny analysis results."""
    gene_order: GeneOrderAnalysis
    new_overlaps: List[OverlapConflict]
    
    @property
    def new_overlap_count(self) -> int:
        return len(self.new_overlaps)


@dataclass
class MappingStatistics:
    """Comprehensive mapping statistics."""
    
    # Input counts
    total_genes: int = 0
    total_transcripts: int = 0
    total_exons: int = 0
    total_cds: int = 0
    
    # Mapping success counts (pre-refinement)
    genes_mapped: int = 0
    genes_partial: int = 0
    genes_unmapped: int = 0
    
    transcripts_mapped: int = 0
    transcripts_partial: int = 0
    transcripts_unmapped: int = 0
    
    exons_mapped: int = 0
    exons_unmapped: int = 0
    
    cds_mapped: int = 0
    cds_unmapped: int = 0
    
    # Post-refinement output counts
    genes_in_output: int = 0
    transcripts_in_output: int = 0
    genes_removed_by_refinement: int = 0
    
    # Validation
    splice_sites_checked: int = 0
    splice_sites_valid: int = 0
    start_codons_checked: int = 0
    start_codons_valid: int = 0
    stop_codons_checked: int = 0
    stop_codons_valid: int = 0
    
    # Identity
    mean_mapping_identity: float = 0.0
    min_mapping_identity: float = 0.0
    max_mapping_identity: float = 0.0
    
    # Performance
    total_time_seconds: float = 0.0
    synteny_time_seconds: float = 0.0
    mapping_time_seconds: float = 0.0
    validation_time_seconds: float = 0.0
    
    # Biotype statistics (populated in calculate_statistics)
    # Major biotype group stats (coding, long_noncoding, small_noncoding, pseudogene, other)
    biotype_group_stats: Dict[str, BiotypeStat] = field(default_factory=dict)
    
    # Individual biotype stats (protein_coding, lncRNA, miRNA, etc.)
    biotype_stats: Dict[str, BiotypeStat] = field(default_factory=dict)

    
    @property
    def gene_mapping_rate(self) -> float:
        return self.genes_mapped / max(self.total_genes, 1)
    
    @property
    def transcript_mapping_rate(self) -> float:
        return self.transcripts_mapped / max(self.total_transcripts, 1)
    
    @property
    def exon_mapping_rate(self) -> float:
        return self.exons_mapped / max(self.total_exons, 1)
    
    @property
    def splice_site_validation_rate(self) -> float:
        return self.splice_sites_valid / max(self.splice_sites_checked, 1)
    
    @property
    def start_codon_validation_rate(self) -> float:
        return self.start_codons_valid / max(self.start_codons_checked, 1)
    
    @property
    def stop_codon_validation_rate(self) -> float:
        return self.stop_codons_valid / max(self.stop_codons_checked, 1)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary with computed rates."""
        d = asdict(self)
        d["gene_mapping_rate"] = self.gene_mapping_rate
        d["transcript_mapping_rate"] = self.transcript_mapping_rate
        d["exon_mapping_rate"] = self.exon_mapping_rate
        d["splice_site_validation_rate"] = self.splice_site_validation_rate
        d["start_codon_validation_rate"] = self.start_codon_validation_rate
        d["stop_codon_validation_rate"] = self.stop_codon_validation_rate
        
        # Properly serialize biotype stats
        d["biotype_group_stats"] = {
            k: v.to_dict() for k, v in self.biotype_group_stats.items()
        }
        d["biotype_stats"] = {
            k: v.to_dict() for k, v in self.biotype_stats.items()
        }
        return d


def calculate_statistics(
    original_genes: Dict[str, Gene],
    mapping_results: List[GeneMappingResult],
    validation_results: Optional[Dict[str, GeneValidation]] = None
) -> MappingStatistics:
    """Calculate comprehensive mapping statistics.
    
    Args:
        original_genes: Original reference genes
        mapping_results: Results from feature mapping
        validation_results: Optional validation results
        
    Returns:
        MappingStatistics object
    """
    stats = MappingStatistics()
    
    # Initialize biotype tracking
    # Track counts per biotype group and individual biotype
    biotype_group_counts = defaultdict(lambda: {'total_genes': 0, 'mapped_genes': 0, 
                                                  'total_transcripts': 0, 'mapped_transcripts': 0})
    biotype_counts = defaultdict(lambda: {'total_genes': 0, 'mapped_genes': 0,
                                           'total_transcripts': 0, 'mapped_transcripts': 0})
    
    # Count original features and collect biotypes
    stats.total_genes = len(original_genes)
    for gene in original_genes.values():
        stats.total_transcripts += len(gene.transcripts)
        for tx in gene.transcripts:
            stats.total_exons += len(tx.exons)
            stats.total_cds += len(tx.cds_list)
        
        # Get biotype from gene or from attributes
        biotype = gene.biotype
        if biotype is None and gene.attributes:
            biotype = gene.attributes.get('gene_type') or gene.attributes.get('biotype')
        
        group = get_biotype_group(biotype)
        biotype_key = biotype or 'unknown'
        
        # Track totals
        biotype_group_counts[group]['total_genes'] += 1
        biotype_group_counts[group]['total_transcripts'] += len(gene.transcripts)
        biotype_counts[biotype_key]['total_genes'] += 1
        biotype_counts[biotype_key]['total_transcripts'] += len(gene.transcripts)
    
    # Count mapping results
    identities = []
    
    for result in mapping_results:
        # Get biotype for this gene
        gene = result.original
        biotype = gene.biotype
        if biotype is None and gene.attributes:
            biotype = gene.attributes.get('gene_type') or gene.attributes.get('biotype')
        group = get_biotype_group(biotype)
        biotype_key = biotype or 'unknown'
        
        # Track gene mapping by status
        if result.status == MappingStatus.MAPPED:
            stats.genes_mapped += 1
            biotype_group_counts[group]['mapped_genes'] += 1
            biotype_counts[biotype_key]['mapped_genes'] += 1
        elif result.status == MappingStatus.PARTIAL:
            stats.genes_partial += 1
            # Partial also counts as mapped for biotype stats
            biotype_group_counts[group]['mapped_genes'] += 1
            biotype_counts[biotype_key]['mapped_genes'] += 1
        else:
            stats.genes_unmapped += 1
        
        for tx_result in result.transcript_results:
            if tx_result.status == MappingStatus.MAPPED:
                stats.transcripts_mapped += 1
                biotype_group_counts[group]['mapped_transcripts'] += 1
                biotype_counts[biotype_key]['mapped_transcripts'] += 1
            elif tx_result.status == MappingStatus.PARTIAL:
                stats.transcripts_partial += 1
                biotype_group_counts[group]['mapped_transcripts'] += 1
                biotype_counts[biotype_key]['mapped_transcripts'] += 1
            else:
                stats.transcripts_unmapped += 1
            
            # Count exons and CDS
            for exon_result in tx_result.exon_results:
                if exon_result.status == MappingStatus.MAPPED:
                    stats.exons_mapped += 1
                    if exon_result.mapped and exon_result.mapped.mapping_identity:
                        identities.append(exon_result.mapped.mapping_identity)
                else:
                    stats.exons_unmapped += 1
            
            for cds_result in tx_result.cds_results:
                if cds_result.status == MappingStatus.MAPPED:
                    stats.cds_mapped += 1
                else:
                    stats.cds_unmapped += 1
    
    # Identity statistics
    if identities:
        stats.mean_mapping_identity = sum(identities) / len(identities)
        stats.min_mapping_identity = min(identities)
        stats.max_mapping_identity = max(identities)
    
    # Validation statistics
    if validation_results:
        for gene_val in validation_results.values():
            for tx_val in gene_val.transcript_results:
                stats.splice_sites_checked += tx_val.splice_sites_checked
                stats.splice_sites_valid += tx_val.splice_sites_valid
                
                if tx_val.start_codon_result:
                    stats.start_codons_checked += 1
                    if tx_val.start_codon_valid:
                        stats.start_codons_valid += 1
                
                if tx_val.stop_codon_result:
                    stats.stop_codons_checked += 1
                    if tx_val.stop_codon_valid:
                        stats.stop_codons_valid += 1
    
    # Convert biotype counts to BiotypeStat objects
    for group, counts in biotype_group_counts.items():
        stats.biotype_group_stats[group] = BiotypeStat(
            name=BIOTYPE_GROUP_NAMES.get(group, group),
            total_genes=counts['total_genes'],
            mapped_genes=counts['mapped_genes'],
            total_transcripts=counts['total_transcripts'],
            mapped_transcripts=counts['mapped_transcripts']
        )
    
    for biotype, counts in biotype_counts.items():
        stats.biotype_stats[biotype] = BiotypeStat(
            name=biotype,
            total_genes=counts['total_genes'],
            mapped_genes=counts['mapped_genes'],
            total_transcripts=counts['total_transcripts'],
            mapped_transcripts=counts['mapped_transcripts']
        )
    
    return stats



def generate_report(
    stats: MappingStatistics,
    refinement_report: Optional[Dict] = None,
    output_path: Optional[Union[str, Path]] = None
) -> Dict:
    """Generate a comprehensive report.
    
    Args:
        stats: Mapping statistics
        refinement_report: Optional refinement stage report
        output_path: Optional path to save JSON report
        
    Returns:
        Report dictionary
    """
    report = {
        "timestamp": datetime.now().isoformat(),
        "summary": {
            "genes": {
                "total": stats.total_genes,
                "mapped": stats.genes_mapped,
                "partial": stats.genes_partial,
                "unmapped": stats.genes_unmapped,
                "mapping_rate": f"{stats.gene_mapping_rate:.2%}"
            },
            "transcripts": {
                "total": stats.total_transcripts,
                "mapped": stats.transcripts_mapped,
                "partial": stats.transcripts_partial,
                "unmapped": stats.transcripts_unmapped,
                "mapping_rate": f"{stats.transcript_mapping_rate:.2%}"
            },
            "exons": {
                "total": stats.total_exons,
                "mapped": stats.exons_mapped,
                "unmapped": stats.exons_unmapped,
                "mapping_rate": f"{stats.exon_mapping_rate:.2%}"
            }
        },
        "validation": {
            "splice_sites": {
                "checked": stats.splice_sites_checked,
                "valid": stats.splice_sites_valid,
                "rate": f"{stats.splice_site_validation_rate:.2%}"
            },
            "start_codons": {
                "checked": stats.start_codons_checked,
                "valid": stats.start_codons_valid,
                "rate": f"{stats.start_codon_validation_rate:.2%}"
            },
            "stop_codons": {
                "checked": stats.stop_codons_checked,
                "valid": stats.stop_codons_valid,
                "rate": f"{stats.stop_codon_validation_rate:.2%}"
            }
        },
        "identity": {
            "mean": stats.mean_mapping_identity,
            "min": stats.min_mapping_identity,
            "max": stats.max_mapping_identity
        },
        "performance": {
            "total_time_seconds": stats.total_time_seconds,
            "synteny_time_seconds": stats.synteny_time_seconds,
            "mapping_time_seconds": stats.mapping_time_seconds,
            "validation_time_seconds": stats.validation_time_seconds
        },
        "biotype_group_stats": {
            k: v.to_dict() for k, v in stats.biotype_group_stats.items()
        },
        "biotype_stats": {
            k: v.to_dict() for k, v in stats.biotype_stats.items()
        }
    }
    
    if refinement_report:
        report["refinement"] = {
            "conflicts_resolved": len(refinement_report.get("conflicts", [])),
            "copy_number_expansions": sum(
                1 for c in refinement_report.get("copy_number_changes", [])
                if c.get("change_type") == "expansion"
            ),
            "copy_number_contractions": sum(
                1 for c in refinement_report.get("copy_number_changes", [])
                if c.get("change_type") == "contraction"
            ),
            "unmapped_regions": len(refinement_report.get("unmapped_regions", []))
        }
        
        # Include details
        report["details"] = refinement_report
    
    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {output_path}")
    
    return report


def print_summary(stats: MappingStatistics):
    """Print a summary to console."""
    print("\n" + "="*60)
    print("PANGENOME MAPPING SUMMARY")
    print("="*60)
    
    # Calculate totals including partial mappings
    genes_total_mapped = stats.genes_mapped + stats.genes_partial
    tx_total_mapped = stats.transcripts_mapped + stats.transcripts_partial
    
    genes_rate = genes_total_mapped / max(stats.total_genes, 1)
    tx_rate = tx_total_mapped / max(stats.total_transcripts, 1)
    
    print(f"\n{'Feature':<20} {'Total':>10} {'Mapped':>10} {'Rate':>10}")
    print("-"*50)
    
    # Genes
    if stats.genes_partial > 0:
        complete_rate = stats.genes_mapped / max(stats.total_genes, 1) * 100
        partial_rate = stats.genes_partial / max(stats.total_genes, 1) * 100
        rate_str = f"{genes_rate:>6.1%} [{complete_rate:.1f}% + {partial_rate:.1f}%]"
        print(f"{'Genes':<20} {stats.total_genes:>10} {genes_total_mapped:>10} {rate_str}")
    else:
        print(f"{'Genes':<20} {stats.total_genes:>10} {stats.genes_mapped:>10} {stats.gene_mapping_rate:>10.1%}")
    
    # Transcripts
    if stats.transcripts_partial > 0:
        complete_rate = stats.transcripts_mapped / max(stats.total_transcripts, 1) * 100
        partial_rate = stats.transcripts_partial / max(stats.total_transcripts, 1) * 100
        rate_str = f"{tx_rate:>6.1%} [{complete_rate:.1f}% + {partial_rate:.1f}%]"
        print(f"{'Transcripts':<20} {stats.total_transcripts:>10} {tx_total_mapped:>10} {rate_str}")
    else:
        print(f"{'Transcripts':<20} {stats.total_transcripts:>10} {stats.transcripts_mapped:>10} {stats.transcript_mapping_rate:>10.1%}")
    
    # Exons
    print(f"{'Exons':<20} {stats.total_exons:>10} {stats.exons_mapped:>10} {stats.exon_mapping_rate:>10.1%}")
    
    print(f"\n{'Validation':<20} {'Checked':>10} {'Valid':>10} {'Rate':>10}")
    print("-"*50)
    print(f"{'Splice sites':<20} {stats.splice_sites_checked:>10} {stats.splice_sites_valid:>10} {stats.splice_site_validation_rate:>10.1%}")
    print(f"{'Start codons':<20} {stats.start_codons_checked:>10} {stats.start_codons_valid:>10} {stats.start_codon_validation_rate:>10.1%}")
    print(f"{'Stop codons':<20} {stats.stop_codons_checked:>10} {stats.stop_codons_valid:>10} {stats.stop_codon_validation_rate:>10.1%}")
    
    print(f"\nMapping identity: {stats.mean_mapping_identity:.4f} (min: {stats.min_mapping_identity:.4f}, max: {stats.max_mapping_identity:.4f})")
    print(f"Total time: {stats.total_time_seconds:.1f}s")
    print("="*60 + "\n")
    
    # Print biotype group statistics if available
    if stats.biotype_group_stats:
        print("\n" + "="*60)
        print("MAPPING BY BIOTYPE GROUP")
        print("="*60)
        print(f"\n{'Biotype Group':<25} {'Genes':<20} {'Transcripts':<20}")
        print("-"*65)
        
        # Sort by gene count descending
        sorted_groups = sorted(
            stats.biotype_group_stats.items(),
            key=lambda x: x[1].total_genes,
            reverse=True
        )
        
        for group_key, stat in sorted_groups:
            gene_str = f"{stat.mapped_genes}/{stat.total_genes} ({stat.gene_mapping_rate:.1%})"
            tx_str = f"{stat.mapped_transcripts}/{stat.total_transcripts} ({stat.transcript_mapping_rate:.1%})"
            print(f"{stat.name:<25} {gene_str:<20} {tx_str:<20}")
        
        print("="*60 + "\n")
    
    # Print individual biotype statistics (top 15)
    if stats.biotype_stats:
        print("\n" + "="*60)
        print("MAPPING BY BIOTYPE (Top 15)")
        print("="*60)
        print(f"\n{'Biotype':<35} {'Genes':<18} {'Transcripts':<18}")
        print("-"*71)
        
        # Sort by gene count descending
        sorted_biotypes = sorted(
            stats.biotype_stats.items(),
            key=lambda x: x[1].total_genes,
            reverse=True
        )[:15]  # Show top 15
        
        for biotype_key, stat in sorted_biotypes:
            gene_str = f"{stat.mapped_genes}/{stat.total_genes} ({stat.gene_mapping_rate:.1%})"
            tx_str = f"{stat.mapped_transcripts}/{stat.total_transcripts} ({stat.transcript_mapping_rate:.1%})"
            print(f"{stat.name:<35} {gene_str:<18} {tx_str:<18}")
        
        # Show how many more biotypes there are
        remaining = len(stats.biotype_stats) - 15
        if remaining > 0:
            print(f"\n... and {remaining} more biotypes (see JSON report for full details)")
        
        print("="*60 + "\n")



def analyze_gene_order(
    original_genes: Dict[str, Gene],
    mapped_genes: Dict[str, Gene]
) -> GeneOrderAnalysis:
    """Analyze conservation of gene order between reference and target.
    
    Compares adjacent gene pairs on each chromosome to check if their
    relative order is preserved after mapping.
    
    Args:
        original_genes: Original reference genes (keyed by gene_id)
        mapped_genes: Mapped target genes (keyed by mapped_gene_id)
        
    Returns:
        GeneOrderAnalysis with order conservation statistics
    """
    analysis = GeneOrderAnalysis()
    
    # Build mapping from original ID to mapped gene
    orig_to_mapped: Dict[str, Gene] = {}
    for mapped_id, mapped_gene in mapped_genes.items():
        orig_id = mapped_gene.mapped_from
        if orig_id:
            orig_to_mapped[orig_id] = mapped_gene
    
    # Group original genes by chromosome, sorted by position
    genes_by_chr: Dict[str, List[Gene]] = defaultdict(list)
    for gene in original_genes.values():
        genes_by_chr[gene.seq_region].append(gene)
    
    for chr_genes in genes_by_chr.values():
        chr_genes.sort(key=lambda g: g.start)
    
    # Analyze adjacent pairs
    for chrom, genes in genes_by_chr.items():
        for i in range(len(genes) - 1):
            gene1 = genes[i]
            gene2 = genes[i + 1]
            
            analysis.total_pairs += 1
            
            # Check if both genes are mapped
            mapped1 = orig_to_mapped.get(gene1.feature_id)
            mapped2 = orig_to_mapped.get(gene2.feature_id)
            
            if mapped1 is None or mapped2 is None:
                analysis.unmappable_pairs += 1
                continue
            
            # Check if still on same chromosome
            if mapped1.seq_region != mapped2.seq_region:
                analysis.different_chromosomes += 1
                continue
            
            # Compare order: in reference, gene1 comes before gene2
            ref_order = gene1.start < gene2.start  # Should always be True
            target_order = mapped1.start < mapped2.start
            
            if ref_order == target_order:
                analysis.order_preserved += 1
            else:
                analysis.order_inverted += 1
                analysis.inversion_details.append({
                    "gene1_id": gene1.feature_id,
                    "gene1_name": gene1.gene_name,
                    "gene2_id": gene2.feature_id,
                    "gene2_name": gene2.gene_name,
                    "ref_positions": f"{gene1.start}-{gene1.end} < {gene2.start}-{gene2.end}",
                    "target_positions": f"{mapped1.start}-{mapped1.end} vs {mapped2.start}-{mapped2.end}"
                })
    
    logger.info(
        f"Gene order analysis: {analysis.order_preserved}/{analysis.total_pairs - analysis.unmappable_pairs} "
        f"pairs preserved ({analysis.order_conservation_rate:.1%}), "
        f"{analysis.order_inverted} inversions"
    )
    
    return analysis


def find_new_overlaps(
    original_genes: Dict[str, Gene],
    mapped_genes: Dict[str, Gene]
) -> List[OverlapConflict]:
    """Find genes that didn't overlap in reference but now overlap in target.
    
    Args:
        original_genes: Original reference genes
        mapped_genes: Mapped target genes
        
    Returns:
        List of OverlapConflict describing new overlaps
    """
    conflicts = []
    
    # Build mapping from original ID to mapped gene
    orig_to_mapped: Dict[str, Gene] = {}
    for mapped_id, mapped_gene in mapped_genes.items():
        orig_id = mapped_gene.mapped_from
        if orig_id:
            orig_to_mapped[orig_id] = mapped_gene
    
    # Group original genes by chromosome
    genes_by_chr: Dict[str, List[Gene]] = defaultdict(list)
    for gene in original_genes.values():
        if gene.feature_id in orig_to_mapped:  # Only include mapped genes
            genes_by_chr[gene.seq_region].append(gene)
    
    for chr_genes in genes_by_chr.values():
        chr_genes.sort(key=lambda g: g.start)
    
    # For each chromosome, check adjacent/nearby gene pairs
    for chrom, genes in genes_by_chr.items():
        for i in range(len(genes)):
            gene1 = genes[i]
            mapped1 = orig_to_mapped.get(gene1.feature_id)
            if not mapped1:
                continue
            
            # Check against genes within a window (adjacent genes + nearby)
            for j in range(i + 1, min(i + 50, len(genes))):  # Check up to 50 genes ahead
                gene2 = genes[j]
                mapped2 = orig_to_mapped.get(gene2.feature_id)
                if not mapped2:
                    continue
                
                # Skip if on different chromosomes in target
                if mapped1.seq_region != mapped2.seq_region:
                    continue
                
                # Check original overlap (should be non-overlapping)
                ref_overlap = max(0, min(gene1.end, gene2.end) - max(gene1.start, gene2.start) + 1)
                
                # Check target overlap
                target_overlap = max(0, min(mapped1.end, mapped2.end) - max(mapped1.start, mapped2.start) + 1)
                
                # If originally non-overlapping but now overlapping
                if ref_overlap == 0 and target_overlap > 0:
                    conflicts.append(OverlapConflict(
                        gene1_id=gene1.feature_id,
                        gene1_name=gene1.gene_name,
                        gene2_id=gene2.feature_id,
                        gene2_name=gene2.gene_name,
                        ref_gene1_start=gene1.start,
                        ref_gene1_end=gene1.end,
                        ref_gene2_start=gene2.start,
                        ref_gene2_end=gene2.end,
                        target_gene1_start=mapped1.start,
                        target_gene1_end=mapped1.end,
                        target_gene2_start=mapped2.start,
                        target_gene2_end=mapped2.end,
                        overlap_bp=target_overlap
                    ))
    
    logger.info(f"Found {len(conflicts)} new overlap conflicts")
    return conflicts


def analyze_synteny(
    original_genes: Dict[str, Gene],
    mapped_genes: Dict[str, Gene]
) -> SyntenyAnalysis:
    """Perform complete synteny analysis.
    
    Args:
        original_genes: Original reference genes
        mapped_genes: Mapped target genes
        
    Returns:
        SyntenyAnalysis with gene order and overlap conflict results
    """
    gene_order = analyze_gene_order(original_genes, mapped_genes)
    new_overlaps = find_new_overlaps(original_genes, mapped_genes)
    
    return SyntenyAnalysis(
        gene_order=gene_order,
        new_overlaps=new_overlaps
    )


def print_synteny_summary(analysis: SyntenyAnalysis):
    """Print synteny analysis summary."""
    go = analysis.gene_order
    
    print("\n" + "="*60)
    print("SYNTENY ANALYSIS")
    print("="*60)
    
    print(f"\nGene Order Conservation:")
    print(f"  Total adjacent pairs: {go.total_pairs}")
    print(f"  Order preserved:      {go.order_preserved} ({go.order_conservation_rate:.1%})")
    print(f"  Order inverted:       {go.order_inverted}")
    print(f"  Different chromosomes: {go.different_chromosomes}")
    print(f"  Unmappable pairs:     {go.unmappable_pairs}")
    
    print(f"\nNew Overlap Conflicts:")
    print(f"  Gene pairs now overlapping: {analysis.new_overlap_count}")
    
    if analysis.new_overlaps[:5]:
        print(f"\n  Top conflicts by overlap size:")
        sorted_conflicts = sorted(analysis.new_overlaps, key=lambda c: c.overlap_bp, reverse=True)
        for c in sorted_conflicts[:5]:
            print(f"    {c.gene1_name or c.gene1_id} / {c.gene2_name or c.gene2_id}: {c.overlap_bp}bp overlap")
    
    print("="*60 + "\n")


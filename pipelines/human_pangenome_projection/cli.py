#!/usr/bin/env python3
"""Command-line interface for the pangenome mapping system."""

import logging
import sys
from pathlib import Path
from typing import Optional, Set

import click

from pipeline import run_pipeline, PangenomeMappingPipeline
from src.config import MappingConfig


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


@click.group()
@click.version_option(version="0.1.0")
def cli():
    """Pangenome Annotation Mapping System.
    
    Maps genome annotations from a reference assembly to a target assembly
    using a three-stage approach: synteny detection, feature mapping,
    and conflict resolution.
    """
    pass


@cli.command()
@click.option("--ref-fasta", "-r", required=True, type=click.Path(exists=True),
              help="Reference genome FASTA file")
@click.option("--ref-gff", "-g", required=True, type=click.Path(exists=True),
              help="Reference annotation GFF3 file")
@click.option("--target-fasta", "-t", required=True, type=click.Path(exists=True),
              help="Target genome FASTA file")
@click.option("--output-gff", "-o", required=True, type=click.Path(),
              help="Output GFF3 file for mapped annotations")
@click.option("--output-stats", "-s", type=click.Path(),
              help="Output JSON file for mapping statistics")
@click.option("--chromosomes", "-c", multiple=True,
              help="Limit to specific chromosomes (can be specified multiple times)")
@click.option("--threads", "-j", default=4, type=int,
              help="Number of threads for alignment (default: 4)")
@click.option("--preset-profile", default="none", type=click.Choice(["none", "ultra-close", "less-close"]),
              help="Production parameter profile (default: none)")
@click.option("--min-identity", default=0.95, type=float,
              help="Minimum alignment identity for syntenic blocks (default: 0.95)")
@click.option("--min-block-length", default=10000, type=int,
              help="Minimum syntenic block length in bp (default: 10000)")
@click.option("--min-mapq", default=10, type=int,
              help="Minimum MAPQ for synteny alignments (default: 10)")
@click.option("--cnv-max-total-copies", default=3, type=int,
              help="Max total copies to emit per gene when expansion loci exist (default: 3)")
@click.option("--cnv-min-expansion-coverage", default=0.60, type=float,
              help="Minimum reference coverage for expansion loci (default: 0.60)")
@click.option("--cnv-min-score-ratio", default=0.85, type=float,
              help="Minimum expansion locus score ratio vs primary locus (default: 0.85)")
@click.option("--cnv-max-locus-overlap", default=0.30, type=float,
              help="Max reciprocal overlap allowed between selected expansion loci (default: 0.30)")
@click.option("--cnv-group-min-reciprocal-overlap", default=0.50, type=float,
              help="Min reciprocal overlap to treat expansion candidates as same target locus (default: 0.50)")
@click.option("--cnv-ambiguity-distance-bp", default=10000, type=int,
              help="Distance tolerance in bp for ambiguous expansion-source assignment (default: 10000)")
@click.option("--cnv-ambiguity-score-delta", default=0.01, type=float,
              help="Score tolerance for ambiguous expansion-source assignment (default: 0.01)")
@click.option("--disable-paralog-reassignment", is_flag=True,
              help="Disable synteny-aware reassignment of overlapping paralog clusters")
@click.option("--paralog-overlap-threshold", default=0.50, type=float,
              help="Reciprocal overlap threshold for paralog conflict groups (default: 0.50)")
@click.option("--paralog-max-locus-candidates", default=8, type=int,
              help="Maximum alternate loci tested per conflicting paralog (default: 8)")
@click.option("--paralog-max-rounds", default=2, type=int,
              help="Maximum reassignment rounds for paralog conflict resolution (default: 2)")
@click.option("--disable-missed-locus-recovery", is_flag=True,
              help="Disable second-pass local recovery for unmapped/partial loci")
@click.option("--missed-locus-max-genes", default=1000, type=int,
              help="Max genes considered in second-pass missed-locus recovery (default: 1000)")
@click.option("--missed-locus-search-padding-bp", default=200000, type=int,
              help="Reference/target padding around nearby synteny for missed-locus search (default: 200000)")
@click.option("--missed-locus-max-window-bp", default=800000, type=int,
              help="Maximum target window size for missed-locus search (default: 800000)")
@click.option("--missed-locus-min-identity", default=0.85, type=float,
              help="Minimum identity for accepting missed-locus local recovery (default: 0.85)")
@click.option("--missed-locus-min-coverage", default=0.60, type=float,
              help="Minimum coverage for accepting missed-locus local recovery (default: 0.60)")
@click.option("--disable-boundary-refinement", is_flag=True,
              help="Disable local base-level edge refinement for uncertain projections")
@click.option("--boundary-refine-min-coverage", default=0.98, type=float,
              help="Coverage threshold to trigger boundary refinement (default: 0.98)")
@click.option("--boundary-refine-window-bp", default=500, type=int,
              help="Window around projected edge for local refinement (default: 500)")
@click.option("--boundary-refine-anchor-bp", default=80, type=int,
              help="Reference anchor length for edge refinement (default: 80)")
@click.option("--boundary-refine-max-unaligned-bp", default=8, type=int,
              help="Max unclipped boundary bases allowed in refinement alignment (default: 8)")
@click.option("--disable-validation-rescue", is_flag=True,
              help="Disable validation-guided remapping rescue for uncertain genes")
@click.option("--validation-rescue-window-bp", default=1200, type=int,
              help="Window around projected edge during validation rescue remap (default: 1200)")
@click.option("--validation-rescue-anchor-bp", default=120, type=int,
              help="Anchor length used during validation rescue remap (default: 120)")
@click.option("--validation-rescue-max-genes", default=2000, type=int,
              help="Maximum genes to remap during validation rescue (default: 2000)")
@click.option("--disable-internal-realign-rescue", is_flag=True,
              help="Disable targeted internal realignment rescue for indel/strand continuity")
@click.option("--internal-realign-max-genes", default=3000, type=int,
              help="Maximum genes considered for targeted internal realignment rescue (default: 3000)")
@click.option("--internal-realign-ref-flank-bp", default=3000, type=int,
              help="Reference flank around gene for internal realignment windows (default: 3000)")
@click.option("--internal-realign-target-flank-bp", default=10000, type=int,
              help="Target flank around mapped locus for internal realignment windows (default: 10000)")
@click.option("--internal-realign-min-path-coverage", default=0.98, type=float,
              help="Minimum path coverage for accepting minimap2 internal realignment (default: 0.98)")
@click.option("--internal-realign-max-internal-gap-bp", default=6, type=int,
              help="Maximum uncovered internal reference gap for minimap2 path acceptance (default: 6)")
@click.option("--internal-realign-trigger-indel-jump-bp", default=6, type=int,
              help="Exon-offset jump threshold (bp) used to trigger internal realignment rescue (default: 6)")
@click.option("--internal-realign-fallback-backend", default="edlib", type=click.Choice(["edlib"]),
              help="Fallback backend when minimap2 path is fragmented (default: edlib)")
@click.option("--disable-internal-realign-accept-only-if-improved", is_flag=True,
              help="Allow internal realignment rescue updates even without validation-score improvement")
@click.option("--disable-splice-shift-rescue", is_flag=True,
              help="Disable splice micro-shift rescue for non-canonical introns")
@click.option("--splice-shift-max-bp", default=6, type=int,
              help="Maximum per-boundary splice rescue shift (default: 6)")
@click.option("--disable-codon-frame-rescue", is_flag=True,
              help="Disable CDS boundary codon/frame rescue for coding transcripts")
@click.option("--codon-rescue-max-bp", default=6, type=int,
              help="Maximum CDS boundary shift for codon/frame rescue (default: 6)")
@click.option("--disable-protein-qc", is_flag=True,
              help="Disable protein-level QC against reference coding transcripts")
@click.option("--protein-qc-min-identity", default=0.90, type=float,
              help="Protein identity threshold for QC pass classification (default: 0.90)")
@click.option("--protein-qc-min-coverage", default=0.90, type=float,
              help="Protein coverage threshold for QC pass classification (default: 0.90)")
@click.option("--allow-secondary/--primary-only", default=False,
              help="Include secondary synteny alignments (default: primary only)")
@click.option("--keep-temp", is_flag=True,
              help="Keep temporary files (PAF alignment)")
@click.option("--verbose", "-v", is_flag=True,
              help="Enable verbose logging")
def map(
    ref_fasta: str,
    ref_gff: str,
    target_fasta: str,
    output_gff: str,
    output_stats: Optional[str],
    chromosomes: tuple,
    threads: int,
    preset_profile: str,
    min_identity: float,
    min_block_length: int,
    min_mapq: int,
    cnv_max_total_copies: int,
    cnv_min_expansion_coverage: float,
    cnv_min_score_ratio: float,
    cnv_max_locus_overlap: float,
    cnv_group_min_reciprocal_overlap: float,
    cnv_ambiguity_distance_bp: int,
    cnv_ambiguity_score_delta: float,
    disable_paralog_reassignment: bool,
    paralog_overlap_threshold: float,
    paralog_max_locus_candidates: int,
    paralog_max_rounds: int,
    disable_missed_locus_recovery: bool,
    missed_locus_max_genes: int,
    missed_locus_search_padding_bp: int,
    missed_locus_max_window_bp: int,
    missed_locus_min_identity: float,
    missed_locus_min_coverage: float,
    disable_boundary_refinement: bool,
    boundary_refine_min_coverage: float,
    boundary_refine_window_bp: int,
    boundary_refine_anchor_bp: int,
    boundary_refine_max_unaligned_bp: int,
    disable_validation_rescue: bool,
    validation_rescue_window_bp: int,
    validation_rescue_anchor_bp: int,
    validation_rescue_max_genes: int,
    disable_internal_realign_rescue: bool,
    internal_realign_max_genes: int,
    internal_realign_ref_flank_bp: int,
    internal_realign_target_flank_bp: int,
    internal_realign_min_path_coverage: float,
    internal_realign_max_internal_gap_bp: int,
    internal_realign_trigger_indel_jump_bp: int,
    internal_realign_fallback_backend: str,
    disable_internal_realign_accept_only_if_improved: bool,
    disable_splice_shift_rescue: bool,
    splice_shift_max_bp: int,
    disable_codon_frame_rescue: bool,
    codon_rescue_max_bp: int,
    disable_protein_qc: bool,
    protein_qc_min_identity: float,
    protein_qc_min_coverage: float,
    allow_secondary: bool,
    keep_temp: bool,
    verbose: bool
):
    """Run the complete mapping pipeline.
    
    Example:
        python cli.py map \
            --ref-fasta GRCh38.fa \
            --ref-gff gencode.v44.gff3 \
            --target-fasta CHM13v2.fa \
            --output-gff chm13.mapped.gff3 \
            --output-stats stats.json \
            --chromosomes chr20 \
            --threads 8
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Convert chromosomes tuple to set
    chr_set = set(chromosomes) if chromosomes else None

    if preset_profile == "ultra-close":
        min_identity = max(min_identity, 0.98)
        cnv_min_expansion_coverage = max(cnv_min_expansion_coverage, 0.70)
        cnv_min_score_ratio = max(cnv_min_score_ratio, 0.90)
        paralog_overlap_threshold = max(paralog_overlap_threshold, 0.55)
        missed_locus_min_identity = max(missed_locus_min_identity, 0.90)
        missed_locus_min_coverage = max(missed_locus_min_coverage, 0.70)
    elif preset_profile == "less-close":
        min_identity = min(min_identity, 0.93)
        cnv_min_expansion_coverage = min(cnv_min_expansion_coverage, 0.55)
        cnv_min_score_ratio = min(cnv_min_score_ratio, 0.80)
        boundary_refine_window_bp = max(boundary_refine_window_bp, 800)
        validation_rescue_window_bp = max(validation_rescue_window_bp, 1600)
        missed_locus_search_padding_bp = max(missed_locus_search_padding_bp, 300000)
    
    config = MappingConfig(
        ref_fasta=Path(ref_fasta),
        ref_gff=Path(ref_gff),
        target_fasta=Path(target_fasta),
        output_gff=Path(output_gff),
        output_stats=Path(output_stats) if output_stats else None,
        chromosomes=chr_set,
        threads=threads,
        min_identity=min_identity,
        min_block_length=min_block_length,
        min_mapq=min_mapq,
        primary_only=not allow_secondary,
        cnv_max_total_copies=cnv_max_total_copies,
        cnv_min_expansion_coverage=cnv_min_expansion_coverage,
        cnv_min_score_ratio=cnv_min_score_ratio,
        cnv_max_locus_overlap=cnv_max_locus_overlap,
        cnv_group_min_reciprocal_overlap=cnv_group_min_reciprocal_overlap,
        cnv_ambiguity_distance_bp=cnv_ambiguity_distance_bp,
        cnv_ambiguity_score_delta=cnv_ambiguity_score_delta,
        enable_paralog_reassignment=not disable_paralog_reassignment,
        paralog_overlap_threshold=paralog_overlap_threshold,
        paralog_max_locus_candidates=paralog_max_locus_candidates,
        paralog_max_rounds=paralog_max_rounds,
        enable_missed_locus_recovery=not disable_missed_locus_recovery,
        missed_locus_max_genes=missed_locus_max_genes,
        missed_locus_search_padding_bp=missed_locus_search_padding_bp,
        missed_locus_max_window_bp=missed_locus_max_window_bp,
        missed_locus_min_identity=missed_locus_min_identity,
        missed_locus_min_coverage=missed_locus_min_coverage,
        enable_boundary_refinement=not disable_boundary_refinement,
        boundary_refine_min_coverage=boundary_refine_min_coverage,
        boundary_refine_window_bp=boundary_refine_window_bp,
        boundary_refine_anchor_bp=boundary_refine_anchor_bp,
        boundary_refine_max_unaligned_bp=boundary_refine_max_unaligned_bp,
        enable_validation_rescue=not disable_validation_rescue,
        validation_rescue_window_bp=validation_rescue_window_bp,
        validation_rescue_anchor_bp=validation_rescue_anchor_bp,
        validation_rescue_max_genes=validation_rescue_max_genes,
        enable_internal_realign_rescue=not disable_internal_realign_rescue,
        internal_realign_max_genes=internal_realign_max_genes,
        internal_realign_ref_flank_bp=internal_realign_ref_flank_bp,
        internal_realign_target_flank_bp=internal_realign_target_flank_bp,
        internal_realign_min_path_coverage=internal_realign_min_path_coverage,
        internal_realign_max_internal_gap_bp=internal_realign_max_internal_gap_bp,
        internal_realign_trigger_indel_jump_bp=internal_realign_trigger_indel_jump_bp,
        internal_realign_fallback_backend=internal_realign_fallback_backend,
        internal_realign_accept_only_if_improved=(
            not disable_internal_realign_accept_only_if_improved
        ),
        enable_splice_shift_rescue=not disable_splice_shift_rescue,
        splice_shift_max_bp=splice_shift_max_bp,
        enable_codon_frame_rescue=not disable_codon_frame_rescue,
        codon_rescue_max_bp=codon_rescue_max_bp,
        enable_protein_qc=not disable_protein_qc,
        protein_qc_min_identity=protein_qc_min_identity,
        protein_qc_min_coverage=protein_qc_min_coverage,
        keep_temp=keep_temp
    )
    
    try:
        pipeline = PangenomeMappingPipeline(config)
        stats = pipeline.run()
        
        click.echo(f"\nMapping complete!")
        click.echo(f"  Genes in output:    {stats.genes_in_output}/{stats.total_genes} ({stats.genes_in_output/stats.total_genes:.1%})")
        click.echo(f"  Transcripts output: {stats.transcripts_in_output}/{stats.total_transcripts} ({stats.transcripts_in_output/stats.total_transcripts:.1%})")
        if stats.genes_removed_by_refinement > 0:
            click.echo(f"  (Removed by conflict resolution: {stats.genes_removed_by_refinement})")
        click.echo(f"  Output: {output_gff}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


@cli.command()
@click.option("--ref-fasta", "-r", required=True, type=click.Path(exists=True),
              help="Reference genome FASTA file")
@click.option("--target-fasta", "-t", required=True, type=click.Path(exists=True),
              help="Target genome FASTA file")
@click.option("--output", "-o", required=True, type=click.Path(),
              help="Output PAF alignment file")
@click.option("--threads", "-j", default=4, type=int,
              help="Number of threads")
@click.option("--preset", default="asm5", type=click.Choice(["asm5", "asm10", "asm20"]),
              help="minimap2 preset (asm5=0.1%% div, asm10=1%% div)")
def synteny(
    ref_fasta: str,
    target_fasta: str,
    output: str,
    threads: int,
    preset: str
):
    """Run synteny detection only.
    
    Produces a PAF file with whole-genome alignments that can be used
    for feature projection.
    
    Example:
        python cli.py synteny \
            --ref-fasta GRCh38.fa \
            --target-fasta CHM13v2.fa \
            --output alignment.paf \
            --threads 8
    """
    from stages.synteny import run_minimap2, build_syntenic_blocks, get_synteny_statistics
    from src.models import SyntenicMap
    
    click.echo(f"Running synteny detection...")
    click.echo(f"  Reference: {ref_fasta}")
    click.echo(f"  Target: {target_fasta}")
    
    success = run_minimap2(ref_fasta, target_fasta, output, threads=threads, preset=preset)
    
    if success:
        click.echo(f"  Output: {output}")
        
        # Show statistics
        blocks = build_syntenic_blocks(output)
        synteny_map = SyntenicMap(blocks=blocks)
        stats = get_synteny_statistics(synteny_map)
        
        click.echo(f"\nSynteny Statistics:")
        click.echo(f"  Blocks: {stats['total_blocks']}")
        click.echo(f"  Mean identity: {stats.get('mean_identity', 0):.4f}")
        click.echo(f"  Inversions: {stats.get('inversions', 0)}")
    else:
        click.echo("Synteny detection failed!", err=True)
        sys.exit(1)


@cli.command()
@click.option("--gff", "-g", required=True, type=click.Path(exists=True),
              help="Mapped GFF3 file to validate")
@click.option("--target-fasta", "-t", required=True, type=click.Path(exists=True),
              help="Target genome FASTA file")
@click.option("--output", "-o", type=click.Path(),
              help="Output JSON validation report")
def validate(gff: str, target_fasta: str, output: Optional[str]):
    """Validate a mapped GFF3 file.
    
    Checks splice sites, start/stop codons, and structural integrity.
    
    Example:
        python cli.py validate \
            --gff mapped.gff3 \
            --target-fasta CHM13v2.fa \
            --output validation.json
    """
    import json
    from src.gff3_parser import parse_gff3
    from src.fasta_handler import FastaHandler
    from validation.structural import validate_mapped_genes
    
    click.echo(f"Validating {gff}...")
    
    genes = parse_gff3(gff)
    click.echo(f"  Loaded {len(genes)} genes")
    
    with FastaHandler(target_fasta) as fasta:
        results = validate_mapped_genes(genes, fasta)
    
    valid_count = sum(1 for r in results.values() if r.is_valid)
    total = len(results)
    
    click.echo(f"\nValidation Results:")
    click.echo(f"  Genes valid: {valid_count}/{total} ({valid_count/total:.1%})")
    
    # Count transcript-level stats
    tx_valid = sum(
        sum(1 for t in r.transcript_results if t.is_valid)
        for r in results.values()
    )
    tx_total = sum(len(r.transcript_results) for r in results.values())
    click.echo(f"  Transcripts valid: {tx_valid}/{tx_total} ({tx_valid/tx_total:.1%})")
    
    if output:
        report = {
            "genes_valid": valid_count,
            "genes_total": total,
            "transcripts_valid": tx_valid,
            "transcripts_total": tx_total,
            "details": {
                gid: {
                    "is_valid": r.is_valid,
                    "transcripts": [
                        {
                            "id": t.transcript_id,
                            "is_valid": t.is_valid,
                            "splice_sites_valid": f"{t.splice_sites_valid}/{t.splice_sites_checked}",
                            "start_codon_valid": t.start_codon_valid,
                            "stop_codon_valid": t.stop_codon_valid,
                            "errors": t.errors
                        }
                        for t in r.transcript_results
                    ]
                }
                for gid, r in results.items()
            }
        }
        with open(output, "w") as f:
            json.dump(report, f, indent=2)
        click.echo(f"  Report saved to {output}")


@cli.command()
@click.option("--stats", "-s", required=True, type=click.Path(exists=True),
              help="Statistics JSON file from mapping")
def summary(stats: str):
    """Display summary of mapping results.
    
    Example:
        python cli.py summary --stats stats.json
    """
    import json
    
    with open(stats) as f:
        data = json.load(f)
    
    click.echo("\n" + "="*60)
    click.echo("PANGENOME MAPPING SUMMARY")
    click.echo("="*60)
    
    if "summary" in data:
        s = data["summary"]
        click.echo(f"\nGenes:       {s['genes']['mapped']}/{s['genes']['total']} mapped ({s['genes']['mapping_rate']})")
        click.echo(f"Transcripts: {s['transcripts']['mapped']}/{s['transcripts']['total']} mapped ({s['transcripts']['mapping_rate']})")
        click.echo(f"Exons:       {s['exons']['mapped']}/{s['exons']['total']} mapped ({s['exons']['mapping_rate']})")
    
    if "validation" in data:
        v = data["validation"]
        click.echo(f"\nSplice sites: {v['splice_sites']['valid']}/{v['splice_sites']['checked']} valid ({v['splice_sites']['rate']})")
        click.echo(f"Start codons: {v['start_codons']['valid']}/{v['start_codons']['checked']} valid ({v['start_codons']['rate']})")
        click.echo(f"Stop codons:  {v['stop_codons']['valid']}/{v['stop_codons']['checked']} valid ({v['stop_codons']['rate']})")
    
    if "refinement" in data:
        r = data["refinement"]
        click.echo(f"\nConflicts resolved: {r.get('conflicts_resolved', 0)}")
        click.echo(f"Copy number expansions: {r.get('copy_number_expansions', 0)}")
        click.echo(f"Copy number contractions: {r.get('copy_number_contractions', 0)}")
    
    if "performance" in data:
        p = data["performance"]
        click.echo(f"\nTotal time: {p.get('total_time_seconds', 0):.1f}s")
    
    click.echo("="*60 + "\n")


if __name__ == "__main__":
    cli()

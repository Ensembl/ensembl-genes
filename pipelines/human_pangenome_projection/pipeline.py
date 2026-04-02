"""Main pipeline orchestration for pangenome mapping."""

import copy
import json
import logging
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from src.models import (
    Gene,
    GenomicInterval,
    Strand,
    SyntenicBlock,
    SyntenicMap,
    Transcript,
)
from src.config import (
    MappingConfig,
    MappingStatus,
    CANONICAL_DONOR_SITES,
    CANONICAL_ACCEPTOR_SITES,
)
from src.gff3_parser import parse_gff3
from src.gff3_writer import write_genes_to_gff3
from src.fasta_handler import FastaHandler

from stages.synteny import create_syntenic_map, get_synteny_statistics
from stages.mapping import (
    map_genes,
    GeneMappingResult,
    resolve_paralog_clusters_synteny_aware,
)
from stages.internal_realign import targeted_internal_realign
from stages.refinement import refine_mappings
from stages.gap_filling import fill_synteny_gaps
from stages.sex_chromosome import detect_sex_chromosomes, SexChromosomeMap

from validation.structural import (
    validate_mapped_genes, GeneValidation, StructuralValidator, TranscriptValidation
)
from validation.protein import GeneProteinQC, run_protein_qc
from validation.statistics import (
    MappingStatistics, calculate_statistics, generate_report, print_summary,
    analyze_synteny, print_synteny_summary, SyntenyAnalysis
)

logger = logging.getLogger(__name__)


class PangenomeMappingPipeline:
    """Main pipeline for pangenome annotation mapping."""
    
    def __init__(self, config: MappingConfig):
        """Initialize pipeline with configuration.
        
        Args:
            config: Pipeline configuration
        """
        self.config = config
        
        # Results storage
        self.original_genes: Dict[str, Gene] = {}
        self.syntenic_map: Optional[SyntenicMap] = None
        self.mapping_results: List[GeneMappingResult] = []
        self.mapped_genes: Dict[str, Gene] = {}
        self.validation_results: Dict[str, GeneValidation] = {}
        self.initial_validation_results: Dict[str, GeneValidation] = {}
        self.protein_qc_results: Dict[str, GeneProteinQC] = {}
        self.protein_qc_report: Dict = {}
        self.audit_traces: Dict = {}
        self.refinement_report: Dict = {}
        self.paralog_resolution_report: Dict = {}
        self.stats: Optional[MappingStatistics] = None
        self.synteny_analysis: Optional[SyntenyAnalysis] = None
        self.pre_refinement_genes: Dict[str, Gene] = {}  # Genes before refinement filtering
        
        # Timing
        self.timings: Dict[str, float] = {}
        
        # Sex chromosome map (initialized during synteny stage)
        self.sex_chrom_map: Optional[SexChromosomeMap] = None
    
    def run(self) -> MappingStatistics:
        """Run the complete mapping pipeline.
        
        Returns:
            MappingStatistics with results
        """
        start_time = time.time()
        
        logger.info("="*60)
        logger.info("Starting Pangenome Mapping Pipeline")
        logger.info("="*60)
        
        # Stage 0: Load input data
        self._load_inputs()
        
        # Stage 1: Synteny detection
        self._run_synteny_detection()
        
        # Sex chromosome detection (after synteny, before mapping)
        self._detect_sex_chromosomes()
        
        # Stage 2: Feature mapping
        self._run_feature_mapping()
        
        # Validation
        self._run_validation()
        
        # Stage 3: Refinement
        self._run_refinement()
        
        # Calculate statistics
        self._calculate_statistics()
        
        # Synteny analysis (gene order, new overlaps)
        self._run_synteny_analysis()
        
        # Final timing must be set before outputs so reports include total runtime
        total_time = time.time() - start_time
        self.stats.total_time_seconds = total_time
        
        # Generate outputs
        self._generate_outputs()
        
        logger.info(f"Pipeline completed in {total_time:.1f} seconds")
        
        return self.stats
    
    def _load_inputs(self):
        """Load reference GFF3 and validate inputs."""
        logger.info("Loading input data...")
        
        # Validate input files exist
        if not self.config.ref_fasta.exists():
            raise FileNotFoundError(f"Reference FASTA not found: {self.config.ref_fasta}")
        if not self.config.ref_gff.exists():
            raise FileNotFoundError(f"Reference GFF3 not found: {self.config.ref_gff}")
        if not self.config.target_fasta.exists():
            raise FileNotFoundError(f"Target FASTA not found: {self.config.target_fasta}")
        
        # Parse reference GFF3
        self.original_genes = parse_gff3(
            self.config.ref_gff,
            chromosomes=self.config.chromosomes,
            biotypes=self.config.biotypes
        )
        
        logger.info(f"Loaded {len(self.original_genes)} genes from reference")
    
    def _run_synteny_detection(self):
        """Run Stage 1: Synteny detection."""
        logger.info("Stage 1: Synteny Detection")
        start = time.time()
        
        # Use output directory for temp files so they persist for caching
        output_dir = self.config.output_gff.parent if self.config.output_gff else None
        paf_output = output_dir / f"synteny_{self.config.ref_fasta.stem}_{self.config.target_fasta.stem}.paf" if output_dir else None
        
        self.syntenic_map = create_syntenic_map(
            ref_fasta=self.config.ref_fasta,
            target_fasta=self.config.target_fasta,
            threads=self.config.threads,
            min_identity=self.config.min_identity,
            min_block_length=self.config.min_block_length,
            min_mapq=self.config.min_mapq,
            primary_only=self.config.primary_only,
            chromosomes=self.config.chromosomes,
            temp_dir=output_dir,
            keep_paf=True,  # Keep PAF for potential reruns
            paf_output=paf_output
        )
        
        self.timings["synteny"] = time.time() - start
        
        # Log synteny statistics
        synteny_stats = get_synteny_statistics(self.syntenic_map)
        logger.info(f"Created {synteny_stats['total_blocks']} syntenic blocks")
        logger.info(f"Mean identity: {synteny_stats.get('mean_identity', 0):.4f}")
        logger.info(f"Total reference coverage: {synteny_stats.get('total_ref_coverage_bp', 0):,} bp")
        
        # Fill gaps for genes in regions without synteny coverage
        self.syntenic_map, gaps_filled = fill_synteny_gaps(
            self.original_genes,
            self.syntenic_map,
            self.config.ref_fasta,
            self.config.target_fasta,
            threads=self.config.threads,
            temp_dir=output_dir
        )
        if gaps_filled > 0:
            logger.info(f"Filled {gaps_filled} synteny gaps via local alignment")
    
    def _detect_sex_chromosomes(self):
        """Detect sex chromosomes in target genome for mapping constraints."""
        logger.info("Detecting sex chromosomes in target genome...")
        start = time.time()
        
        self.sex_chrom_map = detect_sex_chromosomes(
            self.config.target_fasta,
            threads=self.config.threads
        )
        
        self.timings["sex_detection"] = time.time() - start
        
        if self.sex_chrom_map.has_x or self.sex_chrom_map.has_y:
            logger.info(f"Sex chromosome detection: X={self.sex_chrom_map.has_x}, Y={self.sex_chrom_map.has_y}")
            logger.info(f"  Detection method: {self.sex_chrom_map.detection_method}")
        else:
            logger.info("No sex chromosomes detected - all mappings allowed")
    
    def _run_feature_mapping(self):
        """Run Stage 2: Feature mapping."""
        logger.info("Stage 2: Feature Mapping")
        start = time.time()
        
        def progress(current, total):
            if current % 500 == 0:
                logger.info(f"  Mapped {current}/{total} genes...")
        
        with FastaHandler(self.config.ref_fasta) as ref_fasta, FastaHandler(self.config.target_fasta) as target_fasta:
            self.mapped_genes, self.mapping_results = map_genes(
                self.original_genes,
                self.syntenic_map,
                progress_callback=progress,
                sex_chrom_map=self.sex_chrom_map,
                cnv_max_total_copies=self.config.cnv_max_total_copies,
                cnv_min_expansion_coverage=self.config.cnv_min_expansion_coverage,
                cnv_min_score_ratio=self.config.cnv_min_score_ratio,
                cnv_max_locus_overlap=self.config.cnv_max_locus_overlap,
                cnv_group_min_reciprocal_overlap=self.config.cnv_group_min_reciprocal_overlap,
                cnv_ambiguity_distance_bp=self.config.cnv_ambiguity_distance_bp,
                cnv_ambiguity_score_delta=self.config.cnv_ambiguity_score_delta,
                ref_fasta_handler=ref_fasta,
                target_fasta_handler=target_fasta,
                enable_boundary_refinement=self.config.enable_boundary_refinement,
                boundary_refine_min_coverage=self.config.boundary_refine_min_coverage,
                boundary_refine_window_bp=self.config.boundary_refine_window_bp,
                boundary_refine_anchor_bp=self.config.boundary_refine_anchor_bp,
                boundary_refine_max_unaligned_bp=self.config.boundary_refine_max_unaligned_bp,
            )

            if self.config.enable_paralog_reassignment:
                self.mapped_genes, self.paralog_resolution_report = (
                    resolve_paralog_clusters_synteny_aware(
                        self.mapped_genes,
                        self.original_genes,
                        self.syntenic_map,
                        sex_chrom_map=self.sex_chrom_map,
                        ref_fasta_handler=ref_fasta,
                        target_fasta_handler=target_fasta,
                        enable_boundary_refinement=self.config.enable_boundary_refinement,
                        boundary_refine_min_coverage=self.config.boundary_refine_min_coverage,
                        boundary_refine_window_bp=self.config.boundary_refine_window_bp,
                        boundary_refine_anchor_bp=self.config.boundary_refine_anchor_bp,
                        boundary_refine_max_unaligned_bp=self.config.boundary_refine_max_unaligned_bp,
                        overlap_threshold=self.config.paralog_overlap_threshold,
                        max_locus_candidates=self.config.paralog_max_locus_candidates,
                        max_rounds=self.config.paralog_max_rounds,
                    )
                )
                if self.paralog_resolution_report.get("groups_examined", 0) > 0:
                    logger.info(
                        "Paralog resolver: groups=%s remapped=%s conflicts %s -> %s",
                        self.paralog_resolution_report.get("groups_examined", 0),
                        self.paralog_resolution_report.get("genes_reassigned", 0),
                        self.paralog_resolution_report.get("conflicts_before", 0),
                        self.paralog_resolution_report.get("conflicts_after", 0),
                    )
            if self.config.enable_missed_locus_recovery:
                recovered = self._run_missed_locus_recovery(ref_fasta, target_fasta)
                if recovered > 0:
                    logger.info(f"Missed-locus recovery added {recovered} gene(s)")
        
        self.timings["mapping"] = time.time() - start
        logger.info(f"Mapped {len(self.mapped_genes)} genes")

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Reverse complement for DNA sequences."""
        table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq.translate(table)[::-1]

    @classmethod
    def _find_anchor_interval_in_window(
        cls,
        ref_seq: str,
        target_window: str,
        window_start: int,
        expected_center: int,
        anchor_len: int
    ) -> Optional[Tuple[int, int, bool]]:
        """Find likely locus interval by anchor placement in one target window."""
        if len(ref_seq) < anchor_len * 2:
            return None

        best = None

        def consider(query_seq: str, is_inverted: bool):
            nonlocal best
            left_anchor = query_seq[:anchor_len]
            right_anchor = query_seq[-anchor_len:]
            left_pos = target_window.find(left_anchor)
            right_pos = target_window.rfind(right_anchor)
            if left_pos < 0 or right_pos < 0 or right_pos < left_pos:
                return
            cand_start = window_start + left_pos
            cand_end = window_start + right_pos + anchor_len - 1
            cand_center = (cand_start + cand_end) // 2
            penalty = abs((cand_end - cand_start + 1) - len(ref_seq)) + abs(cand_center - expected_center)
            if best is None or penalty < best[0]:
                best = (penalty, cand_start, cand_end, is_inverted)

        consider(ref_seq, False)
        consider(cls._reverse_complement(ref_seq), True)

        if best is None:
            return None
        return best[1], best[2], best[3]

    def _recover_unmapped_gene_local_search(
        self,
        gene: Gene,
        ref_fasta: FastaHandler,
        target_fasta: FastaHandler,
        mapper
    ) -> Optional[Gene]:
        """Recover one unmapped gene via anchor-based local search in syntenic gaps."""
        search_padding = self.config.missed_locus_search_padding_bp
        max_window = self.config.missed_locus_max_window_bp
        nearby_blocks = self.syntenic_map.find_blocks_for_ref_region(
            gene.seq_region,
            max(1, gene.start - search_padding),
            gene.end + search_padding,
        )
        if not nearby_blocks:
            return None

        ref_seq = ref_fasta.fetch(gene.seq_region, gene.start, gene.end)
        if not ref_seq or len(ref_seq) < 24:
            return None

        anchor_len = min(120, max(24, len(ref_seq) // 4))

        block_groups: Dict[Tuple[str, bool], List[SyntenicBlock]] = {}
        for block in nearby_blocks:
            key = (block.target_interval.seq_region, block.is_inverted)
            block_groups.setdefault(key, []).append(block)

        for (target_chr, _), blocks in sorted(
            block_groups.items(),
            key=lambda item: len(item[1]),
            reverse=True,
        ):
            target_len = target_fasta.get_length(target_chr)
            if target_len <= 0:
                continue

            raw_start = min(b.target_interval.start for b in blocks)
            raw_end = max(b.target_interval.end for b in blocks)
            expected_center = (raw_start + raw_end) // 2
            window_start = max(1, raw_start - search_padding)
            window_end = min(target_len, raw_end + search_padding)
            if window_end - window_start + 1 > max_window:
                half = max_window // 2
                window_start = max(1, expected_center - half)
                window_end = min(target_len, expected_center + half)
            if window_start >= window_end:
                continue

            target_window = target_fasta.fetch(target_chr, window_start, window_end)
            interval = self._find_anchor_interval_in_window(
                ref_seq,
                target_window,
                window_start,
                expected_center,
                anchor_len=anchor_len,
            )
            if interval is None:
                continue
            cand_start, cand_end, is_inverted = interval
            if cand_start >= cand_end:
                continue

            strand = Strand.MINUS if is_inverted else Strand.PLUS
            oriented_target_seq = target_fasta.fetch(
                target_chr,
                cand_start,
                cand_end,
                strand=strand,
            )
            compare_len = min(len(ref_seq), len(oriented_target_seq))
            if compare_len <= 0:
                continue
            matches = sum(
                1 for i in range(compare_len)
                if ref_seq[i] == oriented_target_seq[i]
            )
            identity = matches / compare_len
            coverage = compare_len / len(ref_seq)
            if identity < self.config.missed_locus_min_identity:
                continue
            if coverage < self.config.missed_locus_min_coverage:
                continue

            synthetic_block = SyntenicBlock(
                ref_interval=GenomicInterval(
                    seq_region=gene.seq_region,
                    start=gene.start,
                    end=gene.end,
                    strand=gene.strand,
                ),
                target_interval=GenomicInterval(
                    seq_region=target_chr,
                    start=cand_start,
                    end=cand_end,
                    strand=strand,
                ),
                identity=identity,
                alignment_length=compare_len,
                matches=matches,
                mismatches=compare_len - matches,
                block_id=f"missed_locus_{gene.feature_id}",
            )
            recovery_result = mapper._map_gene_single_copy(
                gene,
                mapped_gene_id=f"mapped_{gene.feature_id}",
                candidate_blocks=[synthetic_block],
                transcript_id_suffix="",
            )
            if recovery_result.mapped is None:
                continue

            recovered_gene = recovery_result.mapped
            recovered_gene.attributes["missed_locus_recovered"] = "true"
            recovered_gene.attributes["missed_locus_recovery_strategy"] = "anchor_local_search"
            recovered_gene.attributes["missed_locus_identity"] = f"{identity:.4f}"
            recovered_gene.attributes["missed_locus_coverage"] = f"{coverage:.4f}"
            return recovered_gene

        return None

    def _run_missed_locus_recovery(
        self,
        ref_fasta: FastaHandler,
        target_fasta: FastaHandler
    ) -> int:
        """Second-pass recovery for unmapped/low-confidence genes near syntenic gaps."""
        from stages.mapping import FeatureMapper

        candidates = [
            r for r in self.mapping_results
            if r.status != MappingStatus.MAPPED
        ]
        if not candidates:
            return 0
        if self.config.missed_locus_max_genes > 0:
            candidates = candidates[: self.config.missed_locus_max_genes]

        mapper = FeatureMapper(
            self.syntenic_map,
            sex_chrom_map=self.sex_chrom_map,
            cnv_max_total_copies=self.config.cnv_max_total_copies,
            cnv_min_expansion_coverage=self.config.cnv_min_expansion_coverage,
            cnv_min_score_ratio=self.config.cnv_min_score_ratio,
            cnv_max_locus_overlap=self.config.cnv_max_locus_overlap,
            ref_fasta_handler=ref_fasta,
            target_fasta_handler=target_fasta,
            enable_boundary_refinement=self.config.enable_boundary_refinement,
            boundary_refine_min_coverage=self.config.boundary_refine_min_coverage,
            boundary_refine_window_bp=self.config.boundary_refine_window_bp,
            boundary_refine_anchor_bp=self.config.boundary_refine_anchor_bp,
            boundary_refine_max_unaligned_bp=self.config.boundary_refine_max_unaligned_bp,
        )

        recovered = 0
        for mapping_result in candidates:
            original_gene = mapping_result.original
            mapped_id = f"mapped_{original_gene.feature_id}"
            if mapped_id in self.mapped_genes:
                continue

            rescued_gene = self._recover_unmapped_gene_local_search(
                original_gene,
                ref_fasta,
                target_fasta,
                mapper,
            )
            if rescued_gene is None:
                continue

            self.mapped_genes[mapped_id] = rescued_gene
            mapping_result.mapped = rescued_gene
            mapping_result.status = rescued_gene.mapping_status or MappingStatus.PARTIAL
            mapping_result.errors.append(
                "Recovered by second-pass missed-locus local search"
            )
            recovered += 1

        return recovered
    
    def _run_validation(self):
        """Run structural validation."""
        logger.info("Validating mapped features...")
        start = time.time()
        
        with FastaHandler(self.config.target_fasta) as fasta:
            self.validation_results = validate_mapped_genes(
                self.mapped_genes, fasta
            )
            self.initial_validation_results = copy.deepcopy(self.validation_results)
            validation_updated = False
            if self.config.enable_validation_rescue:
                rescued = self._run_validation_guided_rescue(fasta)
                if rescued > 0:
                    validation_updated = True
                    logger.info(f"Validation rescue updated {rescued} gene(s)")
            if self.config.enable_internal_realign_rescue:
                realign_rescued = self._run_internal_realign_rescue(fasta)
                if realign_rescued > 0:
                    validation_updated = True
                    logger.info(
                        f"Internal realignment rescue updated {realign_rescued} gene(s)"
                    )
            if self.config.enable_splice_shift_rescue:
                splice_rescued = self._run_splice_shift_rescue(fasta)
                if splice_rescued > 0:
                    validation_updated = True
                    logger.info(
                        f"Splice shift rescue updated {splice_rescued} gene(s)"
                    )
            if self.config.enable_codon_frame_rescue:
                codon_rescued = self._run_codon_frame_rescue(fasta)
                if codon_rescued > 0:
                    validation_updated = True
                    logger.info(
                        f"Codon/frame rescue updated {codon_rescued} gene(s)"
                    )
            if validation_updated:
                self.validation_results = validate_mapped_genes(
                    self.mapped_genes, fasta
                )
                logger.info("Recomputed validation metrics after rescue updates")

            if self.config.enable_protein_qc:
                with FastaHandler(self.config.ref_fasta) as ref_fasta:
                    self.protein_qc_results = run_protein_qc(
                        self.original_genes,
                        self.mapped_genes,
                        ref_fasta,
                        fasta,
                        min_identity=self.config.protein_qc_min_identity,
                        min_coverage=self.config.protein_qc_min_coverage,
                    )
                self._annotate_protein_qc_confidence()
        
        self.timings["validation"] = time.time() - start
        
        valid_count = sum(1 for v in self.validation_results.values() if v.is_valid)
        logger.info(f"Validation: {valid_count}/{len(self.validation_results)} genes fully valid")

    @staticmethod
    def _classify_gene_confidence(
        structural_valid: bool,
        protein_qc: Optional[GeneProteinQC]
    ) -> str:
        """Combine structural/protein signals into a confidence class."""
        if protein_qc is None:
            return "high" if structural_valid else "low"
        if structural_valid and protein_qc.status == "ok":
            return "high"
        if not structural_valid and protein_qc.status == "ok":
            return "method_error_suspected"
        if structural_valid and protein_qc.status in {"diverged", "internal_stop"}:
            return "biological_divergence_suspected"
        if not structural_valid and protein_qc.status in {"diverged", "internal_stop"}:
            return "low"
        return "medium"

    def _annotate_protein_qc_confidence(self):
        """Attach protein QC metrics and confidence classes to mapped features."""
        genes_checked = len(self.protein_qc_results)
        genes_ok = sum(1 for q in self.protein_qc_results.values() if q.status == "ok")
        genes_diverged = sum(
            1 for q in self.protein_qc_results.values() if q.status == "diverged"
        )
        genes_internal_stop = sum(
            1 for q in self.protein_qc_results.values() if q.status == "internal_stop"
        )

        for mapped_id, gene in self.mapped_genes.items():
            gene_validation = self.validation_results.get(mapped_id)
            protein_qc = self.protein_qc_results.get(mapped_id)
            confidence = self._classify_gene_confidence(
                gene_validation.is_valid if gene_validation else False,
                protein_qc,
            )
            gene.attributes["projection_confidence_class"] = confidence
            if protein_qc is None:
                continue
            gene.attributes["protein_qc_status"] = protein_qc.status
            gene.attributes["protein_qc_mean_identity"] = f"{protein_qc.mean_identity:.4f}"
            gene.attributes["protein_qc_mean_coverage"] = f"{protein_qc.mean_coverage:.4f}"
            gene.attributes["protein_qc_internal_stop_transcripts"] = str(
                protein_qc.transcripts_with_internal_stop
            )
            tx_by_id = {t.feature_id: t for t in gene.transcripts}
            for tx_result in protein_qc.transcript_results:
                tx = tx_by_id.get(tx_result.transcript_id)
                if tx is None:
                    continue
                tx.attributes["protein_qc_status"] = tx_result.status
                tx.attributes["protein_qc_identity"] = f"{tx_result.identity:.4f}"
                tx.attributes["protein_qc_coverage"] = f"{tx_result.coverage:.4f}"
                if tx_result.internal_stop:
                    tx.attributes["protein_qc_internal_stop"] = "true"

        self.protein_qc_report = {
            "genes_checked": genes_checked,
            "genes_ok": genes_ok,
            "genes_diverged": genes_diverged,
            "genes_internal_stop": genes_internal_stop,
            "genes": {
                gene_id: qc.to_dict()
                for gene_id, qc in self.protein_qc_results.items()
            },
        }
        logger.info(
            "Protein QC: checked=%s ok=%s diverged=%s internal_stop=%s",
            genes_checked,
            genes_ok,
            genes_diverged,
            genes_internal_stop,
        )

    @staticmethod
    def _gene_has_uncertain_projection(gene: Gene) -> bool:
        """Heuristic for method-induced boundary uncertainty."""
        if gene.mapping_status and gene.mapping_status != "mapped":
            return True
        
        def has_uncertain_attrs(feature) -> bool:
            attrs = feature.attributes or {}
            if attrs.get("projection_clamped_5prime") == "true":
                return True
            if attrs.get("projection_clamped_3prime") == "true":
                return True
            cov = attrs.get("projection_coverage")
            if cov is not None:
                try:
                    if float(cov) < 1.0:
                        return True
                except (TypeError, ValueError):
                    pass
            return False
        
        for tx in gene.transcripts:
            if tx.mapping_status and tx.mapping_status != "mapped":
                return True
            if has_uncertain_attrs(tx):
                return True
            for exon in tx.exons:
                if has_uncertain_attrs(exon):
                    return True
            for cds in tx.cds_list:
                if has_uncertain_attrs(cds):
                    return True
            for utr in tx.utrs:
                if has_uncertain_attrs(utr):
                    return True
        return False

    @staticmethod
    def _validation_score(gene_validation: GeneValidation) -> Tuple[float, float, float, float, float, int]:
        """Lexicographic score for comparing two structural validation results."""
        tx_checked = max(1, gene_validation.transcripts_checked)
        tx_score = gene_validation.transcripts_valid / tx_checked
        
        splice_checked = sum(t.splice_sites_checked for t in gene_validation.transcript_results)
        splice_valid = sum(t.splice_sites_valid for t in gene_validation.transcript_results)
        splice_score = splice_valid / max(1, splice_checked)
        
        start_checked = sum(1 for t in gene_validation.transcript_results if t.start_codon_result)
        start_valid = sum(
            1 for t in gene_validation.transcript_results
            if t.start_codon_result and t.start_codon_valid
        )
        start_score = start_valid / max(1, start_checked)
        
        stop_checked = sum(1 for t in gene_validation.transcript_results if t.stop_codon_result)
        stop_valid = sum(
            1 for t in gene_validation.transcript_results
            if t.stop_codon_result and t.stop_codon_valid
        )
        stop_score = stop_valid / max(1, stop_checked)
        
        frame_checked = sum(1 for t in gene_validation.transcript_results if t.cds_length > 0)
        frame_valid = sum(
            1 for t in gene_validation.transcript_results
            if t.cds_length == 0 or t.cds_length_valid
        )
        frame_score = frame_valid / max(1, frame_checked)
        
        total_errors = sum(len(t.errors) for t in gene_validation.transcript_results)
        return (
            tx_score,
            splice_score,
            start_score,
            stop_score,
            frame_score,
            -total_errors,
        )

    def _run_validation_guided_rescue(self, target_fasta: FastaHandler) -> int:
        """Remap uncertain/invalid genes with stricter edge refinement and keep improvements."""
        if self.config.validation_rescue_max_genes == 0:
            return 0
        
        candidate_pairs: List[Tuple[str, str]] = []
        for mapped_id, mapped_gene in self.mapped_genes.items():
            if mapped_gene.attributes.get("cnv_copy_role") == "expansion":
                continue
            gene_val = self.validation_results.get(mapped_id)
            if not gene_val or gene_val.is_valid:
                continue
            orig_id = mapped_gene.mapped_from
            if not orig_id or orig_id not in self.original_genes:
                continue
            if not self._gene_has_uncertain_projection(mapped_gene):
                continue
            candidate_pairs.append((mapped_id, orig_id))
        
        if not candidate_pairs:
            return 0
        
        if len(candidate_pairs) > self.config.validation_rescue_max_genes:
            candidate_pairs = candidate_pairs[:self.config.validation_rescue_max_genes]
        
        logger.info(
            f"Validation rescue: attempting remap for {len(candidate_pairs)} uncertain/invalid gene(s)"
        )
        
        rescued = 0
        validator = StructuralValidator(target_fasta)
        with FastaHandler(self.config.ref_fasta) as ref_fasta:
            from stages.mapping import FeatureMapper
            
            rescue_mapper = FeatureMapper(
                self.syntenic_map,
                sex_chrom_map=self.sex_chrom_map,
                cnv_max_total_copies=self.config.cnv_max_total_copies,
                cnv_min_expansion_coverage=self.config.cnv_min_expansion_coverage,
                cnv_min_score_ratio=self.config.cnv_min_score_ratio,
                cnv_max_locus_overlap=self.config.cnv_max_locus_overlap,
                ref_fasta_handler=ref_fasta,
                target_fasta_handler=target_fasta,
                enable_boundary_refinement=True,
                boundary_refine_min_coverage=1.0,
                boundary_refine_window_bp=max(
                    self.config.boundary_refine_window_bp,
                    self.config.validation_rescue_window_bp
                ),
                boundary_refine_anchor_bp=max(
                    self.config.boundary_refine_anchor_bp,
                    self.config.validation_rescue_anchor_bp
                ),
                boundary_refine_max_unaligned_bp=self.config.boundary_refine_max_unaligned_bp,
            )
            
            for mapped_id, orig_id in candidate_pairs:
                old_gene = self.mapped_genes.get(mapped_id)
                old_val = self.validation_results.get(mapped_id)
                if old_gene is None or old_val is None:
                    continue
                
                remap = rescue_mapper.map_gene(self.original_genes[orig_id])
                if remap.mapped is None:
                    continue
                
                new_gene = remap.mapped
                new_val = validator.validate_gene(new_gene)
                
                old_score = self._validation_score(old_val)
                new_score = self._validation_score(new_val)
                if new_score > old_score:
                    new_gene.feature_id = mapped_id
                    new_gene.attributes["rescue_validation_remap"] = "true"
                    self.mapped_genes[mapped_id] = new_gene
                    self.validation_results[mapped_id] = new_val
                    rescued += 1
        
        return rescued

    @staticmethod
    def _safe_attr_float(attrs: Dict, key: str, default: float = 0.0) -> float:
        value = attrs.get(key)
        if value is None:
            return default
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    def _has_sharp_protein_inconsistency(self, mapped_gene: Gene) -> bool:
        """Detect likely projection errors from nucleotide/protein disagreement."""
        attrs = mapped_gene.attributes or {}
        nucleotide_identity = mapped_gene.mapping_identity or self._safe_attr_float(
            attrs, "mapping_identity", default=0.0
        )
        protein_identity = self._safe_attr_float(
            attrs, "protein_qc_mean_identity", default=1.0
        )
        protein_coverage = self._safe_attr_float(
            attrs, "protein_qc_mean_coverage", default=1.0
        )
        return (
            nucleotide_identity >= 0.98
            and (protein_identity < 0.40 or protein_coverage < 0.60)
        )

    def _has_exon_offset_discontinuity(
        self,
        original_gene: Gene,
        mapped_gene: Gene,
        jump_bp: int,
    ) -> bool:
        """Detect abrupt exon-offset jumps suggestive of missed internal indels."""
        if jump_bp <= 0:
            return False

        original_txs = {tx.feature_id: tx for tx in original_gene.transcripts}
        for mapped_tx in mapped_gene.transcripts:
            orig_tx_id = mapped_tx.mapped_from
            if not orig_tx_id or orig_tx_id not in original_txs:
                continue
            orig_tx = original_txs[orig_tx_id]

            mapped_by_source = {
                exon.mapped_from: exon for exon in mapped_tx.exons if exon.mapped_from
            }
            ordered_orig = sorted(orig_tx.exons, key=lambda e: e.start)
            start_deltas: List[int] = []
            end_deltas: List[int] = []
            for orig_exon in ordered_orig:
                mapped_exon = mapped_by_source.get(orig_exon.feature_id)
                if mapped_exon is None:
                    continue
                start_deltas.append(mapped_exon.start - orig_exon.start)
                end_deltas.append(mapped_exon.end - orig_exon.end)

            if len(start_deltas) < 2:
                continue

            for idx in range(1, len(start_deltas)):
                if abs(start_deltas[idx] - start_deltas[idx - 1]) >= jump_bp:
                    return True
                if abs(end_deltas[idx] - end_deltas[idx - 1]) >= jump_bp:
                    return True

        return False

    def _run_internal_realign_rescue(
        self,
        target_fasta: FastaHandler,
    ) -> int:
        """Targeted internal realignment rescue for large-indel continuity issues."""
        if self.config.internal_realign_max_genes == 0:
            return 0

        candidate_pairs: List[Tuple[str, str]] = []
        trigger_counts: Counter = Counter()
        valid_uncertain_skipped = 0
        for mapped_id, mapped_gene in self.mapped_genes.items():
            if mapped_gene.attributes.get("cnv_copy_role") == "expansion":
                continue
            orig_id = mapped_gene.mapped_from
            if not orig_id or orig_id not in self.original_genes:
                continue

            old_val = self.validation_results.get(mapped_id)
            is_invalid = old_val is not None and not old_val.is_valid
            uncertain_projection = self._gene_has_uncertain_projection(mapped_gene)
            indel_jump = self._has_exon_offset_discontinuity(
                self.original_genes[orig_id],
                mapped_gene,
                self.config.internal_realign_trigger_indel_jump_bp,
            )
            protein_inconsistency = self._has_sharp_protein_inconsistency(mapped_gene)
            orientation_ambiguous = (
                (mapped_gene.attributes or {}).get("orientation_confidence") == "ambiguous"
            )
            should_trigger = False
            if is_invalid:
                trigger_counts["invalid"] += 1
                if uncertain_projection:
                    trigger_counts["invalid_uncertain"] += 1
                should_trigger = True
            if indel_jump:
                trigger_counts["indel_jump"] += 1
                should_trigger = True
            if protein_inconsistency:
                trigger_counts["protein_inconsistency"] += 1
                should_trigger = True
            if orientation_ambiguous:
                trigger_counts["orientation_ambiguous"] += 1
                should_trigger = True

            if should_trigger:
                candidate_pairs.append((mapped_id, orig_id))
            elif uncertain_projection and not is_invalid:
                # Keep targeted mode conservative and avoid spending rescue on
                # already-valid genes that only have soft uncertainty markers.
                valid_uncertain_skipped += 1

        if not candidate_pairs:
            return 0
        if len(candidate_pairs) > self.config.internal_realign_max_genes:
            candidate_pairs = candidate_pairs[: self.config.internal_realign_max_genes]

        logger.info(
            "Internal realignment rescue: attempting targeted remap for %s gene(s)",
            len(candidate_pairs),
        )
        if trigger_counts:
            logger.info(
                "Internal realignment rescue triggers: invalid=%s invalid_uncertain=%s indel_jump=%s protein_inconsistency=%s orientation_ambiguous=%s valid_uncertain_skipped=%s",
                trigger_counts.get("invalid", 0),
                trigger_counts.get("invalid_uncertain", 0),
                trigger_counts.get("indel_jump", 0),
                trigger_counts.get("protein_inconsistency", 0),
                trigger_counts.get("orientation_ambiguous", 0),
                valid_uncertain_skipped,
            )
        total_candidates = len(candidate_pairs)
        progress_started = time.time()
        if total_candidates >= 2000:
            progress_interval = 100
        elif total_candidates >= 500:
            progress_interval = 50
        else:
            progress_interval = 25

        from stages.mapping import FeatureMapper

        validator = StructuralValidator(target_fasta)
        rescued = 0
        backend_counts = {"minimap2_cs": 0, "edlib": 0, "pairwise": 0}
        unresolved = 0
        not_improved = 0
        remap_failed = 0
        filtered_path = 0
        skipped_inputs = 0
        align_seconds = 0.0
        remap_seconds = 0.0
        validate_seconds = 0.0
        with FastaHandler(self.config.ref_fasta) as ref_fasta:
            rescue_mapper = FeatureMapper(
                self.syntenic_map,
                sex_chrom_map=self.sex_chrom_map,
                cnv_max_total_copies=self.config.cnv_max_total_copies,
                cnv_min_expansion_coverage=self.config.cnv_min_expansion_coverage,
                cnv_min_score_ratio=self.config.cnv_min_score_ratio,
                cnv_max_locus_overlap=self.config.cnv_max_locus_overlap,
                ref_fasta_handler=ref_fasta,
                target_fasta_handler=target_fasta,
                enable_boundary_refinement=False,
                boundary_refine_min_coverage=self.config.boundary_refine_min_coverage,
                boundary_refine_window_bp=self.config.boundary_refine_window_bp,
                boundary_refine_anchor_bp=self.config.boundary_refine_anchor_bp,
                boundary_refine_max_unaligned_bp=self.config.boundary_refine_max_unaligned_bp,
            )

            for idx, (mapped_id, orig_id) in enumerate(candidate_pairs, start=1):
                old_gene = self.mapped_genes.get(mapped_id)
                reference_gene = self.original_genes.get(orig_id)
                if old_gene is None or reference_gene is None:
                    skipped_inputs += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue

                ref_len = ref_fasta.get_length(reference_gene.seq_region)
                target_len = target_fasta.get_length(old_gene.seq_region)
                if ref_len <= 0 or target_len <= 0:
                    skipped_inputs += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue

                ref_window_start = max(
                    1,
                    reference_gene.start - self.config.internal_realign_ref_flank_bp,
                )
                ref_window_end = min(
                    ref_len,
                    reference_gene.end + self.config.internal_realign_ref_flank_bp,
                )
                target_window_start = max(
                    1,
                    old_gene.start - self.config.internal_realign_target_flank_bp,
                )
                target_window_end = min(
                    target_len,
                    old_gene.end + self.config.internal_realign_target_flank_bp,
                )
                if ref_window_start >= ref_window_end:
                    skipped_inputs += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue
                if target_window_start >= target_window_end:
                    skipped_inputs += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue

                ref_seq = ref_fasta.fetch(
                    reference_gene.seq_region,
                    ref_window_start,
                    ref_window_end,
                )
                target_seq = target_fasta.fetch(
                    old_gene.seq_region,
                    target_window_start,
                    target_window_end,
                )
                if not ref_seq or not target_seq:
                    skipped_inputs += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue

                # Per-gene windows are usually small; high thread counts can add
                # overhead for repeated short minimap2 invocations.
                local_window_span = max(
                    ref_window_end - ref_window_start + 1,
                    target_window_end - target_window_start + 1,
                )
                if self.config.threads <= 2:
                    local_threads = max(1, self.config.threads)
                elif local_window_span < 500_000:
                    local_threads = 1
                elif local_window_span < 2_000_000:
                    local_threads = 2
                else:
                    local_threads = min(4, self.config.threads)

                t_align = time.time()
                realigned = targeted_internal_realign(
                    ref_seq=ref_seq,
                    target_seq=target_seq,
                    ref_chr=reference_gene.seq_region,
                    target_chr=old_gene.seq_region,
                    gene_start=reference_gene.start,
                    gene_end=reference_gene.end,
                    ref_window_start=ref_window_start,
                    ref_window_end=ref_window_end,
                    target_window_start=target_window_start,
                    target_window_end=target_window_end,
                    threads=local_threads,
                    min_path_coverage=self.config.internal_realign_min_path_coverage,
                    max_internal_gap_bp=self.config.internal_realign_max_internal_gap_bp,
                    fallback_backend=self.config.internal_realign_fallback_backend,
                    temp_dir=self.config.temp_dir,
                )
                align_seconds += time.time() - t_align
                if realigned is None or not realigned.blocks:
                    unresolved += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue
                if realigned.backend in backend_counts:
                    backend_counts[realigned.backend] += 1
                if (
                    realigned.path_coverage < self.config.internal_realign_min_path_coverage
                    or realigned.internal_gap_bp > self.config.internal_realign_max_internal_gap_bp
                ):
                    filtered_path += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue

                t_remap = time.time()
                remap = rescue_mapper._map_gene_single_copy(
                    reference_gene,
                    mapped_gene_id=mapped_id,
                    candidate_blocks=realigned.blocks,
                    transcript_id_suffix="",
                )
                remap_seconds += time.time() - t_remap
                if remap.mapped is None:
                    remap_failed += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue

                candidate_gene = remap.mapped
                candidate_gene.feature_id = mapped_id
                candidate_gene.attributes["rescue_internal_realign"] = "true"
                candidate_gene.attributes["internal_realign_backend"] = realigned.backend
                candidate_gene.attributes["internal_realign_path_coverage"] = (
                    f"{realigned.path_coverage:.4f}"
                )
                candidate_gene.attributes["internal_realign_internal_gap_bp"] = str(
                    realigned.internal_gap_bp
                )
                for transcript in candidate_gene.transcripts:
                    transcript.attributes["rescue_internal_realign"] = "true"
                    transcript.attributes["internal_realign_backend"] = realigned.backend
                    transcript.attributes["internal_realign_path_coverage"] = (
                        f"{realigned.path_coverage:.4f}"
                    )
                    transcript.attributes["internal_realign_internal_gap_bp"] = str(
                        realigned.internal_gap_bp
                    )

                old_val = self.validation_results.get(mapped_id)
                t_validate = time.time()
                new_val = validator.validate_gene(candidate_gene)
                validate_seconds += time.time() - t_validate
                improved = True
                if old_val is not None:
                    improved = self._validation_score(new_val) > self._validation_score(old_val)
                    if not improved:
                        old_jump = self._has_exon_offset_discontinuity(
                            reference_gene,
                            old_gene,
                            self.config.internal_realign_trigger_indel_jump_bp,
                        )
                        new_jump = self._has_exon_offset_discontinuity(
                            reference_gene,
                            candidate_gene,
                            self.config.internal_realign_trigger_indel_jump_bp,
                        )
                        if old_jump and not new_jump:
                            improved = True
                        else:
                            old_orient = old_gene.attributes.get(
                                "orientation_confidence", "ambiguous"
                            )
                            new_orient = candidate_gene.attributes.get(
                                "orientation_confidence", "ambiguous"
                            )
                            if old_orient == "ambiguous" and new_orient == "high":
                                improved = True

                if self.config.internal_realign_accept_only_if_improved and not improved:
                    not_improved += 1
                    if idx % progress_interval == 0 or idx == total_candidates:
                        elapsed = max(0.001, time.time() - progress_started)
                        rate = idx / elapsed
                        remaining = max(0, total_candidates - idx)
                        eta_seconds = int(remaining / rate) if rate > 0 else -1
                        logger.info(
                            "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                            idx,
                            total_candidates,
                            (idx / total_candidates) * 100.0,
                            rescued,
                            unresolved,
                            filtered_path,
                            remap_failed,
                            not_improved,
                            backend_counts["minimap2_cs"],
                            backend_counts["edlib"],
                            backend_counts["pairwise"],
                            rate,
                            eta_seconds,
                        )
                    continue

                self.mapped_genes[mapped_id] = candidate_gene
                self.validation_results[mapped_id] = new_val
                rescued += 1
                if idx % progress_interval == 0 or idx == total_candidates:
                    elapsed = max(0.001, time.time() - progress_started)
                    rate = idx / elapsed
                    remaining = max(0, total_candidates - idx)
                    eta_seconds = int(remaining / rate) if rate > 0 else -1
                    logger.info(
                        "Internal realignment rescue progress: %s/%s (%.1f%%) updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] rate=%.2f genes/s eta=%ss",
                        idx,
                        total_candidates,
                        (idx / total_candidates) * 100.0,
                        rescued,
                        unresolved,
                        filtered_path,
                        remap_failed,
                        not_improved,
                        backend_counts["minimap2_cs"],
                        backend_counts["edlib"],
                        backend_counts["pairwise"],
                        rate,
                        eta_seconds,
                    )

        elapsed_total = max(0.001, time.time() - progress_started)
        logger.info(
            "Internal realignment rescue complete: processed=%s updated=%s unresolved=%s filtered=%s remap_failed=%s not_improved=%s skipped_inputs=%s backends[minimap2_cs=%s, edlib=%s, pairwise=%s] elapsed=%.1fs",
            total_candidates,
            rescued,
            unresolved,
            filtered_path,
            remap_failed,
            not_improved,
            skipped_inputs,
            backend_counts["minimap2_cs"],
            backend_counts["edlib"],
            backend_counts["pairwise"],
            elapsed_total,
        )
        core_seconds = align_seconds + remap_seconds + validate_seconds
        if core_seconds > 0:
            logger.info(
                "Internal realignment rescue timing breakdown: align=%.1fs (%.1f%%) remap=%.1fs (%.1f%%) validate=%.1fs (%.1f%%)",
                align_seconds,
                (align_seconds / core_seconds) * 100.0,
                remap_seconds,
                (remap_seconds / core_seconds) * 100.0,
                validate_seconds,
                (validate_seconds / core_seconds) * 100.0,
            )

        return rescued

    @staticmethod
    def _is_canonical_splice_for_positions(
        fasta: FastaHandler,
        seq_region: str,
        exon1_end: int,
        exon2_start: int,
        strand: str
    ) -> bool:
        """Check canonical donor/acceptor for a candidate intron boundary."""
        if strand == "+":
            donor = fasta.fetch(seq_region, exon1_end + 1, exon1_end + 2, strand=Strand.PLUS)
            acceptor = fasta.fetch(seq_region, exon2_start - 2, exon2_start - 1, strand=Strand.PLUS)
        else:
            donor = fasta.fetch(seq_region, exon2_start - 2, exon2_start - 1, strand=Strand.MINUS)
            acceptor = fasta.fetch(seq_region, exon1_end + 1, exon1_end + 2, strand=Strand.MINUS)
        return donor in CANONICAL_DONOR_SITES and acceptor in CANONICAL_ACCEPTOR_SITES

    @classmethod
    def _propose_transcript_splice_shift(
        cls,
        transcript: Transcript,
        tx_validation: TranscriptValidation,
        fasta: FastaHandler,
        max_shift_bp: int
    ) -> bool:
        """Apply minimal exon-boundary shifts that restore canonical splice sites."""
        if len(transcript.exons) < 2:
            return False
        
        sorted_exons = sorted(transcript.exons, key=lambda e: e.start)
        original_bounds = {
            exon.feature_id: (exon.start, exon.end)
            for exon in transcript.exons
        }
        changed = False
        
        for splice in tx_validation.splice_site_results:
            if splice.is_valid:
                continue
            intron_idx = splice.intron_index
            if intron_idx < 0 or intron_idx >= len(sorted_exons) - 1:
                continue
            
            left_exon = sorted_exons[intron_idx]
            right_exon = sorted_exons[intron_idx + 1]
            best = None
            
            for delta_end in range(-max_shift_bp, max_shift_bp + 1):
                for delta_start in range(-max_shift_bp, max_shift_bp + 1):
                    if delta_end == 0 and delta_start == 0:
                        continue
                    new_end = left_exon.end + delta_end
                    new_start = right_exon.start + delta_start
                    
                    if new_end < left_exon.start:
                        continue
                    if new_start > right_exon.end:
                        continue
                    if new_end + 1 >= new_start:
                        continue
                    if not cls._is_canonical_splice_for_positions(
                        fasta,
                        left_exon.seq_region,
                        new_end,
                        new_start,
                        str(transcript.strand),
                    ):
                        continue
                    
                    penalty = (
                        abs(delta_end) + abs(delta_start),
                        abs(delta_end - delta_start),
                        abs(delta_end),
                        abs(delta_start),
                    )
                    if best is None or penalty < best[0]:
                        best = (penalty, new_end, new_start)
            
            if best is None:
                continue
            
            _, new_end, new_start = best
            left_exon.interval.end = new_end
            right_exon.interval.start = new_start
            changed = True
        
        if not changed:
            return False
        
        # Update transcript interval from shifted exons.
        transcript.interval.start = min(e.start for e in transcript.exons)
        transcript.interval.end = max(e.end for e in transcript.exons)
        
        # Propagate shifts to CDS/UTR boundaries when they share exon boundaries.
        for exon in transcript.exons:
            old_start, old_end = original_bounds.get(exon.feature_id, (exon.start, exon.end))
            if exon.start == old_start and exon.end == old_end:
                continue
            for child in list(transcript.cds_list) + list(transcript.utrs):
                child_old_start, child_old_end = child.start, child.end
                if child.start == old_start:
                    child.interval.start = exon.start
                if child.end == old_end:
                    child.interval.end = exon.end
                # Guard against accidental inversion of child intervals.
                if child.start > child.end:
                    child.interval.start = child_old_start
                    child.interval.end = child_old_end
        
        return True

    def _run_splice_shift_rescue(self, target_fasta: FastaHandler) -> int:
        """Rescue non-canonical splice sites via small exon-boundary micro-shifts."""
        if self.config.splice_shift_max_bp <= 0:
            return 0
        
        validator = StructuralValidator(target_fasta)
        rescued = 0
        
        for mapped_id, mapped_gene in list(self.mapped_genes.items()):
            if mapped_gene.attributes.get("cnv_copy_role") == "expansion":
                continue
            old_val = self.validation_results.get(mapped_id)
            if old_val is None or old_val.is_valid:
                continue
            if not self._gene_has_uncertain_projection(mapped_gene):
                continue
            
            tx_val_by_id = {
                t.transcript_id: t for t in old_val.transcript_results
            }
            candidate_gene = copy.deepcopy(mapped_gene)
            changed = False
            for transcript in candidate_gene.transcripts:
                tx_val = tx_val_by_id.get(transcript.feature_id)
                if tx_val is None:
                    continue
                if tx_val.splice_sites_checked == tx_val.splice_sites_valid:
                    continue
                if self._propose_transcript_splice_shift(
                    transcript,
                    tx_val,
                    target_fasta,
                    max_shift_bp=self.config.splice_shift_max_bp,
                ):
                    changed = True
            
            if not changed:
                continue
            
            candidate_gene.interval.start = min(t.start for t in candidate_gene.transcripts)
            candidate_gene.interval.end = max(t.end for t in candidate_gene.transcripts)
            new_val = validator.validate_gene(candidate_gene)
            
            if self._validation_score(new_val) > self._validation_score(old_val):
                candidate_gene.attributes["rescue_splice_shift"] = "true"
                self.mapped_genes[mapped_id] = candidate_gene
                self.validation_results[mapped_id] = new_val
                rescued += 1
        
        return rescued

    @staticmethod
    def _transcript_validation_score(tx_validation: TranscriptValidation) -> Tuple[int, int, int, float, int]:
        """Score one transcript validation result for rescue comparison."""
        splice_score = (
            tx_validation.splice_sites_valid / max(1, tx_validation.splice_sites_checked)
        )
        return (
            1 if (tx_validation.start_codon_result is None or tx_validation.start_codon_valid) else 0,
            1 if (tx_validation.stop_codon_result is None or tx_validation.stop_codon_valid) else 0,
            1 if tx_validation.cds_length_valid else 0,
            splice_score,
            -len(tx_validation.errors),
        )

    @staticmethod
    def _apply_cds_boundary_shift(
        transcript: Transcript,
        cds_index: int,
        boundary_kind: str,
        delta: int
    ) -> bool:
        """Shift one CDS boundary while keeping it in exon bounds."""
        if delta == 0:
            return True
        if cds_index < 0 or cds_index >= len(transcript.cds_list):
            return False

        cds = transcript.cds_list[cds_index]
        if boundary_kind == "start":
            old_value = cds.start
            new_value = old_value + delta
            if new_value > cds.end:
                return False
        else:
            old_value = cds.end
            new_value = old_value + delta
            if new_value < cds.start:
                return False

        # Clamp to containing exon when possible.
        exon_limit = None
        for exon in transcript.exons:
            if exon.start <= old_value <= exon.end:
                exon_limit = (exon.start, exon.end)
                break
        if exon_limit is None:
            exon_limit = (transcript.start, transcript.end)

        if new_value < exon_limit[0] or new_value > exon_limit[1]:
            return False

        if boundary_kind == "start":
            cds.interval.start = new_value
        else:
            cds.interval.end = new_value

        # If UTR boundaries were coupled to CDS edge, keep them synchronized.
        for utr in transcript.utrs:
            if utr.start == old_value:
                utr.interval.start = new_value
            if utr.end == old_value:
                utr.interval.end = new_value

        return True

    @classmethod
    def _propose_transcript_codon_frame_shift(
        cls,
        transcript: Transcript,
        tx_validation: TranscriptValidation,
        validator: StructuralValidator,
        max_shift_bp: int
    ) -> bool:
        """Apply minimal CDS edge shifts that improve start/stop/frame validity."""
        if max_shift_bp <= 0 or not transcript.cds_list:
            return False
        if (
            tx_validation.start_codon_valid
            and tx_validation.stop_codon_valid
            and tx_validation.cds_length_valid
        ):
            return False

        sorted_indices = sorted(
            range(len(transcript.cds_list)),
            key=lambda i: transcript.cds_list[i].start
        )
        if transcript.strand == Strand.PLUS:
            start_idx = sorted_indices[0]
            stop_idx = sorted_indices[-1]
            start_boundary = "start"
            stop_boundary = "end"
        else:
            start_idx = sorted_indices[-1]
            stop_idx = sorted_indices[0]
            start_boundary = "end"
            stop_boundary = "start"

        start_deltas = [0] if tx_validation.start_codon_valid else list(
            range(-max_shift_bp, max_shift_bp + 1)
        )
        stop_deltas = [0] if tx_validation.stop_codon_valid else list(
            range(-max_shift_bp, max_shift_bp + 1)
        )
        if tx_validation.cds_length_valid:
            # Prefer codon repairs first when frame is already valid.
            pass

        baseline_score = cls._transcript_validation_score(tx_validation)
        best_score = baseline_score
        best_penalty = None
        best_transcript = None

        for delta_start in start_deltas:
            for delta_stop in stop_deltas:
                if delta_start == 0 and delta_stop == 0:
                    continue

                candidate = copy.deepcopy(transcript)
                if not cls._apply_cds_boundary_shift(
                    candidate, start_idx, start_boundary, delta_start
                ):
                    continue
                if not cls._apply_cds_boundary_shift(
                    candidate, stop_idx, stop_boundary, delta_stop
                ):
                    continue

                candidate.interval.start = min(e.start for e in candidate.exons)
                candidate.interval.end = max(e.end for e in candidate.exons)
                candidate_validation = validator.validate_transcript(candidate)
                candidate_score = cls._transcript_validation_score(candidate_validation)
                if candidate_score < best_score:
                    continue

                penalty = (
                    abs(delta_start) + abs(delta_stop),
                    abs(delta_start),
                    abs(delta_stop),
                )
                if candidate_score > best_score or best_penalty is None or penalty < best_penalty:
                    best_score = candidate_score
                    best_penalty = penalty
                    best_transcript = candidate

        if best_transcript is None or best_score <= baseline_score:
            return False

        transcript.interval = best_transcript.interval
        transcript.exons = best_transcript.exons
        transcript.cds_list = best_transcript.cds_list
        transcript.utrs = best_transcript.utrs
        return True

    def _run_codon_frame_rescue(self, target_fasta: FastaHandler) -> int:
        """Rescue start/stop/frame failures with micro-shifts at CDS boundaries."""
        if self.config.codon_rescue_max_bp <= 0:
            return 0

        validator = StructuralValidator(target_fasta)
        rescued = 0

        for mapped_id, mapped_gene in list(self.mapped_genes.items()):
            if mapped_gene.attributes.get("cnv_copy_role") == "expansion":
                continue
            old_val = self.validation_results.get(mapped_id)
            if old_val is None or old_val.is_valid:
                continue
            if not self._gene_has_uncertain_projection(mapped_gene):
                continue

            tx_val_by_id = {
                t.transcript_id: t for t in old_val.transcript_results
            }
            candidate_gene = copy.deepcopy(mapped_gene)
            changed = False
            for transcript in candidate_gene.transcripts:
                tx_val = tx_val_by_id.get(transcript.feature_id)
                if tx_val is None or not transcript.is_coding:
                    continue
                if self._propose_transcript_codon_frame_shift(
                    transcript,
                    tx_val,
                    validator,
                    max_shift_bp=self.config.codon_rescue_max_bp,
                ):
                    changed = True

            if not changed:
                continue

            candidate_gene.interval.start = min(t.start for t in candidate_gene.transcripts)
            candidate_gene.interval.end = max(t.end for t in candidate_gene.transcripts)
            new_val = validator.validate_gene(candidate_gene)

            if self._validation_score(new_val) > self._validation_score(old_val):
                candidate_gene.attributes["rescue_codon_frame"] = "true"
                self.mapped_genes[mapped_id] = candidate_gene
                self.validation_results[mapped_id] = new_val
                rescued += 1

        return rescued

    @staticmethod
    def _validation_snapshot(gene_validation: Optional[GeneValidation]) -> Dict:
        """Compact validation snapshot for audit records."""
        if gene_validation is None:
            return {}
        splice_checked = sum(
            t.splice_sites_checked for t in gene_validation.transcript_results
        )
        splice_valid = sum(
            t.splice_sites_valid for t in gene_validation.transcript_results
        )
        return {
            "is_valid": gene_validation.is_valid,
            "transcripts_checked": gene_validation.transcripts_checked,
            "transcripts_valid": gene_validation.transcripts_valid,
            "splice_checked": splice_checked,
            "splice_valid": splice_valid,
        }

    @staticmethod
    def _projection_snapshot(gene: Gene) -> Dict:
        """Projection uncertainty indicators derived from transcript attributes."""
        coverages = []
        clamped_5prime = 0
        clamped_3prime = 0
        boundary_refined = 0
        for tx in gene.transcripts:
            attrs = tx.attributes or {}
            cov = attrs.get("projection_coverage")
            if cov is not None:
                try:
                    coverages.append(float(cov))
                except (TypeError, ValueError):
                    pass
            if attrs.get("projection_clamped_5prime") == "true":
                clamped_5prime += 1
            if attrs.get("projection_clamped_3prime") == "true":
                clamped_3prime += 1
            if (
                attrs.get("projection_boundary_refined_5prime") == "true"
                or attrs.get("projection_boundary_refined_3prime") == "true"
            ):
                boundary_refined += 1
        return {
            "min_coverage": min(coverages) if coverages else None,
            "max_coverage": max(coverages) if coverages else None,
            "clamped_5prime_transcripts": clamped_5prime,
            "clamped_3prime_transcripts": clamped_3prime,
            "boundary_refined_transcripts": boundary_refined,
        }

    def _build_gene_audit_traces(self) -> Dict:
        """Per-gene audit trace including projection, rescue, and validation deltas."""
        traces = {}
        for result in self.mapping_results:
            ref_id = result.original.feature_id
            mapped_id = f"mapped_{ref_id}"
            pre_gene = self.pre_refinement_genes.get(mapped_id)
            final_gene = self.mapped_genes.get(mapped_id)
            gene_for_trace = final_gene or pre_gene

            actions = []
            if gene_for_trace is not None:
                attrs = gene_for_trace.attributes or {}
                if attrs.get("missed_locus_recovered") == "true":
                    actions.append("missed_locus_recovery")
                if attrs.get("paralog_reassigned") == "true":
                    actions.append("paralog_reassignment")
                if attrs.get("rescue_validation_remap") == "true":
                    actions.append("validation_remap_rescue")
                if attrs.get("rescue_internal_realign") == "true":
                    actions.append("internal_realign_rescue")
                if attrs.get("rescue_splice_shift") == "true":
                    actions.append("splice_shift_rescue")
                if attrs.get("rescue_codon_frame") == "true":
                    actions.append("codon_frame_rescue")

            initial_validation = self._validation_snapshot(
                self.initial_validation_results.get(mapped_id)
            )
            final_validation = self._validation_snapshot(
                self.validation_results.get(mapped_id)
            )
            protein_qc = self.protein_qc_results.get(mapped_id)

            trace = {
                "reference_gene_id": ref_id,
                "mapped_gene_id": mapped_id if gene_for_trace is not None else None,
                "in_output": final_gene is not None,
                "mapping_status_initial": result.status,
                "mapping_status_final": (
                    final_gene.mapping_status if final_gene is not None else "removed"
                ),
                "projection": (
                    self._projection_snapshot(gene_for_trace)
                    if gene_for_trace is not None else {}
                ),
                "rescue_actions": actions,
                "validation_initial": initial_validation,
                "validation_final": final_validation,
                "validation_improved": (
                    bool(initial_validation)
                    and bool(final_validation)
                    and (
                        final_validation.get("transcripts_valid", 0)
                        > initial_validation.get("transcripts_valid", 0)
                    )
                ),
                "protein_qc_status": protein_qc.status if protein_qc else None,
                "confidence_class": (
                    (gene_for_trace.attributes or {}).get("projection_confidence_class")
                    if gene_for_trace is not None else None
                ),
            }
            traces[ref_id] = trace
        return traces
    
    def _run_refinement(self):
        """Run Stage 3: Refinement."""
        logger.info("Stage 3: Refinement")
        start = time.time()
        
        # Save pre-refinement genes for synteny analysis
        self.pre_refinement_genes = dict(self.mapped_genes)
        
        # Get target chromosome lengths
        with FastaHandler(self.config.target_fasta) as fasta:
            target_chr_lengths = {
                seq: fasta.get_length(seq) for seq in fasta.sequences
            }
        
        self.mapped_genes, self.refinement_report = refine_mappings(
            self.original_genes,
            self.mapped_genes,
            target_chr_lengths,
            validation_results=self.validation_results,
        )
        if self.paralog_resolution_report:
            self.refinement_report["paralog_reassignment"] = self.paralog_resolution_report
        if self.protein_qc_report:
            self.refinement_report["protein_qc"] = self.protein_qc_report
        self.audit_traces = self._build_gene_audit_traces()
        self.refinement_report["gene_audit_traces"] = self.audit_traces
        
        self.timings["refinement"] = time.time() - start
        logger.info(f"After refinement: {len(self.mapped_genes)} genes (from {len(self.pre_refinement_genes)} mapped)")
    
    def _run_synteny_analysis(self):
        """Run synteny analysis (gene order, new overlaps).
        
        Uses pre-refinement mapped genes to give accurate picture of
        gene order conservation across all successfully mapped genes.
        """
        logger.info("Running synteny analysis...")
        # Use pre-refinement genes for synteny analysis - gives more complete picture
        genes_for_analysis = self.pre_refinement_genes if self.pre_refinement_genes else self.mapped_genes
        self.synteny_analysis = analyze_synteny(
            self.original_genes,
            genes_for_analysis
        )
    
    def _calculate_statistics(self):
        """Calculate mapping statistics."""
        self.stats = calculate_statistics(
            self.original_genes,
            self.mapping_results,
            self.validation_results
        )
        
        self.stats.synteny_time_seconds = self.timings.get("synteny", 0)
        self.stats.mapping_time_seconds = self.timings.get("mapping", 0)
        self.stats.validation_time_seconds = self.timings.get("validation", 0)
        
        # Set post-refinement output counts
        self.stats.genes_in_output = len(self.mapped_genes)
        self.stats.transcripts_in_output = sum(
            len(g.transcripts) for g in self.mapped_genes.values()
        )
        self.stats.genes_removed_by_refinement = len(self.pre_refinement_genes) - len(self.mapped_genes)
    
    def _generate_outputs(self):
        """Generate output files."""
        logger.info("Generating outputs...")
        
        # Write GFF3
        genes_list = list(self.mapped_genes.values())
        write_genes_to_gff3(
            genes_list,
            self.config.output_gff,
            include_provenance=True
        )
        logger.info(f"Wrote {len(genes_list)} genes to {self.config.output_gff}")
        
        # Write statistics
        if self.config.output_stats:
            generate_report(
                self.stats,
                self.refinement_report,
                self.config.output_stats
            )
            logger.info(f"Wrote statistics to {self.config.output_stats}")
            
            # Write removed genes report
            if self.refinement_report.get("conflict_filtered_genes"):
                removed_file = self.config.output_stats.with_suffix('.removed.tsv')
                with open(removed_file, 'w') as f:
                    f.write("gene_id\tgene_name\treason\tmapped_position\n")
                    for entry in self.refinement_report["conflict_filtered_genes"]:
                        f.write(f"{entry['gene_id']}\t{entry.get('gene_name', '')}\t{entry['reason']}\t{entry['mapped_position']}\n")
                logger.info(f"Wrote {len(self.refinement_report['conflict_filtered_genes'])} removed genes to {removed_file}")
            if self.audit_traces:
                audit_file = self.config.output_stats.with_suffix('.audit.json')
                with open(audit_file, 'w') as f:
                    json.dump(self.audit_traces, f, indent=2)
                logger.info(f"Wrote {len(self.audit_traces)} per-gene audit traces to {audit_file}")
        
        # Print summary
        print_summary(self.stats)
        
        # Print synteny analysis
        if self.synteny_analysis:
            print_synteny_summary(self.synteny_analysis)


def run_pipeline(
    ref_fasta: str,
    ref_gff: str,
    target_fasta: str,
    output_gff: str,
    output_stats: Optional[str] = None,
    chromosomes: Optional[Set[str]] = None,
    threads: int = 4,
    min_identity: float = 0.95,
    min_mapq: int = 10,
    primary_only: bool = True,
    cnv_max_total_copies: int = 3,
    cnv_min_expansion_coverage: float = 0.60,
    cnv_min_score_ratio: float = 0.85,
    cnv_max_locus_overlap: float = 0.30,
    cnv_group_min_reciprocal_overlap: float = 0.50,
    cnv_ambiguity_distance_bp: int = 10000,
    cnv_ambiguity_score_delta: float = 0.01,
    enable_paralog_reassignment: bool = True,
    paralog_overlap_threshold: float = 0.50,
    paralog_max_locus_candidates: int = 8,
    paralog_max_rounds: int = 2,
    enable_missed_locus_recovery: bool = True,
    missed_locus_max_genes: int = 1000,
    missed_locus_search_padding_bp: int = 200000,
    missed_locus_max_window_bp: int = 800000,
    missed_locus_min_identity: float = 0.85,
    missed_locus_min_coverage: float = 0.60,
    enable_boundary_refinement: bool = True,
    boundary_refine_min_coverage: float = 0.98,
    boundary_refine_window_bp: int = 500,
    boundary_refine_anchor_bp: int = 80,
    boundary_refine_max_unaligned_bp: int = 8,
    enable_validation_rescue: bool = True,
    validation_rescue_window_bp: int = 1200,
    validation_rescue_anchor_bp: int = 120,
    validation_rescue_max_genes: int = 2000,
    enable_internal_realign_rescue: bool = True,
    internal_realign_max_genes: int = 3000,
    internal_realign_ref_flank_bp: int = 3000,
    internal_realign_target_flank_bp: int = 10000,
    internal_realign_min_path_coverage: float = 0.98,
    internal_realign_max_internal_gap_bp: int = 6,
    internal_realign_trigger_indel_jump_bp: int = 6,
    internal_realign_fallback_backend: str = "edlib",
    internal_realign_accept_only_if_improved: bool = True,
    enable_splice_shift_rescue: bool = True,
    splice_shift_max_bp: int = 6,
    enable_codon_frame_rescue: bool = True,
    codon_rescue_max_bp: int = 6,
    enable_protein_qc: bool = True,
    protein_qc_min_identity: float = 0.90,
    protein_qc_min_coverage: float = 0.90,
) -> MappingStatistics:
    """Convenience function to run the complete pipeline.
    
    Args:
        ref_fasta: Path to reference FASTA
        ref_gff: Path to reference GFF3
        target_fasta: Path to target FASTA
        output_gff: Path for output GFF3
        output_stats: Optional path for statistics JSON
        chromosomes: Optional set of chromosomes to process
        threads: Number of threads
        min_identity: Minimum synteny identity
        min_mapq: Minimum synteny MAPQ
        primary_only: Keep only primary synteny alignments
        cnv_max_total_copies: Maximum total copies per mapped gene family candidate
        cnv_min_expansion_coverage: Minimum reference coverage for expansion loci
        cnv_min_score_ratio: Minimum score ratio vs primary locus for expansions
        cnv_max_locus_overlap: Maximum reciprocal overlap allowed between selected loci
        cnv_group_min_reciprocal_overlap: Minimum overlap to group expansion candidates
        cnv_ambiguity_distance_bp: Distance tolerance for ambiguous source assignment
        cnv_ambiguity_score_delta: Score tolerance for ambiguous source assignment
        enable_paralog_reassignment: Enable synteny-aware remapping of overlapping paralog clusters
        paralog_overlap_threshold: Reciprocal-overlap threshold for paralog conflict grouping
        paralog_max_locus_candidates: Maximum alternate loci evaluated per conflicting gene
        paralog_max_rounds: Maximum reassignment rounds for unresolved clusters
        enable_missed_locus_recovery: Enable second-pass local recovery for unmapped/partial genes
        missed_locus_max_genes: Maximum genes considered in missed-locus recovery
        missed_locus_search_padding_bp: Ref/target padding around nearby synteny for local search
        missed_locus_max_window_bp: Maximum target window size for local search
        missed_locus_min_identity: Minimum identity to accept missed-locus recovery
        missed_locus_min_coverage: Minimum coverage to accept missed-locus recovery
        enable_boundary_refinement: Enable local boundary realignment for uncertain edges
        boundary_refine_min_coverage: Coverage threshold for refinement trigger
        boundary_refine_window_bp: Window around projected edge for local realignment
        boundary_refine_anchor_bp: Reference anchor length for boundary realignment
        boundary_refine_max_unaligned_bp: Maximum unclipped boundary bases allowed
        enable_validation_rescue: Enable validation-guided remapping rescue
        validation_rescue_window_bp: Rescue remap window around uncertain edges
        validation_rescue_anchor_bp: Rescue remap anchor length
        validation_rescue_max_genes: Maximum number of genes to rescue-remap
        enable_internal_realign_rescue: Enable targeted internal realignment rescue
        internal_realign_max_genes: Maximum genes considered for internal realignment rescue
        internal_realign_ref_flank_bp: Reference flank around gene used in rescue windows
        internal_realign_target_flank_bp: Target flank around mapped locus used in rescue windows
        internal_realign_min_path_coverage: Minimum path coverage required for rescue path acceptance
        internal_realign_max_internal_gap_bp: Maximum uncovered internal gap for rescue path acceptance
        internal_realign_trigger_indel_jump_bp: Exon-offset jump threshold triggering rescue
        internal_realign_fallback_backend: Fallback backend used when minimap2 path is fragmented
        internal_realign_accept_only_if_improved: Keep rescue update only if validation score improves
        enable_splice_shift_rescue: Enable splice micro-shift rescue
        splice_shift_max_bp: Maximum per-boundary shift for splice rescue
        enable_codon_frame_rescue: Enable codon/frame rescue via CDS boundary micro-shifts
        codon_rescue_max_bp: Maximum CDS boundary shift considered in codon/frame rescue
        enable_protein_qc: Enable protein-level QC against reference coding transcripts
        protein_qc_min_identity: Identity threshold for classifying protein QC as OK
        protein_qc_min_coverage: Coverage threshold for classifying protein QC as OK
        
    Returns:
        MappingStatistics
    """
    config = MappingConfig(
        ref_fasta=Path(ref_fasta),
        ref_gff=Path(ref_gff),
        target_fasta=Path(target_fasta),
        output_gff=Path(output_gff),
        output_stats=Path(output_stats) if output_stats else None,
        chromosomes=chromosomes,
        threads=threads,
        min_identity=min_identity,
        min_mapq=min_mapq,
        primary_only=primary_only,
        cnv_max_total_copies=cnv_max_total_copies,
        cnv_min_expansion_coverage=cnv_min_expansion_coverage,
        cnv_min_score_ratio=cnv_min_score_ratio,
        cnv_max_locus_overlap=cnv_max_locus_overlap,
        cnv_group_min_reciprocal_overlap=cnv_group_min_reciprocal_overlap,
        cnv_ambiguity_distance_bp=cnv_ambiguity_distance_bp,
        cnv_ambiguity_score_delta=cnv_ambiguity_score_delta,
        enable_paralog_reassignment=enable_paralog_reassignment,
        paralog_overlap_threshold=paralog_overlap_threshold,
        paralog_max_locus_candidates=paralog_max_locus_candidates,
        paralog_max_rounds=paralog_max_rounds,
        enable_missed_locus_recovery=enable_missed_locus_recovery,
        missed_locus_max_genes=missed_locus_max_genes,
        missed_locus_search_padding_bp=missed_locus_search_padding_bp,
        missed_locus_max_window_bp=missed_locus_max_window_bp,
        missed_locus_min_identity=missed_locus_min_identity,
        missed_locus_min_coverage=missed_locus_min_coverage,
        enable_boundary_refinement=enable_boundary_refinement,
        boundary_refine_min_coverage=boundary_refine_min_coverage,
        boundary_refine_window_bp=boundary_refine_window_bp,
        boundary_refine_anchor_bp=boundary_refine_anchor_bp,
        boundary_refine_max_unaligned_bp=boundary_refine_max_unaligned_bp,
        enable_validation_rescue=enable_validation_rescue,
        validation_rescue_window_bp=validation_rescue_window_bp,
        validation_rescue_anchor_bp=validation_rescue_anchor_bp,
        validation_rescue_max_genes=validation_rescue_max_genes,
        enable_internal_realign_rescue=enable_internal_realign_rescue,
        internal_realign_max_genes=internal_realign_max_genes,
        internal_realign_ref_flank_bp=internal_realign_ref_flank_bp,
        internal_realign_target_flank_bp=internal_realign_target_flank_bp,
        internal_realign_min_path_coverage=internal_realign_min_path_coverage,
        internal_realign_max_internal_gap_bp=internal_realign_max_internal_gap_bp,
        internal_realign_trigger_indel_jump_bp=internal_realign_trigger_indel_jump_bp,
        internal_realign_fallback_backend=internal_realign_fallback_backend,
        internal_realign_accept_only_if_improved=internal_realign_accept_only_if_improved,
        enable_splice_shift_rescue=enable_splice_shift_rescue,
        splice_shift_max_bp=splice_shift_max_bp,
        enable_codon_frame_rescue=enable_codon_frame_rescue,
        codon_rescue_max_bp=codon_rescue_max_bp,
        enable_protein_qc=enable_protein_qc,
        protein_qc_min_identity=protein_qc_min_identity,
        protein_qc_min_coverage=protein_qc_min_coverage,
    )
    
    pipeline = PangenomeMappingPipeline(config)
    return pipeline.run()

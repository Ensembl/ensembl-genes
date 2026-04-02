"""Stage 2: Feature mapping from reference to target genome.

Maps genes, transcripts, exons, and CDS features using syntenic blocks
and coordinate projection.
"""

import copy
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from src.models import (
    CDS, Exon, Gene, GenomicFeature, GenomicInterval, MappingResult,
    Strand, SyntenicBlock, SyntenicMap, Transcript, UTR
)
from src.config import (
    CANONICAL_ACCEPTOR_SITES,
    CANONICAL_DONOR_SITES,
    MappingStatus,
)

from .coordinate_projection import (
    CoordinateProjector,
    project_feature_coordinates,
    select_synteny_path_blocks,
)

logger = logging.getLogger(__name__)


@dataclass
class FeatureMappingResult:
    """Result of mapping a single feature."""
    original: GenomicFeature
    mapped: Optional[GenomicFeature] = None
    status: str = MappingStatus.UNMAPPED
    source_blocks: List[SyntenicBlock] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)


@dataclass
class TranscriptMappingResult:
    """Result of mapping a transcript with all its children."""
    original: Transcript
    mapped: Optional[Transcript] = None
    status: str = MappingStatus.UNMAPPED
    
    # Detailed results for sub-features
    exon_results: List[FeatureMappingResult] = field(default_factory=list)
    cds_results: List[FeatureMappingResult] = field(default_factory=list)
    utr_results: List[FeatureMappingResult] = field(default_factory=list)
    
    # Validation
    all_exons_mapped: bool = False
    all_cds_mapped: bool = False
    exon_order_preserved: bool = False
    
    errors: List[str] = field(default_factory=list)


@dataclass
class GeneMappingResult:
    """Result of mapping a gene with all its transcripts."""
    original: Gene
    mapped: Optional[Gene] = None
    status: str = MappingStatus.UNMAPPED
    
    # Transcript results
    transcript_results: List[TranscriptMappingResult] = field(default_factory=list)
    
    # Additional projected copies (CNV expansions)
    additional_mapped: List[Gene] = field(default_factory=list)
    
    # Summary
    transcripts_mapped: int = 0
    transcripts_partial: int = 0
    transcripts_unmapped: int = 0
    
    errors: List[str] = field(default_factory=list)


@dataclass
class GeneLocusCandidate:
    """Candidate target locus for mapping a reference gene."""
    locus_id: str
    target_chr: str
    is_inverted: bool
    target_start: int
    target_end: int
    blocks: List[SyntenicBlock]
    ref_coverage: float
    weighted_identity: float
    score: float


@dataclass
class ExpansionCandidate:
    """Potential extra copy projected for a reference gene."""
    gene: Gene
    owner_result: GeneMappingResult


@dataclass
class OrientationCandidateEvidence:
    """Evidence summary used for deterministic orientation arbitration."""
    orientation_score: int
    splice_support: int
    neighborhood_score: int
    confidence: str
    source: str
    reconciled: bool


class FeatureMapper:
    """Maps features from reference to target using syntenic blocks."""
    
    def __init__(
        self,
        syntenic_map: SyntenicMap,
        sex_chrom_map=None,
        cnv_max_total_copies: int = 3,
        cnv_min_expansion_coverage: float = 0.60,
        cnv_min_score_ratio: float = 0.85,
        cnv_max_locus_overlap: float = 0.30,
        ref_fasta_handler=None,
        target_fasta_handler=None,
        enable_boundary_refinement: bool = True,
        boundary_refine_min_coverage: float = 0.98,
        boundary_refine_window_bp: int = 500,
        boundary_refine_anchor_bp: int = 80,
        boundary_refine_max_unaligned_bp: int = 8,
    ):
        """Initialize mapper with syntenic map.
        
        Args:
            syntenic_map: Pre-computed syntenic map between genomes
            sex_chrom_map: Optional SexChromosomeMap for X/Y constraint enforcement
        """
        self.syntenic_map = syntenic_map
        self.sex_chrom_map = sex_chrom_map
        self.cnv_max_total_copies = cnv_max_total_copies
        self.cnv_min_expansion_coverage = cnv_min_expansion_coverage
        self.cnv_min_score_ratio = cnv_min_score_ratio
        self.cnv_max_locus_overlap = cnv_max_locus_overlap
        self.ref_fasta_handler = ref_fasta_handler
        self.target_fasta_handler = target_fasta_handler
        self.enable_boundary_refinement = enable_boundary_refinement
        self.boundary_refine_min_coverage = boundary_refine_min_coverage
        self.boundary_refine_window_bp = boundary_refine_window_bp
        self.boundary_refine_anchor_bp = boundary_refine_anchor_bp
        self.boundary_refine_max_unaligned_bp = boundary_refine_max_unaligned_bp

    @staticmethod
    def _flip_strand(strand: Strand) -> Strand:
        """Return the opposite strand for stranded features."""
        if strand == Strand.PLUS:
            return Strand.MINUS
        if strand == Strand.MINUS:
            return Strand.PLUS
        return strand
    
    @staticmethod
    def _interval_union_length(intervals: List[Tuple[int, int]]) -> int:
        """Length of union for closed intervals."""
        if not intervals:
            return 0
        
        merged = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1] + 1:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        return sum(end - start + 1 for start, end in merged)
    
    def _collect_gene_locus_candidates(
        self,
        gene: Gene,
        max_target_gap: int = 100000
    ) -> List[GeneLocusCandidate]:
        """Collect and rank candidate loci for a reference gene.
        
        Loci are clustered by target chromosome/orientation and proximity.
        """
        candidate_blocks = self.syntenic_map.find_blocks_for_ref_region(
            gene.seq_region,
            gene.start,
            gene.end
        )
        
        # Apply sex chromosome constraints if configured
        if self.sex_chrom_map and candidate_blocks:
            candidate_blocks = [
                b for b in candidate_blocks
                if self.sex_chrom_map.should_allow_mapping(
                    gene.seq_region,
                    b.target_interval.seq_region,
                    gene.start,
                    gene.end
                )
            ]
        
        block_rows = []
        for block in candidate_blocks:
            overlap_start = max(gene.start, block.ref_interval.start)
            overlap_end = min(gene.end, block.ref_interval.end)
            if overlap_end < overlap_start:
                continue
            overlap_len = overlap_end - overlap_start + 1
            block_rows.append(
                (
                    block.target_interval.seq_region,
                    block.is_inverted,
                    block.target_interval.start,
                    block.target_interval.end,
                    overlap_start,
                    overlap_end,
                    overlap_len,
                    block,
                )
            )
        
        if not block_rows:
            return []
        
        # Sort blocks by target locus then position and cluster nearby intervals.
        block_rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
        
        clusters = []
        current = None
        for row in block_rows:
            target_chr, is_inverted, t_start, t_end, r_start, r_end, overlap_len, block = row
            if (
                current is None
                or current["target_chr"] != target_chr
                or current["is_inverted"] != is_inverted
                or t_start > current["target_end"] + max_target_gap
            ):
                current = {
                    "target_chr": target_chr,
                    "is_inverted": is_inverted,
                    "target_start": t_start,
                    "target_end": t_end,
                    "blocks": [],
                    "ref_intervals": [],
                    "identity_weighted": 0.0,
                    "overlap_bp": 0,
                }
                clusters.append(current)
            else:
                current["target_end"] = max(current["target_end"], t_end)
            
            current["blocks"].append(block)
            current["ref_intervals"].append((r_start, r_end))
            current["identity_weighted"] += overlap_len * block.identity
            current["overlap_bp"] += overlap_len
        
        candidates: List[GeneLocusCandidate] = []
        for cluster in clusters:
            ref_cov_bp = self._interval_union_length(cluster["ref_intervals"])
            ref_cov = ref_cov_bp / max(gene.length, 1)
            weighted_identity = (
                cluster["identity_weighted"] / max(cluster["overlap_bp"], 1)
            )
            # Synteny-first ranking with identity as secondary signal.
            score = (0.75 * ref_cov) + (0.25 * weighted_identity)
            
            locus_id = (
                f"{cluster['target_chr']}:{cluster['target_start']}-{cluster['target_end']}:"
                f"{'-' if cluster['is_inverted'] else '+'}"
            )
            candidates.append(
                GeneLocusCandidate(
                    locus_id=locus_id,
                    target_chr=cluster["target_chr"],
                    is_inverted=cluster["is_inverted"],
                    target_start=cluster["target_start"],
                    target_end=cluster["target_end"],
                    blocks=sorted(
                        cluster["blocks"],
                        key=lambda b: (b.ref_interval.start, b.ref_interval.end)
                    ),
                    ref_coverage=ref_cov,
                    weighted_identity=weighted_identity,
                    score=score,
                )
            )
        
        candidates.sort(
            key=lambda c: (c.score, c.ref_coverage, c.weighted_identity),
            reverse=True
        )
        return candidates
    
    def _map_single_feature(
        self,
        feature: GenomicFeature,
        new_id: Optional[str] = None,
        candidate_blocks: Optional[List[SyntenicBlock]] = None
    ) -> FeatureMappingResult:
        """Map a single feature's coordinates to target.
        
        Args:
            feature: The feature to map
            new_id: Optional new ID for mapped feature
            
        Returns:
            FeatureMappingResult
        """
        result = FeatureMappingResult(original=feature)
        
        # Find overlapping syntenic blocks
        if candidate_blocks is None:
            blocks = self.syntenic_map.find_blocks_for_ref_region(
                feature.seq_region,
                feature.start,
                feature.end
            )
        else:
            blocks = [
                b for b in candidate_blocks
                if b.ref_interval.seq_region == feature.seq_region
                and b.ref_interval.start <= feature.end
                and feature.start <= b.ref_interval.end
            ]
        
        # Filter blocks by sex chromosome constraints
        if self.sex_chrom_map and blocks:
            blocks = [
                b for b in blocks
                if self.sex_chrom_map.should_allow_mapping(
                    feature.seq_region,
                    b.target_interval.seq_region,
                    feature.start,
                    feature.end
                )
            ]
        
        if not blocks:
            result.errors.append(f"No syntenic blocks found for {feature.seq_region}:{feature.start}-{feature.end}")
            return result
        
        result.source_blocks = blocks
        
        # Project coordinates
        projection = project_feature_coordinates(
            feature.seq_region,
            feature.start,
            feature.end,
            blocks,
            ref_fasta=self.ref_fasta_handler,
            target_fasta=self.target_fasta_handler,
            enable_boundary_refinement=self.enable_boundary_refinement,
            boundary_refine_min_coverage=self.boundary_refine_min_coverage,
            boundary_refine_window_bp=self.boundary_refine_window_bp,
            boundary_refine_anchor_bp=self.boundary_refine_anchor_bp,
            boundary_refine_max_unaligned_bp=self.boundary_refine_max_unaligned_bp,
        )
        
        if projection is None:
            result.errors.append("Coordinate projection failed")
            return result
        
        target_chr, target_start, target_end, is_inverted, proj_info = projection
        
        # Track partial coverage and cross-chromosome mappings
        if proj_info.get('coverage', 1.0) < 1.0:
            result.errors.append(f"Partial coverage: {proj_info['coverage']*100:.1f}%")
        if proj_info.get('clamped_5prime'):
            result.errors.append("5' end clamped to block boundary")
        if proj_info.get('clamped_3prime'):
            result.errors.append("3' end clamped to block boundary")
        if proj_info.get('cross_chromosome'):
            result.errors.append(f"Cross-chromosome mapping: {feature.seq_region} -> {target_chr}")
        if proj_info.get('boundary_refined_5prime'):
            result.errors.append("5' boundary refined via local realignment")
        if proj_info.get('boundary_refined_3prime'):
            result.errors.append("3' boundary refined via local realignment")
        
        # Determine target strand
        if is_inverted:
            if feature.strand == Strand.PLUS:
                target_strand = Strand.MINUS
            elif feature.strand == Strand.MINUS:
                target_strand = Strand.PLUS
            else:
                target_strand = feature.strand
        else:
            target_strand = feature.strand
        
        # Create mapped feature (shallow copy with new coordinates)
        mapped = copy.copy(feature)
        mapped.interval = GenomicInterval(
            seq_region=target_chr,
            start=target_start,
            end=target_end,
            strand=target_strand
        )
        mapped.feature_id = new_id or f"mapped_{feature.feature_id}"
        mapped.mapped_from = feature.feature_id
        mapped.attributes["projection_coverage"] = f"{proj_info.get('coverage', 1.0):.4f}"
        if proj_info.get("clamped_5prime"):
            mapped.attributes["projection_clamped_5prime"] = "true"
        if proj_info.get("clamped_3prime"):
            mapped.attributes["projection_clamped_3prime"] = "true"
        if proj_info.get("boundary_refined_5prime"):
            mapped.attributes["projection_boundary_refined_5prime"] = "true"
        if proj_info.get("boundary_refined_3prime"):
            mapped.attributes["projection_boundary_refined_3prime"] = "true"
        if proj_info.get("boundary_refinement_triggered"):
            mapped.attributes["projection_boundary_refinement_triggered"] = "true"
        if proj_info.get("path_blocks") is not None:
            mapped.attributes["projection_path_blocks"] = str(proj_info["path_blocks"])

        # Set status based on coverage
        if proj_info.get('coverage', 1.0) < 0.5:
            mapped.mapping_status = MappingStatus.PARTIAL
            result.status = MappingStatus.PARTIAL
        else:
            mapped.mapping_status = MappingStatus.MAPPED
            result.status = MappingStatus.MAPPED
        
        # Calculate identity from best block
        best_block = max(blocks, key=lambda b: b.identity)
        mapped.mapping_identity = best_block.identity
        
        result.mapped = mapped
        
        return result
    
    def map_exon(
        self,
        exon: Exon,
        transcript_id: str,
        candidate_blocks: Optional[List[SyntenicBlock]] = None
    ) -> FeatureMappingResult:
        """Map an exon to target.
        
        Args:
            exon: The exon to map
            transcript_id: Parent transcript ID for new ID generation
            
        Returns:
            FeatureMappingResult
        """
        result = self._map_single_feature(
            exon,
            new_id=f"{transcript_id}:exon:{exon.exon_number or 'unknown'}",
            candidate_blocks=candidate_blocks
        )
        
        if result.mapped:
            # Preserve exon-specific attributes
            mapped_exon = Exon(
                feature_id=result.mapped.feature_id,
                feature_type="exon",
                interval=result.mapped.interval,
                source=result.mapped.source,
                score=result.mapped.score,
                phase=result.mapped.phase,
                attributes=dict(result.mapped.attributes),
                exon_number=exon.exon_number,
                mapped_from=exon.feature_id,
                mapping_status=result.status,
                mapping_identity=result.mapped.mapping_identity
            )
            result.mapped = mapped_exon
        
        return result
    
    def map_cds(
        self,
        cds: CDS,
        transcript_id: str,
        candidate_blocks: Optional[List[SyntenicBlock]] = None
    ) -> FeatureMappingResult:
        """Map a CDS to target.
        
        Args:
            cds: The CDS to map
            transcript_id: Parent transcript ID
            
        Returns:
            FeatureMappingResult
        """
        result = self._map_single_feature(
            cds,
            new_id=f"{transcript_id}:cds",
            candidate_blocks=candidate_blocks
        )
        
        if result.mapped:
            mapped_cds = CDS(
                feature_id=result.mapped.feature_id,
                feature_type="CDS",
                interval=result.mapped.interval,
                source=result.mapped.source,
                score=result.mapped.score,
                phase=cds.phase,  # Preserve phase
                attributes=dict(result.mapped.attributes),
                mapped_from=cds.feature_id,
                mapping_status=result.status,
                mapping_identity=result.mapped.mapping_identity
            )
            result.mapped = mapped_cds
        
        return result
    
    def map_utr(
        self,
        utr: UTR,
        transcript_id: str,
        candidate_blocks: Optional[List[SyntenicBlock]] = None
    ) -> FeatureMappingResult:
        """Map a UTR to target.
        
        Args:
            utr: The UTR to map
            transcript_id: Parent transcript ID
            
        Returns:
            FeatureMappingResult
        """
        result = self._map_single_feature(
            utr,
            new_id=f"{transcript_id}:utr",
            candidate_blocks=candidate_blocks
        )
        
        if result.mapped:
            mapped_utr = UTR(
                feature_id=result.mapped.feature_id,
                feature_type=utr.feature_type,  # Preserve 5'/3' UTR type
                interval=result.mapped.interval,
                source=result.mapped.source,
                score=result.mapped.score,
                phase=result.mapped.phase,
                attributes=dict(result.mapped.attributes),
                mapped_from=utr.feature_id,
                mapping_status=result.status,
                mapping_identity=result.mapped.mapping_identity
            )
            result.mapped = mapped_utr
        
        return result
    
    def map_transcript(
        self,
        transcript: Transcript,
        gene_id: str,
        candidate_blocks: Optional[List[SyntenicBlock]] = None,
        expected_target_strand: Optional[Strand] = None,
        id_suffix: str = ""
    ) -> TranscriptMappingResult:
        """Map a transcript with all its sub-features.
        
        Args:
            transcript: The transcript to map
            gene_id: Parent gene ID
            
        Returns:
            TranscriptMappingResult
        """
        result = TranscriptMappingResult(original=transcript)
        mapped_transcript_id = f"mapped_{transcript.feature_id}{id_suffix}"

        transcript_candidate_blocks = candidate_blocks
        if candidate_blocks:
            path_blocks = select_synteny_path_blocks(
                transcript.seq_region,
                transcript.start,
                transcript.end,
                candidate_blocks,
            )
            if path_blocks:
                transcript_candidate_blocks = path_blocks
                if len(path_blocks) < len(candidate_blocks):
                    result.errors.append(
                        f"Restricted transcript projection to synteny path blocks "
                        f"({len(path_blocks)}/{len(candidate_blocks)})"
                    )
        
        # Map exons
        mapped_exons = []
        for exon in transcript.exons:
            exon_result = self.map_exon(
                exon,
                mapped_transcript_id,
                candidate_blocks=transcript_candidate_blocks
            )
            result.exon_results.append(exon_result)
            if exon_result.mapped:
                mapped_exons.append(exon_result.mapped)
        
        # Map CDS
        mapped_cds = []
        for cds in transcript.cds_list:
            cds_result = self.map_cds(
                cds,
                mapped_transcript_id,
                candidate_blocks=transcript_candidate_blocks
            )
            result.cds_results.append(cds_result)
            if cds_result.mapped:
                mapped_cds.append(cds_result.mapped)
        
        # Map UTRs
        mapped_utrs = []
        for utr in transcript.utrs:
            utr_result = self.map_utr(
                utr,
                mapped_transcript_id,
                candidate_blocks=transcript_candidate_blocks
            )
            result.utr_results.append(utr_result)
            if utr_result.mapped:
                mapped_utrs.append(utr_result.mapped)

        # ---- Proactive exon-offset continuity correction ----
        # When an indel between syntenic blocks shifts downstream exons, catch
        # the artefact here by comparing the (mapped − original) delta between
        # consecutive exons.  If the delta jumps by ≥ jump_threshold, correct
        # every downstream feature by the offset change.
        #
        # IMPORTANT: Only apply when the exon was projected through a clamped
        # boundary or multi-block path — these are the cases where indels
        # between blocks are not accounted for by the cs-tag offset map.
        # When both exons are projected through the same block with a cs-tag,
        # the delta change is correct and should NOT be corrected.
        if len(mapped_exons) >= 2 and len(transcript.exons) >= 2:
            jump_threshold = 3  # bp

            orig_by_id = {e.feature_id: e for e in transcript.exons}
            ordered_orig = sorted(transcript.exons, key=lambda e: e.start)

            # Build list of (original_exon, mapped_exon) pairs in original order
            paired = []
            for oe in ordered_orig:
                me = next(
                    (m for m in mapped_exons if m.mapped_from == oe.feature_id),
                    None,
                )
                if me is not None:
                    paired.append((oe, me))

            if len(paired) >= 2:
                start_deltas = [me.start - oe.start for oe, me in paired]
                end_deltas = [me.end - oe.end for oe, me in paired]

                # Find the first discontinuity and compute correction
                correction = 0
                corrections_applied = 0
                for idx in range(1, len(start_deltas)):
                    _, me_cur = paired[idx]
                    me_attrs = me_cur.attributes or {}

                    # Only consider correction when the projection shows signs
                    # of boundary clamping or multi-block stitching — these are
                    # the cases where the cs-tag offset map could not account
                    # for indels between blocks.
                    has_clamping = (
                        me_attrs.get("projection_clamped_5prime") == "true"
                        or me_attrs.get("projection_clamped_3prime") == "true"
                    )
                    multi_block = False
                    try:
                        multi_block = int(me_attrs.get("projection_path_blocks", "1")) > 1
                    except (TypeError, ValueError):
                        pass
                    partial_coverage = False
                    try:
                        partial_coverage = float(me_attrs.get("projection_coverage", "1.0")) < 1.0
                    except (TypeError, ValueError):
                        pass

                    if not (has_clamping or multi_block or partial_coverage):
                        continue

                    jump_start = start_deltas[idx] - start_deltas[idx - 1]
                    jump_end = end_deltas[idx] - end_deltas[idx - 1]

                    # Use the dominant jump direction
                    jump = jump_start if abs(jump_start) >= abs(jump_end) else jump_end

                    if abs(jump) >= jump_threshold:
                        correction += jump

                    if correction != 0:
                        me_cur.interval = GenomicInterval(
                            seq_region=me_cur.seq_region,
                            start=me_cur.start - correction,
                            end=me_cur.end - correction,
                            strand=me_cur.strand,
                        )
                        me_cur.attributes["exon_offset_corrected"] = str(correction)
                        corrections_applied += 1

                # Apply same correction to CDS and UTR features that overlap
                # corrected exon ranges.
                if corrections_applied > 0:
                    result.errors.append(
                        f"Exon-offset continuity correction applied to "
                        f"{corrections_applied} exon(s)"
                    )
                    corrected_exon_ids = {
                        me.mapped_from
                        for _, me in paired
                        if me.attributes.get("exon_offset_corrected")
                    }

                    # Re-derive correction for each CDS/UTR by matching its
                    # original position to the nearest corrected exon.
                    for cds in mapped_cds:
                        cds_correction = self._find_subfeature_correction(
                            cds, paired, correction, orig_by_id
                        )
                        if cds_correction != 0:
                            cds.interval = GenomicInterval(
                                seq_region=cds.seq_region,
                                start=cds.start - cds_correction,
                                end=cds.end - cds_correction,
                                strand=cds.strand,
                            )
                            cds.attributes["exon_offset_corrected"] = str(cds_correction)

                    for utr in mapped_utrs:
                        utr_correction = self._find_subfeature_correction(
                            utr, paired, correction, orig_by_id
                        )
                        if utr_correction != 0:
                            utr.interval = GenomicInterval(
                                seq_region=utr.seq_region,
                                start=utr.start - utr_correction,
                                end=utr.end - utr_correction,
                                strand=utr.strand,
                            )
                            utr.attributes["exon_offset_corrected"] = str(utr_correction)



        # Enforce single-locus transcript projection.
        # If exons land on multiple loci, keep only the dominant locus and mark partial.
        dominant_locus = None
        if mapped_exons:
            locus_counts = Counter((e.seq_region, e.strand) for e in mapped_exons)
            dominant_locus, dominant_count = max(locus_counts.items(), key=lambda item: item[1])
            if len(locus_counts) > 1:
                removed = len(mapped_exons) - dominant_count
                result.errors.append(
                    f"Exons mapped to multiple loci ({len(locus_counts)}); "
                    f"keeping dominant locus {dominant_locus[0]} {dominant_locus[1]} ({dominant_count} exons), "
                    f"dropping {removed}"
                )
                mapped_exons = [
                    e for e in mapped_exons
                    if (e.seq_region, e.strand) == dominant_locus
                ]
                mapped_cds = [
                    c for c in mapped_cds
                    if (c.seq_region, c.strand) == dominant_locus
                ]
                mapped_utrs = [
                    u for u in mapped_utrs
                    if (u.seq_region, u.strand) == dominant_locus
                ]
        
        # Check mapping completeness
        result.all_exons_mapped = len(mapped_exons) == len(transcript.exons)
        result.all_cds_mapped = len(mapped_cds) == len(transcript.cds_list)
        
        # If no exons mapped, transcript unmapped
        if not mapped_exons:
            result.status = MappingStatus.UNMAPPED
            result.errors.append("No exons could be mapped")
            return result
        
        # Check exon order preservation
        result.exon_order_preserved = self._check_exon_order(
            transcript.exons, mapped_exons
        )
        
        if not result.exon_order_preserved:
            result.errors.append("Exon order not preserved in mapping")

        (
            target_strand,
            strand_corrected,
            strand_source,
            strand_confidence,
        ) = self._reconcile_transcript_strand(
            transcript,
            mapped_exons,
            mapped_cds,
            mapped_utrs,
            expected_target_strand=expected_target_strand,
        )
        if strand_corrected:
            result.errors.append(
                "Transcript strand corrected using exon-order/splice evidence"
            )
        
        # Determine transcript span from mapped exons
        target_chr = mapped_exons[0].seq_region
        target_start = min(e.start for e in mapped_exons)
        target_end = max(e.end for e in mapped_exons)
        
        # Create mapped transcript
        mapped_transcript = Transcript(
            feature_id=mapped_transcript_id,
            feature_type=transcript.feature_type,
            interval=GenomicInterval(
                seq_region=target_chr,
                start=target_start,
                end=target_end,
                strand=target_strand
            ),
            source=transcript.source,
            score=transcript.score,
            attributes=dict(transcript.attributes),
            gene_id=gene_id,
            transcript_name=transcript.transcript_name,
            biotype=transcript.biotype,
            exons=mapped_exons,
            cds_list=mapped_cds,
            utrs=mapped_utrs,
            mapped_from=transcript.feature_id,
            mapping_status=MappingStatus.PARTIAL
        )
        
        # Set mapping identity as average of exon identities
        identities = [e.mapping_identity for e in mapped_exons if e.mapping_identity]
        if identities:
            mapped_transcript.mapping_identity = sum(identities) / len(identities)
        mapped_transcript.attributes["orientation_evidence_source"] = strand_source
        mapped_transcript.attributes["orientation_confidence"] = strand_confidence
        mapped_transcript.attributes["orientation_reconciled"] = (
            "true" if strand_corrected else "false"
        )
        
        result.mapped = mapped_transcript
        
        # Determine status
        if result.all_exons_mapped and result.all_cds_mapped and result.exon_order_preserved:
            result.status = MappingStatus.MAPPED
        else:
            result.status = MappingStatus.PARTIAL
        mapped_transcript.mapping_status = result.status
        
        return result

    @staticmethod
    def _find_subfeature_correction(
        subfeature: GenomicFeature,
        paired_exons: list,
        final_correction: int,
        orig_by_id: dict,
    ) -> int:
        """Find the exon-offset correction applicable to a CDS or UTR.

        Matches the sub-feature's original reference position to the nearest
        corrected exon and returns that exon's correction value.
        """
        if not subfeature.mapped_from:
            return 0

        # Try to find the original position of this sub-feature
        # by checking which original exon it overlapped.
        sub_start = subfeature.start
        sub_end = subfeature.end

        best_correction = 0
        best_overlap = 0
        for orig_exon, mapped_exon in paired_exons:
            corr_str = mapped_exon.attributes.get("exon_offset_corrected")
            if corr_str is None:
                continue
            exon_corr = int(corr_str)
            # Check overlap between sub-feature and mapped exon
            # (before correction was applied, the exon was at start+corr, end+corr)
            ov_start = max(sub_start, mapped_exon.start)
            ov_end = min(sub_end, mapped_exon.end)
            if ov_end >= ov_start:
                overlap = ov_end - ov_start + 1
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_correction = exon_corr

        return best_correction

    def _check_exon_order(
        self,
        original_exons: List[Exon],
        mapped_exons: List[Exon]
    ) -> bool:
        """Check if exon order is preserved after mapping.
        
        Args:
            original_exons: Original exons
            mapped_exons: Mapped exons
            
        Returns:
            True if order preserved
        """
        if len(mapped_exons) <= 1:
            return True
        
        # Get original order
        orig_sorted = sorted(original_exons, key=lambda e: e.start)
        
        # For mapped, check if positions are monotonic
        mapped_sorted = sorted(mapped_exons, key=lambda e: e.start)
        
        # Map original exon IDs to their positions
        orig_positions = {e.feature_id: i for i, e in enumerate(orig_sorted)}
        
        # Check if mapped exons maintain relative order
        mapped_orig_positions = []
        for me in mapped_sorted:
            orig_id = me.mapped_from
            if orig_id in orig_positions:
                mapped_orig_positions.append(orig_positions[orig_id])
        
        # Check monotonicity
        if len(mapped_orig_positions) <= 1:
            return True
        
        is_increasing = all(
            mapped_orig_positions[i] < mapped_orig_positions[i+1]
            for i in range(len(mapped_orig_positions) - 1)
        )
        is_decreasing = all(
            mapped_orig_positions[i] > mapped_orig_positions[i+1]
            for i in range(len(mapped_orig_positions) - 1)
        )
        
        return is_increasing or is_decreasing

    @staticmethod
    def _infer_inversion_from_exon_order(
        original_exons: List[Exon],
        mapped_exons: List[Exon],
    ) -> Optional[bool]:
        """Infer inversion from exon-order monotonicity.

        Returns:
            True for inverted order, False for forward order, None if ambiguous.
        """
        if len(mapped_exons) <= 1:
            return None

        orig_sorted = sorted(original_exons, key=lambda e: e.start)
        orig_positions = {e.feature_id: i for i, e in enumerate(orig_sorted)}
        mapped_sorted = sorted(mapped_exons, key=lambda e: e.start)

        mapped_orig_positions = [
            orig_positions[me.mapped_from]
            for me in mapped_sorted
            if me.mapped_from in orig_positions
        ]
        if len(mapped_orig_positions) <= 1:
            return None

        increasing = 0
        decreasing = 0
        for i in range(len(mapped_orig_positions) - 1):
            cur = mapped_orig_positions[i]
            nxt = mapped_orig_positions[i + 1]
            if cur < nxt:
                increasing += 1
            elif cur > nxt:
                decreasing += 1

        if increasing == 0 and decreasing == 0:
            return None
        if increasing == decreasing:
            return None
        return decreasing > increasing

    def _count_canonical_splice_sites(
        self,
        mapped_exons: List[Exon],
        strand: Strand,
    ) -> Optional[Tuple[int, int]]:
        """Count canonical donor/acceptor pairs for one strand hypothesis."""
        if self.target_fasta_handler is None or len(mapped_exons) < 2:
            return None
        if strand not in {Strand.PLUS, Strand.MINUS}:
            return None

        sorted_exons = sorted(mapped_exons, key=lambda e: e.start)
        canonical = 0
        total = 0

        for i in range(len(sorted_exons) - 1):
            exon1 = sorted_exons[i]
            exon2 = sorted_exons[i + 1]
            try:
                if strand == Strand.PLUS:
                    donor = self.target_fasta_handler.fetch(
                        exon1.seq_region, exon1.end + 1, exon1.end + 2, Strand.PLUS
                    )
                    acceptor = self.target_fasta_handler.fetch(
                        exon2.seq_region, exon2.start - 2, exon2.start - 1, Strand.PLUS
                    )
                else:
                    donor = self.target_fasta_handler.fetch(
                        exon2.seq_region, exon2.start - 2, exon2.start - 1, Strand.MINUS
                    )
                    acceptor = self.target_fasta_handler.fetch(
                        exon1.seq_region, exon1.end + 1, exon1.end + 2, Strand.MINUS
                    )
            except Exception:
                continue

            if len(donor) != 2 or len(acceptor) != 2:
                continue

            total += 1
            if (
                donor in CANONICAL_DONOR_SITES
                and acceptor in CANONICAL_ACCEPTOR_SITES
            ):
                canonical += 1

        if total == 0:
            return None
        return canonical, total

    def _strand_from_splice_motifs(
        self,
        mapped_exons: List[Exon],
        current_strand: Strand,
    ) -> Optional[Strand]:
        """Pick strand with stronger canonical splice support if decisive."""
        if current_strand not in {Strand.PLUS, Strand.MINUS}:
            return None

        current_score = self._count_canonical_splice_sites(mapped_exons, current_strand)
        flipped = self._flip_strand(current_strand)
        flipped_score = self._count_canonical_splice_sites(mapped_exons, flipped)
        if current_score is None or flipped_score is None:
            return None

        current_ok, total = current_score
        flipped_ok, flipped_total = flipped_score
        if total != flipped_total or total == 0:
            return None

        # Require a clear improvement to avoid noisy flips.
        min_flipped_support = 1 if total == 1 else max(2, int(total * 0.5))
        if (
            flipped_ok > current_ok
            and flipped_ok >= min_flipped_support
            and current_ok <= int(total * 0.25)
        ):
            return flipped
        return None

    def _reconcile_transcript_strand(
        self,
        original_transcript: Transcript,
        mapped_exons: List[Exon],
        mapped_cds: List[CDS],
        mapped_utrs: List[UTR],
        expected_target_strand: Optional[Strand] = None,
    ) -> Tuple[Strand, bool, str, str]:
        """Correct transcript/subfeature strand using order and splice evidence."""
        target_strand = mapped_exons[0].strand
        desired_strand: Optional[Strand] = None
        evidence_sources: List[str] = []

        if expected_target_strand in {Strand.PLUS, Strand.MINUS}:
            desired_strand = expected_target_strand
            evidence_sources.append("neighborhood")
        elif original_transcript.strand != Strand.UNSTRANDED:
            inferred_inversion = self._infer_inversion_from_exon_order(
                original_transcript.exons,
                mapped_exons,
            )
            if inferred_inversion is not None:
                desired_strand = (
                    self._flip_strand(original_transcript.strand)
                    if inferred_inversion
                    else original_transcript.strand
                )
                evidence_sources.append("order")

        splice_strand = self._strand_from_splice_motifs(mapped_exons, target_strand)
        if splice_strand is not None:
            desired_strand = splice_strand
            evidence_sources.append("splice")

        if desired_strand is None or desired_strand == target_strand:
            source = "mixed" if len(set(evidence_sources)) > 1 else (
                evidence_sources[-1] if evidence_sources else "neighborhood"
            )
            confidence = "high" if "splice" in evidence_sources else "ambiguous"
            return target_strand, False, source, confidence

        for exon in mapped_exons:
            exon.interval.strand = desired_strand
        for cds in mapped_cds:
            cds.interval.strand = desired_strand
        for utr in mapped_utrs:
            utr.interval.strand = desired_strand
        source = "mixed" if len(set(evidence_sources)) > 1 else (
            evidence_sources[-1] if evidence_sources else "neighborhood"
        )
        confidence = "high" if "splice" in evidence_sources else "ambiguous"
        return desired_strand, True, source, confidence
    
    @staticmethod
    def _reciprocal_overlap(
        start1: int,
        end1: int,
        start2: int,
        end2: int
    ) -> float:
        """Reciprocal overlap between two closed intervals."""
        ov_start = max(start1, start2)
        ov_end = min(end1, end2)
        if ov_end < ov_start:
            return 0.0
        ov_len = ov_end - ov_start + 1
        len1 = max(1, end1 - start1 + 1)
        len2 = max(1, end2 - start2 + 1)
        return ov_len / min(len1, len2)

    @staticmethod
    def _gene_exon_interval_set(gene: Gene) -> Set[Tuple[int, int]]:
        """Unique exon intervals across all transcripts."""
        intervals: Set[Tuple[int, int]] = set()
        for transcript in gene.transcripts:
            for exon in transcript.exons:
                intervals.add((exon.start, exon.end))
        return intervals

    def _is_redundant_expansion_copy(
        self,
        candidate: Gene,
        existing: List[Gene],
        gene_overlap_threshold: float = 0.98,
        exon_overlap_threshold: float = 0.95,
    ) -> bool:
        """True when an expansion effectively stacks on an existing mapped copy."""
        candidate_exons = self._gene_exon_interval_set(candidate)
        for chosen in existing:
            if candidate.seq_region != chosen.seq_region:
                continue

            ro = self._reciprocal_overlap(
                candidate.start,
                candidate.end,
                chosen.start,
                chosen.end,
            )
            if ro < gene_overlap_threshold:
                continue

            chosen_exons = self._gene_exon_interval_set(chosen)
            if candidate_exons and chosen_exons:
                shared = len(candidate_exons & chosen_exons)
                min_count = min(len(candidate_exons), len(chosen_exons))
                if min_count > 0 and (shared / min_count) >= exon_overlap_threshold:
                    return True

            # Near-identical locus span without exon-level disagreement is redundant.
            return True
        return False

    def _expected_target_strand_for_locus(
        self,
        reference_strand: Strand,
        locus_is_inverted: bool,
    ) -> Strand:
        """Expected mapped strand for a locus orientation hypothesis."""
        if reference_strand == Strand.UNSTRANDED:
            return Strand.UNSTRANDED
        if locus_is_inverted:
            return self._flip_strand(reference_strand)
        return reference_strand

    def _score_primary_gene_candidate(
        self,
        reference_gene: Gene,
        locus: GeneLocusCandidate,
        candidate_result: GeneMappingResult,
    ) -> Tuple[Tuple[int, int, int, float, float, float], OrientationCandidateEvidence]:
        """Rank candidate primary loci with deterministic orientation-first evidence."""
        mapped_gene = candidate_result.mapped
        if mapped_gene is None:
            return (
                (-1, -10**9, -10**9, 0.0, 0.0, 0.0),
                OrientationCandidateEvidence(
                    orientation_score=-1,
                    splice_support=-10**9,
                    neighborhood_score=-10**9,
                    confidence="ambiguous",
                    source="neighborhood",
                    reconciled=False,
                ),
            )

        expected_strand = self._expected_target_strand_for_locus(
            reference_gene.strand,
            locus.is_inverted,
        )
        tx_total = max(1, len(mapped_gene.transcripts))
        orientation_matches = 0
        if expected_strand in {Strand.PLUS, Strand.MINUS}:
            orientation_matches = sum(
                1 for transcript in mapped_gene.transcripts
                if transcript.strand == expected_strand
            )
        orientation_score = int(round(1000 * (orientation_matches / tx_total)))

        canonical = 0
        total = 0
        for transcript in mapped_gene.transcripts:
            strand_for_scoring = transcript.strand
            if expected_strand in {Strand.PLUS, Strand.MINUS}:
                strand_for_scoring = expected_strand
            score = self._count_canonical_splice_sites(
                transcript.exons,
                strand_for_scoring,
            )
            if score is None:
                continue
            tx_canonical, tx_total = score
            canonical += tx_canonical
            total += tx_total

        mapped_tx = candidate_result.transcripts_mapped
        partial_tx = candidate_result.transcripts_partial
        total_transcripts = max(
            1,
            mapped_tx + partial_tx + candidate_result.transcripts_unmapped,
        )
        completeness = (mapped_tx + (0.5 * partial_tx)) / total_transcripts
        noncanonical = total - canonical
        identity = mapped_gene.mapping_identity or 0.0

        # Neighbor/synteny proxy: prefer candidates that preserve reference span and
        # transcript completeness at the same time.
        neighborhood_score = int(round(1000 * min(locus.ref_coverage, completeness)))
        evidence_sources: List[str] = ["neighborhood"]
        if total > 0 and canonical > 0 and noncanonical <= max(0, total // 3):
            evidence_sources.append("splice")
        if orientation_score >= 800:
            evidence_sources.append("order")
        evidence_source = (
            "mixed" if len(set(evidence_sources)) > 1 else evidence_sources[0]
        )
        confidence = (
            "high"
            if (
                orientation_score >= 800
                or (total > 0 and canonical >= max(1, int(total * 0.6)))
            )
            else "ambiguous"
        )

        score = (
            orientation_score,
            canonical,
            neighborhood_score,
            completeness,
            locus.score,
            identity,
        )
        evidence = OrientationCandidateEvidence(
            orientation_score=orientation_score,
            splice_support=canonical,
            neighborhood_score=neighborhood_score,
            confidence=confidence,
            source=evidence_source,
            reconciled=False,
        )
        return score, evidence
    
    def _select_loci_for_expansion(
        self,
        candidates: List[GeneLocusCandidate],
        max_total_copies: int = 3,
        min_expansion_coverage: float = 0.60,
        min_score_ratio: float = 0.85,
        max_locus_overlap: float = 0.30
    ) -> List[GeneLocusCandidate]:
        """Choose loci for primary + expansion copies with parsimonious rules."""
        if not candidates:
            return []
        
        primary = candidates[0]
        selected = [primary]
        if max_total_copies <= 1:
            return selected
        
        for cand in candidates[1:]:
            if len(selected) >= max_total_copies:
                break
            if cand.ref_coverage < min_expansion_coverage:
                continue
            if primary.score > 0 and cand.score < (primary.score * min_score_ratio):
                continue
            
            overlaps_existing = False
            for chosen in selected:
                if cand.target_chr != chosen.target_chr:
                    continue
                if cand.is_inverted != chosen.is_inverted:
                    continue
                ro = self._reciprocal_overlap(
                    cand.target_start,
                    cand.target_end,
                    chosen.target_start,
                    chosen.target_end
                )
                if ro > max_locus_overlap:
                    overlaps_existing = True
                    break
            if overlaps_existing:
                continue
            selected.append(cand)
        
        return selected
    
    def _map_gene_single_copy(
        self,
        gene: Gene,
        mapped_gene_id: str,
        candidate_blocks: Optional[List[SyntenicBlock]],
        transcript_id_suffix: str
    ) -> GeneMappingResult:
        """Map one gene copy constrained to a candidate locus."""
        result = GeneMappingResult(original=gene)

        expected_target_strands: Dict[str, Strand] = {}
        if candidate_blocks:
            orientation_set = {block.is_inverted for block in candidate_blocks}
            if len(orientation_set) == 1:
                locus_inverted = next(iter(orientation_set))
                for transcript in gene.transcripts:
                    expected = transcript.strand
                    if transcript.strand in {Strand.PLUS, Strand.MINUS} and locus_inverted:
                        expected = self._flip_strand(transcript.strand)
                    expected_target_strands[transcript.feature_id] = expected
        
        # Map all transcripts
        mapped_transcripts = []
        for transcript in gene.transcripts:
            tx_result = self.map_transcript(
                transcript,
                mapped_gene_id,
                candidate_blocks=candidate_blocks,
                expected_target_strand=expected_target_strands.get(transcript.feature_id),
                id_suffix=transcript_id_suffix
            )
            result.transcript_results.append(tx_result)
            
            if tx_result.status == MappingStatus.MAPPED:
                result.transcripts_mapped += 1
                if tx_result.mapped:
                    mapped_transcripts.append(tx_result.mapped)
            elif tx_result.status == MappingStatus.PARTIAL:
                result.transcripts_partial += 1
                if tx_result.mapped:
                    mapped_transcripts.append(tx_result.mapped)
            else:
                result.transcripts_unmapped += 1
        
        # If no transcripts mapped, gene unmapped
        if not mapped_transcripts:
            result.status = MappingStatus.UNMAPPED
            result.errors.append("No transcripts could be mapped")
            return result
        
        # Determine gene span from mapped transcripts
        target_chr = mapped_transcripts[0].seq_region
        target_start = min(t.start for t in mapped_transcripts)
        target_end = max(t.end for t in mapped_transcripts)
        target_strand = mapped_transcripts[0].strand
        
        # VALIDATION: Check if mapped size is coherent with reference size
        # Reject mappings where target span is unreasonably large (>10x reference)
        # This catches cases where exons scatter across the genome
        ref_size = gene.end - gene.start + 1
        target_size = target_end - target_start + 1
        max_expansion = 10  # Allow up to 10x size increase
        
        if target_size > ref_size * max_expansion:
            result.status = MappingStatus.UNMAPPED
            result.errors.append(
                f"Mapped size ({target_size:,}bp) exceeds {max_expansion}x reference "
                f"size ({ref_size:,}bp) - likely scattered alignment"
            )
            return result
        
        # Also check if transcripts are on different chromosomes
        all_chrs = set(t.seq_region for t in mapped_transcripts)
        if len(all_chrs) > 1:
            result.status = MappingStatus.UNMAPPED
            result.errors.append(
                f"Transcripts mapped to multiple chromosomes: {all_chrs}"
            )
            return result
        
        all_strands = set(t.strand for t in mapped_transcripts)
        if len(all_strands) > 1:
            result.status = MappingStatus.UNMAPPED
            result.errors.append(
                f"Transcripts mapped to multiple strands: {all_strands}"
            )
            return result
        
        # Create mapped gene
        mapped_gene = Gene(
            feature_id=mapped_gene_id,
            feature_type=gene.feature_type,
            interval=GenomicInterval(
                seq_region=target_chr,
                start=target_start,
                end=target_end,
                strand=target_strand
            ),
            source=gene.source,
            score=gene.score,
            attributes=dict(gene.attributes),
            gene_name=gene.gene_name,
            biotype=gene.biotype,
            description=gene.description,
            transcripts=mapped_transcripts,
            mapped_from=gene.feature_id,
            mapping_status=MappingStatus.PARTIAL
        )
        
        # Set mapping identity as average of transcript identities
        identities = [t.mapping_identity for t in mapped_transcripts if t.mapping_identity]
        if identities:
            mapped_gene.mapping_identity = sum(identities) / len(identities)
        
        result.mapped = mapped_gene
        
        # Determine overall status
        if result.transcripts_unmapped == 0 and result.transcripts_partial == 0:
            result.status = MappingStatus.MAPPED
        elif result.transcripts_mapped > 0 or result.transcripts_partial > 0:
            result.status = MappingStatus.PARTIAL
        else:
            result.status = MappingStatus.UNMAPPED
        mapped_gene.mapping_status = result.status
        
        return result
    
    def map_gene(self, gene: Gene) -> GeneMappingResult:
        """Map a gene with all its transcripts.
        
        Args:
            gene: The gene to map
            
        Returns:
            GeneMappingResult
        """
        locus_candidates = self._collect_gene_locus_candidates(gene)
        if not locus_candidates:
            result = GeneMappingResult(original=gene)
            result.status = MappingStatus.UNMAPPED
            result.errors.append(
                f"No syntenic locus candidates for {gene.seq_region}:{gene.start}-{gene.end}"
            )
            return result
        
        selected_loci = self._select_loci_for_expansion(
            locus_candidates,
            max_total_copies=self.cnv_max_total_copies,
            min_expansion_coverage=self.cnv_min_expansion_coverage,
            min_score_ratio=self.cnv_min_score_ratio,
            max_locus_overlap=self.cnv_max_locus_overlap
        )
        if not selected_loci:
            selected_loci = [locus_candidates[0]]
        
        # Evaluate candidate loci and choose primary using splice/structure evidence.
        primary_idx = None
        primary_result = None
        primary_locus = None
        primary_orientation_evidence: Optional[OrientationCandidateEvidence] = None
        candidate_primary_results: List[
            Tuple[
                Tuple[int, int, int, float, float, float],
                OrientationCandidateEvidence,
                int,
                GeneLocusCandidate,
                GeneMappingResult,
            ]
        ] = []
        for idx, locus in enumerate(selected_loci, start=1):
            copy_result = self._map_gene_single_copy(
                gene,
                mapped_gene_id=f"mapped_{gene.feature_id}",
                candidate_blocks=locus.blocks,
                transcript_id_suffix=""
            )
            if copy_result.mapped:
                score, evidence = self._score_primary_gene_candidate(
                    gene,
                    locus,
                    copy_result,
                )
                candidate_primary_results.append((score, evidence, idx, locus, copy_result))

        if candidate_primary_results:
            candidate_primary_results.sort(key=lambda item: item[0], reverse=True)
            (
                _,
                primary_orientation_evidence,
                primary_idx,
                primary_locus,
                primary_result,
            ) = candidate_primary_results[0]
            if primary_idx != 1:
                primary_result.errors.append(
                    f"Primary locus switched from rank 1 to rank {primary_idx} "
                    "using orientation/splice support"
                )
        
        if primary_result is None or primary_result.mapped is None:
            # Fall back to unconstrained single-copy mapping.
            fallback = self._map_gene_single_copy(
                gene,
                mapped_gene_id=f"mapped_{gene.feature_id}",
                candidate_blocks=None,
                transcript_id_suffix=""
            )
            if fallback.mapped:
                fallback.mapped.attributes["cnv_copy_index"] = "1"
                fallback.mapped.attributes["cnv_copy_role"] = "primary"
                fallback.mapped.attributes["cnv_assignment_strategy"] = "fallback_unconstrained"
                fallback.mapped.attributes["cnv_assignment_confidence"] = "ambiguous"
                if fallback.mapped.gene_name:
                    fallback.mapped.attributes["cnv_reference_gene_name"] = fallback.mapped.gene_name
                    fallback.mapped.attributes["cnv_label"] = fallback.mapped.gene_name
            return fallback
        
        assert primary_locus is not None
        primary_gene = primary_result.mapped
        primary_gene.attributes["cnv_copy_index"] = "1"
        primary_gene.attributes["cnv_copy_role"] = "primary"
        primary_gene.attributes["cnv_locus_id"] = primary_locus.locus_id
        primary_gene.attributes["cnv_locus_score"] = f"{primary_locus.score:.4f}"
        primary_gene.attributes["cnv_locus_coverage"] = f"{primary_locus.ref_coverage:.4f}"
        primary_gene.attributes["cnv_assignment_strategy"] = "synteny_first_parsimony"
        primary_gene.attributes["cnv_assignment_confidence"] = "high"
        if primary_orientation_evidence is not None:
            primary_orientation_evidence.reconciled = primary_idx != 1
            primary_gene.attributes["orientation_confidence"] = (
                primary_orientation_evidence.confidence
            )
            primary_gene.attributes["orientation_evidence_source"] = (
                primary_orientation_evidence.source
            )
            primary_gene.attributes["orientation_reconciled"] = (
                "true" if primary_orientation_evidence.reconciled else "false"
            )
            for transcript in primary_gene.transcripts:
                transcript.attributes["orientation_confidence"] = (
                    primary_orientation_evidence.confidence
                )
                transcript.attributes["orientation_evidence_source"] = (
                    primary_orientation_evidence.source
                )
                transcript.attributes["orientation_reconciled"] = (
                    "true" if primary_orientation_evidence.reconciled else "false"
                )
        if primary_gene.gene_name:
            primary_gene.attributes["cnv_reference_gene_name"] = primary_gene.gene_name
            primary_gene.attributes["cnv_label"] = primary_gene.gene_name
        
        # Add extra mapped copies for expansion loci.
        additional_mapped: List[Gene] = []
        copy_number = 2
        for idx, locus in enumerate(selected_loci, start=1):
            if idx == primary_idx:
                continue
            
            suffix = f"__copy{copy_number}"
            mapped_gene_id = f"mapped_{gene.feature_id}{suffix}"
            copy_result = self._map_gene_single_copy(
                gene,
                mapped_gene_id=mapped_gene_id,
                candidate_blocks=locus.blocks,
                transcript_id_suffix=suffix
            )
            if (
                not copy_result.mapped
                or copy_result.status != MappingStatus.MAPPED
            ):
                continue
            
            copied_gene = copy_result.mapped
            if self._is_redundant_expansion_copy(
                copied_gene,
                [primary_gene] + additional_mapped,
            ):
                primary_result.errors.append(
                    "Suppressed redundant expansion copy overlapping an existing mapped copy "
                    f"at {copied_gene.seq_region}:{copied_gene.start}-{copied_gene.end}"
                )
                continue

            copied_gene.attributes["cnv_copy_index"] = str(copy_number)
            copied_gene.attributes["cnv_copy_role"] = "expansion"
            copied_gene.attributes["cnv_locus_id"] = locus.locus_id
            copied_gene.attributes["cnv_locus_score"] = f"{locus.score:.4f}"
            copied_gene.attributes["cnv_locus_coverage"] = f"{locus.ref_coverage:.4f}"
            copied_gene.attributes["cnv_assignment_strategy"] = "synteny_first_parsimony"
            copied_gene.attributes["cnv_assignment_confidence"] = (
                "high" if locus.score >= (primary_locus.score * 0.95) else "ambiguous"
            )
            copied_gene.attributes["cnv_source_inference"] = (
                f"candidate_from:{gene.feature_id}"
            )
            if copied_gene.gene_name:
                copied_gene.attributes["cnv_reference_gene_name"] = copied_gene.gene_name
                cnv_label = f"{copied_gene.gene_name}_copy{copy_number}"
                copied_gene.attributes["cnv_label"] = cnv_label
                copied_gene.attributes["Name"] = cnv_label
                copied_gene.attributes["gene_name"] = cnv_label
                copied_gene.gene_name = cnv_label
            
            additional_mapped.append(copied_gene)
            copy_number += 1
        
        primary_result.additional_mapped = additional_mapped
        
        return primary_result


def _interval_distance(
    start1: int,
    end1: int,
    start2: int,
    end2: int
) -> int:
    """Distance between two closed intervals (0 if overlapping)."""
    if end1 < start2:
        return start2 - end1
    if end2 < start1:
        return start1 - end2
    return 0


def _gene_reciprocal_overlap(gene1: Gene, gene2: Gene) -> float:
    """Reciprocal overlap between two genes on the same chromosome."""
    if gene1.seq_region != gene2.seq_region:
        return 0.0
    ov_start = max(gene1.start, gene2.start)
    ov_end = min(gene1.end, gene2.end)
    if ov_end < ov_start:
        return 0.0
    ov_len = ov_end - ov_start + 1
    min_len = max(1, min(gene1.length, gene2.length))
    return ov_len / min_len


def _safe_attr_float(gene: Gene, key: str, default: float = 0.0) -> float:
    """Read a float attribute from gene metadata."""
    value = gene.attributes.get(key)
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _score_expansion_candidate(
    candidate: ExpansionCandidate,
    primary_by_source_ref: Dict[str, Gene],
    source_context_by_ref: Dict[str, Dict[str, Optional[Gene]]]
) -> Tuple[int, int, int, int, float, float, int]:
    """Score one expansion candidate for cross-gene locus arbitration.
    
    Returns:
        Tuple ordered for descending comparison:
        (
            same_chrom_as_source,
            flank_support,
            breakpoint_support,
            -distance_to_source,
            locus_score,
            mapping_identity,
            distance_to_source
        )
    """
    gene = candidate.gene
    source_ref_id = gene.mapped_from
    source_primary = primary_by_source_ref.get(source_ref_id) if source_ref_id else None
    
    same_chr = 0
    distance = 10**12
    if source_primary and source_primary.seq_region == gene.seq_region:
        same_chr = 1
        distance = _interval_distance(
            gene.start, gene.end,
            source_primary.start, source_primary.end
        )

    flank_support = 0
    breakpoint_support = 0
    if source_ref_id and source_ref_id in source_context_by_ref:
        ctx = source_context_by_ref[source_ref_id]
        left_gene = ctx.get("left")
        right_gene = ctx.get("right")
        if left_gene and left_gene.seq_region == gene.seq_region and left_gene.end <= gene.start:
            flank_support += 1
        if right_gene and right_gene.seq_region == gene.seq_region and gene.end <= right_gene.start:
            flank_support += 1
        if (
            left_gene
            and right_gene
            and left_gene.seq_region == gene.seq_region
            and right_gene.seq_region == gene.seq_region
            and left_gene.end <= gene.start <= right_gene.start
        ):
            breakpoint_support = 1
    
    locus_score = _safe_attr_float(gene, "cnv_locus_score", default=0.0)
    identity = gene.mapping_identity or 0.0
    
    return (
        same_chr,
        flank_support,
        breakpoint_support,
        -distance,
        locus_score,
        identity,
        distance,
    )


def _is_confident_expansion_choice(
    best_score: Tuple[int, int, int, int, float, float, int],
    second_score: Optional[Tuple[int, int, int, int, float, float, int]],
    distance_tolerance_bp: int = 10000,
    locus_score_tolerance: float = 0.01
) -> bool:
    """Decide whether chosen expansion source assignment is high confidence."""
    (
        best_same_chr,
        best_flank_support,
        best_breakpoint_support,
        _,
        best_locus,
        _,
        best_distance,
    ) = best_score
    
    if best_same_chr == 0:
        return False
    if second_score is None:
        return best_flank_support > 0 or best_breakpoint_support > 0
    
    (
        second_same_chr,
        second_flank_support,
        second_breakpoint_support,
        _,
        second_locus,
        _,
        second_distance,
    ) = second_score
    if second_same_chr == 0:
        return True

    if best_flank_support > second_flank_support:
        return True
    if best_breakpoint_support > second_breakpoint_support:
        return True
    
    # Ambiguous when competing sources have very similar local support.
    if (
        second_flank_support >= best_flank_support
        and second_breakpoint_support >= best_breakpoint_support
        and
        second_distance <= best_distance + distance_tolerance_bp
        and second_locus >= best_locus - locus_score_tolerance
    ):
        return False
    
    return True


def _group_overlapping_expansion_candidates(
    candidates: List[ExpansionCandidate],
    min_reciprocal_overlap: float = 0.50
) -> List[List[ExpansionCandidate]]:
    """Group expansion candidates that represent the same target locus."""
    groups: List[List[ExpansionCandidate]] = []
    
    for candidate in sorted(
        candidates,
        key=lambda c: (c.gene.seq_region, c.gene.start, c.gene.end)
    ):
        placed = False
        for group in groups:
            if candidate.gene.seq_region != group[0].gene.seq_region:
                continue
            if any(
                _gene_reciprocal_overlap(candidate.gene, other.gene) >= min_reciprocal_overlap
                for other in group
            ):
                group.append(candidate)
                placed = True
                break
        if not placed:
            groups.append([candidate])
    
    return groups


def _build_source_synteny_context(
    sorted_reference_genes: List[Gene],
    primary_by_source_ref: Dict[str, Gene]
) -> Dict[str, Dict[str, Optional[Gene]]]:
    """Nearest mapped primary neighbors for each reference gene (left/right)."""
    context: Dict[str, Dict[str, Optional[Gene]]] = {}
    ordered = sorted(
        sorted_reference_genes,
        key=lambda g: (g.seq_region, g.start, g.feature_id),
    )
    for idx, ref_gene in enumerate(ordered):
        left = None
        for j in range(idx - 1, -1, -1):
            left_candidate = primary_by_source_ref.get(ordered[j].feature_id)
            if left_candidate is not None:
                left = left_candidate
                break
        right = None
        for j in range(idx + 1, len(ordered)):
            right_candidate = primary_by_source_ref.get(ordered[j].feature_id)
            if right_candidate is not None:
                right = right_candidate
                break
        context[ref_gene.feature_id] = {"left": left, "right": right}
    return context


def _canonical_cnv_base_name(gene: Gene) -> str:
    """Base symbol used for deterministic CNV labels."""
    base = (
        gene.attributes.get("cnv_reference_gene_name")
        or gene.gene_name
        or gene.mapped_from
        or gene.feature_id
    )
    if isinstance(base, str):
        if base.startswith("Name="):
            base = base[len("Name="):]
        if base.startswith("gene:"):
            base = base[len("gene:"):]
    return base


def _assign_deterministic_cnv_labels(mapped_genes: Dict[str, Gene]):
    """Apply deterministic cluster-level CNV numbering and ambiguity tags."""
    genes_by_source: Dict[str, List[Gene]] = defaultdict(list)
    for gene in mapped_genes.values():
        source_ref = gene.mapped_from or gene.feature_id
        genes_by_source[source_ref].append(gene)

    # Stable per-source copy index (primary first, then target order).
    for source_ref, source_genes in genes_by_source.items():
        ordered_source = sorted(
            source_genes,
            key=lambda g: (
                0 if g.attributes.get("cnv_copy_role") == "primary" else 1,
                g.seq_region,
                g.start,
                g.end,
                g.feature_id,
            ),
        )
        for idx, gene in enumerate(ordered_source, start=1):
            gene.attributes["cnv_copy_index"] = str(idx)

    cluster_members: Dict[str, List[Gene]] = defaultdict(list)
    for gene in mapped_genes.values():
        sources: Set[str] = set()
        if gene.mapped_from:
            sources.add(gene.mapped_from)
        candidate_sources = gene.attributes.get("cnv_candidate_sources")
        if candidate_sources:
            for source_id in candidate_sources.split(","):
                source_id = source_id.strip()
                if source_id:
                    sources.add(source_id)
        if not sources:
            sources.add(gene.feature_id)
        cluster_id = ",".join(sorted(sources))
        gene.attributes["cnv_cluster_id"] = cluster_id
        if len(sources) > 1:
            gene.attributes["cnv_cluster_sources"] = cluster_id
        cluster_members[cluster_id].append(gene)

    for cluster_genes in cluster_members.values():
        ordered_cluster = sorted(
            cluster_genes,
            key=lambda g: (g.seq_region, g.start, g.end, g.feature_id),
        )
        primary = [g for g in ordered_cluster if g.attributes.get("cnv_copy_role") == "primary"]
        non_primary = [g for g in ordered_cluster if g.attributes.get("cnv_copy_role") != "primary"]
        ordered_cluster = primary + non_primary
        for cluster_idx, gene in enumerate(ordered_cluster, start=1):
            copy_idx = int(gene.attributes.get("cnv_copy_index", "1"))

            gene.attributes["cnv_cluster_copy_number"] = str(cluster_idx)
            gene.attributes["cnv_labeling_strategy"] = "deterministic_cluster_order"

            base_name = _canonical_cnv_base_name(gene)
            if copy_idx <= 1:
                label = base_name
            else:
                label = f"{base_name}_copy{copy_idx}"

            if gene.attributes.get("cnv_assignment_confidence") == "ambiguous":
                gene.attributes["cnv_ambiguity_tag"] = "source_ambiguous"
                if copy_idx > 1:
                    label = f"{label}_ambig"

            gene.attributes["cnv_label"] = label
            gene.attributes["Name"] = label
            gene.attributes["gene_name"] = label
            gene.gene_name = label


def map_genes(
    genes: Dict[str, Gene],
    syntenic_map: SyntenicMap,
    progress_callback=None,
    sex_chrom_map=None,
    cnv_max_total_copies: int = 3,
    cnv_min_expansion_coverage: float = 0.60,
    cnv_min_score_ratio: float = 0.85,
    cnv_max_locus_overlap: float = 0.30,
    cnv_group_min_reciprocal_overlap: float = 0.50,
    cnv_ambiguity_distance_bp: int = 10000,
    cnv_ambiguity_score_delta: float = 0.01,
    ref_fasta_handler=None,
    target_fasta_handler=None,
    enable_boundary_refinement: bool = True,
    boundary_refine_min_coverage: float = 0.98,
    boundary_refine_window_bp: int = 500,
    boundary_refine_anchor_bp: int = 80,
    boundary_refine_max_unaligned_bp: int = 8,
) -> Tuple[Dict[str, Gene], List[GeneMappingResult]]:
    """Map all genes to target genome with CNV-aware expansion reconciliation."""
    mapper = FeatureMapper(
        syntenic_map,
        sex_chrom_map=sex_chrom_map,
        cnv_max_total_copies=cnv_max_total_copies,
        cnv_min_expansion_coverage=cnv_min_expansion_coverage,
        cnv_min_score_ratio=cnv_min_score_ratio,
        cnv_max_locus_overlap=cnv_max_locus_overlap,
        ref_fasta_handler=ref_fasta_handler,
        target_fasta_handler=target_fasta_handler,
        enable_boundary_refinement=enable_boundary_refinement,
        boundary_refine_min_coverage=boundary_refine_min_coverage,
        boundary_refine_window_bp=boundary_refine_window_bp,
        boundary_refine_anchor_bp=boundary_refine_anchor_bp,
        boundary_refine_max_unaligned_bp=boundary_refine_max_unaligned_bp,
    )
    
    # Sort genes by reference position for order-preserving mapping
    sorted_genes = sorted(
        genes.values(),
        key=lambda g: (g.seq_region, g.start)
    )
    
    mapped_genes: Dict[str, Gene] = {}
    results: List[GeneMappingResult] = []
    pending_expansions: List[ExpansionCandidate] = []
    
    for i, gene in enumerate(sorted_genes):
        result = mapper.map_gene(gene)
        results.append(result)
        
        if result.mapped:
            mapped_genes[result.mapped.feature_id] = result.mapped
            for extra_gene in result.additional_mapped:
                pending_expansions.append(
                    ExpansionCandidate(gene=extra_gene, owner_result=result)
                )
        
        if progress_callback and (i + 1) % 500 == 0:
            progress_callback(i + 1, len(genes))
    
    # Reconcile expansion candidates so one extra locus maps to one best source gene.
    accepted_expansion_ids: Set[str] = set()
    if pending_expansions:
        primary_by_source_ref: Dict[str, Gene] = {}
        for primary_gene in mapped_genes.values():
            source_ref = primary_gene.mapped_from
            if source_ref:
                primary_by_source_ref[source_ref] = primary_gene
        source_context_by_ref = _build_source_synteny_context(
            sorted_genes,
            primary_by_source_ref,
        )
        
        accepted_expansions: List[Gene] = []
        for group in _group_overlapping_expansion_candidates(
            pending_expansions,
            min_reciprocal_overlap=cnv_group_min_reciprocal_overlap
        ):
            scored: List[
                Tuple[
                    Tuple[int, int, int, int, float, float, int],
                    ExpansionCandidate
                ]
            ] = [
                (
                    _score_expansion_candidate(
                        candidate,
                        primary_by_source_ref,
                        source_context_by_ref,
                    ),
                    candidate,
                )
                for candidate in group
            ]
            scored.sort(key=lambda item: item[0][:-1], reverse=True)
            
            best_score, best_candidate = scored[0]
            second_score = scored[1][0] if len(scored) > 1 else None
            confident = _is_confident_expansion_choice(
                best_score,
                second_score,
                distance_tolerance_bp=cnv_ambiguity_distance_bp,
                locus_score_tolerance=cnv_ambiguity_score_delta
            )
            
            chosen_gene = best_candidate.gene
            source_ref_id = chosen_gene.mapped_from or best_candidate.owner_result.original.feature_id
            chosen_gene.attributes["cnv_source_inference"] = (
                f"duplicated_from:{source_ref_id}"
                if confident else
                f"ambiguous_best_source:{source_ref_id}"
            )
            chosen_gene.attributes["cnv_assignment_confidence"] = (
                "high" if confident else "ambiguous"
            )
            if best_score[0] == 1 and best_score[6] < 10**12:
                chosen_gene.attributes["cnv_source_distance_bp"] = str(best_score[6])
            chosen_gene.attributes["cnv_flank_support"] = str(best_score[1])
            chosen_gene.attributes["cnv_breakpoint_support"] = str(best_score[2])
            if len(group) > 1:
                candidate_sources = sorted({
                    c.gene.mapped_from or c.owner_result.original.feature_id
                    for c in group
                })
                chosen_gene.attributes["cnv_candidate_sources"] = ",".join(candidate_sources)
            if not confident:
                chosen_gene.attributes["cnv_ambiguity_tag"] = "source_ambiguous"
            
            accepted_expansions.append(chosen_gene)
            accepted_expansion_ids.add(chosen_gene.feature_id)
            
            for _, candidate in scored[1:]:
                candidate.owner_result.errors.append(
                    f"Suppressed overlapping expansion candidate at "
                    f"{candidate.gene.seq_region}:{candidate.gene.start}-{candidate.gene.end}"
                )
            if not confident:
                best_candidate.owner_result.errors.append(
                    "Expansion source assignment is ambiguous; kept best-scoring candidate"
                )
        
        for expansion_gene in accepted_expansions:
            mapped_genes[expansion_gene.feature_id] = expansion_gene
    
    # Trim per-gene expansion lists to accepted copies after arbitration.
    for result in results:
        if not result.additional_mapped:
            continue
        kept = [
            gene for gene in result.additional_mapped
            if gene.feature_id in accepted_expansion_ids
        ]
        dropped = len(result.additional_mapped) - len(kept)
        result.additional_mapped = kept
        if dropped > 0:
            result.errors.append(
                f"Suppressed {dropped} expansion copy/copies after cross-gene arbitration"
            )

    _assign_deterministic_cnv_labels(mapped_genes)
    
    logger.info(f"Mapped {len(mapped_genes)} genes from {len(genes)} input genes")
    
    return mapped_genes, results


def _reference_gene_for_mapped(
    gene: Gene,
    original_genes: Dict[str, Gene]
) -> Optional[Gene]:
    """Resolve mapped gene back to the originating reference gene."""
    ref_id = gene.mapped_from
    if ref_id and ref_id in original_genes:
        return original_genes[ref_id]
    if gene.feature_id in original_genes:
        return original_genes[gene.feature_id]
    return None


def _genes_overlap_in_reference(
    gene1: Gene,
    gene2: Gene,
    original_genes: Dict[str, Gene]
) -> bool:
    """Whether two mapped genes overlapped in the source reference annotation."""
    ref1 = _reference_gene_for_mapped(gene1, original_genes)
    ref2 = _reference_gene_for_mapped(gene2, original_genes)
    if ref1 is None or ref2 is None:
        return False
    if ref1.seq_region != ref2.seq_region:
        return False
    return ref1.start <= ref2.end and ref2.start <= ref1.end


def _is_disallowed_target_overlap(
    gene1: Gene,
    gene2: Gene,
    original_genes: Dict[str, Gene],
    overlap_threshold: float
) -> bool:
    """True when two genes overlap in target but should not overlap by reference model."""
    if gene1.feature_id == gene2.feature_id:
        return False
    if gene1.seq_region != gene2.seq_region:
        return False
    if gene1.strand != gene2.strand:
        return False
    if _gene_reciprocal_overlap(gene1, gene2) < overlap_threshold:
        return False
    if _genes_overlap_in_reference(gene1, gene2, original_genes):
        return False
    return True


def _split_disallowed_overlap_components(
    group: List[Gene],
    original_genes: Dict[str, Gene],
    overlap_threshold: float
) -> List[List[Gene]]:
    """Split one physical overlap group into connected disallowed-overlap components."""
    if len(group) < 2:
        return []

    index_by_gene = {g.feature_id: i for i, g in enumerate(group)}
    adjacency: Dict[int, Set[int]] = defaultdict(set)

    for i in range(len(group)):
        for j in range(i + 1, len(group)):
            g1 = group[i]
            g2 = group[j]
            if _is_disallowed_target_overlap(g1, g2, original_genes, overlap_threshold):
                adjacency[i].add(j)
                adjacency[j].add(i)

    components: List[List[Gene]] = []
    visited: Set[int] = set()

    for i in range(len(group)):
        if i in visited or not adjacency.get(i):
            continue
        stack = [i]
        comp_indices = []
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            comp_indices.append(cur)
            for nxt in adjacency.get(cur, set()):
                if nxt not in visited:
                    stack.append(nxt)

        if len(comp_indices) > 1:
            comp = [group[idx] for idx in sorted(comp_indices)]
            components.append(comp)

    return components


def _find_disallowed_overlap_groups(
    mapped_genes: Dict[str, Gene],
    original_genes: Dict[str, Gene],
    overlap_threshold: float
) -> List[List[Gene]]:
    """Find groups of mapped genes that overlap in target but not in reference."""
    groups: List[List[Gene]] = []
    genes_by_chr: Dict[str, List[Gene]] = defaultdict(list)
    for gene in mapped_genes.values():
        genes_by_chr[gene.seq_region].append(gene)

    for genes in genes_by_chr.values():
        if len(genes) < 2:
            continue
        sorted_genes = sorted(genes, key=lambda g: g.start)

        current_group: List[Gene] = []
        current_end = -1
        for gene in sorted_genes:
            if not current_group or gene.start <= current_end:
                current_group.append(gene)
                current_end = max(current_end, gene.end)
            else:
                groups.extend(
                    _split_disallowed_overlap_components(
                        current_group,
                        original_genes,
                        overlap_threshold,
                    )
                )
                current_group = [gene]
                current_end = gene.end

        if current_group:
            groups.extend(
                _split_disallowed_overlap_components(
                    current_group,
                    original_genes,
                    overlap_threshold,
                )
            )

    return groups


def _count_disallowed_overlaps(
    mapped_genes: Dict[str, Gene],
    original_genes: Dict[str, Gene],
    overlap_threshold: float
) -> int:
    """Count pairwise disallowed overlaps for reporting."""
    count = 0
    genes = list(mapped_genes.values())
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            if _is_disallowed_target_overlap(
                genes[i], genes[j], original_genes, overlap_threshold
            ):
                count += 1
    return count


def _reference_sort_key(gene: Gene, original_genes: Dict[str, Gene]) -> Tuple[str, int, str]:
    """Stable sort key by reference position, falling back to mapped position."""
    ref_gene = _reference_gene_for_mapped(gene, original_genes)
    if ref_gene is not None:
        return (ref_gene.seq_region, ref_gene.start, gene.feature_id)
    return (gene.seq_region, gene.start, gene.feature_id)


def _score_paralog_candidate(
    candidate: Gene,
    current_gene: Gene,
    fixed_genes: List[Gene],
    prior_group_genes: List[Gene],
    original_genes: Dict[str, Gene],
    overlap_threshold: float
) -> Tuple[int, int, int, int, float, float, float]:
    """Score candidate placements for deterministic paralog reassignment."""
    conflict_count = 0
    for other in fixed_genes:
        if _is_disallowed_target_overlap(
            candidate, other, original_genes, overlap_threshold
        ):
            conflict_count += 1

    order_violations = 0
    for prior in prior_group_genes:
        if (
            candidate.seq_region == prior.seq_region
            and candidate.start < prior.start
            and not _genes_overlap_in_reference(candidate, prior, original_genes)
        ):
            order_violations += 1

    same_chr_with_prev = 0
    distance_to_prev = 10**12
    if prior_group_genes:
        prev = prior_group_genes[-1]
        if candidate.seq_region == prev.seq_region:
            same_chr_with_prev = 1
            distance_to_prev = abs(candidate.start - prev.start)

    locus_cov = _safe_attr_float(candidate, "cnv_locus_coverage", default=1.0)
    locus_score = _safe_attr_float(candidate, "cnv_locus_score", default=0.0)
    identity = candidate.mapping_identity or 0.0

    same_chr_current = 1 if candidate.seq_region == current_gene.seq_region else 0
    distance_from_current = abs(candidate.start - current_gene.start) if same_chr_current else 10**12

    # Higher is better for lexicographic tuple.
    return (
        -conflict_count,
        -order_violations,
        same_chr_with_prev,
        -distance_to_prev,
        locus_cov,
        max(locus_score, identity),
        -distance_from_current,
    )


def resolve_paralog_clusters_synteny_aware(
    mapped_genes: Dict[str, Gene],
    original_genes: Dict[str, Gene],
    syntenic_map: SyntenicMap,
    sex_chrom_map=None,
    ref_fasta_handler=None,
    target_fasta_handler=None,
    enable_boundary_refinement: bool = True,
    boundary_refine_min_coverage: float = 0.98,
    boundary_refine_window_bp: int = 500,
    boundary_refine_anchor_bp: int = 80,
    boundary_refine_max_unaligned_bp: int = 8,
    overlap_threshold: float = 0.50,
    max_locus_candidates: int = 8,
    max_rounds: int = 2
) -> Tuple[Dict[str, Gene], Dict]:
    """Reassign overlapping paralog clusters to alternative syntenic loci.

    This operates on already-mapped genes and remaps only the conflicting
    members of each cluster, preserving reference ordering where possible.
    """
    if not mapped_genes:
        return mapped_genes, {
            "groups_examined": 0,
            "genes_examined": 0,
            "genes_reassigned": 0,
            "rounds_run": 0,
            "conflicts_before": 0,
            "conflicts_after": 0,
        }

    mapper = FeatureMapper(
        syntenic_map,
        sex_chrom_map=sex_chrom_map,
        ref_fasta_handler=ref_fasta_handler,
        target_fasta_handler=target_fasta_handler,
        enable_boundary_refinement=enable_boundary_refinement,
        boundary_refine_min_coverage=boundary_refine_min_coverage,
        boundary_refine_window_bp=boundary_refine_window_bp,
        boundary_refine_anchor_bp=boundary_refine_anchor_bp,
        boundary_refine_max_unaligned_bp=boundary_refine_max_unaligned_bp,
    )

    resolved = dict(mapped_genes)
    total_groups_examined = 0
    total_genes_examined = 0
    total_reassigned = 0
    rounds_run = 0
    max_locus_candidates = max(1, max_locus_candidates)
    max_rounds = max(1, max_rounds)

    conflicts_before = _count_disallowed_overlaps(
        resolved, original_genes, overlap_threshold
    )

    for _ in range(max_rounds):
        groups = _find_disallowed_overlap_groups(
            resolved,
            original_genes,
            overlap_threshold,
        )
        if not groups:
            break

        rounds_run += 1
        round_reassigned = 0
        total_groups_examined += len(groups)
        total_genes_examined += sum(len(g) for g in groups)

        for group in groups:
            group_ids = {g.feature_id for g in group}
            fixed_outside_group = [
                g for gid, g in resolved.items()
                if gid not in group_ids
            ]

            ordered_group = sorted(
                [resolved[g.feature_id] for g in group if g.feature_id in resolved],
                key=lambda g: _reference_sort_key(g, original_genes),
            )
            if len(ordered_group) < 2:
                continue

            assigned_genes = [ordered_group[0]]
            next_assignments: Dict[str, Gene] = {
                ordered_group[0].feature_id: ordered_group[0]
            }

            for current in ordered_group[1:]:
                ref_gene = _reference_gene_for_mapped(current, original_genes)
                if ref_gene is None:
                    next_assignments[current.feature_id] = current
                    assigned_genes.append(current)
                    continue

                candidate_pool: List[Tuple[str, Gene]] = [("existing", current)]
                locus_candidates = mapper._collect_gene_locus_candidates(ref_gene)
                for locus in locus_candidates[:max_locus_candidates]:
                    remap_result = mapper._map_gene_single_copy(
                        ref_gene,
                        mapped_gene_id=current.feature_id,
                        candidate_blocks=locus.blocks,
                        transcript_id_suffix="",
                    )
                    if remap_result.mapped is None:
                        continue
                    cand = remap_result.mapped
                    merged_attrs = dict(current.attributes)
                    merged_attrs.update(cand.attributes)
                    cand.attributes = merged_attrs
                    cand.gene_name = current.gene_name
                    if current.gene_name:
                        cand.attributes["gene_name"] = current.gene_name
                    candidate_pool.append((f"locus:{locus.locus_id}", cand))

                best_label = "existing"
                best_gene = current
                best_score = _score_paralog_candidate(
                    current,
                    current,
                    fixed_outside_group + assigned_genes,
                    assigned_genes,
                    original_genes,
                    overlap_threshold,
                )
                for label, candidate in candidate_pool[1:]:
                    score = _score_paralog_candidate(
                        candidate,
                        current,
                        fixed_outside_group + assigned_genes,
                        assigned_genes,
                        original_genes,
                        overlap_threshold,
                    )
                    if score > best_score:
                        best_score = score
                        best_gene = candidate
                        best_label = label

                if best_label != "existing":
                    round_reassigned += 1
                    best_gene.attributes["paralog_reassigned"] = "true"
                    best_gene.attributes["paralog_reassignment_from"] = (
                        f"{current.seq_region}:{current.start}-{current.end}"
                    )
                    best_gene.attributes["paralog_reassignment_strategy"] = (
                        "synteny_ordered_alternative_locus"
                    )
                    best_gene.attributes["paralog_reassignment_choice"] = best_label

                next_assignments[current.feature_id] = best_gene
                assigned_genes.append(best_gene)

            resolved.update(next_assignments)

        total_reassigned += round_reassigned
        if round_reassigned == 0:
            break

    conflicts_after = _count_disallowed_overlaps(
        resolved, original_genes, overlap_threshold
    )

    report = {
        "groups_examined": total_groups_examined,
        "genes_examined": total_genes_examined,
        "genes_reassigned": total_reassigned,
        "rounds_run": rounds_run,
        "conflicts_before": conflicts_before,
        "conflicts_after": conflicts_after,
    }
    return resolved, report


def reassign_overlapping_clusters(
    mapped_genes: Dict[str, Gene],
    original_genes: Dict[str, Gene],
    syntenic_map: "SyntenicMap"
) -> Dict[str, Gene]:
    """Backward-compatible wrapper for synteny-aware paralog reassignment."""
    reassigned, _ = resolve_paralog_clusters_synteny_aware(
        mapped_genes,
        original_genes,
        syntenic_map,
    )
    return reassigned

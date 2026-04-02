"""Coordinate projection using minimap2 cs tag.

Parses the cs (difference string) tag to exactly project coordinates
from reference to target genome.
"""

import re
from dataclasses import dataclass
from typing import List, Optional, Tuple

from src.models import SyntenicBlock


@dataclass
class CSOperation:
    """A single operation from a cs tag."""
    op_type: str  # ':', '*', '+', '-', '~'
    length: int   # Number of bases affected
    ref_bases: str = ""  # Reference bases (for substitutions/deletions)
    query_bases: str = ""  # Query bases (for substitutions/insertions)


def parse_cs_tag(cs_tag: str) -> List[CSOperation]:
    """Parse a minimap2 cs tag into operations.
    
    cs tag format (short form):
        :N      - N bases identical
        *XY     - substitution (ref X -> query Y)
        +SEQ    - insertion in query
        -SEQ    - deletion in query
        ~nnXnn  - intron (for splice alignments)
    
    cs tag format (long form with --cs=long):
        =SEQ    - identical bases
        *XY     - substitution
        +SEQ    - insertion
        -SEQ    - deletion
        ~nnXnn  - intron
    
    Args:
        cs_tag: The cs tag string
        
    Returns:
        List of CSOperation objects
    """
    operations = []
    
    # Pattern for long-form cs tag (used with --cs=long)
    # =ACGT (match), *ac (sub), +acgt (ins), -acgt (del)
    pattern = re.compile(r'([:=])(\d+|[ACGTNacgtn]+)|(\*)([ACGTNacgtn])([ACGTNacgtn])|([+-])([ACGTNacgtn]+)')
    
    pos = 0
    while pos < len(cs_tag):
        # Match operation
        match_type = None
        
        # Check for match/identity
        if cs_tag[pos] == ':':
            # Short form: :N (N identical bases)
            end = pos + 1
            while end < len(cs_tag) and cs_tag[end].isdigit():
                end += 1
            length = int(cs_tag[pos+1:end])
            operations.append(CSOperation(':', length))
            pos = end
            continue
        
        elif cs_tag[pos] == '=':
            # Long form: =SEQ (identical sequence)
            end = pos + 1
            while end < len(cs_tag) and cs_tag[end] in 'ACGTNacgtn':
                end += 1
            seq = cs_tag[pos+1:end]
            operations.append(CSOperation(':', len(seq), seq.upper(), seq.upper()))
            pos = end
            continue
        
        elif cs_tag[pos] == '*':
            # Substitution: *XY
            if pos + 2 < len(cs_tag):
                ref_base = cs_tag[pos+1].upper()
                query_base = cs_tag[pos+2].upper()
                operations.append(CSOperation('*', 1, ref_base, query_base))
                pos += 3
            else:
                pos += 1
            continue
        
        elif cs_tag[pos] == '+':
            # Insertion: +SEQ
            end = pos + 1
            while end < len(cs_tag) and cs_tag[end] in 'ACGTNacgtn':
                end += 1
            seq = cs_tag[pos+1:end]
            operations.append(CSOperation('+', len(seq), '', seq.upper()))
            pos = end
            continue
        
        elif cs_tag[pos] == '-':
            # Deletion: -SEQ
            end = pos + 1
            while end < len(cs_tag) and cs_tag[end] in 'ACGTNacgtn':
                end += 1
            seq = cs_tag[pos+1:end]
            operations.append(CSOperation('-', len(seq), seq.upper(), ''))
            pos = end
            continue
        
        elif cs_tag[pos] == '~':
            # Intron (splice alignment): ~nnXnn
            end = pos + 1
            while end < len(cs_tag) and (cs_tag[end].isdigit() or cs_tag[end] in 'gtag'):
                end += 1
            # For now, treat introns as large deletions
            intron_match = re.match(r'~(\d+)[a-z]+(\d+)', cs_tag[pos:end])
            if intron_match:
                intron_len = int(intron_match.group(1))
                operations.append(CSOperation('~', intron_len))
            pos = end
            continue
        
        else:
            # Unknown character, skip
            pos += 1
    
    return operations


class CoordinateProjector:
    """Projects coordinates between reference and target genomes using cs tags."""
    
    def __init__(self, block: SyntenicBlock):
        """Initialize projector with a syntenic block.
        
        Args:
            block: The syntenic block containing alignment info
        """
        self.block = block
        # Use pre-cached offset map if available, otherwise parse
        if block.cached_offset_map is not None:
            self._offset_changes = block.cached_offset_map
            self._has_cs_tag = True
        elif block.cs_tag:
            # Fallback: parse cs tag (should not happen if build_index was called)
            from src.models import SyntenicMap
            self._offset_changes = SyntenicMap._parse_cs_to_offsets(block)
            self._has_cs_tag = True
        else:
            self._offset_changes = []
            self._has_cs_tag = False
    
    def _get_offset_at(self, ref_pos: int) -> int:
        """Get the cumulative offset at a reference position.
        
        Uses binary search for efficiency.
        """
        if not self._offset_changes:
            return 0
        
        # Binary search for largest ref_pos <= query
        lo, hi = 0, len(self._offset_changes)
        while lo < hi:
            mid = (lo + hi) // 2
            if self._offset_changes[mid][0] <= ref_pos:
                lo = mid + 1
            else:
                hi = mid
        
        if lo == 0:
            return 0
        return self._offset_changes[lo - 1][1]
    
    def project_position(self, ref_pos: int) -> Optional[int]:
        """Project a single position from reference to target.
        
        Args:
            ref_pos: Position in reference (1-based)
            
        Returns:
            Position in target (1-based), or None if position is in a deletion
        """
        # Check if position is within block
        if not self.block.ref_interval.contains(ref_pos):
            return None
        
        if self._has_cs_tag:
            # Calculate target position using cumulative offset
            relative_ref = ref_pos - self.block.ref_interval.start
            offset = self._get_offset_at(ref_pos)
            
            if self.block.is_inverted:
                # For inversions, target genomic coordinates run opposite to
                # alignment/query progression. Apply cs-derived offset with
                # opposite sign versus forward orientation.
                target_pos = self.block.target_interval.end - relative_ref - offset
            else:
                target_pos = self.block.target_interval.start + relative_ref + offset
            
            return target_pos
        
        # Fallback: linear interpolation without cs tag
        return self._linear_project(ref_pos)
    
    def _linear_project(self, ref_pos: int) -> int:
        """Simple linear projection without cs tag.
        
        Args:
            ref_pos: Reference position
            
        Returns:
            Target position (approximate)
        """
        ref_interval = self.block.ref_interval
        target_interval = self.block.target_interval
        
        # Calculate relative position
        relative_pos = (ref_pos - ref_interval.start) / max(ref_interval.length - 1, 1)
        
        if self.block.is_inverted:
            # Inverted: end of ref maps to start of target
            target_pos = target_interval.end - int(relative_pos * (target_interval.length - 1))
        else:
            # Normal: start of ref maps to start of target
            target_pos = target_interval.start + int(relative_pos * (target_interval.length - 1))
        
        return target_pos
    
    def project_interval(
        self,
        ref_start: int,
        ref_end: int
    ) -> Optional[Tuple[int, int]]:
        """Project an interval from reference to target.
        
        Args:
            ref_start: Start position in reference (1-based)
            ref_end: End position in reference (1-based)
            
        Returns:
            Tuple of (target_start, target_end), or None if interval can't be projected
        """
        target_start = self.project_position(ref_start)
        target_end = self.project_position(ref_end)
        
        if target_start is None or target_end is None:
            return None
        
        # Handle inversions
        if target_start > target_end:
            target_start, target_end = target_end, target_start
        
        return (target_start, target_end)


def _reverse_complement(seq: str) -> str:
    """Reverse-complement an uppercase/lowercase DNA sequence."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _boundary_refinement_needed(
    coverage: float,
    clamped_5prime: bool,
    clamped_3prime: bool,
    min_coverage: float
) -> bool:
    """Whether feature boundaries should be locally refined."""
    return clamped_5prime or clamped_3prime or coverage < min_coverage


def _boundary_anchor_from_feature(
    ref_start: int,
    ref_end: int,
    anchor_len: int,
    is_five_prime: bool
) -> Optional[Tuple[int, int, int]]:
    """Get reference anchor coordinates and boundary index in anchor."""
    if ref_start > ref_end:
        return None
    if anchor_len < 8:
        anchor_len = 8
    
    if is_five_prime:
        anchor_start = ref_start
        anchor_end = min(ref_end, ref_start + anchor_len - 1)
        boundary_idx = 0
    else:
        anchor_start = max(ref_start, ref_end - anchor_len + 1)
        anchor_end = ref_end
        boundary_idx = anchor_end - anchor_start
    
    if anchor_start > anchor_end:
        return None
    return anchor_start, anchor_end, boundary_idx


def _map_query_index_to_target_index(
    query_index: int,
    aligned_query,
    aligned_target,
    query_len: int,
    target_len: int,
    max_boundary_unaligned: int
) -> Optional[int]:
    """Map query index to target index from a local alignment."""
    if len(aligned_query) == 0 or len(aligned_target) == 0:
        return None
    
    for (qs, qe), (ts, te) in zip(aligned_query, aligned_target):
        qs_i, qe_i = int(qs), int(qe)
        ts_i, te_i = int(ts), int(te)
        if qs_i <= query_index < qe_i:
            return ts_i + (query_index - qs_i)
    
    # Allow small unclipped boundary offsets for local alignment.
    first_q_start = int(aligned_query[0][0])
    first_t_start = int(aligned_target[0][0])
    last_q_end = int(aligned_query[-1][1])
    last_t_end = int(aligned_target[-1][1])
    
    if query_index == 0:
        if first_q_start <= max_boundary_unaligned:
            return max(0, first_t_start - first_q_start)
        return None
    
    if query_index == query_len - 1:
        tail = query_len - last_q_end
        if tail <= max_boundary_unaligned:
            return min(target_len - 1, last_t_end + tail - 1)
        return None
    
    return None


def _refine_boundary_with_local_alignment(
    ref_fasta,
    target_fasta,
    ref_seq_region: str,
    ref_start: int,
    ref_end: int,
    boundary_ref_pos: int,
    target_seq_region: str,
    boundary_target_pos: int,
    is_inverted: bool,
    is_five_prime: bool,
    window_bp: int,
    anchor_len: int,
    max_boundary_unaligned: int
) -> Optional[int]:
    """Refine one projected boundary using local base-level sequence alignment."""
    anchor_info = _boundary_anchor_from_feature(
        ref_start, ref_end, anchor_len=anchor_len, is_five_prime=is_five_prime
    )
    if anchor_info is None:
        return None
    
    anchor_ref_start, anchor_ref_end, boundary_idx = anchor_info
    ref_anchor = ref_fasta.fetch(ref_seq_region, anchor_ref_start, anchor_ref_end)
    if not ref_anchor:
        return None
    
    query = _reverse_complement(ref_anchor) if is_inverted else ref_anchor
    query_len = len(query)
    if query_len < 8:
        return None
    
    if is_inverted:
        boundary_idx = query_len - 1 - boundary_idx
    
    target_chr_len = target_fasta.get_length(target_seq_region)
    if target_chr_len <= 0:
        return None
    
    window_start = max(1, boundary_target_pos - window_bp)
    window_end = min(target_chr_len, boundary_target_pos + window_bp)
    if window_start > window_end:
        return None
    
    target_window = target_fasta.fetch(target_seq_region, window_start, window_end)
    if not target_window:
        return None
    
    # Fast path: exact-anchor occurrences near projected boundary.
    candidate_positions = []
    search_from = 0
    while True:
        hit = target_window.find(query, search_from)
        if hit == -1:
            break
        boundary_pos = window_start + hit + boundary_idx
        candidate_positions.append(boundary_pos)
        search_from = hit + 1
    
    if candidate_positions:
        return min(candidate_positions, key=lambda p: abs(p - boundary_target_pos))
    
    # Fallback: local alignment around uncertain edge.
    try:
        from Bio.Align import PairwiseAligner
    except Exception:
        return None
    
    aligner = PairwiseAligner(mode="local")
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    
    alignments = aligner.align(query, target_window)
    if not alignments:
        return None
    
    best = alignments[0]
    mapped_idx = _map_query_index_to_target_index(
        boundary_idx,
        best.aligned[0],
        best.aligned[1],
        query_len=query_len,
        target_len=len(target_window),
        max_boundary_unaligned=max_boundary_unaligned,
    )
    if mapped_idx is None:
        return None
    
    return window_start + mapped_idx


@dataclass
class ProjectionSegment:
    """A projected sub-interval of the feature through one block."""
    block: SyntenicBlock
    ref_start: int
    ref_end: int
    target_start: int
    target_end: int

    @property
    def identity(self) -> float:
        return self.block.identity

    @property
    def target_seq_region(self) -> str:
        return self.block.target_interval.seq_region

    @property
    def is_inverted(self) -> bool:
        return self.block.is_inverted

    @property
    def ref_length(self) -> int:
        return self.ref_end - self.ref_start + 1


def _build_projection_segments(
    ref_seq_region: str,
    ref_start: int,
    ref_end: int,
    blocks: List[SyntenicBlock],
    require_same_chromosome: bool
) -> List[ProjectionSegment]:
    """Build candidate projected segments for a reference interval."""
    segments: List[ProjectionSegment] = []

    for block in blocks:
        if block.ref_interval.seq_region != ref_seq_region:
            continue
        if block.ref_interval.start > ref_end or ref_start > block.ref_interval.end:
            continue
        if require_same_chromosome and block.target_interval.seq_region != ref_seq_region:
            continue

        seg_ref_start = max(ref_start, block.ref_interval.start)
        seg_ref_end = min(ref_end, block.ref_interval.end)
        if seg_ref_start > seg_ref_end:
            continue

        projector = CoordinateProjector(block)
        projected = projector.project_interval(seg_ref_start, seg_ref_end)
        if projected is None:
            continue

        segments.append(
            ProjectionSegment(
                block=block,
                ref_start=seg_ref_start,
                ref_end=seg_ref_end,
                target_start=projected[0],
                target_end=projected[1],
            )
        )

    return segments


def _path_compatible(prev: ProjectionSegment, cur: ProjectionSegment) -> bool:
    """Check if two projected segments can be chained in a syntenic path."""
    if prev.target_seq_region != cur.target_seq_region:
        return False
    if prev.is_inverted != cur.is_inverted:
        return False
    if cur.ref_start <= prev.ref_start:
        return False

    # Allow small local backtracking in repetitive regions, but reject
    # strong order violations.
    max_backtrack = 5000
    if not prev.is_inverted:
        return cur.target_start + max_backtrack >= prev.target_start
    return cur.target_end - max_backtrack <= prev.target_end


def _transition_penalty(prev: ProjectionSegment, cur: ProjectionSegment) -> float:
    """Penalty for chaining two projected segments."""
    ref_gap = max(0, cur.ref_start - prev.ref_end - 1)
    ref_overlap = max(0, prev.ref_end - cur.ref_start + 1)

    if not prev.is_inverted:
        target_gap = max(0, cur.target_start - prev.target_end - 1)
        target_overlap = max(0, prev.target_end - cur.target_start + 1)
    else:
        target_gap = max(0, prev.target_start - cur.target_end - 1)
        target_overlap = max(0, cur.target_end - prev.target_start + 1)

    gap_mismatch = abs(ref_gap - target_gap)
    return (
        0.8 * ref_overlap +
        0.4 * target_overlap +
        0.02 * gap_mismatch
    )


def _select_best_synteny_path(segments: List[ProjectionSegment]) -> List[ProjectionSegment]:
    """Select best multi-block path using DP over synteny-compatible segments."""
    if not segments:
        return []

    segments = sorted(
        segments,
        key=lambda s: (s.ref_start, s.ref_end, -s.identity)
    )
    n = len(segments)

    dp = [0.0] * n
    prev_idx = [-1] * n

    for i in range(n):
        seg = segments[i]
        base_score = seg.ref_length * seg.identity
        dp[i] = base_score

        for j in range(i):
            prev = segments[j]
            if not _path_compatible(prev, seg):
                continue

            score = dp[j] + base_score - _transition_penalty(prev, seg)
            if score > dp[i]:
                dp[i] = score
                prev_idx[i] = j

    best_i = max(range(n), key=lambda i: dp[i])
    path = []
    cur = best_i
    while cur != -1:
        path.append(segments[cur])
        cur = prev_idx[cur]
    path.reverse()
    return path


def _coverage_from_path(path: List[ProjectionSegment], ref_start: int, ref_end: int) -> float:
    """Fraction of feature bases covered by path segments."""
    if not path:
        return 0.0

    intervals = sorted((s.ref_start, s.ref_end) for s in path)
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    covered = sum(end - start + 1 for start, end in merged)
    total = max(1, ref_end - ref_start + 1)
    return covered / total


def _project_boundary_from_path(
    ref_pos: int,
    path: List[ProjectionSegment],
    use_first_fallback: bool
) -> Optional[int]:
    """Project one boundary coordinate using the selected path.

    When *ref_pos* falls inside a covered segment, the cs-tag offset map is
    used directly.  When it falls in a **gap** between two consecutive path
    segments the function interpolates through the gap, accounting for the
    net indel (target_gap − ref_gap) so that downstream coordinates are
    shifted by the correct amount.
    """
    if not path:
        return None

    # --- 1. Direct hit inside a segment ---
    for seg in path:
        if seg.ref_start <= ref_pos <= seg.ref_end:
            projector = CoordinateProjector(seg.block)
            projected = projector.project_position(ref_pos)
            if projected is not None:
                return projected
            # Position is within segment span but lands on a deleted base;
            # fall through to edge/interpolation below.
            break

    # --- 2. Position falls in a gap between two consecutive segments ---
    for i in range(len(path) - 1):
        seg_before = path[i]
        seg_after = path[i + 1]

        # Check whether ref_pos is in the gap between these two segments.
        if seg_before.ref_end < ref_pos < seg_after.ref_start:
            # Project the edges of the flanking segments.
            proj_before = CoordinateProjector(seg_before.block)
            proj_after = CoordinateProjector(seg_after.block)

            target_before = proj_before.project_position(seg_before.ref_end)
            target_after = proj_after.project_position(seg_after.ref_start)
            if target_before is None or target_after is None:
                break  # can't interpolate — fall through to clamp

            ref_gap = seg_after.ref_start - seg_before.ref_end  # always > 0
            frac = (ref_pos - seg_before.ref_end) / ref_gap

            if not seg_before.is_inverted:
                target_gap = target_after - target_before
            else:
                target_gap = target_before - target_after

            interpolated_offset = int(round(frac * target_gap))

            if not seg_before.is_inverted:
                return target_before + interpolated_offset
            else:
                return target_before - interpolated_offset

    # --- 3. Position is outside the path span (before first / after last) ---
    # Clamp to the nearest path edge.
    if ref_pos < path[0].ref_start:
        selected = path[0]
        edge_ref = selected.ref_start
    elif ref_pos > path[-1].ref_end:
        selected = path[-1]
        edge_ref = selected.ref_end
    else:
        # Inside overall span but not in any segment or gap (shouldn't
        # happen with a well-formed path, but be defensive).
        selected = path[0] if use_first_fallback else path[-1]
        edge_ref = min(max(ref_pos, selected.ref_start), selected.ref_end)

    projector = CoordinateProjector(selected.block)
    projected = projector.project_position(edge_ref)
    if projected is not None:
        return projected

    # Final fallback to segment edge
    fallback_ref = selected.ref_start if use_first_fallback else selected.ref_end
    return projector.project_position(fallback_ref)


def project_coordinate_through_blocks(
    ref_seq_region: str,
    ref_pos: int,
    blocks: List[SyntenicBlock]
) -> Optional[Tuple[str, int]]:
    """Project a coordinate through any matching syntenic block.
    
    Args:
        ref_seq_region: Reference chromosome/scaffold
        ref_pos: Reference position (1-based)
        blocks: List of syntenic blocks
        
    Returns:
        Tuple of (target_seq_region, target_pos), or None if not in any block
    """
    for block in blocks:
        if block.ref_interval.seq_region != ref_seq_region:
            continue
        
        if block.ref_interval.contains(ref_pos):
            projector = CoordinateProjector(block)
            target_pos = projector.project_position(ref_pos)
            
            if target_pos is not None:
                return (block.target_interval.seq_region, target_pos)
    
    return None


def select_synteny_path_blocks(
    ref_seq_region: str,
    ref_start: int,
    ref_end: int,
    blocks: List[SyntenicBlock],
    require_same_chromosome: bool = False,
) -> List[SyntenicBlock]:
    """Select the best synteny-consistent block path for one reference interval."""
    segments = _build_projection_segments(
        ref_seq_region,
        ref_start,
        ref_end,
        blocks,
        require_same_chromosome=require_same_chromosome,
    )
    if not segments:
        return []

    path = _select_best_synteny_path(segments)
    if not path:
        return []
    return [seg.block for seg in path]


def project_feature_coordinates(
    ref_seq_region: str,
    ref_start: int,
    ref_end: int,
    blocks: List[SyntenicBlock],
    require_same_chromosome: bool = False,
    ref_fasta=None,
    target_fasta=None,
    enable_boundary_refinement: bool = True,
    boundary_refine_min_coverage: float = 0.98,
    boundary_refine_window_bp: int = 500,
    boundary_refine_anchor_bp: int = 80,
    boundary_refine_max_unaligned_bp: int = 8,
) -> Optional[Tuple[str, int, int, bool, dict]]:
    """Project feature coordinates to target genome.
    
    Args:
        ref_seq_region: Reference chromosome/scaffold
        ref_start: Reference start (1-based)
        ref_end: Reference end (1-based)
        blocks: List of syntenic blocks
        require_same_chromosome: If True, only use blocks that map to a target
            chromosome matching the reference chromosome name
        
    Returns:
        Tuple of (target_seq_region, target_start, target_end, is_inverted, info_dict),
        or None if feature can't be projected.
        
        info_dict contains:
            - 'coverage': fraction of feature covered by block (0.0-1.0)
            - 'clamped_5prime': True if 5' end was outside block and clamped
            - 'clamped_3prime': True if 3' end was outside block and clamped
            - 'cross_chromosome': True if target is different chromosome
            - 'boundary_refinement_triggered': refinement attempted for uncertain edges
            - 'boundary_refined_5prime': 5' boundary changed by local refinement
            - 'boundary_refined_3prime': 3' boundary changed by local refinement
    """
    # Build candidate projected segments
    segments = _build_projection_segments(
        ref_seq_region,
        ref_start,
        ref_end,
        blocks,
        require_same_chromosome=require_same_chromosome,
    )
    
    if not segments:
        return None
    
    overlapping = [s.block for s in segments]

    # Find all blocks that fully contain this feature
    containing_blocks = [
        block for block in overlapping
        if block.ref_interval.start <= ref_start and ref_end <= block.ref_interval.end
    ]
    
    # If feature is contained in one or more blocks, use the best by identity
    if containing_blocks:
        # Sort by identity (highest first) to pick the best mapping
        best_block = max(containing_blocks, key=lambda b: b.identity)
        projector = CoordinateProjector(best_block)
        result = projector.project_interval(ref_start, ref_end)
        if result:
            info = {
                'coverage': 1.0,
                'clamped_5prime': False,
                'clamped_3prime': False,
                'cross_chromosome': best_block.target_interval.seq_region != ref_seq_region
            }
            return (
                best_block.target_interval.seq_region,
                result[0],
                result[1],
                best_block.is_inverted,
                info
            )
    
    # Multi-block path projection for features spanning block boundaries.
    path = _select_best_synteny_path(segments)
    if not path:
        return None

    projected_start = _project_boundary_from_path(ref_start, path, use_first_fallback=True)
    projected_end = _project_boundary_from_path(ref_end, path, use_first_fallback=False)
    if projected_start is None or projected_end is None:
        return None

    coverage = _coverage_from_path(path, ref_start, ref_end)
    clamped_5prime = not any(s.ref_start <= ref_start <= s.ref_end for s in path)
    clamped_3prime = not any(s.ref_start <= ref_end <= s.ref_end for s in path)
    path_chr = path[0].target_seq_region
    path_inverted = path[0].is_inverted
    
    refined_5prime = False
    refined_3prime = False
    refinement_triggered = False
    refinement_possible = (
        enable_boundary_refinement
        and ref_fasta is not None
        and target_fasta is not None
    )
    if refinement_possible and _boundary_refinement_needed(
        coverage,
        clamped_5prime,
        clamped_3prime,
        min_coverage=boundary_refine_min_coverage,
    ):
        refinement_triggered = True
        
        refined_start = _refine_boundary_with_local_alignment(
            ref_fasta=ref_fasta,
            target_fasta=target_fasta,
            ref_seq_region=ref_seq_region,
            ref_start=ref_start,
            ref_end=ref_end,
            boundary_ref_pos=ref_start,
            target_seq_region=path_chr,
            boundary_target_pos=projected_start,
            is_inverted=path_inverted,
            is_five_prime=True,
            window_bp=boundary_refine_window_bp,
            anchor_len=boundary_refine_anchor_bp,
            max_boundary_unaligned=boundary_refine_max_unaligned_bp,
        )
        if refined_start is not None and refined_start != projected_start:
            projected_start = refined_start
            refined_5prime = True
        
        refined_end = _refine_boundary_with_local_alignment(
            ref_fasta=ref_fasta,
            target_fasta=target_fasta,
            ref_seq_region=ref_seq_region,
            ref_start=ref_start,
            ref_end=ref_end,
            boundary_ref_pos=ref_end,
            target_seq_region=path_chr,
            boundary_target_pos=projected_end,
            is_inverted=path_inverted,
            is_five_prime=False,
            window_bp=boundary_refine_window_bp,
            anchor_len=boundary_refine_anchor_bp,
            max_boundary_unaligned=boundary_refine_max_unaligned_bp,
        )
        if refined_end is not None and refined_end != projected_end:
            projected_end = refined_end
            refined_3prime = True
    
    target_start = min(projected_start, projected_end)
    target_end = max(projected_start, projected_end)

    info = {
        'coverage': coverage,
        'clamped_5prime': clamped_5prime,
        'clamped_3prime': clamped_3prime,
        'cross_chromosome': path_chr != ref_seq_region,
        'path_blocks': len(path),
        'path_ref_span': sum(s.ref_length for s in path),
        'boundary_refinement_triggered': refinement_triggered,
        'boundary_refined_5prime': refined_5prime,
        'boundary_refined_3prime': refined_3prime,
    }
    
    return (
        path_chr,
        target_start,
        target_end,
        path_inverted,
        info
    )

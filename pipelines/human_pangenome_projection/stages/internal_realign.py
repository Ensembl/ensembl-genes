"""Targeted internal realignment rescue for large-indel continuity issues."""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from src.models import GenomicInterval, Strand, SyntenicBlock
from stages.coordinate_projection import select_synteny_path_blocks
from stages.synteny import PAFRecord

logger = logging.getLogger(__name__)


@dataclass
class InternalRealignResult:
    """Result from one targeted internal realignment attempt."""

    blocks: List[SyntenicBlock]
    backend: str
    path_coverage: float
    internal_gap_bp: int
    reason: Optional[str] = None


def _reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def _merge_intervals(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(s, e) for s, e in merged]


def path_coverage_and_gap(
    path_blocks: Sequence[SyntenicBlock],
    ref_start: int,
    ref_end: int,
) -> Tuple[float, int]:
    """Coverage and maximum uncovered internal gap on reference span."""
    if ref_end < ref_start:
        return 0.0, 0

    overlaps: List[Tuple[int, int]] = []
    for block in path_blocks:
        ov_start = max(ref_start, block.ref_interval.start)
        ov_end = min(ref_end, block.ref_interval.end)
        if ov_end >= ov_start:
            overlaps.append((ov_start, ov_end))

    merged = _merge_intervals(overlaps)
    if not merged:
        return 0.0, max(0, ref_end - ref_start + 1)

    covered = sum(end - start + 1 for start, end in merged)
    total = max(1, ref_end - ref_start + 1)
    coverage = covered / total

    max_gap = 0
    cur = ref_start
    for start, end in merged:
        if start > cur:
            max_gap = max(max_gap, start - cur)
        cur = max(cur, end + 1)
    if cur <= ref_end:
        max_gap = max(max_gap, ref_end - cur + 1)

    return coverage, max_gap


def _write_single_fasta(path: Path, name: str, seq: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f">{name}\n")
        for i in range(0, len(seq), 120):
            handle.write(seq[i : i + 120] + "\n")


def _run_targeted_minimap2(
    ref_seq: str,
    target_seq: str,
    threads: int,
    temp_dir: Optional[Path] = None,
) -> List[PAFRecord]:
    """Align one reference window to one target window with minimap2."""
    if not ref_seq or not target_seq:
        return []

    with tempfile.TemporaryDirectory(
        prefix="internal_realign_",
        dir=str(temp_dir) if temp_dir else None,
    ) as tmpdir:
        tmp = Path(tmpdir)
        ref_fa = tmp / "ref_window.fa"
        target_fa = tmp / "target_window.fa"
        _write_single_fasta(ref_fa, "ref_window", ref_seq)
        _write_single_fasta(target_fa, "target_window", target_seq)

        cmd = [
            "minimap2",
            "-c",
            "--cs=long",
            "-x",
            "asm5",
            "-N",
            "8",
            "-t",
            str(max(1, threads)),
            str(target_fa),
            str(ref_fa),
        ]

        try:
            run = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )
        except FileNotFoundError:
            logger.warning("minimap2 not available for internal realignment")
            return []

        if run.returncode != 0:
            stderr = (run.stderr or "").strip().splitlines()
            detail = stderr[0] if stderr else "unknown minimap2 error"
            logger.warning("Internal minimap2 failed: %s", detail)
            return []

        records: List[PAFRecord] = []
        for line in run.stdout.splitlines():
            if not line.strip():
                continue
            rec = PAFRecord.from_line(line)
            if rec is not None:
                records.append(rec)
        return records


def _record_to_local_block(
    record: PAFRecord,
    ref_chr: str,
    target_chr: str,
    ref_window_start: int,
    target_window_start: int,
) -> Optional[SyntenicBlock]:
    """Convert local-window PAF record to absolute-coordinate syntenic block."""
    ref_start = ref_window_start + record.query_start
    ref_end = ref_window_start + record.query_end - 1
    target_start = target_window_start + record.target_start
    target_end = target_window_start + record.target_end - 1

    if ref_start > ref_end or target_start > target_end:
        return None

    return SyntenicBlock(
        ref_interval=GenomicInterval(
            seq_region=ref_chr,
            start=ref_start,
            end=ref_end,
            strand=Strand.PLUS,
        ),
        target_interval=GenomicInterval(
            seq_region=target_chr,
            start=target_start,
            end=target_end,
            strand=Strand.from_string(record.strand),
        ),
        identity=record.identity,
        alignment_length=record.alignment_length,
        matches=record.matches,
        mismatches=max(0, record.alignment_length - record.matches),
        cs_tag=record.cs_tag,
        cigar=record.cigar,
        block_id=f"internal_mm2_{ref_start}_{ref_end}_{target_start}_{target_end}",
    )


def _select_minimap2_path(
    records: Sequence[PAFRecord],
    ref_chr: str,
    target_chr: str,
    ref_window_start: int,
    target_window_start: int,
    gene_start: int,
    gene_end: int,
    min_path_coverage: float,
    max_internal_gap_bp: int,
) -> Optional[InternalRealignResult]:
    blocks = []
    for rec in records:
        block = _record_to_local_block(
            rec,
            ref_chr=ref_chr,
            target_chr=target_chr,
            ref_window_start=ref_window_start,
            target_window_start=target_window_start,
        )
        if block is None:
            continue
        blocks.append(block)

    if not blocks:
        return None

    path = select_synteny_path_blocks(
        ref_chr,
        gene_start,
        gene_end,
        blocks,
    )
    if not path:
        return None

    target_chrs = {b.target_interval.seq_region for b in path}
    orientations = {b.is_inverted for b in path}
    if len(target_chrs) != 1 or len(orientations) != 1:
        return InternalRealignResult(
            blocks=[],
            backend="minimap2_cs",
            path_coverage=0.0,
            internal_gap_bp=max_internal_gap_bp + 1,
            reason="inconsistent_path_orientation_or_chromosome",
        )

    coverage, max_gap = path_coverage_and_gap(path, gene_start, gene_end)
    if coverage < min_path_coverage:
        return InternalRealignResult(
            blocks=[],
            backend="minimap2_cs",
            path_coverage=coverage,
            internal_gap_bp=max_gap,
            reason="insufficient_path_coverage",
        )
    if max_gap > max_internal_gap_bp:
        return InternalRealignResult(
            blocks=[],
            backend="minimap2_cs",
            path_coverage=coverage,
            internal_gap_bp=max_gap,
            reason="internal_gap_too_large",
        )

    return InternalRealignResult(
        blocks=path,
        backend="minimap2_cs",
        path_coverage=coverage,
        internal_gap_bp=max_gap,
    )


def _parse_cigar_ops(cigar: str) -> List[Tuple[int, str]]:
    return [(int(length), op) for length, op in re.findall(r"(\d+)([=XIDM])", cigar or "")]


def _compress_offset_map(changes: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not changes:
        return []
    changes.sort(key=lambda item: item[0])
    compressed: List[Tuple[int, int]] = []
    for pos, off in changes:
        if not compressed:
            compressed.append((pos, off))
            continue
        last_pos, last_off = compressed[-1]
        if pos == last_pos:
            compressed[-1] = (pos, off)
            continue
        if off == last_off:
            continue
        compressed.append((pos, off))
    return compressed


def _offset_map_from_cigar(
    cigar: str,
    ref_start: int,
    initial_offset: int,
) -> Tuple[List[Tuple[int, int]], int, int, int]:
    """Build sparse offset map from extended CIGAR.

    Returns:
        (offset_changes, query_consumed, target_consumed, max_query_gap)
    """
    ref_pos = ref_start
    offset = initial_offset
    q_consumed = 0
    t_consumed = 0
    max_query_gap = 0
    changes: List[Tuple[int, int]] = [(ref_pos, offset)]

    for length, op in _parse_cigar_ops(cigar):
        if op in {"=", "X", "M"}:
            ref_pos += length
            q_consumed += length
            t_consumed += length
            continue
        if op == "I":
            offset += length
            t_consumed += length
            changes.append((ref_pos, offset))
            continue
        if op == "D":
            ref_pos += length
            q_consumed += length
            offset -= length
            max_query_gap = max(max_query_gap, length)
            changes.append((ref_pos, offset))
            continue

    return _compress_offset_map(changes), q_consumed, t_consumed, max_query_gap


def _edlib_block(
    ref_seq: str,
    target_seq: str,
    ref_chr: str,
    target_chr: str,
    ref_window_start: int,
    ref_window_end: int,
    target_window_start: int,
    target_window_end: int,
) -> Optional[InternalRealignResult]:
    """Fallback internal realignment using edlib (or local safe fallback)."""
    try:
        import edlib  # type: ignore
    except Exception:
        return None

    best = None
    for is_rev in (False, True):
        oriented_target = _reverse_complement(target_seq) if is_rev else target_seq
        result = edlib.align(ref_seq, oriented_target, mode="HW", task="path")
        edit_distance = int(result.get("editDistance", -1))
        locations = result.get("locations") or []
        cigar = result.get("cigar")
        if edit_distance < 0 or not locations or not cigar:
            continue

        loc_start, _loc_end = locations[0]
        offset_map, q_used, t_used, max_q_gap = _offset_map_from_cigar(
            cigar,
            ref_start=ref_window_start,
            initial_offset=int(loc_start),
        )
        if q_used <= 0:
            continue

        identity = max(0.0, 1.0 - (edit_distance / max(1, q_used)))
        strand = Strand.MINUS if is_rev else Strand.PLUS
        block = SyntenicBlock(
            ref_interval=GenomicInterval(
                seq_region=ref_chr,
                start=ref_window_start,
                end=ref_window_end,
                strand=Strand.PLUS,
            ),
            target_interval=GenomicInterval(
                seq_region=target_chr,
                start=target_window_start,
                end=target_window_end,
                strand=strand,
            ),
            identity=identity,
            alignment_length=max(q_used, t_used),
            matches=max(0, max(q_used, t_used) - edit_distance),
            mismatches=max(0, edit_distance),
            block_id="internal_edlib",
        )
        block.cached_offset_map = offset_map

        score = (identity, -edit_distance)
        candidate = InternalRealignResult(
            blocks=[block],
            backend="edlib",
            path_coverage=q_used / max(1, len(ref_seq)),
            internal_gap_bp=max_q_gap,
        )
        if best is None or score > best[0]:
            best = (score, candidate)

    return None if best is None else best[1]


def _pairwise_block(
    ref_seq: str,
    target_seq: str,
    ref_chr: str,
    target_chr: str,
    ref_window_start: int,
    ref_window_end: int,
    target_window_start: int,
    target_window_end: int,
) -> Optional[InternalRealignResult]:
    """Safe fallback when edlib is unavailable."""
    try:
        from Bio import pairwise2
    except Exception:
        return None

    best = None
    for is_rev in (False, True):
        oriented_target = _reverse_complement(target_seq) if is_rev else target_seq
        alignments = pairwise2.align.globalms(
            ref_seq,
            oriented_target,
            2,
            -3,
            -6,
            -1,
            one_alignment_only=True,
            penalize_end_gaps=(False, False),
        )
        if not alignments:
            continue

        aln = alignments[0]
        aligned_ref = aln.seqA
        aligned_target = aln.seqB

        ref_pos = ref_window_start
        offset = 0
        matches = 0
        q_used = 0
        t_used = 0
        del_run = 0
        max_q_gap = 0
        changes: List[Tuple[int, int]] = [(ref_pos, offset)]

        for rb, tb in zip(aligned_ref, aligned_target):
            if rb != "-" and tb != "-":
                if rb == tb:
                    matches += 1
                q_used += 1
                t_used += 1
                ref_pos += 1
                del_run = 0
                continue
            if rb == "-" and tb != "-":
                offset += 1
                t_used += 1
                changes.append((ref_pos, offset))
                del_run = 0
                continue
            if rb != "-" and tb == "-":
                q_used += 1
                ref_pos += 1
                offset -= 1
                changes.append((ref_pos, offset))
                del_run += 1
                max_q_gap = max(max_q_gap, del_run)

        if q_used <= 0:
            continue

        identity = matches / max(1, q_used)
        block = SyntenicBlock(
            ref_interval=GenomicInterval(
                seq_region=ref_chr,
                start=ref_window_start,
                end=ref_window_end,
                strand=Strand.PLUS,
            ),
            target_interval=GenomicInterval(
                seq_region=target_chr,
                start=target_window_start,
                end=target_window_end,
                strand=Strand.MINUS if is_rev else Strand.PLUS,
            ),
            identity=identity,
            alignment_length=max(q_used, t_used),
            matches=matches,
            mismatches=max(0, max(q_used, t_used) - matches),
            block_id="internal_pairwise",
        )
        block.cached_offset_map = _compress_offset_map(changes)

        candidate = InternalRealignResult(
            blocks=[block],
            backend="pairwise",
            path_coverage=q_used / max(1, len(ref_seq)),
            internal_gap_bp=max_q_gap,
        )

        score = (identity, aln.score)
        if best is None or score > best[0]:
            best = (score, candidate)

    return None if best is None else best[1]


def targeted_internal_realign(
    ref_seq: str,
    target_seq: str,
    ref_chr: str,
    target_chr: str,
    gene_start: int,
    gene_end: int,
    ref_window_start: int,
    ref_window_end: int,
    target_window_start: int,
    target_window_end: int,
    threads: int,
    min_path_coverage: float,
    max_internal_gap_bp: int,
    fallback_backend: str = "edlib",
    temp_dir: Optional[Path] = None,
) -> Optional[InternalRealignResult]:
    """Attempt targeted minimap2+cs first, then fallback backend."""
    mm2_records = _run_targeted_minimap2(
        ref_seq=ref_seq,
        target_seq=target_seq,
        threads=threads,
        temp_dir=temp_dir,
    )

    if mm2_records:
        mm2_path = _select_minimap2_path(
            mm2_records,
            ref_chr=ref_chr,
            target_chr=target_chr,
            ref_window_start=ref_window_start,
            target_window_start=target_window_start,
            gene_start=gene_start,
            gene_end=gene_end,
            min_path_coverage=min_path_coverage,
            max_internal_gap_bp=max_internal_gap_bp,
        )
        if mm2_path is not None and mm2_path.blocks:
            return mm2_path

    if fallback_backend == "edlib":
        edlib_result = _edlib_block(
            ref_seq=ref_seq,
            target_seq=target_seq,
            ref_chr=ref_chr,
            target_chr=target_chr,
            ref_window_start=ref_window_start,
            ref_window_end=ref_window_end,
            target_window_start=target_window_start,
            target_window_end=target_window_end,
        )
        if edlib_result is not None:
            return edlib_result

    # Keep pipeline functional even when edlib is unavailable.
    fallback = _pairwise_block(
        ref_seq=ref_seq,
        target_seq=target_seq,
        ref_chr=ref_chr,
        target_chr=target_chr,
        ref_window_start=ref_window_start,
        ref_window_end=ref_window_end,
        target_window_start=target_window_start,
        target_window_end=target_window_end,
    )
    return fallback

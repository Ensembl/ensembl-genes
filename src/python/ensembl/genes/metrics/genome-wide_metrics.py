#!/usr/bin/env python3

import argparse
from collections import Counter
from dataclasses import dataclass
import gzip
import heapq
import re
import sys


LOWERCASE_RUN = re.compile(rb"[a-z]+")
OPPOSITE_STRAND = {"+": "-", "-": "+"}


@dataclass(frozen=True)
class FastaMetrics:
    genome_size: int
    total_masked: int = 0
    masked_annotated: int = 0


def open_maybe_gzip(path, mode):
    """Open plain or gzip-compressed files using the requested mode."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_gff(gff_path, feature_type):
    """Parse GFF/GTF and return feature coordinates as 0-based half-open intervals."""
    features = []

    with open_maybe_gzip(gff_path, "rt") as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue

            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 7 or cols[2] != feature_type:
                continue

            start = int(cols[3]) - 1
            end = int(cols[4])
            if end <= start:
                continue
            features.append((cols[0], start, end, cols[6]))

    return features


def genome_size_from_fasta(fasta_path):
    """Calculate genome size from a FASTA file, including soft-masked bases."""
    return scan_fasta(fasta_path).genome_size


def merge_feature_intervals(features_sorted):
    """Merge feature intervals by chromosome and return merged intervals plus covered bp."""
    merged_by_chrom = {}
    covered_bp = 0
    current_chrom = None
    current_start = None
    current_end = None

    for chrom, start, end, _strand in features_sorted:
        if chrom != current_chrom:
            if current_chrom is not None:
                merged_by_chrom.setdefault(current_chrom, []).append((current_start, current_end))
                covered_bp += current_end - current_start
            current_chrom = chrom
            current_start = start
            current_end = end
        elif start <= current_end:
            if end > current_end:
                current_end = end
        else:
            merged_by_chrom.setdefault(current_chrom, []).append((current_start, current_end))
            covered_bp += current_end - current_start
            current_start = start
            current_end = end

    if current_chrom is not None:
        merged_by_chrom.setdefault(current_chrom, []).append((current_start, current_end))
        covered_bp += current_end - current_start

    return merged_by_chrom, covered_bp


def overlap_length_from_cursor(intervals, cursor, start, end):
    """Return bp overlap with sorted, merged intervals and an advanced cursor."""
    while cursor < len(intervals) and intervals[cursor][1] <= start:
        cursor += 1

    total = 0
    idx = cursor
    while idx < len(intervals) and intervals[idx][0] < end:
        interval_start, interval_end = intervals[idx]
        overlap_start = max(start, interval_start)
        overlap_end = min(end, interval_end)
        if overlap_end > overlap_start:
            total += overlap_end - overlap_start
        idx += 1

    while cursor < len(intervals) and intervals[cursor][1] <= end:
        cursor += 1

    return total, cursor


def scan_fasta(fasta_path, merged_features_by_chrom=None):
    """Scan FASTA once for genome size and, when intervals are supplied, masking metrics.

    masked_annotated intentionally matches the original bedtools command:
    each full soft-masked run is counted if it overlaps any feature.
    """
    genome_size = 0
    total_masked = 0
    masked_annotated = 0
    calculate_masking = merged_features_by_chrom is not None

    chrom = None
    pos = 0
    in_mask = False
    mask_start = 0
    mask_end = 0
    interval_cursor = 0

    def finish_mask(mask_chrom, start, end):
        nonlocal total_masked, masked_annotated, interval_cursor

        if end <= start:
            return

        mask_len = end - start
        total_masked += mask_len

        intervals = merged_features_by_chrom.get(mask_chrom, ())
        if not intervals:
            return

        overlap_bp, interval_cursor = overlap_length_from_cursor(intervals, interval_cursor, start, end)
        if overlap_bp:
            masked_annotated += mask_len

    with open_maybe_gzip(fasta_path, "rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(b">"):
                if calculate_masking and in_mask:
                    finish_mask(chrom, mask_start, mask_end)
                    in_mask = False

                chrom = raw_line[1:].split(None, 1)[0].decode()
                pos = 0
                interval_cursor = 0
                continue

            seq = raw_line.strip()
            if not seq:
                continue

            line_start = pos
            line_end = line_start + len(seq)
            genome_size += len(seq)

            if calculate_masking:
                saw_mask_run = False
                for match in LOWERCASE_RUN.finditer(seq):
                    saw_mask_run = True
                    run_start = line_start + match.start()
                    run_end = line_start + match.end()

                    if in_mask and run_start == mask_end:
                        mask_end = run_end
                    else:
                        if in_mask:
                            finish_mask(chrom, mask_start, mask_end)
                        mask_start = run_start
                        mask_end = run_end
                        in_mask = True

                if not saw_mask_run and in_mask:
                    finish_mask(chrom, mask_start, mask_end)
                    in_mask = False
                elif in_mask and mask_end < line_end:
                    finish_mask(chrom, mask_start, mask_end)
                    in_mask = False

            pos = line_end

    if calculate_masking and in_mask:
        finish_mask(chrom, mask_start, mask_end)

    return FastaMetrics(genome_size, total_masked, masked_annotated)


def gene_density(num_features, genome_size):
    """Calculate feature density as features per megabase."""
    return (num_features / genome_size) * 1e6


def pct_genome_covered(covered_bp, genome_size):
    """Calculate percentage of the genome covered by merged features."""
    return (covered_bp / genome_size) * 100


def pct_masked_annotated(masked_annotated, total_masked):
    """Calculate percentage of soft-masked bases in runs overlapping annotations."""
    if total_masked == 0:
        return None
    return (masked_annotated / total_masked) * 100


def count_stranded_overlap_pairs(features_sorted):
    """Count same-strand and antisense overlapping feature pairs in one sweep."""
    same_strand = 0
    antisense = 0
    current_chrom = None
    active_by_end = []
    active_strands = Counter()

    for chrom, start, end, strand in features_sorted:
        if chrom != current_chrom:
            current_chrom = chrom
            active_by_end = []
            active_strands.clear()

        while active_by_end and active_by_end[0][0] <= start:
            _active_end, active_strand = heapq.heappop(active_by_end)
            active_strands[active_strand] -= 1

        same_strand += active_strands[strand]
        opposite = OPPOSITE_STRAND.get(strand)
        if opposite:
            antisense += active_strands[opposite]

        heapq.heappush(active_by_end, (end, strand))
        active_strands[strand] += 1

    return same_strand, antisense


def positive_int(value):
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def main():
    parser = argparse.ArgumentParser(description="Calculate genome-wide metrics for Ensembl features.")
    parser.add_argument("gff", help="Path to the input GFF/GTF file containing annotations")
    parser.add_argument(
        "--genome-size",
        type=positive_int,
        help="Size of the genome in base pairs. Not needed if --fasta is provided",
    )
    parser.add_argument("--fasta", help="Path to the soft-masked genome FASTA file.")
    parser.add_argument(
        "--feature-type",
        default="gene",
        help="Feature type to consider for counting (default: gene)",
    )
    args = parser.parse_args()

    if args.genome_size is None and args.fasta is None:
        parser.error("You must provide either --genome-size or --fasta to determine the genome size.")

    features = parse_gff(args.gff, args.feature_type)
    if not features:
        print(f"[ERROR] No {args.feature_type} features found in GFF/GTF.", file=sys.stderr)
        sys.exit(1)

    features_sorted = sorted(features, key=lambda item: (item[0], item[1], item[2], item[3]))
    merged_features_by_chrom, covered_bp = merge_feature_intervals(features_sorted)
    same_strand_pairs, antisense_pairs = count_stranded_overlap_pairs(features_sorted)

    fasta_metrics = None
    if args.fasta:
        fasta_metrics = scan_fasta(args.fasta, merged_features_by_chrom)
        genome_size = args.genome_size if args.genome_size is not None else fasta_metrics.genome_size
    else:
        genome_size = args.genome_size

    print(f"{'Metric':<45} {'Value'}")
    print("-" * 65)
    print(f"{'Genome size (bp)':<45} {genome_size}")
    print(f"{f'{args.feature_type} count':<45} {len(features)}")
    print(f"{f'{args.feature_type} density ({args.feature_type}s/Mb)':<45} {gene_density(len(features), genome_size):.4f}")
    print(f"{'% genome covered by ' + args.feature_type:<45} {pct_genome_covered(covered_bp, genome_size):.4f}")

    if fasta_metrics:
        if fasta_metrics.total_masked == 0:
            print(f"{'% masked regions annotated':<45} N/A (no soft-masking detected)")
        else:
            masked_pct = pct_masked_annotated(fasta_metrics.masked_annotated, fasta_metrics.total_masked)
            masked_genome_pct = (fasta_metrics.total_masked / genome_size) * 100
            print(f"{'Total soft-masked bases':<45} {fasta_metrics.total_masked} ({masked_genome_pct:.2f}%)")
            print(f"{'% masked regions annotated':<45} {masked_pct:.4f}")
    else:
        print(f"{'% masked regions annotated':<45} skipped (no --fasta provided)")

    print(f"{f'Overlapping {args.feature_type} pairs (same strand)':<45} {same_strand_pairs}")
    print(f"{f'Antisense {args.feature_type} pairs (opposite strand)':<45} {antisense_pairs}")


if __name__ == "__main__":
    main()


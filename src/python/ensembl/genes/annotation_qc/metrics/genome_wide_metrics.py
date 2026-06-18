#!/usr/bin/env python3
"""Calculate genome-wide annotation QC metrics from GFF/GTF and FASTA files."""

# pulint: disable=line-too-long,missing-module-docstring,missing-class-docstring
# pylint: disable=too-many-arguments,too-many-branches,too-many-statements
# pylint: disable=too-many-locals,missing-class-docstring,missing-function-docstring

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
import gzip
import heapq
import os
import re
import sys


LOWERCASE_RUN = re.compile(rb"[a-z]+")
FILENAME_TAG_RE = re.compile(r"[^A-Za-z0-9_.-]+")
GTF_GFF_COLUMN_COUNT = 9
OPPOSITE_STRAND = {"+": "-", "-": "+"}
CANONICAL_TAGS = {"ensembl_canonical", "canonical", "mane_select"}
TRUE_VALUES = {"1", "true", "yes"}
CSV_HEADERS = ("metrics_name", "metrics_value")
ANNOTATION_SUFFIXES = (".gff3", ".gff", ".gtf")
BIOTYPE_KEYS = (
    "biotype",
    "gene_biotype",
    "transcript_biotype",
    "gene_type",
    "transcript_type",
)


@dataclass(frozen=True)
class FastaMetrics:
    genome_size: int
    total_masked: int = 0
    masked_annotated: int = 0


@dataclass(frozen=True)
class Feature:
    chrom: str
    start: int
    end: int
    strand: str
    feature_id: str
    gene_id: str = ""
    transcript_id: str = ""


def open_maybe_gzip(path, mode, encoding="utf-8", newline=None):
    """Open plain or gzip-compressed files using the requested mode."""
    if str(path).endswith(".gz"):
        if "b" in mode:
            return gzip.open(path, mode)  # pylint: disable=unspecified-encoding
        return gzip.open(path, mode, encoding=encoding, newline=newline)

    if "b" in mode:
        return open(path, mode)  # pylint: disable=unspecified-encoding
    return open(path, mode, encoding=encoding, newline=newline)


def warn(message):
    """Print a non-fatal warning to stderr."""
    print(f"[WARNING] {message}", file=sys.stderr)


def sanitize_filename_tag(value, default="all"):
    """Return a filesystem-safe filename tag."""
    tag = FILENAME_TAG_RE.sub("_", str(value)).strip("._-")
    return tag or default


def annotation_base_name(annotation_path):
    """Return annotation filename without GFF/GTF and optional gzip suffix."""
    name = os.path.basename(str(annotation_path))
    if name.lower().endswith(".gz"):
        name = name[:-3]

    lower_name = name.lower()
    for suffix in ANNOTATION_SUFFIXES:
        if lower_name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    return sanitize_filename_tag(name, default="annotation")


def biotype_filename_tag(biotypes):
    """Return the filename tag for the selected biotypes."""
    if not biotypes:
        return "all"

    return "-".join(
        sanitize_filename_tag(biotype, default="biotype")
        for biotype in sorted(biotypes)
    )


def metrics_csv_filename(annotation_path, biotypes, canonical_only):
    """Build the automatic CSV metrics output filename."""
    canonical_tag = "canonical" if canonical_only else "all"
    filename_parts = (
        annotation_base_name(annotation_path),
        biotype_filename_tag(biotypes),
        canonical_tag,
        "genome",
        "csv",
    )
    return ".".join(filename_parts)


def add_metric(metric_rows, name, value):
    """Append a metric row using the terminal display value."""
    metric_rows.append((name, str(value)))


def print_metric_rows(metric_rows):
    """Print collected metric rows to stdout."""
    print(f"{'Metric':<45} {'Value'}")
    print("-" * 65)
    for metric_name, metric_value in metric_rows:
        print(f"{metric_name:<45} {metric_value}")


def write_metric_rows_csv(metric_rows, csv_path):
    """Write collected metric rows to CSV."""
    with open(csv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADERS)
        writer.writerows(metric_rows)


def split_annotation_line(line):
    """Split a GFF/GTF row while preserving spaces inside the attributes column."""
    stripped = line.rstrip("\r\n")
    if "\t" in stripped:
        return stripped.split("\t")
    return stripped.split(None, GTF_GFF_COLUMN_COUNT - 1)


def parse_attributes(attr_text):
    """Parse GTF/GFF3 attributes into a small dictionary."""
    attrs = {}

    for raw_field in attr_text.strip().rstrip(";").split(";"):
        field = raw_field.strip()
        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            key = key.strip()
            value = value.strip().strip('"')
            values = [
                part.strip().strip('"') for part in value.split(",") if part.strip()
            ]
        else:
            parts = field.split(None, 1)
            key = parts[0]
            value = parts[1].strip().strip('"') if len(parts) > 1 else ""
            values = [value]

        if key.lower() == "tag":
            attrs.setdefault("tag", []).extend(values)
        elif len(values) == 1:
            attrs[key] = values[0]
        else:
            attrs[key] = values

    return attrs


def attr_first(attrs, keys):
    """Return the first non-empty attribute value for any key in keys."""
    for key in keys:
        value = attrs.get(key)
        if isinstance(value, list):
            if value:
                return value[0]
        elif value:
            return value
    return ""


def is_canonical_attributes(attrs):
    """Detect common Ensembl GTF canonical transcript annotations."""
    tags = attrs.get("tag", [])
    if isinstance(tags, str):
        tags = [tags]
    if any(tag.strip().strip('"').lower() in CANONICAL_TAGS for tag in tags):
        return True

    for key in ("transcript_is_canonical", "is_canonical"):
        value = attrs.get(key, "")
        if isinstance(value, str) and value.strip().strip('"').lower() in TRUE_VALUES:
            return True

    return False


def feature_id_from_attrs(attrs, feature_type, fallback):
    """Choose a useful report ID from GTF/GFF3 attributes."""
    if feature_type == "gene":
        keys = ("gene_id", "ID", "Name", "gene_name")
    elif feature_type == "transcript":
        keys = ("transcript_id", "ID", "Name")
    else:
        keys = ("ID", "transcript_id", "gene_id", "Name")

    return attr_first(attrs, keys) or fallback


def parse_biotype_filter(values):
    """Normalize one or more comma-separated --biotype values."""
    if not values:
        return None

    biotypes = set()
    for value in values:
        for biotype in value.split(","):
            biotype = biotype.strip()
            if biotype:
                biotypes.add(biotype.lower())

    return biotypes or None


def feature_matches_biotype(attrs, allowed_biotypes):
    """Return True when any known biotype attribute matches the requested filter."""
    if not allowed_biotypes:
        return True

    for key in BIOTYPE_KEYS:
        value = attrs.get(key)
        if isinstance(value, list):
            if any(item.lower() in allowed_biotypes for item in value):
                return True
        elif isinstance(value, str) and value.lower() in allowed_biotypes:
            return True

    return False


def collect_canonical_ids(gtf_path):
    """Collect transcript/gene IDs marked as canonical in a GTF-like annotation."""
    transcript_ids = set()
    gene_ids = set()

    with open_maybe_gzip(gtf_path, "rt") as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue

            cols = split_annotation_line(line)
            if len(cols) < 9:
                continue

            attrs = parse_attributes(cols[8])
            if not is_canonical_attributes(attrs):
                continue

            transcript_id = attr_first(attrs, ("transcript_id", "ID"))
            gene_id = attr_first(attrs, ("gene_id", "Parent"))
            if transcript_id:
                transcript_ids.add(transcript_id)
            if gene_id:
                gene_ids.add(gene_id)

    return transcript_ids, gene_ids


def parse_gff(
    gff_path,
    feature_type,
    canonical_only=False,
    canonical_transcript_ids=None,
    canonical_gene_ids=None,
    biotypes=None,
):
    """Parse GFF/GTF and return feature coordinates as 0-based half-open intervals."""
    features = []
    canonical_transcript_ids = canonical_transcript_ids or set()
    canonical_gene_ids = canonical_gene_ids or set()

    with open_maybe_gzip(gff_path, "rt") as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue

            cols = split_annotation_line(line)
            if len(cols) < 7 or cols[2] != feature_type:
                continue

            attrs = parse_attributes(cols[8]) if len(cols) >= 9 else {}
            transcript_id = attr_first(attrs, ("transcript_id", "ID"))
            gene_id = attr_first(attrs, ("gene_id", "Parent"))

            if biotypes and not feature_matches_biotype(attrs, biotypes):
                continue

            if canonical_only:
                feature_is_canonical = is_canonical_attributes(attrs)
                if transcript_id:
                    keep_feature = (
                        feature_is_canonical
                        or transcript_id in canonical_transcript_ids
                    )
                elif feature_type == "gene" and gene_id:
                    keep_feature = feature_is_canonical or gene_id in canonical_gene_ids
                else:
                    keep_feature = feature_is_canonical
                if not keep_feature:
                    continue

            start = int(cols[3]) - 1
            end = int(cols[4])
            if end <= start:
                continue

            feature_id = feature_id_from_attrs(
                attrs, feature_type, f"{feature_type}_{len(features) + 1}"
            )
            features.append(
                Feature(
                    cols[0], start, end, cols[6], feature_id, gene_id, transcript_id
                )
            )

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

    for feature in features_sorted:
        if feature.chrom != current_chrom:
            if current_chrom is not None:
                merged_by_chrom.setdefault(current_chrom, []).append(
                    (current_start, current_end)
                )
                covered_bp += current_end - current_start
            current_chrom = feature.chrom
            current_start = feature.start
            current_end = feature.end
        elif feature.start <= current_end:
            current_end = max(current_end, feature.end)
        else:
            merged_by_chrom.setdefault(current_chrom, []).append(
                (current_start, current_end)
            )
            covered_bp += current_end - current_start
            current_start = feature.start
            current_end = feature.end

    if current_chrom is not None:
        merged_by_chrom.setdefault(current_chrom, []).append(
            (current_start, current_end)
        )
        covered_bp += current_end - current_start

    return merged_by_chrom, covered_bp


def group_features_by_chrom(features_sorted):
    """Group sorted features by chromosome for overlap reporting."""
    features_by_chrom = {}
    for feature in features_sorted:
        features_by_chrom.setdefault(feature.chrom, []).append(feature)
    return features_by_chrom


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


def scan_fasta(
    fasta_path,
    merged_features_by_chrom=None,
    masked_overlap_output=None,
    report_features_by_chrom=None,
):
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
    report_interval_cursor = 0
    overlap_out = None

    def finish_mask(mask_chrom, start, end):
        nonlocal total_masked, masked_annotated, interval_cursor

        if end <= start:
            return

        mask_len = end - start
        total_masked += mask_len

        intervals = merged_features_by_chrom.get(mask_chrom, ())
        if not intervals:
            return

        overlap_bp, interval_cursor = overlap_length_from_cursor(
            intervals, interval_cursor, start, end
        )
        if overlap_bp:
            masked_annotated += mask_len

    def write_mask_overlaps(mask_chrom, start, end):
        nonlocal report_interval_cursor

        if overlap_out is None or report_features_by_chrom is None:
            return

        features = report_features_by_chrom.get(mask_chrom, ())
        while (
            report_interval_cursor < len(features)
            and features[report_interval_cursor].end <= start
        ):
            report_interval_cursor += 1

        idx = report_interval_cursor
        while idx < len(features) and features[idx].start < end:
            feature = features[idx]
            overlap_start = max(start, feature.start)
            overlap_end = min(end, feature.end)
            if overlap_end > overlap_start:
                overlap_out.write(
                    f"{mask_chrom}\t{overlap_start}\t{overlap_end}\t"
                    f"{start}\t{end}\t{feature.feature_id}\t"
                    f"{feature.start}\t{feature.end}\t{feature.strand}\t"
                    f"{overlap_end - overlap_start}\n"
                )
            idx += 1

        while (
            report_interval_cursor < len(features)
            and features[report_interval_cursor].end <= end
        ):
            report_interval_cursor += 1

    if masked_overlap_output:
        overlap_out = open_maybe_gzip(masked_overlap_output, "wt")
        overlap_out.write(
            "#chrom\toverlap_start\toverlap_end\tmask_start\tmask_end\t"
            "feature_id\tfeature_start\tfeature_end\tfeature_strand\toverlap_bp\n"
        )

    try:
        with open_maybe_gzip(fasta_path, "rb") as handle:
            for raw_line in handle:
                if raw_line.startswith(b">"):
                    if calculate_masking and in_mask:
                        finish_mask(chrom, mask_start, mask_end)
                        write_mask_overlaps(chrom, mask_start, mask_end)
                        in_mask = False

                    chrom = raw_line[1:].split(None, 1)[0].decode()
                    pos = 0
                    interval_cursor = 0
                    report_interval_cursor = 0
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
                                write_mask_overlaps(chrom, mask_start, mask_end)
                            mask_start = run_start
                            mask_end = run_end
                            in_mask = True

                    if not saw_mask_run and in_mask:
                        finish_mask(chrom, mask_start, mask_end)
                        write_mask_overlaps(chrom, mask_start, mask_end)
                        in_mask = False
                    elif in_mask and mask_end < line_end:
                        finish_mask(chrom, mask_start, mask_end)
                        write_mask_overlaps(chrom, mask_start, mask_end)
                        in_mask = False

                pos = line_end

        if calculate_masking and in_mask:
            finish_mask(chrom, mask_start, mask_end)
            write_mask_overlaps(chrom, mask_start, mask_end)
    finally:
        if overlap_out is not None:
            overlap_out.close()

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


def count_stranded_overlap_pairs(features_sorted, feature_overlap_output=None):
    """Count same-strand and antisense overlapping feature pairs in one sweep."""
    if feature_overlap_output:
        return count_and_write_feature_overlap_pairs(
            features_sorted, feature_overlap_output
        )

    same_strand = 0
    antisense = 0
    current_chrom = None
    active_by_end = []
    active_strands = Counter()

    for feature in features_sorted:
        if feature.chrom != current_chrom:
            current_chrom = feature.chrom
            active_by_end = []
            active_strands.clear()

        while active_by_end and active_by_end[0][0] <= feature.start:
            _active_end, active_strand = heapq.heappop(active_by_end)
            active_strands[active_strand] -= 1

        same_strand += active_strands[feature.strand]
        opposite = OPPOSITE_STRAND.get(feature.strand)
        if opposite:
            antisense += active_strands[opposite]

        heapq.heappush(active_by_end, (feature.end, feature.strand))
        active_strands[feature.strand] += 1

    return same_strand, antisense


def count_and_write_feature_overlap_pairs(features_sorted, feature_overlap_output):
    """Count overlaps and write a BED-like TSV describing every overlapping pair."""
    same_strand = 0
    antisense = 0
    current_chrom = None
    active_by_end = []
    active_records = []
    active_ids = set()
    next_active_id = 0

    with open_maybe_gzip(feature_overlap_output, "wt") as out:
        out.write(
            "#chrom\toverlap_start\toverlap_end\tfeature1_id\tfeature1_start\tfeature1_end\t"
            "feature1_strand\tfeature2_id\tfeature2_start\tfeature2_end\tfeature2_strand\t"
            "relation\toverlap_bp\n"
        )

        for feature in features_sorted:
            if feature.chrom != current_chrom:
                current_chrom = feature.chrom
                active_by_end = []
                active_records = []
                active_ids = set()

            while active_by_end and active_by_end[0][0] <= feature.start:
                _active_end, active_id = heapq.heappop(active_by_end)
                active_ids.discard(active_id)

            for active_id, other in active_records:
                if active_id not in active_ids:
                    continue

                overlap_start = max(other.start, feature.start)
                overlap_end = min(other.end, feature.end)
                if overlap_end <= overlap_start:
                    continue
                overlap_bp = overlap_end - overlap_start

                if other.strand == feature.strand:
                    relation = "same_strand"
                    same_strand += 1
                elif OPPOSITE_STRAND.get(feature.strand) == other.strand:
                    relation = "antisense"
                    antisense += 1
                else:
                    relation = "other_strand"

                out.write(
                    f"{feature.chrom}\t{overlap_start}\t{overlap_end}\t"
                    f"{other.feature_id}\t{other.start}\t{other.end}\t{other.strand}\t"
                    f"{feature.feature_id}\t{feature.start}\t{feature.end}\t{feature.strand}\t"
                    f"{relation}\t{overlap_bp}\n"
                )

            active_id = next_active_id
            next_active_id += 1
            active_ids.add(active_id)
            active_records.append((active_id, feature))
            heapq.heappush(active_by_end, (feature.end, active_id))

            if len(active_records) > (len(active_ids) * 2) + 1000:
                active_records = [
                    (active_id, item)
                    for active_id, item in active_records
                    if active_id in active_ids
                ]

    return same_strand, antisense


def positive_int(value):
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def main():
    """Parse command-line arguments and print genome-wide QC metrics."""
    parser = argparse.ArgumentParser(
        description="Calculate genome-wide metrics for Ensembl features."
    )
    parser.add_argument(
        "gff", help="Path to the input GFF/GTF file containing annotations"
    )
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
    parser.add_argument(
        "--canonical-only",
        action="store_true",
        help="For GTF input, keep only features belonging to canonical transcripts.",
    )
    parser.add_argument(
        "--biotype",
        action="append",
        help="Keep only features with this biotype. Can be repeated or comma-separated.",
    )
    parser.add_argument(
        "--masked-overlaps-output",
        help="Write a TSV of soft-masked/feature overlap coordinates. Requires --fasta.",
    )
    parser.add_argument(
        "--feature-overlaps-output",
        help="Write a TSV of overlapping feature pairs and their overlap coordinates.",
    )
    parser.add_argument(
        "--csv-output",
        action="store_true",
        help="Write metric rows to an automatically named CSV file.",
    )
    args = parser.parse_args()

    if args.genome_size is None and args.fasta is None:
        parser.error(
            "You must provide either --genome-size or --fasta to determine the genome size."
        )
    if args.masked_overlaps_output and args.fasta is None:
        parser.error("--masked-overlaps-output requires --fasta.")

    biotypes = parse_biotype_filter(args.biotype)

    canonical_transcript_ids = set()
    canonical_gene_ids = set()
    canonical_filter_applied = args.canonical_only
    if args.canonical_only:
        canonical_transcript_ids, canonical_gene_ids = collect_canonical_ids(args.gff)
        if not canonical_transcript_ids and not canonical_gene_ids:
            warn(
                "--canonical-only was requested, but no canonical annotations were detected; using all features."
            )
            canonical_filter_applied = False

    features = parse_gff(
        args.gff,
        args.feature_type,
        canonical_only=canonical_filter_applied,
        canonical_transcript_ids=canonical_transcript_ids,
        canonical_gene_ids=canonical_gene_ids,
        biotypes=biotypes,
    )
    if canonical_filter_applied and not features:
        warn(
            f"--canonical-only removed all {args.feature_type} features; "
            f"using all {args.feature_type} features instead."
        )
        canonical_filter_applied = False
        features = parse_gff(args.gff, args.feature_type, biotypes=biotypes)
    if not features:
        print(
            f"[ERROR] No {args.feature_type} features found in GFF/GTF.",
            file=sys.stderr,
        )
        sys.exit(1)

    features_sorted = sorted(
        features,
        key=lambda item: (
            item.chrom,
            item.start,
            item.end,
            item.strand,
            item.feature_id,
        ),
    )
    merged_features_by_chrom, covered_bp = merge_feature_intervals(features_sorted)
    same_strand_pairs, antisense_pairs = count_stranded_overlap_pairs(
        features_sorted,
        feature_overlap_output=args.feature_overlaps_output,
    )

    fasta_metrics = None
    if args.fasta:
        report_features_by_chrom = (
            group_features_by_chrom(features_sorted)
            if args.masked_overlaps_output
            else None
        )
        fasta_metrics = scan_fasta(
            args.fasta,
            merged_features_by_chrom,
            masked_overlap_output=args.masked_overlaps_output,
            report_features_by_chrom=report_features_by_chrom,
        )
        genome_size = (
            args.genome_size
            if args.genome_size is not None
            else fasta_metrics.genome_size
        )
    else:
        genome_size = args.genome_size

    metric_rows = []
    add_metric(metric_rows, "Genome size (bp)", genome_size)
    add_metric(metric_rows, f"{args.feature_type} count", len(features))
    if args.canonical_only:
        filter_status = (
            "applied"
            if canonical_filter_applied
            else "not applied (using all features)"
        )
        add_metric(metric_rows, "Canonical-only filter", filter_status)
    if biotypes:
        add_metric(metric_rows, "Biotype filter", ", ".join(sorted(biotypes)))

    density_name = f"{args.feature_type} density ({args.feature_type}s/Mb)"
    add_metric(
        metric_rows,
        density_name,
        f"{gene_density(len(features), genome_size):.4f}",
    )
    add_metric(
        metric_rows,
        f"% genome covered by {args.feature_type}",
        f"{pct_genome_covered(covered_bp, genome_size):.4f}",
    )

    if fasta_metrics:
        if fasta_metrics.total_masked == 0:
            add_metric(
                metric_rows,
                "% masked regions annotated",
                "N/A (no soft-masking detected)",
            )
        else:
            masked_pct = pct_masked_annotated(
                fasta_metrics.masked_annotated, fasta_metrics.total_masked
            )
            masked_genome_pct = (fasta_metrics.total_masked / genome_size) * 100
            add_metric(
                metric_rows,
                "Total soft-masked bases",
                f"{fasta_metrics.total_masked} ({masked_genome_pct:.2f}%)",
            )
            add_metric(metric_rows, "% masked regions annotated", f"{masked_pct:.4f}")
    else:
        add_metric(
            metric_rows,
            "% masked regions annotated",
            "skipped (no --fasta provided)",
        )

    add_metric(
        metric_rows,
        f"Overlapping {args.feature_type} pairs (same strand)",
        same_strand_pairs,
    )
    add_metric(
        metric_rows,
        f"Antisense {args.feature_type} pairs (opposite strand)",
        antisense_pairs,
    )

    print_metric_rows(metric_rows)
    if args.csv_output:
        write_metric_rows_csv(
            metric_rows,
            metrics_csv_filename(args.gff, biotypes, args.canonical_only),
        )


if __name__ == "__main__":
    main()

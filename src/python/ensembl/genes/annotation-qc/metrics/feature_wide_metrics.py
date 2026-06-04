#!/usr/bin/env python3
"""Calculate feature-level QC metrics for GFF3/GTF genome annotations."""

# pylint: disable=missing-function-docstring,missing-class-docstring

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import re
import sys
from urllib.parse import unquote


VALID_PHASES = {"0", "1", "2"}
PHASE_LABELS = {"0": "3n", "1": "3n+1", "2": "3n+2"}
GTF_ATTRIBUTE = re.compile(r'(\S+)\s+"([^"]*)"')
UTR_LABELS = {"five_prime": "5' UTR", "three_prime": "3' UTR"}
BIOTYPE_ATTRIBUTE_KEYS = (
    "transcript_biotype",
    "transcript_type",
    "biotype",
    "gene_biotype",
    "gene_type",
)
CANONICAL_ATTRIBUTE_KEYS = (
    "canonical",
    "is_canonical",
    "transcript_is_canonical",
)
CANONICAL_TAGS = {"canonical", "ensembl_canonical", "mane_select"}
TRANSCRIPT_FEATURE_TYPES = {
    "transcript",
    "mrna",
    "rrna",
    "trna",
    "ncrna",
    "lncrna",
    "snorna",
    "snrna",
}
FIVE_PRIME_UTR_TYPES = {
    "five_prime_utr",
    "five_prime_untranslated_region",
    "5_prime_utr",
    "5_utr",
    "5utr",
}
THREE_PRIME_UTR_TYPES = {
    "three_prime_utr",
    "three_prime_untranslated_region",
    "3_prime_utr",
    "3_utr",
    "3utr",
}
INTRON_LENGTH_BINS = (
    (1, 19),
    (20, 99),
    (100, 499),
    (500, 999),
    (1000, 4999),
    (5000, 9999),
    (10000, 49999),
    (50000, None),
)
EXON_LENGTH_BINS = (
    (1, 9),
    (10, 24),
    (25, 49),
    (50, 99),
    (100, 199),
    (200, 499),
    (500, 999),
    (1000, None),
)


@dataclass(frozen=True)
class Exon:
    chrom: str
    start: int
    end: int
    strand: str
    phase: str
    index: int


@dataclass(frozen=True)
class FeatureRecord:
    chrom: str
    start: int
    end: int
    strand: str
    phase: str
    transcript_ids: tuple
    utr_type: str | None = None


@dataclass(frozen=True)
class FlaggedFeature:
    flag_type: str
    transcript_id: str
    feature_type: str
    chrom: str
    start: int
    end: int
    strand: str
    length: int
    details: str


@dataclass(frozen=True)
class AnnotationMetrics:
    annotation_format: str
    exons: list
    exon_phase_counts: Counter
    exon_phase_missing: int
    cds_phase_counts: Counter
    cds_phase_missing: int
    exons_by_transcript: dict
    cds_lengths_by_transcript: dict
    utrs_by_transcript: dict
    cds_records_without_phase: list
    filters_applied: bool
    selected_transcript_count: int


@dataclass(frozen=True)
class NumericSummary:
    count: int
    mean: float
    median: float
    minimum: float
    maximum: float


def open_maybe_gzip(path, mode, encoding="utf-8"):
    """Open plain or gzip-compressed files using the requested mode."""
    if str(path).endswith(".gz"):
        if "b" in mode:
            return gzip.open(path, mode)
        return gzip.open(path, mode, encoding=encoding)

    if "b" in mode:
        return open(path, mode)
    return open(path, mode, encoding=encoding)


def positive_int(value):
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def gc_bin_size(value):
    parsed = positive_int(value)
    if parsed > 100:
        raise argparse.ArgumentTypeError("must be between 1 and 100")
    return parsed


def normalize_annotation_format(value):
    if value == "gff":
        return "gff3"
    return value


def detect_annotation_format(annotation_path):
    """Detect GFF3 vs GTF from header/attributes, falling back to filename."""
    with open_maybe_gzip(annotation_path, "rt") as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith("##gff-version"):
                return "gff3"
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue

            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                continue

            attributes = cols[8]
            if GTF_ATTRIBUTE.search(attributes):
                return "gtf"
            if "=" in attributes:
                return "gff3"

    lower_path = str(annotation_path).lower()
    if lower_path.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    return "gff3"


def parse_gff3_attributes(attributes):
    parsed = {}
    for item in attributes.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item, ""
        parsed[unquote(key)] = unquote(value)
    return parsed


def add_attribute(parsed, key, value):
    """Keep repeated GTF attributes such as tag by storing comma-separated values."""
    if key in parsed and parsed[key]:
        parsed[key] = f"{parsed[key]},{value}"
    else:
        parsed[key] = value


def parse_gtf_attributes(attributes):
    parsed = {}
    for key, value in GTF_ATTRIBUTE.findall(attributes):
        add_attribute(parsed, key, value)
    if parsed:
        return parsed

    for item in attributes.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if len(parts) == 2:
            add_attribute(parsed, parts[0], parts[1].strip('"'))
    return parsed


def parse_attributes(attributes, annotation_format):
    if annotation_format == "gtf":
        return parse_gtf_attributes(attributes)
    return parse_gff3_attributes(attributes)


def transcript_ids_from_attributes(attributes, annotation_format):
    """Return transcript IDs for transcript-level child features.

    GFF3 Parent may contain several comma-separated transcript IDs, so the same
    exon/CDS/UTR can contribute to multiple transcript-level calculations.
    """
    if annotation_format == "gtf":
        transcript_id = attributes.get("transcript_id")
        return [transcript_id] if transcript_id else []

    parents = attributes.get("Parent")
    if not parents:
        return []
    return [parent for parent in parents.split(",") if parent]


def transcript_id_from_transcript_feature(attributes, annotation_format):
    if annotation_format == "gtf":
        return attributes.get("transcript_id")
    return attributes.get("ID")


def split_attribute_values(value):
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def collect_biotypes(attributes):
    biotypes = set()
    for key in BIOTYPE_ATTRIBUTE_KEYS:
        for value in split_attribute_values(attributes.get(key)):
            biotypes.add(value.lower())
    return biotypes


def has_canonical_tag(attributes):
    for key in ("tag", "transcript_tag"):
        for value in split_attribute_values(attributes.get(key)):
            if value.lower() in CANONICAL_TAGS:
                return True

    for key in CANONICAL_ATTRIBUTE_KEYS:
        for value in split_attribute_values(attributes.get(key)):
            if value.lower() in {"1", "true", "yes", *CANONICAL_TAGS}:
                return True

    return False


def parse_filter_values(values):
    filters = set()
    for value in values or ():
        for item in value.split(","):
            item = item.strip()
            if item:
                filters.add(item.lower())
    return filters


def format_transcript_ids(transcript_ids):
    if not transcript_ids:
        return "."
    return ",".join(sorted(transcript_ids))


def utr_type_from_feature(feature_type):
    """Return normalized UTR type for common GFF3/GTF 5'/3' UTR feature names."""
    normalized = feature_type.lower().replace("-", "_").replace("'", "")
    if normalized in FIVE_PRIME_UTR_TYPES:
        return "five_prime"
    if normalized in THREE_PRIME_UTR_TYPES:
        return "three_prime"
    return None


def feature_passes_filters(transcript_ids, selected_transcripts):
    if selected_transcripts is None:
        return True
    return bool(set(transcript_ids) & selected_transcripts)


def filter_transcript_ids(transcript_ids, selected_transcripts):
    if selected_transcripts is None:
        return transcript_ids
    return tuple(
        transcript_id
        for transcript_id in transcript_ids
        if transcript_id in selected_transcripts
    )


def select_transcripts(
    all_transcript_ids,
    transcript_biotypes,
    canonical_transcripts,
    include_biotypes,
    exclude_biotypes,
    canonical_only,
):
    filters_applied = bool(include_biotypes or exclude_biotypes or canonical_only)
    if not filters_applied:
        return None

    selected = set()
    for transcript_id in all_transcript_ids:
        biotypes = transcript_biotypes.get(transcript_id, set())
        if include_biotypes and not biotypes & include_biotypes:
            continue
        if exclude_biotypes and (biotypes & exclude_biotypes):
            continue
        if canonical_only and transcript_id not in canonical_transcripts:
            continue
        selected.add(transcript_id)

    return selected


def parse_annotation(
    annotation_path,
    requested_format,
    include_biotypes=None,
    exclude_biotypes=None,
    canonical_only=False,
):
    """Parse annotation features needed for feature-level metrics."""
    annotation_format = (
        detect_annotation_format(annotation_path)
        if requested_format == "auto"
        else normalize_annotation_format(requested_format)
    )
    include_biotypes = include_biotypes or set()
    exclude_biotypes = exclude_biotypes or set()
    filters_applied = bool(include_biotypes or exclude_biotypes or canonical_only)
    all_transcript_ids = set()
    transcript_biotypes = defaultdict(set)
    canonical_transcripts = set()
    exon_records = []
    cds_records = []
    utr_records = []

    with open_maybe_gzip(annotation_path, "rt") as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue

            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 8:
                continue

            feature_type = cols[2]
            feature_type_lower = feature_type.lower()
            utr_type = utr_type_from_feature(feature_type)
            is_metric_feature = (
                feature_type_lower in {"exon", "cds"} or utr_type is not None
            )
            is_transcript_feature = feature_type_lower in TRANSCRIPT_FEATURE_TYPES
            if not is_metric_feature and not is_transcript_feature:
                continue

            attributes = parse_attributes(
                cols[8] if len(cols) > 8 else "", annotation_format
            )
            if is_transcript_feature:
                transcript_id = transcript_id_from_transcript_feature(
                    attributes, annotation_format
                )
                if transcript_id:
                    all_transcript_ids.add(transcript_id)
                    transcript_biotypes[transcript_id].update(
                        collect_biotypes(attributes)
                    )
                    if has_canonical_tag(attributes):
                        canonical_transcripts.add(transcript_id)

            if not is_metric_feature:
                continue

            try:
                start = int(cols[3]) - 1
                end = int(cols[4])
            except ValueError:
                continue
            if end <= start:
                continue

            transcript_ids = tuple(
                transcript_ids_from_attributes(attributes, annotation_format)
            )
            all_transcript_ids.update(transcript_ids)
            child_biotypes = collect_biotypes(attributes)
            child_is_canonical = has_canonical_tag(attributes)
            for transcript_id in transcript_ids:
                transcript_biotypes[transcript_id].update(child_biotypes)
                if child_is_canonical:
                    canonical_transcripts.add(transcript_id)

            phase = cols[7] if cols[7] in VALID_PHASES else None
            record = FeatureRecord(
                cols[0], start, end, cols[6], phase, transcript_ids, utr_type
            )

            if feature_type_lower == "exon":
                exon_records.append(record)
            elif feature_type_lower == "cds":
                cds_records.append(record)
            else:
                utr_records.append(record)

    selected_transcripts = select_transcripts(
        all_transcript_ids,
        transcript_biotypes,
        canonical_transcripts,
        include_biotypes,
        exclude_biotypes,
        canonical_only,
    )
    exons = []
    exon_phase_counts = Counter()
    cds_phase_counts = Counter()
    exon_phase_missing = 0
    cds_phase_missing = 0
    exons_by_transcript = defaultdict(list)
    cds_lengths_by_transcript = defaultdict(int)
    cds_records_without_phase = []
    utrs_by_transcript = {
        "five_prime": defaultdict(list),
        "three_prime": defaultdict(list),
    }

    for record in exon_records:
        if not feature_passes_filters(record.transcript_ids, selected_transcripts):
            continue

        exon = Exon(
            record.chrom,
            record.start,
            record.end,
            record.strand,
            record.phase,
            len(exons),
        )
        exons.append(exon)
        if record.phase is None:
            exon_phase_missing += 1
        else:
            exon_phase_counts[record.phase] += 1

        for transcript_id in filter_transcript_ids(
            record.transcript_ids, selected_transcripts
        ):
            exons_by_transcript[transcript_id].append(
                (record.chrom, record.start, record.end, record.strand)
            )

    for record in cds_records:
        if not feature_passes_filters(record.transcript_ids, selected_transcripts):
            continue

        filtered_transcript_ids = filter_transcript_ids(
            record.transcript_ids, selected_transcripts
        )
        if record.phase is None:
            cds_phase_missing += 1
            cds_records_without_phase.append(
                FlaggedFeature(
                    "cds_without_usable_phase",
                    format_transcript_ids(filtered_transcript_ids),
                    "CDS",
                    record.chrom,
                    record.start,
                    record.end,
                    record.strand,
                    record.end - record.start,
                    "phase is not 0, 1, or 2",
                )
            )
        else:
            cds_phase_counts[record.phase] += 1

        length = record.end - record.start
        for transcript_id in filtered_transcript_ids:
            cds_lengths_by_transcript[transcript_id] += length

    for record in utr_records:
        if not feature_passes_filters(record.transcript_ids, selected_transcripts):
            continue

        for transcript_id in filter_transcript_ids(
            record.transcript_ids, selected_transcripts
        ):
            utrs_by_transcript[record.utr_type][transcript_id].append(
                (record.chrom, record.start, record.end, record.strand)
            )

    return AnnotationMetrics(
        annotation_format,
        exons,
        exon_phase_counts,
        exon_phase_missing,
        cds_phase_counts,
        cds_phase_missing,
        dict(exons_by_transcript),
        dict(cds_lengths_by_transcript),
        {
            utr_type: dict(transcripts)
            for utr_type, transcripts in utrs_by_transcript.items()
        },
        cds_records_without_phase,
        filters_applied,
        len(selected_transcripts)
        if selected_transcripts is not None
        else len(all_transcript_ids),
    )


def summarize_numbers(values):
    if not values:
        return None

    values_sorted = sorted(values)
    count = len(values_sorted)
    midpoint = count // 2
    if count % 2:
        median = values_sorted[midpoint]
    else:
        median = (values_sorted[midpoint - 1] + values_sorted[midpoint]) / 2

    return NumericSummary(
        count,
        sum(values_sorted) / count,
        median,
        values_sorted[0],
        values_sorted[-1],
    )


def merge_exon_intervals(exons):
    """Merge overlaps before deriving gaps so duplicated/overlapping exons do not create negative introns."""
    exons_sorted = sorted(exons)
    merged = []

    for start, end in exons_sorted:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end

    return merged


def derive_intron_lengths(exons_by_transcript):
    """Derive intron lengths from gaps between consecutive exons per transcript/locus."""
    intron_lengths = []

    for transcript_exons in exons_by_transcript.values():
        exons_by_locus = defaultdict(list)
        for chrom, start, end, strand in transcript_exons:
            exons_by_locus[(chrom, strand)].append((start, end))

        for locus_exons in exons_by_locus.values():
            merged_exons = merge_exon_intervals(locus_exons)
            for idx in range(1, len(merged_exons)):
                intron_start = merged_exons[idx - 1][1]
                intron_end = merged_exons[idx][0]
                if intron_end > intron_start:
                    intron_lengths.append(intron_end - intron_start)

    return intron_lengths


def build_exon_intervals_by_chrom(exons):
    intervals_by_chrom = defaultdict(list)
    for exon in exons:
        intervals_by_chrom[exon.chrom].append((exon.start, exon.end, exon.index))

    for chrom in intervals_by_chrom:
        intervals_by_chrom[chrom].sort()

    return dict(intervals_by_chrom)


def scan_fasta_for_exon_gc(fasta_path, exons):
    """Scan FASTA once and accumulate GC/base counts for every exon interval."""
    intervals_by_chrom = build_exon_intervals_by_chrom(exons)
    gc_counts = [0] * len(exons)
    base_counts = [0] * len(exons)

    chrom = None
    pos = 0
    intervals = ()
    cursor = 0
    seen_chroms = set()

    with open_maybe_gzip(fasta_path, "rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(b">"):
                chrom = raw_line[1:].split(None, 1)[0].decode()
                pos = 0
                intervals = intervals_by_chrom.get(chrom, ())
                cursor = 0
                seen_chroms.add(chrom)
                continue

            seq = raw_line.strip()
            if not seq:
                continue

            line_start = pos
            line_end = line_start + len(seq)

            while cursor < len(intervals) and intervals[cursor][1] <= line_start:
                cursor += 1

            idx = cursor
            while idx < len(intervals) and intervals[idx][0] < line_end:
                exon_start, exon_end, exon_index = intervals[idx]
                overlap_start = max(line_start, exon_start)
                overlap_end = min(line_end, exon_end)

                if overlap_end > overlap_start:
                    seq_start = overlap_start - line_start
                    seq_end = overlap_end - line_start
                    chunk = seq[seq_start:seq_end].upper()
                    gc_counts[exon_index] += chunk.count(b"G") + chunk.count(b"C")
                    base_counts[exon_index] += (
                        chunk.count(b"A")
                        + chunk.count(b"C")
                        + chunk.count(b"G")
                        + chunk.count(b"T")
                    )

                idx += 1

            pos = line_end

    missing_chroms = sorted(set(intervals_by_chrom) - seen_chroms)
    return gc_counts, base_counts, missing_chroms


def exon_gc_percentages(gc_counts, base_counts):
    gc_percentages = []
    skipped = 0

    for gc_count, base_count in zip(gc_counts, base_counts):
        if base_count == 0:
            skipped += 1
            continue
        gc_percentages.append((gc_count / base_count) * 100)

    return gc_percentages, skipped


def gc_histogram(values, bin_size):
    num_bins = (100 + bin_size - 1) // bin_size
    counts = [0] * num_bins

    for value in values:
        bin_index = min(int(value // bin_size), num_bins - 1)
        counts[bin_index] += 1

    rows = []
    for idx, count in enumerate(counts):
        start = idx * bin_size
        end = min(100, start + bin_size)
        if end == 100:
            label = f"{start}-100%"
        else:
            label = f"{start}-<{end}%"
        rows.append((label, count))

    return rows


def intron_length_histogram(values):
    return length_histogram(values, INTRON_LENGTH_BINS)


def exon_length_histogram(values):
    return length_histogram(values, EXON_LENGTH_BINS)


def length_histogram(values, bins):
    counts = []
    for minimum, maximum in bins:
        if maximum is None:
            count = sum(1 for value in values if value >= minimum)
            label = f">={minimum} bp"
        else:
            count = sum(1 for value in values if minimum <= value <= maximum)
            label = f"{minimum}-{maximum} bp"
        counts.append((label, count))
    return counts


def frame_consistency(cds_lengths_by_transcript):
    consistent = 0
    inconsistent = 0

    for cds_length in cds_lengths_by_transcript.values():
        if cds_length % 3 == 0:
            consistent += 1
        else:
            inconsistent += 1

    return consistent, inconsistent


def collect_flagged_features(annotation_metrics, min_intron_length, min_exon_length):
    flagged = []

    for transcript_id, cds_length in sorted(
        annotation_metrics.cds_lengths_by_transcript.items()
    ):
        remainder = cds_length % 3
        if remainder:
            flagged.append(
                FlaggedFeature(
                    "frame_inconsistent_transcript",
                    transcript_id,
                    "transcript",
                    ".",
                    None,
                    None,
                    ".",
                    cds_length,
                    f"CDS length modulo 3 = {remainder}",
                )
            )

    for transcript_id, transcript_exons in sorted(
        annotation_metrics.exons_by_transcript.items()
    ):
        exons_by_locus = defaultdict(list)
        for chrom, start, end, strand in transcript_exons:
            length = end - start
            if length < min_exon_length:
                flagged.append(
                    FlaggedFeature(
                        "short_exon",
                        transcript_id,
                        "exon",
                        chrom,
                        start,
                        end,
                        strand,
                        length,
                        f"exon length < {min_exon_length} bp",
                    )
                )
            exons_by_locus[(chrom, strand)].append((start, end))

        for (chrom, strand), locus_exons in exons_by_locus.items():
            merged_exons = merge_exon_intervals(locus_exons)
            for idx in range(1, len(merged_exons)):
                intron_start = merged_exons[idx - 1][1]
                intron_end = merged_exons[idx][0]
                length = intron_end - intron_start
                if 0 < length < min_intron_length:
                    flagged.append(
                        FlaggedFeature(
                            "short_intron",
                            transcript_id,
                            "intron",
                            chrom,
                            intron_start,
                            intron_end,
                            strand,
                            length,
                            f"intron length < {min_intron_length} bp",
                        )
                    )

    flagged.extend(annotation_metrics.cds_records_without_phase)
    return flagged


def format_output_start(start):
    if start is None:
        return "."
    return start + 1


def format_output_end(end):
    if end is None:
        return "."
    return end


def write_flagged_features(flagged_features, output_path):
    with open_maybe_gzip(output_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "flag_type",
                "transcript_id",
                "feature_type",
                "seqid",
                "start",
                "end",
                "strand",
                "length",
                "details",
            )
        )
        for feature in flagged_features:
            writer.writerow(
                (
                    feature.flag_type,
                    feature.transcript_id,
                    feature.feature_type,
                    feature.chrom,
                    format_output_start(feature.start),
                    format_output_end(feature.end),
                    feature.strand,
                    feature.length,
                    feature.details,
                )
            )


def summarize_utr_type(utr_segments_by_transcript, coding_transcripts):
    """Summarize one UTR class over coding transcripts only."""
    lengths = []
    exon_counts = []
    transcripts_with_utr = set()

    for transcript_id in coding_transcripts:
        utr_exons = utr_segments_by_transcript.get(transcript_id, ())
        if not utr_exons:
            continue

        transcripts_with_utr.add(transcript_id)
        exon_counts.append(len(utr_exons))
        lengths.append(sum(end - start for _chrom, start, end, _strand in utr_exons))

    return (
        transcripts_with_utr,
        summarize_numbers(lengths),
        summarize_numbers(exon_counts),
    )


def utr_presence_counts(
    five_prime_transcripts, three_prime_transcripts, coding_transcripts
):
    both = five_prime_transcripts & three_prime_transcripts
    neither = coding_transcripts - five_prime_transcripts - three_prime_transcripts

    return {
        "Coding transcripts with 5' UTR": len(five_prime_transcripts),
        "Coding transcripts with 3' UTR": len(three_prime_transcripts),
        "Coding transcripts with both UTRs": len(both),
        "Coding transcripts with neither UTR": len(neither),
    }


def format_count_pct(count, total):
    if total == 0:
        return f"{count} (N/A)"
    return f"{count} ({(count / total) * 100:.2f}%)"


def print_metric(name, value):
    print(f"{name:<45} {value}")


def print_summary_metrics(prefix, summary, suffix="", include_count=True):
    if summary is None:
        if include_count:
            print_metric(f"{prefix} count", 0)
        print_metric(f"{prefix} mean{suffix}", "N/A")
        print_metric(f"{prefix} median{suffix}", "N/A")
        print_metric(f"{prefix} min{suffix}", "N/A")
        print_metric(f"{prefix} max{suffix}", "N/A")
        return

    if include_count:
        print_metric(f"{prefix} count", summary.count)
    print_metric(f"{prefix} mean{suffix}", f"{summary.mean:.2f}")
    print_metric(f"{prefix} median{suffix}", f"{summary.median:.2f}")
    print_metric(f"{prefix} min{suffix}", f"{summary.minimum:.2f}")
    print_metric(f"{prefix} max{suffix}", f"{summary.maximum:.2f}")


def format_filter_list(values):
    if not values:
        return "none"
    return ", ".join(sorted(values))


def print_utr_metrics(annotation_metrics, coding_transcripts):
    (
        five_prime_transcripts,
        five_prime_length_summary,
        five_prime_count_summary,
    ) = summarize_utr_type(
        annotation_metrics.utrs_by_transcript.get("five_prime", {}),
        coding_transcripts,
    )
    (
        three_prime_transcripts,
        three_prime_length_summary,
        three_prime_count_summary,
    ) = summarize_utr_type(
        annotation_metrics.utrs_by_transcript.get("three_prime", {}),
        coding_transcripts,
    )
    presence_counts = utr_presence_counts(
        five_prime_transcripts, three_prime_transcripts, coding_transcripts
    )
    coding_transcript_count = len(coding_transcripts)

    for label, count in presence_counts.items():
        print_metric(label, format_count_pct(count, coding_transcript_count))

    for label, length_summary, count_summary in (
        (UTR_LABELS["five_prime"], five_prime_length_summary, five_prime_count_summary),
        (
            UTR_LABELS["three_prime"],
            three_prime_length_summary,
            three_prime_count_summary,
        ),
    ):
        print_summary_metrics(
            f"{label} length per transcript",
            length_summary,
            " (bp)",
            include_count=False,
        )
        print_summary_metrics(
            f"{label} exon count per transcript", count_summary, include_count=False
        )


def main():
    parser = argparse.ArgumentParser(
        description="Calculate feature-level QC metrics for genome annotations."
    )
    parser.add_argument(
        "gff", help="Path to the input GFF3/GTF file containing annotations"
    )
    parser.add_argument(
        "--fasta", required=True, help="Path to the soft-masked genome FASTA file."
    )
    parser.add_argument(
        "--format",
        "--annotation-format",
        choices=("auto", "gff3", "gff", "gtf"),
        default="auto",
        dest="annotation_format",
        help="Annotation format (default: auto)",
    )
    parser.add_argument(
        "--gc-bin-size",
        type=gc_bin_size,
        default=10,
        help="GC histogram bin size in percentage points (default: 10)",
    )
    parser.add_argument(
        "--biotype",
        action="append",
        help="Only include transcripts with one of these biotypes. Can be comma-separated or repeated.",
    )
    parser.add_argument(
        "--exclude-biotype",
        action="append",
        help="Exclude transcripts with one of these biotypes. Can be comma-separated or repeated.",
    )
    parser.add_argument(
        "--canonical-only",
        action="store_true",
        help="Only include transcripts tagged as canonical, Ensembl_canonical, or MANE_Select.",
    )
    parser.add_argument(
        "--flagged-features",
        help="Write a TSV of potentially complicated features/transcripts to this path.",
    )
    parser.add_argument(
        "--min-intron-length",
        type=positive_int,
        default=20,
        help="Flag introns shorter than this many bp in --flagged-features output (default: 20).",
    )
    parser.add_argument(
        "--min-exon-length",
        type=positive_int,
        default=3,
        help="Flag exons shorter than this many bp in --flagged-features output (default: 3).",
    )
    args = parser.parse_args()

    include_biotypes = parse_filter_values(args.biotype)
    exclude_biotypes = parse_filter_values(args.exclude_biotype)
    annotation_metrics = parse_annotation(
        args.gff,
        args.annotation_format,
        include_biotypes,
        exclude_biotypes,
        args.canonical_only,
    )
    if not annotation_metrics.exons:
        print("[ERROR] No exon features found in GFF3/GTF.", file=sys.stderr)
        sys.exit(1)

    gc_counts, base_counts, missing_chroms = scan_fasta_for_exon_gc(
        args.fasta, annotation_metrics.exons
    )
    gc_percentages, gc_skipped = exon_gc_percentages(gc_counts, base_counts)
    exon_gc_summary = summarize_numbers(gc_percentages)
    exon_lengths = [exon.end - exon.start for exon in annotation_metrics.exons]
    exon_length_summary = summarize_numbers(exon_lengths)
    intron_lengths = derive_intron_lengths(annotation_metrics.exons_by_transcript)
    intron_summary = summarize_numbers(intron_lengths)
    frame_consistent, frame_inconsistent = frame_consistency(
        annotation_metrics.cds_lengths_by_transcript
    )
    flagged_features = []

    if args.flagged_features:
        flagged_features = collect_flagged_features(
            annotation_metrics,
            args.min_intron_length,
            args.min_exon_length,
        )
        write_flagged_features(flagged_features, args.flagged_features)

    if annotation_metrics.exon_phase_counts:
        phase_counts = annotation_metrics.exon_phase_counts
        phase_total = sum(phase_counts.values())
        phase_source = "exon phase column"
        phase_prefix = "Exon"
        phase_record_label = "Exon records"
        phase_missing = annotation_metrics.exon_phase_missing
    else:
        # Exon phase is commonly absent; CDS phase is a practical coding-exon proxy.
        phase_counts = annotation_metrics.cds_phase_counts
        phase_total = sum(phase_counts.values())
        phase_source = "CDS phase column"
        phase_prefix = "CDS"
        phase_record_label = "CDS records"
        phase_missing = annotation_metrics.cds_phase_missing

    coding_transcripts = len(annotation_metrics.cds_lengths_by_transcript)
    frame_total = frame_consistent + frame_inconsistent

    print(f"{'Metric':<45} {'Value'}")
    print("-" * 65)
    print_metric("Annotation format", annotation_metrics.annotation_format)
    if annotation_metrics.filters_applied:
        print_metric("Included biotypes", format_filter_list(include_biotypes))
        print_metric("Excluded biotypes", format_filter_list(exclude_biotypes))
        print_metric(
            "Canonical-only filter", "enabled" if args.canonical_only else "disabled"
        )
        print_metric(
            "Transcripts passing filters", annotation_metrics.selected_transcript_count
        )

    print("-" * 65)
    print_metric("Exon count", len(annotation_metrics.exons))
    print_metric(
        "Transcripts with exon records", len(annotation_metrics.exons_by_transcript)
    )
    print_metric("Coding transcripts with CDS records", coding_transcripts)

    print_metric("Phase source", phase_source if phase_total else "N/A")
    print_metric(f"{phase_record_label} with usable phase", phase_total)
    for phase, label in PHASE_LABELS.items():
        print_metric(
            f"{phase_prefix} phase {label}",
            format_count_pct(phase_counts[phase], phase_total),
        )
    print_metric(
        f"{phase_record_label} without usable phase",
        phase_missing if phase_total else "N/A",
    )

    print_summary_metrics(
        "Exon length", exon_length_summary, " (bp)", include_count=False
    )
    for label, count in exon_length_histogram(exon_lengths):
        print_metric(f"Exon length {label}", format_count_pct(count, len(exon_lengths)))

    print_summary_metrics("Exon GC", exon_gc_summary, " (%)")
    print_metric("Exons skipped for GC", gc_skipped)
    if missing_chroms:
        print_metric("FASTA contigs missing annotation exons", len(missing_chroms))
    for label, count in gc_histogram(gc_percentages, args.gc_bin_size):
        print_metric(f"Exon GC {label}", format_count_pct(count, len(gc_percentages)))

    print_summary_metrics("Intron length", intron_summary, " (bp)")
    for label, count in intron_length_histogram(intron_lengths):
        print_metric(
            f"Intron length {label}", format_count_pct(count, len(intron_lengths))
        )

    print_utr_metrics(
        annotation_metrics, set(annotation_metrics.cds_lengths_by_transcript)
    )

    print_metric(
        "Frame-consistent coding transcripts",
        format_count_pct(frame_consistent, frame_total),
    )
    print_metric(
        "Frame-inconsistent coding transcripts",
        format_count_pct(frame_inconsistent, frame_total),
    )

    print("-" * 65)
    if args.flagged_features:
        print_metric("Flagged features output", args.flagged_features)
        print_metric("Flagged records written", len(flagged_features))
        print_metric("Short intron threshold (bp)", f"< {args.min_intron_length}")
        print_metric("Short exon threshold (bp)", f"< {args.min_exon_length}")


if __name__ == "__main__":
    main()

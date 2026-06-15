#!/usr/bin/env python3
"""
Find BUSCOs complete in a reference annotation but not complete in a test
annotation, then compare their reference GFF3 structure against the test GFF3.

This is a two-annotation version of find_missing_busco_reason.py. It has no
evidence/transcriptomic/protein-layer stack. Optionally, one additional GFF3
layer can be supplied and checked for the same missing BUSCO structures.
"""

import argparse
import csv
import re
from collections import Counter, defaultdict
from urllib.parse import unquote


COMPLETE_STATUSES = {"Complete", "Duplicated"}


def parse_busco_table(path):
    """Parse BUSCO full_table.tsv into BUSCO_ID -> BUSCO info."""
    buscos = {}

    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            busco_id = parts[0]
            sequence = parts[2] if len(parts) > 2 and parts[2] != "-" else None

            buscos[busco_id] = {
                "status": parts[1],
                "sequence": sequence,
                "score": parts[3] if len(parts) > 3 else None,
                "length": parts[4] if len(parts) > 4 else None,
            }

    return buscos


def parse_gff3_attributes(attr_string):
    """Parse the GFF3 attributes column into key -> list(values)."""
    attrs = defaultdict(list)

    for item in attr_string.strip().split(";"):
        if not item or "=" not in item:
            continue

        key, value = item.split("=", 1)
        key = unquote(key.strip())
        value = unquote(value.strip())
        attrs[key].append(value)

    return attrs


def split_multi_value(value):
    for part in value.split(","):
        part = part.strip()
        if part:
            yield part


def id_aliases(value):
    """Return lookup aliases for IDs that may differ by prefixes/version suffixes."""
    if value is None:
        return set()

    value = str(value).strip().strip('"')
    if not value:
        return set()

    aliases = {value}

    first_token = value.split()[0]
    aliases.add(first_token)

    if ":" in first_token:
        aliases.add(first_token.rsplit(":", 1)[-1])

    for alias in list(aliases):
        if "." in alias:
            base, suffix = alias.rsplit(".", 1)
            if suffix.isdigit():
                aliases.add(base)

    return {alias for alias in aliases if alias}


def add_id_map(id_map, value, transcript_id):
    for alias in id_aliases(value):
        id_map[alias] = transcript_id


def transcript_id_from_attrs(attrs):
    for key in ("ID", "transcript_id"):
        if attrs.get(key):
            return attrs[key][0]
    return None


def parent_ids_from_attrs(attrs):
    parents = []

    for value in attrs.get("Parent", []):
        parents.extend(split_multi_value(value))

    if parents:
        return parents

    if attrs.get("transcript_id"):
        return [attrs["transcript_id"][0]]

    return []


def parse_gff3_annotation(path, cds_only=False):
    """
    Parse GFF3 into transcript structures and ID mappings.

    Returns:
        transcripts: transcript_id -> structure dict
        id_map: protein/transcript/feature aliases -> transcript_id
    """
    transcripts = defaultdict(lambda: {
        "chrom": None,
        "strand": None,
        "start": None,
        "end": None,
        "exons": [],
        "cds_exons": [],
        "feature_ids": set(),
    })
    id_map = {}

    with open(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _source, feature, start, end, _score, strand, _phase, attr_string = parts

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue

            attrs = parse_gff3_attributes(attr_string)
            feature_lower = feature.lower()

            transcript_features = {
                "mrna",
                "transcript",
                "rna",
                "primary_transcript",
                "processed_transcript",
            }

            if feature_lower in transcript_features:
                transcript_ids = [transcript_id_from_attrs(attrs)]
            elif feature_lower in {"exon", "cds"}:
                transcript_ids = parent_ids_from_attrs(attrs)
            else:
                transcript_ids = []

            transcript_ids = [tid for tid in transcript_ids if tid]
            if not transcript_ids:
                continue

            for transcript_id in transcript_ids:
                transcript = transcripts[transcript_id]
                transcript["chrom"] = chrom
                transcript["strand"] = strand

                if feature_lower == "exon":
                    transcript["exons"].append((start, end))
                elif feature_lower == "cds":
                    transcript["cds_exons"].append((start, end))

                for key in ("ID", "Name", "protein_id", "transcript_id", "Target"):
                    for value in attrs.get(key, []):
                        for part in split_multi_value(value):
                            transcript["feature_ids"].add(part.split()[0])
                            add_id_map(id_map, part, transcript_id)

                add_id_map(id_map, transcript_id, transcript_id)

    for transcript_id in list(transcripts.keys()):
        transcript = transcripts[transcript_id]

        if cds_only:
            structure = sorted(transcript["cds_exons"])
        else:
            structure = sorted(transcript["exons"] or transcript["cds_exons"])

        if not structure:
            del transcripts[transcript_id]
            continue

        transcript["exons"] = structure
        transcript["start"] = min(start for start, _end in structure)
        transcript["end"] = max(end for _start, end in structure)

    return transcripts, id_map


def map_buscos_to_reference(buscos, transcripts, id_map):
    """Attach reference transcript coordinates/structure to BUSCO records."""
    for busco_id, info in buscos.items():
        transcript_id = None

        for alias in id_aliases(info.get("sequence")):
            if alias in id_map:
                transcript_id = id_map[alias]
                break

        if transcript_id and transcript_id in transcripts:
            transcript = transcripts[transcript_id]
            info["ref_transcript"] = transcript_id
            info["chrom"] = transcript["chrom"]
            info["start"] = transcript["start"]
            info["end"] = transcript["end"]
            info["strand"] = transcript["strand"]
            info["ref_exons"] = sorted(transcript["exons"])
        else:
            info["ref_transcript"] = None
            info["chrom"] = None
            info["start"] = None
            info["end"] = None
            info["strand"] = None
            info["ref_exons"] = []

    return buscos


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def clip_exons(exons, region_start, region_end):
    clipped = []

    for exon_start, exon_end in exons:
        clipped_start = max(exon_start, region_start)
        clipped_end = min(exon_end, region_end)
        if clipped_start <= clipped_end:
            clipped.append((clipped_start, clipped_end))

    return sorted(clipped)


def exonic_length(exons):
    return sum(end - start + 1 for start, end in exons)


def intron_count(exons):
    return max(0, len(exons) - 1)


def exon_span(exons):
    if not exons:
        return None
    return min(start for start, _end in exons), max(end for _start, end in exons)


def contained_by(container, interval):
    return container[0] <= interval[0] and interval[1] <= container[1]


def compatible_short_geometry(test_exons, ref_exons):
    """
    Return True when test_exons are an exact boundary-compatible subset of
    ref_exons.

    This allows terminal truncation, but it does not allow shifted internal
    exon boundaries, exon splits inside a reference exon, merged reference
    exons, skipped internal exons, or exonic sequence across reference introns.
    """
    if not test_exons or tuple(test_exons) == tuple(ref_exons):
        return False

    mapped = []
    for test_exon in test_exons:
        matching_ref_indexes = [
            index
            for index, ref_exon in enumerate(ref_exons)
            if contained_by(ref_exon, test_exon)
        ]

        if len(matching_ref_indexes) != 1:
            return False

        mapped.append((matching_ref_indexes[0], test_exon))

    indexes = [index for index, _test_exon in mapped]
    if indexes != sorted(indexes) or len(indexes) != len(set(indexes)):
        return False

    if len(indexes) > 1:
        expected = list(range(indexes[0], indexes[-1] + 1))
        if indexes != expected:
            return False

    for position, (ref_index, test_exon) in enumerate(mapped):
        ref_start, ref_end = ref_exons[ref_index]
        test_start, test_end = test_exon

        is_first = position == 0
        is_last = position == len(mapped) - 1

        if not is_first and test_start != ref_start:
            return False
        if not is_last and test_end != ref_end:
            return False
        if not is_first and not is_last and test_exon != (ref_start, ref_end):
            return False

    return True


def classify_candidate(raw_exons, clipped_exons, ref_exons):
    raw_exons = tuple(sorted(raw_exons))
    clipped_exons = tuple(sorted(clipped_exons))
    ref_exons = tuple(sorted(ref_exons))

    if tuple(raw_exons) == tuple(ref_exons):
        return "PRESENT_FULL"

    if tuple(clipped_exons) == tuple(ref_exons):
        return "PRESENT_LONG"

    ref_span = exon_span(ref_exons)
    raw_span = exon_span(raw_exons)
    if raw_span == ref_span and intron_count(raw_exons) != intron_count(ref_exons):
        return "PRESENT_INTRONS"

    if compatible_short_geometry(tuple(clipped_exons), tuple(ref_exons)):
        if intron_count(clipped_exons) != intron_count(ref_exons):
            return "PRESENT_SHORT_INTRONS"
        return "PRESENT_SHORT"

    return None


def status_priority(status):
    if status == "PRESENT_FULL":
        return 5
    if status == "PRESENT_LONG":
        return 4
    if status == "PRESENT_SHORT":
        return 3
    if status in {"PRESENT_INTRONS", "PRESENT_SHORT_INTRONS"}:
        return 2
    return 0


def classify_region(chrom, start, end, ref_exons, test_transcripts):
    """
    Compare reference BUSCO exon/CDS structure against test transcripts.

    Geometry rules:
        PRESENT_FULL requires the complete exon list to match the reference
            exactly, including all exon boundaries and therefore all introns.
        PRESENT_INTRONS requires the same outer start/end as the reference but
            a different number of introns.
        PRESENT_SHORT requires an exact boundary-compatible truncation of the
            reference exon geometry at the beginning and/or end, with the same
            number of introns as the reference.
        PRESENT_SHORT_INTRONS uses the same short/truncated geometry rule, but
            with more or fewer introns than the reference.
        PRESENT_LONG requires the reference exon list to match exactly after
            clipping the candidate to the reference BUSCO span, with extra
            exonic sequence outside that span on either side.
        ABSENT means no transcript overlaps the reference BUSCO region.
        OVERLAP means at least one transcript overlaps the region, but none
            match the requested PRESENT_* categories.
    """
    if chrom is None or start is None or end is None:
        return "UNMAPPED_REFERENCE", 0, 0, 0, ""

    ref_exons = tuple(sorted(ref_exons))
    ref_len = exonic_length(ref_exons)
    if ref_len == 0:
        return "NO_REFERENCE_STRUCTURE", 0, 0, 0, ""

    ref_introns = max(0, len(ref_exons) - 1)
    candidates = []
    overlaps = []

    for transcript_id, transcript in test_transcripts.items():
        if transcript["chrom"] != chrom:
            continue

        clipped = clip_exons(transcript["exons"], start, end)
        if not clipped:
            continue

        introns = max(0, len(clipped) - 1)
        exon_overlap = 0

        for ref_start, ref_end in ref_exons:
            for test_start, test_end in clipped:
                exon_overlap += overlap(ref_start, ref_end, test_start, test_end)

        coverage = round(exon_overlap / ref_len * 100, 1)
        overlaps.append((transcript_id, coverage, introns, exon_overlap))

        status = classify_candidate(tuple(transcript["exons"]), tuple(clipped), ref_exons)
        if status is None:
            continue
        candidates.append((transcript_id, status, coverage, introns, exon_overlap))

    if not candidates:
        if overlaps:
            transcript_id, coverage, introns, _exon_overlap = max(
                overlaps,
                key=lambda item: (item[1], item[3], -abs(item[2] - ref_introns)),
            )
            return "OVERLAP", coverage, introns, len(overlaps), transcript_id
        return "ABSENT", 0, 0, 0, ""

    best = max(
        candidates,
        key=lambda item: (status_priority(item[1]), item[2], -abs(item[3] - ref_introns)),
    )
    transcript_id, status, coverage, introns, _exon_overlap = best

    return status, coverage, introns, len(overlaps), transcript_id


def safe_column_prefix(name):
    prefix = re.sub(r"[^A-Za-z0-9]+", "_", name.strip()).strip("_").lower()
    return prefix or "layer"


def add_structure_classification(row, prefix, chrom, start, end, ref_exons, transcripts):
    status, coverage, introns, hits, transcript_names = classify_region(
        chrom,
        start,
        end,
        ref_exons,
        transcripts,
    )

    row[f"{prefix}_status"] = status
    row[f"{prefix}_coverage_pct"] = coverage
    row[f"{prefix}_introns"] = introns
    row[f"{prefix}_hits"] = hits
    row[f"{prefix}_transcripts"] = transcript_names


def find_missing_buscos(ref_buscos, test_buscos, strict_missing=False):
    complete_ref = {
        busco_id
        for busco_id, info in ref_buscos.items()
        if info["status"] in COMPLETE_STATUSES
    }

    if strict_missing:
        return sorted(
            busco_id
            for busco_id in complete_ref
            if test_buscos.get(busco_id, {}).get("status") == "Missing"
        )

    complete_test = {
        busco_id
        for busco_id, info in test_buscos.items()
        if info["status"] in COMPLETE_STATUSES
    }
    return sorted(complete_ref - complete_test)


def build_rows(ref_buscos, test_buscos, missing_buscos, test_transcripts, extra_layers=None):
    rows = []
    extra_layers = extra_layers or []

    for busco_id in missing_buscos:
        ref_info = ref_buscos[busco_id]
        test_info = test_buscos.get(busco_id, {})
        chrom = ref_info.get("chrom")
        start = ref_info.get("start")
        end = ref_info.get("end")
        ref_exons = clip_exons(ref_info.get("ref_exons", []), start, end) if chrom else []

        row = {
            "BUSCO_ID": busco_id,
            "Reference_BUSCO_status": ref_info.get("status"),
            "Test_BUSCO_status": test_info.get("status", "Not found"),
            "Reference_sequence": ref_info.get("sequence"),
            "Test_sequence": test_info.get("sequence"),
            "Reference_Transcript": ref_info.get("ref_transcript"),
            "Chrom": chrom,
            "Start": start,
            "End": end,
            "Strand": ref_info.get("strand"),
            "Reference_exon_count": len(ref_exons),
            "Reference_exonic_length": exonic_length(ref_exons),
        }

        add_structure_classification(
            row,
            "test",
            chrom,
            start,
            end,
            ref_exons,
            test_transcripts,
        )

        for layer_name, layer_transcripts in extra_layers:
            add_structure_classification(
                row,
                layer_name,
                chrom,
                start,
                end,
                ref_exons,
                layer_transcripts,
            )

        rows.append(row)

    return rows


def structure_fieldnames(prefix):
    return [
        f"{prefix}_status",
        f"{prefix}_coverage_pct",
        f"{prefix}_introns",
        f"{prefix}_hits",
        f"{prefix}_transcripts",
    ]


def write_rows(path, rows, extra_layer_names=None):
    extra_layer_names = extra_layer_names or []

    fieldnames = [
        "BUSCO_ID",
        "Reference_BUSCO_status",
        "Test_BUSCO_status",
        "Reference_sequence",
        "Test_sequence",
        "Reference_Transcript",
        "Chrom",
        "Start",
        "End",
        "Strand",
        "Reference_exon_count",
        "Reference_exonic_length",
    ]
    fieldnames.extend(structure_fieldnames("test"))

    for layer_name in extra_layer_names:
        fieldnames.extend(structure_fieldnames(layer_name))

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def compare_missing_busco_structures(
    busco_ref,
    busco_test,
    ref_gff3,
    test_gff3,
    output,
    cds_only=False,
    strict_missing=False,
    layer_gff3=None,
    layer_name="provenance",
):
    print(f"Parsing reference BUSCO table: {busco_ref}")
    ref_buscos = parse_busco_table(busco_ref)

    print(f"Parsing test BUSCO table: {busco_test}")
    test_buscos = parse_busco_table(busco_test)

    missing_buscos = find_missing_buscos(ref_buscos, test_buscos, strict_missing=strict_missing)
    print(f"BUSCOs complete in reference but not complete in test: {len(missing_buscos)}")

    print(f"Parsing reference GFF3: {ref_gff3}")
    ref_transcripts, ref_id_map = parse_gff3_annotation(ref_gff3, cds_only=cds_only)
    ref_buscos = map_buscos_to_reference(ref_buscos, ref_transcripts, ref_id_map)

    print(f"Parsing test GFF3: {test_gff3}")
    test_transcripts, _test_id_map = parse_gff3_annotation(test_gff3, cds_only=cds_only)

    extra_layers = []
    extra_layer_names = []
    if layer_gff3:
        layer_prefix = safe_column_prefix(layer_name)
        print(f"Parsing {layer_prefix} GFF3 layer: {layer_gff3}")
        layer_transcripts, _layer_id_map = parse_gff3_annotation(layer_gff3, cds_only=cds_only)
        extra_layers.append((layer_prefix, layer_transcripts))
        extra_layer_names.append(layer_prefix)

    rows = build_rows(
        ref_buscos,
        test_buscos,
        missing_buscos,
        test_transcripts,
        extra_layers=extra_layers,
    )
    write_rows(output, rows, extra_layer_names=extra_layer_names)

    print(f"\nWritten {output} ({len(rows)} BUSCOs)")

    status_counts = Counter(row["test_status"] for row in rows)
    if status_counts:
        print("\n=== TEST STATUS SUMMARY ===")
        for status, count in status_counts.most_common():
            print(f"{status}\t{count}")

    for layer_prefix in extra_layer_names:
        layer_counts = Counter(row[f"{layer_prefix}_status"] for row in rows)
        if layer_counts:
            print(f"\n=== {layer_prefix.upper()} STATUS SUMMARY ===")
            for status, count in layer_counts.most_common():
                print(f"{status}\t{count}")

    unmapped = sum(1 for row in rows if row["test_status"] == "UNMAPPED_REFERENCE")
    if unmapped:
        print(
            "\nWARNING: some reference BUSCO sequences could not be mapped to "
            "reference GFF3 transcripts. Check Reference_sequence and "
            "Reference_Transcript in the output."
        )

    return rows


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Find BUSCOs complete in a reference BUSCO run but not complete in a "
            "test BUSCO run, then classify their structure in a test GFF3."
        )
    )
    parser.add_argument("--busco_ref", required=True, help="Reference BUSCO full_table.tsv")
    parser.add_argument("--busco_test", required=True, help="Test BUSCO full_table.tsv")
    parser.add_argument("--ref_gff3", "--ref_annotation", required=True, help="Reference GFF3")
    parser.add_argument("--test_gff3", "--test_annotation", required=True, help="Test GFF3")
    parser.add_argument(
        "--layer_gff3",
        "--provenance_gff3",
        "--layer_db",
        help=(
            "Optional additional GFF3 layer to check for the same missing BUSCO "
            "structures, for example the provenance layer DB."
        ),
    )
    parser.add_argument(
        "--layer_name",
        default="provenance",
        help="Column prefix for --layer_gff3 results",
    )
    parser.add_argument(
        "--out",
        default="missing_busco_annotation_presence.csv",
        help="Output CSV table",
    )
    parser.add_argument(
        "--cds_only",
        action="store_true",
        default=False,
        help=(
            "Use CDS structure instead of exon structure for both reference and test. "
            "Useful when UTRs differ or when GFF3 exon features are absent."
        ),
    )
    parser.add_argument(
        "--strict_missing",
        action="store_true",
        default=False,
        help=(
            "Only include BUSCOs with test status exactly Missing. Default includes "
            "anything complete in reference but not Complete/Duplicated in test."
        ),
    )

    args = parser.parse_args()

    compare_missing_busco_structures(
        busco_ref=args.busco_ref,
        busco_test=args.busco_test,
        ref_gff3=args.ref_gff3,
        test_gff3=args.test_gff3,
        output=args.out,
        cds_only=args.cds_only,
        strict_missing=args.strict_missing,
        layer_gff3=args.layer_gff3,
        layer_name=args.layer_name,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
BUSCO Missing Annotation Scanner

For BUSCOs complete in reference but missing in test:
    • map BUSCO → reference transcript + coordinates
    • check multiple GTF/GFF files
    • detect presence + introns
    • output summary table

Output:
    missing_busco_annotation_presence.csv
"""

import argparse
from pathlib import Path
from collections import defaultdict
import pandas as pd


# ---------------------------------------------------
# BUSCO parsing
# ---------------------------------------------------

def parse_busco_table(busco_file):
    buscos = {}
    with open(busco_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            buscos[parts[0]] = {
                "status": parts[1],
                "sequence": parts[2] if len(parts) > 2 else None
            }
    return buscos


# ---------------------------------------------------
# Attribute parser (GTF + GFF3)
# ---------------------------------------------------

def parse_attributes(attr_string):
    attrs = {}

    if "=" in attr_string:  # GFF3
        for a in attr_string.split(";"):
            if "=" in a:
                k, v = a.split("=", 1)
                attrs[k.strip()] = v.strip()
    else:  # GTF
        for a in attr_string.split(";"):
            a = a.strip()
            if not a:
                continue
            k, v = a.split(" ", 1)
            attrs[k.strip()] = v.strip('"')

    return attrs


# ---------------------------------------------------
# Parse annotation → transcript models
# ---------------------------------------------------

def parse_annotation(annotation_file):

    transcripts = defaultdict(lambda: {
        "chrom": None,
        "strand": None,
        "start": 10**12,
        "end": 0,
        "exons": [],
        "name": None
    })

    id_map = {}

    with open(annotation_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)

            attr = parse_attributes(attrs)

            tid = (
                attr.get("transcript_id")
                or attr.get("ID")
                or attr.get("protein_id")
            )

            pid = attr.get("protein_id")

            if not tid:
                continue

            t = transcripts[tid]

            t["chrom"] = chrom
            t["strand"] = strand
            t["start"] = min(t["start"], start)
            t["end"] = max(t["end"], end)
            t["name"] = tid

            if feature.lower() == "exon":
                t["exons"].append((start, end))

            id_map[tid] = tid
            if pid:
                id_map[pid] = tid

    return transcripts, id_map


# ---------------------------------------------------
# Map BUSCO → reference transcript + coords
# ---------------------------------------------------

def map_buscos_to_reference(buscos, transcripts, id_map):

    for bid, info in buscos.items():
        seq = info["sequence"]

        if seq in id_map:
            tid = id_map[seq]
            t = transcripts[tid]

            info["ref_transcript"] = tid
            info["chrom"] = t["chrom"]
            info["start"] = t["start"]
            info["end"] = t["end"]
            info["strand"] = t["strand"]
        else:
            info["ref_transcript"] = None
            info["chrom"] = info["start"] = info["end"] = info["strand"] = None

    return buscos


# ---------------------------------------------------
# Presence classification
# ---------------------------------------------------

def classify_presence(chrom, start, end, transcripts):

    for t in transcripts.values():

        if t["chrom"] != chrom:
            continue

        if not (t["end"] < start or t["start"] > end):

            if len(t["exons"]) <= 1:
                return "PRESENT_FULL"
            else:
                return "PRESENT_INTRONS"

    return "ABSENT"


# ---------------------------------------------------
# Main analysis
# ---------------------------------------------------

def analyze(buscos_ref, buscos_test, reference_anno, other_annos, names):

    complete_ref = {b for b, i in buscos_ref.items() if i["status"] in ["Complete", "Duplicated"]}
    complete_test = {b for b, i in buscos_test.items() if i["status"] in ["Complete", "Duplicated"]}

    missing = sorted(complete_ref - complete_test)

    print(f"Missing BUSCOs to investigate: {len(missing)}")

    ref_tx, ref_map = parse_annotation(reference_anno)
    buscos_ref = map_buscos_to_reference(buscos_ref, ref_tx, ref_map)

    parsed_annos = [parse_annotation(f)[0] for f in other_annos]

    rows = []

    for bid in missing:

        info = buscos_ref[bid]

        chrom = info["chrom"]
        start = info["start"]
        end = info["end"]

        row = {
            "BUSCO_ID": bid,
            "Reference_Transcript": info["ref_transcript"],  # NEW COLUMN
            "Chrom": chrom,
            "Start": start,
            "End": end
        }

        for name, transcripts in zip(names, parsed_annos):
            status = "ABSENT"
            if chrom:
                status = classify_presence(chrom, start, end, transcripts)
            row[name] = status

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv("missing_busco_annotation_presence.csv", index=False)
    print("✓ missing_busco_annotation_presence.csv written")


# ---------------------------------------------------
# CLI
# ---------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--busco_ref", required=True)
    parser.add_argument("--busco_test", required=True)

    parser.add_argument("--ref_annotation", required=True)

    parser.add_argument("--transcriptomic", required=True)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--protein_sel", required=True)
    parser.add_argument("--uniprot", required=True)
    parser.add_argument("--uniprot_sel", required=True)
    parser.add_argument("--test_annotation", required=True)

    args = parser.parse_args()

    buscos_ref = parse_busco_table(args.busco_ref)
    buscos_test = parse_busco_table(args.busco_test)

    files = [
        args.transcriptomic,
        args.protein,
        args.protein_sel,
        args.uniprot,
        args.uniprot_sel,
        args.test_annotation
    ]

    names = [
        "transcriptomic",
        "protein",
        "protein_sel",
        "uniprot",
        "uniprot_sel",
        "test"
    ]

    analyze(
        buscos_ref,
        buscos_test,
        args.ref_annotation,
        files,
        names
    )
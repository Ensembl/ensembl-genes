#!/usr/bin/env python3
"""
BUSCO Missing Diagnostic Analysis (Coordinate-based)
Compare BUSCO results from funannotate (canonical) and anno (GenBlastG) alignments.
Detect missing BUSCOs in funannotate and potential intron insertion/truncation using GFF3 coordinates.
"""

import argparse
from collections import defaultdict
import pandas as pd

# -----------------------------
# BUSCO parsing
# -----------------------------
def parse_busco_table(busco_file):
    """Parse BUSCO full_table.tsv"""
    buscos = {}
    with open(busco_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            bid = parts[0]
            status = parts[1]
            sequence = parts[2] if len(parts) > 2 else None
            buscos[bid] = {'status': status, 'sequence': sequence}
    return buscos

# -----------------------------
# GFF parsing
# -----------------------------
def parse_gff_cds(gff_file):
    """
    Parse GFF3 into transcripts with CDS coordinates.
    Returns a dict keyed by transcript ID with chromosome, strand, CDS list, total CDS length.
    """
    transcripts = defaultdict(lambda: {"chrom": None, "strand": None, "cds_list": [], "cds_len": 0})
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)
            if feature.lower() != 'cds':
                continue

            # Extract transcript/protein ID
            feature_id = None
            for key in ['ID', 'Name', 'protein_id', 'transcript_id']:
                if f'{key}=' in attributes:
                    for attr in attributes.split(';'):
                        if attr.startswith(f'{key}='):
                            feature_id = attr.split('=')[1]
                            break
                    if feature_id:
                        break
            if not feature_id:
                continue

            transcripts[feature_id]['chrom'] = chrom
            transcripts[feature_id]['strand'] = strand
            transcripts[feature_id]['cds_list'].append((start, end))
            transcripts[feature_id]['cds_len'] += (end - start + 1)
    return transcripts

def get_overlap_region(a_start, a_end, b_start, b_end):
    """Return length of overlap between two regions"""
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    return max(0, overlap_end - overlap_start + 1)

def find_best_overlap(fun_tx, anno_tx_dict):
    """
    Given a funannotate transcript and anno transcripts, find the best overlapping transcript
    on the same chromosome and strand.
    """
    best_id = None
    best_overlap = 0
    for tid, info in anno_tx_dict.items():
        if info['chrom'] != fun_tx['chrom'] or info['strand'] != fun_tx['strand']:
            continue
        total_overlap = 0
        for f_start, f_end in fun_tx['cds_list']:
            for a_start, a_end in info['cds_list']:
                total_overlap += get_overlap_region(f_start, f_end, a_start, a_end)
        if total_overlap > best_overlap:
            best_overlap = total_overlap
            best_id = tid
    return best_id, best_overlap

# -----------------------------
# Analysis function
# -----------------------------
def analyze_missing_buscos(busco_fun, busco_anno, gff_fun, gff_anno):
    """
    Compare BUSCOs missing in funannotate vs anno.
    """
    fun_tx = parse_gff_cds(gff_fun)
    anno_tx = parse_gff_cds(gff_anno)

    # BUSCOs complete in anno but missing in funannotate
    complete_anno = {bid for bid, info in busco_anno.items() if info['status'] in ['Complete', 'Duplicated']}
    complete_fun = {bid for bid, info in busco_fun.items() if info['status'] in ['Complete', 'Duplicated']}
    missing_buscos = complete_anno - complete_fun

    print(f"Found {len(missing_buscos)} BUSCOs complete in anno but missing in funannotate.")

    diagnostics = []
    for bid in sorted(missing_buscos):
        seq_id_fun = busco_fun.get(bid, {}).get('sequence')
        seq_id_anno = busco_anno[bid]['sequence']

        fun_transcript = fun_tx.get(seq_id_fun, None)
        if not fun_transcript:
            fun_transcript = {'chrom': None, 'strand': None, 'cds_list': [], 'cds_len': 0}

        rescued_id, best_overlap = find_best_overlap(fun_transcript, anno_tx)
        rescued_tx = anno_tx.get(rescued_id, None)
        if not rescued_tx:
            rescued_tx = {'chrom': None, 'strand': None, 'cds_list': [], 'cds_len': 0}

        cds_count_fun = len(fun_transcript['cds_list'])
        cds_count_rescued = len(rescued_tx['cds_list'])
        cds_len_fun = fun_transcript['cds_len']
        cds_len_rescued = rescued_tx['cds_len']

        reasons = []
        if not rescued_tx['cds_list']:
            reasons.append("No transcript found in anno")
        elif cds_count_fun > cds_count_rescued:
            reasons.append("Funannotate transcript has more CDS blocks → possible inserted intron")
        elif cds_len_fun < cds_len_rescued * 0.8:
            reasons.append("Funannotate CDS shorter → truncation")
        else:
            reasons.append("Unknown / minor differences")

        diagnostics.append({
            'BUSCO_ID': bid,
            'Funannotate_transcript': seq_id_fun,
            'Anno_transcript': rescued_id,
            'Funannotate_CDS_count': cds_count_fun,
            'Anno_CDS_count': cds_count_rescued,
            'Funannotate_CDS_length': cds_len_fun,
            'Anno_CDS_length': cds_len_rescued,
            'Likely_reason': "; ".join(reasons)
        })

    df = pd.DataFrame(diagnostics)
    df.to_csv("busco_fun_vs_anno_diagnostics.csv", index=False)
    print("✓ Diagnostic table saved to busco_fun_vs_anno_diagnostics.csv")
    return df

# -----------------------------
# CLI
# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare BUSCO missing in funannotate vs anno using GFF coordinates")
    parser.add_argument("--busco_fun", required=True, help="BUSCO full_table.tsv for funannotate proteins")
    parser.add_argument("--busco_anno", required=True, help="BUSCO full_table.tsv for anno (GenBlastG) proteins")
    parser.add_argument("--gff_fun", required=True, help="GFF3 of funannotate protein alignments")
    parser.add_argument("--gff_anno", required=True, help="GFF3 of anno (GenBlastG) protein alignments")
    args = parser.parse_args()

    busco_fun = parse_busco_table(args.busco_fun)
    busco_anno = parse_busco_table(args.busco_anno)

    analyze_missing_buscos(busco_fun, busco_anno, args.gff_fun, args.gff_anno)
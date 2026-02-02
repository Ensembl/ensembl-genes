#!/usr/bin/env python3
"""
BUSCO Missing Diagnostic Analysis (Protein ID Mapping)
Identify which conserved BUSCOs are missing in an annotation and why,
by mapping BUSCO protein IDs to annotation GFF features.
"""

import argparse
from collections import defaultdict
from pathlib import Path
import pandas as pd

# -----------------------------
# Helper Functions
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
            score = float(parts[3]) if len(parts) > 3 else None
            length = int(parts[4]) if len(parts) > 4 else None

            buscos[bid] = {
                'status': status,
                'sequence': sequence,
                'score': score,
                'length': length,
                'chrom': None,
                'start': None,
                'end': None,
                'strand': None
            }
    return buscos


def parse_gff(gff_file):
    """Parse GFF3 into dictionary keyed by protein/transcript ID"""
    gff_data = defaultdict(list)
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)

            # Extract ID or Name from attributes
            feature_id = None
            for key in ['ID', 'Name', 'protein_id', 'transcript_id']:
                if f'{key}=' in attributes:
                    for attr in attributes.split(';'):
                        if attr.startswith(f'{key}='):
                            feature_id = attr.split('=')[1]
                            break
                    if feature_id:
                        break
            if feature_id:
                gff_data[feature_id].append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'feature_type': feature
                })
    return gff_data


def load_genome_fasta(fasta_file):
    """Load genome sequences as dict"""
    sequences = {}
    current_seq = None
    current_data = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    sequences[current_seq] = ''.join(current_data).upper()
                current_seq = line[1:].strip().split()[0]
                current_data = []
            else:
                current_data.append(line.strip())
        if current_seq:
            sequences[current_seq] = ''.join(current_data).upper()
    return sequences


def region_stats(sequence, start, end):
    """Calculate basic stats for a genomic region"""
    if not sequence or start is None or end is None:
        return None, None
    region_seq = sequence[start-1:end]  # GFF is 1-based
    length = len(region_seq)
    n_fraction = region_seq.count('N') / length * 100 if length > 0 else 0
    gc_count = region_seq.count('G') + region_seq.count('C')
    atgc_count = sum(region_seq.count(b) for b in 'ATGC')
    gc_content = (gc_count / atgc_count * 100) if atgc_count > 0 else None
    return gc_content, n_fraction


def map_buscos_to_gff(buscos, gff_data):
    """Fill in chromosome/start/end/strand from GFF based on Sequence ID"""
    for bid, info in buscos.items():
        seq_id = info['sequence']
        if seq_id and seq_id in gff_data:
            # Use the first feature if multiple
            feat = gff_data[seq_id][0]
            info['chrom'] = feat['chrom']
            info['start'] = feat['start']
            info['end'] = feat['end']
            info['strand'] = feat['strand']
    return buscos


def analyze_missing_buscos(buscos_ref, buscos_test, gff_test_data, genome_seq=None, name_ref="Funannotate", name_test="EnsemblAnno"):
    """Analyze missing BUSCOs and generate diagnostics"""
    complete_ref = {bid for bid, info in buscos_ref.items() if info.get('status') in ['Complete', 'Duplicated']}
    complete_test = {bid for bid, info in buscos_test.items() if info.get('status') in ['Complete', 'Duplicated']}
    missing_buscos = complete_ref - complete_test

    print(f"Found {len(missing_buscos)} BUSCOs complete in {name_ref} but missing in {name_test}")

    diagnostics = []

    for bid in sorted(missing_buscos):
        info = buscos_ref[bid]

        # Safely get keys
        chrom = info.get('chrom')
        start = info.get('start')
        end = info.get('end')
        strand = info.get('strand')

        # Check if any gene model exists in test annotation
        feature_present = False
        seq_id = info.get('sequence')
        if seq_id and seq_id in gff_test_data:
            feature_present = True

        # Genome stats
        gc_content, n_fraction = None, None
        if genome_seq and chrom in genome_seq and start and end:
            gc_content, n_fraction = region_stats(genome_seq[chrom], start, end)

        # Determine likely reason
        reasons = []
        if not feature_present:
            reasons.append("No gene model in test annotation")
        if n_fraction is not None and n_fraction > 20:
            reasons.append(f"High N content ({n_fraction:.1f}%)")
        if gc_content is not None and (gc_content < 30 or gc_content > 70):
            reasons.append(f"Extreme GC ({gc_content:.1f}%)")
        if not reasons:
            reasons.append("Unknown / minor differences")

        diagnostics.append({
            'BUSCO_ID': bid,
            'Protein_ID': seq_id,
            'Chromosome': chrom,
            'Start': start,
            'End': end,
            'Strand': strand,
            'AnnotationFeaturePresent': feature_present,
            'GC_Content': gc_content,
            'N_Content': n_fraction,
            'LikelyReason': "; ".join(reasons)
        })

    df = pd.DataFrame(diagnostics)
    df.to_csv(f'missing_busco_diagnostics_{name_test}.csv', index=False)
    print(f"\n✓ Diagnostic table saved to missing_busco_diagnostics_{name_test}.csv")


# -----------------------------
# Command-line Interface
# -----------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagnose missing BUSCOs by mapping protein IDs to annotation GFFs")
    parser.add_argument("--busco_ref", required=True, help="Reference BUSCO full_table.tsv (Funannotate)")
    parser.add_argument("--busco_test", required=True, help="Test BUSCO full_table.tsv (EnsemblAnno)")
    parser.add_argument("--gff_ref", required=True, help="GFF3 of reference annotation (Funannotate)")
    parser.add_argument("--gff_test", required=True, help="GFF3 of test annotation (EnsemblAnno)")
    parser.add_argument("--genome", help="Genome FASTA (optional for GC/N analysis)")
    parser.add_argument("--name_ref", default="Funannotate", help="Reference annotation name")
    parser.add_argument("--name_test", default="EnsemblAnno", help="Test annotation name")

    args = parser.parse_args()

    print("Parsing BUSCO tables...")
    buscos_ref = parse_busco_table(args.busco_ref)
    buscos_test = parse_busco_table(args.busco_test)

    print("Parsing GFFs...")
    gff_ref_data = parse_gff(args.gff_ref)
    gff_test_data = parse_gff(args.gff_test)

    print("Mapping BUSCOs to reference GFF coordinates...")
    buscos_ref = map_buscos_to_gff(buscos_ref, gff_ref_data)

    genome_seq = None
    if args.genome:
        print("Loading genome FASTA...")
        genome_seq = load_genome_fasta(args.genome)

    analyze_missing_buscos(buscos_ref, buscos_test, gff_test_data, genome_seq, args.name_ref, args.name_test)
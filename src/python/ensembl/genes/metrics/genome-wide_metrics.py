#!/usr/bin/env python3

import argparse
import sys
import tempfile
import subprocess
import os
import datetime

# --- Parsers and utility --- #

def genome_size_from_fasta(fasta_path):
    """Calculate the genome size from a soft-masked FASTA file.
    This function reads the FASTA file, counts the total number of bases,
    and returns the genome size as an integer. Soft-masked bases (lowercase) are included
    in the count.

    Args:
        fasta_path (str): Path to the input FASTA file.
    Returns:
        int: The total genome size in base pairs.
    """
    genome_size = 0
    with open(fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            if not line.startswith('>'):
                genome_size += len(line.strip())
    return genome_size

def parse_gff(gff_path, feature_type):
    """Parse the GFF file and extract gene coordinates.
    This function reads the GFF file, extracts gene features, and returns a list of gene annotations.
    Alternatively, if feature_type is specified, it will extract features of that type.

    Args:
        gff_path (str): Path to the input GFF file.
        feature_type (str): The type of feature to extract (default: gene).
    Returns:
        list: A list of gene coordinates.
    """
    features = []
    with open (gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 7 or cols[2] != feature_type:
                continue
            chrom = cols[0]
            start = int(cols[3]) - 1  # Convert to 0-based
            end = int(cols[4])  # GFF is 1-based, end is inclusive
            strand = cols[6]
            features.append((chrom, start, end, strand))
    return features

def write_bed(features, bed_path, feature_type):
    """Write coordinates features to a BED file.
    This function takes the list of coordinates features and writes them to a BED file with the format:
    chrom, start, end, name, score, strand. The score field is set to '.'.

    Args:
        features (list): A list of coordinate features.
        bed_path (str): Path to the output BED file.
        feature_type (str): The type of feature to write.
    """
    features_sorted = sorted(features, key=lambda x: (x[0], x[1]))
    with open(bed_path, 'w') as bed_file:
        for i, (chrom, start, end, strand) in enumerate(features_sorted):
            name = f"{feature_type}_{i+1}"
            bed_file.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")

def run(cmd):
    """Run a shell command and return stdout as string."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout

# --- Metrics calculations --- #

def gene_density(features, genome_size):
    """Calculate gene density as the number of genes per megabase (Mb).
    This function takes the list of gene features and the genome size, and calculates
    the gene density using the formula: (number of genes / genome size in Mb).

    Args:
        features (list): A list of coordinate features.
        genome_size (int): The total genome size in base pairs.
    Returns:
        float: The gene density in genes per megabase.
    """
    num_features = len(features)
    return (num_features / genome_size) * 1e6

def pct_genome_covered(bed_path, genome_size):
    """Calculate the percentage of the genome covered by features.
    This function uses bedtools to calculate the total length of the genome covered by the features
    in the BED file, and then calculates the percentage of the genome that is covered.

    Args:
        bed_path (str): Path to the input BED file containing feature coordinates.
        genome_size (int): The total genome size in base pairs.
    Returns:
        float: The percentage of the genome covered by features.
    """
    merged = run(f"bedtools merge -i {bed_path}")
    total = 0
    for line in merged.strip().split('\n'):
        if line:
            cols = line.split('\t')
            total += int(cols[2]) - int(cols[1])
    return (total / genome_size) * 100

def pct_masked_annotated(bed_path, fasta_path):
    """Calculate the percentage of masked regions that are annotated by features.
    This function uses bedtools to calculate the total length of the genome that is soft-masked,
    and the total length of the genome that is both soft-masked and annotated by features.
    It then calculates the percentage of masked regions that are annotated.

    Args:
        bed_path (str): Path to the input BED file containing feature coordinates.
        fasta_path (str): Path to the input FASTA file containing the soft-masked genome.
    Returns:
        tuple: A tuple containing the percentage of masked regions that are annotated and the total length of masked regions.
    """
    masked_bed = tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w')
    masked_bed_path = masked_bed.name
    masked_bed.close()

    chrom = None
    pos = 0
    in_mask = False
    mask_start = 0
    total_masked = 0

    with open(fasta_path) as fa, open(masked_bed_path, 'w') as out:
        for line in fa:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if in_mask:
                    out.write(f"{chrom}\t{mask_start}\t{pos}\n")
                    total_masked += pos - mask_start
                chrom = line[1:].split()[0]
                pos = 0
                in_mask = False
            else:
                for base in line.strip():
                    if base.islower() and not in_mask:
                        mask_start = pos
                        in_mask = True
                    elif not base.islower() and in_mask:
                        out.write(f"{chrom}\t{mask_start}\t{pos}\n")
                        total_masked += pos - mask_start
                        in_mask = False
                    pos += 1
        if in_mask:
            out.write(f"{chrom}\t{mask_start}\t{pos}\n")
            total_masked += pos - mask_start
    print(datetime.datetime.now())

    if total_masked == 0:
        os.unlink(masked_bed_path)
        return None, None

    intersect = run(
        f"bedtools intersect -a {masked_bed_path} -b {bed_path} -u | "
        f"awk '{{sum += $3 - $2}} END {{print sum}}'"
    )
    os.unlink(masked_bed_path)

    masked_annotated = int(intersect.strip()) if intersect.strip() else 0
    print(datetime.datetime.now())
    return (masked_annotated / total_masked) * 100, total_masked

def overlapping_features(bed_path):
    """Calculate the number of overlapping feature pairs on the same strand.
    This function uses bedtools to find overlapping features in the BED file that are on the same strand,
    and counts the number of unique pairs of overlapping features.

    Args:
        bed_path (str): Path to the input BED file containing feature coordinates.
    Returns:
        int: The number of overlapping feature pairs on the same strand.
    """
    out = run(
        f"bedtools intersect -a {bed_path} -b {bed_path} -s -wa -wb | "
        f"awk '$4 != $10'"  # filter self-hits by name
    )
    lines = [l for l in out.strip().split('\n') if l]
    # Each pair appears twice (A vs B and B vs A), divide by 2
    return len(lines) // 2

def antisense_features(bed_path):
    """Calculate the number of antisense feature pairs on opposite strands.
    This function uses bedtools to find overlapping features in the BED file that are on opposite strands,
    and counts the number of unique pairs of antisense features.

    Args:
        bed_path (str): Path to the input BED file containing feature coordinates.
    Returns:
        int: The number of antisense feature pairs on opposite strands.
    """
    out = run(
        f"bedtools intersect -a {bed_path} -b {bed_path} -S -wa -wb | "
        f"awk '$4 != $10'"  # filter self-hits by name
    )
    lines = [l for l in out.strip().split('\n') if l]
    return len(lines) // 2

def main():
    parser = argparse.ArgumentParser(description="Calculate genome-wide metrics for Ensembl features.")
    parser.add_argument('gff', help='Path to the input GFF file containing gene annotations')
    parser.add_argument('--genome-size', type=int, help='Size of the genome in base pairs. Not needed if --fasta is provided')
    parser.add_argument('--fasta', help='Path to the soft-masked genome FASTA file.')
    parser.add_argument('--feature-type', default='gene', help='Feature type to consider for gene counting (default: gene)')
    args = parser.parse_args()

    if args.genome_size:
        genome_size = args.genome_size
    elif args.fasta:
        genome_size = genome_size_from_fasta(args.fasta)
    else:
        parser.error('You must provide either --genome-size or --fasta to determine the genome size.')

    features = parse_gff(args.gff, args.feature_type)
    if not features:
        print("[ERROR] No gene features found in GFF3.", file=sys.stderr)
        sys.exit(1)

    tmp_bed = tempfile.NamedTemporaryFile(suffix='.bed', delete=False)
    tmp_bed_path = tmp_bed.name
    tmp_bed.close()
    write_bed(features, tmp_bed_path, args.feature_type)

    print(f"{'Metric':<45} {'Value'}")
    print('-' * 65)
    print(f"{'Genome size (bp)':<45} {genome_size}")
    print(f"{f'{args.feature_type} count':<45} {len(features)}")
    print(f"{f'{args.feature_type} density ({args.feature_type}s/Mb)':<45} {gene_density(features, genome_size):.4f}")
    print(f"{'% genome covered by ' + args.feature_type:<45} {pct_genome_covered(tmp_bed_path, genome_size):.4f}")
    if args.fasta:
        pct, total_masked = pct_masked_annotated(tmp_bed_path, args.fasta)
        if pct is None:
            print(f"{'% masked regions annotated':<45} N/A (no soft-masking detected)")
        else:
            print(f"{'Total soft-masked bases':<45} {total_masked} ({(total_masked/genome_size)*100:.2f}%)")
            print(f"{'% masked regions annotated':<45} {pct:.4f}")
    else:
        print(f"{'% masked regions annotated':<45} skipped (no --fasta provided)")
    print(f"{f'Overlapping {args.feature_type} pairs (same strand)':<45} {overlapping_features(tmp_bed_path)}")
    print(f"{f'Antisense {args.feature_type} pairs (opposite strand)':<45} {antisense_features(tmp_bed_path)}")

if __name__ == "__main__":
    main()
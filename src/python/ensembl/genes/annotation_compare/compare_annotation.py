"""
Genome Annotation Comparison Tool
Compares two genome annotations (GFF3/GTF format) and identifies differences
"""
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse


def run_gffcompare(
    gff_ref: Path,
    gff_test: Path,
    gffcompare_bin: Path,
    output_path: Path,
    ref_name: str,
    test_name: str,
) -> None:
    """Run gffcompare on two GFF/GTF files."""

    output_path.mkdir(parents=True, exist_ok=True)

    # Run reference vs test

    prefix1 = output_path / f"{ref_name}_vs_{test_name}"

    cmd1 = [
        str(gffcompare_bin),
        "-r", str(gff_ref),
        str(gff_test),
        "-o", str(prefix1),
    ]

    subprocess.run(cmd1, check=True)

    # Run test vs reference

    prefix2 = output_path / f"{test_name}_vs_{ref_name}"

    cmd2 = [
	    str(gffcompare_bin),
	    "-r", str(gff_ref),
	    str(gff_test),
	    "-o", str(prefix2),
    ]

    subprocess.run(cmd2, check=True)

def parse_annotation(annotation_file: Path):
    """Dispatch parser based on file extension."""
    ext = annotation_file.suffix.lower()

    if ext in [".gff", ".gff3"]:
        return parse_gff(annotation_file)
    elif ext == ".gtf":
        return parse_gtf(annotation_file)
    else:
        raise ValueError(f"Unsupported file type: {ext}")

def parse_gtf(gtf_file):
    """Parse GTF and build structure without requiring CDS."""
    genes = {}
    transcripts = defaultdict(list)
    exons = defaultdict(list)
    cds = defaultdict(list)

    gene_bounds = defaultdict(lambda: [10**12, 0])  # min, max

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)

            # Parse GTF attributes
            attr_dict = {}
            for attr in attributes.split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                key, val = attr.split(" ", 1)
                attr_dict[key] = val.strip('"')

            gene_id = attr_dict.get("gene_id")
            transcript_id = attr_dict.get("transcript_id")

            if not gene_id:
                continue

            # Track gene bounds from all features
            gene_bounds[gene_id][0] = min(gene_bounds[gene_id][0], start)
            gene_bounds[gene_id][1] = max(gene_bounds[gene_id][1], end)

            if feature == "transcript":
                transcripts[gene_id].append({
                    "id": transcript_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": end - start + 1
                })

            elif feature == "exon":
                exons[gene_id].append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": end - start + 1
                })

            elif feature == "CDS":
                cds[gene_id].append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": end - start + 1
                })

    # Create gene entries from collected bounds
    for gene_id, (start, end) in gene_bounds.items():
        genes[gene_id] = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "length": end - start + 1
        }

    return genes, transcripts, exons, cds


def parse_gff(gff_file):
	"""Parse GFF3/GTF file and extract gene/transcript information"""
	genes = {}
	transcripts = defaultdict(list)
	exons = defaultdict(list)
	cds = defaultdict(list)

	with open(gff_file, 'r') as f:
		for line in f:
			if line.startswith('#'):
				continue

			parts = line.strip().split('\t')
			if len(parts) < 9:
				continue

			chrom, source, feature, start, end, score, strand, phase, attributes = parts
			start, end = int(start), int(end)

			# Parse attributes
			attr_dict = {}
			for attr in attributes.split(';'):
				if '=' in attr:
					key, val = attr.split('=', 1)
					attr_dict[key] = val
				elif ' ' in attr:  # GTF format
					parts_attr = attr.strip().split(' ', 1)
					if len(parts_attr) == 2:
						key, val = parts_attr
						attr_dict[key] = val.strip('"')

			# Get gene/transcript IDs
			gene_id = attr_dict.get('gene_id') or attr_dict.get('ID') or attr_dict.get('Parent', '').split(',')[0]
			transcript_id = attr_dict.get('transcript_id') or attr_dict.get('ID')
			parent = attr_dict.get('Parent', '').split(',')[0] if 'Parent' in attr_dict else gene_id

			if feature in ['gene', 'Gene']:
				genes[gene_id] = {
					'chrom': chrom,
					'start': start,
					'end': end,
					'strand': strand,
					'length': end - start + 1
				}
			elif feature in ['mRNA', 'transcript', 'Transcript']:
				transcripts[parent].append({
					'id': transcript_id,
					'chrom': chrom,
					'start': start,
					'end': end,
					'strand': strand,
					'length': end - start + 1
				})
			elif feature in ['exon', 'Exon']:
				exons[parent].append({
					'chrom': chrom,
					'start': start,
					'end': end,
					'strand': strand,
					'length': end - start + 1
				})
			elif feature == 'CDS':
				cds[parent].append({
					'chrom': chrom,
					'start': start,
					'end': end,
					'strand': strand,
					'length': end - start + 1
				})

	return genes, transcripts, exons, cds


def calculate_statistics(genes, transcripts, exons, cds):
	"""Calculate summary statistics for an annotation"""
	stats = {}

	stats['total_genes'] = len(genes)
	stats['total_transcripts'] = sum(len(t) for t in transcripts.values())
	stats['total_exons'] = sum(len(e) for e in exons.values())
	stats['total_cds'] = sum(len(c) for c in cds.values())

	# Gene lengths
	gene_lengths = [g['length'] for g in genes.values()]
	if gene_lengths:
		stats['mean_gene_length'] = np.mean(gene_lengths)
		stats['median_gene_length'] = np.median(gene_lengths)

	# Exons per gene
	exons_per_gene = [len(e) for e in exons.values() if len(e) > 0]
	if exons_per_gene:
		stats['mean_exons_per_gene'] = np.mean(exons_per_gene)
		stats['median_exons_per_gene'] = np.median(exons_per_gene)

	# Transcripts per gene
	transcripts_per_gene = [len(t) for t in transcripts.values() if len(t) > 0]
	if transcripts_per_gene:
		stats['mean_transcripts_per_gene'] = np.mean(transcripts_per_gene)

	# CDS lengths
	cds_lengths = []
	for gene_cds in cds.values():
		total_cds_len = sum(c['length'] for c in gene_cds)
		if total_cds_len > 0:
			cds_lengths.append(total_cds_len)

	if cds_lengths:
		stats['mean_cds_length'] = np.mean(cds_lengths)
		stats['median_cds_length'] = np.median(cds_lengths)

	return stats


def compare_annotations(file1, file2, name1="Funannotate", name2="EnsemblAnno"):
	"""Main comparison function"""

	print(f"Parsing {name1} annotation...")
	genes1, trans1, exons1, cds1 = parse_annotation(Path(file1))

	print(f"Parsing {name2} annotation...")
	genes2, trans2, exons2, cds2 = parse_annotation(Path(file2))

	print("\n" + "=" * 60)
	print("ANNOTATION STATISTICS")
	print("=" * 60)

	stats1 = calculate_statistics(genes1, trans1, exons1, cds1)
	stats2 = calculate_statistics(genes2, trans2, exons2, cds2)

	# Print comparison table
	print(f"\n{'Metric':<30} {name1:>15} {name2:>15}")
	print("-" * 60)
	for key in stats1:
		val1 = stats1.get(key, 0)
		val2 = stats2.get(key, 0)
		if isinstance(val1, float):
			print(f"{key:<30} {val1:>15.1f} {val2:>15.1f}")
		else:
			print(f"{key:<30} {val1:>15} {val2:>15}")

	# Identify unique and shared genes by genomic position
	print("\n" + "=" * 60)
	print("GENOMIC OVERLAP ANALYSIS")
	print("=" * 60)

	# Create position-based gene sets
	def get_gene_positions(genes):
		return {(g['chrom'], g['start'], g['end'], g['strand']) for g in genes.values()}

	pos1 = get_gene_positions(genes1)
	pos2 = get_gene_positions(genes2)

	unique1 = pos1 - pos2
	unique2 = pos2 - pos1
	shared = pos1 & pos2

	print(f"\nGenes unique to {name1}: {len(unique1)}")
	print(f"Genes unique to {name2}: {len(unique2)}")
	print(f"Genes in both (same position): {len(shared)}")
	print(f"\nPercentage of {name1} genes missing in {name2}: {len(unique1) / len(pos1) * 100:.1f}%")
	print(f"Percentage of {name2} genes not in {name1}: {len(unique2) / len(pos2) * 100:.1f}%")

	# Create visualizations
	create_plots(stats1, stats2, genes1, genes2, name1, name2)

	return stats1, stats2, unique1, unique2, shared


def create_plots(stats1, stats2, genes1, genes2, name1, name2):
	"""Create comparison plots"""

	fig, axes = plt.subplots(2, 2, figsize=(14, 10))

	# Plot 1: Gene count comparison
	metrics = ['total_genes', 'total_transcripts', 'total_exons', 'total_cds']
	labels = ['Genes', 'Transcripts', 'Exons', 'CDS']
	values1 = [stats1.get(m, 0) for m in metrics]
	values2 = [stats2.get(m, 0) for m in metrics]

	x = np.arange(len(labels))
	width = 0.35

	axes[0, 0].bar(x - width / 2, values1, width, label=name1, color='steelblue')
	axes[0, 0].bar(x + width / 2, values2, width, label=name2, color='coral')
	axes[0, 0].set_ylabel('Count')
	axes[0, 0].set_title('Feature Counts Comparison')
	axes[0, 0].set_xticks(x)
	axes[0, 0].set_xticklabels(labels)
	axes[0, 0].legend()
	axes[0, 0].grid(axis='y', alpha=0.3)

	# Plot 2: Gene length distributions
	lengths1 = [g['length'] for g in genes1.values()]
	lengths2 = [g['length'] for g in genes2.values()]

	axes[0, 1].hist(lengths1, bins=50, alpha=0.5, label=name1, color='steelblue', density=True)
	axes[0, 1].hist(lengths2, bins=50, alpha=0.5, label=name2, color='coral', density=True)
	axes[0, 1].set_xlabel('Gene Length (bp)')
	axes[0, 1].set_ylabel('Density')
	axes[0, 1].set_title('Gene Length Distribution')
	axes[0, 1].legend()
	axes[0, 1].set_xlim(0, np.percentile(lengths1 + lengths2, 95))

	# Plot 3: Average metrics comparison
	avg_metrics = ['mean_gene_length', 'mean_exons_per_gene', 'mean_cds_length']
	avg_labels = ['Avg Gene\nLength', 'Avg Exons\nper Gene', 'Avg CDS\nLength']
	avg_values1 = [stats1.get(m, 0) for m in avg_metrics]
	avg_values2 = [stats2.get(m, 0) for m in avg_metrics]

	x = np.arange(len(avg_labels))
	axes[1, 0].bar(x - width / 2, avg_values1, width, label=name1, color='steelblue')
	axes[1, 0].bar(x + width / 2, avg_values2, width, label=name2, color='coral')
	axes[1, 0].set_ylabel('Value')
	axes[1, 0].set_title('Average Metrics Comparison')
	axes[1, 0].set_xticks(x)
	axes[1, 0].set_xticklabels(avg_labels)
	axes[1, 0].legend()
	axes[1, 0].grid(axis='y', alpha=0.3)

	# Plot 4: Chromosome distribution
	chrom_counts1 = defaultdict(int)
	chrom_counts2 = defaultdict(int)

	for g in genes1.values():
		chrom_counts1[g['chrom']] += 1
	for g in genes2.values():
		chrom_counts2[g['chrom']] += 1

	# Get common chromosomes
	all_chroms = sorted(set(list(chrom_counts1.keys()) + list(chrom_counts2.keys())))[:20]  # Top 20

	vals1 = [chrom_counts1[c] for c in all_chroms]
	vals2 = [chrom_counts2[c] for c in all_chroms]

	x = np.arange(len(all_chroms))
	axes[1, 1].bar(x - width / 2, vals1, width, label=name1, color='steelblue')
	axes[1, 1].bar(x + width / 2, vals2, width, label=name2, color='coral')
	axes[1, 1].set_ylabel('Gene Count')
	axes[1, 1].set_title('Genes per Chromosome (Top 20)')
	axes[1, 1].set_xticks(x)
	axes[1, 1].set_xticklabels(all_chroms, rotation=45, ha='right')
	axes[1, 1].legend()
	axes[1, 1].grid(axis='y', alpha=0.3)

	plt.tight_layout()
	plt.savefig('annotation_comparison.png', dpi=300, bbox_inches='tight')
	print("\n✓ Plots saved to 'annotation_comparison.png'")

def main(gff_ref= Path, gff_test= Path, ref_name = str, test_name = str, gffcompare_bin = Path, output_path = Path):

	run_gffcompare(Path(gff_ref), Path(gff_test), Path(gffcompare_bin), Path(output_path), str(ref_name), str(test_name))

	compare_annotations(gff_ref, gff_test, str(ref_name), str(test_name))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Compare two genome annotations for different tools."
	)

	parser.add_argument(
		"--gff_ref",
		required=True,
		help="Path to the reference GFF3 file"
	)
	parser.add_argument(
		"--gff_test",
		required=True,
		help="Path to the test GFF3 file"
	)
	parser.add_argument(
		"--ref_name",
		required=True,
		help="Name of first annotation tool"
	)
	parser.add_argument(
		"--test_name",
		required=True,
		help="Name of second annotation tool"
	)
	parser.add_argument(
		"--gffcompare_bin",
		required=True,
		help="Path tho gffcompare"
	)
	parser.add_argument(
		"--output_path",
		required=True,
		help="Path where to save the results"
	)

	args = parser.parse_args()

	main(args.gff_ref, args.gff_test, args.ref_name, args.test_name, args.gffcompare_bin, args.output_path)
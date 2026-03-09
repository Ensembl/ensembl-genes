"""
Script to parse GFF3 and GTF files using PyRanges1.

Usage:
	python parse_annotation.py --file_path annotations.gff3
"""

import argparse
from pathlib import Path
import pyranges1 as pr


def parse_gff3(file_path: str):
	"""
	Parse a GFF3 file using PyRanges1.
	Args:
		file_path: Path to the GFF3 file
	Returns:
		PyRanges object
	"""
	gff3_file = pr.read_gff3(file_path)
	print(gff3_file)

	return gff3_file


def parse_gtf(file_path: str):
	"""
	Parse a GTF file using PyRanges1.
	Args:
		file_path: Path to the GTF file
	Returns:
		PyRanges object
	"""

	gtf_file = pr.read_gtf(file_path)
	print(gtf_file)

	return gtf_file


def main(file_path: str):
	"""
	Parse a GTF or GFF3 file using PyRanges1. Main function and entry to the script
	Args:
		file_path: Path to the GTF file
	Returns:
		PyRanges object
	"""

	path = Path(file_path)
	suffix = path.suffix.lower()

	if suffix == ".gff3":
		print(f"Parsing GFF3... ({file_path})")
		data = parse_gff3(file_path)
	elif suffix == ".gtf":
		print(f"Parsing GTF... ({file_path})")
		data = parse_gtf(file_path)
	else:
		raise ValueError("Unsupported file type. Use .gff3 or .gtf")

	return data


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Parse GFF3 or GTF with PyRanges1")
	parser.add_argument("--file_path", required=True, help="Path to input annotation file")

	args = parser.parse_args()
	main(args.file_path)

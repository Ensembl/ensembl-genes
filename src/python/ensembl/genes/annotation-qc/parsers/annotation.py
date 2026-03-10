
import pyranges1 as pr
from pathlib import Path


def parse_gff3(file_path: str):
	"""
	Parse a GFF3 file using PyRanges1.
	Args:
		file_path: Path to the GFF3 file
	Returns:
		PyRanges object
	"""
	gff3_file = pr.read_gff3(file_path)
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
	return gtf_file


def parse_annotation(file_path: str):
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

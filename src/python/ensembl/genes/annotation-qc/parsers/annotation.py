"""
Script to parse GFF3 and GTF files using PyRanges1.

Usage:
	python parse_annotation.py --file_path annotations.gff3
"""

import argparse
from pathlib import Path
import pyranges1 as pr
import re

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


def fix_gtf_attributes(attr: str) -> str:
    """Fix semicolons inside quoted attributes for GTF."""
    # Replace semicolons inside quotes with commas
    return re.sub(r'\"([^\"]*);([^\"]*)\"', r'"\1,\2"', attr).strip()

def parse_gtf(file_path: str) -> pr.PyRanges:
    path = Path(file_path)
    fixed_file = path.parent / f"fixed_{path.name}"

    skipped_lines = 0
    with open(file_path) as infile, open(fixed_file, "w") as outfile:
        for i, line in enumerate(infile, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                print(f"Line {i} skipped (not 9 columns): {line}")
                skipped_lines += 1
                continue

            # Strip whitespace in the attribute column
            attr = parts[8].strip()
            if not attr:
                print(f"Line {i} skipped (empty Attribute column): {line}")
                skipped_lines += 1
                continue

            # Fix semicolons inside quotes
            parts[8] = fix_gtf_attributes(attr)

            # Validate key-value pairs
            try:
                kvs = [kv.strip() for kv in parts[8].split(";") if kv.strip()]
                for kv in kvs:
                    k, v = kv.split(None, 1)
            except Exception as e:
                print(f"Line {i} skipped (invalid key-value pairs): {line}")
                print(f"  Error: {e}")
                skipped_lines += 1
                continue

            outfile.write("\t".join(parts) + "\n")

    print(f"Finished fixing GTF. Skipped {skipped_lines} invalid lines.")
    gtf_file = pr.read_gtf(str(fixed_file), ignore_bad=False)
    print(f"Loaded fixed GTF: {fixed_file}")
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

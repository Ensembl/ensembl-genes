from pyfaidx import Fasta


def parse_fasta(file_path: str) -> Fasta:
	"""
	Parse a FASTA file using pyfaidx.
	Args:
		file_path: Path to the FASTA file
	Returns:
		Fasta object (supports indexed access by sequence name and slicing)
	"""
	return Fasta(file_path)

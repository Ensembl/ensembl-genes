
"""
Parse FASTA files using pyfaidx.

The main entry point is parse_fasta, which returns a Fasta object that supports
indexed access by sequence name and slicing.
"""


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

"""
Parse FASTA files using pyfaidx.

The main entry point is parse_fasta, which returns a Fasta object that supports
indexed access by sequence name and slicing.
"""

import gzip

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


def parse_sequence_lengths(file_path: str) -> dict[str, int]:
    """
    Return the length of every sequence in a FASTA file.

    Streams the file (plain or gzip) instead of using pyfaidx, so no .fai index
    is written beside the genome. Names are the first word of each header, as
    in pyfaidx.
    Args:
            file_path: Path to the FASTA file
    Returns:
            dict mapping sequence name to length
    """
    with open(file_path, "rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    lengths: dict[str, int] = {}
    name = None
    with gzip.open(file_path, "rt") if is_gzip else open(file_path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split(maxsplit=1)[0]
                if name in lengths:
                    raise ValueError(
                        f"Duplicate FASTA sequence name '{name}' in {file_path}"
                    )
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths

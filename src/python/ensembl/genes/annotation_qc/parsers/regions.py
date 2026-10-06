"""
Parse genomic region strings and region files.

Regions are written as 'seqname', 'seqname:start-end' or 'seqname:pos', with
one-based inclusive coordinates (samtools convention). They are returned as
zero-based half-open intervals so they can be compared directly with the
comparison schema; start/end are None for a whole sequence.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Region:
    """A sequence or a zero-based half-open interval on it."""

    seqname: str
    start: int | None = None
    end: int | None = None

    def __str__(self) -> str:
        if self.start is None:
            return self.seqname
        return f"{self.seqname}:{self.start + 1}-{self.end}"


def parse_region(text: str) -> Region:
    """
    Parse 'seqname', 'seqname:start-end' or 'seqname:pos' (one-based inclusive).
    Args:
            text: Region string
    Returns:
            Region with zero-based half-open coordinates
    """
    text = text.strip()
    seqname, sep, coords = text.rpartition(":")
    if not sep or not coords.replace(",", "").replace("-", "").isdigit():
        return Region(text)
    coords = coords.replace(",", "")
    first, _, last = coords.partition("-")
    start, end = int(first), int(last or first)
    if start < 1 or end < start:
        raise ValueError(f"Invalid region '{text}': expected 1 <= start <= end")
    return Region(seqname, start - 1, end)


def load_regions_file(file_path: str) -> list[Region]:
    """
    Read one region per line; blank lines and '#' comments are ignored.
    Args:
            file_path: Path to the regions file
    Returns:
            list of Region
    """
    with open(file_path, "rt") as handle:
        return [
            parse_region(line)
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]

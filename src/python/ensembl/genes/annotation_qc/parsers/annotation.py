"""
Parse annotation files using PyRanges1.

The main entry point is parse_annotation, which detects whether the input file
is GFF3 or GTF and calls the appropriate parser.

parse_annotation_for_comparison is a separate adapter for pairwise comparison.
It detects the format from file content (including the Tiberius GFF/GTF hybrid)
and returns the normalised comparison schema described in
parsers/annotation_normalise.py. It does not change what parse_annotation
returns.
"""

import gzip
import os
import re
import tempfile
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import pyranges1 as pr

from ensembl.genes.annotation_qc.parsers.annotation_normalise import (
    normalise_gff3,
    normalise_gtf,
    normalise_tiberius,
)


def is_gzip_file(file_path: str) -> bool:
    """True when the file content is gzip (or BGZF) compressed, whatever its name."""
    with open(file_path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


@contextmanager
def _reader_path(file_path: str):
    """
    Yield a path that PyRanges1 reads with the right decompression.

    PyRanges1 (through pandas) infers compression from the file suffix, so a
    gzip/BGZF file named e.g. *.gff3.bgz fails and a plain file named *.gz is
    misread. When content and suffix disagree, a temporary symbolic link whose
    name matches the content is read instead. The file itself is never copied
    or changed.
    """
    compressed = is_gzip_file(file_path)
    if compressed == str(file_path).lower().endswith(".gz"):
        yield file_path
        return
    with tempfile.TemporaryDirectory(prefix="annotation_qc_") as tmp:
        link = Path(tmp) / ("annotation.gz" if compressed else "annotation.txt")
        os.symlink(Path(file_path).resolve(), link)
        yield str(link)


def parse_gff3(file_path: str):
    """
    Parse a GFF3 file using PyRanges1.
    Args:
            file_path: Path to the GFF3 file, plain or gzip/BGZF (detected from content)
    Returns:
            PyRanges object
    """
    with _reader_path(file_path) as path:
        gff3_file = pr.read_gff3(path)
    return gff3_file


def parse_gtf(file_path: str, duplicate_attr: bool = False):
    """
    Parse a GTF file using PyRanges1.
    Args:
            file_path: Path to the GTF file, plain or gzip/BGZF (detected from content)
            duplicate_attr: Keep every value of a repeated attribute (e.g. several
                    'tag' entries) instead of only the last one
    Returns:
            PyRanges object
    """

    with _reader_path(file_path) as path:
        gtf_file = pr.read_gtf(path, duplicate_attr=duplicate_attr)
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
    suffixes = [s.lower() for s in path.suffixes]
    real_suffix = (
        suffixes[-2]
        if suffixes and suffixes[-1] == ".gz" and len(suffixes) >= 2
        else suffixes[-1]
    )

    if real_suffix == ".gff3":
        print(f"Parsing GFF3... ({file_path})")
        data = parse_gff3(file_path)
    elif real_suffix == ".gtf":
        print(f"Parsing GTF... ({file_path})")
        data = parse_gtf(file_path)
    else:
        raise ValueError("Unsupported file type. Use .gff3, .gtf, .gff3.gz, or .gtf.gz")

    return data


ANNOTATION_FORMATS = ("gff3", "gtf", "tiberius")

_GFF3_ATTR = re.compile(r"ID=|Parent=")
_GTF_ATTR = re.compile(r'(gene_id|transcript_id)\s+"')
_BARE_ATTR = re.compile(r"^\S+$")


def _open_text(file_path: str):
    if is_gzip_file(file_path):
        return gzip.open(file_path, "rt")
    return open(file_path, "rt", encoding="utf-8")


def detect_annotation_format(file_path: str, max_records: int = 200) -> str:
    """
    Detect GFF3, GTF or Tiberius hybrid format from file content.

    The Tiberius hybrid has bare identifiers in column 9 of gene/transcript rows
    and GTF-style attributes on exon/CDS rows. File extensions are not used.
    Args:
            file_path: Annotation path, optionally gzip-compressed
            max_records: Number of feature records to inspect
    Returns:
            One of "gff3", "gtf" or "tiberius"
    """
    gff3_score = gtf_score = bare_parents = gtf_children = records = 0
    with _open_text(file_path) as handle:
        for line in handle:
            if line.startswith("#"):
                if "gff-version" in line and "3" in line:
                    gff3_score += 5
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature, attributes = parts[2], parts[8].strip()
            records += 1
            is_gff3 = bool(_GFF3_ATTR.search(attributes))
            is_gtf = bool(_GTF_ATTR.search(attributes))
            gff3_score += is_gff3
            gtf_score += is_gtf
            if feature in ("gene", "transcript") and _BARE_ATTR.match(attributes):
                bare_parents += 1
            if feature in ("exon", "CDS") and is_gtf and not is_gff3:
                gtf_children += 1
            if records >= max_records:
                break

    if records == 0:
        raise ValueError(f"No GFF/GTF feature records found in {file_path}")
    if bare_parents >= 2 and gtf_children >= 2:
        return "tiberius"
    if gff3_score > gtf_score:
        return "gff3"
    if gtf_score > 0:
        return "gtf"
    return "gff3"


def _read_raw_columns(file_path: str) -> pd.DataFrame:
    """Read seqid, type, start, end and column 9 verbatim for every feature record."""
    records = []
    with _open_text(file_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            records.append((parts[0], parts[2], int(parts[3]), int(parts[4]), parts[8]))
    return pd.DataFrame(
        records, columns=["Chromosome", "Feature", "gff_start", "End", "attributes"]
    )


def parse_annotation_for_comparison(
    file_path: str, format_hint: str = "auto"
) -> pd.DataFrame:
    """
    Parse a GFF3, GTF or Tiberius hybrid file into the comparison schema.

    Coordinates are converted once, here, to zero-based half-open intervals
    (PyRanges convention: Start = GFF start - 1, End = GFF end). Genes,
    transcripts, exons and CDS are linked through explicit parentage; exons or
    CDS with several parents become one row per parent transcript. Unresolvable
    relationships raise ValueError. Parse counts are stored in
    ``df.attrs["parse_diagnostics"]``.
    Args:
            file_path: Annotation path, optionally gzip-compressed
            format_hint: "auto" (detect from content), "gff3", "gtf" or "tiberius"
    Returns:
            pandas DataFrame with COMPARISON_COLUMNS
    """
    fmt = detect_annotation_format(file_path) if format_hint == "auto" else format_hint
    if fmt not in ANNOTATION_FORMATS:
        raise ValueError(f"Unsupported annotation format '{fmt}'")

    if fmt == "gff3":
        features = normalise_gff3(pd.DataFrame(parse_gff3(file_path)))
    elif fmt == "gtf":
        features = normalise_gtf(
            pd.DataFrame(parse_gtf(file_path, duplicate_attr=True))
        )
    else:
        features = normalise_tiberius(
            pd.DataFrame(parse_gtf(file_path, duplicate_attr=True)),
            _read_raw_columns(file_path),
        )

    features.attrs["parse_diagnostics"]["format"] = fmt
    return features

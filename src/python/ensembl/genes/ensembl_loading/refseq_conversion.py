"""Converters for RefSeq FASTA and GFF3 files."""

from __future__ import annotations

import gzip
import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Collection, Iterator, TextIO

LOGGER = logging.getLogger(__name__)


@contextmanager
def open_text_maybe_gzip(path: str | Path) -> Iterator[TextIO]:
    """Open plain-text or gzip-compressed files for text reading."""

    input_path = Path(path)
    if input_path.suffix == ".gz":
        with gzip.open(input_path, "rt") as handle:
            yield handle
    else:
        with input_path.open("r") as handle:
            yield handle


def load_refseq_name_map(assembly_report_path: str | Path) -> dict[str, str]:
    """Load RefSeq accession to Ensembl-style seq_region name mappings."""

    accession_to_name: dict[str, str] = {}
    report_path = Path(assembly_report_path)
    with report_path.open("r") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue

            columns = line.strip().split("\t")
            if len(columns) < 7:
                continue

            sequence_name = columns[0]
            assigned_molecule = columns[2]
            refseq_accession = columns[6]
            if not refseq_accession or refseq_accession == "na":
                continue

            accession_to_name[refseq_accession] = (
                assigned_molecule if assigned_molecule != "na" else sequence_name
            )

    LOGGER.info(
        "Loaded %s sequence name mappings from %s",
        len(accession_to_name),
        report_path,
    )
    return accession_to_name


def default_gff_output_path(gff_path: str | Path) -> Path:
    """Return the default converted GFF3 output path."""

    path = Path(gff_path)
    name = path.name
    for suffix in (".gff.gz", ".gff3.gz", ".gff", ".gff3"):
        if name.endswith(suffix):
            return path.with_name(f"{name[: -len(suffix)]}_ensembl.gff3")
    return path.with_name(f"{path.stem}_ensembl.gff3")


def convert_fna_headers(
    fna_in_path: str | Path,
    assembly_report_path: str | Path,
    fna_out_path: str | Path,
    logger: logging.Logger | None = None,
) -> Path:
    """Convert RefSeq FASTA headers to Ensembl-style seq_region names.

    Parameters
    ----------
    fna_in_path
        Input genomic FASTA path.
    assembly_report_path
        NCBI assembly report containing RefSeq accession mappings.
    fna_out_path
        Output FASTA path with converted headers.
    logger
        Optional logger used for progress reporting.
    """

    log = logger or LOGGER
    accession_to_name = load_refseq_name_map(assembly_report_path)
    input_path = Path(fna_in_path)
    output_path = Path(fna_out_path)
    converted_headers = 0
    missing_headers = 0

    with (
        open_text_maybe_gzip(input_path) as input_handle,
        output_path.open("w") as output_handle,
    ):
        for line in input_handle:
            if line.startswith(">"):
                accession = line[1:].split()[0]
                name = accession_to_name.get(accession)
                if name is None:
                    name = accession
                    missing_headers += 1
                else:
                    converted_headers += 1
                output_handle.write(f">{name}\n")
            else:
                output_handle.write(line)

    log.info(
        "Wrote converted FASTA to %s (%s mapped headers, %s unmapped headers)",
        output_path,
        converted_headers,
        missing_headers,
    )
    return output_path


def convert_gff_to_ensembl(
    gff_path: str | Path,
    assembly_report_path: str | Path,
    output_path: str | Path | None = None,
    chrom_filter: Collection[str] | None = None,
    logger: logging.Logger | None = None,
) -> Path:
    """Convert a RefSeq-style GFF3 file to Ensembl-style seq_region names.

    Sequence IDs in feature rows and ``##sequence-region`` directives are
    rewritten using the NCBI assembly report. When ``chrom_filter`` is set, only
    matching original RefSeq accessions or converted names are written.
    """

    log = logger or LOGGER
    input_path = Path(gff_path)
    output = (
        Path(output_path)
        if output_path is not None
        else default_gff_output_path(input_path)
    )
    allowed_sequences = set(chrom_filter) if chrom_filter else None
    refseq_to_name = load_refseq_name_map(assembly_report_path)
    converted_features = 0
    skipped_features = 0
    unmapped_features = 0

    with (
        open_text_maybe_gzip(input_path) as input_handle,
        output.open("w") as output_handle,
    ):
        for line in input_handle:
            if line.startswith("#"):
                if line.lower().startswith("##sequence-region"):
                    parts = line.strip().split()
                    if len(parts) >= 3 and parts[1] in refseq_to_name:
                        parts[1] = refseq_to_name[parts[1]]
                        line = " ".join(parts) + "\n"
                output_handle.write(line)
                continue

            columns = line.rstrip().split("\t")
            if len(columns) < 9:
                skipped_features += 1
                continue

            seq_id = columns[0]
            mapped_seq_id = refseq_to_name.get(seq_id)
            if allowed_sequences is not None and not (
                seq_id in allowed_sequences or mapped_seq_id in allowed_sequences
            ):
                skipped_features += 1
                continue

            if mapped_seq_id is None:
                unmapped_features += 1
            else:
                columns[0] = mapped_seq_id
                converted_features += 1

            output_handle.write("\t".join(columns) + "\n")

    if unmapped_features:
        log.warning(
            "%s GFF3 features had no assembly-report mapping", unmapped_features
        )
    log.info(
        "Wrote converted GFF3 to %s (%s converted features, %s skipped)",
        output,
        converted_features,
        skipped_features,
    )
    return output

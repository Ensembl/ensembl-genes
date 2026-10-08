"""Converters for RefSeq FASTA and GFF3 files."""

from __future__ import annotations

import gzip
import logging
from collections.abc import Collection, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

# pylint: disable=too-many-locals


LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class _RepeatMaskerRecord:  # pylint: disable=too-many-instance-attributes
    """Relevant fields from one whitespace-delimited RepeatMasker row."""

    score: str
    sequence_name: str
    genomic_start: str
    genomic_end: str
    strand: str
    repeat_name: str
    repeat_class: str
    consensus_start: str
    consensus_end: str


def _parse_repeatmasker_row(fields: list[str]) -> _RepeatMaskerRecord | None:
    """Parse a RepeatMasker row and normalize its strand-specific coordinates."""

    if fields and fields[-1] == "*":
        fields = fields[:-1]
    if len(fields) < 15:
        return None

    (
        score,
        _percent_divergence,
        _percent_deletions,
        _percent_insertions,
        sequence_name,
        genomic_start,
        genomic_end,
        _query_left,
        strand_code,
        repeat_name,
        repeat_class,
        consensus_start_forward,
        consensus_end,
        consensus_start_reverse,
        _repeat_id,
        *_extra,
    ) = fields

    strand = "-" if strand_code == "C" else "+"
    consensus_start = (
        consensus_start_reverse if strand == "-" else consensus_start_forward
    )
    return _RepeatMaskerRecord(
        score=score,
        sequence_name=sequence_name,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
        strand=strand,
        repeat_name=repeat_name,
        repeat_class=repeat_class,
        consensus_start=consensus_start,
        consensus_end=consensus_end,
    )


def _get_repeat_type(repeat_class: str) -> str:
    """Resolve an Ensembl repeat type without importing ensembl-tools eagerly."""

    try:
        from ensembl.tools.anno.repeat_annotation.repeatmasker import (  # pylint: disable=import-outside-toplevel
            get_repeat_type,
        )
    except ImportError as error:  # pragma: no cover - environment dependent
        raise ImportError(
            "RepeatMasker conversion requires the installed ensembl-tools package"
        ) from error
    return get_repeat_type(repeat_class)


def default_repeatmasker_output_path(repeatmasker_path: str | Path) -> Path:
    """Return the default GTF path for a RepeatMasker output file."""

    path = Path(repeatmasker_path)
    name = path.name
    for suffix in (".out.gz", ".out"):
        if name.endswith(suffix):
            return path.with_name(f"{name[:-len(suffix)]}_repeatmasker.gtf")
    return path.with_name(f"{path.stem}_repeatmasker.gtf")


@contextmanager
def open_text_maybe_gzip(path: str | Path) -> Iterator[TextIO]:
    """Open plain-text or gzip-compressed files for text reading."""

    input_path = Path(path)
    if input_path.suffix == ".gz":
        with gzip.open(input_path, "rt", encoding="utf-8") as handle:
            yield handle
    else:
        with input_path.open("r", encoding="utf-8") as handle:
            yield handle


def load_refseq_name_map(assembly_report_path: str | Path) -> dict[str, str]:
    """Load RefSeq accession to Ensembl-style seq_region name mappings."""

    return load_assembly_report_name_maps(assembly_report_path)["RefSeq_genomic"]


def load_assembly_report_name_maps(
    assembly_report_path: str | Path,
) -> dict[str, dict[str, str]]:
    """Load RefSeq and GenBank accession to seq_region name mappings."""

    accession_maps: dict[str, dict[str, str]] = {
        "RefSeq_genomic": {},
        "INSDC": {},
    }
    report_path = Path(assembly_report_path)
    with report_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue

            columns = line.strip().split("\t")
            if len(columns) < 7:
                continue

            sequence_name = columns[0]
            assigned_molecule = columns[2]
            sequence_name_for_mapping = (
                assigned_molecule
                if columns[1] == "assembled-molecule" and assigned_molecule != "na"
                else sequence_name
            )

            for index, db_name in ((4, "INSDC"), (6, "RefSeq_genomic")):
                accession = columns[index]
                if accession and accession != "na":
                    accession_maps[db_name][accession] = sequence_name_for_mapping

    LOGGER.info(
        "Loaded %s RefSeq and %s GenBank sequence name mappings from %s",
        len(accession_maps["RefSeq_genomic"]),
        len(accession_maps["INSDC"]),
        report_path,
    )
    return accession_maps


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
        output_path.open("w", encoding="utf-8") as output_handle,
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
        output.open("w", encoding="utf-8") as output_handle,
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


def convert_repeatmasker_to_gtf(
    repeatmasker_path: str | Path,
    assembly_report_path: str | Path,
    output_path: str | Path | None = None,
    logger: logging.Logger | None = None,
) -> Path:
    """Convert NCBI RepeatMasker ``*_rm.out`` output to single-line GTF."""

    log = logger or LOGGER
    output = (
        Path(output_path)
        if output_path is not None
        else default_repeatmasker_output_path(repeatmasker_path)
    )
    refseq_to_name = load_refseq_name_map(assembly_report_path)
    repeat_ids: dict[str, int] = {}
    converted_records = 0
    unmapped_sequences = 0

    with (
        open_text_maybe_gzip(repeatmasker_path) as input_handle,
        output.open("w", encoding="utf-8") as output_handle,
    ):
        for line_number, raw_line in enumerate(input_handle, start=1):
            if not raw_line.strip() or not raw_line.lstrip()[:1].isdigit():
                continue
            record = _parse_repeatmasker_row(raw_line.split())
            if record is None:
                log.debug("Skipping short RepeatMasker row %s", line_number)
                continue

            source_seq_name = record.sequence_name
            seq_name = refseq_to_name.get(source_seq_name, source_seq_name)
            if seq_name == source_seq_name and source_seq_name not in refseq_to_name:
                unmapped_sequences += 1
            repeat_ids[seq_name] = repeat_ids.get(seq_name, 0) + 1
            repeat_type = _get_repeat_type(record.repeat_class)
            attributes = (
                f"repeat_id {repeat_ids[seq_name]}; "
                f'repeat_name "{record.repeat_name}"; '
                f'repeat_class "{record.repeat_class}"; '
                f'repeat_type "{repeat_type}"; '
                f'repeat_start "{record.consensus_start}"; '
                f'repeat_end "{record.consensus_end}"; '
                f'score "{record.score}";'
            )
            output_handle.write(
                f"{seq_name}\tRepeatMasker\trepeat\t{record.genomic_start}"
                f"\t{record.genomic_end}\t.\t{record.strand}\t.\t{attributes}\n"
            )
            converted_records += 1

    if unmapped_sequences:
        log.warning(
            "%s RepeatMasker rows used unmapped sequence names",
            unmapped_sequences,
        )
    log.info("Wrote RepeatMasker GTF to %s (%s records)", output, converted_records)
    return output

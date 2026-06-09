"""Shared constants for RefSeq-to-Ensembl loading utilities."""

from __future__ import annotations

from typing import Final


NCBI_GROUPS: Final[tuple[str, ...]] = (
    "vertebrate_mammalian",
    "vertebrate_other",
    "plant",
    "invertebrate",
    "fungi",
    "protozoa",
)

ASSEMBLY_SUMMARY_MIN_COLUMNS: Final[int] = 20

ANNOTATION_METADATA_FIELDS: Final[tuple[str, ...]] = (
    "group",
    "species_name",
    "taxon_id",
    "assembly_accession",
    "assembly_name",
    "ftp_path",
    "gff3_ftp",
    "fasta_ftp",
    "assembly_report_ftp",
    "gff3_local",
    "fasta_local",
    "assembly_report_local",
    "downloaded",
)

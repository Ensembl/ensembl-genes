"""Typed data structures for NCBI RefSeq assembly discovery and downloads."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AssemblyPaths:
    """FTP and local paths for one RefSeq assembly."""

    assembly_accession: str
    ftp_path: str
    ftp_base: str
    assembly_dir: Path
    gff_url: str
    fasta_url: str
    assembly_report_url: str
    gff_local: Path
    fasta_local: Path
    assembly_report_local: Path


@dataclass(frozen=True)
class AssemblySummaryRecord:
    """Parsed row from an NCBI RefSeq assembly summary file."""

    group: str
    assembly_accession: str
    assembly_name: str
    taxon_id: str
    organism_name: str
    species_name: str
    version_status: str
    ftp_path: str
    paths: AssemblyPaths

    @property
    def downloaded(self) -> bool:
        """Return True when the expected local assembly directory exists."""

        return self.paths.assembly_dir.exists()

    def to_metadata_row(self) -> dict[str, str]:
        """Return a TSV-ready metadata row compatible with the old script."""

        return {
            "group": self.group,
            "species_name": self.species_name,
            "taxon_id": self.taxon_id,
            "assembly_accession": self.assembly_accession,
            "assembly_name": self.assembly_name,
            "ftp_path": self.ftp_path,
            "gff3_ftp": self.paths.gff_url,
            "fasta_ftp": self.paths.fasta_url,
            "assembly_report_ftp": self.paths.assembly_report_url,
            "gff3_local": str(self.paths.gff_local),
            "fasta_local": str(self.paths.fasta_local),
            "assembly_report_local": str(self.paths.assembly_report_local),
            "downloaded": "Downloaded" if self.downloaded else "Not downloaded",
        }

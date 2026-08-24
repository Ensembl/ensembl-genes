"""Assembly-report metadata and shared text-file helpers."""

# Assembly reports contain several independent header fields with different
# fallback rules, so the header parser intentionally has multiple branches.
# pylint: disable=too-many-branches

from __future__ import annotations

import gzip
import logging
import re
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import TextIO

from .ncbi_annotation_report import get_annotation_report_url

LOGGER = logging.getLogger(__name__)
DEFAULT_ASSEMBLY_WEB_ACCESSION_SOURCE = "NCBI"
DEFAULT_ASSEMBLY_WEB_ACCESSION_TYPE = "INSDC Assembly ID"


def species_db_name_token(species_name: str) -> str:
    """Return the species token used in derived Ensembl core DB names."""

    species_parts = species_name.split()
    if len(species_parts) < 2:
        raise ValueError("species_name must contain at least genus and species")
    genus, species = species_parts[:2]
    return "_".join(
        re.sub(r"[^a-z0-9]+", "_", token.lower()).strip("_")
        for token in (genus, species)
    )


def species_url_token(species_name: str) -> str:
    """Return the species token used in the core ``species.url`` metadata."""

    species_parts = species_name.split()
    if len(species_parts) < 2:
        raise ValueError("species_name must contain at least genus and species")
    genus, species = species_parts[:2]
    return f"{genus.capitalize()}_{species.lower()}"


def refseq_accession_db_token(assembly_accession: str) -> str:
    """Return a compact DB-safe token for a RefSeq assembly accession."""

    token = assembly_accession.lower().replace("_", "")
    token = token.replace(".", "v", 1)
    return re.sub(r"[^a-z0-9]+", "_", token).strip("_")


def ncbi_assembly_url(assembly_accession: str) -> str:
    """Return the NCBI assembly page for a GenBank or RefSeq accession."""

    if assembly_accession.upper().startswith(("GCA_", "GCF_")):
        return f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_accession}/"
    return ""


@contextmanager
def open_text_maybe_gzip(path: str | Path) -> Iterator[TextIO]:
    """Open plain-text or gzip-compressed GFF3 files for text reading."""

    input_path = Path(path)
    if input_path.suffix == ".gz":
        with gzip.open(input_path, "rt", encoding="utf-8") as handle:
            yield handle
    else:
        with input_path.open("r", encoding="utf-8") as handle:
            yield handle


def normalize_assembly_report_date(raw_date: str) -> str:
    """Normalize a header date from the NCBI assembly report."""

    value = raw_date.strip()
    if not value:
        return ""
    value = value.replace("/", "-")
    if len(value) >= 10:
        value = value[:10]
    return value


def parse_assembly_report_metadata(
    assembly_report_path: str | Path,
    species_name: str,
    assembly_accession: str,
) -> dict[str, str]:
    """Read assembly-report header metadata needed for core DB bootstrapping."""

    metadata: dict[str, str] = {
        "species.scientific_name": species_name,
        "species.common_name": "",
        "species.display_name": species_name,
        "species.production_name": "",
        "species.url": "",
        "species.taxonomy_id": "",
        "assembly.accession": assembly_accession,
        "assembly.alt_accession": "",
        "assembly.name": assembly_accession,
        "assembly.date": "",
        "assembly.provider_url": ncbi_assembly_url(assembly_accession),
        "assembly.provider": "",
        "genebuild.provider_url": "",
        "annotation.provider_url": (
            "https://www.ncbi.nlm.nih.gov/refseq/"
            if assembly_accession.upper().startswith("GCF_")
            else ""
        ),
        "genebuild.start_date": "",
        "assembly.default": "1",
        "assembly.web_accession_source": DEFAULT_ASSEMBLY_WEB_ACCESSION_SOURCE,
        "assembly.web_accession_type": DEFAULT_ASSEMBLY_WEB_ACCESSION_TYPE,
    }

    metadata["species.production_name"] = (
        f"{species_db_name_token(species_name)}_"
        f"{refseq_accession_db_token(assembly_accession)}"
    )
    metadata["species.url"] = (
        f"{species_url_token(species_name)}_"
        f"{refseq_accession_db_token(assembly_accession)}"
    )

    if assembly_accession.upper().startswith("GCF_"):
        metadata["genebuild.provider_url"] = get_annotation_report_url(
            assembly_accession
        )

    if not str(assembly_report_path):
        return metadata
    report_path = Path(assembly_report_path)

    common_name_pattern = re.compile(r"^(?P<scientific>.+?)\s*\((?P<common>[^()]+)\)$")
    with open_text_maybe_gzip(report_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("#"):
                break
            if ":" not in line:
                continue
            header_key, header_value = line[1:].split(":", 1)
            header_key = header_key.strip().lower()
            header_value = header_value.strip()

            if header_key == "assembly name" and header_value:
                metadata["assembly.name"] = header_value
            elif (
                header_key == "genbank assembly accession"
                and header_value
                and assembly_accession.startswith("GCF")
            ):
                metadata["assembly.alt_accession"] = header_value
                metadata["annotation.provider_name"] = "NCBI RefSeq"
            elif header_key == "date" and header_value:
                normalized_date = normalize_assembly_report_date(header_value)
                metadata["assembly.date"] = normalized_date
                metadata["genebuild.start_date"] = normalized_date
            elif header_key == "taxid" and header_value:
                metadata["species.taxonomy_id"] = header_value
            elif header_key == "submitter" and header_value:
                metadata["assembly.provider_name"] = header_value
            elif header_key == "organism name" and header_value:
                match = common_name_pattern.match(header_value)
                if match:
                    metadata["species.scientific_name"] = match.group(
                        "scientific"
                    ).strip()
                    metadata["species.common_name"] = match.group("common").strip()
                else:
                    metadata["species.scientific_name"] = header_value
            elif header_key == "common name" and header_value:
                metadata["species.common_name"] = header_value

    if metadata["species.common_name"]:
        metadata["species.display_name"] = metadata["species.common_name"]
    else:
        metadata["species.display_name"] = metadata["species.scientific_name"]
    # add species.production_name
    metadata["species.production_name"] = (
        f"{species_db_name_token(metadata['species.scientific_name'])}_"
        f"{refseq_accession_db_token(metadata['assembly.accession'])}"
    )
    metadata["species.url"] = (
        f"{species_url_token(metadata['species.scientific_name'])}_"
        f"{refseq_accession_db_token(metadata['assembly.accession'])}"
    )

    if not metadata["genebuild.start_date"]:
        metadata["genebuild.start_date"] = datetime.now().strftime("%Y-%m-%d")

    return metadata

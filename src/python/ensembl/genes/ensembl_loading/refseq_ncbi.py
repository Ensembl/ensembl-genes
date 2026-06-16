"""NCBI RefSeq discovery and download helpers.

This module contains the reusable parts of listing RefSeq annotation metadata
and downloading the GFF3, genomic FASTA, and assembly report files. It avoids
printing directly and reports progress through the standard logging module.
"""

from __future__ import annotations

import csv
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

import requests  # type: ignore[import]

try:  # Support both package imports and direct same-directory imports.
    from . import refseq_constants as _refseq_constants
    from . import refseq_models as _refseq_models
except ImportError:  # pragma: no cover - used when run beside this file.
    import refseq_constants as _refseq_constants  # type: ignore[import,no-redef]
    import refseq_models as _refseq_models  # type: ignore[import,no-redef]


LOGGER = logging.getLogger(__name__)


def accession_subdir(assembly_accession: str) -> str:
    """Return the NCBI-style three-level accession directory."""

    try:
        numeric_accession = assembly_accession.split("_", 1)[1].split(".", 1)[0]
    except IndexError as exc:
        raise ValueError(
            f"Assembly accession must include an underscore: {assembly_accession}"
        ) from exc
    padded = numeric_accession.zfill(9)
    return "/".join(padded[index : index + 3] for index in range(0, 9, 3))


def build_assembly_paths(
    base_dir: str | Path, assembly_accession: str, ftp_path: str
) -> _refseq_models.AssemblyPaths:
    """Build all local and remote paths for one RefSeq assembly."""

    accession_prefix = assembly_accession.split("_", 1)[0]
    assembly_dir = (
        Path(base_dir)
        / accession_prefix
        / accession_subdir(assembly_accession)
        / assembly_accession
    )
    ftp_base = ftp_path.rstrip("/").split("/")[-1]
    gff_url = f"{ftp_path}/{ftp_base}_genomic.gff.gz"
    fasta_url = f"{ftp_path}/{ftp_base}_genomic.fna.gz"
    assembly_report_url = f"{ftp_path}/{ftp_base}_assembly_report.txt"
    return _refseq_models.AssemblyPaths(
        assembly_accession=assembly_accession,
        ftp_path=ftp_path,
        ftp_base=ftp_base,
        assembly_dir=assembly_dir,
        gff_url=gff_url,
        fasta_url=fasta_url,
        assembly_report_url=assembly_report_url,
        gff_local=assembly_dir / Path(gff_url).name,
        fasta_local=assembly_dir / Path(fasta_url).name,
        assembly_report_local=assembly_dir / Path(assembly_report_url).name,
    )


def parse_assembly_summary(
    summary_text: str,
    group: str,
    base_dir: str | Path = "refseq_data",
    logger: logging.Logger | None = None,
) -> list[_refseq_models.AssemblySummaryRecord]:
    """Parse one NCBI assembly_summary.txt file."""

    log = logger or LOGGER
    records: list[_refseq_models.AssemblySummaryRecord] = []
    for line_number, line in enumerate(summary_text.splitlines(), start=1):
        if line.startswith("#") or not line.strip():
            continue

        columns = line.split("\t")
        if len(columns) < _refseq_constants.ASSEMBLY_SUMMARY_MIN_COLUMNS:
            log.debug(
                "Skipping short assembly summary row %s in %s", line_number, group
            )
            continue

        ftp_path = columns[19]
        if ftp_path == "na":
            log.debug("Skipping assembly %s because ftp_path is 'na'", columns[0])
            continue

        organism_name = columns[7]
        species_name = " ".join(organism_name.split()[:2])
        paths = build_assembly_paths(base_dir, columns[0], ftp_path)
        records.append(
            _refseq_models.AssemblySummaryRecord(
                group=group,
                assembly_accession=columns[0],
                assembly_name=columns[15],
                taxon_id=columns[5],
                organism_name=organism_name,
                species_name=species_name,
                version_status=columns[10],
                ftp_path=ftp_path,
                paths=paths,
            )
        )
    return records


def fetch_assembly_summary(
    group: str,
    base_dir: str | Path = "refseq_data",
    session: requests.Session | None = None,
    timeout: int = 20,
    logger: logging.Logger | None = None,
) -> list[_refseq_models.AssemblySummaryRecord]:
    """Fetch and parse the RefSeq assembly summary for one clade group."""

    log = logger or LOGGER
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{group}/assembly_summary.txt"
    client = session or requests
    log.info("Fetching RefSeq assembly summary for %s", group)
    response = client.get(url, timeout=timeout)
    response.raise_for_status()
    records = parse_assembly_summary(response.text, group, base_dir, logger=log)
    log.info("Parsed %s assemblies for %s", len(records), group)
    return records


def fetch_groups(
    base_dir: str | Path = "refseq_data",
    group_filter: str | None = None,
    session: requests.Session | None = None,
    timeout: int = 20,
    logger: logging.Logger | None = None,
) -> list[_refseq_models.AssemblySummaryRecord]:
    """Fetch assembly summaries for all configured groups or one group."""

    log = logger or LOGGER
    groups = (group_filter,) if group_filter else _refseq_constants.NCBI_GROUPS
    records: list[_refseq_models.AssemblySummaryRecord] = []
    for group in groups:
        try:
            records.extend(
                fetch_assembly_summary(
                    group,
                    base_dir=base_dir,
                    session=session,
                    timeout=timeout,
                    logger=log,
                )
            )
        except requests.RequestException as exc:
            log.error("Failed to fetch RefSeq summary for %s: %s", group, exc)
    return records


def write_annotation_metadata_tsv(
    records: Iterable[_refseq_models.AssemblySummaryRecord], output_tsv: str | Path
) -> Path:
    """Write assembly metadata rows to a tab-separated file."""

    output_path = Path(output_tsv)
    sorted_records = sorted(records, key=lambda row: (row.group, row.species_name))
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=_refseq_constants.ANNOTATION_METADATA_FIELDS,
        )
        writer.writeheader()
        for record in sorted_records:
            writer.writerow(record.to_metadata_row())
    LOGGER.info("Wrote RefSeq annotation metadata to %s", output_path)
    return output_path


def summarize_annotations(
    records: Iterable[_refseq_models.AssemblySummaryRecord],
) -> dict[str, list[str]]:
    """Return clade to species summaries for listing and reporting."""

    summaries: dict[str, list[str]] = {}
    for record in records:
        status = "Downloaded" if record.downloaded else "Not downloaded"
        summaries.setdefault(record.group, []).append(
            f"{record.species_name} [{status}]"
        )

    return {
        group: sorted(species_summaries)
        for group, species_summaries in sorted(summaries.items())
    }


def list_available_annotations(
    base_dir: str | Path = "refseq_data",
    group_filter: str | None = None,
    max_print: int | None = 50,
    return_dict: bool = False,
    output_tsv: str | Path | None = "refseq_annotation_metadata.tsv",
    session: requests.Session | None = None,
    logger: logging.Logger | None = None,
) -> dict[str, list[str]] | None:
    """List RefSeq genome annotations grouped by clade.

    Metadata can optionally be exported to TSV. Summary output is sent through
    logging so callers can decide whether it is shown, captured, or suppressed.
    """

    log = logger or LOGGER
    records = fetch_groups(
        base_dir=base_dir,
        group_filter=group_filter,
        session=session,
        logger=log,
    )
    if output_tsv:
        write_annotation_metadata_tsv(records, output_tsv)

    summaries = summarize_annotations(records)
    for group, species_list in summaries.items():
        log.info("%s: %s species", group, len(species_list))
        shown_species = species_list if max_print is None else species_list[:max_print]
        for species_summary in shown_species:
            log.info(" - %s", species_summary)
        if max_print is not None and len(species_list) > max_print:
            log.info(" ... and %s more species", len(species_list) - max_print)

    return summaries if return_dict else None


def find_annotation_targets(
    base_dir: str | Path = "refseq_data",
    species_name: str | None = None,
    assembly_acc: str | None = None,
    group: str | None = None,
    session: requests.Session | None = None,
    logger: logging.Logger | None = None,
) -> list[_refseq_models.AssemblySummaryRecord]:
    """Locate assembly records by accession, species, or latest group batch."""

    log = logger or LOGGER
    groups = (group,) if group else _refseq_constants.NCBI_GROUPS

    if assembly_acc:
        for group_name in groups:
            try:
                records = fetch_assembly_summary(
                    group_name, base_dir=base_dir, session=session, logger=log
                )
            except requests.RequestException as exc:
                log.error("Failed to fetch RefSeq summary for %s: %s", group_name, exc)
                continue

            for record in records:
                if record.assembly_accession == assembly_acc:
                    return [record]
        raise ValueError(f"Assembly {assembly_acc} not found in RefSeq summaries.")

    if species_name:
        requested = species_name.lower()
        for group_name in groups:
            try:
                records = fetch_assembly_summary(
                    group_name, base_dir=base_dir, session=session, logger=log
                )
            except requests.RequestException as exc:
                log.error("Failed to fetch RefSeq summary for %s: %s", group_name, exc)
                continue

            for record in records:
                if record.organism_name.lower().startswith(requested):
                    return [record]
        raise ValueError(f"Species {species_name} not found in RefSeq summaries.")

    if group:
        return [
            record
            for record in fetch_assembly_summary(
                group, base_dir=base_dir, session=session, logger=log
            )
            if record.version_status == "latest"
        ]

    raise ValueError(
        "Must provide at least one of: assembly_acc, species_name, or group"
    )


def download_file(
    url: str,
    destination: str | Path,
    timeout: int = 30,
    chunk_size: int = 8192,
    logger: logging.Logger | None = None,
) -> Path:
    """Download one URL to destination unless the file already exists."""

    log = logger or LOGGER
    destination_path = Path(destination)
    if destination_path.exists():
        log.info("Already exists: %s", destination_path)
        return destination_path

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    log.info("Downloading %s to %s", url, destination_path)
    response = requests.get(url, stream=True, timeout=timeout)
    response.raise_for_status()
    with destination_path.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk:
                handle.write(chunk)

    if destination_path.stat().st_size == 0:
        raise IOError(f"Downloaded file is empty: {destination_path}")

    return destination_path


def download_assembly(
    record: _refseq_models.AssemblySummaryRecord,
    logger: logging.Logger | None = None,
) -> _refseq_models.AssemblyPaths:
    """Download GFF3, FASTA, and assembly report files for one assembly."""

    log = logger or LOGGER
    paths = record.paths
    paths.assembly_dir.mkdir(parents=True, exist_ok=True)
    download_file(paths.gff_url, paths.gff_local, logger=log)
    download_file(paths.fasta_url, paths.fasta_local, logger=log)
    download_file(paths.assembly_report_url, paths.assembly_report_local, logger=log)
    return paths


def download_annotations(
    base_dir: str | Path = "refseq_data",
    species_name: str | None = None,
    assembly_acc: str | None = None,
    group: str | None = None,
    max_workers: int = 2,
    session: requests.Session | None = None,
    logger: logging.Logger | None = None,
) -> list[_refseq_models.AssemblyPaths]:
    """Download RefSeq annotation files selected by accession, species, or group."""

    log = logger or LOGGER
    targets = find_annotation_targets(
        base_dir=base_dir,
        species_name=species_name,
        assembly_acc=assembly_acc,
        group=group,
        session=session,
        logger=log,
    )
    if not targets:
        log.warning("No RefSeq annotation targets found")
        return []

    worker_count = max(1, min(max_workers, len(targets)))
    downloaded_paths: list[_refseq_models.AssemblyPaths] = []
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        future_to_record = {
            executor.submit(download_assembly, target, log): target
            for target in targets
        }
        for future in as_completed(future_to_record):
            record = future_to_record[future]
            try:
                downloaded_paths.append(future.result())
                log.info(
                    "Downloaded annotation files for %s", record.assembly_accession
                )
            except Exception:
                log.exception(
                    "Failed to download annotation files for %s",
                    record.assembly_accession,
                )
                raise

    return downloaded_paths

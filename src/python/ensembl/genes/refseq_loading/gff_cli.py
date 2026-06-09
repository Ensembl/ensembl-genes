"""Unified command-line wrapper for generic GFF3 and RefSeq loading."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

try:  # Support both package imports and direct same-directory imports.
    from .gff_core_loader import load_gff_features_to_core, load_to_ensembl_core
    from .gff_source_config import available_source_configs, get_source_config
    from .refseq_conversion import convert_fna_headers, convert_gff_to_ensembl
    from .refseq_ncbi import download_annotations, list_available_annotations
except ImportError:  # pragma: no cover - used when run beside this file.
    from gff_core_loader import (  # type: ignore
        load_gff_features_to_core,
        load_to_ensembl_core,
    )
    from gff_source_config import (  # type: ignore
        available_source_configs,
        get_source_config,
    )
    from refseq_conversion import (  # type: ignore
        convert_fna_headers,
        convert_gff_to_ensembl,
    )
    from refseq_ncbi import (  # type: ignore
        download_annotations,
        list_available_annotations,
    )


LOGGER = logging.getLogger(__name__)


def configure_logging(log_level: str, log_file: str | None = None) -> None:
    """Configure standard logging for CLI execution."""

    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        log_path = Path(log_file)
        if log_path.parent != Path("."):
            log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path))

    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=handlers,
    )


def positive_int(raw_value: str) -> int:
    """Parse a positive integer argument."""

    value = int(raw_value)
    if value < 1:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return value


def tcp_port(raw_value: str) -> int:
    """Parse a valid TCP port number."""

    value = positive_int(raw_value)
    if value > 65535:
        raise argparse.ArgumentTypeError("value must be a valid TCP port")
    return value


def default_converted_fna_path(assembly_dir: Path, ftp_base: str) -> Path:
    """Return the default converted FASTA path for one RefSeq assembly."""

    return assembly_dir / f"{ftp_base}_genomic_ensembl.fna"


def default_converted_gff_path(assembly_dir: Path, ftp_base: str) -> Path:
    """Return the default converted GFF3 path for one RefSeq assembly."""

    return assembly_dir / f"{ftp_base}_ensembl.gff3"


def add_source_option(
    parser: argparse.ArgumentParser,
    default: str = "generic",
) -> None:
    """Add the source config selector used by load-capable commands."""

    parser.add_argument(
        "--source",
        default=default,
        choices=available_source_configs(),
        help="Registered GFF source configuration",
    )


def add_existing_core_db_options(parser: argparse.ArgumentParser) -> None:
    """Add DB options for loading features into an existing core DB."""

    parser.add_argument(
        "--db-name",
        "--db_name",
        dest="db_name",
        required=True,
        help="Existing core DB name",
    )
    parser.add_argument(
        "--db-host",
        "--db_host",
        dest="db_host",
        required=True,
        help="MySQL host",
    )
    parser.add_argument(
        "--db-user",
        "--db_user",
        dest="db_user",
        required=True,
        help="MySQL user",
    )
    parser.add_argument(
        "--db-password",
        "--db_password",
        dest="db_password",
        required=True,
        help="MySQL password",
    )
    parser.add_argument(
        "--db-port",
        "--db_port",
        dest="db_port",
        type=tcp_port,
        required=True,
        help="MySQL port",
    )
    parser.add_argument(
        "--coord-system-id",
        type=int,
        help="Explicit coord_system_id; overrides name/version lookup",
    )
    parser.add_argument(
        "--coord-system-name",
        default="primary_assembly",
        help="coord_system.name used to resolve seq_region IDs",
    )
    parser.add_argument(
        "--coord-system-version",
        help="Optional coord_system.version used with --coord-system-name",
    )


def add_create_core_db_options(
    parser: argparse.ArgumentParser,
    required: bool = True,
) -> None:
    """Add DB options for creating/populating a core DB from FASTA and GFF3."""

    required_with_load_core = "" if required else " (required with --load-core)"
    parser.add_argument(
        "--species-name",
        required=required,
        help="Scientific species name used for core DB metadata",
    )
    parser.add_argument(
        "--assembly-accession",
        required=required,
        help="Assembly accession or source assembly identifier",
    )
    parser.add_argument(
        "--db-host",
        "--db_host",
        dest="db_host",
        required=required,
        help=f"MySQL host{required_with_load_core}",
    )
    parser.add_argument(
        "--db-user",
        "--db_user",
        dest="db_user",
        required=required,
        help=f"MySQL user{required_with_load_core}",
    )
    parser.add_argument(
        "--db-password",
        "--db_password",
        dest="db_password",
        required=required,
        help=f"MySQL password{required_with_load_core}",
    )
    parser.add_argument(
        "--db-port",
        "--db_port",
        dest="db_port",
        type=tcp_port,
        required=required,
        help=f"MySQL port{required_with_load_core}",
    )
    parser.add_argument(
        "--schema-sql-path",
        help=(
            "Ensembl core schema SQL path; defaults to bundled "
            "config/core_schema.sql, use an empty string to skip"
        ),
    )


def add_refseq_target_options(parser: argparse.ArgumentParser) -> None:
    """Add RefSeq target selection options."""

    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument(
        "--assembly-acc",
        help="Download one RefSeq assembly accession, for example GCF_000001635.27",
    )
    target_group.add_argument(
        "--species-name",
        help="Download the first matching scientific species name",
    )
    target_group.add_argument(
        "--group",
        help="Download latest assemblies for a RefSeq clade group",
    )


def run_load_features(args: argparse.Namespace) -> int:
    """Load a user-supplied GFF3 into an existing core DB."""

    source_config = get_source_config(args.source)
    summary = load_gff_features_to_core(
        gff_path=args.gff,
        db_name=args.db_name,
        db_host=args.db_host,
        db_user=args.db_user,
        db_password=args.db_password,
        db_port=args.db_port,
        coord_system_id=args.coord_system_id,
        coord_system_name=args.coord_system_name,
        coord_system_version=args.coord_system_version,
        source_config=source_config,
    )
    LOGGER.info(
        "Loaded GFF feature summary: %s genes, %s transcripts, %s CDS transcript groups",
        summary["genes"],
        summary["transcripts"],
        summary["cds_transcript_groups"],
    )
    return 0


def run_create_core(args: argparse.Namespace) -> int:
    """Create/populate a core DB from converted FASTA and converted GFF3."""

    source_config = get_source_config(args.source)
    db_name = load_to_ensembl_core(
        converted_gff_path=args.gff,
        converted_fna_path=args.fna,
        assembly_report_path=args.assembly_report or "",
        species_name=args.species_name,
        assembly_accession=args.assembly_accession,
        db_host=args.db_host,
        db_user=args.db_user,
        db_password=args.db_password,
        db_port=args.db_port,
        schema_sql_path=args.schema_sql_path,
        source_config=source_config,
    )
    LOGGER.info("Loaded core database: %s", db_name)
    return 0


def run_refseq_list(args: argparse.Namespace) -> int:
    """List RefSeq assemblies."""

    list_available_annotations(
        base_dir=args.base_dir,
        group_filter=args.group,
        max_print=args.max_print,
        return_dict=False,
        output_tsv=args.output_tsv,
    )
    return 0


def run_refseq_download(args: argparse.Namespace) -> int:
    """Download RefSeq annotation files."""

    paths = download_annotations(
        base_dir=args.base_dir,
        species_name=args.species_name,
        assembly_acc=args.assembly_acc,
        group=args.group,
        max_workers=args.max_workers,
    )
    for assembly_paths in paths:
        LOGGER.info("Assembly directory: %s", assembly_paths.assembly_dir)
    return 0


def run_refseq_convert_fna(args: argparse.Namespace) -> int:
    """Convert RefSeq FASTA headers."""

    output_path = convert_fna_headers(
        args.fna,
        args.assembly_report,
        args.output,
    )
    LOGGER.info("Converted FASTA: %s", output_path)
    return 0


def run_refseq_convert_gff(args: argparse.Namespace) -> int:
    """Convert RefSeq GFF3 seqids."""

    output_path = convert_gff_to_ensembl(
        args.gff,
        args.assembly_report,
        args.output,
        chrom_filter=args.chrom_filter,
    )
    LOGGER.info("Converted GFF3: %s", output_path)
    return 0


def run_refseq_pipeline(args: argparse.Namespace) -> int:
    """Download plus RefSeq FASTA/GFF3 conversion, optionally followed by DB load."""

    if args.load_core and args.group:
        LOGGER.error(
            "--load-core with --group is not supported; load one assembly at a time"
        )
        return 2

    paths_list = download_annotations(
        base_dir=args.base_dir,
        species_name=args.download_species_name,
        assembly_acc=args.assembly_acc,
        group=args.group,
        max_workers=args.max_workers,
    )
    if not paths_list:
        LOGGER.warning("No assemblies downloaded or found")
        return 1

    if len(paths_list) > 1 and (args.converted_fna or args.converted_gff):
        LOGGER.error("Custom converted output paths can only be used for one assembly")
        return 2

    for assembly_paths in paths_list:
        converted_fna = args.converted_fna or default_converted_fna_path(
            assembly_paths.assembly_dir,
            assembly_paths.ftp_base,
        )
        converted_gff = args.converted_gff or default_converted_gff_path(
            assembly_paths.assembly_dir,
            assembly_paths.ftp_base,
        )

        LOGGER.info("Converting FASTA for %s", assembly_paths.assembly_accession)
        convert_fna_headers(
            assembly_paths.fasta_local,
            assembly_paths.assembly_report_local,
            converted_fna,
        )
        LOGGER.info("Converting GFF3 for %s", assembly_paths.assembly_accession)
        convert_gff_to_ensembl(
            assembly_paths.gff_local,
            assembly_paths.assembly_report_local,
            converted_gff,
            chrom_filter=args.chrom_filter,
        )

        if args.load_core:
            source_config = get_source_config(args.source)
            db_name = load_to_ensembl_core(
                converted_gff_path=converted_gff,
                converted_fna_path=converted_fna,
                assembly_report_path=assembly_paths.assembly_report_local,
                species_name=args.species_name,
                assembly_accession=args.assembly_accession
                or assembly_paths.assembly_accession,
                db_host=args.db_host,
                db_user=args.db_user,
                db_password=args.db_password,
                db_port=args.db_port,
                schema_sql_path=args.schema_sql_path,
                source_config=source_config,
            )
            LOGGER.info("Loaded core database: %s", db_name)

    return 0


def add_load_features_parser(subparsers: argparse._SubParsersAction) -> None:
    """Add the generic existing-core GFF feature loader subcommand."""

    parser = subparsers.add_parser(
        "load-features",
        help="Load a GFF3 into an existing Ensembl core DB",
    )
    parser.add_argument("gff", help="Input GFF3 path, plain text or .gz")
    add_existing_core_db_options(parser)
    add_source_option(parser, default="generic")
    parser.set_defaults(func=run_load_features)


def add_create_core_parser(subparsers: argparse._SubParsersAction) -> None:
    """Add the FASTA plus GFF3 core creation subcommand."""

    parser = subparsers.add_parser(
        "create-core",
        help="Create/populate a core DB from converted FASTA and GFF3",
    )
    parser.add_argument("gff", help="Converted GFF3")
    parser.add_argument("fna", help="Converted FASTA")
    parser.add_argument(
        "assembly_report",
        nargs="?",
        help="Optional source assembly report; retained for RefSeq compatibility",
    )
    add_create_core_db_options(parser)
    add_source_option(parser, default="generic")
    parser.set_defaults(func=run_create_core)


def add_refseq_parser(subparsers: argparse._SubParsersAction) -> None:
    """Add RefSeq-specific subcommands under the unified CLI."""

    parser = subparsers.add_parser(
        "refseq",
        help="RefSeq/NCBI listing, download, conversion, and pipeline commands",
    )
    refseq_subparsers = parser.add_subparsers(dest="refseq_command", required=True)

    list_parser = refseq_subparsers.add_parser(
        "list",
        help="List available RefSeq annotations and optionally write metadata TSV",
    )
    list_parser.add_argument("--base-dir", default="refseq_data")
    list_parser.add_argument("--group", help="Restrict to one RefSeq clade group")
    list_parser.add_argument("--max-print", type=int, default=50)
    list_parser.add_argument(
        "--output-tsv",
        default="refseq_annotation_metadata.tsv",
        help="Metadata TSV output path; use an empty string to disable",
    )
    list_parser.set_defaults(func=run_refseq_list)

    download_parser = refseq_subparsers.add_parser(
        "download",
        help="Download RefSeq GFF3, FASTA, and assembly report files",
    )
    download_parser.add_argument("--base-dir", default="refseq_data")
    add_refseq_target_options(download_parser)
    download_parser.add_argument("--max-workers", type=positive_int, default=2)
    download_parser.set_defaults(func=run_refseq_download)

    fna_parser = refseq_subparsers.add_parser(
        "convert-fna",
        help="Convert RefSeq FASTA headers to Ensembl-style seq_region names",
    )
    fna_parser.add_argument("fna", help="Input RefSeq genomic FASTA")
    fna_parser.add_argument("assembly_report", help="NCBI assembly report")
    fna_parser.add_argument("output", help="Output converted FASTA")
    fna_parser.set_defaults(func=run_refseq_convert_fna)

    gff_parser = refseq_subparsers.add_parser(
        "convert-gff",
        help="Convert RefSeq GFF3 seq IDs to Ensembl-style seq_region names",
    )
    gff_parser.add_argument("gff", help="Input RefSeq GFF3, plain text or .gz")
    gff_parser.add_argument("assembly_report", help="NCBI assembly report")
    gff_parser.add_argument("--output", help="Output converted GFF3")
    gff_parser.add_argument(
        "--chrom-filter",
        action="append",
        help="Keep one sequence by RefSeq accession or converted name; repeatable",
    )
    gff_parser.set_defaults(func=run_refseq_convert_gff)

    pipeline_parser = refseq_subparsers.add_parser(
        "run",
        help="Download and convert RefSeq annotations, with optional core DB load",
    )
    pipeline_parser.add_argument("--base-dir", default="refseq_data")
    pipeline_targets = pipeline_parser.add_mutually_exclusive_group(required=True)
    pipeline_targets.add_argument("--assembly-acc")
    pipeline_targets.add_argument(
        "--download-species-name",
        help="Species name used only for selecting a RefSeq download target",
    )
    pipeline_targets.add_argument("--group")
    pipeline_parser.add_argument("--max-workers", type=positive_int, default=2)
    pipeline_parser.add_argument("--converted-fna", type=Path)
    pipeline_parser.add_argument("--converted-gff", type=Path)
    pipeline_parser.add_argument(
        "--chrom-filter",
        action="append",
        help="Keep one sequence by RefSeq accession or converted name; repeatable",
    )
    pipeline_parser.add_argument(
        "--load-core",
        action="store_true",
        help="Load converted files into a core DB after conversion",
    )
    add_create_core_db_options(pipeline_parser, required=False)
    add_source_option(pipeline_parser, default="refseq")
    pipeline_parser.set_defaults(func=run_refseq_pipeline)


def build_parser() -> argparse.ArgumentParser:
    """Build the unified GFF loader argument parser."""

    parser = argparse.ArgumentParser(
        prog="gff-loader",
        description=(
            "Load GFF3 features into Ensembl core databases, with optional "
            "RefSeq/NCBI download and conversion workflows."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Logging verbosity",
    )
    parser.add_argument(
        "--log-file",
        help="Optional file path for writing the same log output shown on screen",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_load_features_parser(subparsers)
    add_create_core_parser(subparsers)
    add_refseq_parser(subparsers)
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate cross-argument requirements not expressible in argparse."""

    if args.command == "refseq" and args.refseq_command == "run" and args.load_core:
        missing = [
            name
            for name in (
                "species_name",
                "db_host",
                "db_port",
                "db_user",
                "db_password",
            )
            if getattr(args, name) is None
        ]
        if missing:
            parser.error(
                "refseq run --load-core requires: "
                + ", ".join("--" + name.replace("_", "-") for name in missing)
            )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the unified GFF loading CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "output_tsv", None) == "":
        args.output_tsv = None
    validate_args(parser, args)
    configure_logging(args.log_level, args.log_file)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

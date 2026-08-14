"""Load converted GFF3 annotations into an Ensembl-style core database."""

from __future__ import annotations

import logging
from pathlib import Path

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=too-many-branches
# pylint: disable=too-many-return-statements
try:  # Support both package imports and direct same-directory imports.
    from .gff_annotation import (
        prepare_annotation_for_load,
    )
    from .gff_core_database import (
        allocate_numeric_ids,
        allocate_numeric_ids_from_core,
        connect_mysql,
        derive_core_db_name,
        get_coord_system_id,
        get_or_create_analysis,
        initialise_core_tables,
        insert_genes,
        insert_transcripts_and_exons,
        insert_translations,
        load_existing_seq_region_ids,
        load_schema_sql,
        load_seq_regions_from_fna,
        resolve_schema_sql_path,
    )
    from .gff_quality_check import emit_quality_report, run_core_load_quality_check
    from .gff_source_config import GENERIC_GFF_CONFIG, REFSEQ_CONFIG, GffSourceConfig
except ImportError:  # pragma: no cover - used when run beside this file.
    from gff_annotation import (  # type: ignore
        prepare_annotation_for_load,
    )
    from gff_core_database import (  # type: ignore
        allocate_numeric_ids,
        allocate_numeric_ids_from_core,
        connect_mysql,
        derive_core_db_name,
        get_coord_system_id,
        get_or_create_analysis,
        initialise_core_tables,
        insert_genes,
        insert_transcripts_and_exons,
        insert_translations,
        load_existing_seq_region_ids,
        load_schema_sql,
        load_seq_regions_from_fna,
        resolve_schema_sql_path,
    )
    from gff_quality_check import (  # type: ignore
        emit_quality_report,
        run_core_load_quality_check,
    )
    from gff_source_config import (  # type: ignore
        GENERIC_GFF_CONFIG,
        REFSEQ_CONFIG,
        GffSourceConfig,
    )


LOGGER = logging.getLogger(__name__)


def load_to_ensembl_core(
    converted_gff_path: str | Path,
    converted_fna_path: str | Path,
    assembly_report_path: str | Path,
    species_name: str,
    assembly_accession: str,
    db_host: str,
    db_user: str,
    db_password: str,
    db_port: int,
    schema_sql_path: str | Path | None = None,
    logger: logging.Logger | None = None,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
    deduplicate_exons: bool = True,
) -> str:
    """Load converted GFF3/GTF and FASTA files into an Ensembl-style core database.

    The ``assembly_report_path`` argument is retained for API symmetry with the
    original script and for callers that pass the three canonical source files
    together. Sequence names must already be converted in the GFF3 and FASTA.

    Returns
    -------
    str
        The database name that was created or reused.
    """
    log = logger or LOGGER
    db_name = derive_core_db_name(
        species_name,
        assembly_accession,
        source_config=source_config,
    )

    log.info("Connecting to MySQL server %s:%s", db_host, db_port)
    connection = connect_mysql(
        db_host=db_host,
        db_user=db_user,
        db_password=db_password,
        db_port=db_port,
    )
    cursor = connection.cursor()
    try:
        cursor.execute(f"CREATE DATABASE IF NOT EXISTS {db_name}")
        cursor.execute(f"USE {db_name}")
        resolved_schema_sql_path = resolve_schema_sql_path(schema_sql_path)
        if resolved_schema_sql_path:
            log.info("Loading schema SQL from %s", resolved_schema_sql_path)
            load_schema_sql(cursor, resolved_schema_sql_path)
        else:
            log.info("Skipping schema SQL loading")

        coord_system_id, analysis_id = initialise_core_tables(
            cursor,
            species_name,
            assembly_accession,
            assembly_report_path=assembly_report_path,
            source_config=source_config,
        )
        seq_region_ids = load_seq_regions_from_fna(
            converted_fna_path,
            cursor,
            coord_system_id,
            logger=log,
            source_config=source_config,
        )
        annotation = prepare_annotation_for_load(
            converted_gff_path,
            logger=log,
            source_config=source_config,
        )
        gene_id_map, transcript_id_map = allocate_numeric_ids(annotation)
        insert_genes(
            cursor,
            annotation,
            seq_region_ids,
            gene_id_map,
            transcript_id_map,
            analysis_id,
            source_config=source_config,
        )
        exon_id_map, _per_transcript_coord_to_exon_id = insert_transcripts_and_exons(
            cursor,
            annotation,
            seq_region_ids,
            gene_id_map,
            transcript_id_map,
            analysis_id,
            source_config=source_config,
            deduplicate_exons=deduplicate_exons,
        )
        insert_translations(cursor, annotation, transcript_id_map, exon_id_map)
        quality_report = run_core_load_quality_check(
            cursor,
            annotation,
            gene_id_map,
            transcript_id_map,
            exon_id_map,
            source_config,
        )
        if not quality_report.passed:
            emit_quality_report(quality_report, log)
            raise ValueError(
                "Post-load GFF core quality check failed: "
                f"{quality_report.failure_summary()}"
            )
        connection.commit()
        emit_quality_report(quality_report, log)
        log.info("Loaded %s annotation into %s", source_config.name, db_name)
        return db_name
    except Exception:
        connection.rollback()
        log.exception(
            "Failed loading %s annotation into %s; transaction rolled back",
            source_config.name,
            db_name,
        )
        raise
    finally:
        connection.close()


def load_gff_features_to_core(
    gff_path: str | Path,
    db_name: str,
    db_host: str,
    db_user: str,
    db_password: str,
    db_port: int,
    coord_system_id: int | None = None,
    coord_system_name: str = "primary_assembly",
    coord_system_version: str | None = None,
    source_config: GffSourceConfig = GENERIC_GFF_CONFIG,
    logger: logging.Logger | None = None,
    deduplicate_exons: bool = True,
) -> dict[str, int]:
    """Load GFF3/GTF features into an existing Ensembl core database.

    This path is for generic GFF loading. It does not create the database, load
    schema SQL, insert assembly metadata, or load DNA. The target core database
    must already contain matching ``coord_system`` and ``seq_region`` rows for
    the GFF seqids.
    """

    log = logger or LOGGER
    annotation = prepare_annotation_for_load(
        gff_path,
        logger=log,
        source_config=source_config,
    )

    log.info("Connecting to MySQL core database %s on %s:%s", db_name, db_host, db_port)
    connection = connect_mysql(
        db_host=db_host,
        db_user=db_user,
        db_password=db_password,
        db_port=db_port,
        db_name=db_name,
    )
    cursor = connection.cursor()
    try:
        resolved_coord_system_id = get_coord_system_id(
            cursor,
            coord_system_id=coord_system_id,
            coord_system_name=coord_system_name,
            coord_system_version=coord_system_version,
        )
        analysis_id = get_or_create_analysis(cursor, source_config=source_config)
        seq_region_ids = load_existing_seq_region_ids(
            cursor,
            annotation,
            resolved_coord_system_id,
        )
        gene_id_map, transcript_id_map, first_exon_id = allocate_numeric_ids_from_core(
            cursor,
            annotation,
        )
        insert_genes(
            cursor,
            annotation,
            seq_region_ids,
            gene_id_map,
            transcript_id_map,
            analysis_id,
            source_config=source_config,
        )
        exon_id_map, _per_transcript_coord_to_exon_id = insert_transcripts_and_exons(
            cursor,
            annotation,
            seq_region_ids,
            gene_id_map,
            transcript_id_map,
            analysis_id,
            source_config=source_config,
            first_exon_id=first_exon_id,
            deduplicate_exons=deduplicate_exons,
        )
        insert_translations(cursor, annotation, transcript_id_map, exon_id_map)
        quality_report = run_core_load_quality_check(
            cursor,
            annotation,
            gene_id_map,
            transcript_id_map,
            exon_id_map,
            source_config,
        )
        if not quality_report.passed:
            emit_quality_report(quality_report, log)
            raise ValueError(
                "Post-load GFF core quality check failed: "
                f"{quality_report.failure_summary()}"
            )
        connection.commit()
        emit_quality_report(quality_report, log)
        log.info(
            "Loaded %s GFF features into %s: %s genes, %s transcripts",
            source_config.name,
            db_name,
            len(annotation.genes),
            len(annotation.transcripts),
        )
        return {
            "genes": len(annotation.genes),
            "transcripts": len(annotation.transcripts),
            "cds_transcript_groups": len(annotation.cds_segments),
        }
    except Exception:
        connection.rollback()
        log.exception(
            "Failed loading %s GFF features into %s; transaction rolled back",
            source_config.name,
            db_name,
        )
        raise
    finally:
        connection.close()

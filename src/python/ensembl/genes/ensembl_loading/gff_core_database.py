"""Database schema and row helpers for Ensembl core loading."""

# These helpers mirror the Ensembl core schema operations and intentionally
# accept the complete set of values needed for each insert workflow.
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .core_utils.external_db import (
    REFSEQ_GENE_DB_NAME,
    REFSEQ_GENOMIC_DB_NAME,
    get_external_db_id,
    get_refseq_external_db_name,
    is_refseq_accession,
)
from .gff_annotation import DbCursor
from .gff_metadata import (
    parse_assembly_report_metadata,
    refseq_accession_db_token,
    species_db_name_token,
)
from .gff_models import ParsedAnnotation
from .gff_source_config import GENERIC_GFF_CONFIG, REFSEQ_CONFIG, GffSourceConfig

LOGGER = logging.getLogger(__name__)
DEFAULT_CORE_SCHEMA_SQL_PATH = (
    Path(__file__).resolve().parent / "config" / "core_schema.sql"
)
DEFAULT_ENSEMBL_RELEASE = "114"


def load_seq_regions_from_fna(
    fna_path: str | Path,
    cursor: DbCursor,
    coord_system_id: int,
    logger: logging.Logger | None = None,
    seq_region_attributes: Mapping[str, list[tuple[str, str]]] | None = None,
) -> dict[str, int]:
    """Load seq_region and dna rows from a converted FASTA file."""

    log = logger or LOGGER
    seq_region_ids: dict[str, int] = {}
    current_name: str | None = None
    sequence_lines: list[str] = []

    def store_sequence(name: str, sequence: str) -> None:
        cursor.execute(
            "SELECT seq_region_id FROM seq_region WHERE name = %s AND coord_system_id = %s",
            (name, coord_system_id),
        )
        row = cursor.fetchone()
        if row:
            seq_region_id = int(row[0])
            already_exists = True
        else:
            cursor.execute(
                "INSERT INTO seq_region (name, coord_system_id, length) VALUES (%s, %s, %s)",
                (name, coord_system_id, len(sequence)),
            )
            seq_region_id = int(cursor.lastrowid)
            already_exists = False

        if not already_exists:
            cursor.execute(
                "INSERT INTO dna (seq_region_id, sequence) VALUES (%s, %s)",
                (seq_region_id, sequence),
            )
            cursor.execute(
                "INSERT INTO seq_region_attrib (seq_region_id, attrib_type_id, value) "
                "VALUES (%s, %s, '')",
                (seq_region_id, get_or_create_attrib_type(cursor, "toplevel")),
            )

        if name.upper() in {"MT", "CHRM", "CHRMT"}:
            mitochondrial_attributes = [
                ("sequence_location", "mitochondrial_chromosome")
            ]
            mitochondrial_attributes.extend((seq_region_attributes or {}).get(name, []))
            for code, value in mitochondrial_attributes:
                attrib_type_id = get_or_create_attrib_type(cursor, code)
                cursor.execute(
                    "INSERT IGNORE INTO seq_region_attrib "
                    "(seq_region_id, attrib_type_id, value) VALUES (%s,%s,%s)",
                    (seq_region_id, attrib_type_id, value),
                )

        seq_region_ids[name] = seq_region_id

    with Path(fna_path).open("r", encoding="utf-8") as fasta_handle:
        for line in fasta_handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    store_sequence(current_name, "".join(sequence_lines))
                current_name = line[1:].split()[0]
                sequence_lines = []
            else:
                sequence_lines.append(line.upper())

        if current_name is not None:
            store_sequence(current_name, "".join(sequence_lines))

    log.info("Loaded %s seq_regions from %s", len(seq_region_ids), fna_path)
    return seq_region_ids


def insert_refseq_seq_region_synonyms(
    cursor: DbCursor,
    coord_system_id: int,
    accession_to_name: Mapping[str, str],
) -> int:
    """Insert original RefSeq accessions as synonyms for loaded seq_regions."""

    external_db_id = get_external_db_id(cursor, REFSEQ_GENOMIC_DB_NAME)
    inserted = 0
    for accession, sequence_name in accession_to_name.items():
        cursor.execute(
            "SELECT seq_region_id FROM seq_region "
            "WHERE name = %s AND coord_system_id = %s",
            (sequence_name, coord_system_id),
        )
        row = cursor.fetchone()
        if row is None:
            continue

        cursor.execute(
            "INSERT IGNORE INTO seq_region_synonym "
            "(seq_region_id, synonym, external_db_id) VALUES (%s, %s, %s)",
            (int(row[0]), accession, external_db_id),
        )
        inserted += 1

    LOGGER.info("Inserted %s RefSeq seq_region synonyms", inserted)
    return inserted


def replace_mitochondrial_features(cursor: DbCursor) -> None:
    """Remove existing MT features before a RefSeq mitochondrial reload."""

    cursor.execute(
        "DELETE ta FROM translation_attrib ta "
        "JOIN translation tl ON tl.translation_id = ta.translation_id "
        "JOIN transcript t ON t.transcript_id = tl.transcript_id "
        "JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id "
        "WHERE UPPER(sr.name) IN ('MT','CHRM','CHRMT')"
    )
    cursor.execute(
        "DELETE FROM translation WHERE transcript_id IN "
        "(SELECT t.transcript_id FROM transcript t "
        "JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id "
        "WHERE UPPER(sr.name) IN ('MT','CHRM','CHRMT'))"
    )
    cursor.execute(
        "DELETE FROM transcript_attrib WHERE transcript_id IN "
        "(SELECT t.transcript_id FROM transcript t "
        "JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id "
        "WHERE UPPER(sr.name) IN ('MT','CHRM','CHRMT'))"
    )
    cursor.execute(
        "DELETE FROM exon_transcript WHERE transcript_id IN "
        "(SELECT t.transcript_id FROM transcript t "
        "JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id "
        "WHERE UPPER(sr.name) IN ('MT','CHRM','CHRMT'))"
    )
    cursor.execute(
        "DELETE FROM transcript WHERE transcript_id IN "
        "(SELECT transcript_id FROM (SELECT t.transcript_id FROM transcript t "
        "JOIN seq_region sr ON sr.seq_region_id = t.seq_region_id "
        "WHERE UPPER(sr.name) IN ('MT','CHRM','CHRMT')) mt_transcripts)"
    )
    cursor.execute(
        "DELETE FROM gene WHERE seq_region_id IN "
        "(SELECT seq_region_id FROM seq_region "
        "WHERE UPPER(name) IN ('MT','CHRM','CHRMT'))"
    )
    cursor.execute(
        "DELETE e FROM exon e "
        "JOIN seq_region sr ON sr.seq_region_id = e.seq_region_id "
        "LEFT JOIN exon_transcript et ON et.exon_id = e.exon_id "
        "WHERE et.exon_id IS NULL "
        "AND UPPER(sr.name) IN ('MT','CHRM','CHRMT')"
    )


def core_schema_version(schema_sql_path: str | Path | None = None) -> str:
    """Return the bundled core schema version used in derived DB names."""

    schema_path = Path(schema_sql_path or DEFAULT_CORE_SCHEMA_SQL_PATH)
    try:
        schema_sql = schema_path.read_text(encoding="utf-8")
    except OSError:
        return DEFAULT_ENSEMBL_RELEASE

    match = re.search(
        r"['\"]schema_version['\"]\s*,\s*['\"]([^'\"]+)['\"]",
        schema_sql,
    )
    if match:
        return match.group(1)
    return DEFAULT_ENSEMBL_RELEASE


def assembly_version_db_token(assembly_accession: str) -> str:
    """Return the assembly version suffix used after the release number."""

    if "." not in assembly_accession:
        return "1"
    return re.sub(
        r"[^a-z0-9]+",
        "_",
        assembly_accession.rsplit(".", 1)[1].lower(),
    )


def derive_refseq_core_db_name(
    species_name: str,
    assembly_accession: str,
    ensembl_release: str | None = None,
) -> str:
    """Return the core DB name used for RefSeq core creation."""

    species_token = species_db_name_token(species_name)
    accession_token = refseq_accession_db_token(assembly_accession)
    release = ensembl_release or core_schema_version()
    return f"{species_token}_{accession_token}_rs_core_{release}_1"


def derive_core_db_name(
    species_name: str,
    assembly_accession: str,
    source_config: GffSourceConfig = GENERIC_GFF_CONFIG,
) -> str:
    """Return the core DB name for the selected loading source."""

    if source_config.name == REFSEQ_CONFIG.name:
        return derive_refseq_core_db_name(species_name, assembly_accession)

    species_token = species_db_name_token(species_name)
    assembly_token = assembly_accession.split("_", 1)[1].replace(".", "_")
    return f"{species_token}_core_{assembly_token}"


def load_schema_sql(cursor: DbCursor, schema_sql_path: str | Path) -> None:
    """Load a semicolon-delimited Ensembl schema SQL file."""

    raw_schema = Path(schema_sql_path).read_text(encoding="utf-8")
    clean_schema = re.sub(r"/\*\*.*?\*/", "", raw_schema, flags=re.DOTALL)
    for statement in filter(None, map(str.strip, clean_schema.split(";"))):
        cursor.execute(statement)


def resolve_schema_sql_path(schema_sql_path: str | Path | None) -> Path | None:
    """Return the schema SQL path to use for core creation.

    ``None`` means use the bundled default schema. An empty string disables
    schema loading explicitly.
    """

    if schema_sql_path == "":
        return None
    if schema_sql_path is None:
        return DEFAULT_CORE_SCHEMA_SQL_PATH
    return Path(schema_sql_path)


def connect_mysql(
    db_host: str,
    db_user: str,
    db_password: str,
    db_port: int,
    db_name: str | None = None,
) -> Any:
    """Connect to MySQL using the package-declared client dependency."""

    connection_args: dict[str, Any] = {
        "host": db_host,
        "user": db_user,
        "password": db_password,
        "port": db_port,
        "autocommit": False,
    }
    if db_name:
        connection_args["database"] = db_name

    try:
        import pymysql  # pylint: disable=import-outside-toplevel

        return pymysql.connect(**connection_args)
    except ImportError as pymysql_error:
        try:
            import mysql.connector  # pylint: disable=import-outside-toplevel
        except ImportError as mysql_connector_error:  # pragma: no cover
            raise ImportError(
                "A MySQL client is required for core loading. The package "
                "declares pymysql; reinstall the package in this environment "
                "or install pymysql. mysql-connector-python is also supported "
                "as a fallback."
            ) from mysql_connector_error

        LOGGER.warning(
            "pymysql is not installed; falling back to mysql.connector. "
            "Consider reinstalling the package dependencies."
        )
        del pymysql_error
        return mysql.connector.connect(**connection_args)


def get_coord_system_id(
    cursor: DbCursor,
    coord_system_id: int | None = None,
    coord_system_name: str = "primary_assembly",
    coord_system_version: str | None = None,
) -> int:
    """Resolve a coord_system_id from an explicit ID or name/version."""

    if coord_system_id is not None:
        return coord_system_id

    if coord_system_version is None:
        cursor.execute(
            "SELECT coord_system_id FROM coord_system WHERE name = %s "
            "ORDER BY rank LIMIT 1",
            (coord_system_name,),
        )
    else:
        cursor.execute(
            "SELECT coord_system_id FROM coord_system WHERE name = %s "
            "AND version = %s ORDER BY rank LIMIT 1",
            (coord_system_name, coord_system_version),
        )
    row = cursor.fetchone()
    if row is None:
        version_msg = (
            ""
            if coord_system_version is None
            else f" with version '{coord_system_version}'"
        )
        raise ValueError(
            f"Could not find coord_system '{coord_system_name}'{version_msg}"
        )
    return int(row[0])


def get_or_create_analysis(
    cursor: DbCursor,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> int:
    """Return an existing analysis_id or create one for the selected source."""

    cursor.execute(
        "SELECT analysis_id FROM analysis WHERE logic_name = %s LIMIT 1",
        (source_config.analysis_logic_name,),
    )
    row = cursor.fetchone()
    if row is not None:
        return int(row[0])

    cursor.execute(
        "INSERT INTO analysis (logic_name,created,program) VALUES (%s,NOW(),%s)",
        (source_config.analysis_logic_name, source_config.analysis_program),
    )
    return int(cursor.lastrowid)


def required_seq_region_names(annotation: ParsedAnnotation) -> set[str]:
    """Return all seq_region names referenced by parsed annotation records."""

    names = {gene.seq_name for gene in annotation.genes.values()}
    names.update(transcript.seq_name for transcript in annotation.transcripts.values())
    return names


def load_existing_seq_region_ids(
    cursor: DbCursor,
    annotation: ParsedAnnotation,
    coord_system_id: int,
) -> dict[str, int]:
    """Load seq_region IDs for all annotation seqids from an existing core DB."""

    seq_region_ids: dict[str, int] = {}
    missing_names: list[str] = []

    for seq_region_name in sorted(required_seq_region_names(annotation)):
        cursor.execute(
            "SELECT seq_region_id FROM seq_region WHERE name = %s AND coord_system_id = %s",
            (seq_region_name, coord_system_id),
        )
        row = cursor.fetchone()
        if row is None:
            missing_names.append(seq_region_name)
            continue
        seq_region_ids[seq_region_name] = int(row[0])

    if missing_names:
        preview = ", ".join(missing_names[:10])
        suffix = (
            "" if len(missing_names) <= 10 else f", ... ({len(missing_names)} total)"
        )
        raise ValueError(
            "GFF references seq_region names missing from the target core DB: "
            f"{preview}{suffix}"
        )

    return seq_region_ids


def next_table_id(cursor: DbCursor, table_name: str, id_column: str) -> int:
    """Return max(table.id)+1 for tables that require explicit IDs."""

    cursor.execute(f"SELECT COALESCE(MAX({id_column}), 0) + 1 FROM {table_name}")
    row = cursor.fetchone()
    if row is None:
        return 1
    return int(row[0])


def allocate_numeric_ids_from_core(
    cursor: DbCursor,
    annotation: ParsedAnnotation,
) -> tuple[dict[str, int], dict[str, int], int]:
    """Allocate non-conflicting IDs for loading into an existing core DB."""

    next_gene_id = next_table_id(cursor, "gene", "gene_id")
    next_transcript_id = next_table_id(cursor, "transcript", "transcript_id")
    first_exon_id = next_table_id(cursor, "exon", "exon_id")
    gene_id_map = {
        gene_id: index
        for index, gene_id in enumerate(sorted(annotation.genes), start=next_gene_id)
    }
    transcript_id_map = {
        transcript_id: index
        for index, transcript_id in enumerate(
            sorted(annotation.transcripts),
            start=next_transcript_id,
        )
    }
    return gene_id_map, transcript_id_map, first_exon_id


def initialise_core_tables(
    cursor: DbCursor,
    species_name: str,
    assembly_accession: str,
    assembly_report_path: str | Path = "",
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> tuple[int, int]:
    """Insert coord_system, meta, and analysis bootstrap rows."""

    assembly_metadata = parse_assembly_report_metadata(
        assembly_report_path,
        species_name,
        assembly_accession,
    )
    cursor.execute(
        "INSERT INTO coord_system (name,version,rank,attrib) "
        "VALUES ('primary_assembly','',1,'default_version,sequence_level')"
    )
    coord_system_id = int(cursor.lastrowid)
    cursor.executemany(
        "INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,%s,%s)",
        [
            (meta_key, meta_value)
            for meta_key, meta_value in assembly_metadata.items()
            if meta_value != ""
        ]
        + [
            ("genebuild.level", "toplevel"),
            ("transcriptbuild.level", "toplevel"),
            ("exonbuild.level", "toplevel"),
        ],
    )
    cursor.execute(
        "INSERT INTO analysis (logic_name,created,program) VALUES (%s,NOW(),%s)",
        (source_config.analysis_logic_name, source_config.analysis_program),
    )
    analysis_id = int(cursor.lastrowid)
    return coord_system_id, analysis_id


def allocate_numeric_ids(
    annotation: ParsedAnnotation,
) -> tuple[dict[str, int], dict[str, int]]:
    """Allocate deterministic numeric gene and transcript IDs."""

    gene_id_map = {
        gene_id: index
        for index, gene_id in enumerate(sorted(annotation.genes), start=1)
    }
    transcript_id_map = {
        transcript_id: index
        for index, transcript_id in enumerate(sorted(annotation.transcripts), start=1)
    }
    return gene_id_map, transcript_id_map


def insert_genes(
    cursor: DbCursor,
    annotation: ParsedAnnotation,
    seq_region_ids: Mapping[str, int],
    gene_id_map: Mapping[str, int],
    transcript_id_map: Mapping[str, int],
    analysis_id: int,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> None:
    """Insert gene rows into the core database."""

    first_transcript_by_gene: dict[str, str] = {}
    for transcript_id, transcript in annotation.transcripts.items():
        first_transcript_by_gene.setdefault(transcript.gene_id, transcript_id)

    for gene_id, gene in annotation.genes.items():
        seq_region_id = seq_region_ids[gene.seq_name]
        first_transcript = first_transcript_by_gene.get(gene_id)
        canonical_transcript_id = (
            transcript_id_map.get(first_transcript) if first_transcript else None
        )
        cursor.execute(
            """INSERT INTO gene
               (gene_id, seq_region_id, seq_region_start, seq_region_end,
                seq_region_strand, biotype, analysis_id, stable_id,
                display_xref_id, source, canonical_transcript_id, description)
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,NULL,%s,%s,%s)""",
            (
                gene_id_map[gene_id],
                seq_region_id,
                gene.start,
                gene.end,
                gene.strand,
                gene.biotype,
                analysis_id,
                gene.stable_id,
                source_config.source_label,
                canonical_transcript_id,
                gene.description,
            ),
        )
        if source_config.name != REFSEQ_CONFIG.name:
            continue
        external_db_id = get_external_db_id(cursor, REFSEQ_GENE_DB_NAME)
        gene_db_id = gene_id_map[gene_id]
        xref_accession = gene.xref_geneid or gene.stable_id
        cursor.execute(
            """INSERT INTO xref
               (external_db_id, dbprimary_acc, display_label, version,
                description, info_type, info_text)
               VALUES (%s,%s,%s,NULL,NULL,'DIRECT','')""",
            (
                external_db_id,
                xref_accession,
                gene.name,
            ),
        )
        xref_id = cursor.lastrowid
        cursor.execute(
            """INSERT INTO object_xref
               (ensembl_id, ensembl_object_type, xref_id, linkage_annotation,
                analysis_id)
               VALUES (%s,'Gene',%s,NULL,%s)""",
            (gene_db_id, xref_id, analysis_id),
        )
        cursor.execute(
            "UPDATE gene SET display_xref_id = %s WHERE gene_id = %s",
            (xref_id, gene_db_id),
        )


def get_or_create_attrib_type(cursor: DbCursor, code: str) -> int:
    """Return the core attribute type ID for an imported annotation code."""

    cursor.execute(
        "SELECT attrib_type_id FROM attrib_type WHERE code = %s LIMIT 1",
        (code,),
    )
    row = cursor.fetchone()
    if row is not None:
        return int(row[0])
    cursor.execute(
        "INSERT INTO attrib_type (code, name, description) VALUES (%s,%s,%s)",
        (code, code, "Imported from RefSeq annotation exceptions"),
    )
    return cursor.lastrowid


def insert_attributes(
    cursor: DbCursor,
    table: str,
    object_column: str,
    object_id: int,
    attributes: list[tuple[str, str]],
) -> None:
    """Insert core attributes for one transcript or translation."""

    for code, value in attributes:
        attrib_type_id = get_or_create_attrib_type(cursor, code)
        cursor.execute(
            f"INSERT INTO {table} ({object_column}, attrib_type_id, value) "
            "VALUES (%s,%s,%s)",
            (object_id, attrib_type_id, value),
        )


def insert_refseq_object_xref(
    cursor: DbCursor,
    accession: str,
    object_id: int,
    object_type: str,
    analysis_id: int | None,
) -> int:
    """Insert and link one RefSeq xref for a core object."""
    db_name = get_refseq_external_db_name(accession, object_type)
    external_db_id = get_external_db_id(cursor, db_name)
    cursor.execute(
        """INSERT INTO xref
           (external_db_id, dbprimary_acc, display_label, version,
            description, info_type, info_text)
           VALUES (%s,%s,%s,NULL,NULL,'DIRECT','')""",
        (external_db_id, accession, accession),
    )
    xref_id = cursor.lastrowid
    cursor.execute(
        """INSERT INTO object_xref
           (ensembl_id, ensembl_object_type, xref_id, linkage_annotation,
            analysis_id)
           VALUES (%s,%s,%s,NULL,%s)""",
        (object_id, object_type, xref_id, analysis_id),
    )
    return xref_id


def _reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA string."""

    return sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def frameshift_transcript_attributes(
    cursor: DbCursor,
    transcript: Any,
    cds_segments: list[Any],
    seq_region_id: int,
) -> list[tuple[str, str]]:
    """Create transcript SeqEdits for RefSeq ribosomal slippage intervals."""

    if not transcript.frameshift_events:
        return []
    strand = transcript.strand
    ordered_cds = sorted(
        cds_segments,
        key=lambda cds: cds.end if strand == -1 else cds.start,
        reverse=(strand == -1),
    )
    transcript_coordinates: dict[int, int] = {}
    transcript_position = 1
    for exon in sorted(
        transcript.exons,
        key=lambda item: item.end if strand == -1 else item.start,
        reverse=(strand == -1),
    ):
        coordinates = (
            range(exon.end, exon.start - 1, -1)
            if strand == -1
            else range(exon.start, exon.end + 1)
        )
        for coordinate in coordinates:
            transcript_coordinates[coordinate] = transcript_position
            transcript_position += 1

    cursor.execute(
        "SELECT sequence FROM dna WHERE seq_region_id = %s", (seq_region_id,)
    )
    row = cursor.fetchone()
    dna = str(row[0]) if row else ""
    attributes: list[tuple[str, str]] = []
    for _exception, note in transcript.frameshift_events:
        match = re.search(r"([+-])(\d+)", note)
        direction = match.group(1) if match else "-"
        length = int(match.group(2)) if match else 1
        for current, following in zip(ordered_cds, ordered_cds[1:]):
            if strand == 1:
                gap_start = current.end + 1
                gap_end = following.start - 1
                overlap_start = following.start
                overlap_end = current.end
            else:
                gap_start = following.end + 1
                gap_end = current.start - 1
                overlap_start = current.start
                overlap_end = following.end

            if direction == "+" and gap_end - gap_start + 1 == length:
                positions = [
                    transcript_coordinates.get(coordinate)
                    for coordinate in range(gap_start, gap_end + 1)
                ]
                mapped_positions = [
                    position for position in positions if position is not None
                ]
                if len(mapped_positions) == len(positions):
                    attributes.append(
                        (
                            "_rna_edit",
                            f"{min(mapped_positions)} {max(mapped_positions)} ",
                        )
                    )
            elif direction == "-" and overlap_end - overlap_start + 1 == length:
                position = transcript_coordinates.get(overlap_start)
                if position is not None and dna:
                    sequence = dna[overlap_start - 1 : overlap_end]
                    if strand == -1:
                        sequence = _reverse_complement(sequence)
                    attributes.append(
                        ("_rna_edit", f"{position} {position} {sequence}")
                    )
    return attributes


def insert_transcripts_and_exons(
    cursor: DbCursor,
    annotation: ParsedAnnotation,
    seq_region_ids: Mapping[str, int],
    gene_id_map: Mapping[str, int],
    transcript_id_map: Mapping[str, int],
    analysis_id: int,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
    first_exon_id: int = 1,
    deduplicate_exons: bool = True,
) -> tuple[dict[tuple[str, int], int], dict[str, dict[tuple[int, int], int]]]:
    """Insert transcripts, exons, and exon_transcript links."""

    exon_id_map: dict[tuple[str, int], int] = {}
    per_transcript_coord_to_exon_id: dict[str, dict[tuple[int, int], int]] = {}
    shared_exon_ids: dict[tuple[str, int, int, int, int, int, int, str | None], int] = (
        {}
    )
    next_exon_id = first_exon_id

    for transcript_id, transcript in annotation.transcripts.items():
        seq_region_id = seq_region_ids[transcript.seq_name]
        cursor.execute(
            """INSERT INTO transcript
               (transcript_id, gene_id, seq_region_id, seq_region_start, seq_region_end,
                seq_region_strand, biotype, analysis_id, stable_id, display_xref_id, source)
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,NULL,%s)""",
            (
                transcript_id_map[transcript_id],
                gene_id_map[transcript.gene_id],
                seq_region_id,
                transcript.start,
                transcript.end,
                transcript.strand,
                transcript.biotype,
                analysis_id,
                transcript.stable_id,
                source_config.source_label,
            ),
        )
        if source_config.name == REFSEQ_CONFIG.name and is_refseq_accession(
            transcript.stable_id, "Transcript"
        ):
            transcript_xref_id = insert_refseq_object_xref(
                cursor,
                transcript.stable_id,
                transcript_id_map[transcript_id],
                "Transcript",
                analysis_id,
            )
            cursor.execute(
                "UPDATE transcript SET display_xref_id = %s WHERE transcript_id = %s",
                (transcript_xref_id, transcript_id_map[transcript_id]),
            )
        frameshift_attributes = frameshift_transcript_attributes(
            cursor,
            transcript,
            annotation.cds_segments.get(transcript_id, []),
            seq_region_id,
        )
        insert_attributes(
            cursor,
            "transcript_attrib",
            "transcript_id",
            transcript_id_map[transcript_id],
            [*transcript.transcript_attributes, *frameshift_attributes],
        )

        sorted_exons = sorted(
            transcript.exons,
            key=lambda exon: exon.start,
            reverse=(transcript.strand == -1),
        )
        coord_to_exon_id: dict[tuple[int, int], int] = {}
        for rank, exon in enumerate(sorted_exons, start=1):
            phase_value = -1 if exon.phase is None else exon.phase
            end_phase_value = -1 if exon.end_phase is None else exon.end_phase
            exon_key = (
                transcript.gene_id,
                seq_region_id,
                exon.start,
                exon.end,
                exon.strand,
                phase_value,
                end_phase_value,
                exon.stable_id,
            )
            exon_id = shared_exon_ids.get(exon_key) if deduplicate_exons else None
            if exon_id is None:
                exon_id = next_exon_id
                next_exon_id += 1
                cursor.execute(
                    """INSERT INTO exon
                       (exon_id, seq_region_id, seq_region_start, seq_region_end,
                        seq_region_strand, phase, end_phase, stable_id)
                       VALUES (%s,%s,%s,%s,%s,%s,%s,%s)""",
                    (
                        exon_id,
                        seq_region_id,
                        exon.start,
                        exon.end,
                        exon.strand,
                        phase_value,
                        end_phase_value,
                        exon.stable_id,
                    ),
                )
                if deduplicate_exons:
                    shared_exon_ids[exon_key] = exon_id
            exon_id_map[(transcript_id, rank)] = exon_id
            coord_to_exon_id[exon.coordinate_key] = exon_id
            cursor.execute(
                "INSERT INTO exon_transcript (exon_id, transcript_id, rank) "
                "VALUES (%s,%s,%s)",
                (
                    exon_id,
                    transcript_id_map[transcript_id],
                    rank,
                ),
            )

        per_transcript_coord_to_exon_id[transcript_id] = coord_to_exon_id

    return exon_id_map, per_transcript_coord_to_exon_id


def insert_translations(
    cursor: DbCursor,
    annotation: ParsedAnnotation,
    transcript_id_map: Mapping[str, int],
    exon_id_map: Mapping[tuple[str, int], int],
    analysis_id: int | None = None,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> None:
    """Insert translation rows and canonical translation links."""

    for transcript_id, cds_list in annotation.cds_segments.items():
        db_transcript_id = transcript_id_map.get(transcript_id)
        if not db_transcript_id or transcript_id not in annotation.transcripts:
            continue

        transcript = annotation.transcripts[transcript_id]
        strand = transcript.strand
        exons = transcript.exons
        if not cds_list or not exons:
            continue

        if strand == 1:
            translation_start_pos = min(cds.start for cds in cds_list)
            translation_end_pos = max(cds.end for cds in cds_list)
        else:
            translation_start_pos = max(cds.end for cds in cds_list)
            translation_end_pos = min(cds.start for cds in cds_list)

        genomic_sorted_exons = sorted(exons, key=lambda exon: exon.start)
        inserted_order = sorted(
            exons,
            key=lambda exon: exon.start,
            reverse=(strand == -1),
        )
        rank_map = {
            exon.coordinate_key: rank
            for rank, exon in enumerate(inserted_order, start=1)
        }

        start_rank = end_rank = None
        start_offset = end_offset = None
        for exon in genomic_sorted_exons:
            if exon.start <= translation_start_pos <= exon.end:
                start_offset = (
                    translation_start_pos - exon.start + 1
                    if strand == 1
                    else exon.end - translation_start_pos + 1
                )
                start_rank = rank_map[exon.coordinate_key]
            if exon.start <= translation_end_pos <= exon.end:
                end_offset = (
                    translation_end_pos - exon.start + 1
                    if strand == 1
                    else exon.end - translation_end_pos + 1
                )
                end_rank = rank_map[exon.coordinate_key]

        if start_rank and end_rank and start_offset and end_offset:
            protein_id = transcript.protein_id or f"{transcript_id}_prot"
            cursor.execute(
                """INSERT INTO translation
                     (transcript_id, start_exon_id, end_exon_id,
                      seq_start, seq_end, stable_id)
                   VALUES (%s,%s,%s,%s,%s,%s)""",
                (
                    db_transcript_id,
                    exon_id_map[(transcript_id, start_rank)],
                    exon_id_map[(transcript_id, end_rank)],
                    start_offset,
                    end_offset,
                    protein_id,
                ),
            )
            translation_id = cursor.lastrowid
            if (
                source_config.name == REFSEQ_CONFIG.name
                and transcript.protein_id
                and is_refseq_accession(transcript.protein_id, "Translation")
            ):
                insert_refseq_object_xref(
                    cursor,
                    transcript.protein_id,
                    translation_id,
                    "Translation",
                    analysis_id,
                )
            cursor.execute(
                "UPDATE transcript SET canonical_translation_id = %s WHERE transcript_id = %s",
                (translation_id, db_transcript_id),
            )
            insert_attributes(
                cursor,
                "translation_attrib",
                "translation_id",
                translation_id,
                transcript.translation_attributes,
            )

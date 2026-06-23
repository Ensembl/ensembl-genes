"""Load anno pipeline single-line repeat/simple GTF features into a core DB."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

try:  # Support both package imports and direct same-directory imports.
    from .gff_core_loader import (
        DbCursor,
        connect_mysql,
        get_coord_system_id,
        open_text_maybe_gzip,
        parse_gtf_attributes,
    )
except ImportError:  # pragma: no cover - used when run beside this file.
    from gff_core_loader import (  # type: ignore
        DbCursor,
        connect_mysql,
        get_coord_system_id,
        open_text_maybe_gzip,
        parse_gtf_attributes,
    )


LOGGER = logging.getLogger(__name__)

REPEAT_TYPES = {
    "dust": "Dust",
    "repeatdetector": "repeatdetector",
    "repeatmask_repbase_human": "repeatmasker",
    "trf": "Tandem repeats",
}
SIMPLE_FEATURE_ANALYSES = frozenset({"cpg", "eponine"})


@dataclass(frozen=True)
class SingleLineFeatureRecord:
    """One single-line feature row from an anno pipeline GTF."""

    seq_name: str
    start: int
    end: int
    strand: int
    attributes: Mapping[str, str]


def available_single_line_feature_analyses() -> tuple[str, ...]:
    """Return analysis names supported by the old single_line_feature loader."""

    return tuple(sorted((*REPEAT_TYPES, *SIMPLE_FEATURE_ANALYSES)))


def parse_single_line_gtf(gtf_path: str | Path) -> list[SingleLineFeatureRecord]:
    """Parse single-line repeat/simple features from an anno pipeline GTF."""

    records: list[SingleLineFeatureRecord] = []
    with open_text_maybe_gzip(gtf_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) != 9:
                continue

            strand = columns[6]
            if strand == "+":
                strand_value = 1
            elif strand == "-":
                strand_value = -1
            else:
                raise ValueError(
                    f"{gtf_path}:{line_number}: unsupported strand '{strand}'"
                )

            records.append(
                SingleLineFeatureRecord(
                    seq_name=columns[0],
                    start=int(columns[3]),
                    end=int(columns[4]),
                    strand=strand_value,
                    attributes=parse_gtf_attributes(columns[8]),
                )
            )
    return records


def load_seq_region_ids_for_feature_names(
    cursor: DbCursor,
    records: list[SingleLineFeatureRecord],
    coord_system_id: int,
) -> dict[str, int]:
    """Resolve seq_region IDs for all feature seq_names."""

    seq_region_ids: dict[str, int] = {}
    for seq_name in sorted({record.seq_name for record in records}):
        cursor.execute(
            "SELECT seq_region_id FROM seq_region WHERE name = %s "
            "AND coord_system_id = %s",
            (seq_name, coord_system_id),
        )
        row = cursor.fetchone()
        if row is None:
            raise ValueError(
                f"Could not find seq_region '{seq_name}' "
                f"in coord_system_id {coord_system_id}"
            )
        seq_region_ids[seq_name] = int(row[0])
    return seq_region_ids


def get_or_create_single_line_analysis(cursor: DbCursor, analysis_name: str) -> int:
    """Return or create the analysis row used by repeat/simple feature loads."""

    cursor.execute(
        "SELECT analysis_id FROM analysis WHERE logic_name = %s LIMIT 1",
        (analysis_name,),
    )
    row = cursor.fetchone()
    if row is not None:
        return int(row[0])

    cursor.execute(
        "INSERT INTO analysis (logic_name,created,module) VALUES (%s,NOW(),%s)",
        (analysis_name, "Anno"),
    )
    return int(cursor.lastrowid)


def get_or_create_repeat_consensus(
    cursor: DbCursor,
    repeat_name: str,
    repeat_class: str,
    repeat_type: str,
    repeat_consensus: str,
) -> int:
    """Return or create a repeat_consensus row."""

    cursor.execute(
        "SELECT repeat_consensus_id FROM repeat_consensus "
        "WHERE repeat_name = %s AND repeat_class = %s "
        "AND repeat_type = %s AND repeat_consensus = %s LIMIT 1",
        (repeat_name, repeat_class, repeat_type, repeat_consensus),
    )
    row = cursor.fetchone()
    if row is not None:
        return int(row[0])

    cursor.execute(
        "INSERT INTO repeat_consensus "
        "(repeat_name, repeat_class, repeat_type, repeat_consensus) "
        "VALUES (%s,%s,%s,%s)",
        (repeat_name, repeat_class, repeat_type, repeat_consensus),
    )
    return int(cursor.lastrowid)


def load_repeat_records(
    cursor: DbCursor,
    records: list[SingleLineFeatureRecord],
    seq_region_ids: Mapping[str, int],
    analysis_id: int,
    analysis_name: str,
) -> tuple[int, int]:
    """Insert repeat_consensus and repeat_feature rows."""

    consensus_ids: dict[tuple[str, str, str, str], int] = {}
    for record in records:
        repeat_name = record.attributes.get("repeat_name", analysis_name)
        repeat_class = record.attributes.get("repeat_class", analysis_name)
        repeat_type = record.attributes.get(
            "repeat_type",
            REPEAT_TYPES.get(analysis_name, analysis_name),
        )
        repeat_consensus = record.attributes.get("repeat_consensus", "N")
        repeat_score = float(record.attributes.get("score", "0") or 0)
        consensus_key = (repeat_name, repeat_class, repeat_type, repeat_consensus)
        repeat_consensus_id = consensus_ids.get(consensus_key)
        if repeat_consensus_id is None:
            repeat_consensus_id = get_or_create_repeat_consensus(
                cursor,
                repeat_name,
                repeat_class,
                repeat_type,
                repeat_consensus,
            )
            consensus_ids[consensus_key] = repeat_consensus_id

        cursor.execute(
            "INSERT INTO repeat_feature "
            "(seq_region_id, seq_region_start, seq_region_end, "
            "seq_region_strand, repeat_start, repeat_end, "
            "repeat_consensus_id, analysis_id, score) "
            "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)",
            (
                seq_region_ids[record.seq_name],
                record.start,
                record.end,
                record.strand,
                1,
                record.end - record.start + 1,
                repeat_consensus_id,
                analysis_id,
                repeat_score,
            ),
        )
    return len(records), len(consensus_ids)


def load_simple_records(
    cursor: DbCursor,
    records: list[SingleLineFeatureRecord],
    seq_region_ids: Mapping[str, int],
    analysis_id: int,
) -> int:
    """Insert simple_feature rows."""

    for record in records:
        cursor.execute(
            "INSERT INTO simple_feature "
            "(seq_region_id, seq_region_start, seq_region_end, "
            "seq_region_strand, display_label, analysis_id, score) "
            "VALUES (%s,%s,%s,%s,%s,%s,%s)",
            (
                seq_region_ids[record.seq_name],
                record.start,
                record.end,
                record.strand,
                "",
                analysis_id,
                0.0,
            ),
        )
    return len(records)


def load_single_line_features_to_core(  # pylint: disable=too-many-arguments,too-many-locals
    gtf_path: str | Path,
    analysis_name: str,
    db_name: str,
    db_host: str,
    db_user: str,
    db_password: str,
    db_port: int,
    coord_system_id: int | None = None,
    coord_system_name: str = "primary_assembly",
    coord_system_version: str | None = None,
    logger: logging.Logger | None = None,
) -> dict[str, int]:
    """Load repeat/simple single-line GTF features into an existing core DB."""

    log = logger or LOGGER
    if analysis_name not in available_single_line_feature_analyses():
        available = ", ".join(available_single_line_feature_analyses())
        raise ValueError(
            f"Unsupported single-line analysis '{analysis_name}'. "
            f"Available analyses: {available}"
        )

    records = parse_single_line_gtf(gtf_path)
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
        analysis_id = get_or_create_single_line_analysis(cursor, analysis_name)
        seq_region_ids = load_seq_region_ids_for_feature_names(
            cursor,
            records,
            resolved_coord_system_id,
        )
        repeat_features = 0
        repeat_consensus = 0
        simple_features = 0
        if analysis_name in REPEAT_TYPES:
            repeat_features, repeat_consensus = load_repeat_records(
                cursor,
                records,
                seq_region_ids,
                analysis_id,
                analysis_name,
            )
        else:
            simple_features = load_simple_records(
                cursor,
                records,
                seq_region_ids,
                analysis_id,
            )

        connection.commit()
        log.info(
            "Loaded %s single-line features into %s: "
            "%s repeat features, %s simple features",
            analysis_name,
            db_name,
            repeat_features,
            simple_features,
        )
        return {
            "repeat_features": repeat_features,
            "repeat_consensus": repeat_consensus,
            "simple_features": simple_features,
        }
    except Exception:
        connection.rollback()
        log.exception(
            "Failed loading %s single-line features into %s; transaction rolled back",
            analysis_name,
            db_name,
        )
        raise
    finally:
        connection.close()

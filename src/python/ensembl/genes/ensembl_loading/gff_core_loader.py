"""Load converted GFF3 annotations into an Ensembl-style core database."""

from __future__ import annotations

import gzip
import logging
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Mapping, Protocol, Sequence, TextIO

try:  # Support both package imports and direct same-directory imports.
    from .gff_models import (
        CdsSegment,
        ExonRecord,
        GeneRecord,
        ParsedAnnotation,
        TranscriptRecord,
    )
    from .gff_quality_check import emit_quality_report, run_core_load_quality_check
    from .gff_source_config import GENERIC_GFF_CONFIG, REFSEQ_CONFIG, GffSourceConfig
except ImportError:  # pragma: no cover - used when run beside this file.
    from gff_models import (  # type: ignore
        CdsSegment,
        ExonRecord,
        GeneRecord,
        ParsedAnnotation,
        TranscriptRecord,
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
DEFAULT_CORE_SCHEMA_SQL_PATH = (
    Path(__file__).resolve().parent / "config" / "core_schema.sql"
)
DEFAULT_ENSEMBL_RELEASE = "114"


@contextmanager
def open_text_maybe_gzip(path: str | Path) -> Iterator[TextIO]:
    """Open plain-text or gzip-compressed GFF3 files for text reading."""

    input_path = Path(path)
    if input_path.suffix == ".gz":
        with gzip.open(input_path, "rt") as handle:
            yield handle
    else:
        with input_path.open("r") as handle:
            yield handle


class DbCursor(Protocol):
    """Small cursor protocol used by the loader and tests."""

    lastrowid: int

    def execute(self, operation: str, params: Any | None = None) -> Any:
        """Execute one SQL statement."""
        ...

    def executemany(self, operation: str, seq_params: Sequence[Any]) -> Any:
        """Execute one SQL statement for many parameter sets."""
        ...

    def fetchone(self) -> tuple[Any, ...] | None:
        """Fetch one result row."""
        ...


def normalize_id(
    raw_id: str,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> str:
    """Strip configured source prefixes to create stable internal IDs."""

    normalized_id = raw_id
    for prefix in source_config.id_prefixes_to_strip:
        normalized_id = normalized_id.replace(prefix, "")
    return normalized_id


def parent_id(
    raw_parent: str,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> str:
    """Return the normalized Parent ID using the configured source rule."""

    return normalize_id(raw_parent, source_config)


def first_existing_attribute(
    attributes: Mapping[str, str],
    attribute_names: Sequence[str],
) -> str | None:
    """Return the first non-empty configured feature attribute value."""

    for attribute_name in attribute_names:
        value = attributes.get(attribute_name)
        if value:
            return value
    return None


def parse_gff3_attributes(raw_attributes: str) -> dict[str, str]:
    """Parse a GFF3 attribute column into a key-value dictionary."""

    return {
        key: value
        for key, value in (
            item.split("=", 1) for item in raw_attributes.split(";") if "=" in item
        )
    }


def parse_gtf_attributes(raw_attributes: str) -> dict[str, str]:
    """Parse a GTF attribute column into a key-value dictionary."""

    attributes: dict[str, str] = {}
    for item in raw_attributes.rstrip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            attributes[item] = "1"
            continue

        key, raw_value = item.split(None, 1)
        value = raw_value.strip()
        if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
            value = value[1:-1]
        attributes[key] = value
    return attributes


def parse_feature_attributes(
    raw_attributes: str,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> dict[str, str]:
    """Parse a feature attribute column using the selected source syntax."""

    if source_config.attribute_format == "gff3":
        return parse_gff3_attributes(raw_attributes)
    if source_config.attribute_format == "gtf":
        return parse_gtf_attributes(raw_attributes)
    raise ValueError(
        f"Unsupported attribute format '{source_config.attribute_format}' "
        f"for source config '{source_config.name}'"
    )


def required_attribute(
    attributes: Mapping[str, str],
    attribute_name: str,
    feature_type: str,
    line_number: int,
    path: Path,
) -> str:
    """Return a required attribute or raise a parse error with context."""

    value = attributes.get(attribute_name)
    if not value:
        raise ValueError(
            f"{path}:{line_number}: {feature_type} row is missing required "
            f"attribute '{attribute_name}'"
        )
    return value


def resolve_biotype(
    feature_type: str,
    attributes: Mapping[str, str],
    feature_id: str,
    logger: logging.Logger | None = None,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> str:
    """Resolve an Ensembl-compatible biotype from a configured GFF3 source."""

    log = logger or LOGGER
    gbkey = attributes.get(source_config.gbkey_attribute, "")
    pseudo_flag = attributes.get(source_config.pseudo_attribute, "false") == "true"
    transcript_biotype = attributes.get(source_config.transcript_biotype_attribute)
    gene_biotype = attributes.get(source_config.gene_biotype_attribute)

    if gbkey.endswith(source_config.segment_gbkey_suffix):
        segment_type = gbkey[0]
        return f"{source_config.segment_biotype_prefix}_{segment_type}_gene"

    if source_config.transcribed_pseudogene_gbkey_token in gbkey:
        return "transcribed_pseudogene"

    if pseudo_flag:
        return "pseudogene"

    if transcript_biotype:
        return transcript_biotype

    if source_config.translation_coords_attribute and attributes.get(
        source_config.translation_coords_attribute
    ):
        return "protein_coding"

    if feature_type in source_config.biotype_transcript_feature_types:
        return source_config.transcript_feature_biotype_map.get(
            feature_type,
            feature_type,
        )

    if feature_type == "exon":
        exon_biotype = source_config.exon_gbkey_biotype_map.get(gbkey)
        if exon_biotype:
            return exon_biotype

    if gene_biotype:
        return gene_biotype

    log.debug(
        "Defaulting biotype to %s for %s %s %s",
        source_config.default_biotype,
        feature_id,
        feature_type,
        gbkey,
    )
    return source_config.default_biotype


def update_gene_from_transcript(
    annotation: ParsedAnnotation,
    gene_id: str,
    seq_name: str,
    start: int,
    end: int,
    strand: int,
    biotype: str,
    logger: logging.Logger,
) -> None:
    """Create or expand a gene record from transcript-level GTF data."""

    existing_gene = annotation.genes.get(gene_id)
    if existing_gene is None:
        annotation.genes[gene_id] = GeneRecord(
            seq_name=seq_name,
            start=start,
            end=end,
            strand=strand,
            biotype=biotype,
            stable_id=gene_id,
            name=gene_id,
        )
        return

    if existing_gene.seq_name != seq_name or existing_gene.strand != strand:
        logger.warning(
            "Gene %s has transcript rows on inconsistent seq_region/strand: "
            "%s:%s and %s:%s",
            gene_id,
            existing_gene.seq_name,
            existing_gene.strand,
            seq_name,
            strand,
        )
        return

    existing_gene.start = min(existing_gene.start, start)
    existing_gene.end = max(existing_gene.end, end)
    if existing_gene.biotype == "not_set" or biotype == "protein_coding":
        existing_gene.biotype = biotype


def parse_translation_coords(raw_translation_coords: str) -> tuple[int, ...]:
    """Parse anno GTF translation_coords into integer components."""

    match = re.fullmatch(
        r"(\d+):(\d+):(\d+):(\d+):(\d+):(\d+)",
        raw_translation_coords,
    )
    if not match:
        raise ValueError(
            f"Could not parse translation_coords value: {raw_translation_coords!r}"
        )
    return tuple(int(value) for value in match.groups())


def cds_segments_from_translation_coords(
    transcript_id: str,
    transcript: TranscriptRecord,
) -> list[CdsSegment]:
    """Build CDS segments from anno GTF transcript translation coordinates."""

    if not transcript.translation_coords:
        return []

    (
        start_exon_start,
        start_exon_end,
        start_exon_offset,
        end_exon_start,
        end_exon_end,
        end_exon_offset,
    ) = parse_translation_coords(transcript.translation_coords)

    exons_in_transcript_order = sorted(
        transcript.exons,
        key=lambda exon: exon.start,
        reverse=(transcript.strand == -1),
    )
    start_index = next(
        (
            index
            for index, exon in enumerate(exons_in_transcript_order)
            if exon.start == start_exon_start and exon.end == start_exon_end
        ),
        None,
    )
    end_index = next(
        (
            index
            for index, exon in enumerate(exons_in_transcript_order)
            if exon.start == end_exon_start and exon.end == end_exon_end
        ),
        None,
    )
    if start_index is None or end_index is None:
        raise ValueError(
            f"Could not match translation_coords to exons for transcript "
            f"{transcript_id}: {transcript.translation_coords}"
        )
    if start_index > end_index:
        raise ValueError(
            f"translation_coords start exon occurs after end exon for "
            f"transcript {transcript_id}: {transcript.translation_coords}"
        )

    cds_segments: list[CdsSegment] = []
    for index, exon in enumerate(
        exons_in_transcript_order[start_index : end_index + 1],
        start=start_index,
    ):
        cds_start = exon.start
        cds_end = exon.end
        if index == start_index:
            if transcript.strand == 1:
                cds_start = exon.start + start_exon_offset - 1
            else:
                cds_end = exon.end - start_exon_offset + 1
        if index == end_index:
            if transcript.strand == 1:
                cds_end = exon.start + end_exon_offset - 1
            else:
                cds_start = exon.end - end_exon_offset + 1

        if not (exon.start <= cds_start <= cds_end <= exon.end):
            raise ValueError(
                f"translation_coords produced CDS outside exon bounds for "
                f"transcript {transcript_id}: {transcript.translation_coords}"
            )
        cds_segments.append(
            CdsSegment(
                start=cds_start,
                end=cds_end,
                strand=transcript.strand,
                phase="0" if not cds_segments else ".",
            )
        )
    return cds_segments


def synthesize_cds_from_translation_coords(
    annotation: ParsedAnnotation,
    logger: logging.Logger,
) -> None:
    """Populate CDS records for GTF transcripts carrying translation_coords."""

    synthesized = 0
    for transcript_id, transcript in annotation.transcripts.items():
        if (
            not transcript.translation_coords
            or transcript_id in annotation.cds_segments
        ):
            continue
        cds_segments = cds_segments_from_translation_coords(transcript_id, transcript)
        if cds_segments:
            annotation.cds_segments[transcript_id] = cds_segments
            synthesized += 1

    if synthesized:
        logger.info("Synthesized CDS segments for %s transcripts", synthesized)


def parse_converted_gff3(
    converted_gff_path: str | Path,
    logger: logging.Logger | None = None,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> ParsedAnnotation:
    """Parse GFF3/GTF into reusable in-memory annotation records."""

    log = logger or LOGGER
    annotation = ParsedAnnotation()
    gff_path = Path(converted_gff_path)

    with open_text_maybe_gzip(gff_path) as gff_handle:
        for line_number, line in enumerate(gff_handle, start=1):
            if line.startswith("#") or not line.strip():
                continue

            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                log.debug("Skipping non-feature row %s in %s", line_number, gff_path)
                continue

            (
                seq_name,
                _source,
                feature_type,
                start_raw,
                end_raw,
                _score,
                strand,
                phase,
                attrs,
            ) = columns
            start = int(start_raw)
            end = int(end_raw)

            strand_value = 1 if strand == "+" else -1
            attributes = parse_feature_attributes(attrs, source_config)

            if feature_type in source_config.parsed_gene_feature_types:
                gene_id = normalize_id(
                    required_attribute(
                        attributes,
                        source_config.gene_id_attribute,
                        feature_type,
                        line_number,
                        gff_path,
                    ),
                    source_config,
                )
                dbxref_geneid = next(
                    (
                        value.split(":", 1)[1]
                        for value in attributes.get("Dbxref", "").split(",")
                        if value.startswith(source_config.gene_xref_prefix)
                    ),
                    None,
                )
                gene_name = (
                    first_existing_attribute(
                        attributes,
                        source_config.gene_name_attributes,
                    )
                    or gene_id
                )
                annotation.genes[gene_id] = GeneRecord(
                    seq_name=seq_name,
                    start=start,
                    end=end,
                    strand=strand_value,
                    biotype=resolve_biotype(
                        feature_type,
                        attributes,
                        gene_id,
                        log,
                        source_config=source_config,
                    ),
                    stable_id=gene_id,
                    xref_geneid=dbxref_geneid,
                    name=gene_name,
                )
                continue

            if feature_type in source_config.parsed_transcript_feature_types:
                transcript_id = normalize_id(
                    required_attribute(
                        attributes,
                        source_config.transcript_id_attribute,
                        feature_type,
                        line_number,
                        gff_path,
                    ),
                    source_config,
                )
                gene_id = parent_id(
                    required_attribute(
                        attributes,
                        source_config.parent_gene_attribute,
                        feature_type,
                        line_number,
                        gff_path,
                    ),
                    source_config,
                )
                stable_id = (
                    first_existing_attribute(
                        attributes,
                        source_config.transcript_stable_id_attributes,
                    )
                    or transcript_id
                )
                annotation.transcripts[transcript_id] = TranscriptRecord(
                    gene_id=gene_id,
                    seq_name=seq_name,
                    start=start,
                    end=end,
                    strand=strand_value,
                    biotype=resolve_biotype(
                        feature_type,
                        attributes,
                        transcript_id,
                        log,
                        source_config=source_config,
                    ),
                    stable_id=stable_id,
                    translation_coords=(
                        attributes.get(source_config.translation_coords_attribute)
                        if source_config.translation_coords_attribute
                        else None
                    ),
                )
                if source_config.transcript_rows_define_genes:
                    update_gene_from_transcript(
                        annotation,
                        gene_id,
                        seq_name,
                        start,
                        end,
                        strand_value,
                        annotation.transcripts[transcript_id].biotype,
                        log,
                    )
                continue

            if feature_type == "exon":
                transcript_id = parent_id(
                    required_attribute(
                        attributes,
                        source_config.exon_parent_attribute,
                        feature_type,
                        line_number,
                        gff_path,
                    ),
                    source_config,
                )
                if transcript_id not in annotation.transcripts:
                    biotype = (
                        annotation.genes[transcript_id].biotype
                        if transcript_id in annotation.genes
                        else resolve_biotype(
                            feature_type,
                            attributes,
                            transcript_id,
                            log,
                            source_config=source_config,
                        )
                    )
                    if transcript_id not in annotation.genes:
                        annotation.genes[transcript_id] = GeneRecord(
                            seq_name=seq_name,
                            start=start,
                            end=end,
                            strand=strand_value,
                            biotype=biotype,
                            stable_id=transcript_id,
                            name=transcript_id,
                        )
                    dummy_transcript_id = f"{transcript_id}_dTx"
                    annotation.transcripts.setdefault(
                        dummy_transcript_id,
                        TranscriptRecord(
                            gene_id=transcript_id,
                            seq_name=seq_name,
                            start=start,
                            end=end,
                            strand=strand_value,
                            biotype=biotype,
                            stable_id=dummy_transcript_id,
                        ),
                    )
                    transcript_id = dummy_transcript_id

                exon_stable_id = first_existing_attribute(
                    attributes,
                    source_config.exon_stable_id_attributes,
                )
                annotation.transcripts[transcript_id].exons.append(
                    ExonRecord(
                        start=start,
                        end=end,
                        strand=strand_value,
                        phase=None,
                        end_phase=None,
                        stable_id=(
                            normalize_id(exon_stable_id, source_config)
                            if exon_stable_id
                            else None
                        ),
                    )
                )
                continue

            if feature_type == "CDS":
                transcript_id = parent_id(
                    required_attribute(
                        attributes,
                        source_config.cds_parent_attribute,
                        feature_type,
                        line_number,
                        gff_path,
                    ),
                    source_config,
                )
                protein_stable_id = first_existing_attribute(
                    attributes,
                    source_config.translation_stable_id_attributes,
                )
                if (
                    protein_stable_id
                    and transcript_id in annotation.transcripts
                    and annotation.transcripts[transcript_id].protein_id is None
                ):
                    annotation.transcripts[transcript_id].protein_id = normalize_id(
                        protein_stable_id,
                        source_config,
                    )
                annotation.cds_segments.setdefault(transcript_id, []).append(
                    CdsSegment(
                        start=start,
                        end=end,
                        strand=strand_value,
                        phase=phase,
                    )
                )

    synthesize_cds_from_translation_coords(annotation, log)
    log.info(
        "Parsed %s genes, %s transcripts, and %s CDS transcript groups from %s",
        len(annotation.genes),
        len(annotation.transcripts),
        len(annotation.cds_segments),
        gff_path,
    )
    return annotation


def reconcile_annotation(
    annotation: ParsedAnnotation,
    logger: logging.Logger | None = None,
) -> ParsedAnnotation:
    """Patch missing gene-transcript-exon links required by the core schema."""

    log = logger or LOGGER
    added_exons = 0
    added_genes = 0
    added_transcripts = 0

    for transcript in annotation.transcripts.values():
        if not transcript.exons:
            transcript.exons.append(
                ExonRecord(
                    start=transcript.start,
                    end=transcript.end,
                    strand=transcript.strand,
                    phase=None,
                    end_phase=None,
                )
            )
            added_exons += 1

    for transcript_id, transcript in list(annotation.transcripts.items()):
        gene_id = transcript.gene_id
        if gene_id not in annotation.genes:
            annotation.genes[gene_id] = GeneRecord(
                seq_name=transcript.seq_name,
                start=transcript.start,
                end=transcript.end,
                strand=transcript.strand,
                biotype=transcript.biotype,
                stable_id=gene_id,
                name=gene_id,
            )
            added_genes += 1

    gene_ids_with_transcripts = {
        transcript.gene_id for transcript in annotation.transcripts.values()
    }
    for gene_id, gene in list(annotation.genes.items()):
        if gene_id in gene_ids_with_transcripts:
            continue

        dummy_transcript_id = f"{gene_id}_dTx"
        annotation.transcripts[dummy_transcript_id] = TranscriptRecord(
            gene_id=gene_id,
            seq_name=gene.seq_name,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            biotype=gene.biotype,
            stable_id=dummy_transcript_id,
            exons=[
                ExonRecord(
                    start=gene.start,
                    end=gene.end,
                    strand=gene.strand,
                    phase=None,
                    end_phase=None,
                )
            ],
        )
        added_transcripts += 1

    log.info(
        "Reconciled annotation (%s genes, %s transcripts, %s synthetic exons)",
        added_genes,
        added_transcripts,
        added_exons,
    )
    return annotation


def compute_exon_phases(annotation: ParsedAnnotation) -> ParsedAnnotation:
    """Compute exon phase and end_phase values from CDS segments."""

    for transcript_id, cds_list in annotation.cds_segments.items():
        transcript = annotation.transcripts.get(transcript_id)
        if transcript is None:
            continue

        strand = transcript.strand
        exons_in_transcript_order = sorted(
            transcript.exons,
            key=lambda exon: exon.start,
            reverse=(strand == -1),
        )
        if strand == 1:
            ordered_cds = sorted(cds_list, key=lambda cds: cds.start)
        else:
            ordered_cds = sorted(cds_list, key=lambda cds: cds.end, reverse=True)

        first_phase: int | None = None
        for cds in ordered_cds:
            if cds.phase is not None and cds.phase != ".":
                try:
                    first_phase = int(cds.phase)
                except ValueError:
                    pass
                break

        for exon in transcript.exons:
            exon.phase = -1
            exon.end_phase = -1

        coding_bases = 0
        first_coding_exon_seen = False
        for exon in exons_in_transcript_order:
            exon_coding_len = 0
            for cds in ordered_cds:
                overlap_start = max(exon.start, cds.start)
                overlap_end = min(exon.end, cds.end)
                if overlap_start <= overlap_end:
                    exon_coding_len += overlap_end - overlap_start + 1

            if exon_coding_len == 0:
                continue

            if not first_coding_exon_seen and first_phase is not None:
                exon.phase = first_phase
            else:
                exon.phase = coding_bases % 3

            coding_bases += exon_coding_len
            exon.end_phase = coding_bases % 3
            first_coding_exon_seen = True

    return annotation


def apply_biotype_overrides(
    annotation: ParsedAnnotation,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> ParsedAnnotation:
    """Apply final gene and transcript biotype overrides."""

    for gene in annotation.genes.values():
        if gene.biotype in source_config.gene_biotype_overrides:
            gene.biotype = source_config.gene_biotype_overrides[gene.biotype]

    for transcript in annotation.transcripts.values():
        gene_biotype = annotation.genes.get(transcript.gene_id)
        if (
            gene_biotype
            and gene_biotype.biotype in source_config.transcript_biotype_overrides
        ):
            transcript.biotype = source_config.transcript_biotype_overrides[
                gene_biotype.biotype
            ]

    return annotation


def prepare_annotation_for_load(
    converted_gff_path: str | Path,
    logger: logging.Logger | None = None,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> ParsedAnnotation:
    """Parse, reconcile, phase, and normalize one converted GFF3 annotation."""

    annotation = parse_converted_gff3(
        converted_gff_path,
        logger=logger,
        source_config=source_config,
    )
    reconcile_annotation(annotation, logger=logger)
    compute_exon_phases(annotation)
    apply_biotype_overrides(annotation, source_config=source_config)
    return annotation


def load_seq_regions_from_fna(
    fna_path: str | Path,
    cursor: DbCursor,
    coord_system_id: int,
    logger: logging.Logger | None = None,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
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
                (seq_region_id, source_config.toplevel_attrib_type_id),
            )

        seq_region_ids[name] = seq_region_id

    with Path(fna_path).open("r") as fasta_handle:
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


def core_schema_version(schema_sql_path: str | Path | None = None) -> str:
    """Return the bundled core schema version used in derived DB names."""

    schema_path = Path(schema_sql_path or DEFAULT_CORE_SCHEMA_SQL_PATH)
    try:
        schema_sql = schema_path.read_text()
    except OSError:
        return DEFAULT_ENSEMBL_RELEASE

    match = re.search(
        r"['\"]schema_version['\"]\s*,\s*['\"]([^'\"]+)['\"]",
        schema_sql,
    )
    if match:
        return match.group(1)
    return DEFAULT_ENSEMBL_RELEASE


def refseq_accession_db_token(assembly_accession: str) -> str:
    """Return a compact DB-safe token for a RefSeq assembly accession."""

    token = assembly_accession.lower().replace("_", "")
    token = token.replace(".", "v", 1)
    return re.sub(r"[^a-z0-9]+", "_", token).strip("_")


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

    raw_schema = Path(schema_sql_path).read_text()
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
        import pymysql  # type: ignore[import]

        return pymysql.connect(**connection_args)
    except ImportError as pymysql_error:
        try:
            import mysql.connector  # type: ignore[import]
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
    source_config: GffSourceConfig = REFSEQ_CONFIG,
) -> tuple[int, int]:
    """Insert coord_system, meta, and analysis bootstrap rows."""

    cursor.execute(
        "INSERT INTO coord_system (name,version,rank,attrib) "
        "VALUES ('primary_assembly','',1,'default_version,sequence_level')"
    )
    coord_system_id = int(cursor.lastrowid)
    cursor.executemany(
        "INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,%s,%s)",
        [
            ("species.scientific_name", species_name),
            ("assembly.accession", assembly_accession),
            ("assembly.name", assembly_accession),
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
                display_xref_id, source, canonical_transcript_id)
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,NULL,%s,%s)""",
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
            ),
        )


def insert_transcripts_and_exons(
    cursor: DbCursor,
    annotation: ParsedAnnotation,
    seq_region_ids: Mapping[str, int],
    gene_id_map: Mapping[str, int],
    transcript_id_map: Mapping[str, int],
    analysis_id: int,
    source_config: GffSourceConfig = REFSEQ_CONFIG,
    first_exon_id: int = 1,
) -> tuple[dict[tuple[str, int], int], dict[str, dict[tuple[int, int], int]]]:
    """Insert transcripts, exons, and exon_transcript links."""

    exon_id_map: dict[tuple[str, int], int] = {}
    per_transcript_coord_to_exon_id: dict[str, dict[tuple[int, int], int]] = {}
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

        sorted_exons = sorted(
            transcript.exons,
            key=lambda exon: exon.start,
            reverse=(transcript.strand == -1),
        )
        coord_to_exon_id: dict[tuple[int, int], int] = {}
        for rank, exon in enumerate(sorted_exons, start=1):
            exon_id = next_exon_id
            next_exon_id += 1
            phase_value = 0 if exon.phase in (None, -1) else exon.phase
            end_phase_value = 0 if exon.end_phase in (None, -1) else exon.end_phase
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
            exon_id_map[(transcript_id, rank)] = exon_id
            coord_to_exon_id[exon.coordinate_key] = exon_id
            cursor.execute(
                "INSERT INTO exon_transcript (exon_id, transcript_id, rank) "
                "VALUES (%s,%s,%s)",
                (exon_id, transcript_id_map[transcript_id], rank),
            )

        per_transcript_coord_to_exon_id[transcript_id] = coord_to_exon_id

    return exon_id_map, per_transcript_coord_to_exon_id


def insert_translations(
    cursor: DbCursor,
    annotation: ParsedAnnotation,
    transcript_id_map: Mapping[str, int],
    exon_id_map: Mapping[tuple[str, int], int],
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
            cursor.execute(
                "UPDATE transcript SET canonical_translation_id = %s WHERE transcript_id = %s",
                (translation_id, db_transcript_id),
            )


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

    del assembly_report_path
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

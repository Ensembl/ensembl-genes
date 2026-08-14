"""Parse and normalize converted GFF3/GTF annotations."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Protocol

try:  # Support both package imports and direct same-directory imports.
    from .gff_metadata import (
        open_text_maybe_gzip,
    )
    from .gff_models import (
        CdsSegment,
        ExonRecord,
        GeneRecord,
        ParsedAnnotation,
        TranscriptRecord,
    )
    from .gff_source_config import REFSEQ_CONFIG, GffSourceConfig
except ImportError:  # pragma: no cover - used when run beside this file.
    from gff_metadata import (  # type: ignore
        open_text_maybe_gzip,
    )
    from gff_models import (  # type: ignore
        CdsSegment,
        ExonRecord,
        GeneRecord,
        ParsedAnnotation,
        TranscriptRecord,
    )
    from gff_source_config import REFSEQ_CONFIG, GffSourceConfig  # type: ignore


LOGGER = logging.getLogger(__name__)


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

    return dict(item.split("=", 1) for item in raw_attributes.split(";") if "=" in item)


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

        if not exon.start <= cds_start <= cds_end <= exon.end:
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

    for _, transcript in list(annotation.transcripts.items()):
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
        coding_exons_seen = 0
        coding_exon_count = sum(
            1
            for exon in exons_in_transcript_order
            if any(
                max(exon.start, cds.start) <= min(exon.end, cds.end)
                for cds in ordered_cds
            )
        )
        for exon in exons_in_transcript_order:
            exon_coding_len = 0
            exon_overlap_start: int | None = None
            exon_overlap_end: int | None = None
            for cds in ordered_cds:
                overlap_start = max(exon.start, cds.start)
                overlap_end = min(exon.end, cds.end)
                if overlap_start <= overlap_end:
                    exon_coding_len += overlap_end - overlap_start + 1
                    exon_overlap_start = (
                        overlap_start
                        if exon_overlap_start is None
                        else min(exon_overlap_start, overlap_start)
                    )
                    exon_overlap_end = (
                        overlap_end
                        if exon_overlap_end is None
                        else max(exon_overlap_end, overlap_end)
                    )

            if exon_coding_len == 0:
                continue

            coding_exons_seen += 1
            five_prime_is_coding = (
                exon_overlap_start == exon.start
                if strand == 1
                else exon_overlap_end == exon.end
            )
            three_prime_is_coding = (
                exon_overlap_end == exon.end
                if strand == 1
                else exon_overlap_start == exon.start
            )

            if coding_exons_seen == 1:
                if five_prime_is_coding:
                    exon.phase = first_phase if first_phase is not None else 0
            else:
                exon.phase = coding_bases % 3

            coding_bases += exon_coding_len
            if coding_exons_seen < coding_exon_count or three_prime_is_coding:
                exon.end_phase = coding_bases % 3

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
        if transcript.biotype in source_config.transcript_biotype_overrides:
            transcript.biotype = source_config.transcript_biotype_overrides[
                transcript.biotype
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

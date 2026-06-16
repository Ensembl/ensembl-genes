"""Generic GFF3 feature models used by the core loader."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ExonRecord:
    """In-memory representation of an exon before database insertion."""

    start: int
    end: int
    strand: int
    phase: int | None = None
    end_phase: int | None = None
    stable_id: str | None = None

    @property
    def coordinate_key(self) -> tuple[int, int]:
        """Return the coordinate key used to map exons back to DB ids."""

        return (self.start, self.end)


@dataclass
class GeneRecord:
    """In-memory representation of a gene parsed from GFF3."""

    seq_name: str
    start: int
    end: int
    strand: int
    biotype: str
    stable_id: str
    name: str
    xref_geneid: str | None = None


@dataclass
class TranscriptRecord:
    """In-memory representation of a transcript and its exons."""

    gene_id: str
    seq_name: str
    start: int
    end: int
    strand: int
    biotype: str
    stable_id: str
    exons: list[ExonRecord] = field(default_factory=list)
    protein_id: str | None = None
    translation_coords: str | None = None


@dataclass(frozen=True)
class CdsSegment:
    """A CDS segment parsed from GFF3."""

    start: int
    end: int
    strand: int
    phase: str | None


@dataclass
class ParsedAnnotation:
    """All feature records needed to load one parsed GFF3 annotation."""

    genes: dict[str, GeneRecord] = field(default_factory=dict)
    transcripts: dict[str, TranscriptRecord] = field(default_factory=dict)
    cds_segments: dict[str, list[CdsSegment]] = field(default_factory=dict)

"""Shared data models for stable-ID decisions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

FEATURE_ORDER = ("gene", "transcript", "translation")
ACTION_ORDER = ("mapped", "missing", "new")

PK_BY_TYPE = {
    "gene": "gene_id",
    "transcript": "transcript_id",
    "translation": "translation_id",
}

TABLE_BY_TYPE = {
    "gene": "gene",
    "transcript": "transcript",
    "translation": "translation",
}


@dataclass(frozen=True)
class Feature:
    stable_id: str
    version: int
    seqid: str
    start: int
    end: int
    strand: str
    parent_stable_id: Optional[str] = None

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class Decision:
    feature_type: str
    action: str
    current_stable_id: Optional[str]
    current_version: int
    old_stable_id: Optional[str]
    old_version: int
    new_stable_id: Optional[str]
    new_version: int
    mapping_session_id: int
    score: float
    reason: str


@dataclass(frozen=True)
class RawFeature:
    stable_id: str
    version: int
    feature_type: str
    seqid: str
    start: int
    end: int
    strand: str
    parent_ids: tuple[str, ...]


@dataclass(frozen=True)
class CdsRow:
    parent_stable_id: str
    parent_version: Optional[int]
    translation_stable_id: str
    translation_version: Optional[int]
    seqid: str
    start: int
    end: int
    strand: str


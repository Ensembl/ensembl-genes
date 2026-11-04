"""Schema and metadata normalization for expression-track handovers."""

from __future__ import annotations

import re
from typing import Any, Literal
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


def value(raw: Any) -> str:
    """Return a stripped string representation for identifiers and labels."""
    return "" if raw is None else str(raw).strip()


def clean_key(raw: str) -> str:
    """Normalize an input metadata key for alias matching."""
    return re.sub(r"\s+", " ", raw.strip().lower().replace("-", " "))


ATTRIBUTE_ALIASES: dict[str, set[str]] = {
    "organism": {"organism", "scientific_name"},
    "tissue": {"tissue", "organism part", "organism_part", "organ part", "anatomy"},
    "cell_type": {"cell type", "cell_type"},
    "cell_line": {"cell line", "cell_line"},
    "sex": {"sex"},
    "breed": {"breed", "breed name", "breed_name"},
    "strain": {"strain", "strain name", "strain_name"},
    "age": {"age"},
    "developmental_stage": {
        "developmental stage",
        "developmental_stage",
        "dev stage",
        "dev_stage",
    },
    "genotype": {"genotype"},
    "treatment": {"treatment"},
    "condition": {"condition"},
    "disease": {"disease", "disease status", "disease_status"},
    "timepoint": {"timepoint", "time point", "time_point"},
    "replicate": {"replicate", "biological replicate", "biological_replicate"},
    "pooling": {"pooling", "pool", "pool identifier", "pool_identifier"},
    "cultivar": {"cultivar"},
    "isolate": {"isolate"},
}


DISPLAY_LABEL_RULES: dict[str, dict[str, str]] = {
    # Keep this table deliberately explicit. Generic title-casing would damage
    # cell-line names and established abbreviations.
    "tissue": {
        "bone marrow": "Bone Marrow",
        "brain": "Brain",
        "brain cns": "Brain CNS",
        "heart": "Heart",
    },
    "genotype": {
        "wild type": "Wild Type",
        "wildtype": "Wild Type",
    },
    "sex": {
        "female": "Female",
        "male": "Male",
        "unknown": "Unknown",
    },
}


def canonical_attribute_name(raw: Any) -> str | None:
    """Return the project field name for a BioSample attribute alias."""
    cleaned = clean_key(value(raw))
    for canonical, aliases in ATTRIBUTE_ALIASES.items():
        if cleaned in {clean_key(alias) for alias in aliases}:
            return canonical
    return None


def canonical_metadata_field_name(raw: Any) -> str:
    """Return a stable snake-case name for a metadata field."""
    attribute_name = canonical_attribute_name(raw)
    if attribute_name:
        return attribute_name
    return re.sub(r"[^a-z0-9]+", "_", value(raw).casefold()).strip("_")


def standardize_metadata_value(field: str, raw: Any) -> Any:
    """Apply only explicit display rules, preserving unknown labels."""
    if isinstance(raw, list):
        return [standardize_metadata_value(field, item) for item in raw]
    if not isinstance(raw, str):
        return raw
    cleaned = re.sub(r"\s+", " ", raw.strip())
    lookup = re.sub(r"[_-]+", " ", cleaned).casefold()
    return DISPLAY_LABEL_RULES.get(field, {}).get(lookup, cleaned)


def presentable_label(field: str, raw: Any) -> Any:
    """Create a UI label from a controlled aggregate label."""
    standardized = standardize_metadata_value(field, raw)
    if not isinstance(standardized, str):
        return standardized
    raw_clean = re.sub(r"\s+", " ", raw.strip()) if isinstance(raw, str) else raw
    if standardized != raw_clean:
        return standardized
    words = re.sub(r"[_-]+", " ", standardized).split()
    if len(words) < 2:
        return standardized[:1].upper() + standardized[1:] if standardized else standardized
    return " ".join(word if word.isupper() else word.capitalize() for word in words)


class MetadataProvenance(BaseModel):
    """Evidence for one metadata value."""

    model_config = ConfigDict(extra="allow")

    source: str
    source_field: str
    source_id: str | None = None
    source_value: Any = None
    normalization: str | None = None


class MetadataField(BaseModel):
    """A metadata value with extensible annotations."""

    model_config = ConfigDict(extra="allow")

    value: Any
    provenance: MetadataProvenance | None = None


class TargetGenome(BaseModel):
    """Genome information embedded in the project payload."""

    model_config = ConfigDict(extra="allow")

    species: str
    assembly: str
    assembly_uuid: UUID


class PayloadProvenance(BaseModel):
    """Relationships from the payload back to source entities."""

    model_config = ConfigDict(extra="allow")

    biosample_accession: str | None = None
    study_accession: str | None = None
    run_accessions: list[str] = Field(default_factory=list)
    aggregate_id: str | None = None
    aggregate_label: str | None = None


class TrackMetadataPayload(BaseModel):
    """Draft schema for the JSON stored in ``track_metadata_payload``."""

    model_config = ConfigDict(extra="allow")

    schema_version: str
    selection_id: str
    selection_level: Literal["run", "tissue"]
    data_type: Literal["rnaseq"]
    display_label: str
    target_genome: TargetGenome
    metadata: dict[str, MetadataField] = Field(default_factory=dict)
    provenance: PayloadProvenance
    selection_hierarchy: list[str] = Field(min_length=1)

    @field_validator("metadata")
    @classmethod
    def validate_metadata(cls, fields: dict[str, MetadataField]) -> dict[str, MetadataField]:
        for field_name, field in fields.items():
            canonical_name = canonical_metadata_field_name(field_name)
            if field_name != canonical_name:
                raise ValueError(
                    f"metadata field {field_name!r} is not canonical; "
                    f"use {canonical_name!r}"
                )
            source_value = field.value
            if field.provenance and field.provenance.source_value is not None:
                source_value = field.provenance.source_value
            standardized = standardize_metadata_value(field_name, source_value)
            if field_name == "tissue":
                standardized = presentable_label(field_name, source_value)
            if field.value != standardized:
                if field.provenance is None or field.provenance.source_value is None:
                    raise ValueError(
                        f"normalized metadata field {field_name!r} requires "
                        "provenance.source_value"
                    )
                field.value = standardized
            if field.provenance and standardized != source_value:
                field.provenance.normalization = "explicit_display_label_rule"
        return fields

    @model_validator(mode="after")
    def validate_display_label(self) -> "TrackMetadataPayload":
        if self.selection_level == "tissue":
            self.display_label = str(presentable_label("tissue", self.display_label))
        return self

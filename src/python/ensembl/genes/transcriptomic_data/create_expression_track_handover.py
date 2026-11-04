#!/usr/bin/env python3
"""Create an assembly-aware expression-track handover.

This command is the boundary between transcriptomic processing manifests and
the expression-track handover API.  It consumes explicit selected-run,
run-metadata, BioSample, product, aggregate, and optional ENA source-BAM
manifests.  It never discovers products by walking directories or derives
aggregate membership from filenames.

Assembly resolution uses the same metadata database URI convention as the
metadata tooling in this repository: ``--metadata-uri`` or ``METADATA_URI``.
The resolver is kept behind :class:`AssemblyResolver` so it can be replaced by
an API-backed adapter without changing handover rendering.
"""

# This module is intentionally a single, explicit CLI boundary for assembling
# handover records from several manifests.  Keep the data-flow checks together;
# the standalone schema validator owns structural validation.
# pylint: disable=too-many-lines,missing-function-docstring,missing-class-docstring
# pylint: disable=too-few-public-methods,broad-exception-caught,redefined-outer-name
# pylint: disable=too-many-arguments,too-many-locals,too-many-return-statements
# pylint: disable=too-many-branches,too-many-statements,consider-using-dict-items

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, NoReturn, Protocol
from urllib.parse import urlparse

from sqlalchemy import create_engine, text as sql_text

from ensembl.genes.transcriptomic_data.expression_track_handover_schema import (
    ATTRIBUTE_ALIASES,
    canonical_attribute_name,
    canonical_metadata_field_name,
    presentable_label,
    standardize_metadata_value,
)


SCHEMA_VERSION = "transcriptomic-selection-0.3.1"
DATASET_TYPE = "expression_tracks"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
MD5_RE = re.compile(r"^[0-9a-fA-F]{32}$")

PRODUCTS: dict[str, dict[str, str]] = {
    "coverage_forward": {
        "track_type": "rnaseq_coverage",
        "format": "BigWig",
        "kind": "bigwig",
        "role": "summary",
        "strand": "forward",
    },
    "coverage_reverse": {
        "track_type": "rnaseq_coverage",
        "format": "BigWig",
        "kind": "bigwig",
        "role": "summary",
        "strand": "reverse",
    },
    "gene_model": {
        "track_type": "gene_models",
        "format": "BigBed",
        "kind": "bigbed",
        "role": "details",
        "slug": "gene_models",
    },
    "gene_models": {
        "track_type": "gene_models",
        "format": "BigBed",
        "kind": "bigbed",
        "role": "details",
    },
    "splice_junction": {
        "track_type": "splice_junctions",
        "format": "BigBed",
        "kind": "bigbed",
        "role": "details",
        "slug": "splice_junctions",
    },
    "splice_junctions": {
        "track_type": "splice_junctions",
        "format": "BigBed",
        "kind": "bigbed",
        "role": "details",
    },
}
DEFAULT_RUN_PRODUCTS = (
    "coverage_forward",
    "coverage_reverse",
    "gene_models",
    "splice_junctions",
)
DEFAULT_AGGREGATE_PRODUCTS = (
    "coverage_forward",
    "coverage_reverse",
    "splice_junctions",
)


class HandoverError(ValueError):
    """A required handover invariant was violated."""


class AssemblyResolver(Protocol):
    """Resolve a GCA accession to one active assembly UUID."""

    def resolve(self, gca_accession: str) -> "AssemblyResolution": ...


@dataclass(frozen=True)
class AssemblyResolution:
    gca_accession: str
    assembly_uuid: str
    linked_genomes: tuple[dict[str, Any], ...]
    evidence: dict[str, Any]


@dataclass
class Diagnostic:
    severity: str
    code: str
    message: str
    context: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "severity": self.severity,
            "code": self.code,
            "message": self.message,
            "context": self.context,
        }


def fail(message: str) -> NoReturn:
    raise HandoverError(message)


class MetadataDbAssemblyResolver:
    """Resolve assembly UUIDs using the repository's metadata DB URI style."""

    def __init__(self, metadata_uri: str | None = None, engine: Any | None = None):
        if engine is not None:
            self.engine = engine
            self.metadata_uri = metadata_uri or "injected-engine"
            return
        uri = metadata_uri or os.getenv("METADATA_URI")
        if not uri:
            fail("metadata DB URI is required via --metadata-uri or METADATA_URI")
        try:
            self.engine = create_engine(uri)
        except (
            Exception
        ) as exc:  # pragma: no cover - driver errors are environment-specific
            fail(f"cannot create metadata DB engine: {exc}")
        self.metadata_uri = uri

    def resolve(self, gca_accession: str) -> AssemblyResolution:
        uuid_query = sql_text(
            """
            SELECT DISTINCT genome.genome_uuid
            FROM genome
            JOIN assembly ON genome.assembly_id = assembly.assembly_id
            JOIN genome_dataset ON genome.genome_id = genome_dataset.genome_id
            JOIN dataset ON genome_dataset.dataset_id = dataset.dataset_id
            WHERE dataset.name = :dataset_name
              AND assembly.accession = :gca_accession
            """
        )
        linked_query = sql_text(
            """
            SELECT DISTINCT linked_genome.genome_uuid, linked_assembly.accession
            FROM genome AS anchor_genome
            JOIN assembly AS anchor_assembly
              ON anchor_genome.assembly_id = anchor_assembly.assembly_id
            JOIN genome AS linked_genome
              ON linked_genome.assembly_id = anchor_genome.assembly_id
            JOIN assembly AS linked_assembly
              ON linked_genome.assembly_id = linked_assembly.assembly_id
            WHERE anchor_genome.genome_uuid = :assembly_uuid
            """
        )
        try:
            with self.engine.connect() as connection:
                rows = connection.execute(
                    uuid_query,
                    {"dataset_name": "genebuild", "gca_accession": gca_accession},
                ).fetchall()
                uuids = sorted({str(row[0]) for row in rows if row[0]})
                if len(uuids) != 1:
                    fail(
                        f"metadata DB returned {len(uuids)} active assembly UUIDs for "
                        f"{gca_accession}; expected exactly one"
                    )
                linked = connection.execute(
                    linked_query, {"assembly_uuid": uuids[0]}
                ).fetchall()
        except HandoverError:
            raise
        except Exception as exc:
            fail(f"metadata DB assembly lookup failed for {gca_accession}: {exc}")

        linked_genomes = tuple(
            {"genome_uuid": str(row[0]), "gca_accession": str(row[1])}
            for row in linked
            if row[0]
        )
        return AssemblyResolution(
            gca_accession=gca_accession,
            assembly_uuid=uuids[0],
            linked_genomes=linked_genomes,
            evidence={
                "metadata_uri": self.metadata_uri,
                "dataset": "genebuild",
                "gca_accession": gca_accession,
            },
        )


def json_load(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read JSON input {path}: {exc}")


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                fail(f"TSV input has no header: {path}")
            return [dict(row) for row in reader]
    except OSError as exc:
        fail(f"cannot read TSV input {path}: {exc}")


def value(value: Any) -> str:
    return "" if value is None else str(value).strip()


def first(row: Mapping[str, Any], *keys: str) -> str:
    for key in keys:
        result = value(row.get(key))
        if result:
            return result
    return ""


def jsonish(raw: Any, default: Any = None) -> Any:
    if raw is None or raw == "":
        return default
    if not isinstance(raw, str):
        return raw
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        return default


def split_values(raw: Any) -> list[str]:
    parsed = jsonish(raw, raw)
    if isinstance(parsed, list):
        return [value(item) for item in parsed if value(item)]
    if isinstance(parsed, dict):
        parsed = parsed.get("run_accessions", [])
        return [value(item) for item in parsed if value(item)]
    text_value = value(parsed)
    return [item.strip() for item in re.split(r"[;,|]", text_value) if item.strip()]


def clean_value(raw: Any) -> Any:
    if isinstance(raw, str):
        return re.sub(r"\s+", " ", raw.strip())
    return raw


def scalar_or_list(items: Iterable[Any]) -> Any:
    distinct: list[Any] = []
    for item in items:
        item = clean_value(item)
        if item not in (None, "") and item not in distinct:
            distinct.append(item)
    if not distinct:
        return None
    return distinct[0] if len(distinct) == 1 else distinct


def metadata_field(
    field: str,
    value: Any,
    source: str,
    source_field: str,
    source_id: str | None = None,
    source_value: Any = None,
) -> dict[str, Any]:
    """Build one extensible metadata field with auditable source information."""
    standardized = standardize_metadata_value(field, value)
    provenance: dict[str, Any] = {
        "source": source,
        "source_field": source_field,
        "source_value": value if source_value is None else source_value,
    }
    if source_id:
        provenance["source_id"] = source_id
    if standardized != (value if source_value is None else source_value):
        provenance["normalization"] = "explicit_display_label_rule"
    return {"value": standardized, "provenance": provenance}


def source_attributes(record: Mapping[str, Any]) -> list[dict[str, Any]]:
    source = record.get("attributes", record.get("metadata", {}))
    if isinstance(source, list):
        return [item for item in source if isinstance(item, dict) and "name" in item]
    if isinstance(source, dict):
        result = []
        for key, raw in source.items():
            result.append(
                {
                    "name": key,
                    "value": (
                        raw.get("value")
                        if isinstance(raw, dict) and "value" in raw
                        else raw
                    ),
                }
            )
        return result
    fail("BioSample attributes/metadata must be a list or object")


def biological_attributes(record: Mapping[str, Any]) -> dict[str, Any]:
    by_name: dict[str, list[Any]] = defaultdict(list)
    for item in source_attributes(record):
        canonical = canonical_attribute_name(item["name"])
        if canonical:
            by_name[canonical].append(item.get("value"))
    result: dict[str, Any] = {}
    for canonical in ATTRIBUTE_ALIASES:
        candidate_values = by_name.get(canonical, [])
        resolved = scalar_or_list(candidate_values)
        if resolved is not None:
            result[canonical] = resolved
    return result


def json_records(path: Path, keys: tuple[str, ...]) -> list[dict[str, Any]]:
    raw = json_load(path)
    if isinstance(raw, list):
        records = raw
    elif isinstance(raw, dict):
        records = None
        for key in keys:
            candidate = raw.get(key)
            if isinstance(candidate, list):
                records = candidate
                break
        if records is None:
            records = [
                dict(record, _key=key)
                for key, record in raw.items()
                if isinstance(record, dict)
            ]
    else:
        records = None
    if not isinstance(records, list) or not all(
        isinstance(item, dict) for item in records
    ):
        fail(f"{path} must contain a list or an object containing one")
    return records


def load_selected(path: Path) -> dict[str, dict[str, str]]:
    selected: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(read_tsv(path), start=2):
        run_id = first(row, "run_accession", "run_id", "entity_id")
        sample_id = first(row, "sample_id", "sample_accession", "biosample_id")
        gca = first(row, "gca_accession", "assembly_accession", "assembly")
        if not run_id or not sample_id or not gca:
            fail(
                f"selected-run row {line_number} requires run, sample, and GCA accessions"
            )
        if run_id in selected:
            fail(f"duplicate selected run: {run_id}")
        selected[run_id] = {
            "run_accession": run_id,
            "sample_accession": sample_id,
            "gca_accession": gca,
            **row,
        }
    if not selected:
        fail(f"selected-run manifest is empty: {path}")
    return selected


def load_run_metadata(path: Path) -> dict[str, dict[str, Any]]:
    runs: dict[str, dict[str, Any]] = {}
    for record in json_records(path, ("runs", "records", "data")):
        run_id = first(record, "run_accession", "run_id", "accession", "_key")
        if not run_id:
            fail(f"run metadata record has no run accession: {path}")
        sample_id = first(record, "sample_accession", "sample_id", "biosample_id")
        study_id = first(record, "study_accession", "study_id", "project_accession")
        previous = runs.get(run_id)
        if previous is not None and (
            sample_id
            and sample_id != previous["sample_accession"]
            or study_id
            and study_id != previous["study_accession"]
        ):
            fail(f"conflicting ENA metadata for run {run_id}")
        runs.setdefault(
            run_id,
            {
                "run_accession": run_id,
                "sample_accession": sample_id,
                "study_accession": study_id,
                "raw": record,
            },
        )
    return runs


def load_biosamples(path: Path) -> dict[str, dict[str, Any]]:
    samples: dict[str, dict[str, Any]] = {}
    for record in json_records(path, ("samples", "biosamples", "records", "data")):
        sample_id = first(
            record, "sample_accession", "biosample_id", "accession", "id", "_key"
        )
        if not sample_id:
            fail(f"BioSample record has no accession: {path}")
        if sample_id in samples:
            fail(f"duplicate BioSample accession: {sample_id}")
        raw_attributes = source_attributes(record)
        attribute_sources: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for item in raw_attributes:
            canonical = canonical_attribute_name(item["name"])
            if canonical:
                attribute_sources[canonical].append(
                    {
                        "field": value(item["name"]),
                        "value": item.get("value"),
                    }
                )
        samples[sample_id] = {
            "sample_accession": sample_id,
            "attributes": biological_attributes(record),
            "attribute_sources": dict(attribute_sources),
        }
    return samples


def load_source_bams(path: Path | None) -> dict[str, dict[str, Any]]:
    if path is None:
        return {}
    records = (
        read_tsv(path)
        if path.suffix.lower() in {".tsv", ".tab"}
        else json_records(path, ("source_bams", "submissions", "records", "data"))
    )
    submissions: dict[str, dict[str, Any]] = {}
    for record in records:
        run_id = first(record, "run_accession", "run_id", "entity_id")
        if not run_id:
            fail(f"source-BAM record has no run accession: {path}")
        if run_id in submissions:
            fail(f"duplicate source-BAM record for run {run_id}")
        submissions[run_id] = dict(record)
    return submissions


def product_config(raw_type: str) -> tuple[str, dict[str, str]]:
    key = value(raw_type).lower()
    if key not in PRODUCTS:
        fail(f"unsupported track product: {raw_type}")
    config = PRODUCTS[key]
    return config.get("slug", key), config


def load_product_rows(path: Path, level: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for line_number, raw in enumerate(read_tsv(path), start=2):
        entity_id = first(raw, "entity_id", "run_accession", "aggregate_id")
        raw_type = first(
            raw, "product", "product_slug", "track_type", "track_product", "type"
        )
        track_file = first(
            raw, "track_file", "track_path", "path", "artifact", "source_path"
        )
        if not entity_id or not raw_type or not track_file:
            fail(
                f"{path} row {line_number} requires entity_id, product type, and track_file"
            )
        product_slug, config = product_config(raw_type)
        key = (entity_id, product_slug)
        if key in seen:
            fail(f"duplicate product row in {path}: {key}")
        seen.add(key)
        source_runs = split_values(
            first(
                raw,
                "source_run_accessions",
                "source_runs",
                "run_accessions",
                "constituent_run_ids",
                "source_entity_ids",
            )
        )
        if level == "run" and not source_runs:
            source_runs = [first(raw, "run_accession", "run_id", "entity_id")]
        metadata = jsonish(first(raw, "metadata_json", "metadata"), {})
        provenance = jsonish(first(raw, "provenance_json", "provenance"), {})
        if not isinstance(metadata, dict) or not isinstance(provenance, dict):
            fail(f"{path} row {line_number} has non-object metadata/provenance JSON")
        source_count_raw = first(raw, "source_count", "source_run_count", "run_count")
        source_count = (
            int(source_count_raw)
            if source_count_raw.isdigit()
            else (len(source_runs) if level == "run" else None)
        )
        rows.append(
            {
                "selection_level": level,
                "entity_id": entity_id,
                "run_accession": first(raw, "run_accession", "run_id")
                or (entity_id if level == "run" else ""),
                "aggregate_id": (
                    first(raw, "aggregate_id", "entity_id")
                    if level == "aggregate"
                    else ""
                ),
                "aggregate_label": (
                    first(
                        raw,
                        "aggregate_label",
                        "tissue",
                        "label",
                        "group_label",
                        "entity_id",
                    )
                    if level == "aggregate"
                    else ""
                ),
                "gca_accession": first(
                    raw, "gca_accession", "assembly_accession", "assembly"
                ),
                "track_type": config["track_type"],
                "product_slug": product_slug,
                "format": config["format"],
                "kind": config["kind"],
                "role": config["role"],
                "alignment_strand": config.get("strand", ""),
                "track_file": track_file,
                "status": first(raw, "status", "product_status").lower() or "complete",
                "sha256": first(raw, "sha256", "sha_256"),
                "md5": first(raw, "md5"),
                "source_runs": source_runs,
                "source_count": source_count,
                "metadata": metadata,
                "provenance": provenance,
                "raw": raw,
            }
        )
    return rows


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        fail(f"cannot hash artifact {path}: {exc}")
    return digest.hexdigest()


def input_hash(path: Path | None) -> str | None:
    if path is None or not path.is_file():
        return None
    return sha256_file(path)


def add_diagnostic(
    diagnostics: list[Diagnostic],
    code: str,
    message: str,
    severity: str = "error",
    **context: Any,
) -> None:
    diagnostics.append(Diagnostic(severity, code, message, context))


def validate_artifact(
    row: dict[str, Any], diagnostics: list[Diagnostic], mode: str
) -> bool:
    if row["status"] != "complete":
        add_diagnostic(
            diagnostics,
            "incomplete_product",
            f"product status is {row['status']}",
            entity_id=row["entity_id"],
            product=row["product_slug"],
        )
        return False
    parsed = urlparse(row["track_file"])
    if not parsed.scheme:
        path = Path(row["track_file"])
        if not path.is_file():
            add_diagnostic(
                diagnostics,
                "missing_artifact",
                f"artifact is not locally accessible: {path}",
                severity="warning" if mode == "draft" and row["sha256"] else "error",
                entity_id=row["entity_id"],
                product=row["product_slug"],
            )
            return mode == "draft" and bool(row["sha256"])
        actual = sha256_file(path)
        if row["sha256"] and row["sha256"].lower() != actual:
            add_diagnostic(
                diagnostics,
                "checksum_mismatch",
                f"sha256 does not match: {path}",
                entity_id=row["entity_id"],
                product=row["product_slug"],
            )
            return False
        row["sha256"] = actual
    elif not row["sha256"]:
        add_diagnostic(
            diagnostics,
            "missing_checksum",
            f"URI artifact has no sha256: {row['track_file']}",
            entity_id=row["entity_id"],
            product=row["product_slug"],
        )
        return False
    if not SHA256_RE.fullmatch(row["sha256"]):
        add_diagnostic(
            diagnostics,
            "invalid_checksum",
            f"invalid sha256: {row['track_file']}",
            entity_id=row["entity_id"],
            product=row["product_slug"],
        )
        return False
    if row["md5"] and not MD5_RE.fullmatch(row["md5"]):
        add_diagnostic(
            diagnostics,
            "invalid_md5",
            f"invalid md5: {row['track_file']}",
            entity_id=row["entity_id"],
            product=row["product_slug"],
        )
        return False
    return True


def slug(raw: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", raw.casefold()).strip("_") or "unknown"


def accession_slug(raw: str) -> str:
    """Keep accession case while removing unsafe name characters."""
    return re.sub(r"[^A-Za-z0-9]+", "_", raw).strip("_") or "unknown"


def build_track_name(species: str, row: Mapping[str, Any]) -> str:
    scope = "run" if row["selection_level"] == "run" else "tissue"
    entity = (
        accession_slug(row["run_accession"])
        if scope == "run"
        else slug(row["aggregate_label"])
    )
    return f"{slug(species)}_pilot_{scope}_{entity}_{row['product_slug']}"


def build_track_description(species: str, gca: str, row: Mapping[str, Any]) -> str:
    if row["selection_level"] == "run":
        detail = f"run {row['run_accession']}"
        if row["alignment_strand"]:
            detail += f" ({row['alignment_strand']} strand)"
        return (
            f"{species} RNA-seq {row['track_type'].replace('_', ' ')} track for "
            f"{detail} on {gca}. Genebuild data-check pilot."
        )
    count = row.get("source_count") or len(row.get("source_runs", []))
    detail = f"tissue group {row['aggregate_label']}"
    if row["alignment_strand"]:
        detail += f" ({row['alignment_strand']} strand)"
    return (
        f"{species} RNA-seq {row['track_type'].replace('_', ' ')} aggregate for "
        f"{detail}, combining {count} selected runs on {gca}. "
        "Genebuild data-check pilot."
    )


def render_record(
    row: Mapping[str, Any],
    species: str,
    assembly: AssemblyResolution,
    label: str,
    team: str,
    user: str,
    hierarchy: list[str],
    run: Mapping[str, Any] | None,
    sample: Mapping[str, Any] | None,
    source_bams: Mapping[str, Mapping[str, Any]],
    dummy_track_source: bool,
) -> dict[str, Any]:
    name = build_track_name(species, row)
    track_file = row["track_file"]

    if row["selection_level"] == "run":
        if run is None or sample is None:
            fail(f"cannot render run record without metadata: {row['run_accession']}")
        assert sample is not None
        metadata: dict[str, dict[str, Any]] = {}
        for field, field_value in sample["attributes"].items():
            sources = sample.get("attribute_sources", {}).get(field, [])
            source = sources[0] if sources else {}
            metadata[field] = metadata_field(
                field,
                field_value,
                "ENA BioSample",
                source.get("field", field),
                sample["sample_accession"],
                source.get("value"),
            )
        if run.get("study_accession"):
            metadata["study"] = metadata_field(
                "study",
                run["study_accession"],
                "ENA run metadata",
                "study_accession",
                row["run_accession"],
            )
        provenance: dict[str, Any] = {
            "biosample_accession": sample["sample_accession"],
            "run_accessions": [row["run_accession"]],
        }
        if run.get("study_accession"):
            provenance["study_accession"] = run["study_accession"]
        provenance.update(row.get("provenance", {}))
        selection_level = "run"
        display_label = row["run_accession"]
    else:
        metadata = {
            "tissue": metadata_field(
                "tissue",
                presentable_label("tissue", row["aggregate_label"]),
                "aggregate_manifest",
                "aggregate_label",
                row["aggregate_id"],
                row["aggregate_label"],
            )
        }
        aggregate_metadata = dict(row.get("metadata", {}))
        # Keep the structured derived tissue field authoritative; raw
        # aggregate metadata may contain a plain ``tissue`` string.
        aggregate_metadata.pop("tissue", None)
        for raw_field, field_value in aggregate_metadata.items():
            field = canonical_metadata_field_name(raw_field)
            metadata[field] = metadata_field(
                field,
                field_value,
                "aggregate_manifest",
                value(raw_field),
                row["aggregate_id"],
            )
        provenance = {
            "aggregate_id": row["aggregate_id"],
            "aggregate_label": row["aggregate_label"],
            "run_accessions": list(row["source_runs"]),
            **row.get("provenance", {}),
        }
        selection_level = "tissue"
        display_label = presentable_label("tissue", row["aggregate_label"])

    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "selection_id": name,
        "selection_level": selection_level,
        "data_type": "rnaseq",
        "display_label": display_label,
        "target_genome": {
            "species": species,
            "assembly": assembly.gca_accession,
            "assembly_uuid": assembly.assembly_uuid,
        },
        "metadata": metadata,
        "provenance": provenance,
        "selection_hierarchy": hierarchy,
    }
    record: dict[str, Any] = {
        "team_responsible": team,
        "user": user,
        "name": name,
        "description": build_track_description(species, assembly.gca_accession, row),
        "label": label,
        "dataset_type": DATASET_TYPE,
        "dataset_source": {"name": track_file, "type": row["kind"]},
        "assembly_uuid": assembly.assembly_uuid,
        "track_files": {row["role"]: track_file},
        "track_metadata_payload": payload,
    }
    if row["selection_level"] == "run":
        submission = source_bams.get(row["run_accession"])
        source_url = first(submission or {}, "url", "ena_url", "file_url")
        if source_url:
            record["track_source"] = {
                "name": "ENA",
                "details": f"Run {row['run_accession']}",
                "url": source_url,
            }
    if dummy_track_source and "track_source" not in record:
        if row["selection_level"] == "run":
            run_accession = row["run_accession"]
            record["track_source"] = {
                "name": "ENA",
                "details": f"Run {run_accession}",
                "url": f"https://www.ebi.ac.uk/ena/browser/view/{run_accession}",
            }
        else:
            aggregate_id = slug(row["aggregate_id"])
            record["track_source"] = {
                "name": "ENA project placeholder",
                "details": (
                    "DRAFT PLACEHOLDER - Source BAM provenance for aggregate "
                    f"{row['aggregate_id']}"
                ),
                "url": f"https://example.org/ena-project/{aggregate_id}",
            }
    return record


def build_handover(
    args: argparse.Namespace,
    resolver: AssemblyResolver | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    selected = load_selected(args.selected_runs)
    runs = load_run_metadata(args.run_metadata)
    samples = load_biosamples(args.biosamples)
    resolver = resolver or MetadataDbAssemblyResolver(args.metadata_uri)
    gcas = {row["gca_accession"] for row in selected.values()}
    if len(gcas) != 1:
        fail(f"selected runs contain multiple GCA accessions: {sorted(gcas)}")
    assembly = resolver.resolve(next(iter(gcas)))
    source_bams = load_source_bams(args.source_bams)
    diagnostics: list[Diagnostic] = []

    for run_id, selected_row in selected.items():
        run = runs.get(run_id)
        if run is None:
            add_diagnostic(
                diagnostics,
                "missing_run_metadata",
                f"no ENA metadata for {run_id}",
                run_accession=run_id,
            )
            continue
        if not run["sample_accession"]:
            add_diagnostic(
                diagnostics,
                "missing_biosample",
                f"ENA metadata has no BioSample for {run_id}",
                run_accession=run_id,
            )
        elif run["sample_accession"] != selected_row["sample_accession"]:
            add_diagnostic(
                diagnostics,
                "sample_conflict",
                f"selected and ENA BioSample differ for {run_id}",
                run_accession=run_id,
                selected=selected_row["sample_accession"],
                ena=run["sample_accession"],
            )
        if selected_row["sample_accession"] not in samples:
            add_diagnostic(
                diagnostics,
                "missing_biosample_metadata",
                f"no BioSample metadata for {selected_row['sample_accession']}",
                run_accession=run_id,
            )

    run_rows = load_product_rows(args.run_manifest, "run")
    aggregate_rows = (
        load_product_rows(args.aggregate_manifest, "aggregate")
        if args.aggregate_manifest
        else []
    )
    all_rows = run_rows + aggregate_rows
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in all_rows:
        grouped[(row["selection_level"], row["entity_id"])].append(row)
        row["_emit"] = validate_artifact(row, diagnostics, args.mode)
        if row["gca_accession"] and row["gca_accession"] != assembly.gca_accession:
            add_diagnostic(
                diagnostics,
                "wrong_assembly",
                f"product uses {row['gca_accession']}, expected {assembly.gca_accession}",
                entity_id=row["entity_id"],
                product=row["product_slug"],
            )
            row["_emit"] = False
        if row["selection_level"] == "run" and row["run_accession"] not in selected:
            add_diagnostic(
                diagnostics,
                "unselected_run",
                f"product run is not selected: {row['run_accession']}",
                run_accession=row["run_accession"],
            )
            row["_emit"] = False
        if row["selection_level"] == "aggregate":
            row["source_runs"] = list(dict.fromkeys(row["source_runs"]))
            if row["source_count"] is not None and row["source_count"] != len(
                row["source_runs"]
            ):
                add_diagnostic(
                    diagnostics,
                    "source_count_mismatch",
                    "aggregate source count differs from source-run list",
                    entity_id=row["entity_id"],
                )
                row["_emit"] = False
            if not row["source_runs"]:
                add_diagnostic(
                    diagnostics,
                    "missing_aggregate_sources",
                    "aggregate has no source-run list",
                    entity_id=row["entity_id"],
                )
                row["_emit"] = False
            if any(run_id not in selected for run_id in row["source_runs"]):
                add_diagnostic(
                    diagnostics,
                    "unselected_aggregate_source",
                    "aggregate contains a run outside the selected set",
                    entity_id=row["entity_id"],
                    source_runs=row["source_runs"],
                )
                row["_emit"] = False

    for (level, entity_id), rows in grouped.items():
        expected = set(args.run_products if level == "run" else args.aggregate_products)
        actual = {row["product_slug"] for row in rows}
        missing = sorted(expected - actual)
        if missing:
            add_diagnostic(
                diagnostics,
                "incomplete_product_matrix",
                f"{level} entity is missing products",
                entity_id=entity_id,
                missing=missing,
            )
    product_run_ids = {row["entity_id"] for row in run_rows}
    for run_id in selected:
        if run_id not in product_run_ids:
            add_diagnostic(
                diagnostics,
                "missing_run_products",
                f"selected run has no products: {run_id}",
                run_accession=run_id,
            )

    if args.mode == "final":
        for run_id in selected:
            submission = source_bams.get(run_id)
            if not submission:
                add_diagnostic(
                    diagnostics,
                    "missing_source_bam_submission",
                    f"no source-BAM submission for {run_id}",
                    run_accession=run_id,
                )
                continue
            status = first(submission, "status", "submission_status").lower()
            source_url = first(submission, "url", "ena_url", "file_url")
            if status not in {"submitted", "complete", "completed", "released"}:
                add_diagnostic(
                    diagnostics,
                    "source_bam_not_submitted",
                    f"source-BAM status is {status!r}",
                    run_accession=run_id,
                )
            if urlparse(source_url).scheme not in {"http", "https"}:
                add_diagnostic(
                    diagnostics,
                    "invalid_source_bam_url",
                    "source-BAM URL is not verified HTTP(S)",
                    run_accession=run_id,
                )
            selected_bam = first(
                selected[run_id], "bam", "source_bam_path", "source_bam_filename"
            )
            submitted_bam = first(
                submission, "source_bam_path", "source_bam_filename", "source_file"
            )
            if (
                selected_bam
                and submitted_bam
                and Path(selected_bam).name != Path(submitted_bam).name
            ):
                add_diagnostic(
                    diagnostics,
                    "source_bam_mismatch",
                    "source-BAM submission does not match selected source",
                    run_accession=run_id,
                    selected=selected_bam,
                    submitted=submitted_bam,
                )

    output: list[dict[str, Any]] = []
    for row in sorted(
        all_rows,
        key=lambda item: (
            item["selection_level"],
            item["entity_id"],
            item["product_slug"],
        ),
    ):
        if not row.get("_emit"):
            continue
        if row["selection_level"] == "run":
            selected_row = selected[row["run_accession"]]
            run = runs.get(row["run_accession"])
            sample = samples.get(selected_row["sample_accession"])
            if run is None or sample is None:
                continue
        else:
            run = sample = None
        output.append(
            render_record(
                row,
                args.species,
                assembly,
                args.label,
                args.team,
                args.user,
                args.selection_hierarchy,
                run,
                sample,
                source_bams,
                args.dummy_track_source,
            )
        )

    names: set[str] = set()
    for record in output:
        if record["name"] in names:
            fail(f"duplicate output name: {record['name']}")
        names.add(record["name"])
    report = {
        "mode": args.mode,
        "assembly": {
            "gca_accession": assembly.gca_accession,
            "assembly_uuid": assembly.assembly_uuid,
            "linked_genomes": list(assembly.linked_genomes),
            "evidence": assembly.evidence,
        },
        "selected_run_count": len(selected),
        "run_entity_count": len({row["entity_id"] for row in run_rows}),
        "tissue_entity_count": len({row["entity_id"] for row in aggregate_rows}),
        "run_record_count": sum(
            record["track_metadata_payload"]["selection_level"] == "run"
            for record in output
        ),
        "tissue_record_count": sum(
            record["track_metadata_payload"]["selection_level"] == "tissue"
            for record in output
        ),
        "output_record_count": len(output),
        "diagnostics": [diagnostic.as_dict() for diagnostic in diagnostics],
        "inputs": {
            str(path): input_hash(path)
            for path in (
                args.selected_runs,
                args.run_metadata,
                args.biosamples,
                args.run_manifest,
                args.aggregate_manifest,
                args.source_bams,
            )
            if path is not None
        },
    }
    return output, report, all_rows


def write_normalized(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "selection_level",
        "entity_id",
        "run_accession",
        "aggregate_id",
        "aggregate_label",
        "gca_accession",
        "track_type",
        "product_slug",
        "format",
        "track_file",
        "status",
        "sha256",
        "md5",
        "source_count",
        "source_run_accessions",
        "metadata_json",
        "provenance_json",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            output = dict(row)
            output["source_run_accessions"] = json.dumps(
                row.get("source_runs", []), sort_keys=True
            )
            output["metadata_json"] = json.dumps(
                row.get("metadata", {}), sort_keys=True
            )
            output["provenance_json"] = json.dumps(
                row.get("provenance", {}), sort_keys=True
            )
            writer.writerow(output)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-runs", type=Path, required=True)
    parser.add_argument("--run-metadata", type=Path, required=True)
    parser.add_argument("--biosamples", type=Path, required=True)
    parser.add_argument("--run-manifest", type=Path, required=True)
    parser.add_argument("--aggregate-manifest", type=Path)
    parser.add_argument("--source-bams", type=Path)
    parser.add_argument(
        "--metadata-uri",
        default=os.getenv("METADATA_URI"),
        help="metadata DB URI; defaults to METADATA_URI",
    )
    parser.add_argument("--mode", choices=("draft", "final"), default="draft")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--normalized-out", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    parser.add_argument("--species", default="Sus scrofa")
    parser.add_argument("--label", default="pig-pilot-20261104-r1")
    parser.add_argument("--team", default="genebuild")
    parser.add_argument("--user", default="jackt@ebi.ac.uk")
    parser.add_argument(
        "--selection-hierarchy",
        default="tissue,biosample,run",
        type=lambda raw: [item.strip() for item in raw.split(",") if item.strip()],
    )
    parser.add_argument(
        "--dummy-track-source",
        action="store_true",
        help="add draft-only example.org track_source placeholders",
    )
    parser.add_argument(
        "--run-products",
        default=",".join(DEFAULT_RUN_PRODUCTS),
        type=lambda raw: [item.strip() for item in raw.split(",") if item.strip()],
    )
    parser.add_argument(
        "--aggregate-products",
        default=",".join(DEFAULT_AGGREGATE_PRODUCTS),
        type=lambda raw: [item.strip() for item in raw.split(",") if item.strip()],
    )
    parser.add_argument(
        "--dataset-attribute", action="append", default=[], metavar="NAME=VALUE"
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.mode == "final" and args.source_bams is None:
        print("--source-bams is required in final mode", file=sys.stderr)
        return 2
    if args.mode == "final" and args.dummy_track_source:
        print("--dummy-track-source cannot be used in final mode", file=sys.stderr)
        return 2
    try:
        output, report, rows = build_handover(args)
        attributes = []
        for item in args.dataset_attribute:
            if "=" not in item:
                fail(f"--dataset-attribute must be NAME=VALUE: {item}")
            name, attr_value = item.split("=", 1)
            if not name or not attr_value:
                fail(f"--dataset-attribute must have non-empty name and value: {item}")
            attributes.append({"name": name, "value": attr_value})
        for record in output:
            if attributes:
                record["dataset_attribute"] = attributes
        write_json(args.out, output)
        write_normalized(args.normalized_out, rows)
        write_json(args.report_out, report)
        print(
            json.dumps(
                {
                    "output": str(args.out),
                    "normalized": str(args.normalized_out),
                    "report": str(args.report_out),
                    "records": len(output),
                    "diagnostics": len(report["diagnostics"]),
                },
                indent=2,
            )
        )
        return 1 if args.mode == "final" and report["diagnostics"] else 0
    except HandoverError as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

"""Validate an expression-track handover JSON file."""

# The outer contract intentionally checks several independent record invariants.
# pylint: disable=too-many-branches

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping
from urllib.parse import urlparse

from ensembl.genes.transcriptomic_data.expression_track_handover_schema import (
    TrackMetadataPayload,
)


DATASET_TYPE = "expression_tracks"


class HandoverValidationError(ValueError):
    """A handover record or file violates the outer contract."""


def text_value(raw: Any) -> str:
    """Return a stripped string for presence and URL checks."""
    return "" if raw is None else str(raw).strip()


def validate_record(record: Mapping[str, Any]) -> None:
    """Validate one outer handover record and its nested payload."""
    name = text_value(record.get("name")) or "<unnamed>"
    required = {
        "team_responsible",
        "user",
        "name",
        "label",
        "dataset_type",
        "dataset_source",
        "assembly_uuid",
        "track_files",
        "track_metadata_payload",
    }
    missing = sorted(required - set(record))
    if missing:
        raise HandoverValidationError(f"{name} is missing outer fields: {missing}")
    if record["dataset_type"] != DATASET_TYPE:
        raise HandoverValidationError(f"{name} has invalid dataset_type")

    try:
        payload = TrackMetadataPayload.model_validate(record["track_metadata_payload"])
    except ValueError as exc:
        raise HandoverValidationError(
            f"{name} has invalid track_metadata_payload: {exc}"
        ) from exc
    if name != payload.selection_id:
        raise HandoverValidationError(
            f"{name} does not match track_metadata_payload.selection_id"
        )
    if str(payload.target_genome.assembly_uuid) != str(record["assembly_uuid"]):
        raise HandoverValidationError(f"{name} has inconsistent assembly UUIDs")

    source = record["dataset_source"]
    if (
        not isinstance(source, dict)
        or not text_value(source.get("name"))
        or source.get("type") not in {"bigbed", "bigwig"}
    ):
        raise HandoverValidationError(f"{name} has invalid dataset_source")

    files = record["track_files"]
    if not isinstance(files, dict) or not files or set(files) - {"details", "summary"}:
        raise HandoverValidationError(f"{name} has invalid track_files")
    if not all(isinstance(path, str) and text_value(path) for path in files.values()):
        raise HandoverValidationError(f"{name} has empty or invalid track file paths")
    if source["type"] == "bigwig" and set(files) != {"summary"}:
        raise HandoverValidationError(f"{name} BigWig must map to summary")
    if source["type"] == "bigbed" and set(files) != {"details"}:
        raise HandoverValidationError(f"{name} BigBed must map to details")
    role = "summary" if source["type"] == "bigwig" else "details"
    if source["name"] != files[role]:
        raise HandoverValidationError(
            f"{name} dataset_source.name must match track_files.{role}"
        )

    track_source = record.get("track_source")
    if track_source is not None and (
        not isinstance(track_source, dict)
        or not all(
            text_value(track_source.get(key)) for key in ("name", "details", "url")
        )
        or urlparse(track_source["url"]).scheme not in {"http", "https"}
    ):
        raise HandoverValidationError(f"{name} has invalid track_source")

    try:
        json.dumps(record)
    except TypeError as exc:
        raise HandoverValidationError(
            f"{name} is not JSON serializable: {exc}"
        ) from exc


def validate_records(records: Any) -> list[str]:
    """Return all validation errors found in a handover record list."""
    if not isinstance(records, list):
        return ["handover JSON must contain a list of records"]
    errors: list[str] = []
    names: set[str] = set()
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            errors.append(f"record {index} is not an object")
            continue
        name = text_value(record.get("name")) or f"record {index}"
        if name in names:
            errors.append(f"{name} is duplicated")
        names.add(name)
        try:
            validate_record(record)
        except HandoverValidationError as exc:
            errors.append(f"record {index} ({name}): {exc}")
    return errors


def load_handover(path: Path) -> Any:
    """Read and decode a handover JSON file."""
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HandoverValidationError(f"cannot read {path}: {exc}") from exc


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the validator."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("handover", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Validate a handover file and print a machine-readable result."""
    args = parse_args(argv)
    try:
        records = load_handover(args.handover)
        errors = validate_records(records)
    except HandoverValidationError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    if errors:
        print(json.dumps({"valid": False, "errors": errors}, indent=2), file=sys.stderr)
        return 1
    print(json.dumps({"valid": True, "records": len(records)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Changelog comparison logic for project YAML files.

Compares a newly generated YAML against a previous version and reports
added, removed, and modified entries in a human-readable format.
"""

import re
from typing import Any

import yaml

# Fields that are internal / audit-only and should never appear in diffs
_IGNORE_PREFIXES = ("__audit_",)


def _normalise_value(val: Any) -> Any:
    """Normalise a value for comparison.

    - Strings: strip whitespace, trailing slashes; normalise date separators.
    - None / missing → empty string for consistent comparison.
    - Everything else: pass through unchanged.
    """
    if val is None:
        return ""
    if isinstance(val, str):
        v = val.strip().rstrip("/")
        # Normalise YYYY-MM ↔ YYYY_MM so date-only formatting diffs are ignored
        v = re.sub(r"(\d{4})[-_](\d{2})", r"\1_\2", v)
        return v
    return val


def _extract_key(doc: dict[str, Any]) -> str | None:
    """Return the comparison key for a YAML document.

    Priority: accession → assembly_accession → species → assembly (name).
    Returns None if nothing usable is found.
    """
    for field in ("accession", "assembly_accession", "species", "assembly"):
        val = doc.get(field)
        if val:
            return str(val).strip()
    return None


def _display_name(doc: dict[str, Any]) -> str:
    """Return a human-readable label for a YAML entry."""
    species = doc.get("species") or doc.get("assembly") or "unknown"
    return str(species)


def _should_ignore(field: str) -> bool:
    """True if the field should be excluded from comparison."""
    for prefix in _IGNORE_PREFIXES:
        if field.startswith(prefix):
            return True
    return False


def _diff_fields(
    old_doc: dict[str, Any], new_doc: dict[str, Any]
) -> list[tuple[str, Any, Any]]:
    """Return a list of (field, old_value, new_value) for fields that differ."""
    all_keys = sorted(set(old_doc.keys()) | set(new_doc.keys()))
    diffs: list[tuple[str, Any, Any]] = []
    for key in all_keys:
        if _should_ignore(key):
            continue
        old_val = _normalise_value(old_doc.get(key))
        new_val = _normalise_value(new_doc.get(key))
        if old_val != new_val:
            # Use original (un-normalised) values for display
            diffs.append((key, old_doc.get(key), new_doc.get(key)))
    return diffs


def load_yaml_as_keyed_dict(yaml_path: str) -> dict[str, dict[str, Any]]:
    """Load a multi-document YAML file and key entries by accession.

    Returns:
        Dict mapping accession (or fallback key) → document dict.
        Entries without a usable key are silently skipped.
    """
    keyed: dict[str, dict[str, Any]] = {}
    with open(yaml_path, "r", encoding="utf-8") as f:
        content = f.read()

    # The project YAMLs are written as a sequence of single-element lists,
    # so we try multi-document loading first, then fall back to a flat list.
    docs: list = []
    for raw_doc in yaml.safe_load_all(content):
        if raw_doc is None:
            continue
        if isinstance(raw_doc, list):
            docs.extend(raw_doc)
        elif isinstance(raw_doc, dict):
            docs.append(raw_doc)

    for doc in docs:
        if not isinstance(doc, dict):
            continue
        key = _extract_key(doc)
        if key:
            keyed[key] = doc
    return keyed


def compare_yamls(
    old_docs: dict[str, dict[str, Any]],
    new_docs: dict[str, dict[str, Any]],
) -> tuple[
    list[tuple[str, str]],
    list[tuple[str, str]],
    list[tuple[str, str, list[tuple[str, Any, Any]]]],
]:
    """Compare old and new YAML documents keyed by accession.

    Returns:
        (added, removed, modified) where:
          added   = [(key, display_name), ...]
          removed = [(key, display_name), ...]
          modified = [(key, display_name, [(field, old, new), ...]), ...]
    """
    old_keys = set(old_docs.keys())
    new_keys = set(new_docs.keys())

    added = sorted([(k, _display_name(new_docs[k])) for k in (new_keys - old_keys)])
    removed = sorted([(k, _display_name(old_docs[k])) for k in (old_keys - new_keys)])

    modified: list[tuple[str, str, list[tuple[str, Any, Any]]]] = []
    for key in sorted(old_keys & new_keys):
        diffs = _diff_fields(old_docs[key], new_docs[key])
        if diffs:
            name = _display_name(new_docs[key])
            modified.append((key, name, diffs))

    return added, removed, modified


def format_changelog(
    added: list[tuple[str, str]],
    removed: list[tuple[str, str]],
    modified: list[tuple[str, str, list[tuple[str, Any, Any]]]],
) -> str:
    """Format the changelog as a human-readable string."""
    lines: list[str] = []
    lines.append("")
    lines.append("=== CHANGELOG SUMMARY ===")
    lines.append("")

    # Added
    lines.append(f"Added ({len(added)}):")
    if added:
        for key, name in added:
            lines.append(f"  - {key} ({name})")
    else:
        lines.append("  (none)")
    lines.append("")

    # Removed
    lines.append(f"Removed ({len(removed)}):")
    if removed:
        for key, name in removed:
            lines.append(f"  - {key} ({name})")
    else:
        lines.append("  (none)")
    lines.append("")

    # Modified
    lines.append(f"Modified ({len(modified)}):")
    if modified:
        for key, name, diffs in modified:
            lines.append(f"  - {key} ({name})")
            for field, old_val, new_val in diffs:
                lines.append(f"    - {field}:")
                lines.append(f"        OLD: {old_val}")
                lines.append(f"        NEW: {new_val}")
    else:
        lines.append("  (none)")
    lines.append("")

    # Totals
    total_changes = len(added) + len(removed) + len(modified)
    if total_changes == 0:
        lines.append("No differences detected.")
    else:
        lines.append(
            f"Total: {len(added)} added, {len(removed)} removed, "
            f"{len(modified)} modified."
        )
    lines.append("")
    return "\n".join(lines)


def write_changelog_tsv(
    tsv_path: str,
    added: list[tuple[str, str]],
    removed: list[tuple[str, str]],
    modified: list[tuple[str, str, list[tuple[str, Any, Any]]]],
) -> None:
    """Write a machine-readable TSV changelog."""
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("status\taccession\tspecies\tfield\told_value\tnew_value\n")
        f.writelines(f"added\t{key}\t{name}\t\t\t\n" for key, name in added)
        f.writelines(f"removed\t{key}\t{name}\t\t\t\n" for key, name in removed)
        for key, name, diffs in modified:
            f.writelines(
                f"modified\t{key}\t{name}\t{field}\t"
                f"{old_val if old_val is not None else ''}\t"
                f"{new_val if new_val is not None else ''}\n"
                for field, old_val, new_val in diffs
            )

"""
Checksums, implementation fingerprints and environment records.

The implementation fingerprint identifies the comparator code that produced a
result. It hashes the *content* of the source files on the pairwise-compare path
(not the Git commit), so uncommitted changes are captured. A Git commit alone was
not sufficient provenance for the existing benchmark, which ran on uncommitted
code; its file hashes were captured separately and can be supplied as a
"code snapshot" (``sha256sum``-style lines) when importing.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent.parent  # .../annotation_qc
PACKAGE = "ensembl.genes.annotation_qc"
# The comparator entry point; every in-package module it imports, directly or
# transitively, determines pairwise-compare results.
FINGERPRINT_ENTRY = "runners/pairwise_compare.py"
FINGERPRINT_PACKAGES = ("pandas", "pyranges1")


def _module_file(module: str) -> str | None:
    """'ensembl.genes.annotation_qc.a.b' -> 'a/b.py' (or 'a/b/__init__.py'); None outside the package."""
    if module != PACKAGE and not module.startswith(PACKAGE + "."):
        return None
    parts = module[len(PACKAGE) + 1 :].split(".") if module != PACKAGE else []
    for rel in ("/".join(parts) + ".py", "/".join(parts + ["__init__.py"])):
        if rel != ".py" and (PACKAGE_DIR / rel).is_file():
            return rel
    return None


def _has_code(rel: str) -> bool:
    """True unless the file holds nothing but a docstring (empty package markers)."""
    body = ast.parse((PACKAGE_DIR / rel).read_text()).body
    if (
        body
        and isinstance(body[0], ast.Expr)
        and isinstance(body[0].value, ast.Constant)
    ):
        body = body[1:]
    return bool(body)


def _imports(rel: str) -> set[str]:
    """In-package files imported by one file (any depth in the file, including lazy imports)."""
    found = set()
    for node in ast.walk(ast.parse((PACKAGE_DIR / rel).read_text())):
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            # "from pkg import mod" may name a submodule or an attribute; try both
            names = [node.module] + [f"{node.module}.{a.name}" for a in node.names]
        else:
            continue
        for name in names:
            target = _module_file(name)
            if target:
                found.add(target)
            # importing a.b.c executes the __init__ of a and a.b as well
            parts = name.split(".")
            for i in range(len(PACKAGE.split(".")), len(parts)):
                init = _module_file(".".join(parts[:i]))
                if init and init.endswith("__init__.py"):
                    found.add(init)
    return found


def comparator_source_files(entry: str = FINGERPRINT_ENTRY) -> tuple[str, ...]:
    """
    Files whose content determines pairwise-compare results (relative to annotation_qc/).

    The transitive closure of in-package imports from the runner, found by parsing
    the sources, so a new comparator module is covered as soon as the runner (or
    anything it imports) imports it. Package ``__init__`` files that contain only a
    docstring are left out: they cannot change results. ``cli.py`` only dispatches
    and is not imported by the runner, so adding commands does not invalidate results.
    """
    seen, todo = set(), [entry]
    while todo:
        rel = todo.pop()
        if rel in seen:
            continue
        seen.add(rel)
        todo.extend(_imports(rel) - seen)
    return tuple(sorted(rel for rel in seen if _has_code(rel)))


FINGERPRINT_FILES = comparator_source_files()


def sha256_file(path: str | Path, block: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(block), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def stable_hash(obj) -> str:
    """SHA-256 of a JSON-serialisable object with sorted keys."""
    return sha256_text(json.dumps(obj, sort_keys=True, default=str))


class ChecksumCache:
    """
    sha256 per file, cached by (absolute path, size, mtime_ns) so large inputs are
    hashed once. A changed size or modification time forces re-hashing.
    """

    def __init__(self, cache_file: Path):
        self.cache_file = Path(cache_file)
        try:
            self.entries = json.loads(self.cache_file.read_text())
        except (OSError, ValueError):
            self.entries = {}

    def record(self, path: str | Path) -> dict:
        path = Path(path).resolve()
        stat = path.stat()
        key = str(path)
        entry = self.entries.get(key)
        if (
            not entry
            or entry["size"] != stat.st_size
            or entry["mtime_ns"] != stat.st_mtime_ns
        ):
            entry = {
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": sha256_file(path),
            }
            self.entries[key] = entry
            self.save()
        return {"path": key, "size": entry["size"], "sha256": entry["sha256"]}

    def sha256(self, path: str | Path) -> str:
        return self.record(path)["sha256"]

    def save(self) -> None:
        self.cache_file.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.cache_file.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.entries, indent=1, sort_keys=True))
        tmp.replace(self.cache_file)


def package_versions() -> dict:
    versions = {"python": platform.python_version()}
    for name in FINGERPRINT_PACKAGES:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def fingerprint_from_hashes(file_hashes: dict[str, str | None], versions: dict) -> str:
    return stable_hash(
        {
            "files": file_hashes,
            "versions": {
                k: versions.get(k) for k in ("python",) + FINGERPRINT_PACKAGES
            },
        }
    )


def current_implementation() -> dict:
    """Content hashes of the comparator files in this checkout, plus versions and Git state."""
    hashes = {}
    for rel in FINGERPRINT_FILES:
        path = PACKAGE_DIR / rel
        hashes[rel] = sha256_file(path) if path.exists() else None
    versions = package_versions()
    return {
        "fingerprint": fingerprint_from_hashes(hashes, versions),
        "files": hashes,
        "versions": versions,
        "git": git_state(),
    }


def implementation_from_snapshot(
    snapshot_files: list[Path], versions: dict
) -> dict | None:
    """
    Fingerprint of a past implementation from ``sha256sum``-style snapshot files.

    Each line is "<sha256>  <path>"; paths are matched on the part after
    ``annotation_qc/``. Returns None when a fingerprint file is missing from the
    snapshot (the implementation cannot then be verified).
    """
    found: dict[str, str] = {}
    for snapshot in snapshot_files:
        for line in Path(snapshot).read_text().splitlines():
            parts = line.split(maxsplit=1)
            if len(parts) != 2 or "annotation_qc/" not in parts[1]:
                continue
            rel = parts[1].strip().split("annotation_qc/", 1)[1]
            found[rel] = parts[0]
    hashes = {rel: found.get(rel) for rel in FINGERPRINT_FILES}
    missing = [rel for rel, digest in hashes.items() if digest is None]
    if missing:
        return {
            "fingerprint": None,
            "files": hashes,
            "versions": versions,
            "missing_from_snapshot": missing,
        }
    return {
        "fingerprint": fingerprint_from_hashes(hashes, versions),
        "files": hashes,
        "versions": versions,
    }


def git_state() -> dict:
    source = str(PACKAGE_DIR)
    try:
        commit = subprocess.run(
            ["git", "-C", source, "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        branch = subprocess.run(
            ["git", "-C", source, "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", source, "status", "--porcelain", "--", "."],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        return {
            "commit": commit,
            "branch": branch,
            "uncommitted_changes_in_annotation_qc": bool(dirty),
        }
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None}


def environment() -> dict:
    return {
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "versions": package_versions(),
        "ensembl_genes": _safe_version("ensembl-genes"),
        "recorded": datetime.now(timezone.utc).isoformat(),
        "cwd": os.getcwd(),
    }


def _safe_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")

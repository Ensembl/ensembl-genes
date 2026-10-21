"""
Implementation provenance of comparator results: which code produced a manifest.

A comparator manifest (contract 1.0) records ``implementation``: a fingerprint
scheme, a fingerprint, the SHA-256 of every annotation_qc file it covers
(``source_files``) and the pandas / pyranges1 / Python versions. A build is
identified by the pair (scheme, fingerprint); fingerprints of different schemes
are never compared with each other.

    annotation_qc_loaded_modules_by_name/1   manifests without a ``scheme`` key;
                                             modules selected by import name, so the
                                             runner loaded by the ``annotation-qc``
                                             entry point is not covered
    annotation_qc_loaded_files_by_path/2     every loaded module file in the package

Identification of one manifest (``identify``):

    verified          (scheme, fingerprint) is in a registry, the fingerprint
                      recomputes from the listed file hashes, and those hashes are
                      the registered build's. The result is attributed to that
                      (historical) build, whatever code the dashboard runs now.
    current_checkout  scheme 2, not registered, but every listed file is
                      byte-identical to the file in this checkout.
    unverified        anything else (no implementation block, unknown scheme,
                      fingerprint inconsistent with its file list, unknown pair,
                      or file hashes that differ from the registered build).

Registries are JSON files (``implementation_registry.json`` in this package plus
any listed in the experiment configuration under ``provenance.registries``). Each
build entry carries its evidence (commit whose git objects reproduce every file
hash, wheels, verification reports). ``register_build`` writes an entry only
after checking every manifest's file hashes against ``git show <commit>:<file>``.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent.parent  # .../annotation_qc
PACKAGE_REPO_PATH = "src/python/ensembl/genes/annotation_qc"
PACKAGED_REGISTRY = Path(__file__).with_name("implementation_registry.json")

SCHEME_BY_NAME = "annotation_qc_loaded_modules_by_name/1"
SCHEME_BY_PATH = "annotation_qc_loaded_files_by_path/2"
SCHEMES = (SCHEME_BY_NAME, SCHEME_BY_PATH)

VERIFIED = "verified"
CURRENT_CHECKOUT = "current_checkout"
UNVERIFIED = "unverified"


class RegistryError(ValueError):
    """A registry file is malformed or contradicts another registry."""


def scheme_of(implementation: dict) -> str:
    """Scheme recorded in a manifest's implementation block (absent = scheme 1)."""
    return implementation.get("scheme") or SCHEME_BY_NAME


def recompute_fingerprint(source_files: dict, versions: dict) -> str:
    """The comparator's fingerprint formula (identical in schemes 1 and 2)."""
    text = json.dumps({"files": source_files, "versions": versions}, sort_keys=True)
    return hashlib.sha256(text.encode()).hexdigest()


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return None


@dataclass
class Registry:
    """Known builds keyed by (scheme, fingerprint)."""

    builds: dict[tuple[str, str], dict] = field(default_factory=dict)
    sources: list[str] = field(default_factory=list)

    @classmethod
    def load(cls, extra: list[Path] | None = None) -> "Registry":
        registry = cls()
        for path in [PACKAGED_REGISTRY, *(extra or [])]:
            registry.add_file(Path(path))
        return registry

    def add_file(self, path: Path) -> None:
        try:
            data = json.loads(Path(path).read_text())
        except (OSError, ValueError) as error:
            raise RegistryError(f"cannot read registry {path}: {error}") from error
        for entry in data.get("builds", []):
            self.add(entry, str(path))
        self.sources.append(str(path))

    def add(self, entry: dict, source: str = "") -> None:
        missing = [
            k
            for k in ("build_id", "scheme", "fingerprint", "source_files", "versions")
            if not entry.get(k)
        ]
        if missing:
            raise RegistryError(f"{source}: build entry lacks {missing}")
        if entry["scheme"] not in SCHEMES:
            raise RegistryError(f"{source}: unknown scheme {entry['scheme']!r}")
        if (
            recompute_fingerprint(entry["source_files"], entry["versions"])
            != entry["fingerprint"]
        ):
            raise RegistryError(
                f"{source}: build {entry['build_id']}: fingerprint does not match its file hashes"
            )
        key = (entry["scheme"], entry["fingerprint"])
        known = self.builds.get(key)
        if known and known["source_files"] != entry["source_files"]:
            raise RegistryError(
                f"{source}: {key} registered twice with different files"
            )
        self.builds.setdefault(key, {**entry, "registry": source})

    def lookup(self, scheme: str, fingerprint: str) -> dict | None:
        return self.builds.get((scheme, fingerprint))

    def public(self) -> list[dict]:
        return [
            {k: v for k, v in entry.items() if k != "source_files"}
            | {"files": len(entry["source_files"])}
            for entry in self.builds.values()
        ]


_CHECKOUT_CACHE: dict[str, str | None] = {}


def checkout_hash(rel: str) -> str | None:
    """SHA-256 of a file of this checkout's annotation_qc package (cached per process)."""
    if rel not in _CHECKOUT_CACHE:
        _CHECKOUT_CACHE[rel] = _sha256(PACKAGE_DIR / rel)
    return _CHECKOUT_CACHE[rel]


def identify(manifest: dict | None, registry: Registry) -> dict:
    """
    Implementation identity of one comparator manifest (see module docstring).

    Returns {status, scheme, fingerprint, build_id?, commit?, reason, ...}.
    """
    impl = (manifest or {}).get("implementation")
    if not impl or not impl.get("fingerprint"):
        return {
            "status": UNVERIFIED,
            "scheme": None,
            "fingerprint": None,
            "reason": "the manifest records no implementation fingerprint (pre-contract result)",
        }
    scheme, fingerprint = scheme_of(impl), impl["fingerprint"]
    out = {
        "scheme": scheme,
        "scheme_recorded": "scheme" in impl,
        "fingerprint": fingerprint,
        "versions": impl.get("versions"),
        "files": len(impl.get("source_files") or {}),
    }
    if scheme not in SCHEMES:
        return {**out, "status": UNVERIFIED, "reason": f"unsupported scheme {scheme}"}
    files, versions = impl.get("source_files") or {}, impl.get("versions") or {}
    if recompute_fingerprint(files, versions) != fingerprint:
        return {
            **out,
            "status": UNVERIFIED,
            "reason": "the fingerprint does not recompute from the file hashes the manifest lists",
        }
    entry = registry.lookup(scheme, fingerprint)
    if entry is not None:
        if entry["source_files"] != files:
            return {
                **out,
                "status": UNVERIFIED,
                "reason": f"file hashes differ from registered build {entry['build_id']}",
            }
        return {
            **out,
            "status": VERIFIED,
            "build_id": entry["build_id"],
            "commit": entry.get("commit"),
            "contract_version": entry.get("contract_version"),
            "definition_version": entry.get("definition_version"),
            "evidence": entry.get("evidence", []),
            "notes": entry.get("notes"),
            "reason": f"registered build {entry['build_id']} ({scheme})",
        }
    if scheme == SCHEME_BY_PATH and files:
        differing = sorted(rel for rel, h in files.items() if checkout_hash(rel) != h)
        if not differing:
            return {
                **out,
                "status": CURRENT_CHECKOUT,
                "reason": "not registered; every listed file is identical to this checkout",
            }
        return {
            **out,
            "status": UNVERIFIED,
            "reason": f"not registered and {len(differing)} file(s) differ from this checkout "
            f"(e.g. {differing[0]})",
        }
    return {
        **out,
        "status": UNVERIFIED,
        "reason": "(scheme, fingerprint) is not in any registry",
    }


def label(identity: dict) -> str:
    """Short display text for an identity."""
    if identity["status"] == VERIFIED:
        return f"verified implementation {identity['build_id']}"
    if identity["status"] == CURRENT_CHECKOUT:
        return "current code (this checkout)"
    return "unverified implementation"


def _git_blob_hash(repo: Path, commit: str, rel: str) -> str | None:
    proc = subprocess.run(
        ["git", "-C", str(repo), "show", f"{commit}:{PACKAGE_REPO_PATH}/{rel}"],
        capture_output=True,
        check=False,
    )
    return hashlib.sha256(proc.stdout).hexdigest() if proc.returncode == 0 else None


def register_build(
    manifests: list[Path],
    repo: Path,
    commit: str,
    build_id: str,
    evidence: list[dict] | None = None,
    notes: str | None = None,
) -> dict:
    """
    A registry entry for the build behind ``manifests``, checked against git.

    Every manifest must record the same (scheme, fingerprint, files, versions),
    its fingerprint must recompute, and every file hash must equal the SHA-256 of
    that file in ``commit``.
    Raises:
            RegistryError naming the first inconsistency
    """
    seen = None
    for path in manifests:
        manifest = json.loads(Path(path).read_text())
        impl = manifest.get("implementation") or {}
        key = (
            scheme_of(impl),
            impl.get("fingerprint"),
            json.dumps(impl.get("source_files"), sort_keys=True),
            json.dumps(impl.get("versions"), sort_keys=True),
            manifest.get("contract_version"),
            manifest.get("definition_version"),
        )
        if seen is None:
            seen = (key, impl, manifest)
        elif key != seen[0]:
            raise RegistryError(f"{path}: implementation differs from {manifests[0]}")
    if seen is None:
        raise RegistryError("no manifests given")
    _, impl, manifest = seen
    if (
        recompute_fingerprint(impl["source_files"], impl["versions"])
        != impl["fingerprint"]
    ):
        raise RegistryError("fingerprint does not recompute from its file hashes")
    differing = [
        rel
        for rel, digest in impl["source_files"].items()
        if _git_blob_hash(repo, commit, rel) != digest
    ]
    if differing:
        raise RegistryError(f"files differ from commit {commit}: {differing}")
    full = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", commit],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    return {
        "build_id": build_id,
        "commit": full,
        "scheme": scheme_of(impl),
        "fingerprint": impl["fingerprint"],
        "versions": impl["versions"],
        "source_files": impl["source_files"],
        "contract_id": manifest.get("contract_id"),
        "contract_version": manifest.get("contract_version"),
        "definition_version": manifest.get("definition_version"),
        "evidence": [
            {
                "type": "git_objects",
                "detail": f"every listed file hash equals the file in commit {full}",
                "manifests_checked": len(manifests),
            },
            *(evidence or []),
        ],
        "notes": notes,
    }

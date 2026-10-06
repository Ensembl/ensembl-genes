"""
Recorded, reproducible input preparation.

Preparation never changes biology: the only steps are

    decompress               gzip/BGZF -> plain (exact bytes)
    retype_transcript_types  column 3 of the listed transcript-level SO types ->
                             ``transcript``, appending ``original_type=<type>`` to
                             column 9. IDs, Parent links, coordinates and every
                             other record are unchanged. Needed while the parser
                             does not list Ensembl ``unconfirmed_transcript`` /
                             ``gene_segment`` (see validation).

Outputs go to ``<output_dir>/prepared/``; sources are only read. Each output has a
``preparation.json`` with the source checksum, recipe and output checksum.
``reuse_existing`` adopts an already prepared file only after re-deriving it from
the source in memory and confirming the checksum is identical.
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path

from ensembl.genes.annotation_qc.benchmark.config import Preparation, Reference, Run
from ensembl.genes.annotation_qc.benchmark.provenance import now, sha256_text
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace, write_json

PREPARATION_VERSION = "1"


class PreparationError(RuntimeError):
    pass


def is_gzip(path: Path) -> bool:
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _open_text(path: Path):
    return gzip.open(path, "rt") if is_gzip(path) else open(path, "rt")


def _retype_lines(lines, types: set[str], to: str):
    for line in lines:
        if not line.startswith("#"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] in types:
                parts[8] = parts[8].rstrip(";") + f";original_type={parts[2]}"
                parts[2] = to
                line = "\t".join(parts) + "\n"
        yield line


def _transform(source: Path, steps: list[dict], sink, counts: dict) -> None:
    """Stream the prepared content of ``source`` into ``sink(bytes)``."""
    names = [s["name"] for s in steps]
    retype = [s for s in steps if s["name"] == "retype_transcript_types"]
    if "decompress" in names and not is_gzip(source):
        counts["decompress"] = "source not compressed; no-op"
    if not retype:
        opener = gzip.open if is_gzip(source) else open
        with opener(source, "rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                sink(chunk)
        return
    types = set(retype[0]["types"])
    to = retype[0].get("to", "transcript")
    retyped = {t: 0 for t in types}
    with _open_text(source) as handle:
        for line in _retype_lines(handle, types, to):
            if "original_type=" in line and not line.startswith("#"):
                original = line.rsplit("original_type=", 1)[1].strip()
                if original in retyped:
                    retyped[original] += 1
            sink(line.encode())
    counts["retyped_rows"] = retyped


def _needs_output(source: Path, preparation: Preparation) -> bool:
    names = [s["name"] for s in preparation.steps]
    if not names:
        return False
    if names == ["decompress"] and not is_gzip(source):
        return False
    return True


def prepare_file(
    ws: Workspace,
    source: Path,
    preparation: Preparation,
    out_dir: Path,
    out_name: str,
    kind: str,
    force=False,
) -> dict:
    """
    Prepare one file. Returns a record with source/output paths and checksums.
    Unchanged sources and recipes are not prepared again.
    """
    source = Path(source)
    if not source.exists():
        raise PreparationError(f"{kind}: source file not found: {source}")
    source_rec = ws.checksums.record(source)
    recipe = preparation.recipe()
    record_path = out_dir / f"{out_name}.preparation.json"
    previous = _read(record_path)
    if not _needs_output(source, preparation):
        record = {
            "kind": kind,
            "source": source_rec,
            "recipe": recipe,
            "output": source_rec,
            "method": "used as-is (no preparation needed)",
            "preparation_version": PREPARATION_VERSION,
            "created": now(),
        }
        write_json(record_path, record)
        return record
    if (
        not force
        and previous
        and previous.get("source", {}).get("sha256") == source_rec["sha256"]
        and previous.get("recipe") == recipe
        and previous.get("preparation_version") == PREPARATION_VERSION
        and Path(previous["output"]["path"]).exists()
        and ws.checksums.sha256(previous["output"]["path"])
        == previous["output"]["sha256"]
    ):
        return previous

    counts: dict = {}
    if preparation.reuse_existing:
        existing = Path(preparation.reuse_existing)
        if not existing.exists():
            raise PreparationError(f"{kind}: reuse_existing file not found: {existing}")
        digest = hashlib.sha256()
        _transform(source, preparation.steps, digest.update, counts)
        expected = digest.hexdigest()
        existing_rec = ws.checksums.record(existing)
        if existing_rec["sha256"] != expected:
            raise PreparationError(
                f"{kind}: reuse_existing {existing} does not match the recipe applied to {source} "
                f"(sha256 {existing_rec['sha256'][:12]}… vs re-derived {expected[:12]}…); not used"
            )
        record = {
            "kind": kind,
            "source": source_rec,
            "recipe": recipe,
            "output": existing_rec,
            "counts": counts,
            "method": "reused existing file; verified by re-deriving its sha256 from the source",
            "preparation_version": PREPARATION_VERSION,
            "created": now(),
        }
        write_json(record_path, record)
        return record

    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / out_name
    tmp = output.with_name(output.name + ".partial")
    with open(tmp, "wb") as handle:
        _transform(source, preparation.steps, handle.write, counts)
    tmp.replace(output)
    record = {
        "kind": kind,
        "source": source_rec,
        "recipe": recipe,
        "output": ws.checksums.record(output),
        "counts": counts,
        "method": "derived",
        "preparation_version": PREPARATION_VERSION,
        "created": now(),
    }
    write_json(record_path, record)
    return record


def sequence_region_names(path: Path) -> list[str]:
    """Names on ``##sequence-region`` lines, in file order (whole file scanned)."""
    names = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 2:
                    names.append(parts[1])
    return names


def prepare_scope(ws: Workspace, reference: Reference, annotation_path: Path) -> dict:
    """Write the scope regions file (whole sequences, one per line) or record 'all'."""
    scope = reference.scope
    out_dir = ws.prepared_dir(reference.id)
    kind = scope.get("type", "all")
    if kind == "all":
        return {"type": "all", "regions_file": None, "sha256": "all", "sequences": None}
    if kind == "regions_file":
        path = Path(scope["path"])
        rec = ws.checksums.record(path)
        sequences = [
            l.strip()
            for l in path.read_text().splitlines()
            if l.strip() and not l.startswith("#")
        ]
        return {
            "type": kind,
            "regions_file": rec["path"],
            "sha256": rec["sha256"],
            "sequences": sequences,
        }
    if kind == "sequences":
        sequences = list(scope["sequences"])
    else:  # reference_sequence_regions
        cached = _read(out_dir / "scope.json")
        ann_sha = ws.checksums.sha256(annotation_path)
        if cached and cached.get("annotation_sha256") == ann_sha:
            sequences = cached["sequences"]
        else:
            sequences = sequence_region_names(annotation_path)
            if not sequences:
                raise PreparationError(
                    f"{reference.id}: scope 'reference_sequence_regions' but the annotation has no ##sequence-region lines"
                )
            write_json(
                out_dir / "scope.json",
                {"annotation_sha256": ann_sha, "sequences": sequences},
            )
    text = "".join(f"{name}\n" for name in sequences)
    path = out_dir / "scope.regions"
    if not path.exists() or path.read_text() != text:
        out_dir.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return {
        "type": kind,
        "regions_file": str(path),
        "sha256": sha256_text(text),
        "sequences": sequences,
    }


def prepare_reference(ws: Workspace, reference: Reference, force=False) -> dict:
    out_dir = ws.prepared_dir(reference.id)
    annotation = prepare_file(
        ws,
        reference.annotation,
        reference.annotation_preparation,
        out_dir,
        "annotation" + _annotation_suffix(reference.annotation),
        f"reference {reference.id} annotation",
        force,
    )
    genome = None
    if reference.genome:
        genome = prepare_file(
            ws,
            reference.genome,
            reference.genome_preparation,
            out_dir,
            "genome.fa",
            f"reference {reference.id} genome",
            force,
        )
    scope = prepare_scope(ws, reference, Path(annotation["output"]["path"]))
    seqmap = (
        ws.checksums.record(reference.seqname_map) if reference.seqname_map else None
    )
    record = {
        "reference_id": reference.id,
        "annotation": annotation,
        "genome": genome,
        "scope": scope,
        "seqname_map": seqmap,
        "prepared": now(),
    }
    write_json(out_dir / "preparation.json", record)
    return record


def prepare_run(ws: Workspace, run: Run, force=False) -> dict:
    out_dir = ws.prepared_run_dir(run.run_id)
    rec = prepare_file(
        ws,
        run.annotation,
        run.preparation,
        out_dir,
        "annotation" + _annotation_suffix(run.annotation),
        f"run {run.run_id} annotation",
        force,
    )
    if rec["method"] != "used as-is (no preparation needed)":
        write_json(out_dir / "preparation.json", rec)
    return rec


def _annotation_suffix(path: Path) -> str:
    suffixes = [
        s.lower()
        for s in Path(path).suffixes
        if s.lower() not in (".gz", ".bgz", ".bgzf")
    ]
    return (
        suffixes[-1]
        if suffixes and suffixes[-1] in (".gff", ".gff3", ".gtf")
        else ".gff3"
    )


def _read(path: Path):
    from ensembl.genes.annotation_qc.benchmark.workspace import read_json

    return read_json(path)


def copy_if_missing(source: Path, target: Path) -> None:
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

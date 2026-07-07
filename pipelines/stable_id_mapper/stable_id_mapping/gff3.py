"""GFF3 parsing for stable-ID decision generation."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .models import CdsRow, Feature, RawFeature

GENE_FEATURE_TYPES = {"gene"}
TRANSCRIPT_FEATURE_TYPES = {"mrna", "transcript"}


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if not attr_text or attr_text == ".":
        return attrs
    for item in attr_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
        else:
            attrs[item] = ""
    return attrs


def split_stable_id(raw: Optional[str]) -> tuple[Optional[str], Optional[int]]:
    if not raw:
        return None, None
    token = raw.strip()
    if not token:
        return None, None
    if ":" in token:
        token = token.split(":", 1)[1]
    if "." in token:
        stable_id, maybe_version = token.rsplit(".", 1)
        if maybe_version.isdigit():
            return stable_id, int(maybe_version)
    return token, None


def parent_ids(raw: Optional[str]) -> tuple[tuple[str, Optional[int]], ...]:
    if not raw:
        return ()
    out: list[tuple[str, Optional[int]]] = []
    for token in raw.split(","):
        stable_id, version = split_stable_id(token)
        if stable_id:
            out.append((stable_id, version))
    return tuple(out)


def parse_version(
    attrs: dict[str, str],
    embedded_version: Optional[int],
    default: int = 1,
) -> int:
    version_text = attrs.get("version")
    if version_text is not None:
        try:
            return int(version_text)
        except ValueError:
            pass
    if embedded_version is not None:
        return embedded_version
    return default


def parse_gff3(path: str | Path) -> dict[str, dict[str, Feature]]:
    features: dict[str, dict[str, Feature]] = {
        "gene": {},
        "transcript": {},
        "translation": {},
    }
    raw_features: list[RawFeature] = []
    cds_rows: list[CdsRow] = []

    with Path(path).open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            seqid, _source, feature_type, start, end, _score, strand, _phase, attrs_text = fields
            attrs = parse_attrs(attrs_text)
            start_i = int(start)
            end_i = int(end)
            feature_type_lc = feature_type.lower()

            parents = parent_ids(attrs.get("Parent"))
            if feature_type_lc == "cds":
                protein_stable_id, protein_version = split_stable_id(attrs.get("protein_id"))
                for parent_stable_id, parent_version in parents:
                    cds_rows.append(
                        CdsRow(
                            parent_stable_id=parent_stable_id,
                            parent_version=parent_version,
                            translation_stable_id=protein_stable_id or parent_stable_id,
                            translation_version=protein_version
                            if protein_stable_id else parent_version,
                            seqid=seqid,
                            start=start_i,
                            end=end_i,
                            strand=strand,
                        )
                    )
                continue

            if feature_type_lc not in GENE_FEATURE_TYPES | TRANSCRIPT_FEATURE_TYPES:
                continue

            stable_id, embedded_version = split_stable_id(attrs.get("ID"))
            if not stable_id:
                continue
            raw_features.append(
                RawFeature(
                    stable_id=stable_id,
                    version=parse_version(attrs, embedded_version),
                    feature_type=feature_type,
                    seqid=seqid,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    parent_ids=tuple(parent for parent, _version in parents),
                )
            )

    gene_ids = {
        record.stable_id
        for record in raw_features
        if not record.parent_ids and record.feature_type.lower() in GENE_FEATURE_TYPES
    }
    for record in raw_features:
        if record.stable_id not in gene_ids:
            continue
        features["gene"][record.stable_id] = Feature(
            stable_id=record.stable_id,
            version=record.version,
            seqid=record.seqid,
            start=record.start,
            end=record.end,
            strand=record.strand,
        )

    transcript_ids: set[str] = set()
    for record in raw_features:
        if record.feature_type.lower() not in TRANSCRIPT_FEATURE_TYPES:
            continue
        parent_gene = next(
            (parent for parent in record.parent_ids if parent in gene_ids),
            None,
        )
        if parent_gene is None:
            continue
        transcript_ids.add(record.stable_id)
        features["transcript"][record.stable_id] = Feature(
            stable_id=record.stable_id,
            version=record.version,
            seqid=record.seqid,
            start=record.start,
            end=record.end,
            strand=record.strand,
            parent_stable_id=parent_gene,
        )

    for cds in cds_rows:
        if cds.parent_stable_id not in transcript_ids:
            continue
        transcript = features["transcript"][cds.parent_stable_id]
        version = (
            cds.translation_version if cds.translation_version is not None else transcript.version
        )
        previous = features["translation"].get(cds.translation_stable_id)
        if previous is None:
            features["translation"][cds.translation_stable_id] = Feature(
                stable_id=cds.translation_stable_id,
                version=version,
                seqid=cds.seqid,
                start=cds.start,
                end=cds.end,
                strand=cds.strand,
                parent_stable_id=cds.parent_stable_id,
            )
        else:
            features["translation"][cds.translation_stable_id] = Feature(
                stable_id=previous.stable_id,
                version=previous.version,
                seqid=previous.seqid,
                start=min(previous.start, cds.start),
                end=max(previous.end, cds.end),
                strand=previous.strand,
                parent_stable_id=previous.parent_stable_id,
            )

    return features


def parse_missing_gene_ids(report_path: str | Path) -> set[str]:
    missing: set[str] = set()
    in_missing_section = False
    with Path(report_path).open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped == "Missing gene IDs:":
                in_missing_section = True
                continue
            if not in_missing_section:
                continue
            if not stripped:
                continue
            if stripped.endswith(":") and not stripped.startswith("gene:"):
                break
            stable_id, _version = split_stable_id(stripped)
            if stable_id:
                missing.add(stable_id)
    return missing


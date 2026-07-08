#!/usr/bin/env python3
"""Pick real stable-ID truth examples from old/new GFF3 files or FTP releases."""

from __future__ import annotations

import argparse
import csv
import gzip
import html.parser
import io
import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from stable_id_mapping.gff3 import parse_attrs, parse_version, split_stable_id

FEATURE_TYPES = {
    "gene": {"gene"},
    "transcript": {"mrna", "transcript"},
}


@dataclass(frozen=True)
class RealFeature:
    stable_id: str
    version: int
    feature_type: str
    seqid: str
    start: int
    end: int
    strand: str
    source: str

    @property
    def location(self) -> str:
        return f"{self.seqid}:{self.start}-{self.end}({self.strand})"


@dataclass(frozen=True)
class SourcePair:
    old_gff: str
    new_gff: str
    old_fasta: str = ""
    new_fasta: str = ""


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, Optional[str]]]) -> None:
        if tag.lower() != "a":
            return
        attrs_dict = dict(attrs)
        href = attrs_dict.get("href")
        if href:
            self.links.append(href)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare real old/new GFF3 annotations and pick stable-ID truth examples. "
            "Use either --old-gff/--new-gff or Ensembl-style --ftp-root/--old-release/"
            "--new-release/--species."
        )
    )
    parser.add_argument("--old-gff", help="Old/reference GFF3 path or URL")
    parser.add_argument("--new-gff", help="New/target GFF3 path or URL")
    parser.add_argument(
        "--ftp-root",
        help="Ensembl-style FTP root, e.g. https://ftp.ensembl.org/pub or a local pub mirror",
    )
    parser.add_argument("--old-release", type=int)
    parser.add_argument("--new-release", type=int)
    parser.add_argument("--species", help="Species directory name, e.g. gallus_gallus")
    parser.add_argument(
        "--feature-type",
        choices=sorted(FEATURE_TYPES),
        default="gene",
        help="Feature class to compare",
    )
    parser.add_argument("--seqid", help="Restrict examples to one seqid/chromosome")
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="Maximum examples per category",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("real_stable_id_truth_examples.tsv"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Optional JSON manifest with resolved source URLs/paths",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    source_pair = resolve_source_pair(args)
    rows = pick_truth_examples(
        source_pair.old_gff,
        source_pair.new_gff,
        feature_type=args.feature_type,
        seqid=args.seqid,
        limit=args.limit,
    )
    write_rows(rows, args.output)
    if args.manifest:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(
                {
                    "old_gff": source_pair.old_gff,
                    "new_gff": source_pair.new_gff,
                    "old_fasta": source_pair.old_fasta,
                    "new_fasta": source_pair.new_fasta,
                    "feature_type": args.feature_type,
                    "seqid": args.seqid,
                    "output": str(args.output),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
    print(f"Wrote {len(rows)} examples: {args.output}")
    for category in sorted({row['category'] for row in rows}):
        count = sum(1 for row in rows if row["category"] == category)
        print(f"  {category}: {count}")


def resolve_source_pair(args: argparse.Namespace) -> SourcePair:
    if args.old_gff and args.new_gff:
        return SourcePair(old_gff=args.old_gff, new_gff=args.new_gff)
    required = {
        "--ftp-root": args.ftp_root,
        "--old-release": args.old_release,
        "--new-release": args.new_release,
        "--species": args.species,
    }
    missing = [name for name, value in required.items() if value in (None, "")]
    if missing:
        raise SystemExit(
            "Provide either --old-gff/--new-gff or all FTP options: "
            + ", ".join(missing)
        )
    return resolve_ensembl_ftp_pair(
        args.ftp_root,
        args.old_release,
        args.new_release,
        args.species,
    )


def resolve_ensembl_ftp_pair(
    ftp_root: str,
    old_release: int,
    new_release: int,
    species: str,
) -> SourcePair:
    old_base = _join_location(ftp_root, f"release-{old_release}")
    new_base = _join_location(ftp_root, f"release-{new_release}")
    return SourcePair(
        old_gff=_find_one(
            _join_location(old_base, "gff3", species),
            suffix=".gff3.gz",
            contains=f".{old_release}.",
        ),
        new_gff=_find_one(
            _join_location(new_base, "gff3", species),
            suffix=".gff3.gz",
            contains=f".{new_release}.",
        ),
        old_fasta=_find_one(
            _join_location(old_base, "fasta", species, "dna"),
            suffix=".dna.toplevel.fa.gz",
            required=False,
        ),
        new_fasta=_find_one(
            _join_location(new_base, "fasta", species, "dna"),
            suffix=".dna.toplevel.fa.gz",
            required=False,
        ),
    )


def pick_truth_examples(
    old_gff: str,
    new_gff: str,
    feature_type: str,
    seqid: Optional[str],
    limit: int,
) -> list[dict[str, str]]:
    old_features = parse_gff_features(old_gff, feature_type=feature_type, seqid=seqid)
    new_features = parse_gff_features(new_gff, feature_type=feature_type, seqid=seqid)
    shared_ids = sorted(set(old_features) & set(new_features))
    old_only_ids = sorted(set(old_features) - set(new_features))
    new_only_ids = sorted(set(new_features) - set(old_features))

    exact: list[tuple[RealFeature, RealFeature]] = []
    moved: list[tuple[RealFeature, RealFeature]] = []
    version_changed: list[tuple[RealFeature, RealFeature]] = []
    for stable_id in shared_ids:
        old = old_features[stable_id]
        new = new_features[stable_id]
        if _same_location(old, new):
            exact.append((old, new))
        else:
            moved.append((old, new))
        if old.version != new.version:
            version_changed.append((old, new))

    rows: list[dict[str, str]] = []
    rows.extend(_shared_rows("shared_same_locus", exact, limit))
    rows.extend(_shared_rows("shared_changed_locus", moved, limit))
    rows.extend(_shared_rows("shared_version_changed", version_changed, limit))
    rows.extend(_single_rows("old_only_missing_candidate", old_features, old_only_ids[:limit], "old"))
    rows.extend(_single_rows("new_only_new_candidate", new_features, new_only_ids[:limit], "new"))
    return rows


def parse_gff_features(
    source: str,
    feature_type: str,
    seqid: Optional[str] = None,
) -> dict[str, RealFeature]:
    allowed_types = FEATURE_TYPES[feature_type]
    features: dict[str, RealFeature] = {}
    with _open_text(source) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            row_seqid, _source, raw_type, start, end, _score, strand, _phase, attrs_text = fields
            if seqid is not None and row_seqid != seqid:
                continue
            if raw_type.lower() not in allowed_types:
                continue
            attrs = parse_attrs(attrs_text)
            stable_id, embedded_version = split_stable_id(attrs.get("ID"))
            if stable_id is None:
                continue
            features[stable_id] = RealFeature(
                stable_id=stable_id,
                version=parse_version(attrs, embedded_version),
                feature_type=feature_type,
                seqid=row_seqid,
                start=int(start),
                end=int(end),
                strand=strand,
                source=source,
            )
    return features


def write_rows(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "category",
        "feature_type",
        "stable_id",
        "old_stable_id",
        "target_stable_id",
        "old_version",
        "target_version",
        "old_location",
        "target_location",
        "note",
        "old_source",
        "new_source",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _shared_rows(
    category: str,
    pairs: list[tuple[RealFeature, RealFeature]],
    limit: int,
) -> list[dict[str, str]]:
    return [
        {
            "category": category,
            "feature_type": old.feature_type,
            "stable_id": old.stable_id,
            "old_stable_id": old.stable_id,
            "target_stable_id": new.stable_id,
            "old_version": str(old.version),
            "target_version": str(new.version),
            "old_location": old.location,
            "target_location": new.location,
            "note": (
                "same stable ID present in both annotations; expected positive mapping"
            ),
            "old_source": old.source,
            "new_source": new.source,
        }
        for old, new in pairs[:limit]
    ]


def _single_rows(
    category: str,
    features: dict[str, RealFeature],
    stable_ids: Iterable[str],
    side: str,
) -> list[dict[str, str]]:
    rows = []
    for stable_id in stable_ids:
        feature = features[stable_id]
        rows.append(
            {
                "category": category,
                "feature_type": feature.feature_type,
                "stable_id": feature.stable_id,
                "old_stable_id": feature.stable_id if side == "old" else "",
                "target_stable_id": feature.stable_id if side == "new" else "",
                "old_version": str(feature.version) if side == "old" else "",
                "target_version": str(feature.version) if side == "new" else "",
                "old_location": feature.location if side == "old" else "",
                "target_location": feature.location if side == "new" else "",
                "note": (
                    "stable ID found only in the old annotation"
                    if side == "old"
                    else "stable ID found only in the new annotation"
                ),
                "old_source": feature.source if side == "old" else "",
                "new_source": feature.source if side == "new" else "",
            }
        )
    return rows


def _same_location(old: RealFeature, new: RealFeature) -> bool:
    return (
        old.seqid == new.seqid
        and old.start == new.start
        and old.end == new.end
        and old.strand == new.strand
    )


def _find_one(
    directory: str,
    suffix: str,
    contains: str = "",
    required: bool = True,
) -> str:
    links = [
        link
        for link in _list_location(directory)
        if link.endswith(suffix)
        and (not contains or contains in link)
        and "abinitio" not in link.lower()
        and "chr." not in link.lower()
    ]
    if not links:
        if required:
            raise FileNotFoundError(f"No {suffix} file found in {directory}")
        return ""
    links.sort()
    return _join_location(directory, links[0])


def _list_location(location: str) -> list[str]:
    if _is_url(location):
        with urllib.request.urlopen(location.rstrip("/") + "/") as response:
            text = response.read().decode("utf-8", errors="replace")
        parser = LinkParser()
        parser.feed(text)
        return [
            link
            for link in parser.links
            if link not in ("../", "./") and not link.startswith("?")
        ]
    path = Path(location)
    return [child.name for child in path.iterdir()]


def _open_text(source: str):
    if _is_url(source):
        response = urllib.request.urlopen(source)
        data = response.read()
        if source.endswith(".gz"):
            return io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(data)))
        return io.StringIO(data.decode("utf-8"))
    path = Path(source)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open(encoding="utf-8")


def _join_location(*parts: str | int) -> str:
    values = [str(part) for part in parts if str(part)]
    if not values:
        return ""
    first = values[0].rstrip("/")
    rest = [part.strip("/") for part in values[1:]]
    if "://" in first:
        scheme, url_rest = first.split("://", 1)
        return scheme + "://" + "/".join([url_rest.strip("/")] + rest)
    return str(Path(first, *rest))


def _is_url(value: str) -> bool:
    return value.startswith(("http://", "https://", "ftp://"))


if __name__ == "__main__":
    main()

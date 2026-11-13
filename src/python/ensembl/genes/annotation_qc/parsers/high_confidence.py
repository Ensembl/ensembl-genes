# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Parse GTF/GFF3, FASTA and evidence sources into lossless normalized tables.
Coordinates follow the framework: Start is zero-based, End is exclusive."""

from __future__ import annotations

import csv
import json
import gzip
import hashlib
import math
import re
from pathlib import Path
from dataclasses import replace
from urllib.parse import unquote
import pandas as pd
from ensembl.genes.annotation_qc.config.high_confidence import (
    COMPLEMENT,
    HOMOLOGY_FIELDS,
)
from ensembl.genes.annotation_qc.parsers._high_confidence_records import (
    Feature,
    Transcript,
    Gene,
    SequenceRegion,
)


def open_text(path: Path):
    """Open plain or gzipped UTF-8 annotation text."""
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8")
    )


def parse_attributes(raw: str) -> dict[str, list[str]]:
    """Read quoted/unquoted GTF and key=value attributes, retaining repeated tags."""
    result: dict[str, list[str]] = {}
    pattern = re.compile(
        r'\s*([^\s;=]+)(?:\s*=\s*|\s+)(?:"((?:\\.|[^"\\])*)"|([^;\s]+))\s*(?:;|$)'
    )
    offset = 0
    while raw[offset:].strip():
        match = pattern.match(raw, offset)
        if not match:
            raise ValueError(f"Invalid GTF attributes near {raw[offset:]!r}")
        value = match[2] if match[2] is not None else match[3]
        result.setdefault(match[1], []).append(re.sub(r'\\(["\\])', r"\1", value))
        offset = match.end()
    return result


def merge_attrs(target: dict, source: dict) -> None:
    """Merge attribute lists while preserving distinct values."""
    for key, values in source.items():
        current = target.setdefault(key, [])
        current.extend(v for v in values if v not in current)


def single_attr(attrs: dict, *keys: str) -> str | None:
    """Resolve one attribute value and reject conflicting aliases."""
    values = {value for key in keys for value in attrs.get(key, [])}
    if len(values) > 1:
        raise ValueError(f"Conflicting attributes {keys}: {sorted(values)}")
    return next(iter(values), None)


def parse_gff3_attributes(raw: str) -> dict[str, list[str]]:
    """Split GFF3 lists before decoding escaped delimiters in identifiers."""
    attrs: dict[str, list[str]] = {}
    if raw.strip() == ".":
        return attrs
    for item in raw.strip().rstrip(";").split(";"):
        key, separator, value = item.partition("=")
        if not separator or not key or not value:
            raise ValueError(f"Invalid GFF3 attribute: {item!r}")
        attrs.setdefault(unquote(key), []).extend(unquote(v) for v in value.split(","))
    return attrs


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-branches,too-many-locals
def normalize_gff3(records, path):
    """Resolve gene -> transcript -> child links without depending on row order."""
    nodes = {}
    for feature, line_number in records:
        if feature.kind in {"gene", "transcript"}:
            identifier = single_attr(feature.attrs, "ID")
            if not identifier:
                raise ValueError(
                    f"{path}:{line_number}: GFF3 {feature.kind} requires ID"
                )
            if identifier in nodes:
                raise ValueError(
                    f"{path}:{line_number}: duplicate GFF3 parent ID {identifier}"
                )
            nodes[identifier] = feature

    def gene_id(identifier):
        parent = nodes.get(identifier)
        if parent is None or parent.kind != "gene":
            raise ValueError(f"Parent {identifier!r} does not reference a GFF3 gene")
        return single_attr(parent.attrs, "gene_id") or identifier

    def transcript_ids(identifier):
        parent = nodes.get(identifier)
        if parent is None or parent.kind != "transcript":
            raise ValueError(
                f"Parent {identifier!r} does not reference a GFF3 transcript/mRNA"
            )
        parents = parent.attrs.get("Parent", [])
        if len(parents) != 1:
            raise ValueError(
                f"Transcript {identifier!r} requires exactly one gene Parent"
            )
        resolved_gene = gene_id(parents[0])
        supplied_gene = single_attr(parent.attrs, "gene_id")
        if supplied_gene and supplied_gene != resolved_gene:
            raise ValueError(
                f"Conflicting gene_id and Parent for transcript {identifier!r}"
            )
        return resolved_gene, single_attr(parent.attrs, "transcript_id") or identifier

    # Keep the established reconstruction/scoring flow together for parity.
    # pylint: disable-next=too-many-nested-blocks
    for feature, line_number in records:
        try:
            attrs = feature.attrs
            if feature.kind == "gene":
                assignments = [(gene_id(single_attr(attrs, "ID")), None)]
            elif feature.kind == "transcript":
                assignments = [transcript_ids(single_attr(attrs, "ID"))]
            else:
                parents = attrs.get("Parent", [])
                if not parents:
                    raise ValueError(
                        f"GFF3 {feature.kind} requires a transcript Parent"
                    )
                assignments = [transcript_ids(parent) for parent in parents]
            for resolved_gene, resolved_tx in dict.fromkeys(assignments):
                updated = dict(attrs, gene_id=[resolved_gene])
                if resolved_tx is not None:
                    for key, value in [
                        ("gene_id", resolved_gene),
                        ("transcript_id", resolved_tx),
                    ]:
                        supplied = single_attr(attrs, key)
                        if supplied and supplied != value:
                            raise ValueError(f"Conflicting {key} and GFF3 Parent")
                    updated["transcript_id"] = [resolved_tx]
                yield replace(feature, attrs=updated), line_number
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: {error}") from error


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-branches,too-many-locals,too-many-statements
def _read_annotation_records(path: Path) -> dict[str, Gene]:
    genes: dict[str, Gene] = {}
    transcript_genes: dict[str, str] = {}
    deferred = []
    gff_records = []

    def add_feature(feature, line_number):
        attrs, kind = feature.attrs, feature.kind
        gene_id = single_attr(attrs, "gene_id")
        tx_id = single_attr(attrs, "transcript_id")
        if not gene_id:
            if tx_id:
                deferred.append((feature, line_number))
                return
            raise ValueError(f"Missing gene_id and transcript_id on {kind} feature")
        gene = genes.setdefault(
            gene_id,
            Gene(gene_id, feature.seqname, feature.strand, feature.start, feature.end),
        )
        if (gene.seqname, gene.strand) != (feature.seqname, feature.strand):
            raise ValueError(
                f"Gene {gene_id} occurs on multiple contigs or strands; "
                "IDs must be unique across loci"
            )
        gene.start, gene.end = min(gene.start, feature.start), max(
            gene.end, feature.end
        )
        merge_attrs(
            gene.attrs,
            {
                key: values
                for key, values in attrs.items()
                if key in {"gene_biotype", "gene_type"}
            },
        )
        if tx_id and transcript_genes.setdefault(tx_id, gene_id) != gene_id:
            raise ValueError(f"Transcript {tx_id} belongs to multiple genes")
        if kind == "gene":
            if gene.declared is not None:
                raise ValueError(f"Duplicate gene feature: {gene_id}")
            gene.declared = (feature.start, feature.end)
            return
        if not tx_id:
            raise ValueError(f"Missing transcript_id for {kind}")
        transcript = gene.transcripts.setdefault(
            tx_id, Transcript(tx_id, gene_id, feature.seqname, feature.strand)
        )
        merge_attrs(transcript.attrs, attrs)
        transcript.features.append(feature)

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.strip() == "##FASTA":
                break
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                columns = line.rstrip("\r\n").split("\t")
                if len(columns) != 9:
                    columns = line.strip().split(None, 8)
                if len(columns) != 9:
                    raise ValueError("Expected nine GTF columns")
                seqname, _, kind, start, end, _, strand, phase, raw_attrs = columns
                kind = {
                    "cds": "CDS",
                    "mrna": "transcript",
                    "selenocysteine": "Selenocysteine",
                    "utr": "UTR",
                }.get(kind.lower(), kind.lower())
                start, end = int(start), int(end)
                if start < 1 or end < start or strand not in {"+", "-"}:
                    raise ValueError(
                        "Invalid interval or strand; require 1 <= start <= end and +/-"
                    )
                if phase not in {".", "0", "1", "2"}:
                    raise ValueError("Phase must be ., 0, 1 or 2")
                bare = raw_attrs.strip().rstrip(";").strip()
                is_gff3 = bool(re.match(r"[^\s;=]+=", bare))
                if is_gff3:
                    attrs = parse_gff3_attributes(raw_attrs)
                elif (
                    bare != "."
                    and kind in {"gene", "transcript"}
                    and re.fullmatch(r'[^\s;="]+', bare)
                ):
                    attrs = {"gene_id" if kind == "gene" else "transcript_id": [bare]}
                else:
                    attrs = parse_attributes(raw_attrs)
                feature = Feature(
                    seqname,
                    kind,
                    start,
                    end,
                    strand,
                    None if phase == "." else int(phase),
                    attrs,
                )
                if is_gff3 and ("ID" in attrs or "Parent" in attrs):
                    gff_records.append((feature, line_number))
                else:
                    add_feature(feature, line_number)
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: {error}") from error
    for feature, line_number in normalize_gff3(gff_records, path):
        try:
            add_feature(feature, line_number)
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: {error}") from error
    # Tiberius/AUGUSTUS-style parent rows may contain only a transcript ID.
    # Resolve its gene through explicit child attributes, independent of row order.
    for feature, line_number in deferred:
        tx_id = single_attr(feature.attrs, "transcript_id")
        assert tx_id is not None
        gene_id = transcript_genes.get(tx_id)
        if gene_id is None:
            raise ValueError(
                f"{path}:{line_number}: cannot link transcript {tx_id} to a gene; "
                "include gene_id on a transcript or child feature"
            )
        add_feature(
            replace(feature, attrs={**feature.attrs, "gene_id": [gene_id]}), line_number
        )
        tx = genes[gene_id].transcripts[tx_id]
        note = {
            "action": "gene_link_from_transcript",
            "details": "Resolved missing gene_id using explicit attributes on related features",
        }
        if note not in tx.normalization:
            tx.normalization.append(note)
    if not genes:
        raise ValueError("No annotation features found")
    for gene in genes.values():
        if gene.declared and gene.declared != (gene.start, gene.end):
            raise ValueError(f"Features outside declared bounds of gene {gene.id}")
        for tx in gene.transcripts.values():
            if len(gene.transcripts) == 1:
                tx.parent_bounds = (gene.start, gene.end)
    return genes


class Genome:
    """Read regularly wrapped plain FASTA using an in-memory byte-offset index."""

    # Keep the established reconstruction/scoring flow together for parity.
    # pylint: disable-next=too-many-branches
    def __init__(self, path: Path):
        if path.suffix == ".gz":
            raise ValueError("Genome must be plain FASTA; decompress it before running")
        self.handle = path.open("rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        try:
            name = None
            length = offset = width = bases = previous_bases = previous_width = 0
            while True:
                position = self.handle.tell()
                line = self.handle.readline()
                if not line or line.startswith(b">"):
                    if name is not None:
                        if not length:
                            raise ValueError(f"Empty FASTA contig: {name}")
                        self.index[name] = (length, offset, bases, width)
                    if not line:
                        break
                    tokens = line[1:].split()
                    if not tokens:
                        raise ValueError("Empty FASTA header")
                    name = tokens[0].decode("utf-8")
                    if name in self.index:
                        raise ValueError(f"Duplicate FASTA ID: {name}")
                    length, offset, bases, width = 0, self.handle.tell(), 0, 0
                    continue
                sequence = line.rstrip(b"\r\n")
                if name is None or not re.fullmatch(
                    b"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence
                ):
                    raise ValueError(f"Invalid FASTA sequence at byte {position}")
                if length == 0:
                    bases, width = len(sequence), len(line)
                elif (
                    previous_bases != bases
                    or previous_width != width
                    or len(sequence) > bases
                ):
                    raise ValueError(
                        f"FASTA lines must have regular wrapping for {name}"
                    )
                length += len(sequence)
                previous_bases, previous_width = len(sequence), len(line)
            if not self.index:
                raise ValueError("Empty genome FASTA")
        except Exception:
            self.handle.close()
            raise

    def close(self) -> None:
        """Close the indexed FASTA file handle."""
        self.handle.close()

    def contains(self, seqname: str, start: int, end: int) -> bool:
        """Check inclusive bounds against the indexed contig length."""
        entry = self.index.get(seqname)
        return entry is not None and 1 <= start <= end <= entry[0]

    def fetch(self, seqname: str, start: int, end: int, strand: str = "+") -> str:
        """Return an inclusive interval in transcriptional orientation."""
        if seqname not in self.index:
            raise ValueError(
                f"GTF contig {seqname!r} absent from genome FASTA; IDs must match"
            )
        length, offset, bases, width = self.index[seqname]
        if not 1 <= start <= end <= length:
            raise ValueError(
                f"Interval outside FASTA contig: {seqname}:{start}-{end} (length {length})"
            )
        first, last = start - 1, end - 1
        byte_start = offset + first // bases * width + first % bases
        byte_end = offset + last // bases * width + last % bases
        self.handle.seek(byte_start)
        sequence = (
            self.handle.read(byte_end - byte_start + 1)
            .replace(b"\n", b"")
            .replace(b"\r", b"")
            .decode()
            .upper()
        )
        return reverse_complement(sequence) if strand == "-" else sequence


def read_table(path: Path, required: list[str]):
    """Read tab-separated records after checking required columns."""
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"{path}: required TSV columns: {', '.join(required)}")
        for number, row in enumerate(reader, 2):
            if None in row or any(row.get(key) in {None, ""} for key in required):
                raise ValueError(f"{path}:{number}: incomplete or extra TSV fields")
            yield row


def read_junctions(path: Path | None):
    """Read unique inclusive junction intervals and read counts."""
    if path is None:
        return None
    result = {}
    for row in read_table(path, ["seqname", "start", "end", "strand", "read_count"]):
        start, end, count = int(row["start"]), int(row["end"]), int(row["read_count"])
        if start < 1 or end < start or count < 0 or row["strand"] not in {"+", "-"}:
            raise ValueError(f"Invalid junction: {row}")
        key = row["seqname"], start, end, row["strand"]
        if key in result:
            raise ValueError(f"Duplicate junction {key}; aggregate counts upstream")
        result[key] = count
    return result


def read_homology(path: Path | None, selected_ids: set[str]) -> dict:
    """Read seven-column legacy DIAMOND helper results."""
    result: dict[str, dict] = {}
    if path is None:
        return result
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 7:
                raise ValueError(
                    f"{path}:{number}: expected seven helper result columns, without header"
                )
            query, hit, *values = fields
            if not hit or query not in selected_ids or query in result:
                raise ValueError(
                    f"{path}:{number}: unknown/duplicate query or empty hit ID: {query}"
                )
            numeric = [float(value) for value in values]
            if not all(math.isfinite(value) for value in numeric):
                raise ValueError(f"{path}:{number}: homology metrics must be finite")
            result[query] = dict(zip(HOMOLOGY_FIELDS, [hit] + numeric))
    return result


def file_sha256(path: Path) -> str:
    """Fingerprint a file using streaming SHA256."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement(sequence):
    """Normalize genomic sequence orientation without interpreting biological quality."""
    return sequence.upper().translate(COMPLEMENT)[::-1]


def parse_high_confidence_annotation(path: Path, input_mode="auto") -> pd.DataFrame:
    """Return standardized feature rows, preserving repeated attributes and hierarchy metadata.

    The existing PyRanges parser does not cover predictor bare-ID attributes,
    order-independent multi-parent normalization, or repeated tags required here.
    This parser shares the framework's zero-based, half-open table boundary.
    """
    if input_mode not in {"auto", "strict"}:
        raise ValueError(f"Invalid input mode: {input_mode}")
    genes = _read_annotation_records(Path(path))
    records = []
    gene_metadata = {}
    transcript_metadata = {}
    for gene in genes.values():
        gene_metadata[gene.id] = {
            "seqname": gene.seqname,
            "strand": gene.strand,
            "start": gene.start,
            "end": gene.end,
            "attrs": gene.attrs,
            "declared": gene.declared,
            "normalization": gene.normalization,
        }
        records.append(
            {
                "Chromosome": gene.seqname,
                "Feature": "gene",
                "Start": gene.start - 1,
                "End": gene.end,
                "Strand": gene.strand,
                "Frame": None,
                "gene_id": gene.id,
                "transcript_id": None,
                "Attributes": gene.attrs,
            }
        )
        for tx in gene.transcripts.values():
            transcript_metadata[(gene.id, tx.id)] = {
                "attrs": tx.attrs,
                "parent_bounds": tx.parent_bounds,
                "normalization": tx.normalization,
            }
            for feature in tx.features:
                records.append(
                    {
                        "Chromosome": feature.seqname,
                        "Feature": feature.kind,
                        "Start": feature.start - 1,
                        "End": feature.end,
                        "Strand": feature.strand,
                        "Frame": feature.phase,
                        "gene_id": gene.id,
                        "transcript_id": tx.id,
                        "Attributes": feature.attrs,
                    }
                )
    frame = pd.DataFrame(
        records,
        columns=[
            "Chromosome",
            "Feature",
            "Start",
            "End",
            "Strand",
            "Frame",
            "gene_id",
            "transcript_id",
            "Attributes",
        ],
        dtype=object,
    )
    frame.attrs["gene_metadata"] = gene_metadata
    frame.attrs["transcript_metadata"] = transcript_metadata
    return frame


def parse_canonical_assignments(path):
    """Read explicit representative assignments without applying selection policy."""
    rows = list(read_table(Path(path), ["gene_id", "transcript_id"])) if path else []
    return pd.DataFrame(rows, columns=["gene_id", "transcript_id"], dtype=object)


def parse_legacy_diamond_alignments(path):
    """Read the helper's nine-column format (different from the generic 13-column parser)."""
    rows = []
    with open_text(Path(path)) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            values = line.split()
            if len(values) != 9:
                raise ValueError(f"{path}:{number}: expected nine DIAMOND fields")
            rows.append(values)
    return pd.DataFrame(
        rows,
        columns=[
            "qseqid",
            "sseqid",
            "qstart",
            "qend",
            "slen",
            "pident",
            "nident",
            "length",
            "gaps",
        ],
        dtype=object,
    )


def read_transcript_region(genome, tx):
    """Perform sequence I/O before metric evaluation; cache only the current transcript region."""
    length = genome.index[tx.seqname][0]
    first = max(1, min(f.start for f in tx.features) - 3)
    last = min(length, max(f.end for f in tx.features) + 3)
    return SequenceRegion(
        tx.seqname, first, genome.fetch(tx.seqname, first, last), length
    )


def parse_confidence_evidence(path):
    """Read stored evidence as a DataFrame without changing unknown/false values."""
    required = [
        "gene_id",
        "transcript_id",
        "confidence",
        "confidence_rule",
        "homology_status",
        "error",
        "score",
    ]
    rows = list(read_table(Path(path), required))
    with open_text(Path(path)) as handle:
        fields = next(csv.reader(handle, delimiter="\t"))
    return pd.DataFrame(rows, columns=fields, dtype=object)


def parse_json_manifest(path):
    """Read run metadata without inferring scoring decisions."""

    return json.loads(Path(path).read_text(encoding="utf-8"))


def parse_generated_gff(path):
    """Read an exported GFF without altering rows; expose gene IDs for validation."""
    records = []
    for line in Path(path).read_text(encoding="utf-8").splitlines(keepends=True):
        columns = line.rstrip("\r\n").split("\t")
        gene_id = None
        if len(columns) == 9 and columns[2] == "gene":
            attrs = parse_gff3_attributes(columns[8])
            gene_id = single_attr(attrs, "gene_id")
            if gene_id is None:
                gene_id = single_attr(attrs, "ID")
                if gene_id and gene_id.startswith("gene:"):
                    gene_id = gene_id[5:]
        records.append({"raw_line": line, "columns": columns, "gene_id": gene_id})
    return pd.DataFrame(
        records, columns=["raw_line", "columns", "gene_id"], dtype=object
    )


def parse_tsv_frame(path):
    """Read a headered TSV with exact row widths and no automatic NA coercion."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"{path}: missing or duplicate header names")
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"{path}: rows do not match the header")
    return pd.DataFrame(rows, columns=fields, dtype=object)


def parse_confidence_run(path, manifest_path=None, introns_path=None):
    """Read evidence and optional companions, keeping provenance separate from metrics."""
    path = Path(path).resolve()
    if manifest_path is not None and not Path(manifest_path).is_file():
        raise ValueError(f"Manifest not found: {manifest_path}")
    if introns_path is not None and not Path(introns_path).is_file():
        raise ValueError(f"Intron table not found: {introns_path}")
    if path.name.endswith(".evidence.tsv"):
        prefix = path.name[: -len(".evidence.tsv")]
        manifest_path = (
            Path(manifest_path)
            if manifest_path
            else path.with_name(prefix + ".manifest.json")
        )
        introns_path = (
            Path(introns_path)
            if introns_path
            else path.with_name(prefix + ".introns.tsv")
        )
    manifest = (
        parse_json_manifest(manifest_path)
        if manifest_path and Path(manifest_path).is_file()
        else None
    )
    introns = (
        parse_tsv_frame(introns_path)
        if introns_path and Path(introns_path).is_file()
        else None
    )
    return {
        "evidence": parse_tsv_frame(path),
        "manifest": manifest,
        "introns": introns,
        "evidence_sha256": file_sha256(path),
        "source": {
            "evidence": str(path),
            "manifest": str(manifest_path) if manifest else None,
        },
    }


def parse_confidence_models(path):
    """Normalize model/evidence TSV paths relative to the configuration file."""
    path = Path(path).resolve()
    frame = parse_tsv_frame(path)
    if not {"model", "evidence"} <= set(frame.columns) or frame.empty:
        raise ValueError(
            "Models TSV requires model/evidence columns and at least one row"
        )
    if frame["model"].duplicated().any() or any(
        not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value)
        for value in frame["model"]
    ):
        raise ValueError(
            "Models must be unique safe labels starting with a letter/digit"
        )
    rows = []
    for record in frame.to_dict("records"):
        for key in ["evidence", "manifest", "introns"]:
            value = record.get(key, "")
            record[key] = (
                (path.parent / value).resolve() if value and value != "." else None
            )
        if record["evidence"] is None:
            raise ValueError(f'Missing evidence for {record["model"]}')
        rows.append(record)
    return pd.DataFrame(rows, dtype=object)

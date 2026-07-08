#!/usr/bin/env python3
"""Diagnose common LiftOn input failures before rerunning a full projection."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from stable_id_mapping.gff3 import parse_attrs

TRANSCRIPT_TYPES = {"mrna", "transcript"}
CHILD_TYPES = {"exon", "cds"}


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    children: list[tuple[str, int, int]] = field(default_factory=list)


@dataclass(frozen=True)
class Issue:
    severity: str
    issue_type: str
    feature_id: str
    detail: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check GFF3/FASTA inputs for likely LiftOn transcript extraction issues."
    )
    parser.add_argument("--ref-fasta", type=Path, required=True)
    parser.add_argument("--ref-gff", type=Path, required=True)
    parser.add_argument(
        "--lifton-transcripts-fa",
        type=Path,
        help="Optional failed LiftOn intermediate transcripts.fa to validate",
    )
    parser.add_argument(
        "--output-tsv",
        type=Path,
        help="Optional TSV with detailed issues",
    )
    parser.add_argument("--example-limit", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fasta_lengths = read_fasta_lengths(args.ref_fasta)
    transcripts, feature_counts, duplicate_ids = read_gff_transcripts(args.ref_gff)
    issues = diagnose_models(transcripts, fasta_lengths)
    if duplicate_ids:
        for feature_id, count in duplicate_ids.most_common(args.example_limit):
            issues.append(
                Issue(
                    severity="warning",
                    issue_type="duplicate_feature_id",
                    feature_id=feature_id,
                    detail=f"{count} rows share this ID",
                )
            )

    if args.lifton_transcripts_fa is not None and args.lifton_transcripts_fa.exists():
        issues.extend(validate_fasta_records(args.lifton_transcripts_fa))

    print_summary(
        fasta_lengths=fasta_lengths,
        transcripts=transcripts,
        feature_counts=feature_counts,
        issues=issues,
        example_limit=args.example_limit,
    )
    if args.output_tsv is not None:
        try:
            write_issues(issues, args.output_tsv)
        except OSError as error:
            print(f"Warning: could not write diagnostics TSV {args.output_tsv}: {error}")

    fatal_count = sum(1 for issue in issues if issue.severity == "fatal")
    if fatal_count:
        raise SystemExit(2)


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: Optional[str] = None
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths.setdefault(current, 0)
                continue
            if current is None:
                continue
            lengths[current] += len(line)
    return lengths


def read_gff_transcripts(
    path: Path,
) -> tuple[dict[str, TranscriptModel], Counter[str], Counter[str]]:
    transcripts: dict[str, TranscriptModel] = {}
    child_rows: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    feature_counts: Counter[str] = Counter()
    seen_ids: Counter[str] = Counter()

    with open_text(path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature_type, start, end, _score, strand, _phase, attrs_text = fields
            feature_type_lc = feature_type.lower()
            feature_counts[feature_type_lc] += 1
            attrs = parse_attrs(attrs_text)
            feature_id = attrs.get("ID")
            if feature_id:
                seen_ids[feature_id] += 1
            start_i = int(start)
            end_i = int(end)

            if feature_type_lc in TRANSCRIPT_TYPES:
                transcript_id = attrs.get("ID") or attrs.get("transcript_id") or attrs.get("Name")
                if not transcript_id:
                    continue
                transcripts[transcript_id] = TranscriptModel(
                    transcript_id=transcript_id,
                    gene_id=(attrs.get("Parent") or attrs.get("gene_id") or "").split(",")[0],
                    seqid=seqid,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                )
                continue

            if feature_type_lc in CHILD_TYPES:
                parent = attrs.get("Parent")
                if not parent:
                    continue
                for parent_id in parent.split(","):
                    child_rows[parent_id].append((feature_type_lc, start_i, end_i))

    for transcript_id, model in transcripts.items():
        model.children.extend(child_rows.get(transcript_id, []))

    duplicate_ids = Counter(
        {feature_id: count for feature_id, count in seen_ids.items() if count > 1}
    )
    return transcripts, feature_counts, duplicate_ids


def diagnose_models(
    transcripts: dict[str, TranscriptModel],
    fasta_lengths: dict[str, int],
) -> list[Issue]:
    issues: list[Issue] = []
    for transcript in transcripts.values():
        seq_len = fasta_lengths.get(transcript.seqid)
        if seq_len is None:
            issues.append(
                Issue(
                    severity="fatal",
                    issue_type="seqid_missing_from_fasta",
                    feature_id=transcript.transcript_id,
                    detail=f"{transcript.seqid} is present in GFF3 but not FASTA",
                )
            )
            continue
        if transcript.start < 1 or transcript.end < transcript.start:
            issues.append(
                Issue(
                    severity="fatal",
                    issue_type="invalid_transcript_coordinates",
                    feature_id=transcript.transcript_id,
                    detail=f"{transcript.seqid}:{transcript.start}-{transcript.end}",
                )
            )
        if transcript.end > seq_len:
            issues.append(
                Issue(
                    severity="fatal",
                    issue_type="transcript_out_of_fasta_bounds",
                    feature_id=transcript.transcript_id,
                    detail=f"{transcript.seqid}:{transcript.start}-{transcript.end} > length {seq_len}",
                )
            )
        if not transcript.children:
            issues.append(
                Issue(
                    severity="warning",
                    issue_type="transcript_without_exon_or_cds_children",
                    feature_id=transcript.transcript_id,
                    detail=f"{transcript.seqid}:{transcript.start}-{transcript.end}",
                )
            )
        for child_type, start, end in transcript.children:
            if start < 1 or end < start:
                issues.append(
                    Issue(
                        severity="fatal",
                        issue_type=f"invalid_{child_type}_coordinates",
                        feature_id=transcript.transcript_id,
                        detail=f"{transcript.seqid}:{start}-{end}",
                    )
                )
            elif end > seq_len:
                issues.append(
                    Issue(
                        severity="fatal",
                        issue_type=f"{child_type}_out_of_fasta_bounds",
                        feature_id=transcript.transcript_id,
                        detail=f"{transcript.seqid}:{start}-{end} > length {seq_len}",
                    )
                )
    return issues


def validate_fasta_records(path: Path) -> list[Issue]:
    issues: list[Issue] = []
    current_header: Optional[str] = None
    current_len = 0
    seen_any_header = False
    line_no = 0
    with open_text(path) as handle:
        for raw_line in handle:
            line_no += 1
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None and current_len == 0:
                    issues.append(
                        Issue(
                            severity="fatal",
                            issue_type="empty_fasta_record",
                            feature_id=current_header,
                            detail=f"{path}:{line_no}",
                        )
                    )
                current_header = line[1:].split()[0] or "<empty header>"
                current_len = 0
                seen_any_header = True
                continue
            if not line.strip():
                continue
            if current_header is None:
                issues.append(
                    Issue(
                        severity="fatal",
                        issue_type="sequence_before_first_fasta_header",
                        feature_id=str(line_no),
                        detail=line[:80],
                    )
                )
                continue
            current_len += len(line.strip())
    if current_header is not None and current_len == 0:
        issues.append(
            Issue(
                severity="fatal",
                issue_type="empty_fasta_record",
                feature_id=current_header,
                detail=f"{path}:EOF",
            )
        )
    if not seen_any_header:
        issues.append(
            Issue(
                severity="fatal",
                issue_type="fasta_has_no_records",
                feature_id=str(path),
                detail="No FASTA headers found",
            )
        )
    return issues


def print_summary(
    fasta_lengths: dict[str, int],
    transcripts: dict[str, TranscriptModel],
    feature_counts: Counter[str],
    issues: list[Issue],
    example_limit: int,
) -> None:
    print(f"FASTA seqids: {len(fasta_lengths)}")
    print(f"FASTA total bp: {sum(fasta_lengths.values())}")
    print(f"GFF3 transcripts: {len(transcripts)}")
    print(
        "GFF3 feature counts: "
        + ", ".join(f"{key}={value}" for key, value in feature_counts.most_common(10))
    )
    issue_counts = Counter(issue.issue_type for issue in issues)
    severity_counts = Counter(issue.severity for issue in issues)
    if severity_counts:
        severity_text = ", ".join(
            f"{key}={value}" for key, value in sorted(severity_counts.items())
        )
    else:
        severity_text = "none"
    print(f"Issue severity counts: {severity_text}")
    if not issue_counts:
        print("No obvious LiftOn input issues found.")
        return
    print("Issue type counts:")
    for issue_type, count in issue_counts.most_common():
        print(f"  {issue_type}: {count}")
    print("Examples:")
    severity_order = {"fatal": 0, "warning": 1}
    ordered_issues = sorted(
        enumerate(issues),
        key=lambda indexed_issue: (
            severity_order.get(indexed_issue[1].severity, 2),
            indexed_issue[0],
        ),
    )
    for _index, issue in ordered_issues[:example_limit]:
        print(
            f"  [{issue.severity}] {issue.issue_type} "
            f"{issue.feature_id}: {issue.detail}"
        )


def write_issues(issues: list[Issue], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["severity", "issue_type", "feature_id", "detail"])
        for issue in issues:
            writer.writerow(
                [issue.severity, issue.issue_type, issue.feature_id, issue.detail]
            )


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open(encoding="utf-8")


if __name__ == "__main__":
    main()

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from diagnose_lifton_inputs import (
    diagnose_models,
    read_fasta_lengths,
    read_gff_transcripts,
    validate_fasta_records,
)


def write_text(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def test_diagnose_lifton_inputs_flags_missing_seqid_and_bad_children(
    tmp_path: Path,
) -> None:
    fasta = tmp_path / "ref.fa"
    gff = tmp_path / "ref.gff3"
    write_text(
        fasta,
        """
>chr1
AAAAAA
        """,
    )
    write_text(
        gff,
        """
##gff-version 3
chr1	test	gene	1	10	.	+	.	ID=gene:ONE
chr1	test	mRNA	1	10	.	+	.	ID=transcript:ONE;Parent=gene:ONE
chr1	test	exon	1	12	.	+	.	ID=exon:ONE;Parent=transcript:ONE
chrMissing	test	mRNA	1	5	.	+	.	ID=transcript:MISSING;Parent=gene:MISSING
        """,
    )

    lengths = read_fasta_lengths(fasta)
    transcripts, _feature_counts, _duplicates = read_gff_transcripts(gff)
    issues = diagnose_models(transcripts, lengths)
    issue_types = {issue.issue_type for issue in issues}

    assert "transcript_out_of_fasta_bounds" in issue_types
    assert "exon_out_of_fasta_bounds" in issue_types
    assert "seqid_missing_from_fasta" in issue_types


def test_validate_lifton_transcripts_fasta_flags_empty_record(tmp_path: Path) -> None:
    transcripts_fa = tmp_path / "transcripts.fa"
    write_text(
        transcripts_fa,
        """
>good
ACTG
>empty
>next
ACTG
        """,
    )

    issues = validate_fasta_records(transcripts_fa)

    assert any(issue.issue_type == "empty_fasta_record" for issue in issues)

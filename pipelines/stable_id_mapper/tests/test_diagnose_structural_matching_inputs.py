from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from diagnose_structural_matching_inputs import summarize_gff3_structure


def write_text(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def test_summarize_gff3_structure_counts_exact_transcript_children(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "annotation.gff3"
    write_text(
        gff,
        """
##gff-version 3
chr1	test	gene	1	100	.	+	.	ID=gene:GENE1
chr1	test	mRNA	1	100	.	+	.	ID=transcript:TX1;Parent=gene:GENE1
chr1	test	exon	1	50	.	+	.	ID=exon:EX1;Parent=transcript:TX1
chr1	test	CDS	10	40	.	+	0	ID=cds:CDS1;Parent=transcript:TX1
        """,
    )

    summary = summarize_gff3_structure(gff)

    assert summary["transcript_count"] == 1
    assert summary["transcripts_with_exact_children"] == 1
    assert summary["exact_child_links"]["exon"] == 1
    assert summary["exact_child_links"]["cds"] == 1


def test_summarize_gff3_structure_reports_gene_parented_children(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "annotation.gff3"
    write_text(
        gff,
        """
##gff-version 3
chr1	test	gene	1	100	.	+	.	ID=gene:GENE1
chr1	test	mRNA	1	100	.	+	.	ID=transcript:TX1;Parent=gene:GENE1
chr1	test	exon	1	50	.	+	.	ID=exon:EX1;Parent=gene:GENE1
chr1	test	CDS	10	40	.	+	0	ID=cds:CDS1;Parent=gene:GENE1
        """,
    )

    summary = summarize_gff3_structure(gff)

    assert summary["transcript_count"] == 1
    assert summary["transcripts_with_exact_children"] == 0
    assert summary["transcripts_with_core_children"] == 0
    assert summary["unmatched_child_parent_links"]["exon"] == 1
    assert summary["unmatched_child_parent_links"]["cds"] == 1
    assert summary["parent_prefix_counts"]["gene"] == 2

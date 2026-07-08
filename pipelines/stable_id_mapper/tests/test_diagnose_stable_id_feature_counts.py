from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from diagnose_stable_id_feature_counts import summarize_feature_counts


def write_text(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def test_summarize_feature_counts_uses_decision_parser_gene_universe(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "annotation.gff3"
    write_text(
        gff,
        """
##gff-version 3
chr1	test	gene	1	100	.	+	.	ID=gene:ENSG0001.2;biotype=protein_coding
chr1	test	mRNA	1	100	.	+	.	ID=transcript:ENST0001.3;Parent=gene:ENSG0001.2
chr1	test	CDS	10	50	.	+	0	ID=CDS:ENSP0001;Parent=transcript:ENST0001.3;protein_id=ENSP0001.4
chr1	test	gene	200	300	.	+	.	ID=gene:ENSG0002.1;biotype=lncRNA
chr1	test	mRNA	200	300	.	+	.	ID=transcript:ENST0002.1;Parent=gene:ENSG0002.1
chr1	test	gene	400	500	.	+	.	ID=gene:CHILD;Parent=gene:ENSG0002.1;biotype=pseudogene
        """,
    )

    summary = summarize_feature_counts(gff)

    assert summary["parsed_counts"] == {
        "gene": 2,
        "transcript": 2,
        "translation": 1,
    }
    assert summary["parentless_gene_rows"] == 2
    assert summary["parented_gene_rows"] == 1
    assert summary["gene_biotypes"]["protein_coding"] == 1
    assert summary["gene_biotypes"]["lncRNA"] == 1
    assert summary["gene_biotypes"]["pseudogene"] == 1

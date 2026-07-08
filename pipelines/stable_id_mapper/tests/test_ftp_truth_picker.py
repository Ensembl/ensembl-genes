from __future__ import annotations

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pick_ftp_stable_id_truth_examples import (
    pick_truth_examples,
    resolve_ensembl_ftp_pair,
    write_rows,
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def write_gzip(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text.strip() + "\n")


OLD_GFF = """
##gff-version 3
chr1	test	gene	100	200	.	+	.	ID=gene:ENSREALG0001.1
chr1	test	gene	300	400	.	+	.	ID=gene:ENSREALG0002.1
chr1	test	gene	500	600	.	+	.	ID=gene:ENSREALG0003.1
chr1	test	mRNA	100	200	.	+	.	ID=transcript:ENSREALT0001.1;Parent=gene:ENSREALG0001.1
"""

NEW_GFF = """
##gff-version 3
chr1	test	gene	100	200	.	+	.	ID=gene:ENSREALG0001.2
chr1	test	gene	330	430	.	+	.	ID=gene:ENSREALG0002.2
chr1	test	gene	700	800	.	+	.	ID=gene:ENSREALG9001.1
chr1	test	mRNA	100	200	.	+	.	ID=transcript:ENSREALT0001.2;Parent=gene:ENSREALG0001.2
"""


def test_pick_truth_examples_from_old_and_new_gff3(tmp_path: Path) -> None:
    old_gff = tmp_path / "old.gff3"
    new_gff = tmp_path / "new.gff3"
    write_text(old_gff, OLD_GFF)
    write_text(new_gff, NEW_GFF)

    rows = pick_truth_examples(
        str(old_gff),
        str(new_gff),
        feature_type="gene",
        seqid=None,
        limit=10,
    )
    output = tmp_path / "truth.tsv"
    write_rows(rows, output)
    text = output.read_text(encoding="utf-8")

    assert "shared_same_locus" in text
    assert "shared_changed_locus" in text
    assert "shared_version_changed" in text
    assert "old_only_missing_candidate" in text
    assert "new_only_new_candidate" in text
    assert "ENSREALG0001" in text
    assert "ENSREALG0003" in text
    assert "ENSREALG9001" in text


def test_resolve_ensembl_style_local_ftp_pair(tmp_path: Path) -> None:
    old_gff = (
        tmp_path
        / "release-1"
        / "gff3"
        / "synthetic_species"
        / "Synthetic_species.Assembly.1.gff3.gz"
    )
    new_gff = (
        tmp_path
        / "release-2"
        / "gff3"
        / "synthetic_species"
        / "Synthetic_species.Assembly.2.gff3.gz"
    )
    old_fasta = (
        tmp_path
        / "release-1"
        / "fasta"
        / "synthetic_species"
        / "dna"
        / "Synthetic_species.Assembly.dna.toplevel.fa.gz"
    )
    new_fasta = (
        tmp_path
        / "release-2"
        / "fasta"
        / "synthetic_species"
        / "dna"
        / "Synthetic_species.Assembly.dna.toplevel.fa.gz"
    )
    write_gzip(old_gff, OLD_GFF)
    write_gzip(new_gff, NEW_GFF)
    write_gzip(old_fasta, ">chr1\nAAAA")
    write_gzip(new_fasta, ">chr1\nAAAA")

    pair = resolve_ensembl_ftp_pair(
        str(tmp_path),
        old_release=1,
        new_release=2,
        species="synthetic_species",
    )

    assert pair.old_gff == str(old_gff)
    assert pair.new_gff == str(new_gff)
    assert pair.old_fasta == str(old_fasta)
    assert pair.new_fasta == str(new_fasta)

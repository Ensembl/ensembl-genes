from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_stable_id_run import audit_run


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def test_audit_run_writes_missing_coordinate_and_new_gene_tables(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    output_dir = tmp_path / "audit"

    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	100	500	.	+	.	ID=gene:ENSG0001.1;biotype=protein_coding
chr1	test	gene	900	1200	.	+	.	ID=gene:ENSG0002.1;biotype=protein_coding
chr1	test	gene	1500	1800	.	-	.	ID=gene:ENSG0003.1;biotype=lncRNA
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	100	500	.	+	.	ID=gene:ENSG9001.1;biotype=protein_coding
chrT	test	gene	900	1200	.	+	.	ID=gene:ENSG9002.1;biotype=protein_coding
chrT	test	gene	2000	2300	.	+	.	ID=gene:ENSG0003.2;biotype=lncRNA
        """,
    )
    write_text(
        run_dir / "stable_id_decisions.tsv",
        """
type	action	current_stable_id	current_version	old_stable_id	old_version	new_stable_id	new_version	mapping_session_id	score	reason
gene	mapped	ENSG9001	1	ENSG0001	1	ENSG0001	2	1	0.91	mapped gene matched a target gene by lifton structural evidence; score_source=lifton_gene_aggregate confidence=high
gene	mapped	ENSG9002	1	ENSG0002	1	ENSG0002	2	1	0.42	mapped gene matched a target gene by coordinate overlap; score_source=coordinate_overlap
gene	missing		0	ENSG0003	1		0	1	0	mapped gene did not match any target gene by lifton evidence or coordinate overlap
gene	new	ENSG0003	2		0	ENSNEWG1	1	1	0	target gene was not claimed by any mapped gene
        """,
    )
    write_tsv(
        run_dir / "matching" / "lifton.gene_locus_comparison.tsv",
        [
            "lifton_gene",
            "target_gene_by_locus",
            "locus_score",
            "overlap_bp",
            "old_gene_coverage",
            "target_gene_coverage",
            "reciprocal_overlap",
            "span_containment",
            "length_ratio",
            "center_distance",
            "lifton_contig",
            "lifton_start",
            "lifton_end",
            "lifton_strand",
            "target_contig",
            "target_start",
            "target_end",
            "target_strand",
            "structure_accepted_target_gene",
            "structure_weighted_score",
            "structure_fraction_of_total",
            "structure_n_transcripts",
            "best_tx_structure_target_gene",
            "best_tx_structure_target_tx",
            "best_tx_structure_score",
            "locus_vs_structure",
        ],
        [
            [
                "gene:ENSG0002.1",
                "gene:ENSG9002.1",
                "0.420000",
                "126",
                "0.420000",
                "0.420000",
                "0.420000",
                "0.420000",
                "1.000000",
                "0.0",
                "chrT",
                "900",
                "1200",
                "+",
                "chrT",
                "900",
                "1200",
                "+",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "no_accepted_structure",
            ],
            [
                "gene:ENSG0003.1",
                "",
                "0.000000",
                "0",
                "0.000000",
                "0.000000",
                "0.000000",
                "0.000000",
                "0.000000",
                "0.0",
                "chrT",
                "1500",
                "1800",
                "-",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "no_locus_candidate",
            ],
        ],
    )

    summary = audit_run(run_dir, ref_gff, target_gff, output_dir)

    assert summary["gene_counts"]["structural_mapped"] == 1
    assert summary["gene_counts"]["coordinate_mapped"] == 1
    assert summary["gene_counts"]["missing"] == 1
    assert summary["gene_counts"]["new"] == 1
    assert summary["locus_rows_loaded"] == 2
    assert summary["new_target_id_in_ref"]["yes"] == 1
    assert summary["new_target_id_in_ref_old_action"]["missing"] == 1
    assert (output_dir / "missing_genes.tsv").exists()
    assert (output_dir / "coordinate_mapped_genes.tsv").exists()
    assert (output_dir / "new_genes.tsv").exists()
    assert "ENSG0003" in (output_dir / "missing_genes.tsv").read_text()
    assert "ENSG9002" in (output_dir / "coordinate_mapped_genes.tsv").read_text()
    assert "ref_same_id_action" in (output_dir / "new_genes.tsv").read_text()
    assert "missing" in (output_dir / "new_genes.tsv").read_text()

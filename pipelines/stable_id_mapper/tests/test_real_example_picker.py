from __future__ import annotations

import sys
import csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pick_real_stable_id_examples import pick_examples, write_examples


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def test_pick_real_examples_from_completed_run_outputs(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    write_text(
        run_dir / "stable_id_decisions.tsv",
        """
type	action	current_stable_id	current_version	old_stable_id	old_version	new_stable_id	new_version	mapping_session_id	score	reason
gene	mapped	ENSG9001	1	ENSG0001	1	ENSG0001	2	1	0.91	mapped gene matched a target gene by lifton structural evidence; score_source=lifton_gene_aggregate confidence=high
gene	mapped	ENSG9002	1	ENSG0002	1	ENSG0002	2	1	0.52	mapped gene matched a target gene by coordinate overlap; score_source=coordinate_overlap
gene	missing		0	ENSG0003	1		0	1	0	mapped gene did not match any target gene by lifton evidence or coordinate overlap
transcript	mapped	ENST9001	1	ENST0001	1	ENST0001	2	1	0.88	mapped transcript matched a target transcript by lifton structural evidence; score_source=lifton_transcript_structure confidence=high
        """,
    )
    locus_header = [
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
    ]
    write_tsv(
        run_dir / "matching" / "lifton.gene_locus_comparison.tsv",
        locus_header,
        [
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
                "chr1",
                "100",
                "500",
                "+",
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
            [
                "gene:ENSG0004.1",
                "gene:ENSG9004.1",
                "0.800000",
                "401",
                "1.000000",
                "0.800000",
                "0.800000",
                "1.000000",
                "0.800000",
                "0.0",
                "chr1",
                "100",
                "500",
                "+",
                "chr1",
                "100",
                "600",
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
                "gene:ENSG0005.1",
                "gene:ENSG9005.1",
                "0.700000",
                "301",
                "0.700000",
                "0.700000",
                "0.700000",
                "0.700000",
                "1.000000",
                "0.0",
                "chr1",
                "100",
                "530",
                "+",
                "chr1",
                "100",
                "530",
                "+",
                "gene:ENSG9015.1",
                "0.600000",
                "1.000000",
                "1",
                "gene:ENSG9015.1",
                "transcript:ENST9015.1",
                "0.600000",
                "different",
            ],
        ],
    )

    examples = pick_examples(run_dir, limit=2)
    output = tmp_path / "examples.tsv"
    write_examples(examples, output)
    text = output.read_text(encoding="utf-8")

    assert "gene_structural_mapped" in text
    assert "gene_coordinate_mapped" in text
    assert "transcript_structural_mapped" in text
    assert "high_locus_no_accepted_structure" in text
    assert "locus_structure_disagree" in text
    assert "missing_projected_gene_no_locus" in text

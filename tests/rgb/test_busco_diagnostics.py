import pandas as pd

from ensembl.genes.rgb.busco_diagnostics import (
    build_busco_diagnostic_audit,
    identify_busco_problems,
    parse_busco_table,
)


def test_busco_tables_identify_genome_protein_regressions(tmp_path):
    genome = tmp_path / "genome_full_table.tsv"
    protein = tmp_path / "protein_full_table.tsv"
    genome.write_text(
        "\n".join(
            [
                "# BUSCO version",
                "123at456\tComplete\tchr1\t100\t500\t+\t99.0\t250",
                "999at456\tFragmented\tchr2\t1000\t1200\t-\t50.0\t90",
                "777at456\tComplete\tchr3\t1000\t1300\t+\t80.0\t120",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protein.write_text(
        "\n".join(
            [
                "123at456\tMissing",
                "999at456\tMissing",
                "777at456\tComplete\tGENE1\t90.0\t120",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    problems = identify_busco_problems(
        parse_busco_table(str(genome), "genome"),
        parse_busco_table(str(protein), "protein"),
    )

    assert set(problems["busco_id"]) == {"123at456", "999at456"}
    assert (
        problems.set_index("busco_id").loc["123at456", "problem_type"]
        == "genome_complete_protein_missing"
    )


def test_busco_diagnostic_classifies_layer_evidence_not_selected(tmp_path):
    genome = pd.DataFrame(
        [
            {
                "busco_id": "123at456",
                "status": "Complete",
                "orthodb_group": "123",
                "taxid": "456",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 500,
                "seq_region_strand": "+",
            }
        ]
    )
    protein = pd.DataFrame(
        [
            {
                "busco_id": "123at456",
                "status": "Missing",
                "orthodb_group": "123",
                "taxid": "456",
            }
        ]
    )
    layer_genes = pd.DataFrame(
        [
            {
                "gene_id": 10,
                "stable_id": "LAYER1",
                "seq_region_name": "chr1",
                "seq_region_start": 120,
                "seq_region_end": 480,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "protein2genome_high",
            }
        ]
    )

    audit, isoforms, summary = build_busco_diagnostic_audit(
        genome,
        protein,
        pd.DataFrame(),
        layer_genes,
        pd.DataFrame(),
        pd.DataFrame(),
    )

    assert isoforms.empty
    assert (
        summary.set_index("diagnostic_class").loc[
            "layer_evidence_not_selected", "count"
        ]
        == 1
    )
    row = audit.iloc[0]
    assert row["diagnostic_class"] == "layer_evidence_not_selected"
    assert row["layer_logic_names"] == "protein2genome_high(1)"
    assert row["suggested_action"] == "rescue_or_reprioritise_layer_candidate"


def test_busco_diagnostic_reports_good_hmmer_isoform(tmp_path):
    hmmer_dir = tmp_path / "hmmer"
    hmmer_dir.mkdir()
    (hmmer_dir / "123at456.out").write_text(
        "ENSDART0001.2 - 300 123at456 - 250 1e-60 200 0 1 1 1e-60 1e-60 200 0 1 240 1 240\n",
        encoding="utf-8",
    )
    genome = pd.DataFrame(
        [
            {
                "busco_id": "123at456",
                "status": "Complete",
                "orthodb_group": "123",
                "taxid": "456",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 500,
                "seq_region_strand": "+",
            }
        ]
    )
    protein = pd.DataFrame(
        [
            {
                "busco_id": "123at456",
                "status": "Fragmented",
                "orthodb_group": "123",
                "taxid": "456",
            }
        ]
    )
    core_genes = pd.DataFrame(
        [
            {
                "gene_id": 1,
                "stable_id": "ENSDARG0001",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 500,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
            }
        ]
    )
    core_tx = pd.DataFrame(
        [{"gene_id": 1, "transcript_id": 101, "stable_id": "ENSDART0001.2"}]
    )
    core_tr = pd.DataFrame([{"transcript_id": 101, "translation_id": 201}])

    audit, isoforms, _summary = build_busco_diagnostic_audit(
        genome,
        protein,
        core_genes,
        pd.DataFrame(),
        core_tx,
        core_tr,
        hmmer_dir=str(hmmer_dir),
    )

    assert audit.iloc[0]["diagnostic_class"] == "alternative_isoform_satisfies_busco"
    assert audit.iloc[0]["good_isoform_count"] == 1
    assert isoforms.iloc[0]["transcript_id"] == "ENSDART0001.2"

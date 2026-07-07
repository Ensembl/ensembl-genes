from __future__ import annotations

import sys
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import run_stable_id_mapping
import stable_id_mapping.workflow as workflow
import lifton_id_mapper
from stable_id_mapping.ids import parse_id_range
from stable_id_mapping.lifton_matching import LiftonMatchConfig, LiftonMatchSummary
from stable_id_mapping.lifton_matching import run_lifton_matching
from stable_id_mapping.reports import write_missing_gene_report
from stable_id_mapping.workflow import SingleSpeciesRunConfig, run_single_species_pipeline


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def test_write_missing_gene_report(tmp_path: Path) -> None:
    ref_gff = tmp_path / "ref.gff3"
    projected_gff = tmp_path / "projected.gff3"
    report = tmp_path / "report.txt"
    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	10	100	.	+	.	ID=gene:ENSG0001.1
chr1	test	mRNA	10	100	.	+	.	ID=transcript:ENST0001.1;Parent=gene:ENSG0001.1
chr1	test	gene	200	300	.	+	.	ID=gene:ENSG0002.1
chr1	test	mRNA	200	300	.	+	.	ID=transcript:ENST0002.1;Parent=gene:ENSG0002.1
        """,
    )
    write_text(
        projected_gff,
        """
##gff-version 3
chrT	test	gene	10	100	.	+	.	ID=gene:ENSG0001.1
chrT	test	mRNA	10	100	.	+	.	ID=transcript:ENST0001.1;Parent=gene:ENSG0001.1
        """,
    )

    missing = write_missing_gene_report(ref_gff, projected_gff, report)

    assert missing == {"ENSG0002"}
    assert "gene:ENSG0002" in report.read_text(encoding="utf-8")


def test_structural_parser_uses_stable_id_feature_universe(tmp_path: Path) -> None:
    gff = tmp_path / "annotation.gff3"
    write_text(
        gff,
        """
##gff-version 3
chrT	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.1
chrT	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.1;Parent=gene:ENSXG0001.1
chrT	test	exon	100	200	.	+	.	ID=exon:one;Parent=transcript:ENSXT0001.1
chrT	test	ncRNA_gene	700	900	.	+	.	ID=gene:ENSXGNC1.1
chrT	test	lnc_RNA	700	900	.	+	.	ID=transcript:ENSXTNC1.1;Parent=gene:ENSXGNC1.1
chrT	test	exon	700	900	.	+	.	ID=exon:nc;Parent=transcript:ENSXTNC1.1
        """,
    )

    annotation = lifton_id_mapper.load_gff3_as_annotation(str(gff), "test")

    assert set(annotation.genes) == {"gene:ENSXG0001.1"}
    assert set(annotation.tx_index) == {"transcript:ENSXT0001.1"}


def test_lifton_matching_uses_cds_when_lifton_output_has_no_exons(
    tmp_path: Path,
) -> None:
    lifton_gff = tmp_path / "lifton.gff3"
    target_gff = tmp_path / "target.gff3"
    out_prefix = tmp_path / "matching" / "lifton"
    write_text(
        lifton_gff,
        """
##gff-version 3
chrT	LiftOn	gene	100	500	.	+	.	ID=gene:ENSXG0001.1
chrT	LiftOn	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.1;Parent=gene:ENSXG0001.1
chrT	LiftOn	CDS	100	200	.	+	0	ID=CDS:ENSXP0001;Parent=transcript:ENSXT0001.1
chrT	LiftOn	CDS	300	500	.	+	0	ID=CDS:ENSXP0001;Parent=transcript:ENSXT0001.1
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	50	550	.	+	.	ID=gene:ENSXG9001.1
chrT	test	mRNA	50	550	.	+	.	ID=transcript:ENSXT9001.1;Parent=gene:ENSXG9001.1
chrT	test	exon	50	250	.	+	.	ID=exon:one;Parent=transcript:ENSXT9001.1
chrT	test	exon	280	550	.	+	.	ID=exon:two;Parent=transcript:ENSXT9001.1
        """,
    )

    summary = run_lifton_matching(
        LiftonMatchConfig(
            lifton_gff=lifton_gff,
            target_gff=target_gff,
            out_prefix=out_prefix,
            min_score=0.60,
        )
    )

    assert summary.transcript_pairs == 1
    assert summary.gene_pairs == 1
    transcript_pairs_text = summary.transcript_pairs_path.read_text(encoding="utf-8")
    assert "query_coverage" in transcript_pairs_text
    assert "span_containment" in transcript_pairs_text
    assert "transcript:ENSXT0001" in transcript_pairs_text
    assert "gene:ENSXG0001" in summary.gene_pairs_path.read_text(encoding="utf-8")
    assert summary.gene_locus_comparison_path is not None
    locus_text = summary.gene_locus_comparison_path.read_text(encoding="utf-8")
    assert "locus_score" in locus_text
    assert "old_gene_coverage" in locus_text
    assert "structure_accepted_target_gene" in locus_text
    assert "gene:ENSXG0001.1\tgene:ENSXG9001.1" in locus_text


def test_single_species_workflow_with_patched_lifton(
    tmp_path: Path,
    monkeypatch,
) -> None:
    ref_fasta = tmp_path / "ref.fa"
    target_fasta = tmp_path / "target.fa"
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    for fasta in (ref_fasta, target_fasta):
        write_text(fasta, ">chr1\nAAAA")
    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	10	100	.	+	.	ID=gene:ENSXG0001.2
chr1	test	mRNA	10	100	.	+	.	ID=transcript:ENSXT0001.3;Parent=gene:ENSXG0001.2
chr1	test	gene	200	300	.	+	.	ID=gene:ENSXG0002.1
chr1	test	mRNA	200	300	.	+	.	ID=transcript:ENSXT0002.1;Parent=gene:ENSXG0002.1
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	10	100	.	+	.	ID=gene:ENSXG9001.1
chrT	test	mRNA	10	100	.	+	.	ID=transcript:ENSXT9001.1;Parent=gene:ENSXG9001.1
chrT	test	gene	400	500	.	+	.	ID=gene:ENSXG9002.1
chrT	test	mRNA	400	500	.	+	.	ID=transcript:ENSXT9002.1;Parent=gene:ENSXG9002.1
        """,
    )

    def fake_run_lifton(config):
        write_text(
            config.output_gff,
            """
##gff-version 3
chrT	test	gene	10	100	.	+	.	ID=gene:ENSXG0001.2
chrT	test	mRNA	10	100	.	+	.	ID=transcript:ENSXT0001.3;Parent=gene:ENSXG0001.2
            """,
        )

    def fake_run_lifton_matching(config):
        write_text(
            config.transcript_pairs_path,
            """
lifton_tx	ref_tx	score	status	intron_sim	jacc_internal	jacc_all	exon_count_sim	boundary_sim	lifton_identity_prior	lifton_gene	ref_gene	contig	strand	lifton_exons	ref_exons
transcript:ENSXT0001	transcript:ENSXT9001	0.880000	confident	1.000000	0.900000	0.900000	1.000000	0.900000	0.000000	gene:ENSXG0001	gene:ENSXG9001	chrT	+	2	2
            """,
        )
        write_text(
            config.gene_pairs_path,
            """
lifton_gene	ref_gene	weighted_score	fraction_of_total	n_transcripts
gene:ENSXG0001	gene:ENSXG9001	0.880000	1.000000	1
            """,
        )
        return LiftonMatchSummary(
            transcript_pairs=1,
            gene_pairs=1,
            transcript_pairs_path=config.transcript_pairs_path,
            gene_pairs_path=config.gene_pairs_path,
        )

    monkeypatch.setattr(workflow, "run_lifton", fake_run_lifton)
    monkeypatch.setattr(workflow, "run_lifton_matching", fake_run_lifton_matching)

    result = run_single_species_pipeline(
        SingleSpeciesRunConfig(
            ref_fasta=ref_fasta,
            ref_gff=ref_gff,
            target_fasta=target_fasta,
            target_gff=target_gff,
            db_name="species_core",
            mapping_session_id=1,
            gene_range=parse_id_range("ENSXG:7000-7999"),
            transcript_range=parse_id_range("ENSXT:7000-7999"),
            translation_range=parse_id_range("ENSXP:7000-7999"),
            output_dir=tmp_path / "out",
            include_translations=False,
            dry_run_sql=True,
        )
    )

    assert result.lifton_gff.exists()
    assert result.missing_report.exists()
    assert result.match_summary.transcript_pairs == 1
    assert result.output_sql.exists()
    assert result.output_score_evidence_tsv.exists()
    assert "USE `species_core`;" in result.output_sql.read_text(encoding="utf-8")
    assert any(
        decision.feature_type == "gene"
        and decision.action == "mapped"
        and decision.new_version == 3
        and decision.score == 0.88
        for decision in result.decisions
    )
    assert any(
        decision.feature_type == "transcript"
        and decision.action == "mapped"
        and decision.score == 0.88
        for decision in result.decisions
    )


def test_top_level_cli_builds_expected_lifton_command(tmp_path: Path) -> None:
    ref_fasta = tmp_path / "ref.fa"
    target_fasta = tmp_path / "target.fa"
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    for path in (ref_fasta, target_fasta, ref_gff, target_gff):
        write_text(path, "")

    args = run_stable_id_mapping.parse_args(
        [
            "--ref-fasta",
            str(ref_fasta),
            "--ref-gff",
            str(ref_gff),
            "--target-fasta",
            str(target_fasta),
            "--target-gff",
            str(target_gff),
            "--db-name",
            "species_core",
            "--mapping-session-id",
            "4",
            "--gene-range",
            "ENSXG:1-10",
            "--transcript-range",
            "ENSXT:1-10",
            "--translation-range",
            "ENSXP:1-10",
            "--output-dir",
            str(tmp_path / "out"),
            "--lifton-threads",
            "12",
        ]
    )
    config = run_stable_id_mapping.config_from_args(args)

    command = run_stable_id_mapping.lifton_command_for_config(
        config,
        prepare_feature_types=True,
    )
    feature_types_file = tmp_path / "out" / "lifton" / "lifton_feature_types.txt"

    assert command[:6] == [
        "lifton",
        "-f",
        str(feature_types_file),
        "-t",
        "12",
        "-g",
    ]
    assert command[-2:] == [str(ref_fasta), str(target_fasta)]
    assert feature_types_file.read_text(encoding="utf-8") == "gene\n"


def test_top_level_cli_loads_rules_config_with_flag_overrides(tmp_path: Path) -> None:
    ref_fasta = tmp_path / "ref.fa"
    target_fasta = tmp_path / "target.fa"
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    rules_path = tmp_path / "rules.json"
    for path in (ref_fasta, target_fasta, ref_gff, target_gff):
        write_text(path, "")
    rules_path.write_text(
        json.dumps(
            {
                "coordinate_overlap": {"min_overlap": 0.2},
                "structural_matching": {
                    "window": 123,
                    "topk": 2,
                    "min_score": 0.41,
                    "good_score": 0.61,
                    "confident_score": 0.81,
                    "gene_fraction": 0.71,
                    "score_weights": {
                        "query_coverage": 2,
                        "span_containment": 1,
                    },
                },
            }
        ),
        encoding="utf-8",
    )

    args = run_stable_id_mapping.parse_args(
        [
            "--ref-fasta",
            str(ref_fasta),
            "--ref-gff",
            str(ref_gff),
            "--target-fasta",
            str(target_fasta),
            "--target-gff",
            str(target_gff),
            "--db-name",
            "species_core",
            "--mapping-session-id",
            "4",
            "--gene-range",
            "ENSXG:1-10",
            "--transcript-range",
            "ENSXT:1-10",
            "--translation-range",
            "ENSXP:1-10",
            "--output-dir",
            str(tmp_path / "out"),
            "--rules-config",
            str(rules_path),
            "--match-min-score",
            "0.5",
        ]
    )
    config = run_stable_id_mapping.config_from_args(args)

    assert config.rules_config == rules_path
    assert config.min_overlap == 0.2
    assert config.match_window == 123
    assert config.match_topk == 2
    assert config.match_min_score == 0.5
    assert config.match_good == 0.61
    assert config.match_confident == 0.81
    assert config.match_gene_fraction == 0.71
    assert config.match_score_weights["query_coverage"] == 2


def test_single_species_workflow_reuses_existing_outputs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    ref_fasta = tmp_path / "ref.fa"
    target_fasta = tmp_path / "target.fa"
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    lifton_gff = tmp_path / "existing" / "projected.gff3"
    transcript_pairs = tmp_path / "existing" / "lifton.transcript_pairs.tsv"
    gene_pairs = tmp_path / "existing" / "lifton.gene_pairs.tsv"
    for fasta in (ref_fasta, target_fasta):
        write_text(fasta, ">chr1\nAAAA")
    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	10	100	.	+	.	ID=gene:ENSXG0001.1
chr1	test	mRNA	10	100	.	+	.	ID=transcript:ENSXT0001.1;Parent=gene:ENSXG0001.1
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	10	100	.	+	.	ID=gene:ENSXG9001.1
chrT	test	mRNA	10	100	.	+	.	ID=transcript:ENSXT9001.1;Parent=gene:ENSXG9001.1
        """,
    )
    write_text(
        lifton_gff,
        """
##gff-version 3
chrT	test	gene	10	100	.	+	.	ID=gene:ENSXG0001.1
chrT	test	mRNA	10	100	.	+	.	ID=transcript:ENSXT0001.1;Parent=gene:ENSXG0001.1
        """,
    )
    write_text(
        transcript_pairs,
        """
lifton_tx	ref_tx	score	status	intron_sim	jacc_internal	jacc_all	exon_count_sim	boundary_sim	lifton_identity_prior	lifton_gene	ref_gene	contig	strand	lifton_exons	ref_exons
transcript:ENSXT0001	transcript:ENSXT9001	0.770000	good	1	1	1	1	1	0	gene:ENSXG0001	gene:ENSXG9001	chrT	+	1	1
        """,
    )
    write_text(
        gene_pairs,
        """
lifton_gene	ref_gene	weighted_score	fraction_of_total	n_transcripts
gene:ENSXG0001	gene:ENSXG9001	0.770000	1.000000	1
        """,
    )

    def fail_run_lifton(_config):
        raise AssertionError("LiftOn should not run when existing_lifton_gff is provided")

    def fail_run_lifton_matching(_config):
        raise AssertionError("Matching should not run when existing pair tables are provided")

    monkeypatch.setattr(workflow, "run_lifton", fail_run_lifton)
    monkeypatch.setattr(workflow, "run_lifton_matching", fail_run_lifton_matching)

    result = run_single_species_pipeline(
        SingleSpeciesRunConfig(
            ref_fasta=ref_fasta,
            ref_gff=ref_gff,
            target_fasta=target_fasta,
            target_gff=target_gff,
            db_name="species_core",
            mapping_session_id=1,
            gene_range=parse_id_range("ENSXG:7000-7999"),
            transcript_range=parse_id_range("ENSXT:7000-7999"),
            translation_range=parse_id_range("ENSXP:7000-7999"),
            output_dir=tmp_path / "out",
            existing_lifton_gff=lifton_gff,
            existing_transcript_pairs=transcript_pairs,
            existing_gene_pairs=gene_pairs,
            include_translations=False,
            dry_run_sql=True,
        )
    )

    assert result.lifton_gff == lifton_gff
    assert result.match_summary.transcript_pairs == 1
    assert result.match_summary.gene_pairs == 1
    assert any(
        decision.feature_type == "gene"
        and decision.action == "mapped"
        and decision.score == 0.77
        for decision in result.decisions
    )


def test_top_level_cli_dry_run_lifton_command_reports_reuse(tmp_path: Path) -> None:
    existing_lifton_gff = tmp_path / "projected.gff3"
    ref_fasta = tmp_path / "ref.fa"
    target_fasta = tmp_path / "target.fa"
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    for path in (existing_lifton_gff, ref_fasta, target_fasta, ref_gff, target_gff):
        write_text(path, "")

    args = run_stable_id_mapping.parse_args(
        [
            "--ref-fasta",
            str(ref_fasta),
            "--ref-gff",
            str(ref_gff),
            "--target-fasta",
            str(target_fasta),
            "--target-gff",
            str(target_gff),
            "--db-name",
            "species_core",
            "--mapping-session-id",
            "4",
            "--gene-range",
            "ENSXG:1-10",
            "--transcript-range",
            "ENSXT:1-10",
            "--translation-range",
            "ENSXP:1-10",
            "--output-dir",
            str(tmp_path / "out"),
            "--existing-lifton-gff",
            str(existing_lifton_gff),
        ]
    )
    config = run_stable_id_mapping.config_from_args(args)

    assert run_stable_id_mapping.lifton_command_for_config(config) == []

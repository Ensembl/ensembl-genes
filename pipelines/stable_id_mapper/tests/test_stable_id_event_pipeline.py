from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stable_id_mapping.config import StableIdEventConfig
from stable_id_mapping.gff3 import parse_gff3
from stable_id_mapping.ids import collect_reserved_ids, make_allocator, parse_id_range
from stable_id_mapping.pipeline import run_pipeline
from stable_id_mapping.scoring import MappingEvidence, ScoreEvidence


def write_text(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def test_parse_gff3_ignores_region_records(tmp_path: Path) -> None:
    gff = tmp_path / "input.gff3"
    write_text(
        gff,
        """
##gff-version 3
chr1	test	region	1	1000	.	.	.	ID=region:chr1
chr1	test	gene	10	100	.	+	.	ID=gene:ENSG0001.3
chr1	test	mRNA	10	100	.	+	.	ID=transcript:ENST0001.7;Parent=gene:ENSG0001.3
chr1	test	CDS	20	80	.	+	0	ID=CDS:ENSP0001;Parent=transcript:ENST0001.7;protein_id=ENSP0001
        """,
    )

    features = parse_gff3(gff)

    assert set(features["gene"]) == {"ENSG0001"}
    assert set(features["transcript"]) == {"ENST0001"}
    assert set(features["translation"]) == {"ENSP0001"}


def test_stable_id_event_pipeline_versions_and_sql(tmp_path: Path) -> None:
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    mapped_gff = tmp_path / "mapped.gff3"
    report = tmp_path / "report.txt"
    output_sql = tmp_path / "out.sql"
    output_tsv = tmp_path / "out.tsv"

    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.4
chr1	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.9;Parent=gene:ENSXG0001.4
chr1	test	CDS	150	450	.	+	0	ID=CDS:ENSXP0001;Parent=transcript:ENSXT0001.9;protein_id=ENSXP0001.2
chr1	test	gene	900	950	.	+	.	ID=gene:ENSXG0002.1
chr1	test	mRNA	900	950	.	+	.	ID=transcript:ENSXT0002.1;Parent=gene:ENSXG0002.1
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	100	500	.	+	.	ID=gene:ENSXG9001.1
chrT	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT9001.1;Parent=gene:ENSXG9001.1
chrT	test	CDS	150	450	.	+	0	ID=CDS:ENSXP9001;Parent=transcript:ENSXT9001.1;protein_id=ENSXP9001.1
chrT	test	gene	1000	1100	.	+	.	ID=gene:ENSXG9002.1
chrT	test	mRNA	1000	1100	.	+	.	ID=transcript:ENSXT9002.1;Parent=gene:ENSXG9002.1
        """,
    )
    write_text(
        mapped_gff,
        """
##gff-version 3
chrT	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.4
chrT	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.9;Parent=gene:ENSXG0001.4
chrT	test	CDS	150	450	.	+	0	ID=CDS:ENSXP0001;Parent=transcript:ENSXT0001.9;protein_id=ENSXP0001.2
        """,
    )
    write_text(
        report,
        """
ID Mapping Report

Missing gene IDs:
  gene:ENSXG0002
        """,
    )

    config = StableIdEventConfig(
        ref_gff=ref_gff,
        target_gff=target_gff,
        mapped_gff=mapped_gff,
        report=report,
        db_name="species_core",
        mapping_session_id=17,
        gene_range=parse_id_range("ENSXG:7000-7999"),
        transcript_range=parse_id_range("ENSXT:7000-7999"),
        translation_range=parse_id_range("ENSXP:7000-7999"),
        output_sql=output_sql,
        output_tsv=output_tsv,
        include_translations=True,
        dry_run=True,
        backup_prefix="stable_id_mapper_backup_test",
    )

    decisions = run_pipeline(config)
    by_key = {
        (decision.feature_type, decision.action, decision.old_stable_id): decision
        for decision in decisions
    }

    assert by_key[("gene", "mapped", "ENSXG0001")].new_version == 5
    assert by_key[("transcript", "mapped", "ENSXT0001")].new_version == 10
    assert by_key[("translation", "mapped", "ENSXP0001")].new_version == 3
    assert by_key[("gene", "missing", "ENSXG0002")].new_stable_id is None
    assert "USE `species_core`;" in output_sql.read_text(encoding="utf-8")
    assert "ROLLBACK;" in output_sql.read_text(encoding="utf-8")
    assert output_tsv.read_text(encoding="utf-8").splitlines()[0].startswith(
        "type\taction\tcurrent_stable_id"
    )


def test_pipeline_maps_by_lifton_evidence_before_coordinate_overlap(
    tmp_path: Path,
) -> None:
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    mapped_gff = tmp_path / "mapped.gff3"
    report = tmp_path / "report.txt"

    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.4
chr1	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.9;Parent=gene:ENSXG0001.4
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	1000	1500	.	+	.	ID=gene:ENSXG9001.1
chrT	test	mRNA	1000	1500	.	+	.	ID=transcript:ENSXT9001.1;Parent=gene:ENSXG9001.1
        """,
    )
    write_text(
        mapped_gff,
        """
##gff-version 3
chrT	test	gene	5000	5500	.	+	.	ID=gene:ENSXG0001.4
chrT	test	mRNA	5000	5500	.	+	.	ID=transcript:ENSXT0001.9;Parent=gene:ENSXG0001.4
        """,
    )
    write_text(report, "Missing gene IDs:\n")

    evidence = ScoreEvidence()
    evidence.add(
        MappingEvidence(
            feature_type="gene",
            old_stable_id="ENSXG0001",
            target_stable_id="ENSXG9001",
            score=0.91,
            source="lifton_gene_aggregate",
            confidence="high",
        )
    )
    evidence.add(
        MappingEvidence(
            feature_type="transcript",
            old_stable_id="ENSXT0001",
            target_stable_id="ENSXT9001",
            score=0.88,
            source="lifton_transcript_structure",
            confidence="high",
        )
    )

    decisions = run_pipeline(
        StableIdEventConfig(
            ref_gff=ref_gff,
            target_gff=target_gff,
            mapped_gff=mapped_gff,
            report=report,
            mapping_session_id=17,
            gene_range=parse_id_range("ENSXG:7000-7999"),
            transcript_range=parse_id_range("ENSXT:7000-7999"),
            translation_range=parse_id_range("ENSXP:7000-7999"),
            output_sql=tmp_path / "out.sql",
            include_translations=False,
            dry_run=True,
            score_evidence=evidence,
        )
    )

    by_key = {
        (decision.feature_type, decision.action, decision.old_stable_id): decision
        for decision in decisions
    }

    assert by_key[("gene", "mapped", "ENSXG0001")].current_stable_id == "ENSXG9001"
    assert by_key[("gene", "mapped", "ENSXG0001")].score == 0.91
    assert (
        by_key[("transcript", "mapped", "ENSXT0001")].current_stable_id
        == "ENSXT9001"
    )
    assert not any(
        decision.feature_type == "gene" and decision.action == "new"
        for decision in decisions
    )
    assert not any(
        decision.feature_type == "transcript" and decision.action == "new"
        for decision in decisions
    )


def test_pipeline_does_not_map_projected_ids_without_a_target_match(
    tmp_path: Path,
) -> None:
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    mapped_gff = tmp_path / "mapped.gff3"
    report = tmp_path / "report.txt"

    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.4
chr1	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.9;Parent=gene:ENSXG0001.4
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	1000	1500	.	+	.	ID=gene:ENSXG9001.1
chrT	test	mRNA	1000	1500	.	+	.	ID=transcript:ENSXT9001.1;Parent=gene:ENSXG9001.1
        """,
    )
    write_text(
        mapped_gff,
        """
##gff-version 3
chrT	test	gene	5000	5500	.	+	.	ID=gene:ENSXG0001.4
chrT	test	mRNA	5000	5500	.	+	.	ID=transcript:ENSXT0001.9;Parent=gene:ENSXG0001.4
        """,
    )
    write_text(report, "Missing gene IDs:\n")

    decisions = run_pipeline(
        StableIdEventConfig(
            ref_gff=ref_gff,
            target_gff=target_gff,
            mapped_gff=mapped_gff,
            report=report,
            mapping_session_id=17,
            gene_range=parse_id_range("ENSXG:7000-7999"),
            transcript_range=parse_id_range("ENSXT:7000-7999"),
            translation_range=parse_id_range("ENSXP:7000-7999"),
            output_sql=tmp_path / "out.sql",
            include_translations=False,
            dry_run=True,
        )
    )

    assert not any(
        decision.feature_type == "gene" and decision.action == "mapped"
        for decision in decisions
    )
    assert any(
        decision.feature_type == "gene"
        and decision.action == "missing"
        and decision.old_stable_id == "ENSXG0001"
        for decision in decisions
    )
    assert any(
        decision.feature_type == "gene"
        and decision.action == "new"
        and decision.current_stable_id == "ENSXG9001"
        for decision in decisions
    )


def test_gene_evidence_conflicts_are_resolved_by_score_before_overlap(
    tmp_path: Path,
) -> None:
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    mapped_gff = tmp_path / "mapped.gff3"
    report = tmp_path / "report.txt"

    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	100	200	.	+	.	ID=gene:ENSXG0001.1
chr1	test	gene	300	400	.	+	.	ID=gene:ENSXG0002.1
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	100	200	.	+	.	ID=gene:ENSXG9001.1
        """,
    )
    write_text(
        mapped_gff,
        """
##gff-version 3
chrT	test	gene	100	200	.	+	.	ID=gene:ENSXG0001.1
chrT	test	gene	100	200	.	+	.	ID=gene:ENSXG0002.1
        """,
    )
    write_text(report, "Missing gene IDs:\n")

    evidence = ScoreEvidence()
    evidence.add(
        MappingEvidence(
            feature_type="gene",
            old_stable_id="ENSXG0001",
            target_stable_id="ENSXG9001",
            score=0.2,
            source="lifton_gene_aggregate",
            confidence="review",
        )
    )
    evidence.add(
        MappingEvidence(
            feature_type="gene",
            old_stable_id="ENSXG0002",
            target_stable_id="ENSXG9001",
            score=0.9,
            source="lifton_gene_aggregate",
            confidence="high",
        )
    )

    decisions = run_pipeline(
        StableIdEventConfig(
            ref_gff=ref_gff,
            target_gff=target_gff,
            mapped_gff=mapped_gff,
            report=report,
            mapping_session_id=17,
            gene_range=parse_id_range("ENSXG:7000-7999"),
            transcript_range=parse_id_range("ENSXT:7000-7999"),
            translation_range=parse_id_range("ENSXP:7000-7999"),
            output_sql=tmp_path / "out.sql",
            include_translations=False,
            dry_run=True,
            score_evidence=evidence,
        )
    )

    mapped = [
        decision
        for decision in decisions
        if decision.feature_type == "gene" and decision.action == "mapped"
    ]

    assert len(mapped) == 1
    assert mapped[0].old_stable_id == "ENSXG0002"
    assert mapped[0].current_stable_id == "ENSXG9001"
    assert mapped[0].score == 0.9
    assert any(
        decision.feature_type == "gene"
        and decision.action == "missing"
        and decision.old_stable_id == "ENSXG0001"
        for decision in decisions
    )


def test_gene_same_stable_id_target_does_not_win_without_evidence(
    tmp_path: Path,
) -> None:
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    mapped_gff = tmp_path / "mapped.gff3"
    report = tmp_path / "report.txt"

    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.4
chr1	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.4;Parent=gene:ENSXG0001.4
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	100	200	.	+	.	ID=gene:ENSXG0001.5
chrT	test	mRNA	100	200	.	+	.	ID=transcript:ENSXT0001.5;Parent=gene:ENSXG0001.5
chrT	test	gene	100	500	.	+	.	ID=gene:ENSXG9001.1
chrT	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT9001.1;Parent=gene:ENSXG9001.1
        """,
    )
    write_text(
        mapped_gff,
        """
##gff-version 3
chrT	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.4
chrT	test	mRNA	100	500	.	+	.	ID=transcript:ENSXT0001.4;Parent=gene:ENSXG0001.4
        """,
    )
    write_text(report, "Missing gene IDs:\n")

    evidence = ScoreEvidence()
    evidence.add(
        MappingEvidence(
            feature_type="gene",
            old_stable_id="ENSXG0001",
            target_stable_id="ENSXG9001",
            score=0.95,
            source="lifton_gene_aggregate",
            confidence="high",
        )
    )

    decisions = run_pipeline(
        StableIdEventConfig(
            ref_gff=ref_gff,
            target_gff=target_gff,
            mapped_gff=mapped_gff,
            report=report,
            mapping_session_id=17,
            gene_range=parse_id_range("ENSXG:7000-7999"),
            transcript_range=parse_id_range("ENSXT:7000-7999"),
            translation_range=parse_id_range("ENSXP:7000-7999"),
            output_sql=tmp_path / "out.sql",
            include_translations=False,
            dry_run=True,
            score_evidence=evidence,
        )
    )

    gene_decisions = [
        decision for decision in decisions if decision.feature_type == "gene"
    ]
    mapped_gene = next(
        decision for decision in gene_decisions if decision.action == "mapped"
    )
    new_same_id_target = next(
        decision
        for decision in gene_decisions
        if decision.action == "new" and decision.current_stable_id == "ENSXG0001"
    )

    assert mapped_gene.old_stable_id == "ENSXG0001"
    assert mapped_gene.current_stable_id == "ENSXG9001"
    assert mapped_gene.new_version == 5
    assert "by lifton structural evidence" in mapped_gene.reason
    assert new_same_id_target.new_stable_id == "ENSXG7000"


def test_lifton_copy_gene_ids_do_not_become_old_stable_id_decisions(
    tmp_path: Path,
) -> None:
    ref_gff = tmp_path / "ref.gff3"
    target_gff = tmp_path / "target.gff3"
    mapped_gff = tmp_path / "mapped.gff3"
    report = tmp_path / "report.txt"

    write_text(
        ref_gff,
        """
##gff-version 3
chr1	test	gene	100	500	.	+	.	ID=gene:ENSXG0001.4
        """,
    )
    write_text(
        target_gff,
        """
##gff-version 3
chrT	test	gene	100	500	.	+	.	ID=gene:ENSXG9001.1
        """,
    )
    write_text(
        mapped_gff,
        """
##gff-version 3
chrT	test	gene	100	500	.	+	.	ID=gene:ENSXG0001_1.4
        """,
    )
    write_text(report, "Missing gene IDs:\n")

    decisions = run_pipeline(
        StableIdEventConfig(
            ref_gff=ref_gff,
            target_gff=target_gff,
            mapped_gff=mapped_gff,
            report=report,
            mapping_session_id=17,
            gene_range=parse_id_range("ENSXG:7000-7999"),
            transcript_range=parse_id_range("ENSXT:7000-7999"),
            translation_range=parse_id_range("ENSXP:7000-7999"),
            output_sql=tmp_path / "out.sql",
            include_translations=False,
            dry_run=True,
        )
    )

    gene_decisions = [
        decision for decision in decisions if decision.feature_type == "gene"
    ]

    assert not any(
        decision.old_stable_id == "ENSXG0001_1" for decision in gene_decisions
    )
    assert any(
        decision.action == "missing" and decision.old_stable_id == "ENSXG0001"
        for decision in gene_decisions
    )
    assert any(
        decision.action == "new" and decision.current_stable_id == "ENSXG9001"
        for decision in gene_decisions
    )


def test_allocator_skips_reserved_ids() -> None:
    reserved_features = {
        "gene": {
            "ENSXG0001": parse_gff3_feature("ENSXG0001"),
            "ENSXG0002": parse_gff3_feature("ENSXG0002"),
        }
    }
    allocator = make_allocator(
        parse_id_range("ENSXG:0001-0003"),
        collect_reserved_ids(reserved_features),
    )

    assert allocator.allocate() == "ENSXG0003"


def parse_gff3_feature(stable_id: str):
    from stable_id_mapping.models import Feature

    return Feature(
        stable_id=stable_id,
        version=1,
        seqid="chr1",
        start=1,
        end=10,
        strand="+",
    )

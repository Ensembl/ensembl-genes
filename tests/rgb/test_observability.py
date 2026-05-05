from __future__ import annotations

import pandas as pd

from ensembl.genes.rgb.observability import (
    audit_evidence_fate,
    audit_expected_presence,
    build_busco_crosswalk,
    build_busco_proxy_calibration,
    build_completeness_profile,
    build_copy_number_audit,
    build_expected_source_profile,
    build_feature_profile,
    build_non_busco_high_confidence_losses,
    build_recommendations,
    build_release_readiness,
    build_same_assembly_structure_audit,
    audit_reference_protein_set,
    build_review_loci,
    build_source_profile,
    expected_tables_from_core_gene_table,
    expected_tables_from_gff3,
    init_expected_template,
    main,
    run_audit,
    validate_inputs,
)


def _genes(rows):
    return pd.DataFrame(rows)


def _transcripts(rows):
    return pd.DataFrame(rows)


def _translations(transcript_ids):
    return pd.DataFrame(
        [{"translation_id": idx + 1, "transcript_id": tid} for idx, tid in enumerate(transcript_ids)]
    )


def test_evidence_fate_flags_coding_orphan_as_rescue_candidate():
    core_genes = _genes(
        [
            {
                "gene_id": 1,
                "stable_id": "CORE1",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "canonical_transcript_id": 101,
                "logic_name": "core",
            }
        ]
    )
    layer_genes = _genes(
        [
            {
                "gene_id": 10,
                "stable_id": "LAYER_REPRESENTED",
                "seq_region_name": "chr1",
                "seq_region_start": 120,
                "seq_region_end": 280,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "canonical_transcript_id": 1001,
                "logic_name": "rnaseq_high",
            },
            {
                "gene_id": 11,
                "stable_id": "LAYER_ORPHAN",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "canonical_transcript_id": 1101,
                "logic_name": "projection_high",
            },
        ]
    )
    core_tx = _transcripts(
        [
            {
                "transcript_id": 101,
                "gene_id": 1,
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "core",
            }
        ]
    )
    layer_tx = _transcripts(
        [
            {
                "transcript_id": 1001,
                "gene_id": 10,
                "seq_region_name": "chr1",
                "seq_region_start": 120,
                "seq_region_end": 280,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "rnaseq_high",
            },
            {
                "transcript_id": 1101,
                "gene_id": 11,
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "projection_high",
            },
        ]
    )

    fate = audit_evidence_fate(
        core_genes=core_genes,
        layer_genes=layer_genes,
        core_tx=core_tx,
        layer_tx=layer_tx,
        core_tr=_translations([101]),
        layer_tr=_translations([1001, 1101]),
    )

    orphan = fate[fate["layer_stable_id"] == "LAYER_ORPHAN"].iloc[0]
    represented = fate[fate["layer_stable_id"] == "LAYER_REPRESENTED"].iloc[0]

    assert orphan["fate_class"] == "orphan_evidence"
    assert orphan["failure_class"] == "no_core_gene_built"
    assert orphan["review_priority"] == "P1"
    assert orphan["suggested_action"] == "try_candidate_rescue"
    assert represented["review_priority"] == "P4"


def test_expected_presence_classifies_clean_missing_and_assembly_limited_cases():
    core_genes = _genes(
        [
            {
                "gene_id": 1,
                "stable_id": "CORE_EXPECTED",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
            }
        ]
    )
    layer_genes = _genes(
        [
            {
                "gene_id": 20,
                "stable_id": "LAYER_MISSING_EXPECTED",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "same_species_projection",
            }
        ]
    )
    expected_genes = pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP_PRESENT",
                "expected_source": "prior_ensembl",
                "reference_stable_id": "OLD1",
                "biotype": "protein_coding",
                "confidence": "high",
            },
            {
                "expected_gene_id": "EXP_MISSING_WITH_EVIDENCE",
                "expected_source": "prior_ensembl",
                "reference_stable_id": "OLD2",
                "biotype": "protein_coding",
                "confidence": "high",
            },
            {
                "expected_gene_id": "EXP_GAP",
                "expected_source": "prior_ensembl",
                "reference_stable_id": "OLD3",
                "biotype": "protein_coding",
                "confidence": "high",
            },
            {
                "expected_gene_id": "EXP_NO_PROJECTION_ROW",
                "expected_source": "prior_ensembl",
                "reference_stable_id": "OLD4",
                "biotype": "protein_coding",
                "confidence": "high",
            },
        ]
    )
    expected_projections = pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP_PRESENT",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "projection_status": "mapped",
                "projection_coverage": 1.0,
            },
            {
                "expected_gene_id": "EXP_MISSING_WITH_EVIDENCE",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "projection_status": "mapped",
                "projection_coverage": 1.0,
            },
            {
                "expected_gene_id": "EXP_GAP",
                "seq_region_name": "chr1",
                "seq_region_start": 900,
                "seq_region_end": 950,
                "seq_region_strand": 1,
                "projection_status": "mapped",
                "projection_coverage": 0.8,
                "assembly_gap_overlap_bp": 20,
            },
        ]
    )

    presence = audit_expected_presence(expected_genes, expected_projections, core_genes, layer_genes)

    by_id = presence.set_index("expected_gene_id")
    assert by_id.loc["EXP_PRESENT", "presence_class"] == "present_clean"
    assert by_id.loc["EXP_PRESENT", "review_priority"] == "P4"
    assert by_id.loc["EXP_MISSING_WITH_EVIDENCE", "presence_class"] == "missing_with_evidence"
    assert by_id.loc["EXP_MISSING_WITH_EVIDENCE", "review_priority"] == "P1"
    assert by_id.loc["EXP_GAP", "presence_class"] == "assembly_limited"
    assert by_id.loc["EXP_GAP", "suggested_action"] == "check_assembly_gap"
    assert by_id.loc["EXP_NO_PROJECTION_ROW", "presence_class"] == "unresolved"
    assert by_id.loc["EXP_NO_PROJECTION_ROW", "failure_class"] == "projection_unmapped"


def test_review_loci_merges_actionable_evidence_and_expected_rows():
    evidence = pd.DataFrame(
        [
            {
                "review_priority": "P1",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "layer_stable_id": "LAYER_ORPHAN",
                "layer_gene_id": 11,
                "failure_class": "no_core_gene_built",
                "suggested_action": "try_candidate_rescue",
                "review_reason": "no_core_gene_built",
            }
        ]
    )
    expected = pd.DataFrame(
        [
            {
                "review_priority": "P2",
                "seq_region_name": "chr2",
                "seq_region_start": 100,
                "seq_region_end": 200,
                "seq_region_strand": -1,
                "expected_gene_id": "EXP_GAP",
                "presence_class": "assembly_limited",
                "suggested_action": "check_assembly_gap",
                "review_reason": "assembly_gap_limited",
            }
        ]
    )

    review = build_review_loci(evidence, expected, top_n=10)

    assert list(review["audit_track"]) == ["evidence_fate", "expected_presence"]
    assert list(review["suggested_action"]) == ["try_candidate_rescue", "check_assembly_gap"]


def test_profiles_capture_sources_busco_copy_number_and_feature_metrics():
    evidence = pd.DataFrame(
        [
            {
                "layer_logic_name": "projection_high",
                "layer_biotype": "protein_coding",
                "layer_coding_transcript_count": 1,
                "fate_class": "orphan_evidence",
                "review_priority": "P1",
                "failure_class": "no_core_gene_built",
                "layer_span_coverage_by_core": 0.0,
            },
            {
                "layer_logic_name": "rnaseq_high",
                "layer_biotype": "protein_coding",
                "layer_coding_transcript_count": 1,
                "fate_class": "represented_exact_or_near_span",
                "review_priority": "P4",
                "failure_class": "represented_span_only",
                "layer_span_coverage_by_core": 1.0,
            },
        ]
    )
    expected = pd.DataFrame(
        [
            {
                "expected_gene_id": "BUSCO_EXPECTED",
                "expected_source": "busco",
                "confidence": "medium",
                "presence_class": "missing_with_evidence",
                "failure_class": "no_core_gene_built",
                "review_priority": "P2",
                "busco_id": "123at456",
                "orthogroup_id": "",
                "expected_copy_number": 1,
                "best_core_stable_id": "",
                "expected_span_coverage_by_core": 0.0,
                "projection_status": "mapped",
                "suggested_action": "try_candidate_rescue",
                "expected_biotype": "protein_coding",
                "best_core_biotype": "",
            },
            {
                "expected_gene_id": "COPY1",
                "expected_source": "prior_ensembl",
                "confidence": "high",
                "presence_class": "present_clean",
                "failure_class": "represented_span_only",
                "review_priority": "P4",
                "busco_id": "",
                "orthogroup_id": "OG_COPY",
                "expected_copy_number": 2,
                "best_core_stable_id": "CORE_COPY1",
                "expected_span_coverage_by_core": 1.0,
                "projection_status": "mapped",
                "suggested_action": "no_action",
                "expected_biotype": "protein_coding",
                "best_core_biotype": "protein_coding",
            },
        ]
    )
    review = pd.DataFrame(
        [
            {
                "review_priority": "P1",
                "audit_track": "evidence_fate",
                "class": "no_core_gene_built",
            },
            {
                "review_priority": "P2",
                "audit_track": "expected_presence",
                "class": "missing_with_evidence",
            },
        ]
    )

    source_profile = build_source_profile(evidence)
    expected_source_profile = build_expected_source_profile(expected)
    busco_crosswalk = build_busco_crosswalk(expected)
    copy_number = build_copy_number_audit(expected)
    feature_profile = build_feature_profile(evidence, expected, review, copy_number, busco_crosswalk)
    recommendations = build_recommendations(
        evidence,
        expected,
        copy_number,
        build_busco_proxy_calibration(expected),
        build_non_busco_high_confidence_losses(expected),
    )

    assert source_profile.iloc[0]["layer_logic_name"] == "projection_high"
    assert expected_source_profile["n_expected_genes"].sum() == 2
    assert busco_crosswalk.iloc[0]["busco_id"] == "123at456"
    copy_row = copy_number[copy_number["copy_group_id"] == "OG_COPY"].iloc[0]
    assert copy_row["copy_number_class"] == "collapsed_copy_number"
    metric_names = set(feature_profile["metric_name"])
    assert "orphan_layer_model_count" in metric_names
    assert "busco_p1_or_p2_issue_count" in metric_names
    assert "candidate_rescue_from_layer" in set(recommendations["recommendation_id"])


def test_completeness_profile_compares_busco_with_non_busco_expected_genes():
    expected = pd.DataFrame(
        [
            {
                "expected_gene_id": "BUSCO1",
                "expected_source": "busco",
                "reference_stable_id": "BUSCO1",
                "symbol": "",
                "expected_biotype": "protein_coding",
                "orthogroup_id": "",
                "busco_id": "123at456",
                "expected_copy_number": 1,
                "confidence": "medium",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "projection_status": "mapped",
                "presence_class": "present_clean",
                "failure_class": "represented_span_only",
                "review_priority": "P4",
                "best_core_stable_id": "CORE1",
                "best_layer_stable_id": "",
                "best_layer_logic_name": "",
                "expected_span_coverage_by_core": 1.0,
                "suggested_action": "no_action",
                "review_reason": "represented",
            },
            {
                "expected_gene_id": "NONBUSCO1",
                "expected_source": "prior_ensembl",
                "reference_stable_id": "OLD1",
                "symbol": "GENE1",
                "expected_biotype": "protein_coding",
                "orthogroup_id": "OG1",
                "busco_id": "",
                "expected_copy_number": 1,
                "confidence": "high",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "projection_status": "mapped",
                "presence_class": "missing_with_evidence",
                "failure_class": "no_core_gene_built",
                "review_priority": "P1",
                "best_core_stable_id": "",
                "best_layer_stable_id": "LAYER1",
                "best_layer_logic_name": "projection_high",
                "expected_span_coverage_by_core": 0.0,
                "suggested_action": "try_candidate_rescue",
                "review_reason": "missing",
            },
        ]
    )

    completeness = build_completeness_profile(expected)
    losses = build_non_busco_high_confidence_losses(expected)
    calibration = build_busco_proxy_calibration(expected)

    panels = set(completeness["panel_id"])
    assert "busco_linked" in panels
    assert "high_confidence_non_busco" in panels
    assert losses.iloc[0]["expected_gene_id"] == "NONBUSCO1"
    delta = calibration[
        calibration["comparison"] == "non_busco_minus_busco_actionable_loss"
    ].iloc[0]
    assert delta["actionable_loss_fraction"] > 0


def test_expected_tables_from_core_gene_table_can_generate_same_coordinate_inputs():
    genes = pd.DataFrame(
        [
            {
                "gene_id": 1,
                "stable_id": "OLD1",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
            },
            {
                "gene_id": 2,
                "stable_id": "OLD2",
                "seq_region_name": "old_unplaced",
                "seq_region_start": 20,
                "seq_region_end": 90,
                "seq_region_strand": -1,
                "biotype": "lncRNA",
            },
        ]
    )

    expected, projections = expected_tables_from_core_gene_table(
        genes,
        expected_source="prior_ensembl",
        confidence="high",
        projection_mode="same_coordinates",
        target_seq_regions=["chr1"],
    )

    assert list(expected["expected_gene_id"]) == ["OLD1", "OLD2"]
    assert list(expected["expected_source"]) == ["prior_ensembl", "prior_ensembl"]
    assert projections.set_index("expected_gene_id").loc["OLD1", "projection_status"] == "mapped"
    assert projections.set_index("expected_gene_id").loc["OLD2", "projection_status"] == "unmapped"


def test_expected_tables_from_gff3_can_generate_expected_genes_and_projections(tmp_path):
    gff3 = tmp_path / "expected.gff3"
    gff3.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tRefSeq\tgene\t100\t300\t.\t+\t.\tID=gene1;Name=GENE1;biotype=protein_coding;orthogroup_id=OG1",
                "chr1\tRefSeq\tmRNA\t100\t300\t.\t+\t.\tID=tx1;Parent=gene1",
                "old_scaffold\tRefSeq\tgene\t500\t700\t.\t-\t.\tID=gene2;Name=GENE2;gene_biotype=lncRNA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    expected, projections = expected_tables_from_gff3(
        str(gff3),
        expected_source="refseq_gff3",
        confidence="high",
        projection_mode="same_coordinates",
        target_seq_regions=["chr1"],
    )

    by_id = expected.set_index("expected_gene_id")
    assert by_id.loc["gene1", "expected_source"] == "refseq_gff3"
    assert by_id.loc["gene1", "symbol"] == "GENE1"
    assert by_id.loc["gene1", "orthogroup_id"] == "OG1"
    assert projections.set_index("expected_gene_id").loc["gene1", "projection_status"] == "mapped"
    assert projections.set_index("expected_gene_id").loc["gene2", "projection_status"] == "unmapped"


def test_same_assembly_structure_audit_classifies_gffcompare_rows():
    tmap = pd.DataFrame(
        [
            {
                "reference_gene_id": "REFG1",
                "reference_transcript_id": "REFT1",
                "class_code": "=",
                "query_gene_id": "QG1",
                "query_transcript_id": "QT1",
                "num_exons": 3,
                "length": 1000,
                "reference_match_length": 1000,
            },
            {
                "reference_gene_id": "REFG2",
                "reference_transcript_id": "REFT2",
                "class_code": "j",
                "query_gene_id": "QG2",
                "query_transcript_id": "QT2",
                "num_exons": 4,
                "length": 900,
                "reference_match_length": 700,
            },
        ]
    )

    audit = build_same_assembly_structure_audit(tmap)

    by_tx = audit.set_index("query_transcript_id")
    assert by_tx.loc["QT1", "structure_match_class"] == "exact_intron_chain"
    assert by_tx.loc["QT1", "review_priority"] == "P4"
    assert by_tx.loc["QT2", "structure_match_class"] == "splice_compatible_partial"
    assert by_tx.loc["QT2", "review_priority"] == "P2"


def test_reference_protein_audit_finds_protein_supported_missing_core_gene():
    expected_proteins = pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP_PROT1",
                "query_protein_id": "PROT1",
                "expected_source": "prior_ensembl_protein",
                "reference_stable_id": "OLD1",
                "biotype": "protein_coding",
                "confidence": "high",
            },
            {
                "expected_gene_id": "EXP_PROT2",
                "query_protein_id": "PROT2",
                "expected_source": "prior_ensembl_protein",
                "reference_stable_id": "OLD2",
                "biotype": "protein_coding",
                "confidence": "high",
            },
        ]
    )
    protein_hits = pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP_PROT1",
                "query_protein_id": "PROT1",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "percent_identity": 0.95,
                "query_coverage": 0.98,
                "target_coverage": 0.90,
                "hit_rank": 1,
            },
            {
                "expected_gene_id": "EXP_PROT2",
                "query_protein_id": "PROT2",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "percent_identity": 0.90,
                "query_coverage": 0.95,
                "target_coverage": 0.85,
                "hit_rank": 1,
            },
        ]
    )
    core_genes = _genes(
        [
            {
                "gene_id": 1,
                "stable_id": "CORE1",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
            }
        ]
    )

    audit = audit_reference_protein_set(
        expected_proteins,
        protein_hits,
        core_genes,
        min_identity=0.30,
        min_query_coverage=0.70,
    )

    by_id = audit.set_index("expected_gene_id")
    assert by_id.loc["EXP_PROT1", "protein_hit_class"] == "protein_supported_built"
    assert by_id.loc["EXP_PROT2", "protein_hit_class"] == "protein_supported_no_core_gene"
    assert by_id.loc["EXP_PROT2", "review_priority"] == "P1"


def test_release_readiness_fails_bad_busco_copy_and_protein_rescue_signals():
    expected = pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP1",
                "confidence": "high",
                "presence_class": "missing_with_evidence",
            }
        ]
    )
    copy_number = pd.DataFrame(
        [
            {
                "copy_group_id": "OG1",
                "copy_number_class": "collapsed_copy_number",
            }
        ]
    )
    reference_protein = pd.DataFrame(
        [
            {
                "protein_hit_class": "protein_supported_no_core_gene",
            }
        ]
    )
    completeness = pd.DataFrame(
        [
            {
                "panel_id": "high_confidence",
                "actionable_loss_fraction": 1.0,
            }
        ]
    )

    readiness = build_release_readiness(
        expected,
        copy_number,
        reference_protein,
        pd.DataFrame(),
        completeness,
        busco_complete_percent=80.0,
        busco_floor_percent=95.0,
        max_high_confidence_actionable_loss_fraction=0.0,
        max_copy_number_issue_fraction=0.0,
        require_expected_genes=True,
    )

    failed = readiness[readiness["status"] == "FAIL"]
    assert "busco_floor" in set(failed["gate_id"])
    assert "copy_number_regression" in set(failed["gate_id"])
    assert "protein_supported_missing_core" in set(failed["gate_id"])


def test_run_audit_writes_framework_outputs(tmp_path):
    run_id = "demo"
    extract_dir = tmp_path / run_id / "extract"
    extract_dir.mkdir(parents=True)

    core_genes = _genes(
        [
            {
                "gene_id": 1,
                "stable_id": "CORE1",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "canonical_transcript_id": 101,
                "logic_name": "core",
            }
        ]
    )
    layer_genes = _genes(
        [
            {
                "gene_id": 10,
                "stable_id": "LAYER_ORPHAN",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "canonical_transcript_id": 1001,
                "logic_name": "projection_high",
            }
        ]
    )
    core_tx = _transcripts(
        [
            {
                "transcript_id": 101,
                "gene_id": 1,
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "core",
            }
        ]
    )
    layer_tx = _transcripts(
        [
            {
                "transcript_id": 1001,
                "gene_id": 10,
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "projection_high",
            }
        ]
    )
    expected_genes = pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP1",
                "expected_source": "prior_ensembl",
                "reference_stable_id": "OLD1",
                "biotype": "protein_coding",
                "confidence": "high",
                "expected_copy_number": 1,
            }
        ]
    )
    expected_projections = pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP1",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "projection_status": "mapped",
                "projection_coverage": 1.0,
            }
        ]
    )

    core_genes.to_csv(extract_dir / "core_genes.tsv", sep="\t", index=False)
    layer_genes.to_csv(extract_dir / "layer_genes.tsv", sep="\t", index=False)
    core_tx.to_csv(extract_dir / "core_transcripts.tsv", sep="\t", index=False)
    layer_tx.to_csv(extract_dir / "layer_transcripts.tsv", sep="\t", index=False)
    _translations([101]).to_csv(extract_dir / "core_translations.tsv", sep="\t", index=False)
    _translations([1001]).to_csv(extract_dir / "layer_translations.tsv", sep="\t", index=False)
    expected_genes_path = tmp_path / "expected_genes.tsv"
    expected_projections_path = tmp_path / "expected_projections.tsv"
    expected_genes.to_csv(expected_genes_path, sep="\t", index=False)
    expected_projections.to_csv(expected_projections_path, sep="\t", index=False)

    args = type(
        "Args",
        (),
        {
            "output_dir": str(tmp_path),
            "run_id": run_id,
            "format": "tsv",
            "expected_genes": str(expected_genes_path),
            "expected_projections": str(expected_projections_path),
            "locus_gap_bp": 5000,
            "busco_complete_percent": 80.0,
            "busco_floor_percent": 95.0,
            "max_high_confidence_actionable_loss_fraction": 0.0,
            "max_copy_number_issue_fraction": 0.0,
            "require_expected_genes": True,
            "top_n": 100,
            "bed_pad_bp": 0,
        },
    )()

    assert run_audit(args) == 0

    obs_dir = tmp_path / run_id / "observability"
    expected_outputs = [
        "audit_loci.tsv",
        "evidence_fate.tsv",
        "expected_gene_presence.tsv",
        "source_profile.tsv",
        "expected_source_profile.tsv",
        "busco_expected_crosswalk.tsv",
        "completeness_profile.tsv",
        "non_busco_high_confidence_losses.tsv",
        "busco_proxy_calibration.tsv",
        "copy_number_audit.tsv",
        "biotype_transition.tsv",
        "recommendations.tsv",
        "release_readiness.tsv",
        "feature_profile.tsv",
        "failure_mode_summary.tsv",
        "review_loci.tsv",
        "review_loci.bed",
        "actionable_summary.md",
    ]
    for name in expected_outputs:
        assert (obs_dir / name).exists()
    assert (obs_dir / "review_beds" / "missing_with_evidence.bed").exists()
    assert (obs_dir / "completeness_beds" / "non_busco_high_confidence_losses.bed").exists()


def test_template_and_validation_helpers(tmp_path):
    template_args = type("Args", (), {"out_dir": str(tmp_path / "templates")})()
    assert init_expected_template(template_args) == 0
    assert (tmp_path / "templates" / "expected_genes.tsv").exists()
    assert (tmp_path / "templates" / "expected_projections.tsv").exists()
    assert (tmp_path / "templates" / "expected_proteins.tsv").exists()
    assert (tmp_path / "templates" / "reference_protein_hits.tsv").exists()

    run_dir = tmp_path / "run"
    extract_dir = run_dir / "extract"
    extract_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "gene_id": 1,
                "seq_region_name": "chr1",
                "seq_region_start": 1,
                "seq_region_end": 10,
                "seq_region_strand": 1,
            }
        ]
    ).to_csv(extract_dir / "core_genes.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "gene_id": 2,
                "seq_region_name": "chr1",
                "seq_region_start": 20,
                "seq_region_end": 30,
                "seq_region_strand": 1,
            }
        ]
    ).to_csv(extract_dir / "layer_genes.tsv", sep="\t", index=False)

    validate_args = type(
        "Args",
        (),
        {
            "run_dir": str(run_dir),
            "output_dir": None,
            "run_id": None,
            "format": "tsv",
            "expected_genes": str(tmp_path / "templates" / "expected_genes.tsv"),
            "expected_projections": str(tmp_path / "templates" / "expected_projections.tsv"),
            "expected_proteins": str(tmp_path / "templates" / "expected_proteins.tsv"),
            "reference_protein_hits": str(tmp_path / "templates" / "reference_protein_hits.tsv"),
            "gffcompare_tmap": None,
        },
    )()
    assert validate_inputs(validate_args) == 0


def test_cli_audit_smoke_with_expected_gff3_protein_hits_and_release_gates(tmp_path):
    run_dir = tmp_path / "cli_demo"
    extract_dir = run_dir / "extract"
    extract_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "gene_id": 1,
                "stable_id": "CORE1",
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 300,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
            }
        ]
    ).to_csv(extract_dir / "core_genes.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "gene_id": 10,
                "stable_id": "LAYER_EXPECTED",
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 700,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "protein_projection",
            }
        ]
    ).to_csv(extract_dir / "layer_genes.tsv", sep="\t", index=False)

    expected_gff3 = tmp_path / "expected.gff3"
    expected_gff3.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t500\t700\t.\t+\t.\tID=EXP_GFF3;Name=GFF3_EXPECTED;biotype=protein_coding\n",
        encoding="utf-8",
    )
    expected_proteins = tmp_path / "expected_proteins.tsv"
    pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP_PROT",
                "query_protein_id": "PROT1",
                "expected_source": "prior_protein",
                "reference_stable_id": "OLD_PROT",
                "biotype": "protein_coding",
                "confidence": "high",
            }
        ]
    ).to_csv(expected_proteins, sep="\t", index=False)
    protein_hits = tmp_path / "reference_protein_hits.tsv"
    pd.DataFrame(
        [
            {
                "expected_gene_id": "EXP_PROT",
                "query_protein_id": "PROT1",
                "seq_region_name": "chr1",
                "seq_region_start": 900,
                "seq_region_end": 1100,
                "seq_region_strand": 1,
                "percent_identity": 0.95,
                "query_coverage": 0.90,
                "hit_rank": 1,
            }
        ]
    ).to_csv(protein_hits, sep="\t", index=False)

    assert (
        main(
            [
                "audit",
                "--run_dir",
                str(run_dir),
                "--expected_gff3",
                str(expected_gff3),
                "--expected_gff3_projection_mode",
                "same_coordinates",
                "--expected_proteins",
                str(expected_proteins),
                "--reference_protein_hits",
                str(protein_hits),
                "--busco_complete_percent",
                "80.0",
                "--require_expected_genes",
                "--top_n",
                "20",
            ]
        )
        == 0
    )

    obs_dir = run_dir / "observability"
    assert (run_dir / "expected" / "expected_genes.tsv").exists()
    assert (run_dir / "expected" / "expected_projections.tsv").exists()
    readiness = pd.read_csv(obs_dir / "release_readiness.tsv", sep="\t")
    assert "busco_floor" in set(readiness[readiness["status"] == "FAIL"]["gate_id"])
    protein_audit = pd.read_csv(obs_dir / "reference_protein_audit.tsv", sep="\t")
    assert protein_audit.iloc[0]["protein_hit_class"] == "protein_supported_no_core_gene"

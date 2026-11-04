"""Tests for expression-track handover generation and validation."""

# Test doubles and small data builders are deliberately compact.
# pylint: disable=missing-class-docstring,missing-function-docstring
# pylint: disable=too-few-public-methods,use-implicit-booleaness-not-comparison

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import ensembl.genes.transcriptomic_data.create_expression_track_handover as handover_module
from ensembl.genes.transcriptomic_data.create_expression_track_handover import (
    AssemblyResolution,
    HandoverError,
    MetadataDbAssemblyResolver,
    build_handover,
    main,
    parse_args,
)
from ensembl.genes.transcriptomic_data.validate_expression_track_handover import (
    validate_records,
)


class FakeAssemblyResolver:
    def resolve(self, gca_accession: str) -> AssemblyResolution:
        return AssemblyResolution(
            gca_accession=gca_accession,
            assembly_uuid="00000000-0000-0000-0000-000000000001",
            linked_genomes=({"genome_uuid": "00000000-0000-0000-0000-000000000001"},),
            evidence={"source": "test"},
        )


class FakeConnection:
    def __init__(self, assembly_rows, linked_rows):
        self.assembly_rows = assembly_rows
        self.linked_rows = linked_rows

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, _params):
        rows = self.assembly_rows if "dataset_name" in str(query) else self.linked_rows
        return type("Result", (), {"fetchall": lambda _self: rows})()


class FakeEngine:
    def __init__(self, assembly_rows, linked_rows):
        self.connection = FakeConnection(assembly_rows, linked_rows)

    def connect(self):
        return self.connection


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_inputs(tmp_path: Path) -> dict[str, Path]:
    selected = tmp_path / "selected.tsv"
    write_tsv(
        selected,
        [
            {
                "run_accession": "RUN1",
                "sample_id": "SAMPLE1",
                "gca_accession": "GCA_000003025.6",
            },
            {
                "run_accession": "RUN2",
                "sample_id": "SAMPLE2",
                "gca_accession": "GCA_000003025.6",
            },
        ],
    )

    run_metadata = tmp_path / "runs.json"
    write_json(
        run_metadata,
        [
            {
                "run_accession": "RUN1",
                "sample_accession": "SAMPLE1",
                "study_accession": "STUDY1",
            },
            {
                "run_accession": "RUN2",
                "sample_accession": "SAMPLE2",
                "study_accession": "STUDY1",
            },
        ],
    )
    biosamples = tmp_path / "biosamples.json"
    write_json(
        biosamples,
        [
            {
                "sample_accession": "SAMPLE1",
                "attributes": [
                    {"name": "organism", "value": "Sus scrofa"},
                    {"name": "tissue", "value": "brain"},
                    {"name": "genotype", "value": "wild type"},
                    {"name": "sex", "value": "female"},
                ],
            },
            {
                "sample_accession": "SAMPLE2",
                "attributes": [
                    {"name": "organism", "value": "Sus scrofa"},
                    {"name": "tissue", "value": "brain"},
                    {"name": "genotype", "value": "wild type"},
                    {"name": "sex", "value": "female"},
                ],
            },
        ],
    )

    artifact = tmp_path / "track.bin"
    artifact.write_bytes(b"test track\n")
    checksum = hashlib.sha256(artifact.read_bytes()).hexdigest()
    run_rows: list[dict[str, object]] = [
        {
            "entity_id": run_id,
            "track_type": product,
            "track_file": str(artifact),
            "status": "complete",
            "sha256": checksum,
            "gca_accession": "GCA_000003025.6",
        }
        for run_id in ("RUN1", "RUN2")
        for product in (
            "coverage_forward",
            "coverage_reverse",
            "gene_models",
            "splice_junctions",
        )
    ]
    run_manifest = tmp_path / "run_manifest.tsv"
    write_tsv(run_manifest, run_rows)

    aggregate_manifest = tmp_path / "aggregate_manifest.tsv"
    aggregate_rows: list[dict[str, object]] = [
        {
            "entity_id": "AGG_BRAIN",
            "aggregate_label": "brain_cns",
            "track_type": product,
            "track_file": str(artifact),
            "status": "complete",
            "sha256": checksum,
            "source_entity_ids": "RUN1;RUN2",
            "source_count": "2",
            "gca_accession": "GCA_000003025.6",
        }
        for product in ("coverage_forward", "coverage_reverse", "splice_junctions")
    ]
    write_tsv(
        aggregate_manifest,
        aggregate_rows,
    )
    return {
        "selected": selected,
        "run_metadata": run_metadata,
        "biosamples": biosamples,
        "run_manifest": run_manifest,
        "aggregate_manifest": aggregate_manifest,
    }


def test_builds_expected_run_and_tissue_records(tmp_path: Path) -> None:
    files = make_inputs(tmp_path)
    args = parse_args(
        [
            "--selected-runs",
            str(files["selected"]),
            "--run-metadata",
            str(files["run_metadata"]),
            "--biosamples",
            str(files["biosamples"]),
            "--run-manifest",
            str(files["run_manifest"]),
            "--aggregate-manifest",
            str(files["aggregate_manifest"]),
            "--mode",
            "draft",
            "--out",
            str(tmp_path / "handover.json"),
            "--normalized-out",
            str(tmp_path / "normalized.tsv"),
            "--report-out",
            str(tmp_path / "report.json"),
        ]
    )

    output, report, _ = build_handover(args, FakeAssemblyResolver())

    assert validate_records(output) == []
    assert len(output) == 11
    assert report["run_record_count"] == 8
    assert report["tissue_record_count"] == 3
    assert len({record["name"] for record in output}) == 11
    assert all(
        record["assembly_uuid"] == "00000000-0000-0000-0000-000000000001"
        for record in output
    )
    aggregate_records = [record for record in output if "_tissue_" in record["name"]]
    assert all(
        "parent_sample_id" not in record["track_metadata_payload"]
        for record in aggregate_records
    )
    run_record = next(
        record for record in output if record["name"].endswith("RUN1_coverage_forward")
    )
    payload = run_record["track_metadata_payload"]
    assert payload["schema_version"] == "transcriptomic-selection-0.3.1"
    assert payload["metadata"]["tissue"]["value"] == "Brain"
    assert payload["metadata"]["genotype"]["value"] == "Wild Type"
    assert payload["metadata"]["sex"]["value"] == "Female"
    assert payload["metadata"]["genotype"]["provenance"]["source_value"] == "wild type"
    assert payload["metadata"]["sex"]["provenance"]["source_value"] == "female"
    assert payload["metadata"]["tissue"]["provenance"] == {
        "source": "ENA BioSample",
        "source_field": "tissue",
        "source_id": "SAMPLE1",
        "source_value": "brain",
        "normalization": "explicit_display_label_rule",
    }
    assert "biosample" not in payload["metadata"]
    assert "hierarchy_level" not in json.dumps(payload)
    assert "parent_sample_id" not in payload
    assert "track_products" not in payload
    assert "source_run_count" not in payload["provenance"]
    assert "source_bam_submissions" not in payload["provenance"]
    assert payload["provenance"]["biosample_accession"] == "SAMPLE1"
    assert payload["target_genome"]["assembly_uuid"] == (
        "00000000-0000-0000-0000-000000000001"
    )
    assert payload["selection_hierarchy"] == ["tissue", "biosample", "run"]
    aggregate_record = next(
        record
        for record in output
        if "tissue_GCA_000003025_6_AGG_BRAIN_coverage_forward" in record["name"]
    )
    aggregate_tissue = aggregate_record["track_metadata_payload"]["metadata"]["tissue"]
    assert aggregate_tissue["value"] == "Brain CNS"
    assert aggregate_tissue["provenance"]["source"] == "aggregate_manifest"
    assert aggregate_record["track_metadata_payload"]["display_label"] == "Brain CNS"
    assert "_" not in aggregate_record["track_metadata_payload"]["display_label"]
    assert {record["name"].rsplit("_", 1)[-1] for record in output} == {
        "forward",
        "reverse",
        "models",
        "junctions",
    }


def test_metadata_db_resolver_returns_uuid_and_linked_genomes() -> None:
    resolver = MetadataDbAssemblyResolver(
        metadata_uri=(
            "mysql+pymysql://ensro:secret@mysql-ens-production-1:4721/"
            "ensembl_genome_metadata?token=also-secret"
        ),
        engine=FakeEngine(
            [("00000000-0000-0000-0000-000000000001",)],
            [
                ("00000000-0000-0000-0000-000000000001", "GCA_000003025.6"),
                ("linked-uuid", "GCA_000003025.6"),
            ],
        ),
    )

    result = resolver.resolve("GCA_000003025.6")

    assert result.assembly_uuid == "00000000-0000-0000-0000-000000000001"
    assert {row["genome_uuid"] for row in result.linked_genomes} == {
        "00000000-0000-0000-0000-000000000001",
        "linked-uuid",
    }
    assert result.evidence["metadata_uri"] == (
        "mysql+pymysql://ensro:***@mysql-ens-production-1:4721/"
        "ensembl_genome_metadata"
    )
    assert "secret" not in json.dumps(result.evidence)


def test_draft_can_emit_placeholder_track_sources(tmp_path: Path) -> None:
    files = make_inputs(tmp_path)
    args = parse_args(
        [
            "--selected-runs",
            str(files["selected"]),
            "--run-metadata",
            str(files["run_metadata"]),
            "--biosamples",
            str(files["biosamples"]),
            "--run-manifest",
            str(files["run_manifest"]),
            "--aggregate-manifest",
            str(files["aggregate_manifest"]),
            "--mode",
            "draft",
            "--dummy-track-source",
            "--out",
            str(tmp_path / "handover.json"),
            "--normalized-out",
            str(tmp_path / "normalized.tsv"),
            "--report-out",
            str(tmp_path / "report.json"),
        ]
    )

    output, _, _ = build_handover(args, FakeAssemblyResolver())

    run_record = next(
        record
        for record in output
        if record["track_metadata_payload"]["selection_level"] == "run"
    )
    assert run_record["track_source"] == {
        "name": "ENA",
        "details": "Run RUN1",
        "url": "https://www.ebi.ac.uk/ena/browser/view/RUN1",
    }
    aggregate_record = next(
        record
        for record in output
        if record["track_metadata_payload"]["selection_level"] == "tissue"
    )
    assert aggregate_record["track_source"]["name"] == "ENA project placeholder"
    assert (
        "Source BAM provenance for aggregate"
        in aggregate_record["track_source"]["details"]
    )
    assert aggregate_record["track_source"]["url"].startswith(
        "https://example.org/ena-project/"
    )


def test_validator_reports_nested_assembly_uuid_mismatch(tmp_path: Path) -> None:
    files = make_inputs(tmp_path)
    args = parse_args(
        [
            "--selected-runs",
            str(files["selected"]),
            "--run-metadata",
            str(files["run_metadata"]),
            "--biosamples",
            str(files["biosamples"]),
            "--run-manifest",
            str(files["run_manifest"]),
            "--aggregate-manifest",
            str(files["aggregate_manifest"]),
            "--mode",
            "draft",
            "--out",
            str(tmp_path / "handover.json"),
            "--normalized-out",
            str(tmp_path / "normalized.tsv"),
            "--report-out",
            str(tmp_path / "report.json"),
        ]
    )
    output, _, _ = build_handover(args, FakeAssemblyResolver())
    output[0]["track_metadata_payload"]["target_genome"][
        "assembly_uuid"
    ] = "00000000-0000-0000-0000-000000000002"

    errors = validate_records(output)

    assert any("inconsistent assembly UUIDs" in error for error in errors)


def test_metadata_db_resolver_rejects_ambiguous_uuid() -> None:
    resolver = MetadataDbAssemblyResolver(
        engine=FakeEngine([("uuid-1",), ("uuid-2",)], []),
    )

    try:
        resolver.resolve("GCA_000003025.6")
    except HandoverError as exc:
        assert "exactly one" in str(exc)
    else:
        raise AssertionError("ambiguous assembly lookup should fail")


def test_validator_rejects_outer_and_nested_inconsistencies(tmp_path: Path) -> None:
    files = make_inputs(tmp_path)
    args = parse_args(
        [
            "--selected-runs",
            str(files["selected"]),
            "--run-metadata",
            str(files["run_metadata"]),
            "--biosamples",
            str(files["biosamples"]),
            "--run-manifest",
            str(files["run_manifest"]),
            "--aggregate-manifest",
            str(files["aggregate_manifest"]),
            "--mode",
            "draft",
            "--out",
            str(tmp_path / "handover.json"),
            "--normalized-out",
            str(tmp_path / "normalized.tsv"),
            "--report-out",
            str(tmp_path / "report.json"),
        ]
    )
    output, _, _ = build_handover(args, FakeAssemblyResolver())

    wrong_name = dict(output[0])
    wrong_name["name"] = "different-name"
    assert any(
        "does not match track_metadata_payload.selection_id" in error
        for error in validate_records([wrong_name])
    )

    wrong_source = dict(output[0])
    wrong_source["track_files"] = {"summary": "different-track.bw"}
    assert any(
        "dataset_source.name must match track_files.summary" in error
        for error in validate_records([wrong_source])
    )


def test_final_failure_writes_report_without_publishing_outputs(
    tmp_path: Path, monkeypatch
) -> None:
    files = make_inputs(tmp_path)
    source_bams = tmp_path / "source_bams.json"
    write_json(source_bams, [])
    output_path = tmp_path / "handover.json"
    normalized_path = tmp_path / "normalized.tsv"
    report_path = tmp_path / "report.json"
    monkeypatch.setattr(
        handover_module,
        "MetadataDbAssemblyResolver",
        lambda _uri: FakeAssemblyResolver(),
    )

    exit_code = main(
        [
            "--selected-runs",
            str(files["selected"]),
            "--run-metadata",
            str(files["run_metadata"]),
            "--biosamples",
            str(files["biosamples"]),
            "--run-manifest",
            str(files["run_manifest"]),
            "--aggregate-manifest",
            str(files["aggregate_manifest"]),
            "--source-bams",
            str(source_bams),
            "--mode",
            "final",
            "--out",
            str(output_path),
            "--normalized-out",
            str(normalized_path),
            "--report-out",
            str(report_path),
        ]
    )

    assert exit_code == 1
    assert not output_path.exists()
    assert not normalized_path.exists()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert any(
        diagnostic["code"] == "missing_source_bam_submission"
        for diagnostic in report["diagnostics"]
    )


def test_generator_rejects_invalid_generated_payload(
    tmp_path: Path, monkeypatch
) -> None:
    files = make_inputs(tmp_path)
    output_path = tmp_path / "handover.json"
    normalized_path = tmp_path / "normalized.tsv"
    report_path = tmp_path / "report.json"
    monkeypatch.setattr(
        handover_module,
        "MetadataDbAssemblyResolver",
        lambda _uri: FakeAssemblyResolver(),
    )

    exit_code = main(
        [
            "--selected-runs",
            str(files["selected"]),
            "--run-metadata",
            str(files["run_metadata"]),
            "--biosamples",
            str(files["biosamples"]),
            "--run-manifest",
            str(files["run_manifest"]),
            "--aggregate-manifest",
            str(files["aggregate_manifest"]),
            "--selection-hierarchy",
            "",
            "--mode",
            "draft",
            "--out",
            str(output_path),
            "--normalized-out",
            str(normalized_path),
            "--report-out",
            str(report_path),
        ]
    )

    assert exit_code == 1
    assert not output_path.exists()
    assert not normalized_path.exists()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert any(
        diagnostic["code"] == "generated_record_invalid"
        for diagnostic in report["diagnostics"]
    )

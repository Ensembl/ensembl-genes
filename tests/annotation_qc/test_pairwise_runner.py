"""End-to-end runs of annotation-qc pairwise-compare on small fixtures."""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess

import pytest

from ensembl.genes.annotation_qc.runners import pairwise_compare

from .conftest import gene_records

MULTI = [(100, 200), (300, 400), (500, 600)]
OUTPUTS = [
    "comparison_summary.json",
    "comparison_summary.tsv",
    "comparison_details.tsv",
    "consensus_transcript_labels.tsv",
    "reference_filter_audit.json",
    "reference_filter_audit.tsv",
    "comparison_manifest.json",
]


def _run(argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    pairwise_compare.register(subparsers)
    args = parser.parse_args(["pairwise-compare", *argv])
    return args.func(args)


@pytest.fixture
def inputs(write_gff3, write_text):
    ref = write_gff3(
        "ref.gff3.gz",
        gene_records("Pf_01", "R1", "+", {"R1.1": MULTI})
        + gene_records(
            "Pf_01",
            "R2",
            "+",
            {"R2.1": [(2000, 2500)]},
            biotype="lncRNA",
            with_cds=False,
        )
        + gene_records("Pf_MIT", "R3", "+", {"R3.1": [(10, 90)]}),
        compress=True,
    )
    query = write_gff3(
        "consensus.gff3",
        gene_records("1", "Q1", "+", {"Q1.1": MULTI})
        + gene_records("1", "Q2", "-", {"Q2.1": [(5000, 5100)]}),
    )
    genome = write_text("genome.fa", ">1\n" + "A" * 6000 + "\n")
    seqmap = write_text("map.tsv", "from_seqname\tto_seqname\nPf_01\t1\n")
    return {"ref": ref, "query": query, "genome": genome, "map": seqmap}


def test_outputs_manifest_and_coordinates(inputs, tmp_path):
    outdir = tmp_path / "out"
    summary = _run(
        [
            "--query",
            inputs["query"],
            "--reference",
            inputs["ref"],
            "--genome",
            inputs["genome"],
            "--seqname-map",
            inputs["map"],
            "--genome-mismatch",
            "exclude",
            "--evaluation-mode",
            "protein_coding",
            "--outdir",
            str(outdir),
        ]
    )
    for name in OUTPUTS:
        assert (outdir / name).is_file(), name

    # Pf_MIT is absent from the genome and excluded; lncRNA R2 removed by the mode.
    assert summary["total_reference_genes"] == 1
    assert summary["reference_classification"] == {"Exact_Match": 1}
    assert summary["consensus_classification"] == {"Exact_Match": 1, "Novel": 1}

    with open(outdir / "comparison_details.tsv") as handle:
        rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert (rows["R1"]["chrom"], rows["R1"]["start"], rows["R1"]["end"]) == (
        "1",
        "100",
        "600",
    )

    manifest = json.loads((outdir / "comparison_manifest.json").read_text())
    expected = hashlib.sha256(open(inputs["query"], "rb").read()).hexdigest()
    assert manifest["inputs"]["query"]["sha256"] == expected
    assert manifest["inputs"]["reference"]["path"] == os.path.realpath(inputs["ref"])
    assert manifest["seqname_checks"]["excluded"] == {"reference": {"Pf_MIT": 1}}
    assert manifest["options"]["evaluation_mode"] == "protein_coding"
    assert manifest["code_version"]["pyranges1"]

    audit = json.loads((outdir / "reference_filter_audit.json").read_text())
    assert (audit["pre_filter_genes"], audit["post_filter_genes"]) == (2, 1)
    assert audit["excluded_seqnames"] == {"Pf_MIT": 1}


def test_no_shared_seqnames_fails_before_writing(inputs, tmp_path):
    outdir = tmp_path / "out"
    with pytest.raises(SystemExit, match="share no sequence names"):
        _run(
            [
                "--query",
                inputs["query"],
                "--reference",
                inputs["ref"],
                "--outdir",
                str(outdir),
            ]
        )
    assert not outdir.exists()


def test_sequence_missing_from_genome_fails_by_default(inputs, tmp_path):
    with pytest.raises(SystemExit, match="Pf_MIT"):
        _run(
            [
                "--query",
                inputs["query"],
                "--reference",
                inputs["ref"],
                "--genome",
                inputs["genome"],
                "--seqname-map",
                inputs["map"],
                "--outdir",
                str(tmp_path / "out"),
            ]
        )


def test_warn_keeps_legacy_denominator(inputs, tmp_path):
    summary = _run(
        [
            "--query",
            inputs["query"],
            "--reference",
            inputs["ref"],
            "--genome",
            inputs["genome"],
            "--seqname-map",
            inputs["map"],
            "--genome-mismatch",
            "warn",
            "--evaluation-mode",
            "protein_coding",
            "--outdir",
            str(tmp_path / "out"),
        ]
    )
    assert summary["reference_classification"] == {"Exact_Match": 1, "Missed": 1}


def test_regions_and_evidence_attribution(inputs, write_text, tmp_path):
    attribution = write_text(
        "ev.tsv", "transcript_id\tevidence\nQ1.1\tScallop\nQX.1\tHelixer\n"
    )
    outdir = tmp_path / "out"
    summary = _run(
        [
            "--query",
            inputs["query"],
            "--reference",
            inputs["ref"],
            "--seqname-map",
            inputs["map"],
            "--region",
            "1:4000-6000",
            "--evidence-attribution",
            attribution,
            "--outdir",
            str(outdir),
        ]
    )
    assert (
        summary["total_reference_genes"] == 0 and summary["total_consensus_genes"] == 1
    )
    assert summary["subset_regions"] == ["1:4000-6000"]
    with open(outdir / "evidence_attribution_labeled.tsv") as handle:
        labels = {
            row["transcript_id"]: row["comparison_label"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert labels == {"Q1.1": "", "QX.1": ""}


def test_query_transcript_selection_canonical(write_gff3, inputs, tmp_path):
    records = gene_records("1", "Q1", "+", {"Q1.1": MULTI[:1], "Q1.2": MULTI})
    records[1] = (
        *records[1][:5],
        records[1][5] + ";tag=Ensembl_canonical",
    )  # Q1.1 mRNA
    query = write_gff3("canonical.gff3", records)
    common = [
        "--query",
        query,
        "--reference",
        inputs["ref"],
        "--seqname-map",
        inputs["map"],
        "--evaluation-mode",
        "protein_coding",
    ]
    every = _run(common + ["--outdir", str(tmp_path / "all")])
    canonical = _run(
        common
        + [
            "--query-transcript-selection",
            "canonical",
            "--outdir",
            str(tmp_path / "canonical"),
        ]
    )
    assert every["sensitivity"]["exact_match_count"] == 1
    assert canonical["sensitivity"]["exact_match_count"] == 0


@pytest.mark.skipif(
    shutil.which("annotation-qc") is None, reason="package not installed"
)
def test_installed_entry_point(inputs, tmp_path):
    result = subprocess.run(
        [
            "annotation-qc",
            "pairwise-compare",
            "--query",
            inputs["query"],
            "--reference",
            inputs["ref"],
            "--seqname-map",
            inputs["map"],
            "--outdir",
            str(tmp_path / "out"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "out" / "comparison_summary.json").is_file()

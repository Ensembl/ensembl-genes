# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Integration checks for the installed high-confidence workflow and framework layers."""
from __future__ import annotations

import argparse
import ast
import importlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.parsers.high_confidence import (
    parse_high_confidence_annotation,
    parse_tsv_frame,
    parse_legacy_diamond_alignments,
)
from ensembl.genes.annotation_qc.metrics.high_confidence import score_diamond_alignments
from ensembl.genes.annotation_qc.runners.high_confidence import add_arguments, run


def make_case(tmp_path, strand="+", spliced=False, canonical=True):
    """Create a small complete CDS, including a separately annotated stop codon."""
    genome = tmp_path / "genome.fa"
    annotation = tmp_path / "annotation.gtf"
    if spliced:
        motif = ("GT", "AG") if canonical else ("AA", "AA")
        sequence = "ATGA" + motif[0] + "C" * 46 + motif[1] + "AATAA"
        features = [
            ("exon", 1, 4, "."),
            ("exon", 55, 59, "."),
            ("CDS", 1, 4, "0"),
            ("CDS", 55, 56, "2"),
            ("stop_codon", 57, 59, "."),
        ]
    elif strand == "-":
        sequence = "TTATTTCAT"
        features = [("exon", 1, 9, "."), ("CDS", 4, 9, "0"), ("stop_codon", 1, 3, ".")]
    else:
        sequence = "ATGAAATAA"
        features = [("exon", 1, 9, "."), ("CDS", 1, 6, "0"), ("stop_codon", 7, 9, ".")]
    genome.write_text(">chr1\n" + sequence + "\n")
    attrs = 'gene_id "g1"; transcript_id "t1"; gene_biotype "protein_coding";'
    annotation.write_text(
        "".join(
            f"chr1\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"
            for kind, start, end, phase in features
        )
    )
    hits = tmp_path / "hits.tsv"
    hits.write_text("t1\tref\t0\t100\t100\t100\t100\n")
    return genome, annotation, hits


def run_case(tmp_path, **kwargs):
    """Run the orchestrator using its public CLI argument schema."""
    genome, annotation, hits = make_case(tmp_path, **kwargs)
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args(
        [
            "--annotation",
            str(annotation),
            "--genome",
            str(genome),
            "--diamond-results",
            str(hits),
            "--output-prefix",
            str(tmp_path / "result"),
            "--no-utr",
        ]
    )
    return args, run


@pytest.mark.parametrize("strand", ["+", "-"])
def test_sequences_and_annotation_coordinates(tmp_path, strand):
    """Reconstruct identical protein/CDS on both strands; keep 1-based exports."""
    args, execute = run_case(tmp_path, strand=strand)
    manifest = execute(args)
    evidence = parse_tsv_frame(tmp_path / "result.evidence.tsv")
    assert evidence.iloc[0]["transcript_cds_complete"] == "1"
    assert evidence.iloc[0]["confidence"] == "high"
    assert (tmp_path / "result.cds.fa").read_text(
        encoding="utf-8"
    ) == ">t1\nATGAAATAA\n"
    assert (tmp_path / "result.prots.fa").read_text(encoding="utf-8") == ">t1\nMK\n"
    assert manifest["confidence_policy_version"] == 2
    assert manifest["high_confidence_requires_intron_support"] is True
    gff = (tmp_path / "result.annotations.gff3").read_text(encoding="utf-8")
    assert "gene_quality=high" in gff
    assert "\tCDS\t1\t9\t" in gff


def test_split_codon_and_required_junction_support(tmp_path):
    """Internal phase bases complete a split codon rather than being trimmed."""
    args, execute = run_case(tmp_path, spliced=True)
    manifest = execute(args)
    frame = parse_tsv_frame(tmp_path / "result.evidence.tsv")
    assert frame.iloc[0]["transcript_intron_support"] == "."
    assert frame.iloc[0]["confidence"] == "medium"
    assert (tmp_path / "result.prots.fa").read_text(encoding="utf-8") == ">t1\nMK\n"
    junctions = tmp_path / "junctions.tsv"
    junctions.write_text("seqname\tstart\tend\tstrand\tread_count\nchr1\t5\t54\t+\t3\n")
    args.junctions = junctions
    args.min_junction_reads = 3
    args.output_prefix = tmp_path / "supported"
    manifest = execute(args)
    assert manifest["confidence_counts"] == {"high": 1, "medium": 0, "low": 0}
    intron = parse_tsv_frame(tmp_path / "supported.introns.tsv").iloc[0]
    assert (intron["start"], intron["end"], intron["length"]) == ("5", "54", "50")
    assert intron["splice_motif"] == "GT-AG"


def test_normalized_dataframe_boundary(tmp_path):
    """Expose common zero-based/half-open columns without returning legacy objects."""
    _, annotation, _ = make_case(tmp_path)
    frame = parse_high_confidence_annotation(annotation)
    assert isinstance(frame, pd.DataFrame)
    assert {
        "Chromosome",
        "Feature",
        "Start",
        "End",
        "Strand",
        "gene_id",
        "transcript_id",
    } <= set(frame.columns)
    exon = frame[frame["Feature"] == "exon"].iloc[0]
    assert (exon["Start"], exon["End"]) == (0, 9)


def test_metric_modules_have_no_file_io():
    """Reject accidental file/report/CLI operations in metric implementations."""
    forbidden = {
        "open",
        "read_text",
        "read_bytes",
        "write_text",
        "write_bytes",
        "read_csv",
        "to_csv",
        "parse_args",
    }
    for module_name in ["high_confidence", "high_confidence_summary"]:
        module = importlib.import_module(
            "ensembl.genes.annotation_qc.metrics." + module_name
        )
        tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                name = (
                    node.func.id
                    if isinstance(node.func, ast.Name)
                    else node.func.attr if isinstance(node.func, ast.Attribute) else ""
                )
                assert name not in forbidden, (module_name, name)


def test_diamond_parity_and_first_hit_ties(tmp_path):
    """Retain query-span target approximation, rounding, and first-hit tie choice."""
    path = tmp_path / "alignments.tsv"
    path.write_text(
        "t1\tfirst\t1\t110\t100\t90\t99\t110\t0\n"
        "t1\tsecond\t1\t110\t100\t90\t99\t110\t0\n"
    )
    results = score_diamond_alignments(
        parse_legacy_diamond_alignments(path), {"t1": 110}
    )
    assert results["t1"]["hit_id"] == "first"
    assert results["t1"]["target_coverage"] == 110
    assert results["t1"]["error"] == 10


def test_failed_coordinates_have_no_success_manifest(tmp_path):
    """Write diagnostics on coordinate mismatch without a misleading completion marker."""
    args, execute = run_case(tmp_path)
    args.gtf.write_text(args.gtf.read_text(encoding="utf-8").replace("\t9\t", "\t19\t"))
    with pytest.raises(ValueError, match="coordinate mismatch"):
        execute(args)
    assert (tmp_path / "result.input_validation.tsv").is_file()
    assert not (tmp_path / "result.manifest.json").exists()


def test_common_cli_and_summary_plots(tmp_path):
    """Exercise command registration, report tables and the optional rendering path."""
    genome, annotation, hits = make_case(tmp_path)
    common = [sys.executable, "-m", "ensembl.genes.annotation_qc.cli"]
    result = subprocess.run(
        common
        + [
            "high-confidence",
            "--annotation",
            str(annotation),
            "--genome",
            str(genome),
            "--diamond-results",
            str(hits),
            "--no-utr",
            "--output-prefix",
            str(tmp_path / "cli"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    env = dict(__import__("os").environ, MPLCONFIGDIR=str(tmp_path / "matplotlib"))
    plot_args = (
        ["--plots"] if importlib.util.find_spec("matplotlib") is not None else []
    )
    result = subprocess.run(
        common
        + [
            "high-confidence-summary",
            "--evidence",
            str(tmp_path / "cli.evidence.tsv"),
            "--label",
            "Method",
            "--species",
            "test_species",
            "--output-dir",
            str(tmp_path / "summary"),
        ]
        + plot_args,
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "summary/flag_distribution.tsv").is_file()
    if plot_args:
        assert (tmp_path / "summary/figures/evidence_flags_heatmap.png").is_file()
        assert (tmp_path / "summary/figures/confidence_sankey.svg").is_file()
    metadata = json.loads(
        (tmp_path / "summary/provenance.json").read_text(encoding="utf-8")
    )
    assert metadata["run_manifest"]["confidence_policy_version"] == 2

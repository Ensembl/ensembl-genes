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
"""Final EAM summary counts genes rather than repeated transcript rows."""

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.metrics.gene_confidence_summary import (
    summarize_gene_confidence,
)
from ensembl.genes.annotation_qc.reports.gene_confidence_summary import (
    write_gene_confidence_summary,
)
from ensembl.genes.annotation_qc.runners.gene_confidence_summary import run


def frame(rows):
    """Construct a minimal scored-evidence table."""
    return pd.DataFrame(rows, columns=["gene_id", "confidence"], dtype=object)


def test_counts_unique_genes_and_exports_exact_schema(tmp_path):
    """Duplicate isoforms do not inflate gene totals or percentages."""
    evidence = frame(
        [
            ("g1", "high"),
            ("g1", "high"),
            ("g2", "medium"),
            ("g3", "low"),
            ("g4", "high"),
        ]
    )
    before = evidence.copy(deep=True)
    metrics = summarize_gene_confidence(evidence)
    assert list(metrics.items()) == [
        ("eam.total_genes", 4),
        ("eam.total_high_confidence_genes", 2),
        ("eam.total_high_confidence_genes_pct", 50.0),
        ("eam.total_medium_confidence_genes", 1),
        ("eam.total_medium_confidence_genes_pct", 25.0),
        ("eam.total_low_confidence_genes", 1),
        ("eam.total_low_confidence_genes_pct", 25.0),
    ]
    pd.testing.assert_frame_equal(evidence, before)
    output = tmp_path / "reports/final.csv"
    write_gene_confidence_summary(metrics, output)
    assert output.read_text() == (
        "metric,value\n"
        "eam.total_genes,4\n"
        "eam.total_high_confidence_genes,2\n"
        "eam.total_high_confidence_genes_pct,50.00\n"
        "eam.total_medium_confidence_genes,1\n"
        "eam.total_medium_confidence_genes_pct,25.00\n"
        "eam.total_low_confidence_genes,1\n"
        "eam.total_low_confidence_genes_pct,25.00\n"
    )


def test_empty_and_missing_category():
    """Empty evidence is explicit and absent confidence classes have zero counts."""
    empty = summarize_gene_confidence(frame([]))
    assert len(empty) == 7
    assert all(value == 0 for value in empty.values())
    result = summarize_gene_confidence(frame([("g1", "high")]))
    assert result["eam.total_high_confidence_genes_pct"] == 100
    assert result["eam.total_medium_confidence_genes"] == 0
    assert result["eam.total_low_confidence_genes_pct"] == 0


def test_rounding_is_based_on_unique_gene_denominator():
    """Each third rounds to 33.33 percent."""
    result = summarize_gene_confidence(
        frame([("g1", "high"), ("g2", "medium"), ("g3", "low")])
    )
    assert [v for key, v in result.items() if key.endswith("_pct")] == [33.33] * 3


@pytest.mark.parametrize(
    "rows",
    [
        [("g1", "high"), ("g1", "low")],
        [(None, "high")],
        [("", "high")],
        [(" ", "high")],
        [(".", "high")],
        [("g1", "unknown")],
        [("g1", None)],
    ],
)
def test_invalid_assignments_rejected(rows):
    """Invalid inputs never silently change totals or select an arbitrary class."""
    with pytest.raises(ValueError):
        summarize_gene_confidence(frame(rows))


def test_missing_columns_rejected():
    """Require the two fields needed to compute gene totals."""
    with pytest.raises(ValueError, match="Missing required"):
        summarize_gene_confidence(pd.DataFrame({"gene_id": ["g1"]}))


def test_input_cannot_be_overwritten(tmp_path):
    """Protect evidence when an output filename is accidentally reused."""
    source = tmp_path / "input.tsv"
    source.write_text("gene_id\tconfidence\ng1\thigh\n")
    before = source.read_bytes()
    with pytest.raises(ValueError, match="overwrite"):
        run(argparse.Namespace(evidence=source, output=source))
    assert source.read_bytes() == before


def test_shared_cli_generates_csv(tmp_path):
    """Exercise registration and the complete parse/calculate/write command."""
    source = tmp_path / "input.evidence.tsv"
    source.write_text("gene_id\tconfidence\ng1\thigh\ng2\tlow\n")
    output = tmp_path / "summary.csv"
    package_root = str(Path(__file__).resolve().parents[4] / "src/python")
    env = dict(
        os.environ,
        PYTHONPATH=package_root + os.pathsep + os.environ.get("PYTHONPATH", ""),
    )
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "ensembl.genes.annotation_qc.cli",
            "high-confidence-gene-summary",
            "--evidence",
            str(source),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert result.returncode == 0, result.stderr
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 7
    assert rows[0] == {"metric": "eam.total_genes", "value": "2"}
    assert rows[2]["value"] == "50.00"

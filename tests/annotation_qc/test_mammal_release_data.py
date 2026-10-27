"""
Integration checks on the real mammalian dashboard data (opt-in).

Skipped unless ANNOTATION_QC_MAMMAL_CONFIGS names a folder of species configurations
(e.g. <data>/qc_configs/mammal_dashboard). Needs the built datasets
(qc_workspace/mammal_<species>/dashboard/dashboard.sqlite); the comparator output
bundles are used when present (metric-by-metric check against metrics.tsv), otherwise
that check is skipped for the bundle.

    ANNOTATION_QC_MAMMAL_CONFIGS=/path/to/qc_configs/mammal_dashboard \
        python -m pytest tests/annotation_qc/test_mammal_release_data.py
"""

import csv
import json
import math
import os
from pathlib import Path

import pytest

from ensembl.genes.annotation_qc.benchmark.config import load_experiment
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
from ensembl.genes.annotation_qc.dashboard import queries

ROOT = os.environ.get("ANNOTATION_QC_MAMMAL_CONFIGS")
CONFIGS = sorted(Path(ROOT).glob("*.json")) if ROOT else []
pytestmark = pytest.mark.skipif(
    not CONFIGS, reason="set ANNOTATION_QC_MAMMAL_CONFIGS to the species configuration folder"
)
FOCUS_TOOL = "Internal predictions"


@pytest.fixture(params=CONFIGS, ids=[p.stem for p in CONFIGS])
def species(request):
    experiment = load_experiment(request.param)
    db = Workspace(experiment).dashboard_db()
    if not db.is_file():
        pytest.fail(f"{request.param.name}: dataset not built: {db}")
    con = queries.connect(db)
    yield experiment, con
    con.close()


def test_one_reference_and_the_configured_runs(species):
    experiment, con = species
    meta = queries.meta(con)
    assert len(experiment.references) == 1 and len(meta["references"]) == 1
    (ref,) = experiment.references.values()
    assert meta["references"][0]["id"] == ref.id
    assert {r["run_id"]: r["tool"] for r in meta["runs"]} == {
        r.run_id: r.tool for r in experiment.runs.values()
    }
    assert {r["run_id"]: r["display_role"] for r in meta["runs"]} == {
        r.run_id: r.display_role for r in experiment.runs.values()
    }


def test_every_result_is_complete_and_compatible(species):
    experiment, con = species
    (ref,) = experiment.references.values()
    results = queries.overview(con, ref.id)["results"]
    assert len(results) == len(experiment.runs) * len(experiment.policies)
    assert {r["status"] for r in results} == {"complete"}, [
        (r["run_id"], r["policy_id"], r["status"], r["reason"]) for r in results if r["status"] != "complete"
    ]
    assert len({r["definition_version"] for r in results}) == 1
    assert len({(r.get("identity") or {}).get("build_id") for r in results}) == 1


def test_internal_predictions_are_the_stop_corrected_comparison(species):
    experiment, con = species
    for run in experiment.runs.values():
        if run.tool == FOCUS_TOOL:
            assert run.comparator_options.get("add_stop_codon"), run.run_id
    (ref,) = experiment.references.values()
    for r in queries.overview(con, ref.id)["results"]:
        if experiment.runs[r["run_id"]].tool == FOCUS_TOOL:
            assert (r.get("stop_codon_harmonisation") or {}), r["run_id"]


def test_served_metrics_equal_metrics_tsv(species):
    experiment, con = species
    (ref,) = experiment.references.values()
    ov = queries.overview(con, ref.id)
    served = {(m["run_id"], m["policy_id"], m["scope_id"], m["metric_id"]): m for m in ov["metrics"]}
    compared = 0
    for r in ov["results"]:
        metrics = Path(r["output_dir"]) / "metrics.tsv"
        if not metrics.is_file():
            continue
        with open(metrics) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                m = served[(r["run_id"], r["policy_id"], row["scope_id"], row["metric_id"])]
                assert m["status"] == row["status"]
                if row["status"] == "ok":
                    assert math.isclose(float(row["value"]), float(m["value"]), rel_tol=1e-12)
                    assert str(m["numerator"]) == row["numerator"]
                compared += 1
    if not compared:
        pytest.skip("comparator output bundles not present (runtime-only data package)")


def test_scopes_offered_are_scopes_reported(species):
    experiment, con = species
    (ref,) = experiment.references.values()
    meta = queries.meta(con)
    offered = {s["id"] for s in meta["references"][0]["scopes"] if s["available"]}
    reported = {m["scope_id"] for m in queries.overview(con, ref.id)["metrics"]}
    assert offered == reported


def test_species_differ_in_their_predictors():
    tools = {p.stem: {r.tool for r in load_experiment(p).runs.values()} for p in CONFIGS}
    if len(tools) > 1:
        assert len({frozenset(t) for t in tools.values()}) > 1, json.dumps({k: sorted(v) for k, v in tools.items()})

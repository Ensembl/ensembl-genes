"""
Release smoke test on the SYNTHETIC handover example
(src/python/ensembl/genes/annotation_qc/examples/handover_mini; no private data).

The complete handover workflow:

    reference + prediction annotations
      -> `annotation-qc pairwise-compare` (the CLI module, as installed)
      -> saved output bundles (contract 1.0)
      -> species configuration (paths relative to the configuration file)
      -> `benchmark import` + `benchmark build` -> dashboard dataset
      -> dashboard loader / HTTP API

and checks that served values equal metrics.tsv, that the project folder can be moved,
that display roles do not change any value and that missing results are reported.
"""

import csv
import json
import math
import shutil
import subprocess
import sys
import threading
import urllib.request
from pathlib import Path

import pytest

import ensembl.genes.annotation_qc as aqc
from ensembl.genes.annotation_qc.benchmark.cli import checkpoint_dataset
from ensembl.genes.annotation_qc.benchmark.config import load_experiment
from ensembl.genes.annotation_qc.benchmark.importer import import_benchmark
from ensembl.genes.annotation_qc.benchmark.store import build_dashboard
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
from ensembl.genes.annotation_qc.dashboard import queries

QUIET = lambda *a, **k: None  # noqa: E731
CLI = [sys.executable, "-m", "ensembl.genes.annotation_qc.cli"]
POLICIES = {"A": ("all", "all"), "B": ("longest_cds", "longest_cds")}
# (run id, variant, query, extra options): as in examples/handover_mini/run_example.sh
RUNS = [
    ("mini_internal", "stop_harmonised", "focus.gff3", ["--add-stop-codon", "--stop-codon-genetic-code", "1"]),
    ("mini_other", "original", "other.gff3", []),
]


def _example() -> Path:
    """The example folder: in the source tree next to the package (tests run from a checkout)."""
    candidates = [
        Path(aqc.__file__).parent / "examples" / "handover_mini",
        Path(__file__).resolve().parents[2] / "src/python/ensembl/genes/annotation_qc/examples/handover_mini",
    ]
    for path in candidates:
        if (path / "configs" / "mini_species.json").is_file():
            return path
    pytest.skip("handover example not found (tests need a source checkout)")


def _compare(root: Path, run: str, variant: str, policy: str, query: str, extra: list) -> Path:
    """One comparison through the CLI entry module, exactly as a user runs it."""
    ref_sel, query_sel = POLICIES[policy]
    out = Path("results") / run / variant / policy
    argv = CLI + [
        "pairwise-compare",
        "--reference", "inputs/reference.gff3", "--query", f"inputs/{query}",
        "--genome", "inputs/genome.fa", "--genome-mismatch", "error",
        "--evaluation-mode", "protein_coding", "--reference-transcript-biotypes", "protein_coding",
        "--reference-transcript-selection", ref_sel, "--query-transcript-selection", query_sel,
        "--analysis-scope", "chromosomes=1", "--export-evaluated-models", "ids",
        *extra, "--outdir", str(out),
    ]  # fmt: skip
    done = subprocess.run(argv, cwd=root, capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr[-2000:]
    return root / out


def _variant_config(root: Path, name: str, roles=True, results=None) -> Path:
    """A copy of the example configuration (same folder, so relative paths stay valid)."""
    config = json.loads((root / "configs" / "mini_species.json").read_text())
    config["experiment"]["id"] = name
    config["output_dir"] = f"../workspace_{name}"
    if not roles:
        for run in config["runs"]:
            run.pop("display_role", None)
    if results is not None:
        config["import"]["results_dirs"] = results
    path = root / "configs" / f"{name}.json"
    path.write_text(json.dumps(config, indent=1))
    return path


def _index(config: Path):
    experiment = load_experiment(config)
    ws = Workspace(experiment)
    import_benchmark(experiment, ws, None, log=QUIET)
    return build_dashboard(experiment, ws, log=QUIET)


def _served(db: Path) -> dict:
    con = queries.connect(db)
    try:
        ov = queries.overview(con, "mini_ref")
        meta = queries.meta(con)
    finally:
        con.close()
    rows = {(m["run_id"], m["policy_id"], m["scope_id"], m["metric_id"]): m for m in ov["metrics"]}
    return {"rows": rows, "results": ov["results"], "meta": meta}


def _number(text):
    return None if text in ("", "NA", None) else float(text)


@pytest.fixture(scope="module")
def project(tmp_path_factory):
    root = tmp_path_factory.mktemp("handover") / "project"
    shutil.copytree(_example(), root, ignore=shutil.ignore_patterns("results", "workspace*"))
    for run, variant, query, extra in RUNS:
        for policy in POLICIES:
            _compare(root, run, variant, policy, query, extra)
    config = root / "configs" / "mini_species.json"
    return {"root": root, "config": config, "db": _index(config)}


def test_comparator_writes_a_contract_bundle(project):
    out = project["root"] / "results" / "mini_internal" / "stop_harmonised" / "B"
    manifest = json.loads((out / "comparison_manifest.json").read_text())
    assert manifest["contract_version"] == "1.0" and manifest["definition_version"]
    assert manifest["options"]["add_stop_codon"] is True
    for name in ("metrics.tsv", "comparison_summary.json", "comparison_details.tsv",
                 "cds_exact_one_to_one_pairs.tsv", "evaluated_transcripts.tsv.gz",
                 "splice_sites_by_position.tsv", "cds_introns_by_length.tsv",
                 "stop_codon_harmonisation.tsv"):  # fmt: skip
        assert (out / name).is_file(), name
        assert name in manifest["outputs"], name
    with open(out / "metrics.tsv") as handle:
        rows = {(r["scope_id"], r["metric_id"]): r for r in csv.DictReader(handle, delimiter="\t")}
    # R1 (stop restored), R2 and R4 exact of R = 4; R3 missed; scope "chromosomes" = chr 1
    assert rows[("all", "reference_genes")]["numerator"] == "4"
    assert rows[("all", "exact_cds_one_to_one_tp")]["numerator"] == "3"
    assert rows[("all", "stop_codon_extended")]["numerator"] == "1"
    assert rows[("chromosomes", "reference_genes")]["numerator"] == "3"


def test_dashboard_serves_metrics_tsv_values(project):
    served = _served(project["db"])
    status = {(r["run_id"], r["policy_id"]): r["status"] for r in served["results"]}
    assert {s for (_, p), s in status.items() if p in POLICIES} == {"complete"}
    assert status[("mini_internal", "C")] == "unavailable"  # no Ensembl_canonical tags
    checked = 0
    for run, variant, _, _ in RUNS:
        for policy in POLICIES:
            out = project["root"] / "results" / run / variant / policy
            with open(out / "metrics.tsv") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    m = served["rows"].pop((run, policy, row["scope_id"], row["metric_id"]))
                    assert m["status"] == row["status"]
                    for field in ("numerator", "denominator"):
                        assert _number(row[field]) == (None if m[field] is None else float(m[field]))
                    if row["status"] == "ok":
                        assert math.isclose(float(row["value"]), float(m["value"]), rel_tol=1e-12)
                    checked += 1
    assert checked > 100
    assert not served["rows"], "dashboard rows without a metrics.tsv row"
    recall = _served(project["db"])["rows"][("mini_internal", "B", "all", "exact_cds_one_to_one_recall")]
    assert (recall["numerator"], recall["denominator"]) == (3, 4)
    scopes = {s["id"]: s["available"] for s in served["meta"]["references"][0]["scopes"]}
    assert scopes == {"all": True, "chromosomes": True}
    roles = {r["run_id"]: r["display_role"] for r in served["meta"]["runs"]}
    assert roles == {"mini_internal": "focus", "mini_other": "primary_comparator"}


def test_build_leaves_one_self_contained_file(project):
    db = project["db"]
    assert checkpoint_dataset(db)
    wal = db.with_name(db.name + "-wal")
    assert not wal.exists() or wal.stat().st_size == 0


def test_project_folder_can_be_moved(project, tmp_path):
    moved = tmp_path / "elsewhere" / "copy"
    shutil.copytree(project["root"], moved)
    experiment = load_experiment(moved / "configs" / "mini_species.json")
    db = Workspace(experiment).dashboard_db()
    assert db.is_file() and moved in db.parents
    original = _served(project["db"])["rows"]
    assert _served(db)["rows"] == original

    from http.server import ThreadingHTTPServer

    from ensembl.genes.annotation_qc.dashboard.server import Datasets, make_handler

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(Datasets({experiment.id: db})))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        base = f"http://127.0.0.1:{httpd.server_address[1]}"
        (entry,) = json.loads(urllib.request.urlopen(base + "/api/experiments").read())
        assert entry["species"] == "mini species (synthetic)"
        ov = json.loads(urllib.request.urlopen(base + "/api/overview?reference=mini_ref").read())
        assert len(ov["metrics"]) == len(original)
    finally:
        httpd.shutdown()


def test_display_roles_change_no_value(project):
    plain = _served(_index(_variant_config(project["root"], "no_roles", roles=False)))
    roled = _served(project["db"])
    assert plain["rows"] == roled["rows"]
    assert {r["run_id"]: r["display_role"] for r in plain["meta"]["runs"]} == {
        "mini_internal": None,
        "mini_other": None,
    }


def test_missing_results_are_reported(project):
    root = project["root"]
    # a run whose results were not supplied: "not run", with no metric rows
    config = _variant_config(root, "partial", results=["../results/mini_internal/stop_harmonised"])
    served = _served(_index(config))
    status = {(r["run_id"], r["policy_id"]): r["status"] for r in served["results"]}
    assert status[("mini_other", "A")] == "not_run" and status[("mini_internal", "A")] == "complete"
    assert not any(k[0] == "mini_other" for k in served["rows"])
    # a results directory that does not exist: a clear error naming it, no traceback
    config = _variant_config(root, "bad_dir", results=["../results/mini_internal/stop_harmonised", "../results/nonexistent"])
    done = subprocess.run(CLI + ["benchmark", "import", "--config", str(config)], capture_output=True, text=True, check=False)
    assert done.returncode != 0
    assert "result directory not found" in done.stderr and "nonexistent" in done.stderr
    assert "Traceback" not in done.stderr


def test_dashboard_refuses_a_missing_dataset(project):
    config = _variant_config(project["root"], "never_built")
    done = subprocess.run(
        CLI + ["dashboard", "--config", str(config), "--port", "0"],
        capture_output=True, text=True, check=False, timeout=60,
    )  # fmt: skip
    assert done.returncode != 0
    assert "Dashboard dataset not built" in done.stderr + done.stdout
    assert "Traceback" not in done.stderr


def test_port_in_use_is_reported(project):
    import socket

    config = project["root"] / "configs" / "mini_species.json"
    with socket.socket() as busy:
        busy.bind(("127.0.0.1", 0))
        busy.listen()
        port = busy.getsockname()[1]
        done = subprocess.run(
            CLI + ["dashboard", "--config", str(config), "--port", str(port)],
            capture_output=True, text=True, check=False, timeout=60,
        )  # fmt: skip
    assert done.returncode != 0
    assert f"Cannot listen on 127.0.0.1:{port}" in done.stderr and "--port" in done.stderr
    assert "Traceback" not in done.stderr

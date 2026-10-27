"""
Release smoke test (SYNTHETIC, small fixtures; no private data, no absolute paths).

The complete handover workflow on a tiny genome:

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
import sqlite3
import subprocess
import sys
import threading
import urllib.request
from pathlib import Path

import pytest

from ensembl.genes.annotation_qc.benchmark.cli import checkpoint_dataset
from ensembl.genes.annotation_qc.benchmark.config import load_experiment
from ensembl.genes.annotation_qc.benchmark.importer import import_benchmark
from ensembl.genes.annotation_qc.benchmark.store import build_dashboard
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
from ensembl.genes.annotation_qc.dashboard import queries

from .conftest import gene_records, gff3_text

QUIET = lambda *a, **k: None  # noqa: E731
POLICIES = {"A": ("all", "all"), "B": ("longest_cds", "longest_cds")}
MULTI = [(101, 200), (301, 400), (501, 600)]
MINUS = [(1101, 1200), (1301, 1400)]


def _ensembl_like(records):
    """Ensembl GFF3 carries the biotype on transcripts too (used by the biotype filter)."""
    return [(*r[:5], r[5] + ";biotype=protein_coding") if r[1] == "mRNA" else r for r in records]


def _write_inputs(root: Path) -> None:
    inputs = root / "inputs"
    inputs.mkdir(parents=True)
    (inputs / "ref.gff3").write_text(
        gff3_text(
            _ensembl_like(
                gene_records("1", "R1", "+", {"R1.1": MULTI, "R1.2": MULTI[:2]})
                + gene_records("1", "R2", "-", {"R2.1": MINUS})
                + gene_records("1", "R3", "+", {"R3.1": [(2001, 2300)]})
                + gene_records("2", "R4", "+", {"R4.1": [(101, 400)]})
            )
        )
    )
    # "focus" predictor: R1 and R2 exact, R3 missed, R4 exact
    (inputs / "focus.gff3").write_text(
        gff3_text(
            gene_records("1", "F1", "+", {"F1.1": MULTI})
            + gene_records("1", "F2", "-", {"F2.1": MINUS})
            + gene_records("2", "F4", "+", {"F4.1": [(101, 400)]})
        )
    )
    # comparison predictor: R1 exact, R2 different splice site, R3 exact
    (inputs / "other.gff3").write_text(
        gff3_text(
            gene_records("1", "O1", "+", {"O1.1": MULTI})
            + gene_records("1", "O2", "-", {"O2.1": [(1101, 1180), (1301, 1400)]})
            + gene_records("1", "O3", "+", {"O3.1": [(2001, 2300)]})
        )
    )
    (inputs / "genome.fa").write_text(">1\n" + "ACGT" * 1000 + "\n>2\n" + "ACGT" * 250 + "\n")


def _compare(root: Path, prediction: str, policy: str) -> Path:
    """One comparison through the CLI entry module, exactly as a user runs it."""
    ref_sel, query_sel = POLICIES[policy]
    out = root / "results" / prediction / "original" / policy
    argv = [
        sys.executable, "-m", "ensembl.genes.annotation_qc.cli", "pairwise-compare",
        "--query", f"inputs/{prediction}.gff3",
        "--reference", "inputs/ref.gff3",
        "--genome", "inputs/genome.fa", "--genome-mismatch", "error",
        "--evaluation-mode", "protein_coding",
        "--reference-transcript-biotypes", "protein_coding",
        "--reference-transcript-selection", ref_sel,
        "--query-transcript-selection", query_sel,
        "--analysis-scope", "chromosomes=1",
        "--export-evaluated-models", "ids",
        "--outdir", str(out.relative_to(root)),
    ]  # fmt: skip
    done = subprocess.run(argv, cwd=root, capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr[-2000:]
    return out


def _config(root: Path, name: str, runs: list[dict], workspace: str, results: list[str]):
    config = {
        "config_version": 2,
        "experiment": {"id": name, "title": f"Synthetic species ({name})", "synthetic": True},
        "output_dir": f"../{workspace}",
        "evaluation": {
            "evaluation_mode": "protein_coding",
            "reference_gene_biotypes": [],
            "reference_transcript_biotypes": ["protein_coding"],
            "genome_mismatch": "error",
        },
        "policies": {
            pid: {"label": pid, "reference_transcript_selection": r, "query_transcript_selection": q}
            for pid, (r, q) in POLICIES.items()
        },
        "references": [
            {
                "id": "synth_ref",
                "species": "synthetic species",
                "assembly": "SYN1",
                "annotation": "../inputs/ref.gff3",
                "genome": "../inputs/genome.fa",
                "scope": {"type": "all", "label": "Whole assembly"},
                "analysis_scopes": [{"id": "chromosomes", "label": "Chromosome 1", "sequences": ["1"]}],
            }
        ],
        "runs": runs,
        "import": {"results_dirs": [f"../results/{r}/original" for r in results]},
    }
    path = root / "configs" / f"{name}.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(config, indent=1))
    return path


def _runs(roles: bool) -> list[dict]:
    return [
        {"run_id": "other", "reference": "synth_ref", "tool": "Other tool", "annotation": "../inputs/other.gff3",
         **({"display_role": "primary_comparator"} if roles else {})},
        {"run_id": "focus", "reference": "synth_ref", "tool": "Internal predictions", "annotation": "../inputs/focus.gff3",
         **({"display_role": "focus"} if roles else {})},
    ]  # fmt: skip


def _index(config: Path):
    experiment = load_experiment(config)
    ws = Workspace(experiment)
    import_benchmark(experiment, ws, None, log=QUIET)
    db = build_dashboard(experiment, ws, log=QUIET)
    return experiment, ws, db


def _served(db: Path) -> dict:
    con = queries.connect(db)
    try:
        ov = queries.overview(con, "synth_ref")
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
    _write_inputs(root)
    for prediction in ("focus", "other"):
        for policy in POLICIES:
            _compare(root, prediction, policy)
    config = _config(root, "synthetic_species", _runs(roles=True), "workspace", ["focus", "other"])
    experiment, ws, db = _index(config)
    return {"root": root, "config": config, "db": db}


def test_comparator_writes_a_contract_bundle(project):
    out = project["root"] / "results" / "focus" / "original" / "B"
    manifest = json.loads((out / "comparison_manifest.json").read_text())
    assert manifest["contract_version"] == "1.0" and manifest["definition_version"]
    for name in ("metrics.tsv", "comparison_summary.json", "comparison_details.tsv",
                 "cds_exact_one_to_one_pairs.tsv", "evaluated_transcripts.tsv.gz",
                 "splice_sites_by_position.tsv", "cds_introns_by_length.tsv"):  # fmt: skip
        assert (out / name).is_file(), name
        assert name in manifest["outputs"], name
    with open(out / "metrics.tsv") as handle:
        rows = {(r["scope_id"], r["metric_id"]): r for r in csv.DictReader(handle, delimiter="\t")}
    # focus: R1, R2, R4 exact of R = 4 evaluated reference genes; scope "chromosomes" = chr 1
    assert rows[("all", "reference_genes")]["numerator"] == "4"
    assert rows[("all", "exact_cds_one_to_one_tp")]["numerator"] == "3"
    assert rows[("chromosomes", "reference_genes")]["numerator"] == "3"


def test_dashboard_serves_metrics_tsv_values(project):
    served = _served(project["db"])
    assert {r["status"] for r in served["results"]} == {"complete"}
    checked = 0
    for out in sorted((project["root"] / "results").glob("*/original/*")):
        run, policy = out.parts[-3], out.parts[-1]
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
    focus_b = _served(project["db"])["rows"][("focus", "B", "all", "exact_cds_one_to_one_recall")]
    assert (focus_b["numerator"], focus_b["denominator"]) == (3, 4)
    assert not served["rows"], "dashboard rows without a metrics.tsv row"
    scopes = {s["id"]: s["available"] for s in served["meta"]["references"][0]["scopes"]}
    assert scopes == {"all": True, "chromosomes": True}
    runs = {r["run_id"]: r for r in served["meta"]["runs"]}
    assert runs["focus"]["display_role"] == "focus" and runs["other"]["display_role"] == "primary_comparator"


def test_build_leaves_one_self_contained_file(project):
    db = project["db"]
    assert checkpoint_dataset(db)
    wal = db.with_name(db.name + "-wal")
    assert not wal.exists() or wal.stat().st_size == 0


def test_project_folder_can_be_moved(project, tmp_path):
    moved = tmp_path / "elsewhere" / "copy"
    shutil.copytree(project["root"], moved)
    config = moved / "configs" / project["config"].name
    experiment = load_experiment(config)
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
        assert entry["species"] == "synthetic species"
        ov = json.loads(urllib.request.urlopen(base + "/api/overview?reference=synth_ref").read())
        assert len(ov["metrics"]) == len(original)
    finally:
        httpd.shutdown()


def test_display_roles_change_no_value(project):
    config = _config(project["root"], "no_roles", _runs(roles=False), "workspace_no_roles", ["focus", "other"])
    _, _, db = _index(config)
    plain, roled = _served(db), _served(project["db"])
    strip = lambda rows: {k: {f: v for f, v in m.items()} for k, m in rows.items()}  # noqa: E731
    assert strip(plain["rows"]) == strip(roled["rows"])
    assert {r["run_id"]: r["display_role"] for r in plain["meta"]["runs"]} == {"focus": None, "other": None}


def test_missing_results_are_reported(project):
    root = project["root"]
    # a run whose results were not supplied: shown as "not run", with no metric rows
    config = _config(root, "partial", _runs(roles=True), "workspace_partial", ["focus"])
    _, _, db = _index(config)
    served = _served(db)
    status = {(r["run_id"], r["policy_id"]): r["status"] for r in served["results"]}
    assert status[("other", "A")] == "not_run" and status[("focus", "A")] == "complete"
    assert not any(k[0] == "other" for k in served["rows"])
    # a results directory that does not exist: a clear error naming it
    config = _config(root, "bad_dir", _runs(roles=True), "workspace_bad", ["focus", "nonexistent"])
    done = subprocess.run(
        [sys.executable, "-m", "ensembl.genes.annotation_qc.cli", "benchmark", "import", "--config", str(config)],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert done.returncode != 0
    assert "result directory not found" in done.stderr and "nonexistent" in done.stderr
    assert "Traceback" not in done.stderr


def test_dashboard_refuses_a_missing_dataset(project, tmp_path):
    config = _config(project["root"], "never_built", _runs(roles=True), "workspace_never_built", ["focus"])
    done = subprocess.run(
        [sys.executable, "-m", "ensembl.genes.annotation_qc.cli", "dashboard", "--config", str(config), "--port", "0"],
        capture_output=True, text=True, check=False, timeout=60,
    )  # fmt: skip
    assert done.returncode != 0
    assert "Dashboard dataset not built" in done.stderr + done.stdout
    assert "Traceback" not in done.stderr

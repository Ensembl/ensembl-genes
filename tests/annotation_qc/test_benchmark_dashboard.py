"""Experiment configuration, caching, import, dataset build and dashboard API (SYNTHETIC fixtures)."""

import hashlib
import json
import shutil
import sqlite3
import threading
import urllib.request
from pathlib import Path

import pytest

from ensembl.genes.annotation_qc.benchmark import runs as runs_mod
from ensembl.genes.annotation_qc.benchmark.config import ConfigError, load_experiment
from ensembl.genes.annotation_qc.benchmark.derive import derive_metrics
from ensembl.genes.annotation_qc.benchmark.examples.synthetic import (
    write_synthetic_demo,
)
from ensembl.genes.annotation_qc.benchmark.importer import import_benchmark
from ensembl.genes.annotation_qc.benchmark.independent.cds_signatures import (
    run_independent,
)
from ensembl.genes.annotation_qc.benchmark.preparation import (
    PreparationError,
    prepare_reference,
    prepare_run,
)
from ensembl.genes.annotation_qc.benchmark.provenance import (
    FINGERPRINT_FILES,
    PACKAGE_DIR,
)
from ensembl.genes.annotation_qc.benchmark.store import build_dashboard
from ensembl.genes.annotation_qc.benchmark.validation import validate_experiment
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
from ensembl.genes.annotation_qc.dashboard import queries

QUIET = lambda *a, **k: None  # noqa: E731


def _statuses(experiment, ws):
    return {
        (t.run.run_id, t.policy.id): t.status for t in runs_mod.plan(experiment, ws)
    }


def _independent(experiment, ws):
    for ref in experiment.references.values():
        prep = prepare_reference(ws, ref)
        canonical = runs_mod.canonical_available(ws, prep)
        policies = {
            p.id: (p.reference_transcript_selection, p.query_transcript_selection)
            for p in experiment.policies.values()
            if runs_mod.policy_available(p, canonical)[0]
        }
        run_independent(
            prep["annotation"]["output"]["path"],
            set(prep["scope"]["sequences"] or []) or None,
            experiment.evaluation,
            {
                r.run_id: prepare_run(ws, r)["output"]["path"]
                for r in experiment.runs_for(ref.id)
            },
            policies,
            ws.root / "independent",
            log=QUIET,
        )


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    root = tmp_path_factory.mktemp("synthetic")
    config = write_synthetic_demo(root)
    experiment = load_experiment(config)
    ws = Workspace(experiment)
    prepared = {r.id: prepare_reference(ws, r) for r in experiment.references.values()}
    report = validate_experiment(experiment, ws, prepared)
    runs_mod.run_missing(experiment, ws, log=QUIET)
    _independent(experiment, ws)
    build_dashboard(experiment, ws, log=QUIET)
    return {
        "root": root,
        "config": config,
        "experiment": experiment,
        "ws": ws,
        "validation": report,
    }


@pytest.fixture()
def con(demo):
    connection = queries.connect(demo["ws"].dashboard_db())
    yield connection
    connection.close()


# --------------------------------------------------------------- configuration
def test_config_identity_rules(tmp_path):
    base = {
        "config_version": 1,
        "experiment": {"id": "x"},
        "output_dir": "out",
        "references": [
            {"id": "r", "species": "s", "assembly": "a", "annotation": "ref.gff3"}
        ],
        "runs": [
            {"run_id": "one", "reference": "r", "tool": "T", "annotation": "a.gff3"},
            {"run_id": "two", "reference": "r", "tool": "T", "annotation": "b.gff3"},
        ],
    }
    path = tmp_path / "c.json"
    path.write_text(json.dumps(base))
    exp = load_experiment(path)
    assert set(exp.runs) == {"one", "two"} and {r.tool for r in exp.runs.values()} == {
        "T"
    }
    assert (
        exp.runs["one"].annotation == (tmp_path / "a.gff3").resolve()
    )  # relative to the config
    bad = json.loads(json.dumps(base))
    bad["runs"][1]["run_id"] = "one"
    bad["policies"] = {"A": {"reference_transcript_selection": "nonsense"}}
    path.write_text(json.dumps(bad))
    with pytest.raises(ConfigError) as error:
        load_experiment(path)
    assert any("duplicate run_id" in p for p in error.value.problems)
    assert any(
        "policies.A.reference_transcript_selection" in p for p in error.value.problems
    )


def test_template_loads():
    from ensembl.genes.annotation_qc.benchmark import config as config_mod

    template = (
        Path(config_mod.__file__).with_name("templates") / "experiment_template.json"
    )
    exp = load_experiment(template)
    assert set(exp.policies) == {"A", "B", "C"} and len(exp.runs) == 2
    assert {r.tool for r in exp.runs.values()} == {"ToolA"}


def test_validation_reports_and_probes(demo):
    report = demo["validation"]
    messages = " ".join(i["message"] for i in report["issues"])
    assert "Parent is not a record" in messages  # broken run flagged before running
    assert "reused on different sequences" in messages  # collision run noted
    assert report["runs"]["SYNTH_toolX_m1"]["outside_scope"] == {
        "sequences": 1,
        "genes": 1,
    }
    statuses = {lim["id"]: lim["status"] for lim in report["parser_limitations"]}
    assert set(statuses) == {
        "bgz_suffix",
        "unlisted_transcript_types",
        "duplicate_child_rows",
        "gtf_inferred_transcript_fusion",
    }


# --------------------------------------------------------------- runs, failures, unavailable policies
def test_results_failures_and_unavailable_policy(demo):
    st = _statuses(demo["experiment"], demo["ws"])
    assert st[("SYNTH_toolX_m1", "A")] == "complete"
    assert (
        st[("SYNTH_toolX_m1_B", "C")] == "unavailable"
    )  # no canonical tags on SYNTH_B
    assert {st[("SYNTH_toolZ_broken", p)] for p in "ABC"} == {"failed"}
    ws = demo["ws"]
    assert not ws.comparison_dir(
        "SYNTH_toolZ_broken", "A"
    ).exists()  # partial output never promoted
    attempts = ws.attempts("SYNTH_toolZ_broken", "A")
    assert (
        attempts
        and attempts[-1]["status"] == "failed"
        and attempts[-1]["exit_code"] != 0
    )
    record = ws.run_record("SYNTH_toolX_m1", "A")
    assert (
        record["status"] == "complete"
        and record["peak_rss_bytes"]
        and record["implementation"]["fingerprint"]
    )


def test_failed_run_resumes_after_fix(tmp_path):
    config = write_synthetic_demo(tmp_path)
    experiment = load_experiment(config)
    ws = Workspace(experiment)
    runs_mod.run_missing(
        experiment, ws, run_ids=["SYNTH_toolZ_broken"], policy_ids=["A"], log=QUIET
    )
    assert _statuses(experiment, ws)[("SYNTH_toolZ_broken", "A")] == "failed"
    broken = tmp_path / "inputs" / "SYNTH_toolZ_broken.gff3"
    broken.write_text(
        "\n".join(
            l
            for l in broken.read_text().splitlines()
            if "SYNTH_missing_transcript" not in l
        )
        + "\n"
    )
    tasks = runs_mod.run_missing(
        experiment, ws, run_ids=["SYNTH_toolZ_broken"], policy_ids=["A"], log=QUIET
    )
    assert [t.status for t in tasks] == [
        "not_run"
    ]  # new input checksum -> new key, not the failed one
    assert _statuses(experiment, ws)[("SYNTH_toolZ_broken", "A")] == "complete"


def test_cache_reuse_and_invalidation(tmp_path):
    config = write_synthetic_demo(tmp_path)
    experiment = load_experiment(config)
    ws = Workspace(experiment)
    keep = ["SYNTH_toolX_m1", "SYNTH_toolX_m2"]
    runs_mod.run_missing(experiment, ws, run_ids=keep, log=QUIET)
    st = _statuses(experiment, ws)
    assert all(st[(r, p)] == "complete" for r in keep for p in "ABC")
    # 1. adding a predictor through configuration only: existing runs untouched
    raw = json.loads(config.read_text())
    raw["runs"].append(
        {
            **raw["runs"][0],
            "run_id": "SYNTH_toolX_m1_copy",
            "provenance": "SYNTHETIC: added later",
        }
    )
    config.write_text(json.dumps(raw))
    experiment = load_experiment(config)
    st = _statuses(experiment, ws)
    assert all(st[(r, p)] == "complete" for r in keep for p in "ABC")
    assert {st[("SYNTH_toolX_m1_copy", p)] for p in "ABC"} == {"not_run"}
    # 2. changed query input: only that run becomes stale
    q = tmp_path / "inputs" / "SYNTH_toolX_model2.gff3"
    q.write_text(q.read_text().replace("2051\t2300", "2052\t2300"))
    st = _statuses(experiment, ws)
    assert {st[("SYNTH_toolX_m2", p)] for p in "ABC"} == {"stale"}
    assert {st[("SYNTH_toolX_m1", p)] for p in "ABC"} == {"complete"}
    # 3. changed evaluation policy: only that policy becomes stale
    raw["policies"]["C"]["query_transcript_selection"] = "all"
    config.write_text(json.dumps(raw))
    st = _statuses(load_experiment(config), ws)
    assert (
        st[("SYNTH_toolX_m1", "C")] == "stale"
        and st[("SYNTH_toolX_m1", "A")] == "complete"
    )
    # 4. changed reference: every run on it becomes stale
    ref = tmp_path / "inputs" / "SYNTH_reference_A.gff3"
    ref.write_text(
        ref.read_text()
        + "chrA1\tSYNTH\tgene\t18001\t18300\t.\t+\t.\tID=gene:SYNTH_G9;biotype=protein_coding\n"
        "chrA1\tSYNTH\tmRNA\t18001\t18300\t.\t+\t.\tID=transcript:SYNTH_T9;Parent=gene:SYNTH_G9;biotype=protein_coding\n"
        "chrA1\tSYNTH\texon\t18001\t18300\t.\t+\t.\tParent=transcript:SYNTH_T9\n"
    )
    st = _statuses(load_experiment(config), ws)
    assert st[("SYNTH_toolX_m1", "A")] == "stale"


# --------------------------------------------------------------- import
def test_import_matches_by_content_and_separates_pilots(demo, tmp_path):
    ws_src = demo["ws"]
    bench = tmp_path / "old_benchmark"
    # a different directory layout than the workspace: names do not matter, content does
    for run_id, policy, name in (
        ("SYNTH_toolX_m1", "A", "viewA/x1"),
        ("SYNTH_toolX_m1", "B", "viewB/x1"),
        ("SYNTH_toolY_hybrid", "A", "viewA/y"),
    ):
        shutil.copytree(ws_src.comparison_dir(run_id, policy), bench / name)
    record = ws_src.run_record("SYNTH_toolX_m1", "A")
    pilot_argv = [a for a in record["argv"] if a] + ["--region", "chrA1"]
    i = pilot_argv.index("--outdir")
    pilot_argv[i + 1] = str(bench / "pilots" / "x1_chrA1")
    import subprocess

    assert subprocess.run(pilot_argv, capture_output=True).returncode == 0
    snapshot = tmp_path / "code.sha256"
    snapshot.write_text(
        "".join(
            f"{hashlib.sha256((PACKAGE_DIR / f).read_bytes()).hexdigest()}  src/annotation_qc/{f}\n"
            for f in FINGERPRINT_FILES
        )
    )
    raw = json.loads(demo["config"].read_text())
    raw["output_dir"] = str(tmp_path / "fresh_workspace")
    raw["import"] = {"benchmark_dir": str(bench), "code_snapshot": [str(snapshot)]}
    cfg = demo["root"] / "import_test.json"
    cfg.write_text(json.dumps(raw))
    experiment = load_experiment(cfg)
    ws = Workspace(experiment)
    report = import_benchmark(experiment, ws, log=QUIET)
    matched = {(m["run_id"], m["policy_id"]): m["status"] for m in report["matched"]}
    assert matched == {
        ("SYNTH_toolX_m1", "A"): "complete",
        ("SYNTH_toolX_m1", "B"): "complete",
        ("SYNTH_toolY_hybrid", "A"): "complete",
    }
    assert [p["dir"] for p in report["pilots"]] == ["pilots/x1_chrA1"]
    st = _statuses(experiment, ws)
    assert (
        st[("SYNTH_toolX_m1", "A")] == "complete"
        and st[("SYNTH_toolX_m1", "C")] == "not_run"
    )
    # without a code snapshot the implementation cannot be verified: imported but stale
    raw["output_dir"] = str(tmp_path / "fresh2")
    raw["import"] = {"benchmark_dir": str(bench)}
    cfg.write_text(json.dumps(raw))
    exp2 = load_experiment(cfg)
    report2 = import_benchmark(exp2, Workspace(exp2), log=QUIET)
    assert {m["status"] for m in report2["matched"]} == {"stale"}


def test_preparation_reuse_is_verified(demo, tmp_path):
    raw = json.loads(demo["config"].read_text())
    prepared = demo["ws"].preparation_record("SYNTH_A")["annotation"]["output"]["path"]
    wrong = tmp_path / "wrong.gff3"
    wrong.write_text(Path(prepared).read_text().replace("SYNTH_G1", "SYNTH_GX"))
    for path, ok in ((prepared, True), (wrong, False)):
        raw["output_dir"] = str(tmp_path / f"ws_{ok}")
        raw["references"][0]["annotation_preparation"]["reuse_existing"] = str(path)
        cfg = demo["root"] / f"reuse_{ok}.json"
        cfg.write_text(json.dumps(raw))
        exp = load_experiment(cfg)
        if ok:
            rec = prepare_reference(Workspace(exp), exp.references["SYNTH_A"])
            assert rec["annotation"]["method"].startswith(
                "reused existing file; verified"
            )
            assert rec["annotation"]["counts"]["retyped_rows"] == {
                "unconfirmed_transcript": 1
            }
        else:
            with pytest.raises(PreparationError):
                prepare_reference(Workspace(exp), exp.references["SYNTH_A"])


# --------------------------------------------------------------- dataset and queries
def test_identifier_collisions_stay_distinct(con):
    rows = con.execute(
        "SELECT gene_id, original_gene_id, chrom, start, end FROM query_outcomes "
        "WHERE run_id='SYNTH_toolY_hybrid' AND policy_id='A' AND original_gene_id='g1' ORDER BY chrom"
    ).fetchall()
    assert [(r["gene_id"], r["chrom"]) for r in rows] == [
        ("chrA1:g1", "chrA1"),
        ("chrA2:g1", "chrA2"),
    ]
    txs = con.execute(
        "SELECT gene_id, transcript_id, chrom FROM transcripts WHERE source='run:SYNTH_toolY_hybrid' ORDER BY chrom, start"
    ).fetchall()
    assert all(
        t["transcript_id"].startswith(t["chrom"] + ":")
        for t in txs
        if t["gene_id"].endswith(":g1")
    )
    names = con.execute(
        "SELECT gene_id FROM ref_genes WHERE reference_id='SYNTH_A' AND name='ALPHA1'"
    ).fetchall()
    assert len(names) == 2  # same display name, two distinct genes


def test_gene_rows_link_to_locus(con):
    table = queries.genes(con, {"reference": "SYNTH_A", "policy": "A", "limit": "100"})
    assert table["total"] == 6
    for row in table["rows"]:
        loc = queries.locus(
            con,
            {
                "reference": "SYNTH_A",
                "chrom": row["chrom"],
                "start": str(row["start"]),
                "end": str(row["end"]),
                "policy": "A",
                "runs": "SYNTH_toolX_m1",
                "isoforms": "evaluated",
                "gene_id": row["gene_id"],
            },
        )
        genes = {g["gene_id"]: g for g in loc["reference_genes"]}
        assert (
            genes[row["gene_id"]]["start"] == row["start"]
            and genes[row["gene_id"]]["end"] == row["end"]
        )
        evaluated = {
            t["transcript_id"]
            for t in loc["reference_transcripts"]
            if t["gene_id"] == row["gene_id"] and t["evaluated"]
        }
        assert len(evaluated) == row["evaluated_transcript_count"]
    detail = queries.gene_detail(con, "SYNTH_A", "A", "SYNTH_G5")
    assert detail["outcomes"]["SYNTH_toolX_m1"]["cds_status"] == "terminal_diff"
    assert (
        detail["outcomes"]["SYNTH_toolX_m1"]["explanation"]["stop_side_diff_bp"] == -50
    )
    assert detail["outcomes"]["SYNTH_toolX_m2"]["counterpart_count"] == 2  # split
    assert queries.locus(
        con, {"reference": "SYNTH_A", "chrom": "chrA1", "start": "1", "end": "9000000"}
    )["window_clamped"]


def test_filters_and_policy_semantics(con):
    q = lambda **p: queries.genes(con, {"reference": "SYNTH_A", "policy": "A", **p})
    assert {
        r["gene_id"] for r in q(focus_run="SYNTH_toolX_m1", status="missed")["rows"]
    } == {"SYNTH_G3"}
    assert {r["gene_id"] for r in q(focus_run="SYNTH_toolX_m2", split="1")["rows"]} == {
        "SYNTH_G5"
    }
    contrast = q(
        focus_run="SYNTH_toolX_m1",
        compare_run="SYNTH_toolX_m2",
        contrast="exact_in_focus",
    )["rows"]
    assert {r["gene_id"] for r in contrast} == {"SYNTH_G2", "SYNTH_G7"}
    assert q(q="beta2")["rows"][0]["gene_id"] == "SYNTH_G2"
    # alternative isoform: exact under A, not under B (longest CDS) for model1
    exact = lambda policy, run: con.execute(
        "SELECT cds_exact FROM ref_outcomes WHERE policy_id=? AND run_id=? AND gene_id='SYNTH_G1'",
        (policy, run),
    ).fetchone()[0]
    assert (
        exact("A", "SYNTH_toolX_m1"),
        exact("B", "SYNTH_toolX_m1"),
        exact("C", "SYNTH_toolX_m1"),
    ) == (1, 0, 1)
    assert (
        exact("A", "SYNTH_toolX_m2"),
        exact("B", "SYNTH_toolX_m2"),
        exact("C", "SYNTH_toolX_m2"),
    ) == (1, 1, 0)
    preds = queries.query_genes(
        con,
        {
            "reference": "SYNTH_A",
            "policy": "A",
            "run": "SYNTH_toolX_m1",
            "classification": "Novel",
        },
    )
    assert {(r["gene_id"], r["novel_category"]) for r in preds["rows"]} == {
        ("SYNTH_x1_ps", "Overlaps_reference_pseudogene"),
        ("SYNTH_x1_new", "Novel_no_reference_overlap"),
    }


def test_metric_rows_and_missing_values(con, demo):
    m = {
        (r["run_id"], r["policy_id"], r["metric_id"]): r
        for r in con.execute("SELECT * FROM metrics")
    }
    exact = m[("SYNTH_toolX_m1", "A", "cds_exact_ref")]
    assert (
        exact["numerator"],
        exact["denominator"],
        exact["status"],
        exact["source"],
    ) == (3, 6, "ok", "comparator")
    assert (
        m[("SYNTH_toolX_m1", "A", "cds_intron_chain_any_pair")]["source"]
        == "independent"
    )
    assert m[("SYNTH_toolX_m1", "A", "cds_exact_ref_multi")]["denominator"] == 5
    assert not any(
        k[0] == "SYNTH_toolZ_broken" for k in m
    )  # failed: no metrics, never zeros
    results = {
        (r["run_id"], r["policy_id"]): r["status"]
        for r in con.execute("SELECT * FROM results")
    }
    assert (
        results[("SYNTH_toolZ_broken", "A")] == "failed"
        and results[("SYNTH_toolX_m1_B", "C")] == "unavailable"
    )
    # summary-only result: per-gene and segment metrics are not_available, not zero
    summary_only = demo["root"] / "summary_only"
    summary_only.mkdir(exist_ok=True)
    shutil.copy(
        demo["ws"].comparison_dir("SYNTH_toolX_m1", "A") / "comparison_summary.json",
        summary_only,
    )
    rows, ref, qry = derive_metrics(summary_only)
    by_id = {r["metric_id"]: r for r in rows}
    assert ref is None and by_id["cds_exact_ref"]["status"] == "ok"
    assert (
        by_id["exact_f1_one_to_one"]["status"] == "not_available"
        and by_id["exact_f1_one_to_one"]["value"] is None
    )


def test_incremental_build_skips_unchanged(demo):
    db = demo["ws"].dashboard_db()
    before = dict(
        sqlite3.connect(db)
        .execute("SELECT component, built FROM components")
        .fetchall()
    )
    build_dashboard(demo["experiment"], demo["ws"], log=QUIET)
    after = dict(
        sqlite3.connect(db)
        .execute("SELECT component, built FROM components")
        .fetchall()
    )
    assert before == after


def test_http_api(demo):
    from http.server import ThreadingHTTPServer

    from ensembl.genes.annotation_qc.dashboard.server import Datasets, make_handler

    httpd = ThreadingHTTPServer(
        ("127.0.0.1", 0),
        make_handler(Datasets({"SYNTHETIC_demo": demo["ws"].dashboard_db()})),
    )
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}"
    try:
        get = lambda path: urllib.request.urlopen(base + path).read().decode()
        assert json.loads(get("/api/experiments"))[0]["synthetic"] is True
        assert "<title>Annotation QC</title>" in get("/")
        tsv = get(
            "/api/genes?reference=SYNTH_A&policy=A&format=csv&runs=SYNTH_toolX_m1,SYNTH_toolZ_broken"
        )
        header, *rows = tsv.strip().split("\n")
        assert len(rows) == 6 and "SYNTH_toolX_m1.cds_status" in header
        assert all(
            line.split("\t")[header.split("\t").index("SYNTH_toolZ_broken.cds_status")]
            == "NA"
            for line in rows
        )
        meta = json.loads(get("/api/meta"))
        assert meta["experiment"]["synthetic"] and "metric_dictionary" in meta
        with pytest.raises(urllib.error.HTTPError) as error:
            urllib.request.urlopen(base + "/static/../server.py")
        assert error.value.code == 404
    finally:
        httpd.shutdown()

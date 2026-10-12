"""Benchmark/dashboard integration with comparator one-to-one metrics and stop-codon harmonisation (SYNTHETIC fixtures)."""

import hashlib
import json
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from ensembl.genes.annotation_qc.benchmark import runs as runs_mod
from ensembl.genes.annotation_qc.benchmark.config import (
    ConfigError,
    comparator_parameters,
    load_experiment,
)
from ensembl.genes.annotation_qc.benchmark.derive import (
    check_pairs,
    derive_metrics,
    one_to_one_metrics,
)
from ensembl.genes.annotation_qc.benchmark.importer import import_benchmark
from ensembl.genes.annotation_qc.benchmark.provenance import (
    FINGERPRINT_FILES,
    PACKAGE_DIR,
    _has_code,
    implementation_from_snapshot,
    package_versions,
)
from ensembl.genes.annotation_qc.benchmark.store import build_dashboard
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
from ensembl.genes.annotation_qc.dashboard import queries

QUIET = lambda *a, **k: None  # noqa: E731

# chr1 (600 bp of C): ATG at 101, CDS 101-130 + 201-233 with TAA at 231-233.
# The reference CDS includes the stop; the CDS-only query ends at 230.
GENOME = "C" * 100 + "ATG" + "GCC" * 9 + "C" * 70 + "GCC" * 10 + "TAA" + "C" * 367
REFERENCE = """##gff-version 3
##sequence-region chr1 1 600
chr1\tSYNTH\tgene\t101\t233\t.\t+\t.\tID=gene:SYNTH_RG1;Name=SYNTH_RG1;biotype=protein_coding;gene_id=SYNTH_RG1
chr1\tSYNTH\tmRNA\t101\t233\t.\t+\t.\tID=transcript:SYNTH_RT1;Parent=gene:SYNTH_RG1;biotype=protein_coding;tag=Ensembl_canonical;transcript_id=SYNTH_RT1
chr1\tSYNTH\texon\t101\t130\t.\t+\t.\tParent=transcript:SYNTH_RT1
chr1\tSYNTH\texon\t201\t233\t.\t+\t.\tParent=transcript:SYNTH_RT1
chr1\tSYNTH\tCDS\t101\t130\t.\t+\t0\tID=CDS:SYNTH_RT1;Parent=transcript:SYNTH_RT1
chr1\tSYNTH\tCDS\t201\t233\t.\t+\t0\tID=CDS:SYNTH_RT1;Parent=transcript:SYNTH_RT1
"""
QUERY = """chr1\tSYNTH_pred\tCDS\t101\t130\t.\t+\t0\tgene_id "SYNTH_q1"; transcript_id "SYNTH_q1.t1";
chr1\tSYNTH_pred\tCDS\t201\t230\t.\t+\t0\tgene_id "SYNTH_q1"; transcript_id "SYNTH_q1.t1";
"""


def _config(root: Path, runs: list[dict], output="workspace", **extra) -> Path:
    inputs = root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    (inputs / "genome.fa").write_text(
        ">chr1\n" + "\n".join(GENOME[i : i + 60] for i in range(0, 600, 60)) + "\n"
    )
    (inputs / "reference.gff3").write_text(REFERENCE)
    (inputs / "query.gtf").write_text(QUERY)
    config = {
        "config_version": 1,
        "experiment": {"id": "SYNTH_stop", "synthetic": True},
        "output_dir": output,
        "evaluation": {
            "evaluation_mode": "protein_coding",
            "reference_transcript_biotypes": ["protein_coding"],
        },
        "policies": {
            "A": {
                "reference_transcript_selection": "all",
                "query_transcript_selection": "all",
            }
        },
        "references": [
            {
                "id": "SYNTH_R",
                "species": "synthetic",
                "assembly": "SYNTH_v1",
                "annotation": "inputs/reference.gff3",
                "genome": "inputs/genome.fa",
            }
        ],
        "runs": runs,
        **extra,
    }
    path = root / f"{output}.json"
    path.write_text(json.dumps(config))
    return path


ORIGINAL = {
    "run_id": "SYNTH_pred_original",
    "reference": "SYNTH_R",
    "tool": "SYNTH predictor",
    "annotation": "inputs/query.gtf",
}
HARMONISED = {
    **ORIGINAL,
    "run_id": "SYNTH_pred_stop",
    "label": "SYNTH predictor — stop harmonised",
    "role": "supplementary",
    "comparator_options": {"add_stop_codon": True, "stop_codon_genetic_code": 1},
}


def _snapshot(path: Path) -> Path:
    path.write_text(
        "".join(
            f"{hashlib.sha256((PACKAGE_DIR / f).read_bytes()).hexdigest()}  src/annotation_qc/{f}\n"
            for f in FINGERPRINT_FILES
        )
    )
    return path


@pytest.fixture(scope="module")
def computed(tmp_path_factory):
    """Both runs computed by `benchmark run`, then built into a dataset."""
    root = tmp_path_factory.mktemp("stop")
    experiment = load_experiment(_config(root, [ORIGINAL, HARMONISED]))
    ws = Workspace(experiment)
    tasks = runs_mod.run_missing(experiment, ws, log=QUIET)
    build_dashboard(experiment, ws, log=QUIET)
    return {"root": root, "experiment": experiment, "ws": ws, "tasks": tasks}


# --------------------------------------------------------------- implementation fingerprint
def test_fingerprint_covers_every_module_the_runner_loads():
    assert {"metrics/pairwise/matching.py", "metrics/pairwise/stop_codons.py"} <= set(
        FINGERPRINT_FILES
    )
    code = (
        "import sys, ensembl.genes.annotation_qc.runners.pairwise_compare, "
        "ensembl.genes.annotation_qc.reports.pairwise_plots; "
        "print('\\n'.join(m.__file__ for n, m in sys.modules.items() "
        "if n.startswith('ensembl.genes.annotation_qc') and getattr(m, '__file__', None)))"
    )
    loaded = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    ).stdout.split()
    rel = {str(Path(f).resolve().relative_to(PACKAGE_DIR)) for f in loaded}
    # every loaded module with code is fingerprinted; the rest are empty package markers
    assert {r for r in rel if _has_code(r)} == set(FINGERPRINT_FILES)
    assert all(r.endswith("__init__.py") for r in rel - set(FINGERPRINT_FILES))
    assert (
        "cli.py" not in FINGERPRINT_FILES
    )  # dispatch only: new commands do not invalidate results


def test_snapshot_without_a_comparator_file_cannot_verify(tmp_path):
    full = _snapshot(tmp_path / "full.sha256")
    assert implementation_from_snapshot([full], package_versions())["fingerprint"]
    partial = tmp_path / "old.sha256"
    partial.write_text(
        "".join(
            l + "\n" for l in full.read_text().splitlines() if "stop_codons" not in l
        )
    )
    old = implementation_from_snapshot([partial], package_versions())
    assert old["fingerprint"] is None and old["missing_from_snapshot"] == [
        "metrics/pairwise/stop_codons.py"
    ]


# --------------------------------------------------------------- configuration and commands
def test_stop_codon_options_parameters_and_command(tmp_path):
    raw = {
        **HARMONISED,
        "comparator_options": {
            "add_stop_codon": True,
            "stop_codon_genetic_code": 11,
            "stop_codon_sequence_code": {"chrM": 2, "chrP": "none"},
        },
    }
    exp = load_experiment(_config(tmp_path, [ORIGINAL, raw]))
    ref, policy = exp.references["SYNTH_R"], exp.policies["A"]
    original = comparator_parameters(exp, ref, exp.runs["SYNTH_pred_original"], policy)
    harmonised = comparator_parameters(exp, ref, exp.runs["SYNTH_pred_stop"], policy)
    # original-coordinate parameters are unchanged by the new option (cache keys stay valid)
    assert not any(k.startswith(("add_stop", "stop_codon")) for k in original)
    assert harmonised["stop_codon_sequence_code"] == ["chrM=2", "chrP=none"]
    argv = runs_mod.build_argv(
        {"query_path": "q", "reference_path": "r", "genome_path": "g"}, harmonised
    )
    assert argv[argv.index("--add-stop-codon") :] == [
        "--add-stop-codon",
        "--stop-codon-genetic-code",
        "11",
        "--stop-codon-sequence-code",
        "chrM=2",
        "--stop-codon-sequence-code",
        "chrP=none",
    ]
    assert "--add-stop-codon" not in runs_mod.build_argv(
        {"query_path": "q", "reference_path": "r"}, original
    )
    assert (
        exp.runs["SYNTH_pred_stop"].display_label == "SYNTH predictor — stop harmonised"
    )
    assert exp.runs["SYNTH_pred_stop"].role == "supplementary"


@pytest.mark.parametrize(
    "options, message",
    [
        ({"stop_codon_genetic_code": 2}, "need add_stop_codon"),
        (
            {"add_stop_codon": True, "stop_codon_genetic_code": 7},
            "Unsupported genetic code 7",
        ),
        (
            {"add_stop_codon": True, "stop_codon_sequence_code": ["chrM"]},
            "Expected SEQNAME=TABLE",
        ),
        ({"add_stop_codon": "yes"}, "must be true or false"),
        ({"stop_codons": True}, "unknown option"),
    ],
)
def test_invalid_stop_codon_options(tmp_path, options, message):
    with pytest.raises(ConfigError) as error:
        load_experiment(
            _config(tmp_path, [{**ORIGINAL, "comparator_options": options}])
        )
    assert any(message in p for p in error.value.problems)


def test_stop_codon_option_needs_a_genome(tmp_path):
    path = _config(tmp_path, [HARMONISED])
    raw = json.loads(path.read_text())
    del raw["references"][0]["genome"]
    path.write_text(json.dumps(raw))
    with pytest.raises(ConfigError) as error:
        load_experiment(path)
    assert any("add_stop_codon needs a genome" in p for p in error.value.problems)


# --------------------------------------------------------------- computed runs, metrics and dataset
def test_original_and_harmonised_are_distinct_results(computed):
    tasks = {t.run.run_id: t for t in computed["tasks"]}
    assert tasks["SYNTH_pred_original"].key != tasks["SYNTH_pred_stop"].key
    assert runs_mod.plan(computed["experiment"], computed["ws"])[0].status == "complete"
    con = queries.connect(computed["ws"].dashboard_db())
    m = {
        (r["run_id"], r["metric_id"]): r
        for r in queries.overview(con, "SYNTH_R")["metrics"]
    }
    assert (
        m[("SYNTH_pred_original", "exact_tp_one_to_one")]["numerator"],
        m[("SYNTH_pred_stop", "exact_tp_one_to_one")]["numerator"],
    ) == (0, 1)
    f1 = m[("SYNTH_pred_stop", "exact_f1_one_to_one")]
    assert (f1["numerator"], f1["denominator"], f1["value"], f1["source"]) == (
        2,
        2,
        1.0,
        "comparator",
    )
    zero = m[("SYNTH_pred_original", "exact_f1_one_to_one")]
    assert (zero["status"], zero["value"]) == ("ok", 0.0)  # a real zero, not missing
    assert (
        m[("SYNTH_pred_original", "stop_codon_extended")]["status"] == "not_applicable"
    )
    assert (
        m[("SYNTH_pred_stop", "stop_codon_extended")]["numerator"],
        m[("SYNTH_pred_stop", "stop_codon_extended")]["denominator"],
    ) == (1, 1)
    locus = m[("SYNTH_pred_stop", "locus_recovery_cds_overlap")]
    assert (locus["numerator"], locus["source"]) == (1, "comparator")


def test_harmonised_models_show_evaluated_coordinates(computed):
    con = queries.connect(computed["ws"].dashboard_db())
    cds = {
        r["source"]: r["cds"]
        for r in con.execute(
            "SELECT source, cds FROM transcripts WHERE source LIKE 'run:%'"
        )
    }
    assert cds == {
        "run:SYNTH_pred_original": "101-130,201-230",
        "run:SYNTH_pred_stop": "101-130,201-233",
    }
    runs = {r["run_id"]: r for r in queries.meta(con)["runs"]}
    assert (
        runs["SYNTH_pred_stop"]["models_source"]["kind"]
        == "comparator_evaluated_stop_harmonised"
    )
    assert runs["SYNTH_pred_original"]["models_source"]["kind"] == "prediction_file"
    loc = queries.locus(
        con,
        {
            "reference": "SYNTH_R",
            "chrom": "chr1",
            "start": "1",
            "end": "600",
            "policy": "A",
            "runs": "SYNTH_pred_stop,SYNTH_pred_original",
        },
    )
    stop = loc["runs"]["SYNTH_pred_stop"]["transcripts"][0]["stop_codon"]
    assert (
        stop["status"],
        stop["reason"],
        stop["next_codon"],
        stop["terminal_cds_end"],
        stop["new_terminal_cds_end"],
    ) == ("extended", "stop_codon_added", "TAA", 230, 233)
    assert loc["runs"]["SYNTH_pred_original"]["transcripts"][0]["stop_codon"] is None
    detail = queries.gene_detail(con, "SYNTH_R", "A", "SYNTH_RG1")["outcomes"]
    assert detail["SYNTH_pred_stop"]["one_to_one_partner"] == "SYNTH_q1"
    assert detail["SYNTH_pred_original"]["one_to_one_partner"] is None
    pairs = queries.pairs(
        con, {"reference": "SYNTH_R", "policy": "A", "run": "SYNTH_pred_stop"}
    )
    assert pairs["total"] == 1 and pairs["rows"][0]["query_gene_id"] == "SYNTH_q1"
    audit = queries.stop_audit(con, {"run": "SYNTH_pred_stop"})
    assert audit["counts"] == [
        {"status": "extended", "reason": "stop_codon_added", "transcripts": 1}
    ]


def test_harmonised_run_without_evaluated_annotation_is_not_substituted(
    computed, tmp_path
):
    root = tmp_path / "copy"
    shutil.copytree(computed["root"], root)
    experiment = load_experiment(root / "workspace.json")
    ws = Workspace(experiment)
    record = ws.run_record("SYNTH_pred_stop", "A")
    out = ws.comparison_dir("SYNTH_pred_stop", "A")
    ws.write_run_record("SYNTH_pred_stop", "A", {**record, "output_dir": str(out)})
    (out / "query_evaluated.stop_harmonised.gtf").unlink()
    task = next(t for t in runs_mod.plan(experiment, ws, run_ids=["SYNTH_pred_stop"]))
    assert (task.status, task.reason) == (
        "stale",
        "outputs missing: query_evaluated.stop_harmonised.gtf",
    )
    build_dashboard(experiment, ws, log=QUIET)
    con = queries.connect(ws.dashboard_db())
    assert (
        con.execute(
            "SELECT COUNT(*) FROM transcripts WHERE source='run:SYNTH_pred_stop'"
        ).fetchone()[0]
        == 0
    )
    assert (
        "original annotation is not shown"
        in queries.meta(con)["model_errors"]["SYNTH_pred_stop"]
    )


# --------------------------------------------------------------- import
def test_import_matches_stop_codon_settings_and_flags_incompatible(computed, tmp_path):
    bench = tmp_path / "bench"
    for run_id, name in (("SYNTH_pred_original", "orig"), ("SYNTH_pred_stop", "harm")):
        shutil.copytree(computed["ws"].comparison_dir(run_id, "A"), bench / name)
    snapshot = _snapshot(tmp_path / "code.sha256")
    imports = {
        "import": {"benchmark_dir": str(bench), "code_snapshot": [str(snapshot)]}
    }
    exp = load_experiment(
        _config(computed["root"], [ORIGINAL, HARMONISED], output="import_ws", **imports)
    )
    report = import_benchmark(exp, Workspace(exp), log=QUIET)
    assert {(m["dir"], m["run_id"], m["status"]) for m in report["matched"]} == {
        ("orig", "SYNTH_pred_original", "complete"),
        ("harm", "SYNTH_pred_stop", "complete"),
    }
    # only the original run configured: the harmonised manifest is not mistaken for it
    exp = load_experiment(
        _config(computed["root"], [ORIGINAL], output="import_ws2", **imports)
    )
    report = import_benchmark(exp, Workspace(exp), log=QUIET)
    assert [m["dir"] for m in report["matched"]] == ["orig"]
    assert (
        report["unmatched"][0]["dir"] == "harm"
        and "stop-codon settings" in report["unmatched"][0]["reason"]
    )
    # a summary contradicting its manifest is imported as stale, with the reason
    summary = bench / "harm" / "comparison_summary.json"
    data = json.loads(summary.read_text())
    data["query_stop_codon_harmonisation"] = {"applied": False, "annotation": "query"}
    summary.write_text(json.dumps(data))
    exp = load_experiment(
        _config(
            computed["root"], [ORIGINAL, HARMONISED], output="import_ws3", **imports
        )
    )
    report = import_benchmark(exp, Workspace(exp), log=QUIET)
    harm = next(m for m in report["matched"] if m["dir"] == "harm")
    assert (
        harm["status"] == "stale" and "records no harmonisation" in harm["difference"]
    )


# --------------------------------------------------------------- missing, undefined, zero, legacy
def test_one_to_one_undefined_missing_and_inconsistent():
    block = lambda tp, r, q, **kw: {
        "cds_exact_one_to_one": {
            "true_positives": tp,
            "reference_genes": r,
            "query_genes": q,
            "unmatched_reference_genes": r - tp,
            "unmatched_query_genes": q - tp,
            "exact_gene_pairs": tp,
            "exact_components": tp,
            "exact_components_with_several_edges": 0,
            "precision": round(tp / q, 4) if q else None,
            "recall": round(tp / r, 4) if r else None,
            "f1": round(2 * tp / (r + q), 4) if r + q else None,
            **kw,
        }
    }
    rows = {m["metric_id"]: m for m in one_to_one_metrics(block(0, 5, 0), 5, 0)}
    assert rows["exact_precision_one_to_one"]["status"] == "undefined"  # Q = 0
    assert (
        rows["exact_recall_one_to_one"]["status"],
        rows["exact_recall_one_to_one"]["value"],
    ) == ("ok", 0.0)
    rows = {m["metric_id"]: m for m in one_to_one_metrics({}, 5, 4)}
    assert {m["status"] for m in rows.values()} == {"not_available"}
    rows = {m["metric_id"]: m for m in one_to_one_metrics(block(2, 5, 4, f1=0.9), 5, 4)}
    assert (
        rows["exact_f1_one_to_one"]["status"] == "inconsistent"
        and "f1 0.9" in rows["exact_f1_one_to_one"]["note"]
    )
    import pandas as pd

    dup = pd.DataFrame({"reference_gene_id": ["a", "b"], "query_gene_id": ["q", "q"]})
    assert check_pairs(dup, block(2, 5, 4))["status"] == "disagrees"


def test_legacy_dataset_is_relabelled_not_rewritten(computed, tmp_path):
    db = tmp_path / "legacy.sqlite"
    shutil.copy(computed["ws"].dashboard_db(), db)
    with sqlite3.connect(db) as raw:
        raw.execute("UPDATE meta SET value='\"1\"' WHERE key='store_version'")
        raw.execute("DELETE FROM metrics WHERE metric_id LIKE 'exact_%one_to_one'")
        raw.execute(
            "UPDATE metrics SET metric_id='exact_tp_one_to_one' WHERE metric_id='historical_exact_tp_best_pair'"
        )
        raw.execute("DROP TABLE exact_pairs")
    con = queries.connect(db)
    m = {
        (r["run_id"], r["metric_id"]): r
        for r in queries.overview(con, "SYNTH_R")["metrics"]
    }
    assert m[("SYNTH_pred_stop", "historical_exact_tp_best_pair")]["numerator"] == 1
    assert m[("SYNTH_pred_stop", "exact_tp_one_to_one")]["status"] == "not_available"
    assert queries.meta(con)["legacy_dataset"]
    with pytest.raises(queries.QueryError):
        queries.pairs(
            con, {"reference": "SYNTH_R", "policy": "A", "run": "SYNTH_pred_stop"}
        )
    assert (
        sqlite3.connect(db)
        .execute("SELECT value FROM meta WHERE key='store_version'")
        .fetchone()[0]
        == '"1"'
    )


def test_derived_rows_for_a_summary_without_newer_blocks(computed, tmp_path):
    out = tmp_path / "old"
    shutil.copytree(computed["ws"].comparison_dir("SYNTH_pred_original", "A"), out)
    data = json.loads((out / "comparison_summary.json").read_text())
    for key in (
        "cds_exact_one_to_one",
        "cds_overlap_locus_recovery",
        "query_stop_codon_harmonisation",
    ):
        data.pop(key)
    (out / "comparison_summary.json").write_text(json.dumps(data))
    rows = {r["metric_id"]: r for r in derive_metrics(out)[0]}
    assert rows["exact_f1_one_to_one"]["status"] == "not_available"
    assert (
        rows["historical_exact_f1_best_pair"]["status"] == "ok"
    )  # kept, separately labelled
    assert rows["stop_codon_extended"]["status"] == "not_available"
    assert (
        rows["locus_recovery_cds_overlap"]["source"],
        rows["locus_recovery_cds_overlap"]["numerator"],
    ) == ("comparator_per_gene", 1)

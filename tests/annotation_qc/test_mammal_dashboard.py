"""
Species-level dashboard behaviour (SYNTHETIC fixtures): display roles, one reference per
experiment, gene-table contrasts and the experiment list. Presentation only: nothing here
may change a stored comparator value.
"""

import argparse
import json
import threading
import urllib.request

import pytest

from ensembl.genes.annotation_qc.benchmark import store
from ensembl.genes.annotation_qc.benchmark.config import ConfigError, load_experiment
from ensembl.genes.annotation_qc.dashboard import cli as dashboard_cli
from ensembl.genes.annotation_qc.dashboard import queries

from .test_benchmark_dashboard import con, demo  # noqa: F401  (module fixtures)


def _config(tmp_path, runs, references=None):
    raw = {
        "config_version": 2,
        "experiment": {"id": "sp"},
        "output_dir": "out",
        "references": references
        or [{"id": "r", "species": "cat (Felis catus)", "assembly": "a", "annotation": "ref.gff3"}],
        "runs": runs,
    }
    path = tmp_path / "c.json"
    path.write_text(json.dumps(raw))
    return path


def _run(run_id, tool, **extra):
    return {"run_id": run_id, "reference": "r", "tool": tool, "annotation": f"{run_id}.gtf", **extra}


def test_display_roles_reserve_colours_and_need_no_fixed_names(tmp_path):
    runs = [
        _run("x_annevo", "ANNEVO"),
        _run("x_tiberius", "Tiberius", display_role="primary_comparator"),
        _run("x_internal", "Internal predictions", display_role="focus"),
        _run("x_helixer", "Helixer"),
    ]
    exp = load_experiment(_config(tmp_path, runs))
    colors = store.run_colors(exp)
    assert colors["x_internal"] == store.ROLE_COLORS["focus"]
    assert colors["x_tiberius"] == store.ROLE_COLORS["primary_comparator"]
    others = {colors["x_annevo"], colors["x_helixer"]}
    assert len(others) == 2 and not others & set(store.ROLE_COLORS.values())
    meta = {r["run_id"]: r for r in store._run_meta(exp)}  # pylint: disable=protected-access
    assert meta["x_internal"]["display_role"] == "focus"
    assert meta["x_annevo"]["display_role"] is None
    # roles belong to whatever runs carry them: a species without the focus tool still loads
    exp = load_experiment(_config(tmp_path, [_run("y_helixer", "Helixer"), runs[1]]))
    assert store.run_colors(exp)["x_tiberius"] == store.ROLE_COLORS["primary_comparator"]


def test_display_role_is_validated_and_not_a_comparator_parameter(tmp_path):
    with pytest.raises(ConfigError) as error:
        load_experiment(_config(tmp_path, [_run("x", "T", display_role="hero")]))
    assert any("display_role" in p for p in error.value.problems)
    from ensembl.genes.annotation_qc.benchmark.config import comparator_parameters

    plain = load_experiment(_config(tmp_path, [_run("x", "T")]))
    roled = load_experiment(_config(tmp_path, [_run("x", "T", display_role="focus")]))
    args = lambda e: comparator_parameters(  # noqa: E731
        e, e.references["r"], e.runs["x"], next(iter(e.policies.values()))
    )
    assert args(plain) == args(roled)


def test_dashboard_command_needs_one_reference_per_experiment(tmp_path):
    refs = [
        {"id": "r", "species": "human", "assembly": "a", "annotation": "h.gff3"},
        {"id": "m", "species": "mouse", "assembly": "b", "annotation": "m.gff3"},
    ]
    path = _config(tmp_path, [_run("x", "T")], refs)
    args = argparse.Namespace(config=[str(path)], host="127.0.0.1", port=0, verbose=False)
    with pytest.raises(SystemExit) as error:
        dashboard_cli.cmd_dashboard(args)
    assert "2 references" in str(error.value) and "one configuration per reference" in str(error.value)


def _exact_sets(con, a, b):  # noqa: F811
    rows = con.execute(
        "SELECT gene_id, run_id, cds_exact FROM ref_outcomes WHERE reference_id='SYNTH_A' "
        "AND policy_id='A' AND run_id IN (?,?)",
        (a, b),
    ).fetchall()
    by = {}
    for gene, run, exact in rows:
        by.setdefault(gene, {})[run] = exact
    both = {g for g, v in by.items() if v.get(a) == 1 and v.get(b) == 1}
    neither = {g for g, v in by.items() if v.get(a) == 0 and v.get(b) == 0}
    return both, neither


def test_gene_contrasts_both_and_neither(con):  # noqa: F811
    a, b = "SYNTH_toolX_m1", "SYNTH_toolX_m2"
    q = lambda c: {  # noqa: E731
        r["gene_id"]
        for r in queries.genes(
            con,
            {"reference": "SYNTH_A", "policy": "A", "focus_run": a, "compare_run": b, "contrast": c},
        )["rows"]
    }
    both, neither = _exact_sets(con, a, b)
    assert q("exact_in_both") == both and both
    assert q("exact_in_neither") == neither and neither
    assert not both & neither
    with pytest.raises(queries.QueryError):
        q("best_overall")


def test_experiment_list_names_the_species(demo):  # noqa: F811
    from http.server import ThreadingHTTPServer

    from ensembl.genes.annotation_qc.dashboard.server import Datasets, make_handler

    httpd = ThreadingHTTPServer(
        ("127.0.0.1", 0), make_handler(Datasets({"SYNTHETIC_demo": demo["ws"].dashboard_db()}))
    )
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    try:
        url = f"http://127.0.0.1:{httpd.server_address[1]}/api/experiments"
        (entry,) = json.loads(urllib.request.urlopen(url).read())
        first = demo["experiment"].references[next(iter(demo["experiment"].references))]
        assert entry["species"] == first.species and entry["assembly"] == first.assembly
    finally:
        httpd.shutdown()


def test_contract_results_carry_only_comparator_metric_rows(con):  # noqa: F811
    """Every metric row of a contract result comes from metrics.tsv: no dashboard-side estimates."""
    sources = {
        source
        for (source,) in con.execute(
            "SELECT DISTINCT m.source FROM metrics m JOIN results r ON r.reference_id=m.reference_id "
            "AND r.run_id=m.run_id AND r.policy_id=m.policy_id WHERE r.status IN ('complete', 'unverified')"
        )
    }
    assert sources and "historical" not in sources
    assert not con.execute(
        "SELECT COUNT(*) FROM metrics m JOIN results r ON r.reference_id=m.reference_id AND "
        "r.run_id=m.run_id AND r.policy_id=m.policy_id WHERE r.status='complete' AND m.metric_id LIKE 'historical_%'"
    ).fetchone()[0]

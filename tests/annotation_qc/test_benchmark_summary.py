"""Sequence scopes, run display order and splice tables read from comparator outputs (SYNTHETIC fixtures)."""

import json
from pathlib import Path

import pytest

from ensembl.genes.annotation_qc.benchmark import runs as runs_mod
from ensembl.genes.annotation_qc.benchmark.config import ConfigError, load_experiment
from ensembl.genes.annotation_qc.benchmark.store import build_dashboard
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
from ensembl.genes.annotation_qc.dashboard import queries

QUIET = lambda *a, **k: None  # noqa: E731
THREE = ((101, 150), (201, 260), (301, 350))


# chr1: SYNTH_G1 (+) three CDS segments. chr2: SYNTH_G2 (-) three segments and SYNTH_G3 (-)
# with no prediction. The query reproduces SYNTH_G1 exactly and SYNTH_G2 with one donor shifted.
def _ref_gene(chrom, strand, gid, segments):
    lo, hi = segments[0][0], segments[-1][1]
    lines = [
        f"{chrom}\tSYNTH\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID=gene:{gid};Name={gid};biotype=protein_coding;gene_id={gid}",
        f"{chrom}\tSYNTH\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID=transcript:{gid}.t;Parent=gene:{gid};biotype=protein_coding;transcript_id={gid}.t",
    ]
    for s, e in segments:
        lines.append(
            f"{chrom}\tSYNTH\texon\t{s}\t{e}\t.\t{strand}\t.\tParent=transcript:{gid}.t"
        )
        lines.append(
            f"{chrom}\tSYNTH\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=CDS:{gid}.t;Parent=transcript:{gid}.t"
        )
    return "\n".join(lines)


def _query_gene(chrom, strand, gid, segments):
    return "\n".join(
        f'{chrom}\tSYNTH_pred\t{feature}\t{s}\t{e}\t.\t{strand}\t{"0" if feature == "CDS" else "."}\tgene_id "{gid}"; transcript_id "{gid}.t1";'
        for s, e in segments
        for feature in ("exon", "CDS")
    )


REFERENCE = (
    "\n".join(
        [
            "##gff-version 3",
            "##sequence-region chr1 1 1000",
            "##sequence-region chr2 1 1000",
            _ref_gene("chr1", "+", "SYNTH_G1", THREE),
            _ref_gene("chr2", "-", "SYNTH_G2", THREE),
            _ref_gene("chr2", "-", "SYNTH_G3", ((601, 650), (701, 760))),
        ]
    )
    + "\n"
)
QUERY = (
    "\n".join(
        [
            _query_gene("chr1", "+", "SYNTH_q1", THREE),
            _query_gene("chr2", "-", "SYNTH_q2", ((101, 150), (205, 260), (301, 350))),
        ]
    )
    + "\n"
)
# a second tool that reproduces only SYNTH_G1
QUERY_B = _query_gene("chr1", "+", "SYNTH_b1", THREE) + "\n"


def _write(root: Path, **extra) -> Path:
    inputs = root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    (inputs / "reference.gff3").write_text(REFERENCE)
    (inputs / "query.gtf").write_text(QUERY)
    (inputs / "query_b.gtf").write_text(QUERY_B)
    (inputs / "query_a2.gtf").write_text("# SYNTHETIC run 2\n" + QUERY)
    config = {
        "config_version": 2,
        "experiment": {"id": "SYNTH_scope", "synthetic": True},
        "output_dir": "workspace",
        "evaluation": {
            "evaluation_mode": "protein_coding",
            "reference_transcript_biotypes": ["protein_coding"],
        },
        "policies": {
            "A": {
                "reference_transcript_selection": "all",
                "query_transcript_selection": "all",
            },
            "B": {
                "reference_transcript_selection": "longest_cds",
                "query_transcript_selection": "longest_cds",
            },
        },
        "references": [
            {
                "id": "SYNTH_R",
                "species": "synthetic",
                "assembly": "SYNTH_v1",
                "annotation": "inputs/reference.gff3",
                "scope": {"type": "all", "label": "Whole synthetic assembly"},
                "analysis_scopes": [
                    {"id": "SYNTH_chr1", "label": "chr1 only", "sequences": ["chr1"]}
                ],
            }
        ],
        "runs": [
            {
                "run_id": "SYNTH_b",
                "reference": "SYNTH_R",
                "tool": "SYNTH B",
                "annotation": "inputs/query_b.gtf",
            },
            {
                "run_id": "SYNTH_a",
                "reference": "SYNTH_R",
                "tool": "SYNTH A",
                "run_label": "SYNTH A run 1",
                "annotation": "inputs/query.gtf",
            },
            {
                "run_id": "SYNTH_a2",
                "reference": "SYNTH_R",
                "tool": "SYNTH A",
                "run_label": "SYNTH A run 2",
                "annotation": "inputs/query_a2.gtf",
            },
        ],
        "display": {"preferred_tools": ["SYNTH B"]},
        **extra,
    }
    path = root / "config.json"
    path.write_text(json.dumps(config))
    return path


# --------------------------------------------------------------- configuration
@pytest.mark.parametrize(
    "change, message",
    [
        (
            lambda c: c["references"][0]["analysis_scopes"].append(
                {"id": "all", "sequences": ["chr1"]}
            ),
            "other than",
        ),
        (
            lambda c: c["references"][0]["analysis_scopes"].append(
                {"id": "SYNTH_x", "sequences": []}
            ),
            "non-empty",
        ),
        (
            lambda c: c["runs"].append({**c["runs"][1], "run_id": "SYNTH_dup"}),
            "could not be told apart",
        ),
        (
            lambda c: c["runs"].append(
                {
                    **c["runs"][1],
                    "run_id": "SYNTH_c",
                    "companion_of": "SYNTH_nope",
                    "role": "supplementary",
                    "comparator_options": {"add_stop_codon": True},
                }
            ),
            "companion_of",
        ),
        (
            lambda c: c["runs"].append(
                {
                    **c["runs"][1],
                    "run_id": "SYNTH_c",
                    "companion_of": "SYNTH_b",
                    "role": "supplementary",
                    "comparator_options": {"add_stop_codon": True},
                }
            ),
            "same annotation file",
        ),
        (
            lambda c: c["display"].update(preferred_tools=["SYNTH nope"]),
            "preferred_tools",
        ),
        (lambda c: c["runs"][0].update(source={"sha256": "xyz"}), "64 hex digits"),
        (lambda c: c.update(config_version=3), "config_version"),
    ],
)
def test_configuration_is_validated(tmp_path, change, message):
    path = _write(tmp_path)
    config = json.loads(path.read_text())
    change(config)
    path.write_text(json.dumps(config))
    with pytest.raises(ConfigError, match=message):
        load_experiment(path)


def test_version_one_configuration_still_loads(tmp_path):
    path = _write(tmp_path)
    config = json.loads(path.read_text())
    config["config_version"] = 1
    config["display"] = {
        "pinned_runs": ["SYNTH_a"],
        "order_metric": "exact_recall_one_to_one",
    }
    path.write_text(json.dumps(config))
    exp = load_experiment(path)
    assert (
        exp.display["order_metric"] == "exact_cds_one_to_one_recall"
    )  # earlier id, same definition
    assert exp.runs["SYNTH_a"].display_label == "SYNTH A run 1"


def test_analysis_scopes_are_passed_to_the_comparator(tmp_path):
    exp = load_experiment(_write(tmp_path))
    ws = Workspace(exp)
    task = next(
        t
        for t in runs_mod.plan(exp, ws)
        if t.run.run_id == "SYNTH_a" and t.policy.id == "B"
    )
    argv = task.argv_template
    assert argv[argv.index("--analysis-scope") + 1] == "SYNTH_chr1=chr1"
    assert argv[argv.index("--export-evaluated-models") + 1] == "ids"
    assert task.key_components["parameters"]["analysis_scopes"] == {
        "SYNTH_chr1": ["chr1"]
    }


# --------------------------------------------------------------- comparator results end to end
@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("scope")
    experiment = load_experiment(_write(root))
    ws = Workspace(experiment)
    runs_mod.run_missing(experiment, ws, log=QUIET)
    build_dashboard(experiment, ws, log=QUIET)
    return {"experiment": experiment, "ws": ws}


def test_scope_rows_come_from_the_comparator(built):
    con = queries.connect(built["ws"].dashboard_db())
    meta = queries.meta(con)
    scopes = {s["id"]: s for s in meta["references"][0]["scopes"]}
    assert (
        scopes["all"]["label"] == "Whole synthetic assembly"
        and scopes["SYNTH_chr1"]["available"]
    )
    runs = {r["run_id"]: r for r in meta["runs"]}
    assert (runs["SYNTH_a"]["label"], runs["SYNTH_a2"]["label"]) == (
        "SYNTH A run 1",
        "SYNTH A run 2",
    )
    assert meta["display"]["preferred_tools"] == ["SYNTH B"]
    rows = {
        (m["scope_id"], m["policy_id"], m["run_id"], m["metric_id"]): m
        for m in queries.overview(con, "SYNTH_R")["metrics"]
    }
    f1 = rows[("SYNTH_chr1", "B", "SYNTH_a", "exact_cds_one_to_one_f1")]
    assert (f1["numerator"], f1["denominator"], f1["value"]) == (2, 2, 1.0)
    whole = rows[("all", "B", "SYNTH_a", "exact_cds_one_to_one_f1")]
    assert (whole["numerator"], whole["denominator"]) == (
        2,
        5,
    )  # R = 3, Q = 2 on the whole assembly
    # every stored value equals the comparator's metrics.tsv row, scope by scope
    out = Path(built["ws"].run_record("SYNTH_a", "B")["output_dir"])
    for line in (out / "metrics.tsv").read_text().splitlines()[1:]:
        scope_id, metric_id, num, den, value, status, reason = line.split("\t")[:7]
        stored = rows[(scope_id, "B", "SYNTH_a", metric_id)]
        assert (stored["value_text"], stored["status"], stored["reason"] or "") == (
            value,
            status,
            reason,
        )
        assert stored["numerator"] == (None if num == "NA" else int(num))
        assert stored["denominator"] == (None if den in ("NA", "") else int(den))


def test_splice_tables_are_the_comparator_tables(built):
    con = queries.connect(built["ws"].dashboard_db())
    b = queries.splice_sites(
        con, {"reference": "SYNTH_R", "policy": "B", "scope": "all"}
    )
    rows = [
        r
        for r in b["lengths"]
        if r["run_id"] == "SYNTH_a" and r["eligible"] + r["missing_locus"]
    ]
    assert [
        (r["eligible"], r["exact"], r["donor"], r["acceptor"], r["missing_locus"])
        for r in rows
    ] == [(4, 3, 3, 4, 1)]
    a = queries.splice_sites(
        con, {"reference": "SYNTH_R", "policy": "A", "scope": "all"}
    )
    assert a["sites"] == [] and a["lengths"] == []  # not computed under policy A
    chr1 = queries.splice_sites(
        con, {"reference": "SYNTH_R", "policy": "B", "scope": "SYNTH_chr1"}
    )
    donors = [
        r
        for r in chr1["sites"]
        if r["run_id"] == "SYNTH_a"
        and r["kind"] == "donor"
        and r["coordinate"] == "coding"
        and r["min_intron_length"] == 0
    ]
    assert sum(r["missing_locus"] + r["assessed"] for r in donors) == 2
    results = {
        r["run_id"]: r
        for r in queries.overview(con, "SYNTH_R")["results"]
        if r["policy_id"] == "A"
    }
    assert (
        results["SYNTH_a"]["analyses"]["splice_position_length"]["status"]
        == "not_applicable"
    )


def test_gene_table_filters_by_scope_and_one_to_one(built):
    con = queries.connect(built["ws"].dashboard_db())
    q = lambda **p: queries.genes(con, {"reference": "SYNTH_R", "policy": "B", **p})
    assert q()["total"] == 3 and q(scope="SYNTH_chr1")["total"] == 1
    assert {r["gene_id"] for r in q(focus_run="SYNTH_a", one_to_one="1")["rows"]} == {
        "SYNTH_G1"
    }
    assert {r["gene_id"] for r in q(focus_run="SYNTH_b", one_to_one="0")["rows"]} == {
        "SYNTH_G2",
        "SYNTH_G3",
    }
    assert {r["gene_id"] for r in q(multi="1")["rows"]} == {
        "SYNTH_G1",
        "SYNTH_G2",
        "SYNTH_G3",
    }
    preds = queries.query_genes(
        con,
        {
            "reference": "SYNTH_R",
            "policy": "B",
            "run": "SYNTH_a",
            "scope": "SYNTH_chr1",
        },
    )
    assert [r["gene_id"] for r in preds["rows"]] == ["SYNTH_q1"]

"""Sequence scopes, run display order and splice-site analysis (SYNTHETIC fixtures)."""

import json
import shutil
import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.benchmark import runs as runs_mod
from ensembl.genes.annotation_qc.benchmark import splice
from ensembl.genes.annotation_qc.benchmark.config import ConfigError, load_experiment
from ensembl.genes.annotation_qc.benchmark.scope import (
    check_whole_assembly,
    mark_inconsistent,
    scoped_metrics,
)
from ensembl.genes.annotation_qc.benchmark.store import build_dashboard
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace
from ensembl.genes.annotation_qc.dashboard import queries

QUIET = lambda *a, **k: None  # noqa: E731
THREE = ((101, 150), (201, 260), (301, 350))  # CDS 250 bp span, 160 coding bases


# --------------------------------------------------------------- splice-site geometry
def test_plus_strand_sites_and_positions():
    first, second = splice.intron_sites("+", THREE)
    assert (first.lo, first.hi, first.donor, first.acceptor, first.length) == (151, 200, 151, 200, 50)
    assert first.position["cds_span"] == {"donor": 50 / 250, "acceptor": 100 / 250}
    assert first.position["coding"] == {"donor": 50 / 160, "acceptor": 50 / 160}
    assert second.position["cds_span"] == {"donor": 160 / 250, "acceptor": 200 / 250}


def test_minus_strand_is_oriented_in_transcription_direction():
    # transcription runs from 350 down to 101: the first intron met is 261-300
    by_lo = {i.lo: i for i in splice.intron_sites("-", THREE)}
    upstream, downstream = by_lo[261], by_lo[151]
    assert (upstream.donor, upstream.acceptor) == (300, 261)
    assert (downstream.donor, downstream.acceptor) == (200, 151)
    assert upstream.position["cds_span"] == {"donor": 50 / 250, "acceptor": 90 / 250}
    assert downstream.position["cds_span"] == {"donor": 150 / 250, "acceptor": 200 / 250}
    assert upstream.position["coding"]["donor"] == 50 / 160  # coding bases 301-350 upstream
    # a mirrored gene gets mirrored positions
    mirrored = tuple(sorted((1001 - e, 1001 - s) for s, e in THREE))
    plus = {round(i.position["cds_span"]["donor"], 6) for i in splice.intron_sites("+", mirrored)}
    minus = {round(i.position["cds_span"]["donor"], 6) for i in splice.intron_sites("-", THREE)}
    assert plus == minus
    assert splice.intron_sites("+", ((1, 90),)) == []  # single segment: no CDS intron


def test_intron_length_bins_cover_every_length():
    edges = splice.INTRON_LENGTH_EDGES
    assert [splice.length_bin(n) for n in (0, 9, 10, 29, 30, 59, 60, 99_999, 100_000, 10**7)] == [
        0, 0, 1, 1, 2, 2, 3, len(edges) - 2, len(edges) - 1, len(edges) - 1
    ]
    assert splice.length_bin_label(0) == "< 10 bp"
    assert splice.length_bin_label(3) == "60 bp–100 bp"
    assert splice.length_bin_label(len(edges) - 1) == "≥ 100 kb"
    assert splice.position_bin(0.999) == 9 and splice.position_bin(0.0) == 0


def _gene(gid, chrom, strand, *transcripts):
    return splice.ReferenceGene(gid, chrom, strand, list(transcripts))


def _cell(result, coordinate, kind, b, min_length=0):
    return next(
        r for r in result["position"]
        if (r["min_intron_length"], r["coordinate"], r["kind"], r["bin"]) == (min_length, coordinate, kind, b)
    )


def test_support_missing_loci_and_empty_bins():
    # query: chr1 exact; chr2 (minus) first intron with a shifted donor (151-204)
    query = splice.SupportIndex.build(
        [("chr1", "+", THREE), ("chr2", "-", ((101, 150), (205, 260), (301, 350))), ("chr3", "+", THREE)]
    )
    genes = [
        _gene("G1", "chr1", "+", THREE),
        _gene("G2", "chr2", "-", THREE),
        _gene("G3", "chr2", "-", ((601, 650), (701, 760))),  # locus not recovered
    ]
    detected = {"G1": True, "G2": True, "G3": False}
    out = splice.analyse(genes, query, detected, {"all": None, "chr1": {"chr1"}}, one_per_gene=True)
    whole = out["all"]
    length = {r["bin"]: r for r in whole["length"]}[splice.length_bin(50)]
    assert (length["eligible"], length["exact"], length["donor"], length["acceptor"], length["missing_locus"]) == (4, 3, 3, 4, 1)
    # donor bin 6 holds G1's second donor (supported) and G2's 151-200 donor at 200 (not supported)
    cell = _cell(whole, "cds_span", "donor", 6)
    assert (cell["assessed"], cell["supported"], cell["missing_locus"]) == (2, 1, 0)
    # empty bins are explicit zeros (shown as not available), never dropped
    empty = _cell(whole, "cds_span", "donor", 0)
    assert (empty["assessed"], empty["supported"], empty["missing_locus"]) == (0, 0, 0)
    assert len(whole["position"]) == len(splice.POSITION_MIN_LENGTHS) * 2 * 2 * splice.POSITION_BINS
    donors = [r for r in whole["position"] if r["kind"] == "donor" and r["coordinate"] == "coding"]
    assert sum(r["missing_locus"] for r in donors if r["min_intron_length"] == 0) == 1
    # every intron here is 50 bp: the "introns >= 60 bp" rows count nothing
    assert not any(r["assessed"] + r["missing_locus"] for r in donors if r["min_intron_length"] == 60)
    # unique-intron totals; the scope applies to query introns too (chr3 is outside both)
    assert (whole["reference_introns"], whole["query_introns"], whole["shared_introns"]) == (5, 6, 3)
    chr1 = out["chr1"]
    assert (chr1["reference_introns"], chr1["query_introns"], chr1["shared_introns"], chr1["genes"]) == (2, 2, 2, 1)
    assert sum(r["missing_locus"] for r in chr1["length"]) == 0


def test_alternative_isoforms_are_deduplicated_or_refused():
    alt = ((101, 150), (201, 260))  # shares intron 151-200 with THREE
    genes = [_gene("G1", "chr1", "+", THREE, alt)]
    out = splice.analyse(genes, splice.SupportIndex(), {"G1": True}, {"all": None}, one_per_gene=False)
    assert out["all"]["reference_introns"] == 2 and "position" not in out["all"]
    with pytest.raises(ValueError, match="one transcript per gene"):
        splice.analyse(genes, splice.SupportIndex(), {"G1": True}, {"all": None}, one_per_gene=True)
    with pytest.raises(ValueError, match="no locus outcome"):
        splice.analyse(genes, splice.SupportIndex(), {}, {"all": None}, one_per_gene=False)


# --------------------------------------------------------------- scoped metrics
def _outcomes():
    ref = pd.DataFrame(
        {
            "gene_id": ["r1", "r2", "r3", "r4"],
            "chrom": ["1", "1", "2", "Y"],
            "cds_overlap": [1.0, 0.5, 1.0, 0.0],
            "cds_exact": [1, 0, 1, 0],
            "cds_intron_chain_match": ["True", "False", "True", "NA"],
            "counterpart_count": [1, 2, 1, 0],
            "multi_segment": [1, 1, 1, 0],
        }
    )
    qry = pd.DataFrame(
        {"gene_id": ["q1", "q2", "q3", "q4", "q5"], "chrom": ["1", "1", "1", "2", "Y"],
         "cds_exact": [1, 0, 0, 1, 0], "counterpart_count": [1, 1, 1, 1, 0]}
    )
    pairs = pd.DataFrame({"chrom": ["1", "2"]})
    return ref, qry, pairs


def test_scope_is_applied_to_reference_query_and_pairs():
    ref, qry, pairs = _outcomes()
    m = {r["metric_id"]: r for r in scoped_metrics(ref, qry, pairs, {"1"})}
    assert (m["reference_coding_genes"]["numerator"], m["predicted_coding_genes"]["numerator"]) == (2, 3)
    assert (m["exact_tp_one_to_one"]["numerator"], m["exact_recall_one_to_one"]["value"]) == (1, 0.5)
    f1 = m["exact_f1_one_to_one"]
    assert (f1["numerator"], f1["denominator"], f1["value"], f1["source"]) == (2, 5, 0.4, "scope_filter")
    assert (m["locus_recovery_cds_overlap"]["numerator"], m["splits"]["numerator"]) == (2, 1)
    assert (m["cds_intron_chain_best_pair"]["numerator"], m["cds_intron_chain_best_pair"]["denominator"]) == (1, 2)
    assert m["ref_cds_introns_recovered"]["status"] == "not_available"
    whole = {r["metric_id"]: r for r in scoped_metrics(ref, qry, pairs, None)}
    assert whole["exact_f1_one_to_one"]["value"] == round(4 / 9, 4)
    empty = {r["metric_id"]: r for r in scoped_metrics(ref, qry, pairs, {"MT"})}
    assert empty["exact_f1_one_to_one"]["status"] == "undefined"
    assert empty["locus_recovery_cds_overlap"]["value"] is None
    no_pairs = {r["metric_id"]: r for r in scoped_metrics(ref, qry, None, {"1"})}
    assert no_pairs["exact_f1_one_to_one"]["status"] == "not_available"


def test_whole_assembly_check_marks_scoped_rows_inconsistent():
    ref, qry, pairs = _outcomes()
    derived = scoped_metrics(ref, qry, pairs, None)
    comparator = [dict(r, source="comparator") for r in derived]
    assert check_whole_assembly(derived, comparator)["status"] == "verified"
    comparator[2] = {**comparator[2], "numerator": 99}
    check = check_whole_assembly(derived, comparator)
    assert check["status"] == "disagrees" and "locus_recovery_cds_overlap.numerator" in check["differences"][0]
    marked = mark_inconsistent(scoped_metrics(ref, qry, pairs, {"1"}), check)
    assert {r["status"] for r in marked} == {"inconsistent", "not_available"}


# --------------------------------------------------------------- configuration
def _write(root: Path, **extra) -> Path:
    inputs = root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    (inputs / "reference.gff3").write_text(REFERENCE)
    (inputs / "query.gtf").write_text(QUERY)
    config = {
        "config_version": 1,
        "experiment": {"id": "SYNTH_scope", "synthetic": True},
        "output_dir": "workspace",
        "evaluation": {"evaluation_mode": "protein_coding", "reference_transcript_biotypes": ["protein_coding"]},
        "policies": {
            "A": {"reference_transcript_selection": "all", "query_transcript_selection": "all"},
            "B": {"reference_transcript_selection": "longest_cds", "query_transcript_selection": "longest_cds"},
        },
        "references": [
            {
                "id": "SYNTH_R",
                "species": "synthetic",
                "assembly": "SYNTH_v1",
                "annotation": "inputs/reference.gff3",
                "scope": {"type": "all", "label": "Whole synthetic assembly"},
                "analysis_scopes": [
                    {"id": "SYNTH_chr1", "label": "chr1 only", "sequences": ["chr1"]},
                    {"id": "SYNTH_absent", "label": "not compared", "sequences": ["chr1", "chrZ"]},
                ],
            }
        ],
        "runs": [
            {"run_id": "SYNTH_b", "reference": "SYNTH_R", "tool": "SYNTH B", "annotation": "inputs/query.gtf"},
            {"run_id": "SYNTH_a", "reference": "SYNTH_R", "tool": "SYNTH A", "annotation": "inputs/query.gtf"},
            {"run_id": "SYNTH_a_suppl", "reference": "SYNTH_R", "tool": "SYNTH A", "label": "SYNTH A supplementary",
             "role": "supplementary", "companion_of": "SYNTH_a", "annotation": "inputs/query.gtf"},
        ],
        "display": {"pinned_runs": ["SYNTH_a"]},
        **extra,
    }
    path = root / "config.json"
    path.write_text(json.dumps(config))
    return path


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda c: c["references"][0]["analysis_scopes"].append({"id": "whole_assembly", "sequences": ["chr1"]}), "other than"),
        (lambda c: c["references"][0]["analysis_scopes"].append({"id": "SYNTH_x", "sequences": []}), "non-empty"),
        (lambda c: c["runs"][2].update(companion_of="SYNTH_nope"), "companion_of"),
        (lambda c: c["runs"][2].update(role="primary"), "needs role 'supplementary'"),
        (lambda c: c["display"].update(pinned_runs=["SYNTH_nope"]), "pinned_runs"),
    ],
)
def test_presentation_config_is_validated(tmp_path, change, message):
    path = _write(tmp_path)
    config = json.loads(path.read_text())
    change(config)
    path.write_text(json.dumps(config))
    with pytest.raises(ConfigError, match=message):
        load_experiment(path)


# --------------------------------------------------------------- comparator results end to end
# chr1: SYNTH_G1 (+) three CDS segments. chr2: SYNTH_G2 (-) three segments and SYNTH_G3 (-)
# with no prediction. The query reproduces SYNTH_G1 exactly and SYNTH_G2 with one donor shifted.
def _ref_gene(chrom, strand, gid, segments):
    lo, hi = segments[0][0], segments[-1][1]
    lines = [
        f"{chrom}\tSYNTH\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID=gene:{gid};Name={gid};biotype=protein_coding;gene_id={gid}",
        f"{chrom}\tSYNTH\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID=transcript:{gid}.t;Parent=gene:{gid};biotype=protein_coding;transcript_id={gid}.t",
    ]
    for s, e in segments:
        lines.append(f"{chrom}\tSYNTH\texon\t{s}\t{e}\t.\t{strand}\t.\tParent=transcript:{gid}.t")
        lines.append(f"{chrom}\tSYNTH\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=CDS:{gid}.t;Parent=transcript:{gid}.t")
    return "\n".join(lines)


def _query_gene(chrom, strand, gid, segments):
    return "\n".join(
        f'{chrom}\tSYNTH_pred\t{feature}\t{s}\t{e}\t.\t{strand}\t{"0" if feature == "CDS" else "."}\tgene_id "{gid}"; transcript_id "{gid}.t1";'
        for s, e in segments
        for feature in ("exon", "CDS")
    )


REFERENCE = "\n".join(
    [
        "##gff-version 3",
        "##sequence-region chr1 1 1000",
        "##sequence-region chr2 1 1000",
        _ref_gene("chr1", "+", "SYNTH_G1", THREE),
        _ref_gene("chr2", "-", "SYNTH_G2", THREE),
        _ref_gene("chr2", "-", "SYNTH_G3", ((601, 650), (701, 760))),
    ]
) + "\n"
QUERY = "\n".join(
    [
        _query_gene("chr1", "+", "SYNTH_q1", THREE),
        _query_gene("chr2", "-", "SYNTH_q2", ((101, 150), (205, 260), (301, 350))),
    ]
) + "\n"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("scope")
    experiment = load_experiment(_write(root))
    ws = Workspace(experiment)
    runs_mod.run_missing(experiment, ws, log=QUIET)
    build_dashboard(experiment, ws, log=QUIET)
    return {"experiment": experiment, "ws": ws}


def test_build_scopes_order_and_analysis(built):
    con = queries.connect(built["ws"].dashboard_db())
    meta = queries.meta(con)
    scopes = {s["id"]: s for s in meta["references"][0]["scopes"]}
    assert scopes["whole_assembly"]["label"] == "Whole synthetic assembly"
    assert scopes["SYNTH_chr1"]["available"] and not scopes["SYNTH_absent"]["available"]
    assert "chrZ" in scopes["SYNTH_absent"]["reason"]
    runs = {r["run_id"]: r for r in meta["runs"]}
    assert runs["SYNTH_a_suppl"]["companion_of"] == "SYNTH_a"
    assert meta["display"]["pinned_runs"] == ["SYNTH_a"]
    ov = queries.overview(con, "SYNTH_R")
    assert {a["status"] for a in ov["analyses"]} == {"ok"}
    record = ov["analyses"][0]["record"]
    assert record["whole_assembly_check"]["status"] == "verified"
    assert {"support_rule", "intron_length_edges", "inputs", "code"} <= set(record)
    assert record["inputs"]["comparison_details.tsv"]
    sm = {(m["scope_id"], m["policy_id"], m["run_id"], m["metric_id"]): m for m in ov["scope_metrics"]}
    f1 = sm[("SYNTH_chr1", "B", "SYNTH_a", "exact_f1_one_to_one")]
    assert (f1["numerator"], f1["denominator"], f1["value"]) == (2, 2, 1.0)
    whole_f1 = next(m for m in ov["metrics"] if (m["run_id"], m["policy_id"], m["metric_id"]) == ("SYNTH_a", "B", "exact_f1_one_to_one"))
    assert (whole_f1["numerator"], whole_f1["denominator"]) == (2, 5)  # R = 3, Q = 2 on the whole assembly
    assert not any(k[0] == "SYNTH_absent" for k in sm)
    # binned splice analysis: policy B only; minus-strand donor shift detected; missing locus separate
    with pytest.raises(KeyError):
        queries.splice_sites(con, {"reference": "SYNTH_R", "scope": "whole_assembly"})
    b = queries.splice_sites(con, {"reference": "SYNTH_R", "policy": "B", "scope": "whole_assembly"})
    rows = [r for r in b["lengths"] if r["run_id"] == "SYNTH_a" and r["eligible"] + r["missing_locus"]]
    assert [(r["eligible"], r["exact"], r["donor"], r["acceptor"], r["missing_locus"]) for r in rows] == [(4, 3, 3, 4, 1)]
    a = queries.splice_sites(con, {"reference": "SYNTH_R", "policy": "A", "scope": "whole_assembly"})
    assert a["sites"] == [] and a["lengths"] == []
    chr1 = queries.splice_sites(con, {"reference": "SYNTH_R", "policy": "B", "scope": "SYNTH_chr1"})
    assert sum(r["missing_locus"] + r["assessed"] for r in chr1["sites"] if r["run_id"] == "SYNTH_a" and r["kind"] == "donor" and r["coordinate"] == "coding" and r["min_intron_length"] == 0) == 2


def test_dataset_without_analysis_tables_still_loads(built, tmp_path):
    db = tmp_path / "v2.sqlite"
    shutil.copy(built["ws"].dashboard_db(), db)
    with sqlite3.connect(db) as raw:
        raw.execute("UPDATE meta SET value='\"2\"' WHERE key='store_version'")
        for table in ("scope_metrics", "splice_sites", "intron_lengths", "analyses"):
            raw.execute(f"DROP TABLE {table}")
        meta = json.loads(raw.execute("SELECT value FROM meta WHERE key='references'").fetchone()[0])
        for ref in meta:
            ref.pop("scopes")
        raw.execute("UPDATE meta SET value=? WHERE key='references'", (json.dumps(meta),))
    con = queries.connect(db)
    ref = queries.meta(con)["references"][0]
    assert [s["id"] for s in ref["scopes"]] == ["whole_assembly"] and not ref["splice_analysis_available"]
    ov = queries.overview(con, "SYNTH_R")
    assert ov["scope_metrics"] == [] and ov["analyses"] == [] and ov["metrics"]
    with pytest.raises(queries.QueryError, match="store version 3"):
        queries.splice_sites(con, {"reference": "SYNTH_R", "policy": "B"})

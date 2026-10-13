"""Output contract 1.0 of annotation-qc pairwise-compare (metrics.tsv and new files)."""

# pylint: disable=redefined-outer-name

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.metrics.pairwise import contract
from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    build_gene_models,
    cds_differences,
)
from ensembl.genes.annotation_qc.parsers.annotation import (
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.runners import pairwise_compare

from .conftest import gene_records

THREE = [(100, 200), (300, 400), (500, 600)]


def _run(argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    pairwise_compare.register(subparsers)
    args = parser.parse_args(["pairwise-compare", *argv])
    return args.func(args)


def _metrics(outdir) -> pd.DataFrame:
    return pd.read_csv(
        outdir / "metrics.tsv", sep="\t", dtype={"reason": str}, keep_default_na=False
    ).set_index(["scope_id", "metric_id"])


def _num(table, metric, scope="all"):
    return int(table.loc[(scope, metric), "numerator"])


def _transcripts(outdir) -> pd.DataFrame:
    return pd.read_csv(
        outdir / "evaluated_transcripts.tsv.gz",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


@pytest.fixture
def two_chromosomes(write_gff3):
    """Reference and query on chromosomes 1 and 2, both strands, a split and a lncRNA."""
    ref = (
        gene_records("1", "R1", "+", {"R1.1": THREE, "R1.2": THREE[:2]})
        + gene_records("1", "R2", "-", {"R2.1": [(2000, 2100), (2200, 2300)]})
        + gene_records("1", "R3", "+", {"R3.1": [(5000, 5600)]})
        + gene_records(
            "1", "NC", "+", {"NC.1": [(9000, 9100)]}, biotype="lncRNA", with_cds=False
        )
        + gene_records("2", "R4", "+", {"R4.1": [(100, 200), (300, 400)]})
        + gene_records("2", "R5", "-", {"R5.1": [(1000, 1100), (1200, 1300)]})
    )
    query = (
        gene_records("1", "Q1", "+", {"Q1.1": THREE})
        # minus strand: same introns, CDS longer at the start end (higher coordinate)
        + gene_records(
            "1",
            "Q2",
            "-",
            {"Q2.1": [(2000, 2100), (2200, 2309)]},
        )
        + gene_records("1", "Q3a", "+", {"Q3a.1": [(5000, 5200)]})
        + gene_records("1", "Q3b", "+", {"Q3b.1": [(5400, 5600)]})
        + gene_records("1", "QN", "+", {"QN.1": [(20000, 20300)]})
        + gene_records("2", "Q4", "+", {"Q4.1": [(100, 200), (300, 400)]})
        + gene_records("2", "Q5", "-", {"Q5.1": [(1000, 1100), (1210, 1300)]})
    )
    return write_gff3("ref.gff3", ref), write_gff3("query.gff3", query)


@pytest.fixture
def compared(two_chromosomes, tmp_path):
    ref, query = two_chromosomes
    outdir = tmp_path / "out"
    summary = _run(
        [
            "--reference",
            ref,
            "--query",
            query,
            "--evaluation-mode",
            "protein_coding",
            "--analysis-scope",
            "chr1=1",
            "--analysis-scope",
            "chr2=2",
            "--export-evaluated-models",
            "full",
            "--outdir",
            str(outdir),
        ]
    )
    return summary, outdir


class TestMetricsTable:
    def test_schema_and_status_semantics(self, compared):
        _, outdir = compared
        table = _metrics(outdir).reset_index()
        assert list(table.columns) == contract.METRIC_COLUMNS
        assert set(table["scope_id"]) == {"all", "chr1", "chr2"}
        for _, rows in table.groupby("scope_id"):
            assert list(rows["metric_id"]) == list(contract.METRICS)
        ok = table["status"] == "ok"
        assert (table.loc[ok, "value"] != "NA").all()
        assert (table.loc[~ok, "value"] == "NA").all()
        assert (table.loc[~ok, "reason"] != "").all()
        assert set(table.loc[~ok, "reason"]) <= set(contract.REASONS)

    def test_values_match_the_summary(self, compared):
        summary, outdir = compared
        table = _metrics(outdir)
        one = summary["cds_exact_one_to_one"]
        assert _num(table, "reference_genes") == summary["total_reference_genes"] == 5
        assert _num(table, "query_genes") == summary["total_consensus_genes"] == 7
        assert (
            _num(table, "locus_recovery_cds_overlap")
            == summary["cds_overlap_locus_recovery"]["recovered_count"]
        )
        assert _num(table, "exact_cds_one_to_one_tp") == one["true_positives"]
        assert (
            _num(table, "cds_intron_chain_best_pair")
            == summary["sensitivity_cds"]["cds_intron_chain_recovered"]
        )
        assert (
            int(table.loc[("all", "cds_intron_chain_best_pair"), "denominator"])
            == summary["sensitivity_cds"]["multi_segment_cds_reference_genes"]
        )
        assert (
            _num(table, "split_reference_genes")
            == summary["split_merge"]["gene_split_count"]
        )

    def test_invariants(self, compared):
        _, outdir = compared
        table = _metrics(outdir)
        for scope in ("all", "chr1", "chr2"):
            r = _num(table, "reference_genes", scope)
            q = _num(table, "query_genes", scope)
            tp = _num(table, "exact_cds_one_to_one_tp", scope)
            f1 = table.loc[(scope, "exact_cds_one_to_one_f1")]
            assert (int(f1["numerator"]), int(f1["denominator"])) == (2 * tp, r + q)
            assert float(f1["value"]) == pytest.approx(2 * tp / (r + q))
            assert (
                _num(table, "cds_coordinate_exact_reference", scope)
                - _num(table, "exact_cds_reference_with_candidate_unmatched", scope)
                == tp
            )
            assert (
                _num(table, "query_matched", scope)
                + _num(table, "query_novel", scope)
                + _num(table, "strand_mismatch_query", scope)
                == q
            )
            assert (
                sum(
                    _num(table, f"cds_outcome_{k}", scope)
                    for k in contract.CDS_OUTCOMES
                )
                == r
            )
            assert _num(table, "cds_intron_chain_any_pair", scope) >= _num(
                table, "cds_intron_chain_best_pair", scope
            )

    def test_scopes_partition_the_counts(self, compared):
        _, outdir = compared
        table = _metrics(outdir).reset_index()
        counts = table[table["denominator"].isna() | (table["status"] == "ok")]
        for metric in (
            "reference_genes",
            "query_genes",
            "locus_recovery_cds_overlap",
            "exact_cds_one_to_one_tp",
            "ref_unique_cds_introns_recovered",
            "split_reference_genes",
        ):
            values = counts[counts["metric_id"] == metric].set_index("scope_id")
            values = values["numerator"].astype(int)
            assert values["all"] == values["chr1"] + values["chr2"]

    def test_scope_equals_a_region_restricted_run(self, two_chromosomes, tmp_path):
        ref, query = two_chromosomes
        common = ["--reference", ref, "--query", query]
        _run(common + ["--analysis-scope", "chr2=2", "--outdir", str(tmp_path / "a")])
        _run(common + ["--region", "2", "--outdir", str(tmp_path / "b")])
        scoped = _metrics(tmp_path / "a").xs("chr2")
        region = _metrics(tmp_path / "b").xs("all")
        for metric in contract.METRICS:
            if metric.endswith("_evaluated"):  # transcripts outside regions excluded
                continue
            assert (
                scoped.loc[metric, ["numerator", "denominator", "status"]].tolist()
                == region.loc[metric, ["numerator", "denominator", "status"]].tolist()
            ), metric

    def test_undefined_is_not_zero(self, write_gff3, tmp_path):
        ref = write_gff3("r.gff3", gene_records("1", "R", "+", {"R.1": [(10, 90)]}))
        query = write_gff3("q.gff3", gene_records("1", "Q", "+", {"Q.1": [(10, 90)]}))
        summary = _run(
            ["--reference", ref, "--query", query, "--outdir", str(tmp_path / "o")]
        )
        table = _metrics(tmp_path / "o")
        row = table.loc[("all", "cds_intron_chain_best_pair")]
        assert (row["status"], row["reason"], row["value"]) == (
            "undefined",
            "zero_denominator",
            "NA",
        )
        assert summary["intron_chain"]["exon_intron_chain_rate"] is None
        assert (
            summary["intron_support"]["cds"]["reference_introns_recovered_rate"] is None
        )
        tsv = (tmp_path / "o" / "comparison_summary.tsv").read_text()
        assert "intron_chain_exon_intron_chain_rate\tNA" in tsv
        assert _num(table, "cds_coordinate_exact_reference") == 1  # a real value


class TestPerGeneColumns:
    def test_minus_strand_descriptors_and_any_pair(self, compared):
        _, outdir = compared
        details = pd.read_csv(
            outdir / "comparison_details.tsv",
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        ref = details[details["source"] == "reference"].set_index("gene_id")
        # R2 (-): the query CDS ends 9 bp further at the high coordinate = start side.
        assert ref.loc["R2", "cds_start_diff_bp"] == "9"
        assert ref.loc["R2", "cds_stop_diff_bp"] == "0"
        assert ref.loc["R2", "cds_intron_chain_any_pair"] == "True"
        assert ref.loc["R2", "cds_outcome"] == "terminal_diff"
        # R5 (-): acceptor moved; the chain differs.
        assert ref.loc["R5", "cds_introns_missing"] == "1"
        assert ref.loc["R5", "cds_introns_extra"] == "1"
        assert ref.loc["R5", "cds_intron_chain_any_pair"] == "False"
        assert ref.loc["R1", "one_to_one_partner"] == "Q1"
        assert ref.loc["R1", "cds_outcome"] == "exact"
        assert ref.loc["R3", "counterpart_count"] == "2"
        assert ref.loc["R3", "cds_intron_chain_any_pair"] == "NA"
        assert ref.loc["R4", "scope_ids"] == "chr2"
        assert ref.loc["R1", "cds_shared_bases"] == "303"

    def test_cds_differences_plus_and_minus(self):
        ref = ((100, 200), (300, 400))
        query = ((90, 200), (300, 403))
        assert cds_differences(ref, query, "+") == (10, 3, 0, 0)
        assert cds_differences(ref, query, "-") == (3, 10, 0, 0)


class TestEvaluatedTranscripts:
    def test_status_and_reasons(self, write_gff3, write_text, tmp_path):
        ref = write_gff3(
            "r.gff3",
            gene_records("1", "R", "+", {"R.1": THREE, "R.2": THREE[:1]})
            + gene_records(
                "1", "NC", "+", {"NC.1": [(800, 900)]}, biotype="lncRNA", with_cds=False
            )
            + gene_records("1", "FAR", "+", {"FAR.1": [(50000, 50100)]})
            + gene_records("MT", "MTG", "+", {"MTG.1": [(10, 90)]}),
        )
        query = write_gff3(
            "q.gff3",
            gene_records("1", "g1", "+", {"g1.t1": THREE})
            + gene_records("2", "g1", "+", {"g1.t1": [(10, 90)]}),
        )
        genome = write_text("g.fa", ">1\n" + "A" * 60000 + "\n>2\n" + "A" * 200 + "\n")
        outdir = tmp_path / "o"
        _run(
            [
                "--reference",
                ref,
                "--query",
                query,
                "--genome",
                genome,
                "--genome-mismatch",
                "exclude",
                "--evaluation-mode",
                "protein_coding",
                "--reference-transcript-selection",
                "longest_cds",
                "--region",
                "1:1-10000",
                "--region",
                "2",
                "--outdir",
                str(outdir),
            ]
        )
        table = _transcripts(outdir).set_index("transcript_id")
        reasons = table["exclusion_reason"].to_dict()
        assert table.loc["R.1", "status"] == "evaluated"
        assert table.loc["R.1", "selected_by"] == "longest_cds"
        assert reasons["R.2"] == "not_selected:longest_cds"
        assert reasons["NC.1"] == "evaluation_mode:protein_coding"
        assert reasons["FAR.1"] == "outside_regions"
        assert reasons["MTG.1"] == "excluded_seqname"
        # g1 is reused on two sequences: comparison ids are namespaced, file ids kept.
        query_rows = table[table["role"] == "query"]
        assert set(query_rows.index) == {"1:g1.t1", "2:g1.t1"}
        assert set(query_rows["original_transcript_id"]) == {"g1.t1"}
        sha = hashlib.sha256(Path(ref).read_bytes()).hexdigest()
        assert set(table.loc[table["role"] == "reference", "input_sha256"]) == {sha}

    def test_models_rebuilt_without_selecting_again(self, compared, two_chromosomes):
        """IDs + the source file reproduce the evaluated models exactly."""
        _, outdir = compared
        table = _transcripts(outdir)
        evaluated = table[table["status"] == "evaluated"]
        for role, path in zip(("reference", "query"), two_chromosomes):
            rows = evaluated[evaluated["role"] == role]
            parsed = parse_annotation_for_comparison(path)
            kept = parsed[
                (parsed["Feature"] == "gene") & parsed["gene_id"].isin(rows["gene_id"])
                | parsed["transcript_id"].isin(rows["transcript_id"])
            ]
            models = build_gene_models(kept)
            for row in rows.itertuples():
                cds = ",".join(
                    f"{s + 1}-{e}" for s, e in models.cds.get(row.transcript_id, ())
                )
                exons = ",".join(
                    f"{s + 1}-{e}" for s, e in models.exons[row.transcript_id]
                )
                assert (cds, exons) == (row.cds_segments, row.exon_segments)
                assert len(models.cds.get(row.transcript_id, ())) == int(
                    row.cds_segment_count
                )


class TestCanonicalPolicy:
    def test_no_tags_is_an_error_not_a_fallback(self, two_chromosomes, tmp_path):
        ref, query = two_chromosomes
        with pytest.raises(SystemExit, match="Ensembl_canonical"):
            _run(
                [
                    "--reference",
                    ref,
                    "--query",
                    query,
                    "--reference-transcript-selection",
                    "canonical",
                    "--outdir",
                    str(tmp_path / "o"),
                ]
            )

    def test_partial_tags_are_labelled_and_counted(self, write_gff3, tmp_path):
        records = gene_records("1", "T", "+", {"T.1": THREE[:1], "T.2": THREE})
        records[1] = (*records[1][:5], records[1][5] + ";tag=Ensembl_canonical")
        records += gene_records("1", "U", "+", {"U.1": [(2000, 2100)]})
        ref = write_gff3("r.gff3", records)
        query = write_gff3("q.gff3", gene_records("1", "Q", "+", {"Q.1": THREE[:1]}))
        outdir = tmp_path / "o"
        _run(
            [
                "--reference",
                ref,
                "--query",
                query,
                "--reference-transcript-selection",
                "canonical",
                "--outdir",
                str(outdir),
            ]
        )
        table = _transcripts(outdir).set_index("transcript_id")
        assert table.loc["T.1", "selected_by"] == "canonical_tag"
        assert table.loc["U.1", "selected_by"] == "canonical_fallback_longest_cds"
        metrics = _metrics(outdir)
        row = metrics.loc[("all", "reference_canonical_fallback_genes")]
        assert (int(row["numerator"]), int(row["denominator"])) == (1, 2)


def test_cds_only_stop_excluded_query(write_text, tmp_path):
    """A gene+CDS-only GTF ending before the stop: harmonised run is labelled and linked."""
    body = "ATG" + "AAA" * 30 + "TAA"
    genome = write_text("g.fa", ">1\n" + "C" * 99 + body + "C" * 300 + "\n")
    end = 99 + 3 + 90  # CDS 100..192 excludes the stop (193..195)
    gtf = (
        # as in such files, the gene line covers the stop codon the CDS omits
        f'1\t.\tgene\t100\t{end + 3}\t.\t+\t.\tgene_id "G";\n'
        f'1\t.\tCDS\t100\t{end}\t.\t+\t0\tgene_id "G"; transcript_id "T";\n'
    )
    query = write_text("q.gtf", gtf)
    ref = write_text(
        "r.gff3",
        "##gff-version 3\n"
        f"1\tt\tgene\t100\t{end + 3}\t.\t+\t.\tID=R;biotype=protein_coding\n"
        f"1\tt\tmRNA\t100\t{end + 3}\t.\t+\t.\tID=R.1;Parent=R\n"
        f"1\tt\texon\t100\t{end + 3}\t.\t+\t.\tParent=R.1\n"
        f"1\tt\tCDS\t100\t{end + 3}\t.\t+\t0\tParent=R.1\n",
    )
    common = ["--reference", ref, "--query", query, "--genome", genome]
    _run(common + ["--outdir", str(tmp_path / "orig")])
    _run(common + ["--add-stop-codon", "--outdir", str(tmp_path / "stop")])
    original, harmonised = _metrics(tmp_path / "orig"), _metrics(tmp_path / "stop")
    assert _num(original, "exact_cds_one_to_one_tp") == 0
    assert _num(harmonised, "exact_cds_one_to_one_tp") == 1
    assert original.loc[("all", "stop_codon_extended"), "status"] == "not_applicable"
    assert _num(harmonised, "stop_codon_extended") == 1
    rows = _transcripts(tmp_path / "stop").set_index("transcript_id")
    assert rows.loc["T", "transformed"] == "stop_codon_added"
    assert rows.loc["T", "transformation_audit"] == "stop_codon_harmonisation.tsv"
    assert rows.loc["T", "record_origin"] == "inferred"
    assert rows.loc["T", "cds_length"] == str(96)
    assert (
        _transcripts(tmp_path / "orig")
        .set_index("transcript_id")
        .loc["T", "transformed"]
        == ""
    )
    manifest = json.loads((tmp_path / "stop" / "comparison_manifest.json").read_text())
    assert manifest["transformations"][0]["id"] == "stop_codon_added"
    assert manifest["inputs"]["query"]["compression"] == "none"


def test_bgzip_annotation_named_bgz(two_chromosomes, tmp_path):
    ref, query = two_chromosomes
    compressed = tmp_path / "ref.gff3.bgz"
    with open(ref, "rb") as src, gzip.open(compressed, "wb") as dst:
        shutil.copyfileobj(src, dst)
    plain = parse_annotation_for_comparison(ref)
    packed = parse_annotation_for_comparison(str(compressed))
    pd.testing.assert_frame_equal(
        plain.reset_index(drop=True), packed.reset_index(drop=True)
    )
    summary = _run(
        [
            "--reference",
            str(compressed),
            "--query",
            query,
            "--outdir",
            str(tmp_path / "o"),
        ]
    )
    assert summary["total_reference_genes"] == 6


class TestSpliceAvailability:
    def test_not_applicable_under_all(self, compared):
        _, outdir = compared
        manifest = json.loads((outdir / "comparison_manifest.json").read_text())
        assert manifest["analyses"]["splice_position_length"] == {
            "status": "not_applicable",
            "reason": "requires_one_reference_transcript_per_gene",
        }
        assert not (outdir / "splice_sites_by_position.tsv").exists()

    def test_written_per_scope_with_one_transcript_per_gene(
        self, two_chromosomes, tmp_path
    ):
        ref, query = two_chromosomes
        outdir = tmp_path / "o"
        _run(
            [
                "--reference",
                ref,
                "--query",
                query,
                "--reference-transcript-selection",
                "longest_cds",
                "--analysis-scope",
                "chr1=1",
                "--outdir",
                str(outdir),
            ]
        )
        lengths = pd.read_csv(outdir / "cds_introns_by_length.tsv", sep="\t")
        positions = pd.read_csv(outdir / "splice_sites_by_position.tsv", sep="\t")
        assert set(lengths["scope_id"]) == set(positions["scope_id"]) == {"all", "chr1"}
        whole = lengths[lengths["scope_id"] == "all"]
        # every reference CDS intron of the selected transcripts is binned once:
        # R1.1 (2), R2.1, R4.1 and R5.1 (1 each); R3.1 has a single segment
        assert int((whole["eligible"] + whole["missing_locus"]).sum()) == 5
        donors = positions[
            (positions["scope_id"] == "all")
            & (positions["kind"] == "donor")
            & (positions["coordinate"] == "cds_span")
            & (positions["min_intron_length"] == 0)
        ]
        assert int((donors["assessed"] + donors["missing_locus"]).sum()) == 5


def test_manifest_contract_and_inventory(
    compared, write_text, two_chromosomes, tmp_path
):
    _, outdir = compared
    manifest = json.loads((outdir / "comparison_manifest.json").read_text())
    assert manifest["contract_version"] == contract.CONTRACT_VERSION
    assert "comparison_manifest.json" not in manifest["outputs"]
    record = manifest["outputs"]["metrics.tsv"]
    digest = hashlib.sha256((outdir / "metrics.tsv").read_bytes()).hexdigest()
    assert record["sha256"] == digest and record["schema"] == "metrics/1.0"
    assert record["rows"] == 3 * len(contract.METRICS)
    assert manifest["implementation"]["fingerprint"]
    assert set(manifest["metric_definitions"]) == set(contract.METRICS)
    ref, query = two_chromosomes
    provenance = write_text("p.json", json.dumps({"source_path": "old_name.gtf"}))
    out = tmp_path / "p"
    _run(
        [
            "--reference",
            ref,
            "--query",
            query,
            "--provenance-json",
            provenance,
            "--outdir",
            str(out),
        ]
    )
    manifest = json.loads((out / "comparison_manifest.json").read_text())
    assert manifest["provenance"] == {"source_path": "old_name.gtf"}


def test_unknown_scope_sequence_is_rejected(two_chromosomes, write_text, tmp_path):
    ref, query = two_chromosomes
    genome = write_text("g.fa", ">1\n" + "A" * 30000 + "\n>2\n" + "A" * 3000 + "\n")
    with pytest.raises(SystemExit, match="not in the genome FASTA"):
        _run(
            [
                "--reference",
                ref,
                "--query",
                query,
                "--genome",
                genome,
                "--analysis-scope",
                "x=1,chr9",
                "--outdir",
                str(tmp_path / "o"),
            ]
        )


def test_transcript_outside_its_gene_record(write_text, tmp_path):
    """A query transcript beyond its gene line: span pairing misses, the CDS graph does not."""
    ref = write_text(
        "r.gff3",
        "##gff-version 3\n"
        "1\tt\tgene\t100\t200\t.\t+\t.\tID=R2;biotype=protein_coding\n"
        "1\tt\tmRNA\t100\t200\t.\t+\t.\tID=R2.1;Parent=R2\n"
        "1\tt\texon\t100\t200\t.\t+\t.\tParent=R2.1\n"
        "1\tt\tCDS\t100\t200\t.\t+\t0\tParent=R2.1\n"
        "1\tt\tgene\t1000\t1100\t.\t+\t.\tID=R1;biotype=protein_coding\n"
        "1\tt\tmRNA\t1000\t1100\t.\t+\t.\tID=R1.1;Parent=R1\n"
        "1\tt\texon\t1000\t1100\t.\t+\t.\tParent=R1.1\n"
        "1\tt\tCDS\t1000\t1100\t.\t+\t0\tParent=R1.1\n",
    )
    query = write_text(
        "q.gff3",
        "##gff-version 3\n"
        "1\tt\tgene\t100\t200\t.\t+\t.\tID=G\n"
        "1\tt\tmRNA\t100\t200\t.\t+\t.\tID=G.t1;Parent=G\n"
        "1\tt\texon\t100\t200\t.\t+\t.\tParent=G.t1\n"
        "1\tt\tCDS\t100\t200\t.\t+\t0\tParent=G.t1\n"
        "1\tt\tmRNA\t1000\t1100\t.\t+\t.\tID=G.t2;Parent=G\n"
        "1\tt\texon\t1000\t1100\t.\t+\t.\tParent=G.t2\n"
        "1\tt\tCDS\t1000\t1100\t.\t+\t0\tParent=G.t2\n",
    )
    summary = _run(
        ["--reference", ref, "--query", query, "--outdir", str(tmp_path / "o")]
    )
    table = _metrics(tmp_path / "o")
    # span pairing: R1 is never compared with G; the identical-CDS graph links both
    assert _num(table, "cds_coordinate_exact_reference") == 1
    assert _num(table, "exact_cds_gene_pairs") == 2
    assert _num(table, "exact_cds_one_to_one_tp") == 1
    assert _num(table, "exact_cds_reference_with_candidate_unmatched") == 1
    one = summary["cds_exact_one_to_one"]
    assert one["reference_genes_with_exact_candidate_unmatched"] == 1
    outside = summary["gene_span_checks"]["query"]["transcripts_outside"]
    assert outside["transcripts_outside_gene_span"] == 1
    assert outside["examples"] == ["G.t2"]

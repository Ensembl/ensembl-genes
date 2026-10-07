"""
Optional stop-codon harmonisation (--add-stop-codon) on synthetic annotations.

The genome is filled with C; codons are placed explicitly. Coordinates in
fixtures are one-based inclusive; the comparison schema is zero-based half-open.
"""

import hashlib
import json

import pytest

from ensembl.genes.annotation_qc.metrics.pairwise.stop_codons import (
    FRAME_CHECKS,
    GeneticCodePolicy,
    harmonise_stop_codons,
    parse_genetic_code_policy,
)
from ensembl.genes.annotation_qc.parsers.annotation import (
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.sequence import iter_sequences

from .test_pairwise_runner import _run


def _rc(text: str) -> str:
    return text.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _genome(length: int, placements: list[tuple[int, str]]) -> str:
    sequence = ["C"] * length
    for position, text in placements:
        sequence[position - 1 : position - 1 + len(text)] = list(text)
    return "".join(sequence)


def _gtf(models: list[tuple]) -> str:
    """CDS-only GTF: (seqid, gene, strand, [(start, end)], gene_start, gene_end)."""
    lines = []
    for seqid, gene, strand, cds, gene_start, gene_end in models:
        attributes = f'gene_id "{gene}"; transcript_id "{gene}.t";'
        lines.append(
            f"{seqid}\tsrc\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t{attributes}"
        )
        for start, end in cds:
            lines.append(
                f"{seqid}\tsrc\tCDS\t{start}\t{end}\t.\t{strand}\t.\t{attributes}"
            )
    return "\n".join(lines) + "\n"


# Sequence "1": one model per case. Transcript-orientation CDS, then what follows.
CHR1 = _genome(
    3000,
    [
        (101, "ATGAAA"),  # P: + strand, two segments, stop TAA after
        (201, "CCCGGG"),
        (207, "TAA"),
        (498, _rc("ATGAAACCC" + "TAG")),  # M: - strand, one segment
        (901, _rc("ATGAAA")),  # MM: - strand, two segments (901-906 first)
        (801, _rc("CCCGGG")),
        (798, _rc("TGA")),
        (1101, "ATGAAACCCTAA"),  # I: stop already in CDS
        (1201, "ATGAAACCCG"),  # D: 10 bp CDS, then TAA
        (1211, "TAA"),
        (1301, "ATGTAACCC"),  # S: internal stop, then TGA
        (1310, "TGA"),
        (1401, "AAACCCGGG"),  # N: no start codon (5'-partial), then TAA
        (1410, "TAA"),
        (1501, "ATGNAACCC"),  # A1: ambiguous base in CDS, then TAA
        (1510, "TAA"),
        (1601, "ATGAAACCC"),  # A2: ambiguous base after CDS
        (1610, "TNA"),
        (1701, "ATGAAACCC"),  # U: followed by CCC (no stop)
        (1801, "ATGAAACCC"),  # G: followed by AGA (stop only in table 2)
        (1810, "AGA"),
    ],
)
CHR1_MODELS = [
    ("1", "P", "+", [(101, 106), (201, 206)], 101, 209),
    ("1", "M", "-", [(501, 509)], 498, 509),
    ("1", "MM", "-", [(801, 806), (901, 906)], 798, 906),
    ("1", "I", "+", [(1101, 1112)], 1101, 1112),
    ("1", "D", "+", [(1201, 1210)], 1201, 1213),
    ("1", "S", "+", [(1301, 1309)], 1301, 1312),
    ("1", "N", "+", [(1401, 1409)], 1401, 1412),
    ("1", "A1", "+", [(1501, 1509)], 1501, 1512),
    ("1", "A2", "+", [(1601, 1609)], 1601, 1612),
    ("1", "U", "+", [(1701, 1709)], 1701, 1712),
    ("1", "G", "+", [(1801, 1809)], 1801, 1812),
]
EDGE = _genome(60, [(52, "ATGAAACCC"), (1, _rc("ATGAAACCC"))])
EDGE_MODELS = [
    ("edge", "END", "+", [(52, 60)], 52, 60),
    ("edge", "BEGIN", "-", [(1, 9)], 1, 9),
]
MITO = _genome(100, [(11, "ATGAAACCC"), (20, "AGA")])
MITO_MODELS = [("MT", "MTG", "+", [(11, 19)], 11, 22)]


@pytest.fixture
def harmonise(write_text):
    def _harmonise(models, sequences, policy=None, text=None):
        path = write_text("query.gtf", text or _gtf(models))
        df = parse_annotation_for_comparison(path)
        out, audit, summary = harmonise_stop_codons(
            df, list(sequences.items()), policy or GeneticCodePolicy()
        )
        return df, out, audit.set_index("original_transcript_id"), summary

    return _harmonise


def _cds(df, tid):
    rows = df[(df["Feature"] == "CDS") & (df["original_transcript_id"] == tid)]
    return sorted(zip(rows["Start"], rows["End"]))


def _outcome(audit, tid):
    return audit.loc[tid, "status"], audit.loc[tid, "reason"]


def test_both_strands_extended_without_phase(harmonise):
    # Phase is "." on every CDS line: the frame comes from the sequence.
    df, out, audit, summary = harmonise(CHR1_MODELS, {"1": CHR1})
    assert _outcome(audit, "P.t") == ("extended", "stop_codon_added")
    assert _cds(out, "P.t") == [(100, 106), (200, 209)]  # only the terminal segment
    assert _outcome(audit, "M.t") == ("extended", "stop_codon_added")
    assert _cds(out, "M.t") == [(497, 509)]
    assert _cds(out, "MM.t") == [(797, 806), (900, 906)]
    assert audit.loc["MM.t", "next_codon"] == "TGA"
    assert audit.loc["MM.t", "last_codon"] == "GGG"
    # Inferred exons follow the CDS; the inferred transcript span grows with it.
    exons = out[(out["Feature"] == "exon") & (out["original_transcript_id"] == "P.t")]
    assert sorted(zip(exons["Start"], exons["End"])) == [(100, 106), (200, 209)]
    tx = out[
        (out["Feature"] == "transcript") & (out["original_transcript_id"] == "P.t")
    ]
    assert (tx["Start"].iloc[0], tx["End"].iloc[0]) == (100, 209)
    # Identities, parents and row order are untouched.
    for column in ("gene_id", "transcript_id", "Parent", "Feature", "record_index"):
        assert out[column].tolist() == df[column].tolist()
    assert summary["extended"] == 3  # P, M, MM (G needs table 2)


@pytest.mark.parametrize(
    "tid, status, reason",
    [
        ("I.t", "unchanged", "stop_already_in_cds"),
        ("D.t", "skipped", "cds_length_not_multiple_of_3"),
        ("S.t", "skipped", "internal_stop_codon"),
        ("N.t", "skipped", "no_start_codon"),
        ("A1.t", "skipped", "ambiguous_bases_in_cds"),
        ("A2.t", "skipped", "ambiguous_bases_after_cds"),
        ("U.t", "unchanged", "no_adjacent_stop"),
        ("G.t", "unchanged", "no_adjacent_stop"),
    ],
)
def test_unchanged_and_skipped_models(harmonise, tid, status, reason):
    _, out, audit, _ = harmonise(CHR1_MODELS, {"1": CHR1})
    assert _outcome(audit, tid) == (status, reason)
    assert audit.loc[tid, "new_terminal_cds_end"] == audit.loc[tid, "terminal_cds_end"]
    # Eligible = frame established; later outcomes can still skip the model.
    assert audit.loc[tid, "eligible"] == (reason not in FRAME_CHECKS)


def test_disrupted_frame_with_adjacent_stop_is_not_extended(harmonise):
    # The three bases after D's CDS are TAA, but the 10-bp CDS has no consistent
    # frame: an adjacent stop triplet alone is not enough.
    df, out, audit, _ = harmonise(CHR1_MODELS, {"1": CHR1})
    assert CHR1[1210:1213] == "TAA"
    assert _cds(out, "D.t") == _cds(df, "D.t")


def test_repeated_application_does_not_extend_again(harmonise, tmp_path):
    df, once, audit, summary = harmonise(CHR1_MODELS, {"1": CHR1})
    twice, audit2, summary2 = harmonise_stop_codons(
        once, [("1", CHR1)], GeneticCodePolicy()
    )[0:3]
    assert summary2["extended"] == 0
    extended = audit.index[audit["status"] == "extended"]
    assert set(audit2.set_index("original_transcript_id").loc[extended, "reason"]) == {
        "stop_already_in_cds"
    }
    assert twice[["Start", "End"]].equals(once[["Start", "End"]])


def test_sequence_boundaries(harmonise):
    _, out, audit, _ = harmonise(EDGE_MODELS, {"edge": EDGE})
    assert _outcome(audit, "END.t") == ("skipped", "extension_out_of_bounds")
    assert _outcome(audit, "BEGIN.t") == ("skipped", "extension_out_of_bounds")


def test_missing_sequence_is_reported(harmonise):
    _, _, audit, summary = harmonise(EDGE_MODELS, {})
    assert set(audit["reason"]) == {"sequence_not_available"}
    assert summary["eligible"] == 0


def test_genetic_code_selection(harmonise):
    _, _, audit, _ = harmonise(
        CHR1_MODELS, {"1": CHR1}, parse_genetic_code_policy(2, None)
    )
    # Table 2: AGA is a stop and TGA is not; P (followed by TAA) still extends.
    assert _outcome(audit, "G.t") == ("extended", "stop_codon_added")
    assert audit.loc["G.t", "genetic_code"] == 2
    assert _outcome(audit, "MM.t") == ("unchanged", "no_adjacent_stop")
    assert _outcome(audit, "S.t") == ("skipped", "internal_stop_codon")  # TAA


def test_organelle_sequence_needs_an_explicit_code(harmonise):
    _, out, audit, summary = harmonise(MITO_MODELS, {"MT": MITO})
    assert _outcome(audit, "MTG.t") == ("skipped", "genetic_code_not_specified")
    assert audit.loc["MTG.t", "genetic_code_source"] == (
        "organelle_sequence_without_override"
    )
    _, out, audit, _ = harmonise(
        MITO_MODELS, {"MT": MITO}, parse_genetic_code_policy(1, ["MT=2"])
    )
    assert _outcome(audit, "MTG.t") == ("extended", "stop_codon_added")
    assert audit.loc["MTG.t", "genetic_code_source"] == "sequence_override"
    _, _, audit, _ = harmonise(
        CHR1_MODELS, {"1": CHR1}, parse_genetic_code_policy(1, ["1=none"])
    )
    assert set(audit["reason"]) == {"genetic_code_not_specified"}


@pytest.mark.parametrize(
    "default, overrides, message",
    [
        (27, None, "Unsupported genetic code 27"),
        (1, ["MT"], "Expected SEQNAME=TABLE"),
        (1, ["MT=two"], "not a number"),
        ("none", None, "cannot be none"),
    ],
)
def test_invalid_genetic_code_options(default, overrides, message):
    with pytest.raises(ValueError, match=message):
        parse_genetic_code_policy(default, overrides)


EXPLICIT = _genome(
    3000,
    [
        (2001, "ATGAAACCC"),  # X1: stop inside the exon's 3' UTR
        (2010, "TAA"),
        (2101, "ATGAAACCC"),  # X2: CDS ends at exon end, next exon follows
        (2110, "TAA"),  # intronic triplet: must not be used
        (2301, "ATGAAACCC"),  # X3: exon ends with the CDS, no further exon
        (2310, "TAA"),
        (2401, "ATGAAACCC"),  # X4: CDS-only, explicit mRNA ends at CDS end
        (2410, "TAA"),
    ],
)


def _explicit_records():
    def model(gene, exons, cds, mrna=None):
        start = min(s for s, _ in exons or cds)
        end = max(e for _, e in exons or cds)
        records = [
            ("1", "gene", start, end + 3, "+", f"ID={gene}"),
            ("1", "mRNA", *(mrna or (start, end)), "+", f"ID={gene}.t;Parent={gene}"),
        ]
        records += [("1", "exon", s, e, "+", f"Parent={gene}.t") for s, e in exons]
        return records + [("1", "CDS", s, e, "+", f"Parent={gene}.t") for s, e in cds]

    return (
        model("X1", [(2001, 2020)], [(2001, 2009)])
        + model("X2", [(2101, 2109), (2201, 2220)], [(2101, 2109)])
        + model("X3", [(2301, 2309)], [(2301, 2309)])
        + model("X4", [], [(2401, 2409)], mrna=(2401, 2409))
    )


def test_explicit_exons_and_split_stops(write_gff3):
    df = parse_annotation_for_comparison(write_gff3("x.gff3", _explicit_records()))
    out, audit, _ = harmonise_stop_codons(df, [("1", EXPLICIT)], GeneticCodePolicy())
    audit = audit.set_index("original_transcript_id")
    assert _outcome(audit, "X1.t") == ("extended", "stop_codon_added")
    assert audit.loc["X1.t", "exon_representation"] == "explicit"
    assert audit.loc["X1.t", "records_changed"] == "CDS"
    assert _outcome(audit, "X2.t") == ("skipped", "split_stop_codon_unsupported")
    assert _outcome(audit, "X3.t") == ("skipped", "extension_beyond_explicit_exons")
    assert _outcome(audit, "X4.t") == (
        "skipped",
        "extension_outside_transcript_record",
    )
    exons = df["Feature"] == "exon"
    assert out.loc[exons, ["Start", "End"]].equals(df.loc[exons, ["Start", "End"]])


def test_shared_inferred_gene_span_covers_every_extension(write_text):
    text = (
        '1\tsrc\tCDS\t101\t106\t.\t+\t.\tgene_id "G"; transcript_id "G.1";\n'
        '1\tsrc\tCDS\t201\t206\t.\t+\t.\tgene_id "G"; transcript_id "G.1";\n'
        '1\tsrc\tCDS\t101\t106\t.\t+\t.\tgene_id "G"; transcript_id "G.2";\n'
        '1\tsrc\tCDS\t201\t206\t.\t+\t.\tgene_id "G"; transcript_id "G.2";\n'
    )
    df = parse_annotation_for_comparison(write_text("g.gtf", text))
    out, audit, _ = harmonise_stop_codons(df, [("1", CHR1)], GeneticCodePolicy())
    assert set(audit["status"]) == {"extended"}
    gene = out[out["Feature"] == "gene"]
    assert (gene["Start"].iloc[0], gene["End"].iloc[0]) == (100, 209)


def test_iter_sequences_streams_requested_names(write_text):
    path = write_text("g.fa.gz", ">1 x\nacgt\nAC\n>2\nGG\n>3\nTT\n", compress=True)
    assert list(iter_sequences(path, {"1", "3"})) == [("1", "ACGTAC"), ("3", "TT")]


# Runner ------------------------------------------------------------------------

REFERENCE = [
    ("1", "gene", 101, 209, "+", "ID=R;biotype=protein_coding"),
    ("1", "mRNA", 101, 209, "+", "ID=R.1;Parent=R;biotype=protein_coding"),
    ("1", "exon", 101, 106, "+", "Parent=R.1"),
    ("1", "exon", 201, 209, "+", "Parent=R.1"),
    ("1", "CDS", 101, 106, "+", "Parent=R.1"),
    ("1", "CDS", 201, 209, "+", "Parent=R.1"),
]


@pytest.fixture
def runner_inputs(write_gff3, write_text):
    return {
        "ref": write_gff3("ref.gff3", REFERENCE),
        "query": write_text("query.gtf", _gtf(CHR1_MODELS)),
        "genome": write_text("genome.fa", ">1\n" + CHR1 + "\n"),
    }


def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def test_runner_default_compares_original_coordinates(runner_inputs, tmp_path):
    out = tmp_path / "plain"
    args = ["--query", runner_inputs["query"], "--reference", runner_inputs["ref"]]
    _run([*args, "--genome", runner_inputs["genome"], "--outdir", str(out)])
    summary = json.loads((out / "comparison_summary.json").read_text())
    assert summary["query_stop_codon_harmonisation"] == {
        "applied": False,
        "annotation": "query",
    }
    assert summary["cds_exact_one_to_one"]["true_positives"] == 0
    assert not (out / "stop_codon_harmonisation.tsv").exists()
    assert not (out / "query_evaluated.stop_harmonised.gtf").exists()


def test_runner_add_stop_codon(runner_inputs, tmp_path, write_text):
    before = {k: _sha(v) for k, v in runner_inputs.items()}
    out = tmp_path / "harmonised"
    _run(
        [
            "--query",
            runner_inputs["query"],
            "--reference",
            runner_inputs["ref"],
            "--genome",
            runner_inputs["genome"],
            "--add-stop-codon",
            "--outdir",
            str(out),
        ]
    )
    assert {k: _sha(v) for k, v in runner_inputs.items()} == before
    summary = json.loads((out / "comparison_summary.json").read_text())
    stop = summary["query_stop_codon_harmonisation"]
    assert (stop["applied"], stop["extended"], stop["transcripts_with_cds"]) == (
        True,
        3,
        11,
    )
    assert summary["cds_exact_one_to_one"]["true_positives"] == 1
    manifest = json.loads((out / "comparison_manifest.json").read_text())
    record = manifest["stop_codon_harmonisation"]
    assert record["policy"]["genetic_code"]["default_table"] == 1
    assert record["outputs"]["evaluated_query"]["sha256"] == _sha(
        out / "query_evaluated.stop_harmonised.gtf"
    )
    assert manifest["options"]["add_stop_codon"] is True
    audit = (out / "stop_codon_harmonisation.tsv").read_text().splitlines()
    assert len(audit) == 12
    row = dict(zip(audit[0].split("\t"), audit[1].split("\t")))
    assert (row["original_transcript_id"], row["terminal_cds_start"]) == ("P.t", "201")
    assert row["new_terminal_cds_end"] == "209"
    # The written GTF reproduces the coordinates that were compared.
    evaluated = parse_annotation_for_comparison(
        str(out / "query_evaluated.stop_harmonised.gtf")
    )
    assert _cds(evaluated, "P.t") == [(100, 106), (200, 209)]
    assert _cds(evaluated, "D.t") == [(1200, 1210)]
    # The reference is never changed.
    with open(out / "comparison_details.tsv") as handle:
        ref_row = [line for line in handle if line.startswith("reference\tR\t")][0]
    assert ref_row.split("\t")[3:5] == ["101", "209"]


@pytest.mark.parametrize(
    "extra, message",
    [
        (["--add-stop-codon"], "--add-stop-codon needs --genome"),
        (["--stop-codon-genetic-code", "2"], "need --add-stop-codon"),
    ],
)
def test_runner_option_errors(runner_inputs, tmp_path, extra, message, capsys):
    with pytest.raises(SystemExit) as error:
        _run(
            [
                "--query",
                runner_inputs["query"],
                "--reference",
                runner_inputs["ref"],
                "--outdir",
                str(tmp_path / "x"),
                *extra,
            ]
        )
    assert message in str(error.value)

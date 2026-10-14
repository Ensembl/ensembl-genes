"""Splice-site position and intron-length binning (metrics/pairwise/splice.py; SYNTHETIC structures)."""

import pytest

from ensembl.genes.annotation_qc.metrics.pairwise import splice

THREE = ((101, 150), (201, 260), (301, 350))  # CDS 250 bp span, 160 coding bases


# --------------------------------------------------------------- splice-site geometry
def test_plus_strand_sites_and_positions():
    first, second = splice.intron_sites("+", THREE)
    assert (first.lo, first.hi, first.donor, first.acceptor, first.length) == (
        151,
        200,
        151,
        200,
        50,
    )
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
    assert downstream.position["cds_span"] == {
        "donor": 150 / 250,
        "acceptor": 200 / 250,
    }
    assert (
        upstream.position["coding"]["donor"] == 50 / 160
    )  # coding bases 301-350 upstream
    # a mirrored gene gets mirrored positions
    mirrored = tuple(sorted((1001 - e, 1001 - s) for s, e in THREE))
    plus = {
        round(i.position["cds_span"]["donor"], 6)
        for i in splice.intron_sites("+", mirrored)
    }
    minus = {
        round(i.position["cds_span"]["donor"], 6)
        for i in splice.intron_sites("-", THREE)
    }
    assert plus == minus
    assert splice.intron_sites("+", ((1, 90),)) == []  # single segment: no CDS intron


def test_intron_length_bins_cover_every_length():
    edges = splice.INTRON_LENGTH_EDGES
    assert [
        splice.length_bin(n) for n in (0, 9, 10, 29, 30, 59, 60, 99_999, 100_000, 10**7)
    ] == [0, 0, 1, 1, 2, 2, 3, len(edges) - 2, len(edges) - 1, len(edges) - 1]
    assert splice.length_bin_label(0) == "< 10 bp"
    assert splice.length_bin_label(3) == "60 bp–100 bp"
    assert splice.length_bin_label(len(edges) - 1) == "≥ 100 kb"
    assert splice.position_bin(0.999) == 9 and splice.position_bin(0.0) == 0


def _gene(gid, chrom, strand, *transcripts):
    return splice.ReferenceGene(gid, chrom, strand, list(transcripts))


def _cell(result, coordinate, kind, b, min_length=0):
    return next(
        r
        for r in result["position"]
        if (r["min_intron_length"], r["coordinate"], r["kind"], r["bin"])
        == (min_length, coordinate, kind, b)
    )


def test_support_missing_loci_and_empty_bins():
    # query: chr1 exact; chr2 (minus) first intron with a shifted donor (151-204)
    query = splice.SupportIndex.build(
        [
            ("chr1", "+", THREE),
            ("chr2", "-", ((101, 150), (205, 260), (301, 350))),
            ("chr3", "+", THREE),
        ]
    )
    genes = [
        _gene("G1", "chr1", "+", THREE),
        _gene("G2", "chr2", "-", THREE),
        _gene("G3", "chr2", "-", ((601, 650), (701, 760))),  # locus not recovered
    ]
    detected = {"G1": True, "G2": True, "G3": False}
    out = splice.analyse(
        genes, query, detected, {"all": None, "chr1": {"chr1"}}, one_per_gene=True
    )
    whole = out["all"]
    length = {r["bin"]: r for r in whole["length"]}[splice.length_bin(50)]
    assert (
        length["eligible"],
        length["exact"],
        length["donor"],
        length["acceptor"],
        length["missing_locus"],
    ) == (4, 3, 3, 4, 1)
    # donor bin 6 holds G1's second donor (supported) and G2's 151-200 donor at 200 (not supported)
    cell = _cell(whole, "cds_span", "donor", 6)
    assert (cell["assessed"], cell["supported"], cell["missing_locus"]) == (2, 1, 0)
    # empty bins are explicit zeros (shown as not available), never dropped
    empty = _cell(whole, "cds_span", "donor", 0)
    assert (empty["assessed"], empty["supported"], empty["missing_locus"]) == (0, 0, 0)
    assert (
        len(whole["position"])
        == len(splice.POSITION_MIN_LENGTHS) * 2 * 2 * splice.POSITION_BINS
    )
    donors = [
        r
        for r in whole["position"]
        if r["kind"] == "donor" and r["coordinate"] == "coding"
    ]
    assert sum(r["missing_locus"] for r in donors if r["min_intron_length"] == 0) == 1
    # every intron here is 50 bp: the "introns >= 60 bp" rows count nothing
    assert not any(
        r["assessed"] + r["missing_locus"]
        for r in donors
        if r["min_intron_length"] == 60
    )
    # unique-intron totals; the scope applies to query introns too (chr3 is outside both)
    assert (
        whole["reference_introns"],
        whole["query_introns"],
        whole["shared_introns"],
    ) == (5, 6, 3)
    chr1 = out["chr1"]
    assert (
        chr1["reference_introns"],
        chr1["query_introns"],
        chr1["shared_introns"],
        chr1["genes"],
    ) == (2, 2, 2, 1)
    assert sum(r["missing_locus"] for r in chr1["length"]) == 0


def test_alternative_isoforms_are_deduplicated_or_refused():
    alt = ((101, 150), (201, 260))  # shares intron 151-200 with THREE
    genes = [_gene("G1", "chr1", "+", THREE, alt)]
    out = splice.analyse(
        genes, splice.SupportIndex(), {"G1": True}, {"all": None}, one_per_gene=False
    )
    assert out["all"]["reference_introns"] == 2 and "position" not in out["all"]
    with pytest.raises(ValueError, match="one transcript per gene"):
        splice.analyse(
            genes, splice.SupportIndex(), {"G1": True}, {"all": None}, one_per_gene=True
        )
    with pytest.raises(ValueError, match="no locus outcome"):
        splice.analyse(
            genes, splice.SupportIndex(), {}, {"all": None}, one_per_gene=False
        )

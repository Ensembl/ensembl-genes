"""
CDS position coordinates of splice sites (metrics/pairwise/splice.py; SYNTHETIC structures).

The dashboard's "genomic CDS" / "coding sequence" toggle switches between the comparator's
two position coordinates; these tests pin down how they differ and that both follow the
direction of transcription.
"""

import pytest

from ensembl.genes.annotation_qc.metrics.pairwise import splice

# Long first intron, short second intron: the two coordinates disagree strongly.
LONG_FIRST = ((1_001, 1_100), (11_001, 11_100), (11_201, 11_300))  # 300 coding bases


@pytest.mark.parametrize("strand", ["+", "-"])
def test_genomic_and_coding_positions_are_distinct(strand):
    introns = sorted(splice.intron_sites(strand, LONG_FIRST), key=lambda i: i.lo)
    span = LONG_FIRST[-1][1] - LONG_FIRST[0][0] + 1
    for intron in introns:
        pos = intron.position
        # coding: one junction, so donor and acceptor share a position
        assert pos["coding"]["donor"] == pos["coding"]["acceptor"]
        # genomic CDS span: the acceptor lies exactly one intron length after the donor
        assert pos["cds_span"]["acceptor"] - pos["cds_span"]["donor"] == pytest.approx(
            intron.length / span
        )
        for coordinate in splice.POSITION_COORDINATES:
            for kind in splice.SITE_KINDS:
                assert 0 < pos[coordinate][kind] < 1
    first_in_transcription = introns[0] if strand == "+" else introns[-1]
    if strand == "+":
        # 100 coding bases (of 300) precede the long intron's donor, and 100 genomic bases
        assert first_in_transcription.position["coding"]["donor"] == pytest.approx(100 / 300)
        assert first_in_transcription.position["cds_span"]["donor"] == pytest.approx(100 / span)
    else:
        # minus strand: transcription starts at 11,300, so the short intron comes first
        assert first_in_transcription.lo == 11_101
        assert first_in_transcription.position["coding"]["donor"] == pytest.approx(100 / 300)
        assert first_in_transcription.position["cds_span"]["donor"] == pytest.approx(100 / span)


def test_positions_follow_transcription_order():
    for strand in "+-":
        introns = splice.intron_sites(strand, LONG_FIRST)
        ordered = sorted(introns, key=lambda i: i.lo, reverse=strand == "-")
        for coordinate in splice.POSITION_COORDINATES:
            donors = [i.position[coordinate]["donor"] for i in ordered]
            assert donors == sorted(donors), (strand, coordinate)


def test_segment_order_in_the_file_does_not_matter():
    shuffled = "11201-11300,1001-1100,11001-11100"
    assert splice.parse_segments(shuffled) == LONG_FIRST
    assert [i.position for i in splice.intron_sites("-", splice.parse_segments(shuffled))] == [
        i.position for i in splice.intron_sites("-", LONG_FIRST)
    ]


def test_one_base_segments_and_bin_edges():
    tiny = ((10, 10), (20, 20))  # one coding base each side of a 9 bp gap
    (intron,) = splice.intron_sites("+", tiny)
    assert (intron.lo, intron.hi, intron.length) == (11, 19, 9)
    assert intron.position["coding"]["donor"] == 0.5
    assert intron.position["cds_span"] == {"donor": 1 / 11, "acceptor": 10 / 11}
    assert [splice.position_bin(v) for v in (1e-9, 0.1, 0.5, 0.999999)] == [0, 1, 5, 9]


def test_binning_uses_each_coordinate_separately():
    gene = splice.ReferenceGene("g", "chr1", "-", [LONG_FIRST])
    support = splice.SupportIndex.build([])  # nothing supported
    result = splice.analyse([gene], support, {"g": True}, {"all": None}, one_per_gene=True)["all"]
    cells = {
        (r["coordinate"], r["kind"], r["bin"]): r["assessed"]
        for r in result["position"]
        if r["min_intron_length"] == 0 and r["assessed"]
    }
    # minus strand, transcription order: the short intron (donor: coding 1/3, genomic
    # 100/span), then the long intron (donor: coding 2/3, genomic 300/span; its 9.9 kb
    # follow the donor). Genomic positions put both donors in the first tenth.
    span = LONG_FIRST[-1][1] - LONG_FIRST[0][0] + 1
    expected = {}
    for i in splice.intron_sites("-", LONG_FIRST):
        for c in splice.POSITION_COORDINATES:
            for k in splice.SITE_KINDS:
                key = (c, k, splice.position_bin(i.position[c][k]))
                expected[key] = expected.get(key, 0) + 1
    assert cells == expected
    assert cells[("coding", "donor", 3)] == 1 and cells[("coding", "donor", 6)] == 1
    assert splice.position_bin(300 / span) == 0
    assert cells[("cds_span", "donor", 0)] == 2

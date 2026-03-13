"""
Tests for the coordinate- and sequence-based metric functions in metrics.features.

GFF input is mocked as plain pandas DataFrames (pyranges objects support the
same boolean indexing, so this works without a real annotation file).

FASTA input is mocked with a minimal in-memory class that mimics the
pyfaidx.Fasta slice interface.
"""

import pandas as pd
import pytest

from metrics.features import (
    compute_cds_utr5_overlap,
    compute_cds_utr5_overlap_summary,
    compute_splice_junction_metrics,
    compute_splice_junction_summary_stats,
    compute_translation_metrics,
    compute_translation_summary_stats,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gff(rows: list[dict]) -> pd.DataFrame:
    """Build a minimal GFF-like DataFrame from a list of row dicts."""
    df = pd.DataFrame(rows)
    for col in ("Feature", "transcript_id", "Start", "End", "Strand", "Chromosome", "Parent", "gene_id"):
        if col not in df.columns:
            df[col] = None
    return df


class _Interval:
    """Mimics a pyfaidx interval: .seq, .reverse.complement.seq"""
    _COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

    def __init__(self, seq: str):
        self.seq = seq

    @property
    def reverse(self) -> "_Interval":
        return _Interval(self.seq[::-1])

    @property
    def complement(self) -> "_Interval":
        return _Interval(self.seq.translate(self._COMP))


class MockFasta:
    """Minimal pyfaidx.Fasta substitute backed by an in-memory dict."""

    def __init__(self, sequences: dict[str, str]):
        self._seqs = sequences

    def __getitem__(self, chrom: str):
        seq = self._seqs[chrom]

        class _Chrom:
            def __getitem__(_, sl):  # noqa: N805
                return _Interval(seq[sl])

        return _Chrom()


# ---------------------------------------------------------------------------
# compute_cds_utr5_overlap
# ---------------------------------------------------------------------------

def _base_gff(tx_id, gene_id, cds_start, cds_end, utr5_start=None, utr5_end=None):
    """Return a minimal GFF with one mRNA + one CDS (+ optional UTR5)."""
    rows = [
        {"Feature": "mRNA",          "transcript_id": tx_id,  "Parent": f"gene:{gene_id}", "Start": 0, "End": 1000},
        {"Feature": "CDS",           "transcript_id": tx_id,  "Parent": f"transcript:{tx_id}", "Start": cds_start,  "End": cds_end},
    ]
    if utr5_start is not None:
        rows.append({"Feature": "five_prime_UTR", "transcript_id": tx_id, "Parent": f"transcript:{tx_id}", "Start": utr5_start, "End": utr5_end})
    return _gff(rows)


def test_overlap_no_utr5_features():
    gff = _base_gff("TX1", "G1", cds_start=100, cds_end=200)
    df = compute_cds_utr5_overlap(gff)
    assert len(df) == 1
    assert df.loc[0, "n_cds_overlapping_utr5"] == 0
    assert not df.loc[0, "has_overlap"]


def test_overlap_adjacent_no_overlap():
    # CDS starts exactly where UTR5 ends — half-open: no overlap
    gff = _base_gff("TX1", "G1", cds_start=100, cds_end=200, utr5_start=0, utr5_end=100)
    df = compute_cds_utr5_overlap(gff)
    assert df.loc[0, "n_cds_overlapping_utr5"] == 0
    assert not df.loc[0, "has_overlap"]


def test_overlap_detected():
    # CDS starts inside the UTR5 region
    gff = _base_gff("TX1", "G1", cds_start=90, cds_end=200, utr5_start=0, utr5_end=100)
    df = compute_cds_utr5_overlap(gff)
    assert df.loc[0, "n_cds_overlapping_utr5"] == 1
    assert df.loc[0, "has_overlap"]


def test_overlap_empty_cds():
    gff = _gff([{"Feature": "five_prime_UTR", "transcript_id": "TX1", "Start": 0, "End": 100}])
    df = compute_cds_utr5_overlap(gff)
    assert df.empty


def test_overlap_gene_id_propagated():
    gff = _base_gff("TX1", "MYGENE", cds_start=90, cds_end=200, utr5_start=0, utr5_end=100)
    df = compute_cds_utr5_overlap(gff)
    assert df.loc[0, "gene_id"] == "MYGENE"


def test_overlap_mixed_transcripts():
    # TX1 overlaps; TX2 does not (CDS starts at UTR5 end exactly)
    rows = [
        {"Feature": "mRNA",            "transcript_id": "TX1", "Parent": "gene:G1", "Start": 0, "End": 1000},
        {"Feature": "CDS",             "transcript_id": "TX1", "Parent": "transcript:TX1", "Start": 90, "End": 200},
        {"Feature": "five_prime_UTR",  "transcript_id": "TX1", "Parent": "transcript:TX1", "Start": 0,  "End": 100},
        {"Feature": "mRNA",            "transcript_id": "TX2", "Parent": "gene:G2", "Start": 0, "End": 1000},
        {"Feature": "CDS",             "transcript_id": "TX2", "Parent": "transcript:TX2", "Start": 100, "End": 200},
        {"Feature": "five_prime_UTR",  "transcript_id": "TX2", "Parent": "transcript:TX2", "Start": 0,   "End": 100},
    ]
    df = compute_cds_utr5_overlap(_gff(rows)).set_index("transcript_id")
    assert df.loc["TX1", "has_overlap"]
    assert not df.loc["TX2", "has_overlap"]


def test_overlap_n_cds_segments_counted():
    # Two CDS segments for the same transcript; only one overlaps UTR5
    rows = [
        {"Feature": "mRNA",           "transcript_id": "TX1", "Parent": "gene:G1", "Start": 0, "End": 1000},
        {"Feature": "CDS",            "transcript_id": "TX1", "Parent": "transcript:TX1", "Start": 90,  "End": 150},
        {"Feature": "CDS",            "transcript_id": "TX1", "Parent": "transcript:TX1", "Start": 200, "End": 300},
        {"Feature": "five_prime_UTR", "transcript_id": "TX1", "Parent": "transcript:TX1", "Start": 0,   "End": 100},
    ]
    df = compute_cds_utr5_overlap(_gff(rows))
    assert df.loc[0, "n_cds_segments"] == 2
    assert df.loc[0, "n_cds_overlapping_utr5"] == 1


def test_overlap_cross_transcript_nte():
    # NTE scenario: TX_long (same gene) has CDS starting upstream of TX_short's CDS,
    # inside the region annotated as five_prime_UTR on TX_short.
    # TX_short: UTR5 [0,100), CDS [100,200)
    # TX_long:  CDS [50,200) — starts inside TX_short's UTR5 → overlap expected
    rows = [
        {"Feature": "mRNA",           "transcript_id": "TX_short", "Parent": "gene:G1", "Start": 0, "End": 1000},
        {"Feature": "CDS",            "transcript_id": "TX_short", "Parent": "transcript:TX_short", "Start": 100, "End": 200},
        {"Feature": "five_prime_UTR", "transcript_id": "TX_short", "Parent": "transcript:TX_short", "Start": 0,   "End": 100},
        {"Feature": "mRNA",           "transcript_id": "TX_long",  "Parent": "gene:G1", "Start": 0, "End": 1000},
        {"Feature": "CDS",            "transcript_id": "TX_long",  "Parent": "transcript:TX_long",  "Start": 50,  "End": 200},
    ]
    df = compute_cds_utr5_overlap(_gff(rows)).set_index("transcript_id")
    # TX_long's CDS [50,200) overlaps TX_short's UTR5 [0,100)
    assert df.loc["TX_long", "has_overlap"]
    # TX_short's CDS [100,200) is adjacent to its own UTR5 [0,100) — no overlap
    assert not df.loc["TX_short", "has_overlap"]


def test_overlap_cross_transcript_no_cross_gene_bleed():
    # UTR5 from a different gene must not trigger overlap in an unrelated gene
    rows = [
        # Gene G1: TX1 has CDS at [50,200) — no UTR5 anywhere in G1
        {"Feature": "mRNA", "transcript_id": "TX1", "Parent": "gene:G1", "Start": 0, "End": 1000},
        {"Feature": "CDS",  "transcript_id": "TX1", "Parent": "transcript:TX1", "Start": 50, "End": 200},
        # Gene G2: TX2 has UTR5 [0,100) that spans the same coordinates
        {"Feature": "mRNA",           "transcript_id": "TX2", "Parent": "gene:G2", "Start": 0, "End": 1000},
        {"Feature": "CDS",            "transcript_id": "TX2", "Parent": "transcript:TX2", "Start": 200, "End": 300},
        {"Feature": "five_prime_UTR", "transcript_id": "TX2", "Parent": "transcript:TX2", "Start": 0,   "End": 100},
    ]
    df = compute_cds_utr5_overlap(_gff(rows)).set_index("transcript_id")
    assert not df.loc["TX1", "has_overlap"]


# ---------------------------------------------------------------------------
# compute_translation_metrics
# ---------------------------------------------------------------------------

def _cds_gff(tx_id, chrom, start, end, strand="+"):
    return _gff([
        {"Feature": "mRNA", "transcript_id": tx_id, "Parent": f"gene:G1", "Chromosome": chrom, "Start": start, "End": end, "Strand": strand},
        {"Feature": "CDS",  "transcript_id": tx_id, "Parent": f"transcript:{tx_id}", "Chromosome": chrom, "Start": start, "End": end, "Strand": strand},
    ])


def test_translation_canonical_start_stop():
    gff = _cds_gff("TX1", "chr1", 0, 9)  # ATG AAA TAA
    fasta = MockFasta({"chr1": "ATGAAATAA"})
    df = compute_translation_metrics(gff, fasta)
    row = df[df["transcript_id"] == "TX1"].iloc[0]
    assert row["start_codon"] == "ATG"
    assert row["start_status"] == "canonical"
    assert row["stop_codon"] == "TAA"
    assert row["stop_status"] == "canonical"
    assert row["frame_ok"]
    assert not row["has_internal_stop"]


def test_translation_near_cognate_start():
    gff = _cds_gff("TX1", "chr1", 0, 9)  # CTG AAA TAA
    fasta = MockFasta({"chr1": "CTGAAATAA"})
    df = compute_translation_metrics(gff, fasta)
    assert df.iloc[0]["start_status"] == "near_cognate"


def test_translation_frame_error():
    gff = _cds_gff("TX1", "chr1", 0, 8)  # 8 nt — not divisible by 3
    fasta = MockFasta({"chr1": "ATGAAATA"})
    df = compute_translation_metrics(gff, fasta)
    assert df.iloc[0]["frame_error"]


def test_translation_internal_stop():
    # ATG TAA AAA TAA — internal stop at codon 1
    seq = "ATGTAAAAA" + "TAA"
    gff = _cds_gff("TX1", "chr1", 0, len(seq))
    fasta = MockFasta({"chr1": seq})
    df = compute_translation_metrics(gff, fasta)
    assert df.iloc[0]["has_internal_stop"]


def test_translation_no_cds_returns_empty():
    gff = _gff([{"Feature": "mRNA", "transcript_id": "TX1", "Start": 0, "End": 100}])
    fasta = MockFasta({"chr1": "A" * 100})
    df = compute_translation_metrics(gff, fasta)
    assert df.empty


def test_translation_minus_strand():
    # Minus-strand CDS: genomic sequence is reverse complement of coding seq.
    # Coding seq: ATG AAA TAA (canonical). RC = TTA TTT CAT → store on - strand.
    coding = "ATGAAATAA"
    _comp = str.maketrans("ACGT", "TGCA")
    genomic = coding[::-1].translate(_comp)  # reverse complement
    gff = _cds_gff("TX1", "chr1", 0, 9, strand="-")
    fasta = MockFasta({"chr1": genomic})
    df = compute_translation_metrics(gff, fasta)
    assert df.iloc[0]["start_status"] == "canonical"
    assert df.iloc[0]["stop_status"] == "canonical"


# ---------------------------------------------------------------------------
# compute_splice_junction_metrics
# ---------------------------------------------------------------------------

def _exon_gff(tx_id, exons: list[tuple[int, int]], chrom="chr1", strand="+"):
    """Build a GFF with exon features for a single transcript."""
    rows = []
    for start, end in exons:
        rows.append({
            "Feature": "exon", "transcript_id": tx_id,
            "Parent": f"transcript:{tx_id}",
            "Chromosome": chrom, "Strand": strand, "Start": start, "End": end,
        })
    return _gff(rows)


def test_splice_junction_canonical_gt_ag():
    # Exon1: [0,100), intron: [100,200), Exon2: [200,300)
    # Place GT at 100 and AG at 198
    seq = "A" * 100 + "GT" + "N" * 96 + "AG" + "A" * 100
    gff = _exon_gff("TX1", [(0, 100), (200, 300)])
    fasta = MockFasta({"chr1": seq})
    df = compute_splice_junction_metrics(gff, fasta)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["donor_dinucleotide"] == "GT"
    assert row["acceptor_dinucleotide"] == "AG"
    assert row["junction_type"] == "GT-AG"
    assert row["is_canonical"]


def test_splice_junction_noncanonical():
    seq = "A" * 100 + "TT" + "N" * 96 + "CC" + "A" * 100
    gff = _exon_gff("TX1", [(0, 100), (200, 300)])
    fasta = MockFasta({"chr1": seq})
    df = compute_splice_junction_metrics(gff, fasta)
    row = df.iloc[0]
    assert row["junction_type"] == "TT-CC"
    assert not row["is_canonical"]


def test_splice_junction_intron_coordinates():
    gff = _exon_gff("TX1", [(0, 100), (200, 300)])
    fasta = MockFasta({"chr1": "A" * 300})
    df = compute_splice_junction_metrics(gff, fasta)
    row = df.iloc[0]
    assert row["intron_start"] == 100
    assert row["intron_end"] == 200
    assert row["chromosome"] == "chr1"
    assert row["strand"] == "+"


def test_splice_junction_intron_numbering():
    # Three exons → two introns; check numbering is 1-based
    seq = "A" * 400
    gff = _exon_gff("TX1", [(0, 100), (150, 250), (300, 400)])
    fasta = MockFasta({"chr1": seq})
    df = compute_splice_junction_metrics(gff, fasta).sort_values("intron_start")
    assert list(df["intron_number"]) == [1, 2]


def test_splice_junction_single_exon_no_introns():
    gff = _exon_gff("TX1", [(0, 300)])
    fasta = MockFasta({"chr1": "A" * 300})
    df = compute_splice_junction_metrics(gff, fasta)
    assert df.empty


def test_splice_junction_no_exons_returns_empty():
    gff = _gff([{"Feature": "mRNA", "transcript_id": "TX1", "Start": 0, "End": 300}])
    fasta = MockFasta({"chr1": "A" * 300})
    df = compute_splice_junction_metrics(gff, fasta)
    assert df.empty


def test_splice_junction_minus_strand_orientation():
    # On - strand the donor is at the genomic RIGHT end of the intron (RC'd).
    # Genomic: ...NN CT (pos 198-200) ... AC NN (pos 100-102) ...
    # After RC: donor = RC("CT") = "AG"... hmm, let me think properly.
    # We want the minus-strand junction to be GT-AG (canonical).
    # donor_nt = RC(fasta[i_end-2 : i_end]) = RC(fasta[198:200])
    # acceptor_nt = RC(fasta[i_start : i_start+2]) = RC(fasta[100:102])
    # For GT-AG: donor="GT", acceptor="AG"
    #   → fasta[198:200] must RC to "GT" → fasta[198:200] = RC("GT") = "AC"
    #   → fasta[100:102] must RC to "AG" → fasta[100:102] = RC("AG") = "CT"
    seq = "A" * 100 + "CT" + "N" * 96 + "AC" + "A" * 100
    gff = _exon_gff("TX1", [(0, 100), (200, 300)], strand="-")
    fasta = MockFasta({"chr1": seq})
    df = compute_splice_junction_metrics(gff, fasta)
    row = df.iloc[0]
    assert row["donor_dinucleotide"] == "GT"
    assert row["acceptor_dinucleotide"] == "AG"
    assert row["is_canonical"]


# ---------------------------------------------------------------------------
# compute_cds_utr5_overlap_summary
# ---------------------------------------------------------------------------

def test_cds_utr5_overlap_summary_empty():
    s = compute_cds_utr5_overlap_summary(pd.DataFrame(columns=["has_overlap"]))
    assert s["Coding transcripts assessed"] == 0


def test_cds_utr5_overlap_summary_no_overlaps():
    df = pd.DataFrame({"has_overlap": [False, False, False]})
    s = compute_cds_utr5_overlap_summary(df)
    assert s["Coding transcripts assessed"] == 3
    assert s["With CDS/5' UTR overlap"] == 0
    assert s["% with CDS/5' UTR overlap"] == 0.0


def test_cds_utr5_overlap_summary_some_overlaps():
    df = pd.DataFrame({"has_overlap": [True, False, True, True]})
    s = compute_cds_utr5_overlap_summary(df)
    assert s["Coding transcripts assessed"] == 4
    assert s["With CDS/5' UTR overlap"] == 3
    assert s["% with CDS/5' UTR overlap"] == 75.0


# ---------------------------------------------------------------------------
# compute_translation_summary_stats
# ---------------------------------------------------------------------------

def _translation_df(**kwargs):
    """Build a minimal translation_metrics DataFrame with given column values."""
    defaults = {
        "start_status": "canonical",
        "stop_status": "canonical",
        "frame_ok": True,
        "frame_error": False,
        "has_internal_stop": False,
    }
    defaults.update(kwargs)
    return pd.DataFrame([defaults])


def test_translation_summary_canonical():
    df = _translation_df()
    s = compute_translation_summary_stats(df)
    assert s["Coding transcripts assessed"] == 1
    assert s["Canonical start (ATG)"] == 1
    assert s["% canonical start"] == 100.0
    assert s["Canonical stop"] == 1
    assert s["Frame errors"] == 0
    assert s["With internal stops"] == 0


def test_translation_summary_near_cognate():
    df = _translation_df(start_status="near_cognate")
    s = compute_translation_summary_stats(df)
    assert s["Near-cognate start"] == 1
    assert s["% near-cognate start"] == 100.0
    assert s["Canonical start (ATG)"] == 0


def test_translation_summary_frame_and_internal_stop():
    rows = [
        {"start_status": "canonical", "stop_status": "canonical",  "frame_ok": True,  "frame_error": False, "has_internal_stop": False},
        {"start_status": "canonical", "stop_status": "canonical",  "frame_ok": False, "frame_error": True,  "has_internal_stop": True},
    ]
    df = pd.DataFrame(rows)
    s = compute_translation_summary_stats(df)
    assert s["Coding transcripts assessed"] == 2
    assert s["Frame errors"] == 1
    assert s["% frame errors"] == 50.0
    assert s["With internal stops"] == 1
    assert s["% with internal stops"] == 50.0


def test_translation_summary_empty():
    df = pd.DataFrame(columns=["start_status", "stop_status", "frame_ok", "frame_error", "has_internal_stop"])
    s = compute_translation_summary_stats(df)
    assert s["Coding transcripts assessed"] == 0


# ---------------------------------------------------------------------------
# compute_splice_junction_summary_stats
# ---------------------------------------------------------------------------

def _splice_df(junction_types: list[str]) -> pd.DataFrame:
    rows = [{"junction_type": jt, "is_canonical": jt in {"GT-AG", "GC-AG", "AT-AC"}} for jt in junction_types]
    return pd.DataFrame(rows)


def test_splice_summary_all_canonical():
    df = _splice_df(["GT-AG", "GT-AG", "GC-AG"])
    s = compute_splice_junction_summary_stats(df)
    assert s["Introns assessed"] == 3
    assert s["Canonical junctions"] == 3
    assert s["% canonical junctions"] == 100.0
    assert s["Noncanonical junctions"] == 0


def test_splice_summary_mixed():
    df = _splice_df(["GT-AG", "GT-AG", "TT-CC", "AT-AC"])
    s = compute_splice_junction_summary_stats(df)
    assert s["Introns assessed"] == 4
    assert s["Canonical junctions"] == 3
    assert s["Noncanonical junctions"] == 1
    assert s["% noncanonical junctions"] == 25.0
    assert s["GT-AG junctions"] == 2
    assert s["AT-AC junctions"] == 1


def test_splice_summary_empty():
    df = pd.DataFrame(columns=["junction_type", "is_canonical"])
    s = compute_splice_junction_summary_stats(df)
    assert s["Introns assessed"] == 0

"""
Tests for metrics.events — pure functions only, no I/O or mocking needed.
"""

import pytest

from metrics.events import (
    assess_splice_junction,
    compute_cds_metrics,
    get_genetic_code,
    STANDARD_CODE,
)


# ---------------------------------------------------------------------------
# get_genetic_code
# ---------------------------------------------------------------------------

def test_get_genetic_code_standard():
    gc = get_genetic_code(1)
    assert gc.table_id == 1
    assert "TAA" in gc.stop_codons
    assert "TAG" in gc.stop_codons
    assert "TGA" in gc.stop_codons


def test_get_genetic_code_vertebrate_mt():
    gc = get_genetic_code(2)
    assert gc.table_id == 2
    # TGA is not a stop in vertebrate mitochondrial code
    assert "TGA" not in gc.stop_codons


def test_get_genetic_code_invalid():
    with pytest.raises(ValueError, match="Unknown NCBI genetic code"):
        get_genetic_code(999)


# ---------------------------------------------------------------------------
# compute_cds_metrics — start codon
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq, expected_status", [
    ("ATGAAATAA", "canonical"),
    ("CTGAAATAA", "near_cognate"),   # position-1 variant
    ("GTGAAATAA", "near_cognate"),
    ("TTGAAATAA", "near_cognate"),
    ("ACGAAATAA", "near_cognate"),   # position-2 variant
    ("ATAAAATAA", "non_cognate"),    # wobble-position variant — NOT near cognate
    ("CCCAAATAA", "non_cognate"),
])
def test_start_codon_status(seq, expected_status):
    assert compute_cds_metrics(seq).start_status == expected_status


def test_start_codon_identity():
    m = compute_cds_metrics("CTGAAATAA")
    assert m.start_codon == "CTG"


# ---------------------------------------------------------------------------
# compute_cds_metrics — stop codon
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq, expected_stop_codon, expected_status", [
    ("ATGAAATAA", "TAA", "canonical"),
    ("ATGAAATAG", "TAG", "canonical"),
    ("ATGAAATGA", "TGA", "canonical"),
    ("ATGAAAAAA", "AAA", "noncanonical"),
])
def test_stop_codon(seq, expected_stop_codon, expected_status):
    m = compute_cds_metrics(seq)
    assert m.stop_codon == expected_stop_codon
    assert m.stop_status == expected_status


# ---------------------------------------------------------------------------
# compute_cds_metrics — missing (too short)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq", ["", "AT", None])
def test_missing_codons(seq):
    m = compute_cds_metrics(seq)
    assert m.start_status == "missing"
    assert m.stop_status == "missing"
    assert m.start_codon == ""
    assert m.stop_codon == ""


# ---------------------------------------------------------------------------
# compute_cds_metrics — frame
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq, frame_ok", [
    ("ATGAAATAA", True),   # 9 nt
    ("ATGAAATA",  False),  # 8 nt
    ("ATGAA",     False),  # 5 nt
    ("",          True),   # 0 nt — divisible by 3
])
def test_frame(seq, frame_ok):
    m = compute_cds_metrics(seq)
    assert m.frame_ok is frame_ok
    assert m.frame_error is not frame_ok


# ---------------------------------------------------------------------------
# compute_cds_metrics — internal stops
# ---------------------------------------------------------------------------

def test_no_internal_stop():
    # ATG AAA TAA — no internal stop
    m = compute_cds_metrics("ATGAAATAA")
    assert m.has_internal_stop is False


def test_internal_stop_present():
    # ATG TAA AAA TAA — TAA at codon 1 is internal
    m = compute_cds_metrics("ATGTAAAAA" + "TAA")
    assert m.has_internal_stop is True


def test_terminal_stop_not_counted_as_internal():
    # Only the last codon is a stop; should not be flagged
    m = compute_cds_metrics("ATGAAATAA")
    assert m.has_internal_stop is False


def test_custom_genetic_code_stop():
    # In vertebrate mt (table 2), TGA is a sense codon, not a stop
    mt = get_genetic_code(2)
    # ATG + TGA (sense in mt) + TAA — no internal stop under mt code
    m = compute_cds_metrics("ATGTGATAA", mt)
    assert m.has_internal_stop is False


# ---------------------------------------------------------------------------
# assess_splice_junction — classification
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq, expected_type, is_canonical", [
    ("GTNNNNAG", "GT-AG", True),   # major U2
    ("GCNNNNAG", "GC-AG", True),   # minor U2
    ("ATNNNNAC", "AT-AC", True),   # U12
    ("TTNNNNCC", "TT-CC", False),  # noncanonical
    ("GGNNNNTT", "GG-TT", False),
])
def test_splice_junction_classification(seq, expected_type, is_canonical):
    m = assess_splice_junction(seq)
    assert m.junction_type == expected_type
    assert m.is_canonical is is_canonical


def test_splice_junction_dinucleotides_recorded():
    m = assess_splice_junction("GTNNNNAG")
    assert m.donor_dinucleotide == "GT"
    assert m.acceptor_dinucleotide == "AG"


def test_splice_junction_noncanonical_dinucleotides_recorded():
    m = assess_splice_junction("TTNNNNCC")
    assert m.donor_dinucleotide == "TT"
    assert m.acceptor_dinucleotide == "CC"


# ---------------------------------------------------------------------------
# assess_splice_junction — too short / empty
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq", ["", "GT", "GTN", None])
def test_splice_junction_too_short(seq):
    m = assess_splice_junction(seq)
    assert m.junction_type == "unknown"
    assert m.is_canonical is False
    assert m.donor_dinucleotide == ""
    assert m.acceptor_dinucleotide == ""


def test_splice_junction_exactly_4_nt():
    # 4-char string: donor = first 2, acceptor = last 2 (same here)
    m = assess_splice_junction("GTAG")
    assert m.donor_dinucleotide == "GT"
    assert m.acceptor_dinucleotide == "AG"
    assert m.is_canonical is True

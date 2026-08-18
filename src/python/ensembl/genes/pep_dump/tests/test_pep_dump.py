#!/usr/bin/env python3
"""
Unit tests for ensembl_pep_dump.

These cover the pure functions only: no database connection is required.
Run with:  python3 -m unittest discover -s tests -v
"""

import io
import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ensembl_pep_dump import (  # noqa: E402
    apply_seq_edit,
    build_pep_header,
    build_spliced_seq,
    build_translateable_seq,
    classify_slices,
    collect_seq_edits,
    coding_region_bounds,
    compute_cdna_coding_end,
    compute_cdna_coding_start,
    filter_grch38_y_par,
    is_start_codon,
    is_ter_codon,
    parse_seq_edit,
    reverse_complement,
    set_table1_start_codons,
    stable_id_version,
    translate_codon,
    translate_sequence,
    translate_transcript,
    write_fasta_record,
    CoreDbError,
)


def exon(exon_id, start, end, strand=1, phase=0, rank=1, end_phase=-1):
    return {
        'exon_id': exon_id, 'seq_region_id': 1,
        'seq_region_start': start, 'seq_region_end': end,
        'seq_region_strand': strand, 'phase': phase, 'end_phase': end_phase,
        'stable_id': 'ENSE{0:011d}'.format(exon_id), 'version': 1,
        'transcript_id': 1, 'rank': rank,
    }


class FakeDB:
    """
    Minimal stand-in for EnsemblCoreDB exposing only what the transcript
    functions call, so the biological logic can be tested without MySQL.
    """

    def __init__(self, sequences=None, transcript_attribs=None,
                 translation_attribs=None, codon_table=1, xrefs=None):
        self.sequences = sequences or {}
        self._transcript_attribs = transcript_attribs or {}
        self._translation_attribs = translation_attribs or {}
        self._codon_table = codon_table
        self._xrefs = xrefs or {}

    def exon_sequence(self, ex):
        seq = self.sequences[ex['exon_id']]
        return reverse_complement(seq) if ex['seq_region_strand'] == -1 else seq

    def transcript_attribs(self, transcript_id, code):
        return self._transcript_attribs.get(code, [])

    def translation_attribs(self, translation_id, code):
        return self._translation_attribs.get(code, [])

    def codon_table_id(self, seq_region_id):
        return self._codon_table

    def display_xref_label(self, xref_id):
        return self._xrefs.get(xref_id)


class TestStableIdVersion(unittest.TestCase):
    def test_version_appended_when_present(self):
        self.assertEqual(stable_id_version("ENSG00000001", 1), "ENSG00000001.1")
        self.assertEqual(stable_id_version("ENSP05220066032", 12),
                         "ENSP05220066032.12")

    def test_version_omitted_when_absent_or_zero(self):
        self.assertEqual(stable_id_version("ENSG00000001", None), "ENSG00000001")
        self.assertEqual(stable_id_version("ENSG00000001", 0), "ENSG00000001")
        self.assertEqual(stable_id_version("ENSG00000001", ""), "ENSG00000001")


class TestReverseComplement(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(reverse_complement("ATCG"), "CGAT")
        self.assertEqual(reverse_complement(""), "")

    def test_case_preserved(self):
        self.assertEqual(reverse_complement("atcg"), "cgat")
        self.assertEqual(reverse_complement("AtCg"), "cGaT")

    def test_iupac_ambiguity_codes(self):
        self.assertEqual(reverse_complement("ATCGN"), "NCGAT")
        self.assertEqual(reverse_complement("RYSWKMBDHVN"), "NBDHVKMWSRY")

    def test_bytes_input(self):
        self.assertEqual(reverse_complement(b"ATCG"), b"CGAT")

    def test_double_reverse_is_identity(self):
        seq = "ACGTRYSWKMBDHVNacgtn"
        self.assertEqual(reverse_complement(reverse_complement(seq)), seq)


class TestFastaWrapping(unittest.TestCase):
    def test_wraps_at_sixty(self):
        buf = io.StringIO()
        write_fasta_record(buf, "HDR", "A" * 130)
        lines = buf.getvalue().split('\n')
        self.assertEqual(lines[0], ">HDR")
        self.assertEqual(len(lines[1]), 60)
        self.assertEqual(len(lines[2]), 60)
        self.assertEqual(len(lines[3]), 10)
        self.assertEqual(lines[4], '')

    def test_exact_multiple_has_no_blank_line(self):
        buf = io.StringIO()
        write_fasta_record(buf, "HDR", "M" * 120)
        self.assertEqual(buf.getvalue(), ">HDR\n" + "M" * 60 + "\n" + "M" * 60 + "\n")

    def test_short_sequence(self):
        buf = io.StringIO()
        write_fasta_record(buf, "HDR", "MKV")
        self.assertEqual(buf.getvalue(), ">HDR\nMKV\n")

    def test_empty_sequence_writes_header_only(self):
        buf = io.StringIO()
        write_fasta_record(buf, "HDR", "")
        self.assertEqual(buf.getvalue(), ">HDR\n")


class TestCodonTranslation(unittest.TestCase):
    def test_standard_table(self):
        self.assertEqual(translate_codon("ATG", 1), "M")
        self.assertEqual(translate_codon("atg", 1), "M")
        self.assertEqual(translate_codon("AUG", 1), "M")
        self.assertEqual(translate_codon("TTT", 1), "F")
        self.assertEqual(translate_codon("TGG", 1), "W")
        for stop in ("TAA", "TAG", "TGA"):
            self.assertEqual(translate_codon(stop, 1), "*")

    def test_vertebrate_mitochondrial_table(self):
        self.assertEqual(translate_codon("ATA", 2), "M")
        self.assertEqual(translate_codon("TGA", 2), "W")
        self.assertEqual(translate_codon("AGA", 2), "*")
        self.assertEqual(translate_codon("AGG", 2), "*")

    def test_ambiguous_codons(self):
        # Every expansion agrees -> the shared residue.
        self.assertEqual(translate_codon("TTR", 1), "L")
        self.assertEqual(translate_codon("CTN", 1), "L")
        self.assertEqual(translate_codon("GGN", 1), "G")
        # Expansions disagree -> X.
        self.assertEqual(translate_codon("RAY", 1), "X")
        self.assertEqual(translate_codon("NNN", 1), "X")
        self.assertEqual(translate_codon("ATN", 1), "X")

    def test_translate_sequence(self):
        self.assertEqual(translate_sequence("ATGCCCTTT", 1), "MPF")
        self.assertEqual(translate_sequence("", 1), "")

    def test_translate_sequence_ignores_trailing_partial_codon(self):
        self.assertEqual(translate_sequence("ATGCCCTT", 1), "MP")
        self.assertEqual(translate_sequence("AT", 1), "")

    def test_phase_padding_codon_translates_to_x(self):
        self.assertEqual(translate_sequence("NNGATG", 1), "XM")


class TestStartAndStopCodons(unittest.TestCase):
    def tearDown(self):
        set_table1_start_codons('ncbi')

    def test_standard_table_ncbi_starts(self):
        for codon in ("ATG", "CTG", "TTG"):
            self.assertTrue(is_start_codon(codon, 1))
        for codon in ("GTG", "ATT", "TTC"):
            self.assertFalse(is_start_codon(codon, 1))

    def test_atg_only_start_set(self):
        set_table1_start_codons('atg-only')
        self.assertTrue(is_start_codon("ATG", 1))
        for codon in ("CTG", "TTG", "GTG"):
            self.assertFalse(is_start_codon(codon, 1))

    def test_start_set_switch_is_reversible(self):
        set_table1_start_codons('atg-only')
        set_table1_start_codons('ncbi')
        self.assertTrue(is_start_codon("CTG", 1))

    def test_start_set_does_not_touch_table_two(self):
        set_table1_start_codons('atg-only')
        self.assertTrue(is_start_codon("GTG", 2))
        self.assertTrue(is_start_codon("ATT", 2))

    def test_unknown_start_set_rejected(self):
        with self.assertRaises(ValueError):
            set_table1_start_codons('nonsense')

    def test_mitochondrial_starts(self):
        for codon in ("ATT", "ATC", "ATA", "ATG", "GTG"):
            self.assertTrue(is_start_codon(codon, 2))

    def test_terminators(self):
        for codon in ("TAA", "TAG", "TGA"):
            self.assertTrue(is_ter_codon(codon, 1))
        self.assertFalse(is_ter_codon("TGG", 1))
        self.assertTrue(is_ter_codon("AGA", 2))
        self.assertFalse(is_ter_codon("TGA", 2))


class TestSeqEdits(unittest.TestCase):
    def test_apply_substitution(self):
        self.assertEqual(apply_seq_edit("ABCDEF", 2, 4, "XYZ"), "AXYZEF")

    def test_apply_deletion(self):
        self.assertEqual(apply_seq_edit("ABCDEF", 2, 4, ""), "AEF")

    def test_apply_single_residue(self):
        self.assertEqual(apply_seq_edit("ABCDEF", 1, 1, "X"), "XBCDEF")

    def test_apply_insertion_before_first_base(self):
        self.assertEqual(apply_seq_edit("ABC", 1, 0, "XY"), "XYABC")

    def test_parse_seq_edit(self):
        self.assertEqual(parse_seq_edit("10 12 UGA"), (10, 12, "UGA"))
        self.assertEqual(parse_seq_edit("10 12"), (10, 12, ""))
        self.assertIsNone(parse_seq_edit("10"))
        self.assertIsNone(parse_seq_edit(""))

    def test_collect_orders_highest_coordinate_first(self):
        edits = collect_seq_edits(["5 5 A", "20 20 B", "1 1 C"])
        self.assertEqual([e[0] for e in edits], [20, 5, 1])

    def test_reverse_order_keeps_later_edits_correct(self):
        seq = "AAAAAAAAAA"
        for start, end, alt in collect_seq_edits(["2 3 XXXX", "8 8 Y"]):
            seq = apply_seq_edit(seq, start, end, alt)
        self.assertEqual(seq, "AXXXXAAAAYAA")

    def test_selenocysteine_edit_applied_to_peptide(self):
        # A UGA stop read through as selenocysteine.
        exons = [exon(1, 1, 12)]
        db = FakeDB(sequences={1: "ATGTGACCCTAA"},
                    translation_attribs={'_selenocysteine': ['2 2 U']})
        translation = {'translation_id': 1, 'seq_start': 1, 'seq_end': 12,
                       'start_exon_id': 1, 'end_exon_id': 1}
        self.assertEqual(
            translate_transcript(db, 1, exons, translation, 1), "MUP")


class TestCdnaCodingCoordinates(unittest.TestCase):
    def setUp(self):
        # Three 10 bp exons -> a 30 bp transcript.
        self.exons = [exon(1, 1, 10, rank=1),
                      exon(2, 21, 30, rank=2),
                      exon(3, 41, 50, rank=3)]

    def test_coding_start_in_first_exon(self):
        db = FakeDB()
        translation = {'translation_id': 1, 'seq_start': 4, 'seq_end': 6,
                       'start_exon_id': 1, 'end_exon_id': 3}
        self.assertEqual(compute_cdna_coding_start(db, 1, self.exons, translation), 4)

    def test_coding_start_in_later_exon_sums_preceding_lengths(self):
        db = FakeDB()
        translation = {'translation_id': 1, 'seq_start': 3, 'seq_end': 6,
                       'start_exon_id': 2, 'end_exon_id': 3}
        self.assertEqual(compute_cdna_coding_start(db, 1, self.exons, translation), 13)

    def test_coding_end_sums_preceding_lengths(self):
        db = FakeDB()
        translation = {'translation_id': 1, 'seq_start': 1, 'seq_end': 6,
                       'start_exon_id': 1, 'end_exon_id': 3}
        self.assertEqual(compute_cdna_coding_end(db, 1, self.exons, translation), 26)

    def test_transl_start_and_end_attributes_override(self):
        db = FakeDB(transcript_attribs={'_transl_start': ['7'],
                                        '_transl_end': ['19']})
        translation = {'translation_id': 1, 'seq_start': 1, 'seq_end': 6,
                       'start_exon_id': 1, 'end_exon_id': 3}
        self.assertEqual(compute_cdna_coding_start(db, 1, self.exons, translation), 7)
        self.assertEqual(compute_cdna_coding_end(db, 1, self.exons, translation), 19)

    def test_upstream_insertion_shifts_both_coordinates(self):
        # A 4 bp insertion replacing 1 bp at position 2 adds 3 bp upstream.
        db = FakeDB(transcript_attribs={'_rna_edit': ['2 2 AAAA']})
        translation = {'translation_id': 1, 'seq_start': 3, 'seq_end': 6,
                       'start_exon_id': 2, 'end_exon_id': 3}
        self.assertEqual(compute_cdna_coding_start(db, 1, self.exons, translation), 16)
        self.assertEqual(compute_cdna_coding_end(db, 1, self.exons, translation), 29)

    def test_downstream_edit_does_not_shift_coding_start(self):
        db = FakeDB(transcript_attribs={'_rna_edit': ['25 25 AAAA']})
        translation = {'translation_id': 1, 'seq_start': 3, 'seq_end': 6,
                       'start_exon_id': 2, 'end_exon_id': 3}
        self.assertEqual(compute_cdna_coding_start(db, 1, self.exons, translation), 13)


class TestSplicingAndTranslateableSeq(unittest.TestCase):
    def test_spliced_seq_concatenates_in_rank_order(self):
        exons = [exon(1, 1, 3, rank=1), exon(2, 11, 13, rank=2)]
        db = FakeDB(sequences={1: "AAA", 2: "CCC"})
        self.assertEqual(build_spliced_seq(db, 1, exons), "AAACCC")

    def test_reverse_strand_exon_is_reverse_complemented(self):
        exons = [exon(1, 1, 4, strand=-1, rank=1)]
        db = FakeDB(sequences={1: "ATCG"})
        self.assertEqual(build_spliced_seq(db, 1, exons), "CGAT")

    def test_rna_edit_applied_to_spliced_seq(self):
        exons = [exon(1, 1, 6, rank=1)]
        db = FakeDB(sequences={1: "AAAAAA"},
                    transcript_attribs={'_rna_edit': ['3 4 GG']})
        self.assertEqual(build_spliced_seq(db, 1, exons), "AAGGAA")

    def test_phase_padding_prepends_ns(self):
        exons = [exon(1, 1, 6, phase=2, rank=1)]
        db = FakeDB(sequences={1: "GATGCC"})
        translation = {'translation_id': 1, 'seq_start': 1, 'seq_end': 6,
                       'start_exon_id': 1, 'end_exon_id': 1}
        self.assertEqual(
            build_translateable_seq(db, 1, exons, translation), "NNGATGCC")

    def test_phase_zero_adds_no_padding(self):
        exons = [exon(1, 1, 6, phase=0, rank=1)]
        db = FakeDB(sequences={1: "ATGCCC"})
        translation = {'translation_id': 1, 'seq_start': 1, 'seq_end': 6,
                       'start_exon_id': 1, 'end_exon_id': 1}
        self.assertEqual(
            build_translateable_seq(db, 1, exons, translation), "ATGCCC")

    def test_no_translation_gives_empty_cds(self):
        self.assertEqual(build_translateable_seq(FakeDB(), 1, [], None), '')


class TestTranslateTranscript(unittest.TestCase):
    def _translation(self, seq_start=1, seq_end=None):
        return {'translation_id': 1, 'seq_start': seq_start,
                'seq_end': seq_end, 'start_exon_id': 1, 'end_exon_id': 1}

    def test_terminal_stop_removed(self):
        exons = [exon(1, 1, 12)]
        db = FakeDB(sequences={1: "ATGCCCTTTTAA"})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 12), 1), "MPF")

    def test_internal_stop_retained(self):
        exons = [exon(1, 1, 12)]
        db = FakeDB(sequences={1: "ATGTAACCCTTT"})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 12), 1), "M*PF")

    def test_incomplete_trailing_codon_trimmed(self):
        exons = [exon(1, 1, 11)]
        db = FakeDB(sequences={1: "ATGCCCTTTTA"})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 11), 1), "MPF")

    def test_initial_met_conversion_for_alternative_start(self):
        # GTG is a start codon under table 2, so V becomes M.
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "GTGCCCTTT"}, codon_table=2)
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "MPF")

    def test_initial_met_conversion_for_ctg_under_ncbi_starts(self):
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "CTGCCCTTT"})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "MPF")

    def test_no_initial_met_conversion_under_atg_only_starts(self):
        set_table1_start_codons('atg-only')
        self.addCleanup(set_table1_start_codons, 'ncbi')
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "CTGCCCTTT"})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "LPF")

    def test_no_conversion_for_non_start_codon(self):
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "GTGCCCTTT"})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "VPF")

    def test_cds_start_nf_does_not_suppress_initial_met_conversion(self):
        # Transcript::translate applies the rewrite regardless of cds_start_NF;
        # the GRCh38 release-116 dump confirms it for 224 records.
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "GTGCCCTTT"}, codon_table=2,
                    transcript_attribs={'cds_start_NF': ['1']})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "MPF")

    def test_initial_met_seq_edit_applied(self):
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "CTGCCCTTT"},
                    translation_attribs={'initial_met': ['1 1 M']})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "MPF")

    def test_amino_acid_substitution_applied(self):
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "ATGCCCTTT"},
                    translation_attribs={'amino_acid_sub': ['3 3 W']})
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "MPW")

    def test_mitochondrial_codon_table_used(self):
        # TGA is a stop under table 1 and W under table 2.
        exons = [exon(1, 1, 9)]
        db = FakeDB(sequences={1: "ATGTGACCC"}, codon_table=2)
        self.assertEqual(
            translate_transcript(db, 1, exons, self._translation(1, 9), 1), "MWP")

    def test_no_translation_returns_none(self):
        self.assertIsNone(translate_transcript(FakeDB(), 1, [], None, 1))

    def test_empty_cds_returns_none(self):
        self.assertIsNone(
            translate_transcript(FakeDB(), 1, [], self._translation(1, 9), 1,
                                 translateable=''))


class TestCodingRegionBounds(unittest.TestCase):
    def test_forward_strand(self):
        exons = [exon(1, 100, 199, rank=1), exon(2, 300, 399, rank=2)]
        translation = {'translation_id': 1, 'seq_start': 11, 'seq_end': 50,
                       'start_exon_id': 1, 'end_exon_id': 2}
        self.assertEqual(coding_region_bounds(exons, translation), (110, 349))

    def test_reverse_strand_returns_low_coordinate_first(self):
        exons = [exon(1, 300, 399, strand=-1, rank=1),
                 exon(2, 100, 199, strand=-1, rank=2)]
        translation = {'translation_id': 1, 'seq_start': 11, 'seq_end': 50,
                       'start_exon_id': 1, 'end_exon_id': 2}
        start, end = coding_region_bounds(exons, translation)
        self.assertEqual((start, end), (150, 389))
        self.assertLess(start, end)

    def test_missing_translation(self):
        self.assertEqual(coding_region_bounds([], None), (None, None))


class TestHeaderGeneration(unittest.TestCase):
    def setUp(self):
        self.exons = [exon(1, 84487828, 84511613, strand=-1, rank=1)]
        self.translation = {'translation_id': 1, 'stable_id': 'ENSP05220066032',
                            'version': 1, 'seq_start': 1, 'seq_end': 23786,
                            'start_exon_id': 1, 'end_exon_id': 1}
        self.transcript = {'transcript_id': 1, 'stable_id': 'ENST05220147416',
                           'version': 1, 'biotype': 'protein_coding',
                           'seq_region_strand': -1}
        self.gene = {'gene_id': 1, 'stable_id': 'ENSG05220038331', 'version': 1,
                     'biotype': 'protein_coding', 'description': 'SSX2IP',
                     'display_xref_id': 7}

    def test_full_header(self):
        db = FakeDB(xrefs={7: 'SSX2IP'})
        self.assertEqual(
            build_pep_header(db, self.transcript, self.gene, self.translation,
                             self.exons, 'T2T-CHM13v2.0', '1'),
            '>ENSP05220066032.1 pep T2T-CHM13v2.0:1:84487828:84511613:-1 '
            'gene:ENSG05220038331.1 transcript:ENST05220147416.1 '
            'gene_biotype:protein_coding transcript_biotype:protein_coding '
            'gene_symbol:SSX2IP description:"SSX2IP"'.lstrip('>')
        )

    def test_gene_symbol_omitted_without_display_xref(self):
        self.gene['display_xref_id'] = None
        db = FakeDB()
        self.assertNotIn('gene_symbol:', build_pep_header(
            db, self.transcript, self.gene, self.translation, self.exons,
            'T2T-CHM13v2.0', '1'))

    def test_description_omitted_when_null(self):
        self.gene['description'] = None
        db = FakeDB(xrefs={7: 'SSX2IP'})
        header = build_pep_header(db, self.transcript, self.gene,
                                  self.translation, self.exons,
                                  'T2T-CHM13v2.0', '1')
        self.assertNotIn('description:', header)
        self.assertTrue(header.endswith('gene_symbol:SSX2IP'))

    def test_description_is_quoted(self):
        db = FakeDB()
        header = build_pep_header(db, self.transcript, self.gene,
                                  self.translation, self.exons, 'GRCh38', '1')
        self.assertIn('description:"SSX2IP"', header)


class TestSliceClassification(unittest.TestCase):
    def _slice(self, name, length, seq_region_id=None):
        return {'name': name, 'length': length, 'start': 1, 'end': length,
                'seq_region_id': seq_region_id if seq_region_id is not None else length}

    def test_sorted_by_descending_length_and_split(self):
        slices = [self._slice('scaffold', 500, 1),
                  self._slice('1', 1000, 2),
                  self._slice('alt', 800, 3),
                  self._slice('2', 900, 4)]
        chromosomes, non_chromosomal, non_reference = classify_slices(
            slices,
            is_reference=lambda s: s['seq_region_id'] != 3,
            is_chromosome=lambda s: s['name'] in ('1', '2'),
        )
        self.assertEqual([s['name'] for s in chromosomes], ['1', '2'])
        self.assertEqual([s['name'] for s in non_chromosomal], ['scaffold'])
        self.assertEqual([s['name'] for s in non_reference], ['alt'])

    def test_equal_lengths_keep_input_order(self):
        slices = [self._slice('a', 100, 1), self._slice('b', 100, 2),
                  self._slice('c', 100, 3)]
        _, non_chromosomal, _ = classify_slices(
            slices, is_reference=lambda s: True, is_chromosome=lambda s: False)
        self.assertEqual([s['name'] for s in non_chromosomal], ['a', 'b', 'c'])

    def test_y_par_filter_skipped_for_assembly_specific_species(self):
        slices = [self._slice('Y', 1000)]
        self.assertEqual(
            filter_grch38_y_par(slices, 'homo_sapiens_gca009914755v4',
                                'T2T-CHM13v2.0'),
            slices)

    def test_y_par_filter_drops_short_y_regions_on_grch38(self):
        slices = [self._slice('Y', 57227415), self._slice('Y', 1000),
                  self._slice('1', 248956422)]
        kept = filter_grch38_y_par(slices, 'homo_sapiens', 'GRCh38')
        self.assertEqual([(s['name'], s['length']) for s in kept],
                         [('Y', 57227415), ('1', 248956422)])

    def test_y_par_filter_raises_for_other_human_assemblies(self):
        with self.assertRaises(CoreDbError):
            filter_grch38_y_par([self._slice('Y', 1000)], 'homo_sapiens', 'GRCh37')


if __name__ == '__main__':
    unittest.main(verbosity=2)

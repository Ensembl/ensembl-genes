"""Parser contract for pairwise comparison: coordinates, formats, parentage, IDs, seqnames."""

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.parsers.annotation import (
    detect_annotation_format,
    parse_annotation,
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.annotation_normalise import COMPARISON_COLUMNS
from ensembl.genes.annotation_qc.parsers.regions import (
    Region,
    load_regions_file,
    parse_region,
)
from ensembl.genes.annotation_qc.parsers.seqnames import (
    apply_seqname_mapping,
    build_seqname_mapping,
    check_seqnames_against_genome,
    load_seqname_map,
)
from ensembl.genes.annotation_qc.parsers.sequence import parse_sequence_lengths

from .conftest import gene_records

ENSEMBL_GFF3 = [
    ("c1", "region", 1, 1000, ".", "ID=region:c1"),
    ("c1", "gene", 10, 100, "+", "ID=gene:G1;biotype=protein_coding;gene_id=G1"),
    (
        "c1",
        "mRNA",
        10,
        100,
        "+",
        "ID=transcript:T1;Parent=gene:G1;biotype=protein_coding;tag=Ensembl_canonical",
    ),
    (
        "c1",
        "mRNA",
        10,
        100,
        "+",
        "ID=transcript:T2;Parent=gene:G1;biotype=nonsense_mediated_decay",
    ),
    ("c1", "exon", 10, 40, "+", "Parent=transcript:T1,transcript:T2"),
    ("c1", "exon", 60, 100, "+", "Parent=transcript:T1"),
    ("c1", "exon", 70, 100, "+", "Parent=transcript:T2"),
    ("c1", "CDS", 12, 40, "+", "ID=CDS:P1;Parent=transcript:T1"),
    ("c1", "CDS", 60, 90, "+", "ID=CDS:P1;Parent=transcript:T1"),
    ("c1", "five_prime_UTR", 10, 11, "+", "Parent=transcript:T1"),
]

GTF = (
    'c1\tt\tgene\t10\t100\t.\t-\t.\tgene_id "G1"; gene_biotype "protein_coding";\n'
    'c1\tt\ttranscript\t10\t100\t.\t-\t.\tgene_id "G1"; transcript_id "T1"; tag "basic"; tag "Ensembl_canonical";\n'
    'c1\tt\texon\t10\t100\t.\t-\t.\tgene_id "G1"; transcript_id "T1";\n'
    'c1\tt\tCDS\t20\t90\t.\t-\t0\tgene_id "G1"; transcript_id "T1";\n'
    'c1\tt\texon\t300\t400\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n'
)

TIBERIUS = (
    "1\tTiberius\tgene\t100\t500\t.\t+\t.\tg1\n"
    "1\tTiberius\ttranscript\t100\t500\t.\t+\t.\tg1.t1\n"
    '1\tTiberius\texon\t100\t200\t.\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";\n'
    '1\tTiberius\tCDS\t100\t200\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    '1\tTiberius\texon\t300\t500\t.\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";\n'
    "2\tTiberius\tgene\t100\t500\t.\t-\t.\tg1\n"
    "2\tTiberius\ttranscript\t100\t500\t.\t-\t.\tg1.t1\n"
    '2\tTiberius\texon\t100\t500\t.\t-\t.\ttranscript_id "g1.t1"; gene_id "g1";\n'
    '2\tTiberius\tCDS\t100\t500\t.\t-\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
)


def _rows(df, feature):
    return df[df["Feature"] == feature]


class TestCoordinatesAndSchema:
    def test_gff3_start_converted_once_to_zero_based_half_open(self, write_gff3):
        df = parse_annotation_for_comparison(write_gff3("a.gff3", ENSEMBL_GFF3))
        gene = _rows(df, "gene").iloc[0]
        assert (gene["Start"], gene["End"]) == (9, 100)  # GFF 10..100 is 91 bp
        assert gene["End"] - gene["Start"] == 91
        assert list(df.columns) == COMPARISON_COLUMNS

    def test_gtf_uses_same_convention(self, write_text):
        df = parse_annotation_for_comparison(write_text("a.gtf", GTF))
        cds = _rows(df, "CDS").iloc[0]
        assert (cds["Start"], cds["End"]) == (19, 90)

    def test_existing_parse_annotation_unchanged(self, write_gff3):
        data = parse_annotation(write_gff3("a.gff3", ENSEMBL_GFF3))
        assert data.__class__.__name__ == "PyRanges"
        assert (
            "Parent" in data.columns
            and data["Parent"].str.startswith("transcript:").any()
        )

    @pytest.mark.parametrize("compress", [False, True])
    def test_compressed_input_matches_plain(self, write_gff3, compress):
        plain = parse_annotation_for_comparison(write_gff3("p.gff3", ENSEMBL_GFF3))
        other = parse_annotation_for_comparison(
            write_gff3(
                "x.gff3.gz" if compress else "x.gff3", ENSEMBL_GFF3, compress=compress
            )
        )
        pd.testing.assert_frame_equal(plain, other)

    def test_format_detected_from_content_not_extension(self, write_text, write_gff3):
        assert (
            detect_annotation_format(
                write_text(
                    "really_gff3.gtf", open(write_gff3("a.gff3", ENSEMBL_GFF3)).read()
                )
            )
            == "gff3"
        )
        assert detect_annotation_format(write_text("really_gtf.gff3", GTF)) == "gtf"
        assert (
            detect_annotation_format(write_text("tib.gtf.gz", TIBERIUS, compress=True))
            == "tiberius"
        )


class TestParentage:
    def test_ensembl_prefixes_stripped_and_multi_parent_exon_split(self, write_gff3):
        df = parse_annotation_for_comparison(write_gff3("a.gff3", ENSEMBL_GFF3))
        exons = _rows(df, "exon")
        assert sorted(exons["transcript_id"]) == ["T1", "T1", "T2", "T2"]
        shared = exons[exons["Start"] == 9]
        assert sorted(shared["Parent"]) == ["T1", "T2"]
        assert set(exons["gene_id"]) == {"G1"}
        assert df.attrs["parse_diagnostics"]["multi_parent_children"] == 1
        assert (
            "five_prime_UTR" in df.attrs["parse_diagnostics"]["dropped_feature_types"]
        )
        assert "region" in df.attrs["parse_diagnostics"]["dropped_feature_types"]

    def test_biotypes_and_tags_propagated(self, write_gff3):
        df = parse_annotation_for_comparison(write_gff3("a.gff3", ENSEMBL_GFF3))
        tx = _rows(df, "transcript").set_index("transcript_id")
        assert tx.loc["T1", "tags"] == "Ensembl_canonical"
        assert tx.loc["T2", "transcript_biotype"] == "nonsense_mediated_decay"
        assert set(df.loc[df["Feature"] != "gene", "gene_biotype"]) == {
            "protein_coding"
        }
        cds = _rows(df, "CDS")
        assert set(cds["transcript_biotype"]) == {"protein_coding"}

    def test_gtf_repeated_tags_kept_and_missing_parents_synthesised(self, write_text):
        df = parse_annotation_for_comparison(write_text("a.gtf", GTF))
        tx = _rows(df, "transcript").set_index("transcript_id")
        assert (
            "Ensembl_canonical" in tx.loc["T1", "tags"]
            and "basic" in tx.loc["T1", "tags"]
        )
        g2 = _rows(df, "gene").set_index("gene_id").loc["G2"]
        assert (g2["Start"], g2["End"], g2["source_feature"]) == (
            299,
            400,
            "inferred_gene",
        )
        assert df.attrs["parse_diagnostics"]["synthesised_genes"] == 1
        assert df.attrs["parse_diagnostics"]["synthesised_transcripts"] == 1

    def test_orphan_exon_rejected(self, write_gff3):
        records = gene_records("c1", "G1", "+", {"T1": [(10, 50)]})
        records.append(("c1", "exon", 60, 70, "+", "Parent=NOPE"))
        with pytest.raises(ValueError, match="Parent is not a transcript"):
            parse_annotation_for_comparison(write_gff3("a.gff3", records))

    def test_transcript_without_gene_rejected(self, write_gff3):
        records = [
            ("c1", "mRNA", 10, 50, "+", "ID=T1;Parent=MISSING"),
            ("c1", "exon", 10, 50, "+", "Parent=T1"),
        ]
        with pytest.raises(ValueError, match="not a gene record"):
            parse_annotation_for_comparison(write_gff3("a.gff3", records))

    def test_duplicate_gene_id_on_one_sequence_rejected(self, write_gff3):
        records = gene_records("c1", "G1", "+", {"T1": [(10, 50)]})
        records += gene_records("c1", "G1", "+", {"T2": [(100, 150)]})
        with pytest.raises(ValueError, match="Duplicate gene"):
            parse_annotation_for_comparison(write_gff3("a.gff3", records))

    def test_gtf_record_without_transcript_id_rejected(self, write_text):
        bad = GTF + 'c1\tt\texon\t500\t600\t.\t+\t.\tgene_id "G3";\n'
        with pytest.raises(ValueError, match="without gene_id/transcript_id"):
            parse_annotation_for_comparison(write_text("a.gtf", bad))


class TestDuplicateIds:
    def test_tiberius_hybrid_ids_recovered_and_namespaced_per_sequence(
        self, write_text
    ):
        df = parse_annotation_for_comparison(write_text("t.gtf", TIBERIUS))
        assert df.attrs["parse_diagnostics"]["format"] == "tiberius"
        genes = _rows(df, "gene").set_index("Chromosome")
        assert (
            genes.loc["1", "gene_id"] == "1:g1" and genes.loc["2", "gene_id"] == "2:g1"
        )
        assert set(genes["original_gene_id"]) == {"g1"}
        exons = _rows(df, "exon")
        assert set(exons.loc[exons["Chromosome"] == "2", "transcript_id"]) == {
            "2:g1.t1"
        }
        assert set(exons["original_transcript_id"]) == {"g1.t1"}
        assert df.attrs["parse_diagnostics"]["namespaced_gene_ids"] == 1

    def test_unique_ids_untouched(self, write_gff3):
        df = parse_annotation_for_comparison(write_gff3("a.gff3", ENSEMBL_GFF3))
        assert (df["gene_id"] == df["original_gene_id"]).all()


class TestSeqnames:
    @pytest.mark.parametrize(
        "text",
        [
            "from_seqname\tto_seqname\nchrA\t1\n# comment\nchrB\t2\n",
            "chrA,1\nchrB,2\n",
            "chrA   1\n\nchrB   2\n",
        ],
    )
    def test_map_formats(self, write_text, text):
        assert load_seqname_map(write_text("m.tsv", text)) == {"chrA": "1", "chrB": "2"}

    def test_conflicting_map_entry_rejected(self, write_text):
        with pytest.raises(ValueError, match="mapped to both"):
            load_seqname_map(write_text("m.tsv", "chrA\t1\nchrA\t2\n"))

    def test_custom_map_overrides_assembly_report(self, write_text):
        report = write_text(
            "r.txt",
            "# Sequence-Name\tRole\tMolecule\tType\tGenBank\n"
            "chr1\tassembled\t1\tChromosome\tCM000001.1\n"
            "chr2\tassembled\t2\tChromosome\tCM000002.1\n",
        )
        custom = write_text("m.tsv", "CM000002.1\tII\n")
        assert build_seqname_mapping(custom, report) == {
            "CM000001.1": "1",
            "CM000002.1": "II",
        }

    def test_mapping_applies_and_rejects_merges(self, write_gff3):
        records = gene_records("chrA", "G1", "+", {"T1": [(10, 50)]})
        records += gene_records("chrB", "G2", "+", {"T2": [(10, 50)]})
        df = parse_annotation_for_comparison(write_gff3("a.gff3", records))
        mapped, report = apply_seqname_mapping(df, {"chrA": "1"})
        assert set(mapped["Chromosome"]) == {"1", "chrB"}
        assert report["unmapped"] == ["chrB"]
        with pytest.raises(ValueError, match="merges distinct sequences"):
            apply_seqname_mapping(df, {"chrA": "1", "chrB": "1"})

    def test_genome_check_reports_missing_and_out_of_bounds(self, write_gff3):
        records = gene_records("1", "G1", "+", {"T1": [(10, 50)]})
        records += gene_records("2", "G2", "+", {"T2": [(10, 150)]})
        records += gene_records("MIT", "G3", "+", {"T3": [(10, 50)]})
        df = parse_annotation_for_comparison(write_gff3("a.gff3", records))
        result = check_seqnames_against_genome(df, {"1": 50, "2": 100})
        assert result["missing_seqnames"] == {"MIT": 1}
        assert result["out_of_bounds_seqnames"] == {"2": 1}
        assert result["out_of_bounds_features"] == 4  # gene, mRNA, exon, CDS

    @pytest.mark.parametrize("compress", [False, True])
    def test_sequence_lengths(self, write_text, compress):
        fasta = write_text(
            "g.fa", ">1 description\nACGT\nAC\n>2\nA\n", compress=compress
        )
        assert parse_sequence_lengths(fasta) == {"1": 6, "2": 1}


class TestRegions:
    def test_region_strings_are_one_based_inclusive(self):
        assert parse_region("chr1:100-200") == Region("chr1", 99, 200)
        assert parse_region("chr1:1,000-2,000") == Region("chr1", 999, 2000)
        assert parse_region("chr1:5") == Region("chr1", 4, 5)
        assert parse_region("chrUn:KI270302.1") == Region("chrUn:KI270302.1")
        assert str(parse_region("chr1:100-200")) == "chr1:100-200"

    def test_regions_file(self, write_text):
        path = write_text("r.txt", "# regions\n1:10-20\n\n2\n")
        assert load_regions_file(path) == [Region("1", 9, 20), Region("2")]

# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""GFF3 hierarchy compatibility and GTF regression checks."""
import tempfile
import unittest
from pathlib import Path

from ensembl.genes.annotation_qc.parsers.high_confidence import (
    parse_high_confidence_annotation,
)
from ensembl.genes.annotation_qc.metrics.high_confidence import (
    _models_from_dataframe,
    select_transcripts,
)


def read_gtf(path):
    """Adapt the public DataFrame parser for hierarchy assertions."""
    frame = parse_high_confidence_annotation(path)
    assert {"Chromosome", "Start", "End", "Feature", "gene_id", "transcript_id"} <= set(
        frame.columns
    )
    return _models_from_dataframe(frame)


class Gff3InputTests(unittest.TestCase):
    """Check predictor GFF3 hierarchies and malformed parent links."""

    def test_supplied_helixer_file(self):
        """Check supplied helixer file."""
        path = Path(__file__).parent / "data" / "helixer.gff"
        if not path.exists():
            self.skipTest("User-supplied Helixer example is not available")
        genes = read_gtf(path)
        self.assertEqual(len(genes), 3)
        self.assertEqual(sum(len(g.transcripts) for g in genes.values()), 3)
        self.assertEqual(
            sum(
                len(t.of_type("exon"))
                for g in genes.values()
                for t in g.transcripts.values()
            ),
            10,
        )
        self.assertEqual(
            sum(
                len(t.of_type("CDS"))
                for g in genes.values()
                for t in g.transcripts.values()
            ),
            10,
        )
        self.assertEqual(len(select_transcripts(genes, {}, "auto")), 3)
        for g in genes.values():
            for tx in g.transcripts.values():
                self.assertTrue(tx.of_type("five_prime_utr"))
                self.assertTrue(tx.of_type("three_prime_utr"))

    def read(self, lines):
        """Parse temporary annotation rows for a hierarchy test."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "annotation.gff"
            path.write_text("\n".join(lines) + "\n")
            return read_gtf(path)

    def row(self, kind, attrs, start=1, end=90, phase="."):
        """Construct one predictor GFF3 feature row."""
        return f"chr1\tHelixer\t{kind}\t{start}\t{end}\t.\t+\t{phase}\t{attrs}"

    def test_helixer_hierarchy_out_of_order(self):
        """Check helixer hierarchy out of order."""
        genes = self.read(
            [
                "##gff-version 3",
                self.row("CDS", "Parent=g1.t1", 11, 70, "0"),
                self.row("exon", "ID=e1;Parent=g1.t1"),
                self.row("mRNA", "ID=g1.t1;Parent=g1"),
                self.row("gene", "ID=g1"),
                self.row("five_prime_UTR", "Parent=g1.t1", 1, 10),
                self.row("three_prime_UTR", "Parent=g1.t1", 71, 90),
            ]
        )
        tx = genes["g1"].transcripts["g1.t1"]
        self.assertEqual(tx.gene_id, "g1")
        self.assertEqual(len(tx.of_type("transcript")), 1)
        self.assertEqual(len(tx.of_type("exon")), 1)
        self.assertEqual(tx.of_type("CDS")[0].phase, 0)
        self.assertEqual(tx.of_type("five_prime_utr")[0].end, 10)

    def test_shared_exon_and_escaped_ids(self):
        """Check shared exon and escaped ids."""
        genes = self.read(
            [
                self.row("gene", "ID=g%3B1"),
                self.row("mRNA", "ID=t%2C1;Parent=g%3B1"),
                self.row("mRNA", "ID=t2;Parent=g%3B1"),
                self.row("exon", "Parent=t%2C1,t2"),
            ]
        )
        self.assertEqual(set(genes["g;1"].transcripts), {"t,1", "t2"})
        for tx in genes["g;1"].transcripts.values():
            self.assertEqual(len(tx.of_type("exon")), 1)

    def test_explicit_ids_and_fasta_tail(self):
        """Check explicit ids and fasta tail."""
        genes = self.read(
            [
                self.row("gene", "ID=gene:g1;gene_id=g1"),
                self.row("mRNA", "ID=transcript:t1;Parent=gene:g1;transcript_id=t1"),
                self.row("exon", "Parent=transcript:t1"),
                "##FASTA",
                ">chr1",
                "ATG",
            ]
        )
        self.assertIn("t1", genes["g1"].transcripts)

    def test_orphan_parent_rejected(self):
        """Check orphan parent rejected."""
        with self.assertRaisesRegex(
            ValueError, "annotation.gff:1:.*does not reference"
        ):
            self.read([self.row("CDS", "Parent=missing", 1, 90, "0")])

    def test_duplicate_parent_and_conflicting_ids(self):
        """Check duplicate parent and conflicting ids."""
        with self.assertRaisesRegex(ValueError, "duplicate GFF3 parent ID"):
            self.read([self.row("gene", "ID=g1"), self.row("gene", "ID=g1")])
        with self.assertRaisesRegex(ValueError, "Conflicting gene_id"):
            self.read(
                [
                    self.row("gene", "ID=g1"),
                    self.row("mRNA", "ID=t1;Parent=g1;gene_id=wrong"),
                ]
            )

    def test_gtf_and_bare_tiberius_ids_unchanged(self):
        """Check gtf and bare tiberius ids unchanged."""
        for gene_attrs, tx_attrs in [
            ('gene_id "g1";', 'gene_id "g1"; transcript_id "t1";'),
            ("g1", "t1"),
        ]:
            genes = self.read(
                [
                    self.row("gene", gene_attrs),
                    self.row("transcript", tx_attrs),
                    self.row("exon", 'gene_id "g1"; transcript_id "t1";'),
                    self.row("CDS", 'gene_id "g1"; transcript_id "t1";', 1, 90, "0"),
                ]
            )
            self.assertEqual(len(genes["g1"].transcripts["t1"].of_type("CDS")), 1)


if __name__ == "__main__":
    unittest.main()

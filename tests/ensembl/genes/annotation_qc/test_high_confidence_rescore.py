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
"""Regression checks for evidence rescoring and output provenance."""
import json
from pathlib import Path
import tempfile
import unittest

from ensembl.genes.annotation_qc.config.high_confidence import FLAGS
from ensembl.genes.annotation_qc.runners.rescore_high_confidence import rescore
from ensembl.genes.annotation_qc.reports.high_confidence_summary import write_tsv
from ensembl.genes.annotation_qc.parsers.high_confidence import parse_tsv_frame


def read_tsv(path):
    """Read summary rows while retaining the original field order."""
    frame = parse_tsv_frame(path)
    return list(frame.columns), frame.to_dict("records")


class RescoreEvidenceTests(unittest.TestCase):
    """Check policy migration without changing original evidence."""

    def test_rescore_updates_table_manifest_gff_and_preserves_source(self):
        """Check rescore updates table manifest gff and preserves source."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "species.evidence.tsv"
            base = dict.fromkeys(FLAGS, "1")
            base.update(
                transcript_five_prime_utr="0",
                transcript_three_prime_utr="0",
                homology_status="hit",
                error="8",
                score="70",
                confidence="high",
                confidence_rule="utr_ignored_error_le_10",
            )
            rows = [
                dict(
                    base,
                    gene_id="g1",
                    transcript_id="t1",
                    transcript_intron_support="0",
                ),
                dict(
                    base,
                    gene_id="g2",
                    transcript_id="t2",
                    transcript_intron_support="0",
                    error="2",
                    score="90",
                ),
                dict(base, gene_id="g3", transcript_id="t3"),
            ]
            write_tsv(source, rows, list(rows[0]))
            original = source.read_bytes()
            (root / "species.manifest.json").write_text(
                json.dumps(
                    {
                        "processed_transcript_count": 3,
                        "options": {"no_utr": True},
                        "inputs": {"genome": {"sha256": "original-genome-hash"}},
                    }
                )
            )
            gff = root / "species.annotations.gff3"
            gff.write_text(
                "##gff-version 3\n"
                + "".join(
                    f"chr1\tx\tgene\t1\t90\t.\t+\t.\tID=gene:{g};gene_id={g};gene_quality=high\n"
                    for g in ["g1", "g2", "g3"]
                )
            )
            introns = root / "species.introns.tsv"
            introns.write_text("transcript_id\n")
            out = root / "rescored"
            manifest = rescore(source, out, True)
            self.assertEqual(source.read_bytes(), original)
            _, updated = read_tsv(out / source.name)
            self.assertEqual(
                [r["confidence"] for r in updated], ["low", "medium", "high"]
            )
            self.assertEqual(manifest["confidence_policy_version"], 2)
            self.assertTrue(manifest["high_confidence_requires_intron_support"])
            self.assertEqual(
                manifest["confidence_counts"], {"high": 1, "medium": 1, "low": 1}
            )
            self.assertEqual(
                manifest["inputs"]["genome"]["sha256"], "original-genome-hash"
            )
            self.assertEqual(manifest["rescoring"]["class_changed_count"], 2)
            self.assertIn("gene_id=g1;gene_quality=low", (out / gff.name).read_text())
            self.assertIn(
                "gene_id=g2;gene_quality=medium", (out / gff.name).read_text()
            )
            self.assertNotIn("gene_quality=low", gff.read_text())
            self.assertTrue((out / introns.name).is_symlink())
            self.assertEqual((out / introns.name).resolve(), introns.resolve())
            with self.assertRaisesRegex(ValueError, "already exist"):
                rescore(source, out, True)


if __name__ == "__main__":
    unittest.main()

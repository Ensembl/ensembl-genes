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
"""Boundary and exhaustive checks for strict RNA-supported high confidence under no_utr."""
import itertools
import unittest
from ensembl.genes.annotation_qc.config.high_confidence import FLAGS
from ensembl.genes.annotation_qc.metrics.high_confidence import confidence_rating


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-return-statements
def legacy(e, h):
    """Reference the original UTR-aware Perl confidence branches."""
    complete, five, three, canonical, regular, support = (e[f] for f in FLAGS)
    if not complete or h is None:
        return "low", "incomplete_cds_or_no_homology"
    error, score = h["error"], h["score"]
    if regular and canonical and five and three and error <= 10:
        return "high", "both_utrs_error_le_10"
    if regular and canonical and error <= 5 and score >= 80:
        return "high", "error_le_5_score_ge_80"
    if regular and canonical and error <= 10 and support:
        return "high", "supported_error_le_10"
    if (regular or canonical) and error <= 10 and support:
        return "medium", "supported_one_intron_check"
    if (regular or canonical) and error <= 5 and score >= 80:
        return "medium", "strong_homology_one_intron_check"
    return "low", "no_high_or_medium_rule"


class ConfidencePolicyTests(unittest.TestCase):
    """Check legacy parity and the support-required no-UTR policy."""

    def evidence(self, **updates):
        """Create complete evidence with specified flag overrides."""
        e = dict.fromkeys(FLAGS, True)
        e.update(updates)
        return e

    def test_all_flag_combinations_and_threshold_boundaries(self):
        """Check all flag combinations and threshold boundaries."""
        for bits in itertools.product([False, True], repeat=5):
            for support in [False, True, None]:
                e = dict(zip(FLAGS, (*bits, support)))
                for h in [None] + [
                    {"error": error, "score": score}
                    for error in [0, 5, 5.01, 8, 10, 10.01]
                    for score in [79.99, 80, 100]
                ]:
                    self.assertEqual(confidence_rating(e, h), legacy(e, h))
                    actual = confidence_rating(e, h, no_utr=True)[0]
                    complete, _, _, canonical, regular, _ = (e[f] for f in FLAGS)
                    expected = "low"
                    if complete and h is not None:
                        strong = h["error"] <= 5 and h["score"] >= 80
                        if (
                            regular
                            and canonical
                            and support is True
                            and h["error"] <= 10
                        ):
                            expected = "high"
                        elif (regular or canonical) and (
                            (support is True and h["error"] <= 10) or strong
                        ):
                            expected = "medium"
                    self.assertEqual(actual, expected)
                    if actual == "high":
                        self.assertIs(support, True)
                    # UTR values have no effect under no_utr.
                    altered = dict(
                        e,
                        transcript_five_prime_utr=not e["transcript_five_prime_utr"],
                        transcript_three_prime_utr=not e["transcript_three_prime_utr"],
                    )
                    self.assertEqual(
                        confidence_rating(e, h, no_utr=True),
                        confidence_rating(altered, h, no_utr=True),
                    )

    def test_strong_homology_without_support_is_medium(self):
        """Check strong homology without support is medium."""
        for support in [False, None]:
            e = self.evidence(transcript_intron_support=support)
            self.assertEqual(
                confidence_rating(e, {"error": 5, "score": 80}, no_utr=True),
                ("medium", "strong_homology_one_intron_check"),
            )

    def test_support_route_remains_high_without_utrs(self):
        """Check support route remains high without utrs."""
        e = self.evidence(
            transcript_five_prime_utr=False, transcript_three_prime_utr=False
        )
        self.assertEqual(
            confidence_rating(e, {"error": 10, "score": 70}, no_utr=True),
            ("high", "supported_error_le_10"),
        )

    def test_missing_support_and_weaker_homology_is_low(self):
        """Check missing support and weaker homology is low."""
        e = self.evidence(transcript_intron_support=None)
        self.assertEqual(
            confidence_rating(e, {"error": 8, "score": 85}, no_utr=True),
            ("low", "no_high_or_medium_rule"),
        )


if __name__ == "__main__":
    unittest.main()

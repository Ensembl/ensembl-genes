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
"""Pure evidence summaries distinguish unknowns and preserve supplied labels."""
from copy import deepcopy

import pandas as pd

from ensembl.genes.annotation_qc.config.high_confidence import FLAGS
from ensembl.genes.annotation_qc.metrics.high_confidence_summary import (
    summarize_confidence,
)
from ensembl.genes.annotation_qc.metrics.high_confidence import (
    select_representative_transcripts,
)
from ensembl.genes.annotation_qc.parsers.high_confidence import (
    parse_high_confidence_annotation,
)


def test_summary_unknown_support_and_missing_homology():
    """Unknown flags are excluded from measured-only rates but included in all-gene denominators."""
    base = dict.fromkeys(FLAGS, "1")
    base.update(transcript_five_prime_utr="0", transcript_three_prime_utr="0")
    rows = [
        dict(
            base,
            gene_id="g1",
            transcript_id="t1",
            transcript_intron_support=".",
            cds_intron_count="1",
            intron_count="1",
            error="8",
            score="85",
            confidence="low",
            confidence_rule="no_high_or_medium_rule",
            homology_status="hit",
        ),
        dict(
            base,
            gene_id="g2",
            transcript_id="t2",
            cds_intron_count="0",
            intron_count="0",
            error=".",
            score=".",
            confidence="low",
            confidence_rule="incomplete_cds_or_no_homology",
            homology_status="no_hit",
        ),
    ]
    frame = pd.DataFrame(rows, dtype=object)
    before = frame.copy(deep=True)
    result = summarize_confidence(
        frame,
        "M",
        manifest={
            "processed_transcript_count": 2,
            "confidence_policy_version": 2,
            "options": {"no_utr": True},
        },
    )
    pd.testing.assert_frame_equal(frame, before)
    states = {
        r["state"]: r
        for r in result["tables"]["flag_distribution"][0]
        if r["group"] == "all" and r["flag"] == "transcript_intron_support"
    }
    assert states["unknown"]["count"] == 1
    assert states["pass"]["percent"] == 50
    assert states["pass"]["pass_percent_measured"] == 100
    assert result["overview"]["low"] == 2
    assert result["overview"]["utr_sensitive_transcripts"] == 0


def test_selection_does_not_mutate_annotation_metadata(tmp_path):
    """The private reconstruction bridge never mutates the public DataFrame's metadata."""
    annotation = tmp_path / "model.gtf"
    annotation.write_text(
        'chr1\tx\tCDS\t1\t9\t.\t+\t0\tgene_id "g"; transcript_id "t";\n'
    )
    frame = parse_high_confidence_annotation(annotation)
    attrs = deepcopy(frame.attrs)
    assignments = pd.DataFrame(columns=["gene_id", "transcript_id"])
    select_representative_transcripts(frame, assignments)
    assert frame.attrs == attrs

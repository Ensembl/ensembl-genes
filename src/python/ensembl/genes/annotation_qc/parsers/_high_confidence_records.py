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
"""Private, lossless reconstruction records; public annotation boundaries use DataFrames."""

from __future__ import annotations

from dataclasses import dataclass, field
from ensembl.genes.annotation_qc.config.high_confidence import COMPLEMENT


@dataclass(frozen=True)
class Feature:
    """Private annotation row using inclusive genomic coordinates."""

    seqname: str
    kind: str
    start: int
    end: int
    strand: str
    phase: int | None
    attrs: dict[str, list[str]]


@dataclass
# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-instance-attributes
class Transcript:
    """Private feature collection preserving transcript metadata."""

    id: str
    gene_id: str
    seqname: str
    strand: str
    features: list[Feature] = field(default_factory=list)
    attrs: dict[str, list[str]] = field(default_factory=dict)
    parent_bounds: tuple[int, int] | None = None
    normalization: list = field(default_factory=list)

    def of_type(self, kind: str) -> list[Feature]:
        """Return transcript features in transcriptional order."""
        return sorted(
            (f for f in self.features if f.kind == kind),
            key=lambda f: f.start,
            reverse=self.strand == "-",
        )


@dataclass
# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-instance-attributes
class Gene:
    """Private gene bounds and associated transcripts."""

    id: str
    seqname: str
    strand: str
    start: int
    end: int
    attrs: dict[str, list[str]] = field(default_factory=dict)
    transcripts: dict[str, Transcript] = field(default_factory=dict)
    declared: tuple[int, int] | None = None
    normalization: list = field(default_factory=list)


class SequenceRegion:
    """In-memory genomic region used by pure reconstruction metrics (1-based API)."""

    def __init__(self, seqname, start, sequence, contig_length):
        self.seqname = seqname
        self.start = start
        self.sequence = sequence
        self.index = {seqname: (contig_length, 0, 0, 0)}

    def contains(self, seqname, start, end):
        """Check inclusive bounds against the indexed contig length."""
        return seqname == self.seqname and 1 <= start <= end <= self.index[seqname][0]

    def fetch(self, seqname, start, end, strand="+"):
        """Return an inclusive interval in transcriptional orientation."""
        if (
            seqname != self.seqname
            or start < self.start
            or end >= self.start + len(self.sequence)
        ):
            raise ValueError(f"Sequence was not prefetched: {seqname}:{start}-{end}")
        result = self.sequence[start - self.start : end - self.start + 1]
        return result.translate(COMPLEMENT)[::-1] if strand == "-" else result

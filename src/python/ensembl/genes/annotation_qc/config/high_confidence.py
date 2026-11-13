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
"""Constants and output schemas for high-confidence QC."""

from __future__ import annotations


COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

DEFAULT_MITOCHONDRIAL_CONTIGS = ("MT", "chrM", "chrMT", "M")

CANONICAL_SPLICES = {"GT-AG", "GC-AG", "AT-AC"}

FLAGS = [
    "transcript_cds_complete",
    "transcript_five_prime_utr",
    "transcript_three_prime_utr",
    "transcript_canonical_splicing",
    "transcript_regular_intron_size",
    "transcript_intron_support",
]

HOMOLOGY_FIELDS = [
    "hit_id",
    "error",
    "perc_id",
    "score",
    "query_coverage",
    "target_coverage",
]

CONFIDENCE_POLICY_VERSION = 2

MIN_INTRON_LENGTH = 50
MAX_INTRON_LENGTH = 100000
HOMOLOGY_ERROR_MAX = 10
STRONG_HOMOLOGY_ERROR_MAX = 5
STRONG_HOMOLOGY_SCORE_MIN = 80
DIAGNOSTIC_FIELDS = [
    "cds_length",
    "cds_length_multiple_of_three",
    "internal_stop_count",
    "ambiguous_amino_acids",
    "translation_status",
    "selenocysteine_count",
    "translation_table",
    "translation_table_source",
    "initiator_methionine_adjusted",
    "trimmed_terminal_bases",
    "exon_structure_source",
    "cds_phase_source",
    "inferred_cds_phase_count",
    "stop_codon_source",
    "utr_annotation_source",
    "intron_count",
    "cds_intron_count",
    "supported_cds_intron_count",
]
EVIDENCE_FIELDS = (
    [
        "gene_id",
        "transcript_id",
        "seqname",
        "start",
        "end",
        "strand",
        "selection_source",
        "biotype_source",
    ]
    + FLAGS
    + DIAGNOSTIC_FIELDS
    + ["homology_status"]
    + HOMOLOGY_FIELDS
    + ["confidence", "confidence_rule"]
)
INTRON_FIELDS = [
    "gene_id",
    "transcript_id",
    "seqname",
    "strand",
    "intron_number",
    "start",
    "end",
    "length",
    "is_cds_intron",
    "splice_motif",
    "canonical",
    "regular_size",
    "read_count",
    "supported",
]
VALIDATION_FIELDS = [
    "gene_id",
    "transcript_id",
    "seqname",
    "feature_type",
    "start",
    "end",
    "strand",
    "contig_length",
    "bases_beyond_end",
    "reason",
    "policy",
]
NORMALIZATION_FIELDS = ["gene_id", "transcript_id", "action", "details"]
EXCEPTION_FIELDS = ["gene_id", "transcript_id", "policy", "reason"]
OUTPUT_SUFFIXES = {
    "evidence": ".evidence.tsv",
    "introns": ".introns.tsv",
    "cds": ".cds.fa",
    "proteins": ".prots.fa",
    "gff": ".annotations.gff3",
    "manifest": ".manifest.json",
    "translation_exceptions": ".translation_exceptions.tsv",
    "normalization": ".normalization.tsv",
    "input_validation": ".input_validation.tsv",
}


# Evidence-summary schema.
STATES = ("pass", "fail", "unknown")
CLASSES = ("high", "medium", "low")
BINARY = FLAGS + ["cds_length_multiple_of_three"]
NUMERIC = [
    "cds_length",
    "intron_count",
    "cds_intron_count",
    "supported_cds_intron_count",
    "internal_stop_count",
    "ambiguous_amino_acids",
    "selenocysteine_count",
    "inferred_cds_phase_count",
    "trimmed_terminal_bases",
    "error",
    "perc_id",
    "score",
    "query_coverage",
    "target_coverage",
]
CATEGORICAL = [
    "confidence",
    "confidence_rule",
    "homology_status",
    "selection_source",
    "biotype_source",
    "translation_status",
    "translation_table",
    "translation_table_source",
    "exon_structure_source",
    "cds_phase_source",
    "stop_codon_source",
    "utr_annotation_source",
    "initiator_methionine_adjusted",
]
MISSING = {"", ".", "NA", "None", "null"}


# Diagram labels and class colours.
RULES = {
    "utr_ignored_error_le_10": "Legacy no_utr rule: both intron checks; error ≤10",
    "both_utrs_error_le_10": "Both UTRs + both intron checks; error ≤10",
    "error_le_5_score_ge_80": "Both intron checks; error ≤5; score ≥80",
    "supported_error_le_5_score_ge_80": "Both intron checks + RNA support; error ≤5; score ≥80",
    "supported_error_le_10": "Both intron checks + RNA support; error ≤10",
    "supported_one_intron_check": "≥1 intron check + RNA support; error ≤10",
    "strong_homology_one_intron_check": "≥1 intron check; error ≤5; score ≥80",
    "no_high_or_medium_rule": "No high/medium rule satisfied",
}
COLORS = {"high": "#269585", "medium": "#E5AF38", "low": "#CE5D63"}

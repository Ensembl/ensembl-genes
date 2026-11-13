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
"""Pure sequence reconstruction, evidence assessment and confidence policy v2. No
files are opened by this module."""

from __future__ import annotations

import math
from copy import deepcopy
from dataclasses import replace
from functools import lru_cache
import pandas as pd
from ensembl.genes.annotation_qc.config.high_confidence import (
    COMPLEMENT,
    DEFAULT_MITOCHONDRIAL_CONTIGS,
    CANONICAL_SPLICES,
    FLAGS,
    HOMOLOGY_FIELDS,
    MIN_INTRON_LENGTH,
    MAX_INTRON_LENGTH,
    HOMOLOGY_ERROR_MAX,
    STRONG_HOMOLOGY_ERROR_MAX,
    STRONG_HOMOLOGY_SCORE_MIN,
)
from ensembl.genes.annotation_qc.parsers._high_confidence_records import (
    Feature,
    Transcript,
    Gene,
    SequenceRegion,
)
from ensembl.genes.annotation_qc.parsers.high_confidence import single_attr


from ensembl.genes.annotation_qc.config.high_confidence import (
    EVIDENCE_FIELDS,
    INTRON_FIELDS,
)


class TranslationExceptionError(ValueError):
    """An explicitly annotated translation event that cannot be reconstructed."""


@lru_cache(maxsize=1)
def biopython_api():
    """Import lazily so help and scoring-only comparison tools remain usable."""
    # Keep help and evidence-only classification usable without loading Bio.
    # pylint: disable=import-outside-toplevel
    try:
        import Bio
        from Bio.Data import CodonTable
        from Bio.Seq import Seq
    except ImportError as error:
        raise ValueError(
            "Biopython is required for translation. Install with "
            "the ensembl-genes project dependencies using the same Python interpreter."
        ) from error
    return Seq, CodonTable, Bio.__version__


def translation_table_for(tx, mitochondrial_contigs=None):
    """Resolve vertebrate codon table from annotation or contig."""
    declared = single_attr(tx.attrs, "transl_table", "codon_table")
    if declared is not None:
        if declared not in {"1", "2"}:
            raise TranslationExceptionError(
                f"{tx.id}: unsupported translation table {declared!r}; "
                "this vertebrate workflow supports tables 1 and 2"
            )
        return int(declared), "gtf_attribute"
    mitochondrial_contigs = (
        DEFAULT_MITOCHONDRIAL_CONTIGS
        if mitochondrial_contigs is None
        else mitochondrial_contigs
    )
    if tx.seqname in mitochondrial_contigs:
        return 2, "mitochondrial_contig"
    return 1, "nuclear_default"


def translate_cds(coding_sequence, table=1, selenocysteine_codons=()):
    """Match Ensembl's core translation steps with Biopython codon tables.

    Trim a trailing incomplete codon, preserve internal stops, normalize a
    recognized initiator to M, recode annotated Sec, and omit a terminal stop.
    cds=False permits partial/diagnostic gene models and recoding exceptions.
    """
    if table not in {1, 2}:
        raise TranslationExceptionError(
            f"Unsupported translation table {table}; expected 1 or 2"
        )
    Seq, codon_tables, _ = biopython_api()
    coding_sequence = coding_sequence.upper()
    trimmed_bases = len(coding_sequence) % 3
    sequence = (
        coding_sequence[: len(coding_sequence) - trimmed_bases]
        if trimmed_bases
        else coding_sequence
    )
    amino_acids = list(
        str(Seq(sequence).translate(table=table, cds=False, to_stop=False))
    )
    start_codons = codon_tables.unambiguous_dna_by_id[table].start_codons
    adjusted_start = bool(
        amino_acids and amino_acids[0] != "M" and sequence[:3] in start_codons
    )
    if amino_acids and sequence[:3] in start_codons:
        amino_acids[0] = "M"
    for index in set(selenocysteine_codons):
        if (
            not 0 <= index < len(amino_acids)
            or sequence[index * 3 : index * 3 + 3] != "TGA"
        ):
            raise TranslationExceptionError(
                "Annotated Selenocysteine must map to a complete in-frame TGA codon"
            )
        amino_acids[index] = "U"
    protein = "".join(amino_acids)
    if protein.endswith("*"):
        protein = protein[:-1]
    return {
        "protein": protein,
        "translation_table": table,
        "initiator_methionine_adjusted": adjusted_start,
        "trimmed_terminal_bases": trimmed_bases,
    }


def reverse_complement(sequence: str) -> str:
    """Reverse and complement a DNA sequence."""
    return sequence.upper().translate(COMPLEMENT)[::-1]


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-branches
def select_transcripts(
    genes: dict[str, Gene], overrides: dict, policy: str, input_mode="auto"
):
    """Select one coding representative per gene under the chosen policy."""
    selected = []
    for gene in sorted(genes.values(), key=lambda g: (g.seqname, g.start, g.id)):
        biotype = single_attr(gene.attrs, "gene_biotype", "gene_type")
        if biotype is None:
            transcript_types = {
                single_attr(t.attrs, "transcript_biotype", "transcript_type")
                for t in gene.transcripts.values()
            }
            transcript_types.discard(None)
            if not transcript_types:
                if input_mode == "strict":
                    raise ValueError(
                        f"Gene {gene.id} has no gene or transcript biotype"
                    )
                has_cds = any(t.of_type("CDS") for t in gene.transcripts.values())
                biotype = "protein_coding" if has_cds else "unclassified"
                biotype_source = (
                    "cds_presence_inferred" if has_cds else "unclassified_without_cds"
                )
                note = {
                    "action": biotype_source,
                    "details": "No biotype was supplied; "
                    + (
                        "CDS features identify a coding model"
                        if has_cds
                        else "no CDS features; omitted from coding-gene assessment"
                    ),
                }
                if note not in gene.normalization:
                    gene.normalization.append(note)
            else:
                biotype = (
                    "protein_coding"
                    if "protein_coding" in transcript_types
                    else "noncoding"
                )
                biotype_source = "transcript_biotype_fallback"
        else:
            biotype_source = "gene_biotype"
        if biotype != "protein_coding":
            continue
        candidates = list(gene.transcripts.values())
        tagged = [
            t for t in candidates if "Ensembl_canonical" in t.attrs.get("tag", [])
        ]
        if gene.id in overrides:
            chosen, source = gene.transcripts[overrides[gene.id]], "canonical_map"
        elif len(tagged) == 1:
            chosen, source = tagged[0], "Ensembl_canonical_tag"
        elif len(tagged) > 1:
            raise ValueError(
                f"Multiple Ensembl_canonical transcripts for {gene.id}; supply a canonical map"
            )
        elif len(candidates) == 1:
            chosen, source = candidates[0], "single_transcript"
        elif policy in {"auto", "longest-cds"} and candidates:
            chosen = min(
                candidates,
                key=lambda t: (
                    -sum(f.end - f.start + 1 for f in t.of_type("CDS")),
                    -sum(f.end - f.start + 1 for f in t.of_type("exon")),
                    t.id,
                ),
            )
            source = "longest_cds_fallback"
        else:
            raise ValueError(
                f"No unambiguous canonical transcript for {gene.id}; "
                "supply --canonical-map or opt into --canonical-policy longest-cds"
            )
        selected.append((gene, chosen, source, biotype_source))
    if not selected:
        raise ValueError("No protein coding genes selected")
    return selected


def merge_spans(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjoining half-open intervals."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(merged[-1][1], end)
        else:
            merged.append((start, end))
    return merged


def resolve_cds_phases(tx, cds, input_mode):
    """Validate CDS phases and infer missing phases consistently."""
    missing = sum(f.phase is None for f in cds)
    if missing and input_mode == "strict":
        raise ValueError(f"{tx.id}: CDS phase is required in strict mode")
    offset, first_phase = 0, None
    for feature in cds:
        if feature.phase is not None:
            first_phase = (feature.phase + offset) % 3
            break
        offset += feature.end - feature.start + 1
    source = "annotated" if not missing else "inferred_from_provided_phases"
    if first_phase is None:
        first_phase, source = 0, "assumed_first_phase_zero"
    result, consumed = [], 0
    for feature in cds:
        expected = (first_phase - consumed) % 3
        if feature.phase is not None and feature.phase != expected:
            raise ValueError(
                f"{tx.id}: inconsistent CDS phase at {feature.start}-{feature.end}; "
                f"expected {expected}, found {feature.phase}"
            )
        result.append(
            replace(feature, phase=expected) if feature.phase is None else feature
        )
        consumed += feature.end - feature.start + 1
    return result, source, missing


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-locals
def prepare_exon_model(tx, genome, cds, table, input_mode):
    """Build a minimal spliced model when only coding/UTR features are present."""
    exons, stops = tx.of_type("exon"), tx.of_type("stop_codon")
    notes = []
    if not exons and input_mode == "strict":
        raise ValueError(f"{tx.id}: selected transcript requires exons in strict mode")
    source = "annotated" if exons else "inferred_from_features"
    stop_source = "explicit_feature" if stops else "not_annotated"
    if not stops:
        _, codon_tables, _ = biopython_api()
        stop_codons = set(codon_tables.unambiguous_dna_by_id[table].stop_codons)
        sequence = "N" * ((3 - cds[0].phase) % 3) + "".join(
            genome.fetch(tx.seqname, f.start, f.end, tx.strand) for f in cds
        )
        if len(sequence) % 3 == 0 and sequence[-3:] in stop_codons:
            stop_source = "within_cds"
        elif (
            not exons
            and input_mode == "auto"
            and len(sequence) % 3 == 0
            and not any(
                f.kind in {"UTR", "five_prime_utr", "three_prime_utr"}
                for f in tx.features
            )
        ):
            headers = tx.of_type("transcript")
            bounds = (
                (headers[0].start, headers[0].end)
                if len(headers) == 1
                else tx.parent_bounds
            )
            terminal = cds[-1]
            candidate = (
                (terminal.end + 1, terminal.end + 3)
                if tx.strand == "+"
                else (terminal.start - 3, terminal.start - 1)
            )
            if (
                bounds
                and bounds[0] <= candidate[0] <= candidate[1] <= bounds[1]
                and genome.contains(tx.seqname, *candidate)
            ):
                codon = genome.fetch(tx.seqname, *candidate, tx.strand)
                if codon in stop_codons:
                    stops = [
                        Feature(
                            tx.seqname, "stop_codon", *candidate, tx.strand, None, {}
                        )
                    ]
                    stop_source = "inferred_from_genome"
                    notes.append(
                        {
                            "action": "stop_codon_inferred",
                            "details": (
                                f"Verified {codon} at {tx.seqname}:"
                                f"{candidate[0]}-{candidate[1]}:{tx.strand} "
                                "immediately after CDS, within declared parent bounds"
                            ),
                        }
                    )
    if not exons:
        structural = [
            f
            for f in tx.features
            if f.kind
            in {
                "CDS",
                "start_codon",
                "stop_codon",
                "UTR",
                "five_prime_utr",
                "three_prime_utr",
                "Selenocysteine",
            }
        ]
        structural += stops
        spans = merge_spans([(f.start, f.end + 1) for f in structural])
        exons = [
            Feature(tx.seqname, "exon", start, end - 1, tx.strand, None, {})
            for start, end in spans
        ]
        exons.sort(key=lambda f: f.start, reverse=tx.strand == "-")
        notes.append(
            {
                "action": "exons_inferred",
                "details": (
                    f"Built {len(exons)} minimal exon blocks from CDS/codon/UTR "
                    "features; unannotated UTR exons and transcript ends cannot be recovered"
                ),
            }
        )
    return exons, stops, source, stop_source, notes


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-arguments,too-many-branches,too-many-locals,too-many-statements
def assess_transcript(
    tx: Transcript,
    genome: SequenceRegion,
    junctions,
    minimum_reads: int,
    mitochondrial_contigs=None,
    input_mode="auto",
):
    """Reconstruct CDS and evaluate translation, UTR and intron evidence."""
    tags = set(tx.attrs.get("tag", []))
    selenocysteines = [f for f in tx.features if f.kind.lower() == "selenocysteine"]
    unsupported = [f"tag={tag}" for tag in sorted(tags & {"stop_codon_readthrough"})]
    for key in ("transl_except", "exception", "_rna_edit"):
        if key in tx.attrs:
            unsupported.append(f"{key}=" + "|".join(tx.attrs[key]))
    seleno_tags = tags & {"selenocysteine", "seleno"}
    if seleno_tags and not selenocysteines:
        unsupported.append(
            "tag="
            + "|".join(sorted(seleno_tags))
            + " without Selenocysteine feature coordinates"
        )
    if unsupported:
        raise TranslationExceptionError(
            f"{tx.id}: translation exception requires an explicit translation adapter: "
            + "; ".join(unsupported)
            + ". Supply the annotation details for an adapter, or use "
            "--translation-exception-policy skip to exclude and report this transcript."
        )
    translation_table, table_source = translation_table_for(tx, mitochondrial_contigs)
    cds = tx.of_type("CDS")
    if not cds:
        raise ValueError(f"{tx.id}: selected transcript requires CDS features")
    cds, phase_source, inferred_phase_count = resolve_cds_phases(tx, cds, input_mode)
    exons, stop_features, exon_source, stop_source, normalization = prepare_exon_model(
        tx, genome, cds, translation_table, input_mode
    )
    if inferred_phase_count:
        normalization.append(
            {
                "action": "cds_phases_inferred",
                "details": (
                    f"Filled {inferred_phase_count} missing CDS phases; source={phase_source}"
                ),
            }
        )
    genomic_exons = sorted(exons, key=lambda f: f.start)
    if any(
        left.end + 1 >= right.start
        for left, right in zip(genomic_exons, genomic_exons[1:])
    ):
        raise ValueError(f"{tx.id}: overlapping, duplicate or directly adjacent exons")
    transcript_features = tx.of_type("transcript")
    if len(transcript_features) > 1:
        raise ValueError(f"{tx.id}: duplicate transcript features")
    if transcript_features:
        header = transcript_features[0]
        if (
            not header.start
            <= genomic_exons[0].start
            <= genomic_exons[-1].end
            <= header.end
        ):
            raise ValueError(f"{tx.id}: exons outside declared transcript bounds")
        if (header.start, header.end) != (
            genomic_exons[0].start,
            genomic_exons[-1].end,
        ):
            if input_mode == "strict":
                raise ValueError(f"{tx.id}: transcript bounds differ from exon bounds")
            normalization.append(
                {
                    "action": "transcript_bounds_not_filled",
                    "details": (
                        "Declared transcript bounds exceed known exon blocks; "
                        "no UTR sequence invented"
                    ),
                }
            )
    offsets, sequence_parts, total = [], [], 0
    for exon in exons:
        offsets.append((exon, total, total + exon.end - exon.start + 1))
        sequence_parts.append(genome.fetch(tx.seqname, exon.start, exon.end, tx.strand))
        total = offsets[-1][2]
    mrna = "".join(sequence_parts)

    def project(feature):
        for exon, offset, _ in offsets:
            if exon.start <= feature.start <= feature.end <= exon.end:
                if tx.strand == "+":
                    return (
                        offset + feature.start - exon.start,
                        offset + feature.end - exon.start + 1,
                    )
                return (
                    offset + exon.end - feature.end,
                    offset + exon.end - feature.start + 1,
                )
        raise ValueError(
            f"{tx.id}: {feature.kind} {feature.start}-{feature.end} is not inside an exon"
        )

    def genomic_spans(start, end):
        result = []
        for exon, a, b in offsets:
            lo, hi = max(start, a), min(end, b)
            if lo >= hi:
                continue
            if tx.strand == "+":
                result.append((exon.start + lo - a, exon.start + hi - a - 1))
            else:
                result.append((exon.end - (hi - a) + 1, exon.end - (lo - a)))
        return result

    cds_spans = [project(f) for f in cds]
    if any(right[0] < left[1] for left, right in zip(cds_spans, cds_spans[1:])):
        raise ValueError(f"{tx.id}: overlapping CDS features")
    first_phase = cds[0].phase
    stop_spans = [project(f) for f in stop_features]
    merged = merge_spans(cds_spans + stop_spans)
    if len(merged) != 1:
        raise ValueError(
            f"{tx.id}: CDS and stop codon do not form one continuous spliced coding span"
        )
    coding_start, coding_end = merged[0]
    if stop_spans:
        if sum(b - a for a, b in stop_spans) != 3 or merge_spans(stop_spans) != [
            (coding_end - 3, coding_end)
        ]:
            raise ValueError(
                f"{tx.id}: stop_codon must cover exactly the final three coding bases"
            )
    padding = "N" * ((3 - first_phase) % 3)
    coding_sequence = padding + mrna[coding_start:coding_end]
    # Project annotated Sec bases onto the spliced CDS. A Sec codon may span
    # multiple GTF rows/exons; nearby distinct Sec codons must remain distinct.
    sec_bases: set[int] = set()
    for feature in selenocysteines:
        start, end = project(feature)
        if not coding_start <= start < end <= coding_end:
            raise TranslationExceptionError(
                f"{tx.id}: Selenocysteine feature lies outside the coding span"
            )
        sec_bases.update(
            range(
                start - coding_start + len(padding), end - coding_start + len(padding)
            )
        )
    sec_codons = sorted({position // 3 for position in sec_bases})
    for codon_index in sec_codons:
        position = codon_index * 3
        if not set(range(position, position + 3)) <= sec_bases:
            raise TranslationExceptionError(
                f"{tx.id}: Selenocysteine annotation does not cover a complete in-frame codon"
            )
        if coding_sequence[position : position + 3] != "TGA":
            raise TranslationExceptionError(
                f"{tx.id}: annotated Selenocysteine codon is "
                f"{coding_sequence[position:position + 3]!r}, expected TGA; "
                "check annotation/genome assembly and phase"
            )
    translation = translate_cds(
        coding_sequence, table=translation_table, selenocysteine_codons=sec_codons
    )
    protein = translation["protein"]
    introns = []
    for (left, _, boundary), (right, _, _) in zip(offsets, offsets[1:]):
        start, end = min(left.end, right.end) + 1, max(left.start, right.start) - 1
        if tx.strand == "+":
            donor = genome.fetch(tx.seqname, start, min(start + 1, end))
            acceptor = genome.fetch(tx.seqname, max(start, end - 1), end)
        else:
            donor = genome.fetch(tx.seqname, max(start, end - 1), end, "-")
            acceptor = genome.fetch(tx.seqname, start, min(start + 1, end), "-")
        motif = donor + "-" + acceptor
        count = (
            None
            if junctions is None
            else junctions.get((tx.seqname, start, end, tx.strand), 0)
        )
        introns.append(
            {
                "start": start,
                "end": end,
                "length": end - start + 1,
                "is_cds_intron": coding_start < boundary < coding_end,
                "splice_motif": motif,
                "canonical": motif in CANONICAL_SPLICES,
                "regular_size": MIN_INTRON_LENGTH
                <= end - start + 1
                <= MAX_INTRON_LENGTH,
                "read_count": count,
                "supported": None if count is None else count >= minimum_reads,
            }
        )
    coding_introns = [i for i in introns if i["is_cds_intron"]]
    supported = (
        True
        if not coding_introns
        else None if junctions is None else all(i["supported"] for i in coding_introns)
    )
    evidence = {
        "transcript_cds_complete": coding_sequence.startswith("ATG")
        and coding_sequence[-3:] in {"TAA", "TAG", "TGA"},
        "transcript_five_prime_utr": coding_start > 0,
        "transcript_three_prime_utr": coding_end < len(mrna),
        "transcript_canonical_splicing": all((i["canonical"] for i in introns)),
        "transcript_regular_intron_size": all((i["regular_size"] for i in introns)),
        "transcript_intron_support": supported,
        "cds_length": len(coding_sequence),
        "cds_length_multiple_of_three": len(coding_sequence) % 3 == 0,
        "internal_stop_count": protein.count("*"),
        "ambiguous_amino_acids": protein.count("X"),
        "translation_status": (
            "selenocysteine_from_gtf"
            if sec_codons
            else "standard" if translation_table == 1 else "vertebrate_mitochondrial"
        ),
        "translation_table": translation_table,
        "translation_table_source": table_source,
        "initiator_methionine_adjusted": translation["initiator_methionine_adjusted"],
        "trimmed_terminal_bases": translation["trimmed_terminal_bases"],
        "selenocysteine_count": len(sec_codons),
        "exon_structure_source": exon_source,
        "cds_phase_source": phase_source,
        "inferred_cds_phase_count": inferred_phase_count,
        "stop_codon_source": stop_source,
        "utr_annotation_source": (
            "exon_geometry"
            if exon_source == "annotated"
            else (
                "explicit_utr_features"
                if any(
                    (
                        f.kind in {"UTR", "five_prime_utr", "three_prime_utr"}
                        for f in tx.features
                    )
                )
                else "unavailable_cds_only"
            )
        ),
        "intron_count": len(introns),
        "cds_intron_count": len(coding_introns),
        "supported_cds_intron_count": (
            None
            if junctions is None and coding_introns
            else sum((i["supported"] is True for i in coding_introns))
        ),
    }
    return {
        "evidence": evidence,
        "cds": coding_sequence,
        "protein": protein,
        "introns": introns,
        "exons": exons,
        "cds_spans": genomic_spans(coding_start, coding_end),
        "first_phase": first_phase,
        "normalization": normalization,
        "five_prime_UTR": genomic_spans(0, coding_start),
        "three_prime_UTR": genomic_spans(coding_end, len(mrna)),
    }


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-return-statements
def confidence_rating(
    evidence: dict, homology: dict | None, no_utr: bool = False
) -> tuple[str, str]:
    """Preserve Perl defaults; no_utr disables the UTR route and requires RNA support for high."""
    complete, five, three, canonical, regular, support = (evidence[f] for f in FLAGS)
    if not complete or homology is None:
        return "low", "incomplete_cds_or_no_homology"
    error, score = homology["error"], homology["score"]
    if (
        # Keep the established reconstruction/scoring flow together for parity.
        # pylint: disable-next=too-many-boolean-expressions
        not no_utr
        and regular
        and canonical
        and five
        and three
        and error <= HOMOLOGY_ERROR_MAX
    ):
        return "high", "both_utrs_error_le_10"
    if (
        # Keep the established reconstruction/scoring flow together for parity.
        # pylint: disable-next=too-many-boolean-expressions
        regular
        and canonical
        and error <= STRONG_HOMOLOGY_ERROR_MAX
        and score >= STRONG_HOMOLOGY_SCORE_MIN
        and (not no_utr or support)
    ):
        return "high", (
            "supported_error_le_5_score_ge_80" if no_utr else "error_le_5_score_ge_80"
        )
    if regular and canonical and error <= HOMOLOGY_ERROR_MAX and support:
        return "high", "supported_error_le_10"
    if (regular or canonical) and error <= HOMOLOGY_ERROR_MAX and support:
        return "medium", "supported_one_intron_check"
    if (
        (regular or canonical)
        and error <= STRONG_HOMOLOGY_ERROR_MAX
        and score >= STRONG_HOMOLOGY_SCORE_MIN
    ):
        return "medium", "strong_homology_one_intron_check"
    return "low", "no_high_or_medium_rule"


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-locals
def validate_model_bounds(selected, genome):
    """Check all selected feature and gene intervals without fetching sequence."""
    issues, invalid = [], set()
    for gene, tx, _, _ in selected:
        intervals = [("gene", gene.seqname, gene.start, gene.end, gene.strand)]
        intervals.extend(
            (f.kind, f.seqname, f.start, f.end, f.strand) for f in tx.features
        )
        seen = set()
        for kind, seqname, start, end, orientation in intervals:
            key = kind, seqname, start, end, orientation
            if key in seen:
                continue
            seen.add(key)
            entry = genome.index.get(seqname)
            if entry is not None and genome.contains(seqname, start, end):
                continue
            length = entry[0] if entry is not None else None
            reason = (
                "contig_missing_from_fasta"
                if entry is None
                else "interval_outside_fasta"
            )
            issues.append(
                {
                    "gene_id": gene.id,
                    "transcript_id": tx.id,
                    "seqname": seqname,
                    "feature_type": kind,
                    "start": start,
                    "end": end,
                    "strand": orientation,
                    "contig_length": length,
                    "bases_beyond_end": (
                        max(0, end - length) if length is not None else None
                    ),
                    "reason": reason,
                }
            )
            invalid.add((gene.id, tx.id))
    return issues, invalid


def _models_from_dataframe(annotation):
    """Private compatibility adapter preserving proven frame/stop/UTR behavior.

    This is not a second public annotation model: parser/runner boundaries use
    DataFrames. Records remain private to the parity-sensitive reconstruction.
    """
    genes = {}
    for gene_id, metadata in annotation.attrs.get("gene_metadata", {}).items():
        metadata = deepcopy(metadata)
        genes[gene_id] = Gene(
            gene_id,
            metadata["seqname"],
            metadata["strand"],
            metadata["start"],
            metadata["end"],
            attrs=metadata["attrs"],
            declared=metadata["declared"],
            normalization=metadata["normalization"],
        )
    for row in annotation.to_dict("records"):
        if row["Feature"] == "gene":
            continue
        gene = genes[row["gene_id"]]
        tx_id = row["transcript_id"]
        if tx_id not in gene.transcripts:
            metadata = deepcopy(annotation.attrs["transcript_metadata"][gene.id, tx_id])
            gene.transcripts[tx_id] = Transcript(
                tx_id,
                gene.id,
                row["Chromosome"],
                row["Strand"],
                attrs=metadata["attrs"],
                parent_bounds=metadata["parent_bounds"],
                normalization=metadata["normalization"],
            )
        phase = row["Frame"]
        phase = None if phase is None or pd.isna(phase) else int(phase)
        gene.transcripts[tx_id].features.append(
            Feature(
                row["Chromosome"],
                row["Feature"],
                int(row["Start"]) + 1,
                int(row["End"]),
                row["Strand"],
                phase,
                deepcopy(row["Attributes"]),
            )
        )
    return genes


def select_representative_transcripts(
    annotation, assignments, policy="auto", input_mode="auto"
):
    """Select one representative per coding gene from normalized DataFrames."""
    genes = _models_from_dataframe(annotation)
    overrides = {}
    for row in assignments.to_dict("records"):
        gene_id, tx_id = row["gene_id"], row["transcript_id"]
        if (
            gene_id in overrides
            or gene_id not in genes
            or tx_id not in genes[gene_id].transcripts
        ):
            raise ValueError(
                f"Invalid or duplicate canonical assignment: {gene_id} -> {tx_id}"
            )
        overrides[gene_id] = tx_id
    return genes, select_transcripts(genes, overrides, policy, input_mode)


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-locals
def score_diamond_alignments(alignments, query_lengths):
    """Reproduce the legacy best-hit arithmetic using parsed alignments only."""
    best, scores = {}, {}
    for number, row in enumerate(alignments.to_dict("records"), 1):
        query, hit = row["qseqid"], row["sseqid"]
        values = [
            float(row[k])
            for k in ["qstart", "qend", "slen", "pident", "nident", "length", "gaps"]
        ]
        if (
            query not in query_lengths
            or not hit
            or not all(math.isfinite(v) for v in values)
        ):
            raise ValueError(
                f"Alignment {number}: unknown query or invalid alignment value"
            )
        start, end, target_length, _, identical, length, gaps = values
        if not (
            1 <= start <= end <= query_lengths[query]
            and target_length > 0
            and identical >= 0
            and 0 <= gaps < length
        ):
            raise ValueError(
                f"Alignment {number}: invalid coordinates, lengths or counts"
            )
        span = end - start + 1
        query_coverage = span / query_lengths[query] * 100
        target_coverage = span / target_length * 100
        error = abs(query_coverage - 100) + abs(target_coverage - 100)
        identity = identical / (length - gaps) * 100
        score = 0.25 * identity + 0.75 * (100 - error)
        if query not in scores or score > scores[query]:
            scores[query] = score
            values = [error, identity, score, query_coverage, target_coverage]
            best[query] = dict(
                zip(HOMOLOGY_FIELDS, [hit] + [float(f"{v:.2f}") for v in values])
            )
    return best


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-locals
def build_confidence_tables(results, homology, no_utr=False, homology_measured=False):
    """Return computed evidence/intron DataFrames and gene output records, without I/O."""

    rows, introns, gff_records = [], [], []
    counts = {"high": 0, "medium": 0, "low": 0}
    for result in results:
        gene, tx = result["gene"], result["tx"]
        hit = homology.get(tx.id)
        confidence, rule = confidence_rating(result["evidence"], hit, no_utr=no_utr)
        counts[confidence] += 1
        # Keyword construction rejects duplicate export-schema keys.
        # pylint: disable-next=use-dict-literal
        row = dict(
            gene_id=gene.id,
            transcript_id=tx.id,
            seqname=gene.seqname,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            selection_source=result["selection_source"],
            biotype_source=result["biotype_source"],
            **result["evidence"],
            homology_status=(
                "not_measured" if not homology_measured else "hit" if hit else "no_hit"
            ),
            **(hit or {key: None for key in HOMOLOGY_FIELDS}),
            confidence=confidence,
            confidence_rule=rule,
        )
        rows.append(row)
        for number, intron in enumerate(result["introns"], 1):
            introns.append(
                # Keyword construction rejects duplicate export-schema keys.
                # pylint: disable-next=use-dict-literal
                dict(
                    gene_id=gene.id,
                    transcript_id=tx.id,
                    seqname=tx.seqname,
                    strand=tx.strand,
                    intron_number=number,
                    **intron,
                )
            )
        record = dict(result, confidence=confidence, confidence_rule=rule)
        consumed = 0
        phases = []
        for start, end in result["cds_spans"]:
            phases.append((result["first_phase"] - consumed) % 3)
            consumed += end - start + 1
        record["gff_cds_phases"] = phases
        gff_records.append(record)
    return {
        "evidence": pd.DataFrame(rows, columns=EVIDENCE_FIELDS, dtype=object),
        "introns": pd.DataFrame(introns, columns=INTRON_FIELDS, dtype=object),
        "gff_records": gff_records,
        "confidence_counts": counts,
    }


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-branches,too-many-locals
def rescore_evidence(frame, no_utr=False, gff=None):
    """Reclassify parsed evidence and return new tables without reading/writing files."""
    required = set(
        FLAGS
        + [
            "gene_id",
            "transcript_id",
            "homology_status",
            "error",
            "score",
            "confidence",
            "confidence_rule",
        ]
    )
    if not required <= set(frame.columns):
        raise ValueError(
            f"Missing evidence columns: {sorted(required-set(frame.columns))}"
        )
    for column in ["gene_id", "transcript_id"]:
        if frame[column].duplicated().any():
            raise ValueError(
                f"Duplicate {column}; expected one selected transcript per gene"
            )
    rows = frame.to_dict("records")
    changes = []
    counts = {"high": 0, "medium": 0, "low": 0}
    old_counts = {"high": 0, "medium": 0, "low": 0}
    for row in rows:
        evidence = {}
        for flag in FLAGS:
            value = row[flag]
            if value not in {"1", "0", ".", "", "NA", "None", "null"}:
                raise ValueError(f"Invalid flag {flag}={value!r}")
            evidence[flag] = True if value == "1" else False if value == "0" else None
        if row["homology_status"] not in {"hit", "no_hit", "not_measured"}:
            raise ValueError(f'Invalid homology status: {row["homology_status"]}')
        hit = None
        if row["homology_status"] == "hit":
            hit = {key: float(row[key]) for key in ["error", "score"]}
            if not all(math.isfinite(value) for value in hit.values()):
                raise ValueError("Homology metrics must be finite")
        old, old_rule = row["confidence"], row["confidence_rule"]
        if old not in old_counts:
            raise ValueError(f"Invalid original confidence {old!r}")
        confidence, rule = confidence_rating(evidence, hit, no_utr=no_utr)
        old_counts[old] += 1
        counts[confidence] += 1
        if (old, old_rule) != (confidence, rule):
            changes.append(
                {
                    "gene_id": row["gene_id"],
                    "transcript_id": row["transcript_id"],
                    "old_confidence": old,
                    "new_confidence": confidence,
                    "old_rule": old_rule,
                    "new_rule": rule,
                    "class_changed": int(old != confidence),
                }
            )
        row.update(confidence=confidence, confidence_rule=rule)
    if gff is not None:
        genes = [
            row["gene_id"]
            for row in gff.to_dict("records")
            if row["gene_id"] is not None
        ]
        if len(set(genes)) != len(genes) or set(genes) != set(frame["gene_id"]):
            raise ValueError(
                "Exported GFF genes do not match the evidence table exactly"
            )
    return {
        "evidence": pd.DataFrame(rows, columns=frame.columns, dtype=object),
        "changes": changes,
        "confidence_counts": counts,
        "old_confidence_counts": old_counts,
        "class_changed_count": sum((r["class_changed"] for r in changes)),
    }


def compute_run_diagnostics(results):
    """Summarize already computed translation diagnostics without I/O."""
    return {
        "selenocysteine_transcript_count": sum(
            (r["evidence"]["selenocysteine_count"] > 0 for r in results)
        ),
        "translation_table_counts": {
            str(table): sum(
                (r["evidence"]["translation_table"] == table for r in results)
            )
            for table in (1, 2)
        },
        "initiator_methionine_adjusted_count": sum(
            (r["evidence"]["initiator_methionine_adjusted"] for r in results)
        ),
    }

"""Derive and assess sequence-aware annotation events.

The public functions in this module take a standardized annotation table (or
PyRanges object) and, when needed, an indexed genome FASTA. They return
DataFrames: one row per intron/junction, transcript translation, or transcript
boundary. Interpretation is conservative when sequence is unavailable or a
motif profile has not been supplied.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from Bio.Data import CodonTable

from ensembl.genes.annotation_qc.parsers.annotation import standardize_annotation

STANDARD_CODE = {
    "table_id": 1,
    "name": CodonTable.unambiguous_dna_by_id[1].names[0],
    "stop_codons": frozenset(CodonTable.unambiguous_dna_by_id[1].stop_codons),
}

NEAR_COGNATE_STARTS = frozenset({"CTG", "GTG", "TTG", "ACG"})
DEFAULT_CANONICAL_JUNCTIONS = frozenset({"GT-AG", "GC-AG", "AT-AC"})
DEFAULT_POLYA_MOTIFS = {"AATAAA": "AATAAA", "ATTAAA": "ATTAAA"}
EMPTY_EVENT_COLUMNS = [
    "event_type",
    "transcript_id",
    "gene_id",
    "chromosome",
    "strand",
]
_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def get_genetic_code(table_id: int) -> dict:
    """Return a plain dictionary describing an NCBI genetic code."""
    try:
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError as exc:
        raise ValueError(f"Unknown NCBI genetic code table: {table_id}") from exc
    return {
        "table_id": table_id,
        "name": table.names[0],
        "stop_codons": frozenset(table.stop_codons),
    }


def _annotation_df(annotation) -> pd.DataFrame:
    """Standardize a PyRanges/DataFrame input for event calculations."""
    if isinstance(annotation, pd.DataFrame) and {
        "gene_id",
        "transcript_id",
        "feature_length",
    }.issubset(annotation.columns):
        return annotation
    return standardize_annotation(annotation)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_RC_TABLE)[::-1]


def _empty_events(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=[*EMPTY_EVENT_COLUMNS, *columns])


def _sequence(fasta, chromosome: str, start: int, end: int, strand: str) -> str:
    """Extract a genomic interval in transcript orientation."""
    if start >= end:
        return ""
    sequence = fasta[chromosome][start:end].seq
    return sequence if strand != "-" else _reverse_complement(sequence)


def _transcript_boundaries(annotation: pd.DataFrame) -> pd.DataFrame:
    """Return one transcript boundary row, using exons when needed."""
    transcript_rows = annotation[
        annotation["Feature"].isin({"transcript", "mRNA"})
    ].copy()
    if not transcript_rows.empty:
        return transcript_rows[
            ["transcript_id", "gene_id", "Chromosome", "Start", "End", "Strand"]
        ].drop_duplicates("transcript_id")

    exon_rows = annotation[annotation["Feature"] == "exon"]
    if exon_rows.empty:
        return pd.DataFrame(
            columns=["transcript_id", "gene_id", "Chromosome", "Start", "End", "Strand"]
        )
    return (
        exon_rows.groupby("transcript_id", dropna=True)
        .agg(
            gene_id=("gene_id", "first"),
            Chromosome=("Chromosome", "first"),
            Start=("Start", "min"),
            End=("End", "max"),
            Strand=("Strand", "first"),
        )
        .reset_index()
    )


def extract_cds_sequences(annotation, fasta) -> dict[str, str]:
    """Extract one spliced, transcript-oriented CDS sequence per transcript."""
    df = _annotation_df(annotation)
    cds = df[(df["Feature"] == "CDS") & df["transcript_id"].notna()]
    sequences = {}
    for transcript_id, group in cds.groupby("transcript_id", dropna=True):
        strand = str(group["Strand"].iloc[0])
        ordered = group.sort_values("Start", ascending=strand != "-")
        parts = []
        for row in ordered.itertuples(index=False):
            try:
                parts.append(
                    _sequence(
                        fasta,
                        str(row.Chromosome),
                        int(row.Start),
                        int(row.End),
                        strand,
                    )
                )
            except KeyError:
                parts.append("")
        sequences[str(transcript_id)] = "".join(parts)
    return sequences


def compute_cds_metrics(sequence: str, genetic_code: dict = STANDARD_CODE) -> dict:
    """Measure frame, start, stop, and internal-stop properties of one CDS."""
    sequence = (sequence or "").upper()
    length = len(sequence)
    start_codon = sequence[:3] if length >= 3 else ""
    stop_codon = sequence[-3:] if length >= 3 else ""

    if length < 3:
        start_status = "missing"
        stop_status = "missing"
    else:
        if start_codon == "ATG":
            start_status = "canonical"
        elif start_codon in NEAR_COGNATE_STARTS:
            start_status = "near_cognate"
        else:
            start_status = "non_cognate"
        stop_status = (
            "canonical" if stop_codon in genetic_code["stop_codons"] else "noncanonical"
        )

    codons = [sequence[index : index + 3] for index in range(0, length - 2, 3)]
    internal_stops = sum(codon in genetic_code["stop_codons"] for codon in codons[:-1])
    return {
        "cds_len_nt": length,
        "start_codon": start_codon,
        "stop_codon": stop_codon,
        "start_status": start_status,
        "stop_status": stop_status,
        "frame_ok": length % 3 == 0,
        "frame_error": length % 3 != 0,
        "internal_stop_count": internal_stops,
        "has_internal_stop": internal_stops > 0,
    }


def _kozak_metrics(annotation: pd.DataFrame, fasta, profile: dict) -> dict[str, dict]:
    """Return transcript-level Kozak windows around the first CDS segment."""
    cds = annotation[
        (annotation["Feature"] == "CDS") & annotation["transcript_id"].notna()
    ]
    sequences = extract_cds_sequences(annotation, fasta)
    kozak = profile.get("kozak", {})
    upstream = int(kozak.get("upstream", 6))
    # The canonical -6..+4 Kozak window contains six upstream bases, the
    # three-base start codon, and one downstream base.
    downstream = int(kozak.get("downstream", 1))
    strong_minus3 = set(kozak.get("strong_minus3", "AG"))
    strong_plus4 = kozak.get("strong_plus4", "G")
    results = {}

    for transcript_id, group in cds.groupby("transcript_id", dropna=True):
        transcript_id = str(transcript_id)
        strand = str(group["Strand"].iloc[0])
        first = group.sort_values("Start", ascending=strand != "-").iloc[0]
        start = int(first.Start)
        end = int(first.End)
        if strand == "+":
            window_start = start - upstream
            window_end = start + 3 + downstream
        else:
            window_start = end - 3 - downstream
            window_end = end + upstream

        try:
            window = _sequence(
                fasta,
                str(first.Chromosome),
                max(0, window_start),
                window_end,
                strand,
            )
        except KeyError:
            window = ""

        cds_metrics = compute_cds_metrics(sequences.get(transcript_id, ""))
        interpretation_window = window.upper()
        full_window = len(window) == upstream + 3 + downstream
        minus3 = (
            interpretation_window[upstream - 3]
            if len(interpretation_window) >= upstream
            else ""
        )
        plus4 = (
            interpretation_window[upstream + 3]
            if len(interpretation_window) > upstream + 3
            else ""
        )
        if cds_metrics["start_status"] == "missing":
            context_class = "not_applicable"
        elif not full_window or "N" in interpretation_window:
            context_class = "unassessable"
        elif minus3 in strong_minus3 and plus4 == strong_plus4:
            context_class = "strong"
        elif minus3 in strong_minus3 or plus4 == strong_plus4:
            context_class = "moderate"
        else:
            context_class = "weak"

        results[transcript_id] = {
            "kozak_window": window,
            "kozak_minus3_base": minus3,
            "kozak_plus4_base": plus4,
            "kozak_context_class": context_class,
            "kozak_context_assessable": full_window
            and "N" not in interpretation_window,
        }
    return results


def compute_translation_metrics(
    annotation, fasta, genetic_code: dict = STANDARD_CODE, profile: dict | None = None
) -> pd.DataFrame:
    """Measure CDS translation events and Kozak context per coding transcript."""
    profile = profile or {}
    sequences = extract_cds_sequences(annotation, fasta)
    columns = [
        "event_type",
        "transcript_id",
        "gene_id",
        "cds_len_nt",
        "start_codon",
        "stop_codon",
        "start_status",
        "stop_status",
        "frame_ok",
        "frame_error",
        "internal_stop_count",
        "has_internal_stop",
        "kozak_window",
        "kozak_minus3_base",
        "kozak_plus4_base",
        "kozak_context_class",
        "kozak_context_assessable",
    ]
    if not sequences:
        return pd.DataFrame(columns=columns)

    df = _annotation_df(annotation)
    kozak = _kozak_metrics(df, fasta, profile)
    gene_by_transcript = (
        df[df["transcript_id"].notna()]
        .drop_duplicates("transcript_id")
        .set_index("transcript_id")["gene_id"]
        .to_dict()
    )
    rows = []
    for transcript_id, sequence in sequences.items():
        rows.append(
            {
                "event_type": "translation",
                "transcript_id": transcript_id,
                "gene_id": gene_by_transcript.get(transcript_id, ""),
                **compute_cds_metrics(sequence, genetic_code),
                **kozak.get(transcript_id, {}),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def derive_introns(annotation) -> pd.DataFrame:
    """Derive transcript-oriented intron events from exon coordinates."""
    df = _annotation_df(annotation)
    exons = df[(df["Feature"] == "exon") & df["transcript_id"].notna()]
    columns = [
        "event_type",
        "transcript_id",
        "gene_id",
        "chromosome",
        "strand",
        "intron_number",
        "intron_start",
        "intron_end",
        "intron_length_bp",
    ]
    if exons.empty:
        return pd.DataFrame(columns=columns)

    intron_parts = []
    for strand, ascending in (("-", False), ("+", True), (".", True)):
        strand_exons = exons[exons["Strand"] == strand]
        if strand_exons.empty:
            continue
        ordered = strand_exons.sort_values(
            ["transcript_id", "Start"], ascending=[True, ascending]
        )
        grouped = ordered.groupby("transcript_id", sort=False)
        previous_start = grouped["Start"].shift()
        previous_end = grouped["End"].shift()
        intron_start = np.minimum(previous_end, ordered["End"])
        intron_end = np.maximum(previous_start, ordered["Start"])
        valid = previous_start.notna() & (intron_end > intron_start)
        if not valid.any():
            continue
        part = pd.DataFrame(
            {
                "event_type": "intron",
                "transcript_id": ordered["transcript_id"].astype(str),
                "gene_id": ordered["gene_id"],
                "chromosome": ordered["Chromosome"],
                "strand": strand,
                "intron_number": grouped.cumcount() + 1,
                "intron_start": intron_start,
                "intron_end": intron_end,
                "intron_length_bp": intron_end - intron_start,
            },
            index=ordered.index,
        )
        part["intron_start"] = part["intron_start"].astype("Int64")
        part["intron_end"] = part["intron_end"].astype("Int64")
        part["intron_length_bp"] = part["intron_length_bp"].astype("Int64")
        intron_parts.append(part.loc[valid])
    if not intron_parts:
        return pd.DataFrame(columns=columns)
    return pd.concat(intron_parts, ignore_index=True)[columns]


def compute_splice_junction_metrics(
    annotation, fasta, canonical_junctions=DEFAULT_CANONICAL_JUNCTIONS
) -> pd.DataFrame:
    """Measure donor/acceptor sequence for every derived intron event."""
    introns = derive_introns(annotation)
    columns = [
        *introns.columns,
        "donor_dinucleotide",
        "acceptor_dinucleotide",
        "junction_type",
        "is_canonical",
        "sequence_assessable",
    ]
    if introns.empty:
        return pd.DataFrame(columns=columns)

    rows = []
    for row in introns.itertuples(index=False):
        try:
            if row.strand == "-":
                donor = _sequence(
                    fasta, row.chromosome, row.intron_end - 2, row.intron_end, "-"
                )
                acceptor = _sequence(
                    fasta, row.chromosome, row.intron_start, row.intron_start + 2, "-"
                )
            else:
                donor = _sequence(
                    fasta, row.chromosome, row.intron_start, row.intron_start + 2, "+"
                )
                acceptor = _sequence(
                    fasta, row.chromosome, row.intron_end - 2, row.intron_end, "+"
                )
            assessable = len(donor) == 2 and len(acceptor) == 2
        except KeyError:
            donor, acceptor, assessable = "", "", False

        junction_type = (
            f"{donor.upper()}-{acceptor.upper()}" if assessable else "unknown"
        )
        rows.append(
            {
                **row._asdict(),
                "donor_dinucleotide": donor.upper(),
                "acceptor_dinucleotide": acceptor.upper(),
                "junction_type": junction_type,
                "is_canonical": junction_type in canonical_junctions,
                "sequence_assessable": assessable,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _scan_boundary_motifs(
    annotation, fasta, profile: dict, boundary: str
) -> pd.DataFrame:
    """Scan transcript-oriented 5' or 3' boundary windows for configured motifs."""
    df = _annotation_df(annotation)
    transcripts = _transcript_boundaries(df)
    section = profile.get(boundary, {})
    motifs = section.get("motifs", {})
    upstream = int(section.get("upstream", 0))
    columns = [
        "event_type",
        "transcript_id",
        "gene_id",
        "chromosome",
        "strand",
        "boundary",
        "window_sequence",
        "motif_hit_count",
        "nearest_motif_distance_bp",
        "nearest_motif",
        "sequence_assessable",
        "status",
    ]
    if transcripts.empty:
        return pd.DataFrame(columns=columns)

    rows = []
    for row in transcripts.itertuples(index=False):
        if boundary == "promoter":
            boundary_coordinate = int(row.Start) if row.Strand != "-" else int(row.End)
        else:
            boundary_coordinate = int(row.End) if row.Strand != "-" else int(row.Start)

        if row.Strand == "-":
            start, end = boundary_coordinate, boundary_coordinate + upstream
        else:
            start, end = boundary_coordinate - upstream, boundary_coordinate
        try:
            sequence = _sequence(fasta, row.Chromosome, max(0, start), end, row.Strand)
            assessable = len(sequence) == upstream and "N" not in sequence.upper()
        except KeyError:
            sequence, assessable = "", False

        hits = []
        for name, pattern in motifs.items():
            for match in re.finditer(pattern, sequence, flags=re.IGNORECASE):
                hits.append((len(sequence) - match.end(), name))
        hits.sort()
        status = (
            "not_configured"
            if not motifs
            else (
                "unassessable"
                if not assessable
                else "sequence_supported_candidate" if hits else "no_candidate"
            )
        )
        rows.append(
            {
                "event_type": f"{boundary}_boundary",
                "transcript_id": row.transcript_id,
                "gene_id": row.gene_id,
                "chromosome": row.Chromosome,
                "strand": row.Strand,
                "boundary": boundary,
                "window_sequence": sequence,
                "motif_hit_count": len(hits),
                "nearest_motif_distance_bp": hits[0][0] if hits else pd.NA,
                "nearest_motif": hits[0][1] if hits else pd.NA,
                "sequence_assessable": assessable,
                "status": status,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def compute_promoter_proximity_metrics(annotation, fasta, profile=None) -> pd.DataFrame:
    """Scan upstream 5' transcript windows for configured promoter-like motifs."""
    return _scan_boundary_motifs(annotation, fasta, profile or {}, "promoter")


def compute_polyA_motif_metrics(annotation, fasta, profile=None) -> pd.DataFrame:
    """Scan upstream 3' transcript windows for configured poly-A motifs."""
    profile = profile or {}
    profile = dict(profile)
    profile.setdefault("polya", {"motifs": DEFAULT_POLYA_MOTIFS, "upstream": 50})
    return _scan_boundary_motifs(annotation, fasta, profile, "polya")

"""
Optional stop-codon harmonisation of a CDS annotation (stop-excluded convention).

Some annotations end the CDS before the stop codon; others, including Ensembl,
include it. Compared as written, every complete model then differs by 3 bp and
is never coordinate-exact. harmonise_stop_codons extends the terminal CDS
segment over the stop codon that immediately follows it, when the model's
reading frame is established from the sequence itself. It aligns conventions; it
does not validate ORFs and repairs nothing else. Inputs and outputs are
comparison-schema DataFrames; sequences are passed in, nothing reads files here.

Genetic code. Each sequence gets an NCBI translation table: the default table
for every sequence, unless a per-sequence override gives another table or
"none" (skip the sequence). Sequences whose name says they are organelle genomes
(ORGANELLE_SEQNAMES, case-insensitive: MT, chrM, mitochondrion, Pt, chloroplast,
...) are not given the default table: without an override their transcripts are
skipped (genetic_code_not_specified).

Eligibility (checked in this order; the first failure is the reason). The frame
is taken to start at the first CDS base, and is accepted only with sequence
evidence for it. CDS phase attributes are not used, so missing phase does not
matter; a 5'-partial model (frame not starting at its first base) fails the
start-codon or length check and is skipped rather than guessed.

    sequence_not_available        the sequence is not in the genome
    genetic_code_not_specified    organelle sequence without an override, or an
                                  override of "none"
    strand_not_defined            strand is not + or -
    overlapping_cds_segments      CDS segments of the transcript overlap
    cds_outside_sequence          a CDS segment lies beyond the sequence
    ambiguous_bases_in_cds        a CDS base is not A, C, G or T
    cds_length_not_multiple_of_3  disrupted or partial frame
    no_start_codon                first codon is not ATG (5'-partial or
                                  non-canonical start: frame not established)
    internal_stop_codon           an in-frame stop before the last codon

Transcripts passing these checks are "eligible". For them:

    unchanged  stop_already_in_cds     last codon is a stop: nothing to do (this
                                       also makes a second application a no-op)
    skipped    split_stop_codon_unsupported
                                       the terminal CDS segment ends within 3 bp
                                       of its exon end and another exon follows:
                                       a stop split by an intron is not added
    skipped    extension_beyond_explicit_exons
                                       explicit exons end before the extension;
                                       exon records are never rewritten
    skipped    cds_outside_explicit_exons
                                       the terminal CDS segment is not inside an
                                       exon record
    skipped    extension_out_of_bounds the 3 bases would pass the sequence end
    skipped    ambiguous_bases_after_cds
                                       the next 3 bases contain a non-ACGT base
    unchanged  no_adjacent_stop        the next 3 bases are not a stop codon
    skipped    extension_outside_transcript_record /
               extension_outside_gene_record
                                       an explicit transcript or gene record does
                                       not contain the extended CDS
    extended   stop_codon_added        terminal CDS segment extended by 3 bp

Only the terminal CDS segment changes. Exons inferred from CDS (CDS-only files,
source_feature inferred_exon_from_CDS) are extended with it, as are transcript
and gene records inferred from their children. Identifiers, parents, other
segments and explicit records are never changed.
"""

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field

import pandas as pd

# NCBI translation tables with unambiguous stop codons. Tables whose stops are
# context-dependent (27, 28, 31) are not supported.
GENETIC_CODE_STOPS = {
    1: ("TAA", "TAG", "TGA"),
    2: ("TAA", "TAG", "AGA", "AGG"),
    3: ("TAA", "TAG"),
    4: ("TAA", "TAG"),
    5: ("TAA", "TAG"),
    6: ("TGA",),
    9: ("TAA", "TAG"),
    10: ("TAA", "TAG"),
    11: ("TAA", "TAG", "TGA"),
    12: ("TAA", "TAG", "TGA"),
    13: ("TAA", "TAG"),
    14: ("TAG",),
    16: ("TAA", "TGA"),
    21: ("TAA", "TAG"),
    22: ("TCA", "TAA", "TGA"),
    23: ("TTA", "TAA", "TAG", "TGA"),
    24: ("TAA", "TAG"),
    25: ("TAA", "TAG"),
    26: ("TAA", "TAG", "TGA"),
    29: ("TGA",),
    30: ("TGA",),
    33: ("TAG",),
}
START_CODON = "ATG"
ORGANELLE_SEQNAMES = frozenset(
    {
        "mt",
        "m",
        "chrm",
        "chrmt",
        "mito",
        "mitochondrion",
        "mitochondrial",
        "pt",
        "chrc",
        "chrpt",
        "chloroplast",
        "plastid",
    }
)
EXTENDED, UNCHANGED, SKIPPED = "extended", "unchanged", "skipped"
FRAME_CHECKS = (
    "sequence_not_available",
    "genetic_code_not_specified",
    "strand_not_defined",
    "overlapping_cds_segments",
    "cds_outside_sequence",
    "ambiguous_bases_in_cds",
    "cds_length_not_multiple_of_3",
    "no_start_codon",
    "internal_stop_codon",
)
AUDIT_COLUMNS = [
    "transcript_id",
    "original_transcript_id",
    "gene_id",
    "original_gene_id",
    "chrom",
    "strand",
    "genetic_code",
    "genetic_code_source",
    "cds_segments",
    "cds_length",
    "eligible",
    "status",
    "reason",
    "last_codon",
    "next_codon",
    "exon_representation",
    "terminal_cds_start",
    "terminal_cds_end",
    "new_terminal_cds_start",
    "new_terminal_cds_end",
    "records_changed",
]
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_BASES = frozenset("ACGT")


@dataclass
class GeneticCodePolicy:
    """NCBI translation table per sequence (None: do not harmonise)."""

    default: int = 1
    overrides: dict = field(default_factory=dict)

    def code_for(self, seqname: str) -> tuple[int | None, str]:
        """Return (table or None, how it was chosen)."""
        if seqname in self.overrides:
            return self.overrides[seqname], "sequence_override"
        if seqname.lower() in ORGANELLE_SEQNAMES:
            return None, "organelle_sequence_without_override"
        return self.default, "default"

    def describe(self) -> dict:
        """Policy as a JSON-ready dict (tables, overrides, stop codons used)."""
        return {
            "default_table": self.default,
            "sequence_overrides": dict(self.overrides),
            "organelle_seqnames_needing_override": sorted(ORGANELLE_SEQNAMES),
            "stop_codons": {
                str(code): list(GENETIC_CODE_STOPS[code])
                for code in sorted(
                    {self.default}
                    | {c for c in self.overrides.values() if c is not None}
                )
            },
        }


def parse_genetic_code_policy(
    default: int, overrides: list[str] | None
) -> GeneticCodePolicy:
    """
    Build a policy from a default table and NAME=TABLE overrides (TABLE may be none).
    Raises:
            ValueError for an unsupported table or malformed override
    """

    def table(value) -> int | None:
        if str(value).lower() == "none":
            return None
        try:
            code = int(value)
        except ValueError as error:
            raise ValueError(f"Genetic code '{value}' is not a number") from error
        if code not in GENETIC_CODE_STOPS:
            raise ValueError(
                f"Unsupported genetic code {code}; supported NCBI tables: "
                f"{sorted(GENETIC_CODE_STOPS)}"
            )
        return code

    default_code = table(default)
    if default_code is None:
        raise ValueError("The default genetic code cannot be none")
    parsed = {}
    for item in overrides or []:
        name, sep, value = item.partition("=")
        if not sep or not name.strip():
            raise ValueError(f"Expected SEQNAME=TABLE, got '{item}'")
        parsed[name.strip()] = table(value.strip())
    return GeneticCodePolicy(default_code, parsed)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass
class _Model:  # pylint: disable=too-many-instance-attributes
    """One query transcript with CDS (zero-based half-open coordinates)."""

    transcript_id: str
    chrom: str
    strand: str
    cds: list  # [(start, end, row index)], sorted by start
    exons: list  # [(start, end, row index)], sorted by start
    exons_inferred: bool
    transcript_row: tuple | None  # (start, end, row index, inferred)
    gene_row: tuple | None


def _frame_reason(  # pylint: disable=too-many-return-statements
    model: _Model, sequence: str | None, stops
) -> tuple[str, dict]:
    """First failed frame check (FRAME_CHECKS order), or "" when the frame is established."""
    info = {"last_codon": "", "cds_length": sum(e - s for s, e, _ in model.cds)}
    if sequence is None:
        return "sequence_not_available", info
    if stops is None:
        return "genetic_code_not_specified", info
    if model.strand not in ("+", "-"):
        return "strand_not_defined", info
    if any(model.cds[i + 1][0] < model.cds[i][1] for i in range(len(model.cds) - 1)):
        return "overlapping_cds_segments", info
    if model.cds[0][0] < 0 or model.cds[-1][1] > len(sequence):
        return "cds_outside_sequence", info
    nucleotides = "".join(sequence[s:e] for s, e, _ in model.cds)
    if model.strand == "-":
        nucleotides = _reverse_complement(nucleotides)
    if not set(nucleotides) <= _BASES:
        return "ambiguous_bases_in_cds", info
    if len(nucleotides) % 3 or not nucleotides:
        return "cds_length_not_multiple_of_3", info
    codons = [nucleotides[i : i + 3] for i in range(0, len(nucleotides), 3)]
    info["last_codon"] = codons[-1]
    if codons[0] != START_CODON:
        return "no_start_codon", info
    if any(codon in stops for codon in codons[:-1]):
        return "internal_stop_codon", info
    return "", info


def _exon_reason(model: _Model, new_start: int, new_end: int) -> str:
    """Why explicit exons cannot hold the extended terminal CDS segment, or ""."""
    if model.exons_inferred:
        return ""
    terminal = model.cds[-1] if model.strand == "+" else model.cds[0]
    holder = [e for e in model.exons if e[0] <= terminal[0] and terminal[1] <= e[1]]
    if not holder:
        return "cds_outside_explicit_exons"
    exon_start, exon_end, _ = holder[0]
    if exon_start <= new_start and new_end <= exon_end:
        return ""
    downstream = (
        any(e[0] >= exon_end for e in model.exons)
        if model.strand == "+"
        else any(e[1] <= exon_start for e in model.exons)
    )
    return (
        "split_stop_codon_unsupported"
        if downstream
        else "extension_beyond_explicit_exons"
    )


def _parent_reason(model: _Model, new_start: int, new_end: int) -> str:
    for record, reason in (
        (model.transcript_row, "extension_outside_transcript_record"),
        (model.gene_row, "extension_outside_gene_record"),
    ):
        if record is None:
            continue
        start, end, _, inferred = record
        if not inferred and (new_start < start or new_end > end):
            return reason
    return ""


def _extension(model: _Model, sequence: str, stops, terminal) -> dict:
    """Outcome for an eligible model: status, reason, next codon, new coordinates."""
    if model.strand == "+":
        new_start, new_end = terminal[0], terminal[1] + 3
    else:
        new_start, new_end = terminal[0] - 3, terminal[1]
    outcome = {"status": SKIPPED, "next_codon": "", "new": (new_start, new_end)}
    reason = _exon_reason(model, new_start, new_end)
    if not reason and (new_start < 0 or new_end > len(sequence)):
        reason = "extension_out_of_bounds"
    if not reason:
        following = (
            sequence[terminal[1] : new_end]
            if model.strand == "+"
            else _reverse_complement(sequence[new_start : terminal[0]])
        )
        outcome["next_codon"] = following
        if not set(following) <= _BASES:
            reason = "ambiguous_bases_after_cds"
        elif following not in stops:
            outcome["status"], reason = UNCHANGED, "no_adjacent_stop"
        else:
            reason = _parent_reason(model, new_start, new_end)
            if not reason:
                outcome["status"], reason = EXTENDED, "stop_codon_added"
    outcome["reason"] = reason
    return outcome


def _edits(model: _Model, terminal, new_start: int, new_end: int) -> tuple:
    """Row edits for an extension: terminal CDS, inferred exon and inferred spans."""
    edits = {terminal[2]: (new_start, new_end)}
    changed = ["CDS"]
    if model.exons_inferred:
        for start, end, row in model.exons:
            if (start, end) == (terminal[0], terminal[1]):
                edits[row] = (new_start, new_end)
                changed.append("inferred_exon")
                break
    for record, label in (
        (model.transcript_row, "transcript"),
        (model.gene_row, "gene"),
    ):
        if record is None:
            continue
        start, end, row, inferred = record
        if inferred and (new_start < start or new_end > end):
            edits[row] = (min(start, new_start), max(end, new_end))
            changed.append(f"inferred_{label}_span")
    return edits, changed


def _evaluate(model: _Model, sequence: str | None, policy: GeneticCodePolicy):
    """Return (audit fields, {row index: (start, end)} edits)."""
    code, code_source = policy.code_for(model.chrom)
    stops = GENETIC_CODE_STOPS[code] if code is not None else None
    reason, info = _frame_reason(model, sequence, stops)
    terminal = model.cds[-1] if model.strand == "+" else model.cds[0]
    audit = {
        "genetic_code": code if code is not None else "",
        "genetic_code_source": code_source,
        "cds_segments": len(model.cds),
        "cds_length": info["cds_length"],
        "eligible": not reason,
        "status": SKIPPED,
        "reason": reason,
        "last_codon": info["last_codon"],
        "next_codon": "",
        "exon_representation": (
            "inferred_from_cds" if model.exons_inferred else "explicit"
        ),
        "terminal_cds_start": terminal[0],
        "terminal_cds_end": terminal[1],
        "new_terminal_cds_start": terminal[0],
        "new_terminal_cds_end": terminal[1],
        "records_changed": "",
    }
    if reason:
        return audit, {}
    if info["last_codon"] in stops:
        audit.update(status=UNCHANGED, reason="stop_already_in_cds")
        return audit, {}
    outcome = _extension(model, sequence, stops, terminal)
    audit.update(
        status=outcome["status"],
        reason=outcome["reason"],
        next_codon=outcome["next_codon"],
    )
    if outcome["status"] != EXTENDED:
        return audit, {}
    edits, changed = _edits(model, terminal, *outcome["new"])
    audit.update(
        new_terminal_cds_start=outcome["new"][0],
        new_terminal_cds_end=outcome["new"][1],
        records_changed=",".join(changed),
    )
    return audit, edits


def _records(rows: pd.DataFrame, key: str) -> dict:
    """key -> (start, end, row index, inferred) for gene or transcript rows."""
    return {
        k: (s, e, i, f.startswith("inferred_"))
        for k, s, e, i, f in zip(
            rows[key], rows["Start"], rows["End"], rows.index, rows["source_feature"]
        )
    }


def _exons_by_transcript(exons: pd.DataFrame) -> dict:
    """transcript_id -> sorted [(start, end, row index, inferred from CDS)]."""
    by_tx: dict = {}
    for tid, start, end, row, source in zip(
        exons["transcript_id"],
        exons["Start"],
        exons["End"],
        exons.index,
        exons["source_feature"],
    ):
        by_tx.setdefault(tid, []).append(
            (start, end, row, source == "inferred_exon_from_CDS")
        )
    return {tid: sorted(rows) for tid, rows in by_tx.items()}


def _models(df: pd.DataFrame) -> dict[str, list[_Model]]:
    """Query transcripts with CDS, grouped by sequence."""
    tx_records = _records(
        df[df["Feature"] == "transcript"].drop_duplicates("transcript_id"),
        "transcript_id",
    )
    gene_records = _records(
        df[df["Feature"] == "gene"].drop_duplicates("gene_id"), "gene_id"
    )
    exons_by_tx = _exons_by_transcript(df[df["Feature"] == "exon"])
    by_chrom: dict[str, list[_Model]] = {}
    for tid, group in df[df["Feature"] == "CDS"].groupby("transcript_id", sort=False):
        tx_exons = exons_by_tx.get(tid, [])
        strands = set(group["Strand"])
        model = _Model(
            transcript_id=tid,
            chrom=str(group["Chromosome"].iloc[0]),
            strand=strands.pop() if len(strands) == 1 else ".",
            cds=sorted(zip(group["Start"], group["End"], group.index)),
            exons=[(s, e, i) for s, e, i, _ in tx_exons],
            exons_inferred=bool(tx_exons) and all(inf for *_, inf in tx_exons),
            transcript_row=tx_records.get(tid),
            gene_row=gene_records.get(group["gene_id"].iloc[0]),
        )
        by_chrom.setdefault(model.chrom, []).append(model)
    return by_chrom


def _audit_table(records: list, df: pd.DataFrame) -> pd.DataFrame:
    """One AUDIT_COLUMNS row per evaluated transcript, sorted by position."""
    cds = df[df["Feature"] == "CDS"].drop_duplicates("transcript_id")
    names = cds.set_index("transcript_id")[
        ["original_transcript_id", "gene_id", "original_gene_id"]
    ].to_dict("index")
    rows = [
        {
            "transcript_id": model.transcript_id,
            **names[model.transcript_id],
            "chrom": model.chrom,
            "strand": model.strand,
            **audit,
        }
        for model, audit in records
    ]
    table = pd.DataFrame(rows, columns=AUDIT_COLUMNS)
    if len(table):
        table = table.sort_values(
            ["chrom", "terminal_cds_start", "transcript_id"], kind="stable"
        ).reset_index(drop=True)
    return table


def _merge_edits(edits: dict, new: dict) -> None:
    """Add row edits; several transcripts may widen the same inferred gene span."""
    for row, (start, end) in new.items():
        if row in edits:
            start, end = min(start, edits[row][0]), max(end, edits[row][1])
        edits[row] = (start, end)


def _apply_edits(df: pd.DataFrame, edits: dict) -> pd.DataFrame:
    out = df.copy()
    for row, (start, end) in edits.items():
        out.at[row, "Start"], out.at[row, "End"] = start, end
    return out


def harmonise_stop_codons(
    df: pd.DataFrame,
    sequences: Iterable[tuple[str, str]],
    policy: GeneticCodePolicy,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Extend terminal CDS segments over an adjacent stop codon (see module docstring).
    Args:
            df: Comparison-schema DataFrame (the query annotation)
            sequences: (name, upper-case sequence) for the sequences df uses; may
                    yield other names, which are ignored
            policy: genetic-code policy
    Returns:
            (harmonised DataFrame, audit DataFrame with AUDIT_COLUMNS in zero-based
            coordinates, one row per transcript with CDS, summary dict)
    """
    by_chrom = _models(df)
    records, edits = [], {}
    for name, sequence in sequences:
        for model in by_chrom.pop(name, []):
            audit, model_edits = _evaluate(model, sequence, policy)
            records.append((model, audit))
            _merge_edits(edits, model_edits)
    for models in by_chrom.values():  # sequences the genome did not provide
        records += [(model, _evaluate(model, None, policy)[0]) for model in models]
    audit_table = _audit_table(records, df)
    return (
        _apply_edits(df, edits),
        audit_table,
        summarise_harmonisation(audit_table, policy),
    )


def summarise_harmonisation(audit: pd.DataFrame, policy: GeneticCodePolicy) -> dict:
    """Counts by status and reason, plus the policy."""
    status = Counter(audit["status"])
    return {
        "applied": True,
        "annotation": "query",
        "transcripts_with_cds": len(audit),
        "eligible": int(audit["eligible"].sum()) if len(audit) else 0,
        "extended": status.get(EXTENDED, 0),
        "unchanged": status.get(UNCHANGED, 0),
        "skipped": status.get(SKIPPED, 0),
        "reasons": {
            f"{s}:{r}": int(n)
            for (s, r), n in sorted(
                Counter(zip(audit["status"], audit["reason"])).items()
            )
        },
        "genetic_code_by_sequence_source": {
            k: int(v) for k, v in sorted(Counter(audit["genetic_code_source"]).items())
        },
        "policy": {
            "genetic_code": policy.describe(),
            "frame_checks_in_order": list(FRAME_CHECKS),
            "start_codon": START_CODON,
            "phase_used": False,
            "split_stop_codons": "not added (skipped)",
            "explicit_exon_records": "never changed; incompatible extensions skipped",
            "inferred_records": (
                "exons inferred from CDS and inferred transcript/gene spans are "
                "extended with the CDS"
            ),
        },
    }

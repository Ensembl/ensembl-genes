"""
Splice-site and intron recovery along reference CDS structures (exploratory).

Tool-agnostic: works on the stored evaluated models of the reference and of one
prediction run under one transcript policy (one-based inclusive CDS segments, as
selected by the comparator's own selection functions). Nothing here reads
annotation files or knows which predictor produced the models.

Definitions

    CDS intron      the gap between consecutive CDS segments of one transcript
                    (one-based inclusive); a single-segment CDS has none. UTR
                    exons are never used.
    donor/acceptor  the 5' / 3' end of a CDS intron in the direction of
                    transcription; on the minus strand the donor is the intron's
                    higher coordinate. A splice site is identified by (sequence,
                    strand, kind, intron-side coordinate).
    support         exact, same sequence and same strand. A reference intron is
                    recovered when an evaluated query transcript has the identical
                    CDS intron; a donor (acceptor) is supported when an evaluated
                    query transcript has a CDS intron with the same donor
                    (acceptor). The query transcript need not belong to the
                    prediction matched to the gene: this is the comparator's
                    unique-intron rule (``intron_support``) extended to single sites.
    detected locus  the comparator recovered the reference gene's locus
                    (same-strand CDS overlap, ``cds_overlap > 0``). Sites of genes
                    whose locus was not recovered are counted separately as
                    missing-locus sites and are excluded from the recovery rates,
                    so the rates describe structure at loci the predictor found.

Binned results need one reference transcript per gene (a policy whose reference
selection is not "all"), so every gene contributes each of its introns once.
Under "all" only the unique-intron totals are computed, deduplicated by
sequence, strand and coordinates exactly as the comparator does.

Position along the reference CDS, in the direction of transcription (0 = first
coding base, 1 = last coding base):

    cds_span  genomic distance from the coding start to the exon-intron junction,
              over the genomic CDS span (introns included)
    coding    coding bases upstream of the junction over the CDS length

A donor junction lies after the last base of its upstream segment, an acceptor
junction before the first base of its downstream segment, so both kinds always
fall strictly inside (0, 1).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping

POSITION_BINS = 10
POSITION_COORDINATES = ("cds_span", "coding")
SITE_KINDS = ("donor", "acceptor")
# Lower edges (bp) of the reference CDS intron-length bins; the last bin is open.
# Fixed, so every run and scope uses the same bins. The first two bins keep gaps
# of a few bases (often frameshift or indel gaps in a CDS rather than spliced
# introns) apart from short introns.
INTRON_LENGTH_EDGES = (
    0,
    10,
    30,
    60,
    100,
    200,
    500,
    1_000,
    2_000,
    5_000,
    10_000,
    20_000,
    50_000,
    100_000,
)
# Position bins are also reported without reference CDS introns shorter than
# this (an intron-length bin edge): such gaps are rarely modelled by predictors
# and would otherwise dominate the unsupported sites in every position bin.
POSITION_MIN_LENGTHS = (0, 60)
assert set(POSITION_MIN_LENGTHS) <= set(INTRON_LENGTH_EDGES)
SUPPORT_RULE = (
    "exact same-sequence, same-strand support by any evaluated query transcript's CDS intron "
    "(identical intron for intron recovery; identical donor or acceptor coordinate for a splice site)"
)
DETECTION_RULE = (
    "reference gene locus recovered by the comparator (same-strand CDS overlap, cds_overlap > 0); "
    "sites of other genes are counted as missing-locus sites and excluded from rates"
)

Segments = tuple[tuple[int, int], ...]


def parse_segments(text: str | None) -> Segments:
    """'start-end,...' (one-based inclusive) -> sorted segments."""
    if not text:
        return ()
    return tuple(
        sorted(tuple(int(x) for x in part.split("-")) for part in text.split(",") if part)
    )


def cds_introns(segments: Segments) -> list[tuple[int, int]]:
    """One-based inclusive (lo, hi) of each CDS intron, in coordinate order."""
    return [(segments[i][1] + 1, segments[i + 1][0] - 1) for i in range(len(segments) - 1)]


def site_coordinates(strand: str, lo: int, hi: int) -> tuple[int, int]:
    """(donor, acceptor) intron-side coordinates of one intron."""
    return (hi, lo) if strand == "-" else (lo, hi)


def length_bin(length: int) -> int:
    index = 0
    for i, edge in enumerate(INTRON_LENGTH_EDGES):
        if length >= edge:
            index = i
    return index


def length_bin_label(index: int) -> str:
    lo = INTRON_LENGTH_EDGES[index]
    if index + 1 == len(INTRON_LENGTH_EDGES):
        return f"≥ {_bp(lo)}"
    hi = INTRON_LENGTH_EDGES[index + 1]
    return f"< {_bp(hi)}" if lo == 0 else f"{_bp(lo)}–{_bp(hi)}"


def _bp(n: int) -> str:
    return f"{n // 1000} kb" if n >= 1000 and n % 1000 == 0 else f"{n} bp"


def position_bin(position: float) -> int:
    return min(int(position * POSITION_BINS), POSITION_BINS - 1)


@dataclass
class Intron:
    lo: int
    hi: int
    donor: int
    acceptor: int
    length: int
    position: dict  # coordinate -> {kind: fraction}


def intron_sites(strand: str, segments: Segments) -> list[Intron]:
    """CDS introns of one transcript with donor/acceptor and normalised positions."""
    if len(segments) < 2:
        return []
    cds_lo, cds_hi = segments[0][0], segments[-1][1]
    span = cds_hi - cds_lo + 1
    lengths = [e - s + 1 for s, e in segments]
    total = sum(lengths)
    out = []
    for i in range(len(segments) - 1):
        lo, hi = segments[i][1] + 1, segments[i + 1][0] - 1
        donor, acceptor = site_coordinates(strand, lo, hi)
        if strand == "-":
            donor_offset = cds_hi - segments[i + 1][0] + 1
            acceptor_offset = cds_hi - segments[i][1]
            coding = sum(lengths[i + 1 :])
        else:
            donor_offset = segments[i][1] - cds_lo + 1
            acceptor_offset = segments[i + 1][0] - cds_lo
            coding = sum(lengths[: i + 1])
        out.append(
            Intron(
                lo=lo,
                hi=hi,
                donor=donor,
                acceptor=acceptor,
                length=hi - lo + 1,
                position={
                    "cds_span": {
                        "donor": donor_offset / span,
                        "acceptor": acceptor_offset / span,
                    },
                    "coding": {"donor": coding / total, "acceptor": coding / total},
                },
            )
        )
    return out


@dataclass
class SupportIndex:
    """Unique CDS introns and splice sites of the evaluated query transcripts."""

    introns: set = field(default_factory=set)  # (chrom, strand, lo, hi)
    donors: set = field(default_factory=set)  # (chrom, strand, coordinate)
    acceptors: set = field(default_factory=set)

    @classmethod
    def build(cls, transcripts: Iterable[tuple[str, str, Segments]]) -> "SupportIndex":
        index = cls()
        for chrom, strand, segments in transcripts:
            for lo, hi in cds_introns(segments):
                donor, acceptor = site_coordinates(strand, lo, hi)
                index.introns.add((chrom, strand, lo, hi))
                index.donors.add((chrom, strand, donor))
                index.acceptors.add((chrom, strand, acceptor))
        return index


@dataclass
class ReferenceGene:
    gene_id: str
    chrom: str
    strand: str
    transcripts: list  # evaluated transcripts' CDS segments under the policy


def _empty_scope(binned: bool) -> dict:
    out = {"unique": {"reference": set(), "query": 0}, "genes": 0, "genes_with_introns": 0,
           "missing_locus_genes_with_introns": 0}
    if binned:
        out["position"] = {
            (m, c, k, b): [0, 0, 0]
            for m in POSITION_MIN_LENGTHS
            for c in POSITION_COORDINATES
            for k in SITE_KINDS
            for b in range(POSITION_BINS)
        }
        out["length"] = {b: [0, 0, 0, 0, 0] for b in range(len(INTRON_LENGTH_EDGES))}
    return out


def analyse(
    genes: Iterable[ReferenceGene],
    query: SupportIndex,
    detected: Mapping[str, bool],
    scopes: Mapping[str, set | None],
    one_per_gene: bool,
) -> dict:
    """
    Splice-site and intron recovery of one result, for each sequence scope.

    Args:
            genes: evaluated reference genes (all compared sequences)
            query: support index of the evaluated query transcripts
            detected: reference gene_id -> locus recovered (every gene must be present)
            scopes: scope id -> set of sequence names (None = every compared sequence);
                    applied to reference genes and to query introns alike
            one_per_gene: True when the policy selects one reference transcript per
                    gene; required for the binned results
    Returns:
            scope id -> {"reference_introns", "query_introns", "shared_introns",
            "genes", "genes_with_introns", "missing_locus_genes_with_introns",
            and, when one_per_gene, "position" rows and "length" rows}
    Raises:
            ValueError when a gene has no detection status, or has several
            transcripts although one_per_gene is set
    """
    acc = {sid: _empty_scope(one_per_gene) for sid in scopes}
    in_scope = {
        sid: (lambda chrom, seqs=seqs: seqs is None or chrom in seqs)
        for sid, seqs in scopes.items()
    }
    for gene in genes:
        if gene.gene_id not in detected:
            raise ValueError(f"reference gene {gene.gene_id} has no locus outcome")
        targets = [sid for sid in scopes if in_scope[sid](gene.chrom)]
        if not targets:
            continue
        if one_per_gene and len(gene.transcripts) > 1:
            raise ValueError(
                f"reference gene {gene.gene_id} has {len(gene.transcripts)} evaluated transcripts; "
                "binned results need one transcript per gene"
            )
        found = bool(detected[gene.gene_id])
        introns_by_tx = [intron_sites(gene.strand, segs) for segs in gene.transcripts]
        has_introns = any(introns_by_tx)
        for sid in targets:
            a = acc[sid]
            a["genes"] += 1
            a["genes_with_introns"] += has_introns
            a["missing_locus_genes_with_introns"] += has_introns and not found
            for introns in introns_by_tx:
                for intron in introns:
                    a["unique"]["reference"].add((gene.chrom, gene.strand, intron.lo, intron.hi))
        if not one_per_gene:
            continue
        for intron in introns_by_tx[0] if introns_by_tx else []:
            key = (gene.chrom, gene.strand)
            exact = (*key, intron.lo, intron.hi) in query.introns
            donor = (*key, intron.donor) in query.donors
            acceptor = (*key, intron.acceptor) in query.acceptors
            lbin = length_bin(intron.length)
            for sid in targets:
                a = acc[sid]
                row = a["length"][lbin]
                if found:
                    row[0] += 1
                    row[1] += exact
                    row[2] += donor
                    row[3] += acceptor
                else:
                    row[4] += 1
                for min_length in POSITION_MIN_LENGTHS:
                    if intron.length < min_length:
                        continue
                    for coordinate in POSITION_COORDINATES:
                        for kind, supported in (("donor", donor), ("acceptor", acceptor)):
                            cell = a["position"][
                                (
                                    min_length,
                                    coordinate,
                                    kind,
                                    position_bin(intron.position[coordinate][kind]),
                                )
                            ]
                            if found:
                                cell[0] += 1
                                cell[1] += supported
                            else:
                                cell[2] += 1
    out = {}
    for sid, seqs in scopes.items():
        a = acc[sid]
        ref_introns = a["unique"]["reference"]
        query_introns = {i for i in query.introns if seqs is None or i[0] in seqs}
        result = {
            "reference_introns": len(ref_introns),
            "query_introns": len(query_introns),
            "shared_introns": len(ref_introns & query_introns),
            "genes": a["genes"],
            "genes_with_introns": a["genes_with_introns"],
            "missing_locus_genes_with_introns": a["missing_locus_genes_with_introns"],
        }
        if one_per_gene:
            result["position"] = [
                {
                    "min_intron_length": m,
                    "coordinate": c,
                    "kind": k,
                    "bin": b,
                    "lo": b / POSITION_BINS,
                    "hi": (b + 1) / POSITION_BINS,
                    "assessed": v[0],
                    "supported": v[1],
                    "missing_locus": v[2],
                }
                for (m, c, k, b), v in a["position"].items()
            ]
            result["length"] = [
                {
                    "bin": b,
                    "lo": INTRON_LENGTH_EDGES[b],
                    "hi": INTRON_LENGTH_EDGES[b + 1] if b + 1 < len(INTRON_LENGTH_EDGES) else None,
                    "label": length_bin_label(b),
                    "eligible": v[0],
                    "exact": v[1],
                    "donor": v[2],
                    "acceptor": v[3],
                    "missing_locus": v[4],
                }
                for b, v in a["length"].items()
            ]
        out[sid] = result
    return out

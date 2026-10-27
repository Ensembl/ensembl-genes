"""
Write the SYNTHETIC handover example inputs (deterministic; run from this folder).

    python make_inputs.py

inputs/genome.fa         chr 1 (3,000 bp) and chr 2 (1,000 bp); every reference CDS is a real
                         ORF (ATG ... TAA, no internal stop), so --add-stop-codon has an effect
inputs/reference.gff3    R1 (+, 3 CDS exons), R2 (-, 2 CDS exons), R3 (+, 1 exon) on chr 1;
                         R4 (+, 1 exon) on chr 2; Ensembl-style biotypes
inputs/focus.gff3        "Internal predictions": R1 with the stop codon left out of the CDS
                         (exact only after --add-stop-codon), R2 and R4 exact, R3 missed
inputs/other.gff3        a second predictor: R1 exact, R2 with a different donor, R3 exact
"""

import random
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")
STOPS = {"TAA", "TAG", "TGA"}
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def orf(length: int, rng: random.Random) -> str:
    assert length % 3 == 0
    return "ATG" + "".join(rng.choice(CODONS) for _ in range(length // 3 - 2)) + "TAA"


def place(seq: list, segments, coding: str, strand: str) -> None:
    """Write the spliced coding sequence into one-based inclusive genomic segments."""
    if strand == "-":
        coding = coding.translate(COMP)[::-1]
    i = 0
    for start, end in segments:
        n = end - start + 1
        seq[start - 1 : end] = list(coding[i : i + n])
        i += n
    for (_, a_end), (b_start, _) in zip(segments, segments[1:]):  # GT...AG introns
        donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
        seq[a_end : a_end + 2] = list(donor)
        seq[b_start - 3 : b_start - 1] = list(acceptor)


def gene(seqid, gid, strand, transcripts, source="example"):
    """transcripts: {tid: (exons, cds)}; Ensembl-style gene/mRNA/exon/CDS rows."""
    starts = [s for ex, _ in transcripts.values() for s, _ in ex]
    ends = [e for ex, _ in transcripts.values() for _, e in ex]
    rows = [(seqid, source, "gene", min(starts), max(ends), strand, f"ID={gid};biotype=protein_coding")]
    for tid, (exons, cds) in transcripts.items():
        rows.append((seqid, source, "mRNA", exons[0][0], exons[-1][1], strand,
                     f"ID={tid};Parent={gid};biotype=protein_coding"))  # fmt: skip
        rows += [(seqid, source, "exon", s, e, strand, f"Parent={tid}") for s, e in exons]
        rows += [(seqid, source, "CDS", s, e, strand, f"Parent={tid}") for s, e in cds]
    return rows


def write_gff3(path: Path, rows) -> None:
    lines = ["##gff-version 3"]
    for seqid, source, kind, start, end, strand, attrs in rows:
        lines.append(f"{seqid}\t{source}\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")
    path.write_text("\n".join(lines) + "\n")


def main():
    rng = random.Random(20261027)
    chr1 = [rng.choice("ACGT") for _ in range(3000)]
    chr2 = [rng.choice("ACGT") for _ in range(1000)]
    r1 = [(101, 200), (301, 400), (501, 600)]  # 300 bp CDS, + strand
    r2 = [(1101, 1200), (1301, 1401)]  # 201 bp CDS, - strand
    r3 = [(2001, 2300)]  # 300 bp, + strand
    r4 = [(101, 400)]  # 300 bp, + strand, chr 2
    place(chr1, r1, orf(300, rng), "+")
    place(chr1, r2, orf(201, rng), "-")
    place(chr1, r3, orf(300, rng), "+")
    place(chr2, r4, orf(300, rng), "+")
    out = Path("inputs")
    out.mkdir(exist_ok=True)
    fasta = lambda name, s: f">{name}\n" + "\n".join("".join(s[i : i + 60]) for i in range(0, len(s), 60))  # noqa: E731
    (out / "genome.fa").write_text(fasta("1", chr1) + "\n" + fasta("2", chr2) + "\n")
    write_gff3(out / "reference.gff3", gene("1", "R1", "+", {"R1.1": (r1, r1)})
               + gene("1", "R2", "-", {"R2.1": (r2, r2)})
               + gene("1", "R3", "+", {"R3.1": (r3, r3)})
               + gene("2", "R4", "+", {"R4.1": (r4, r4)}))  # fmt: skip
    no_stop = [(101, 200), (301, 400), (501, 597)]  # CDS without its stop codon
    write_gff3(out / "focus.gff3", gene("1", "F1", "+", {"F1.1": (r1, no_stop)}, "focus")
               + gene("1", "F2", "-", {"F2.1": (r2, r2)}, "focus")
               + gene("2", "F4", "+", {"F4.1": (r4, r4)}, "focus"))  # fmt: skip
    other_r2 = [(1101, 1180), (1301, 1401)]
    write_gff3(out / "other.gff3", gene("1", "O1", "+", {"O1.1": (r1, r1)}, "other")
               + gene("1", "O2", "-", {"O2.1": (other_r2, other_r2)}, "other")
               + gene("1", "O3", "+", {"O3.1": (r3, r3)}, "other"))  # fmt: skip


if __name__ == "__main__":
    main()

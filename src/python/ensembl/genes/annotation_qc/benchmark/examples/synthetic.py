"""
SYNTHETIC demonstration experiment (not real data).

``write_synthetic_demo(target)`` writes two tiny genomes, two reference annotations
and five prediction runs plus ``synthetic_demo.json``. Every identifier is
prefixed SYNTH and the experiment is flagged ``synthetic: true`` so it can never be
mistaken for a real benchmark. It exercises:

    reference SYNTH_A   canonical tags, an Ensembl TEC ``unconfirmed_transcript``
                        (needs the audited retype step), a pseudogene, a lncRNA and
                        an annotated sequence plus an unannotated one (scope)
    reference SYNTH_B   a second genome without canonical tags (policy C unavailable)
    SYNTH_toolX_m1 / _m2  two runs (models) of the same tool on SYNTH_A
    SYNTH_toolY_hybrid  Tiberius-style GTF reusing gene IDs on two chromosomes
    SYNTH_toolZ_broken  an orphan exon: the comparison fails and is recorded
    SYNTH_toolX_m1_B    the same tool on the second genome

Cases (policy A): exact CDS, alternative-isoform exact (fails under B or C), stop-side
difference, different splice site, < 0.8 partial, missed, split, merge, strand
mismatch, Novel on a pseudogene, Novel with no overlap, prediction outside scope.

Usage: python -m ensembl.genes.annotation_qc.benchmark.examples.synthetic TARGET_DIR
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path


def _fasta(seqs: dict[str, int], seed: int) -> str:
    rng = random.Random(seed)
    out = []
    for name, length in seqs.items():
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        out.append(f">{name} SYNTHETIC")
        out += [seq[i : i + 60] for i in range(0, length, 60)]
    return "\n".join(out) + "\n"


def _gff_gene(
    seq, gid, strand, transcripts, biotype="protein_coding", gene_type="gene", name=None
):
    """transcripts: [(tid, exons, cds or None, tx_type, tags, tx_biotype)]"""
    starts = [s for t in transcripts for s, _ in t[1]]
    ends = [e for t in transcripts for _, e in t[1]]
    rows = [
        f"{seq}\tSYNTH\t{gene_type}\t{min(starts)}\t{max(ends)}\t.\t{strand}\t.\tID=gene:{gid};Name={name or gid};biotype={biotype};gene_id={gid}"
    ]
    for tid, exons, cds, ttype, tags, tbiotype in transcripts:
        tag = f";tag={tags}" if tags else ""
        rows.append(
            f"{seq}\tSYNTH\t{ttype}\t{min(s for s, _ in exons)}\t{max(e for _, e in exons)}\t.\t{strand}\t.\t"
            f"ID=transcript:{tid};Parent=gene:{gid};biotype={tbiotype}{tag};transcript_id={tid}"
        )
        for s, e in exons:
            rows.append(
                f"{seq}\tSYNTH\texon\t{s}\t{e}\t.\t{strand}\t.\tParent=transcript:{tid}"
            )
        phase = 0
        ordered = sorted(cds or [], reverse=(strand == "-"))
        for s, e in ordered:
            rows.append(
                f"{seq}\tSYNTH\tCDS\t{s}\t{e}\t.\t{strand}\t{phase}\tID=CDS:{tid};Parent=transcript:{tid}"
            )
            phase = (3 - ((e - s + 1 - phase) % 3)) % 3
    return rows


def _pred_gff(seq, gid, strand, cds, source="SYNTH_pred"):
    rows = [
        f"{seq}\t{source}\tgene\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\tID={gid};Name={gid}",
        f"{seq}\t{source}\tmRNA\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}",
    ]
    for i, (s, e) in enumerate(cds, 1):
        rows.append(
            f"{seq}\t{source}\texon\t{s}\t{e}\t.\t{strand}\t.\tID={gid}.t1.exon{i};Parent={gid}.t1"
        )
        rows.append(
            f"{seq}\t{source}\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID={gid}.t1.cds{i};Parent={gid}.t1"
        )
    return rows


def _pred_tiberius(seq, gid, strand, cds):
    rows = [
        f"{seq}\tSYNTH_hybrid\tgene\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\t{gid}",
        f"{seq}\tSYNTH_hybrid\ttranscript\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\t{gid}.t1",
    ]
    for s, e in cds:
        attrs = f'transcript_id "{gid}.t1"; gene_id "{gid}";'
        rows.append(f"{seq}\tSYNTH_hybrid\tCDS\t{s}\t{e}\t.\t{strand}\t0\t{attrs}")
        rows.append(f"{seq}\tSYNTH_hybrid\texon\t{s}\t{e}\t.\t{strand}\t0\t{attrs}")
    return rows


def write_synthetic_demo(target: str | Path) -> Path:
    t = Path(target)
    (t / "inputs").mkdir(parents=True, exist_ok=True)
    i = t / "inputs"
    (i / "SYNTH_genome_A.fa").write_text(
        _fasta({"chrA1": 20000, "chrA2": 20000, "chrA_unplaced": 5000}, 1)
    )
    (i / "SYNTH_genome_B.fa").write_text(_fasta({"chrB1": 15000}, 2))

    a = [
        "##gff-version 3",
        "##sequence-region   chrA1 1 20000",
        "##sequence-region   chrA2 1 20000",
        "#!genome-build SYNTHETIC_A",
    ]
    # g1: two coding isoforms; t1 canonical (shorter CDS), t2 longest CDS
    a += _gff_gene(
        "chrA1",
        "SYNTH_G1",
        "+",
        [
            (
                "SYNTH_T1a",
                [(1001, 1200), (1501, 1700), (2001, 2300)],
                [(1051, 1200), (1501, 1700), (2001, 2210)],
                "mRNA",
                "Ensembl_canonical",
                "protein_coding",
            ),
            (
                "SYNTH_T1b",
                [(1001, 1200), (1501, 1700), (1801, 1900), (2001, 2400)],
                [(1051, 1200), (1501, 1700), (1801, 1900), (2001, 2350)],
                "mRNA",
                "",
                "protein_coding",
            ),
            (
                "SYNTH_T1c",
                [(1001, 1200), (1501, 2300)],
                None,
                "lnc_RNA",
                "",
                "retained_intron",
            ),
        ],
        name="ALPHA1",
    )
    a += _gff_gene(
        "chrA1",
        "SYNTH_G2",
        "-",
        [
            (
                "SYNTH_T2",
                [(3001, 3600)],
                [(3051, 3551)],
                "mRNA",
                "Ensembl_canonical",
                "protein_coding",
            )
        ],
        name="BETA2",
    )
    a += _gff_gene(
        "chrA1",
        "SYNTH_G3",
        "+",
        [
            (
                "SYNTH_T3",
                [(5001, 5300), (6001, 6300)],
                [(5101, 5300), (6001, 6200)],
                "mRNA",
                "Ensembl_canonical",
                "protein_coding",
            )
        ],
        name="GAMMA3",
    )
    a += _gff_gene(
        "chrA1",
        "SYNTH_G4",
        "+",
        [
            (
                "SYNTH_T4",
                [(8001, 8900)],
                None,
                "pseudogenic_transcript",
                "",
                "processed_pseudogene",
            )
        ],
        biotype="processed_pseudogene",
        gene_type="pseudogene",
        name="DELTA4P",
    )
    a += _gff_gene(
        "chrA1",
        "SYNTH_G5",
        "+",
        [
            (
                "SYNTH_T5",
                [(10001, 10400), (11001, 11300), (12001, 12500)],
                [(10051, 10400), (11001, 11300), (12001, 12450)],
                "mRNA",
                "Ensembl_canonical",
                "protein_coding",
            )
        ],
        name="EPSILON5",
    )
    a += _gff_gene(
        "chrA1",
        "SYNTH_G6",
        "+",
        [
            (
                "SYNTH_T6",
                [(14001, 14500)],
                None,
                "unconfirmed_transcript",
                "Ensembl_canonical",
                "TEC",
            )
        ],
        biotype="TEC",
        gene_type="ncRNA_gene",
        name="ZETA6",
    )
    a += _gff_gene(
        "chrA2",
        "SYNTH_G7",
        "+",
        [
            (
                "SYNTH_T7",
                [(2001, 2300), (2601, 2900)],
                [(2051, 2300), (2601, 2850)],
                "mRNA",
                "Ensembl_canonical",
                "protein_coding",
            )
        ],
        name="ETA7",
    )
    a += _gff_gene(
        "chrA2",
        "SYNTH_G8",
        "-",
        [
            (
                "SYNTH_T8",
                [(5001, 5200), (5501, 5800)],
                [(5051, 5200), (5501, 5750)],
                "mRNA",
                "Ensembl_canonical",
                "protein_coding",
            )
        ],
        name="ALPHA1",
    )  # repeated display name
    (i / "SYNTH_reference_A.gff3").write_text("\n".join(a) + "\n")

    b = ["##gff-version 3", "#!genome-build SYNTHETIC_B (no canonical tags)"]
    b += _gff_gene(
        "chrB1",
        "SYNTH_H1",
        "+",
        [
            (
                "SYNTH_U1",
                [(1001, 1300), (2001, 2300)],
                [(1101, 1300), (2001, 2200)],
                "mRNA",
                "",
                "protein_coding",
            )
        ],
        name="THETA1",
    )
    b += _gff_gene(
        "chrB1",
        "SYNTH_H2",
        "-",
        [("SYNTH_U2", [(6001, 6600)], [(6101, 6500)], "mRNA", "", "protein_coding")],
        name="IOTA2",
    )
    (i / "SYNTH_reference_B.gff3").write_text("\n".join(b) + "\n")

    x1 = ["##gff-version 3"]
    x1 += _pred_gff(
        "chrA1", "SYNTH_x1_g1", "+", [(1051, 1200), (1501, 1700), (2001, 2210)]
    )  # exact canonical isoform
    x1 += _pred_gff("chrA1", "SYNTH_x1_g2", "-", [(3051, 3551)])  # exact single segment
    x1 += _pred_gff("chrA1", "SYNTH_x1_ps", "+", [(8101, 8700)])  # Novel on pseudogene
    x1 += _pred_gff(
        "chrA1", "SYNTH_x1_g5", "+", [(10051, 10400), (11001, 11300), (12001, 12400)]
    )  # stop side differs
    x1 += _pred_gff("chrA1", "SYNTH_x1_new", "+", [(16001, 16300)])  # Novel, no overlap
    x1 += _pred_gff("chrA2", "SYNTH_x1_g7", "+", [(2051, 2300), (2601, 2850)])  # exact
    x1 += _pred_gff(
        "chrA2", "SYNTH_x1_g8", "-", [(5051, 5200), (5451, 5750)]
    )  # different splice site
    x1 += _pred_gff(
        "chrA_unplaced", "SYNTH_x1_out", "+", [(1001, 1300)]
    )  # outside scope
    (i / "SYNTH_toolX_model1.gff3").write_text("\n".join(x1) + "\n")

    x2 = ["##gff-version 3"]
    x2 += _pred_gff(
        "chrA1",
        "SYNTH_x2_g1",
        "+",
        [(1051, 1200), (1501, 1700), (1801, 1900), (2001, 2350)],
    )  # exact longest isoform
    x2 += _pred_gff(
        "chrA1", "SYNTH_x2_anti", "+", [(3101, 3500)]
    )  # opposite strand of G2
    x2 += _pred_gff(
        "chrA1", "SYNTH_x2_g5a", "+", [(10051, 10400), (11001, 11300)]
    )  # split G5 (part 1)
    x2 += _pred_gff("chrA1", "SYNTH_x2_g5b", "+", [(12001, 12450)])  # split G5 (part 2)
    x2 += _pred_gff("chrA2", "SYNTH_x2_g7", "+", [(2051, 2300)])  # partial < 0.8
    x2 += _pred_gff("chrA2", "SYNTH_x2_g8", "-", [(5051, 5200), (5501, 5750)])  # exact
    (i / "SYNTH_toolX_model2.gff3").write_text("\n".join(x2) + "\n")

    y = []
    y += _pred_tiberius(
        "chrA1", "g1", "+", [(1051, 1200), (1501, 1700), (2001, 2210)]
    )  # exact G1 (ID g1 also used on chrA2)
    y += _pred_tiberius(
        "chrA1",
        "g2",
        "+",
        [(5101, 5300), (6001, 6200), (10051, 10400), (11001, 11300), (12001, 12450)],
    )  # merges G3 + G5
    y += _pred_tiberius(
        "chrA2", "g1", "+", [(2051, 2300), (2601, 2850)]
    )  # exact G7, same ID as chrA1 g1
    y += _pred_tiberius("chrA2", "g2", "-", [(5051, 5200), (5501, 5750)])  # exact G8
    (i / "SYNTH_toolY_hybrid.gtf").write_text("\n".join(y) + "\n")

    z = ["##gff-version 3"] + _pred_gff("chrA1", "SYNTH_z_g1", "+", [(1051, 1200)])
    z.append(
        "chrA1\tSYNTH\texon\t5101\t5300\t.\t+\t.\tParent=SYNTH_missing_transcript"
    )  # orphan exon -> comparator error
    (i / "SYNTH_toolZ_broken.gff3").write_text("\n".join(z) + "\n")

    xb = ["##gff-version 3"] + _pred_gff(
        "chrB1", "SYNTH_xb_h1", "+", [(1101, 1300), (2001, 2200)]
    )
    (i / "SYNTH_toolX_model1_genomeB.gff3").write_text("\n".join(xb) + "\n")

    config = {
        "config_version": 1,
        "experiment": {
            "id": "SYNTHETIC_demo",
            "title": "SYNTHETIC demonstration (not real data)",
            "synthetic": True,
            "description": "Tiny generated genomes and predictions that exercise every dashboard case.",
        },
        "output_dir": "workspace",
        "evaluation": {
            "evaluation_mode": "protein_coding",
            "reference_transcript_biotypes": ["protein_coding"],
            "genome_mismatch": "error",
        },
        "policies": {
            "A": {
                "reference_transcript_selection": "all",
                "query_transcript_selection": "all",
            },
            "B": {
                "reference_transcript_selection": "longest_cds",
                "query_transcript_selection": "longest_cds",
            },
            "C": {
                "reference_transcript_selection": "canonical",
                "query_transcript_selection": "longest_cds",
                "optional": True,
            },
        },
        "references": [
            {
                "id": "SYNTH_A",
                "species": "synthetic species A",
                "assembly": "SYNTH_A_v1",
                "annotation_source": "generated",
                "annotation_release": "synthetic",
                "annotation": "inputs/SYNTH_reference_A.gff3",
                "genome": "inputs/SYNTH_genome_A.fa",
                "annotation_preparation": {
                    "steps": [
                        {
                            "name": "retype_transcript_types",
                            "types": ["unconfirmed_transcript"],
                        }
                    ]
                },
                "scope": {"type": "reference_sequence_regions"},
            },
            {
                "id": "SYNTH_B",
                "species": "synthetic species B",
                "assembly": "SYNTH_B_v1",
                "annotation_source": "generated",
                "annotation_release": "synthetic (no canonical tags)",
                "annotation": "inputs/SYNTH_reference_B.gff3",
                "genome": "inputs/SYNTH_genome_B.fa",
                "scope": {"type": "all"},
            },
        ],
        "runs": [
            {
                "run_id": "SYNTH_toolX_m1",
                "reference": "SYNTH_A",
                "tool": "ToolX",
                "tool_version": "1.0",
                "model": "model1",
                "annotation": "inputs/SYNTH_toolX_model1.gff3",
                "provenance": "SYNTHETIC: generated by synthetic.py",
            },
            {
                "run_id": "SYNTH_toolX_m2",
                "reference": "SYNTH_A",
                "tool": "ToolX",
                "tool_version": "1.0",
                "model": "model2",
                "annotation": "inputs/SYNTH_toolX_model2.gff3",
                "provenance": "SYNTHETIC: second model of the same tool",
            },
            {
                "run_id": "SYNTH_toolY_hybrid",
                "reference": "SYNTH_A",
                "tool": "ToolY",
                "annotation": "inputs/SYNTH_toolY_hybrid.gtf",
                "provenance": "SYNTHETIC: Tiberius-style GTF reusing gene IDs on two chromosomes",
                "caveats": ["version and model unknown (left unknown on purpose)"],
            },
            {
                "run_id": "SYNTH_toolZ_broken",
                "reference": "SYNTH_A",
                "tool": "ToolZ",
                "annotation": "inputs/SYNTH_toolZ_broken.gff3",
                "provenance": "SYNTHETIC: contains an orphan exon so the comparison fails",
            },
            {
                "run_id": "SYNTH_toolX_m1_B",
                "reference": "SYNTH_B",
                "tool": "ToolX",
                "tool_version": "1.0",
                "model": "model1",
                "annotation": "inputs/SYNTH_toolX_model1_genomeB.gff3",
                "provenance": "SYNTHETIC: same tool on a second genome",
            },
        ],
    }
    path = t / "synthetic_demo.json"
    path.write_text(json.dumps(config, indent=2) + "\n")
    return path


if __name__ == "__main__":
    print(write_synthetic_demo(sys.argv[1] if len(sys.argv) > 1 else "synthetic_demo"))

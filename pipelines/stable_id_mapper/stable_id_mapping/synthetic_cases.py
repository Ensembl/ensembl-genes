"""Synthetic stable-ID mapping fixtures with known outcomes."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


@dataclass(frozen=True)
class GeneModel:
    gene_id: str
    transcript_id: str
    protein_id: str
    seqid: str
    gene_start: int
    gene_end: int
    strand: str
    exons: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class SyntheticCase:
    name: str
    description: str
    expected_decisions: dict[str, dict[str, int]]
    case_dir: Path

    @property
    def ref_fasta(self) -> Path:
        return self.case_dir / "ref.fa"

    @property
    def target_fasta(self) -> Path:
        return self.case_dir / "target.fa"

    @property
    def ref_gff(self) -> Path:
        return self.case_dir / "ref.gff3"

    @property
    def target_gff(self) -> Path:
        return self.case_dir / "target.gff3"

    @property
    def lifton_gff(self) -> Path:
        return self.case_dir / "lifton_projected.gff3"

    @property
    def metadata_path(self) -> Path:
        return self.case_dir / "case.json"


def available_case_names() -> tuple[str, ...]:
    return (
        "high_identity",
        "isoform_shift_gene_only",
        "unrelated_empty_lifton",
        "duplicated_competing_locus",
        "split_old_gene",
        "merged_old_genes",
        "strand_mismatch",
        "contig_mismatch",
    )


def write_synthetic_cases(
    output_dir: Path,
    case_names: Optional[Iterable[str]] = None,
) -> list[SyntheticCase]:
    selected = tuple(case_names) if case_names is not None else available_case_names()
    cases = []
    for name in selected:
        cases.append(_write_case_by_name(output_dir, name))
    return cases


def stable_id_mapping_command(
    case: SyntheticCase,
    run_dir: Optional[Path] = None,
    rules_config: str | Path = "pipelines/stable_id_mapper/stable_id_mapping_rules.json",
) -> list[str]:
    run_dir = run_dir or case.case_dir / "run"
    return [
        "python3",
        "pipelines/stable_id_mapper/run_stable_id_mapping.py",
        "--ref-fasta",
        str(case.ref_fasta),
        "--ref-gff",
        str(case.ref_gff),
        "--target-fasta",
        str(case.target_fasta),
        "--target-gff",
        str(case.target_gff),
        "--db-name",
        "synthetic_core_test",
        "--mapping-session-id",
        "1",
        "--gene-range",
        "ENSFAKEG:990001-990999",
        "--transcript-range",
        "ENSFAKET:990001-990999",
        "--translation-range",
        "ENSFAKEP:990001-990999",
        "--output-dir",
        str(run_dir),
        "--existing-lifton-gff",
        str(case.lifton_gff),
        "--rules-config",
        str(rules_config),
    ]


def format_shell_command(command: list[str]) -> str:
    return " ".join(_shell_quote(part) for part in command)


def _write_case_by_name(output_dir: Path, name: str) -> SyntheticCase:
    if name == "high_identity":
        return _write_high_identity(output_dir / name)
    if name == "isoform_shift_gene_only":
        return _write_isoform_shift_gene_only(output_dir / name)
    if name == "unrelated_empty_lifton":
        return _write_unrelated_empty_lifton(output_dir / name)
    if name == "duplicated_competing_locus":
        return _write_duplicated_competing_locus(output_dir / name)
    if name == "split_old_gene":
        return _write_split_old_gene(output_dir / name)
    if name == "merged_old_genes":
        return _write_merged_old_genes(output_dir / name)
    if name == "strand_mismatch":
        return _write_strand_mismatch(output_dir / name)
    if name == "contig_mismatch":
        return _write_contig_mismatch(output_dir / name)
    raise ValueError(f"Unknown synthetic case: {name}")


def _write_high_identity(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000001", "ENSFAKET000001", "ENSFAKEP000001", "chrRef", 100, 500, "+", ((100, 180), (300, 380), (450, 500))),
        GeneModel("ENSFAKEG000002", "ENSFAKET000002", "ENSFAKEP000002", "chrRef", 900, 1300, "+", ((900, 1000), (1100, 1300))),
        GeneModel("ENSFAKEG000003", "ENSFAKET000003", "ENSFAKEP000003", "chrRef", 1700, 2100, "-", ((1700, 1800), (1900, 2100))),
    )
    lifton = tuple(_with_seqid(model, "chrSynthetic") for model in old)
    target = (
        GeneModel("ENSFAKEG900001", "ENSFAKET900001", "ENSFAKEP900001", "chrSynthetic", 100, 500, "+", ((100, 180), (300, 380), (450, 500))),
        GeneModel("ENSFAKEG900002", "ENSFAKET900002", "ENSFAKEP900002", "chrSynthetic", 900, 1300, "+", ((900, 1000), (1100, 1300))),
        GeneModel("ENSFAKEG900003", "ENSFAKET900003", "ENSFAKEP900003", "chrSynthetic", 1700, 2100, "-", ((1700, 1800), (1900, 2100))),
    )
    return _write_case(
        case_dir,
        "high_identity",
        "Three projected old genes exactly match three target genes. Expected: all stable IDs map.",
        old,
        target,
        lifton,
        {
            "gene": {"mapped": 3, "missing": 0, "new": 0},
            "transcript": {"mapped": 3, "missing": 0, "new": 0},
            "translation": {"mapped": 3, "missing": 0, "new": 0},
        },
    )


def _write_isoform_shift_gene_only(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000101", "ENSFAKET000101", "ENSFAKEP000101", "chrRef", 100, 1000, "+", ((100, 200),)),
    )
    lifton = (_with_seqid(old[0], "chrSynthetic"),)
    target = (
        GeneModel("ENSFAKEG900101", "ENSFAKET900101", "ENSFAKEP900101", "chrSynthetic", 100, 1000, "+", ((900, 1000),)),
    )
    return _write_case(
        case_dir,
        "isoform_shift_gene_only",
        (
            "The projected old gene and target gene have the same locus, but their "
            "single transcript models occupy opposite ends of the gene. Expected: "
            "gene maps by coordinate overlap; transcript and translation do not map."
        ),
        old,
        target,
        lifton,
        {
            "gene": {"mapped": 1, "missing": 0, "new": 0},
            "transcript": {"mapped": 0, "missing": 1, "new": 1},
            "translation": {"mapped": 0, "missing": 1, "new": 1},
        },
    )


def _write_unrelated_empty_lifton(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000201", "ENSFAKET000201", "ENSFAKEP000201", "chrRef", 100, 300, "+", ((100, 300),)),
        GeneModel("ENSFAKEG000202", "ENSFAKET000202", "ENSFAKEP000202", "chrRef", 600, 900, "+", ((600, 700), (800, 900))),
        GeneModel("ENSFAKEG000203", "ENSFAKET000203", "ENSFAKEP000203", "chrRef", 1200, 1600, "-", ((1200, 1300), (1500, 1600))),
    )
    target = (
        GeneModel("ENSFAKEG900201", "ENSFAKET900201", "ENSFAKEP900201", "chrTarget", 5000, 5300, "+", ((5000, 5300),)),
        GeneModel("ENSFAKEG900202", "ENSFAKET900202", "ENSFAKEP900202", "chrTarget", 5600, 5900, "+", ((5600, 5700), (5800, 5900))),
        GeneModel("ENSFAKEG900203", "ENSFAKET900203", "ENSFAKEP900203", "chrTarget", 6200, 6600, "-", ((6200, 6300), (6500, 6600))),
    )
    return _write_case(
        case_dir,
        "unrelated_empty_lifton",
        "No LiftOn-projected genes are present. Expected: all old features missing; all target features new.",
        old,
        target,
        (),
        {
            "gene": {"mapped": 0, "missing": 3, "new": 3},
            "transcript": {"mapped": 0, "missing": 3, "new": 3},
            "translation": {"mapped": 0, "missing": 3, "new": 3},
        },
    )


def _write_duplicated_competing_locus(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000301", "ENSFAKET000301", "ENSFAKEP000301", "chrRef", 100, 500, "+", ((100, 150), (320, 500))),
        GeneModel("ENSFAKEG000302", "ENSFAKET000302", "ENSFAKEP000302", "chrRef", 100, 500, "+", ((100, 180), (300, 500))),
    )
    lifton = tuple(_with_seqid(model, "chrSynthetic") for model in old)
    target = (
        GeneModel("ENSFAKEG900301", "ENSFAKET900301", "ENSFAKEP900301", "chrSynthetic", 100, 500, "+", ((100, 180), (300, 500))),
    )
    return _write_case(
        case_dir,
        "duplicated_competing_locus",
        (
            "Two projected old genes compete for one target locus. The better "
            "structural match should claim the target; the other old gene is missing."
        ),
        old,
        target,
        lifton,
        {
            "gene": {"mapped": 1, "missing": 1, "new": 0},
            "transcript": {"mapped": 1, "missing": 1, "new": 0},
            "translation": {"mapped": 1, "missing": 1, "new": 0},
        },
    )


def _write_split_old_gene(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000401", "ENSFAKET000401", "ENSFAKEP000401", "chrRef", 100, 1000, "+", ((100, 300), (700, 1000))),
    )
    lifton = (_with_seqid(old[0], "chrSynthetic"),)
    target = (
        GeneModel("ENSFAKEG900401", "ENSFAKET900401", "ENSFAKEP900401", "chrSynthetic", 100, 300, "+", ((100, 300),)),
        GeneModel("ENSFAKEG900402", "ENSFAKET900402", "ENSFAKEP900402", "chrSynthetic", 700, 1000, "+", ((700, 1000),)),
    )
    return _write_case(
        case_dir,
        "split_old_gene",
        (
            "One old projected gene overlaps two target genes. Current conservative "
            "behavior maps the old stable ID to one representative target and leaves "
            "the other target feature new."
        ),
        old,
        target,
        lifton,
        {
            "gene": {"mapped": 1, "missing": 0, "new": 1},
            "transcript": {"mapped": 1, "missing": 0, "new": 1},
            "translation": {"mapped": 1, "missing": 0, "new": 1},
        },
    )


def _write_merged_old_genes(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000501", "ENSFAKET000501", "ENSFAKEP000501", "chrRef", 100, 300, "+", ((100, 300),)),
        GeneModel("ENSFAKEG000502", "ENSFAKET000502", "ENSFAKEP000502", "chrRef", 700, 1000, "+", ((700, 1000),)),
    )
    lifton = tuple(_with_seqid(model, "chrSynthetic") for model in old)
    target = (
        GeneModel("ENSFAKEG900501", "ENSFAKET900501", "ENSFAKEP900501", "chrSynthetic", 100, 1000, "+", ((100, 300), (700, 1000))),
    )
    return _write_case(
        case_dir,
        "merged_old_genes",
        (
            "Two old projected genes overlap one merged target gene. Current "
            "one-to-one behavior maps one old stable ID and marks the other old gene missing."
        ),
        old,
        target,
        lifton,
        {
            "gene": {"mapped": 1, "missing": 1, "new": 0},
            "transcript": {"mapped": 1, "missing": 1, "new": 0},
            "translation": {"mapped": 1, "missing": 1, "new": 0},
        },
    )


def _write_strand_mismatch(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000601", "ENSFAKET000601", "ENSFAKEP000601", "chrRef", 100, 500, "+", ((100, 500),)),
    )
    lifton = (_with_seqid(old[0], "chrSynthetic"),)
    target = (
        GeneModel("ENSFAKEG900601", "ENSFAKET900601", "ENSFAKEP900601", "chrSynthetic", 100, 500, "-", ((100, 500),)),
    )
    return _write_case(
        case_dir,
        "strand_mismatch",
        "Projected and target genes overlap perfectly in span but are on opposite strands. Expected: no mapping.",
        old,
        target,
        lifton,
        {
            "gene": {"mapped": 0, "missing": 1, "new": 1},
            "transcript": {"mapped": 0, "missing": 1, "new": 1},
            "translation": {"mapped": 0, "missing": 1, "new": 1},
        },
    )


def _write_contig_mismatch(case_dir: Path) -> SyntheticCase:
    old = (
        GeneModel("ENSFAKEG000701", "ENSFAKET000701", "ENSFAKEP000701", "chrRef", 100, 500, "+", ((100, 500),)),
    )
    lifton = (_with_seqid(old[0], "chrSyntheticA"),)
    target = (
        GeneModel("ENSFAKEG900701", "ENSFAKET900701", "ENSFAKEP900701", "chrSyntheticB", 100, 500, "+", ((100, 500),)),
    )
    return _write_case(
        case_dir,
        "contig_mismatch",
        "Projected and target genes overlap numerically but are on different contigs. Expected: no mapping.",
        old,
        target,
        lifton,
        {
            "gene": {"mapped": 0, "missing": 1, "new": 1},
            "transcript": {"mapped": 0, "missing": 1, "new": 1},
            "translation": {"mapped": 0, "missing": 1, "new": 1},
        },
    )


def _write_case(
    case_dir: Path,
    name: str,
    description: str,
    ref_models: tuple[GeneModel, ...],
    target_models: tuple[GeneModel, ...],
    lifton_models: tuple[GeneModel, ...],
    expected_decisions: dict[str, dict[str, int]],
) -> SyntheticCase:
    case_dir.mkdir(parents=True, exist_ok=True)
    case = SyntheticCase(
        name=name,
        description=description,
        expected_decisions=expected_decisions,
        case_dir=case_dir,
    )
    _write_fasta(case.ref_fasta, _seqids(ref_models), _max_end(ref_models) + 1000)
    _write_fasta(
        case.target_fasta,
        _seqids(target_models + lifton_models),
        max(_max_end(target_models), _max_end(lifton_models)) + 1000,
    )
    _write_gff(case.ref_gff, ref_models)
    _write_gff(case.target_gff, target_models)
    _write_gff(case.lifton_gff, lifton_models)
    _write_metadata(case)
    return case


def _write_metadata(case: SyntheticCase) -> None:
    case.metadata_path.write_text(
        json.dumps(
            {
                "name": case.name,
                "description": case.description,
                "expected_decisions": case.expected_decisions,
                "command": stable_id_mapping_command(case),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def _write_fasta(path: Path, seqids: tuple[str, ...], length: int) -> None:
    sequence = ("ACGT" * ((length // 4) + 1))[:length]
    with path.open("w", encoding="utf-8") as handle:
        for seqid in seqids:
            handle.write(f">{seqid}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def _write_gff(path: Path, models: tuple[GeneModel, ...]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for index, model in enumerate(models, start=1):
            tx_start = min(start for start, _end in model.exons)
            tx_end = max(end for _start, end in model.exons)
            gene_id = f"gene:{model.gene_id}.1"
            transcript_id = f"transcript:{model.transcript_id}.1"
            protein_id = f"{model.protein_id}.1"
            handle.write(
                "\t".join(
                    [
                        model.seqid,
                        "synthetic",
                        "gene",
                        str(model.gene_start),
                        str(model.gene_end),
                        ".",
                        model.strand,
                        ".",
                        f"ID={gene_id}",
                    ]
                )
                + "\n"
            )
            handle.write(
                "\t".join(
                    [
                        model.seqid,
                        "synthetic",
                        "mRNA",
                        str(tx_start),
                        str(tx_end),
                        ".",
                        model.strand,
                        ".",
                        f"ID={transcript_id};Parent={gene_id}",
                    ]
                )
                + "\n"
            )
            for exon_index, (start, end) in enumerate(model.exons, start=1):
                exon_id = f"exon:{model.transcript_id}.{exon_index}"
                cds_id = f"CDS:{model.protein_id}.{exon_index}"
                handle.write(
                    "\t".join(
                        [
                            model.seqid,
                            "synthetic",
                            "exon",
                            str(start),
                            str(end),
                            ".",
                            model.strand,
                            ".",
                            f"ID={exon_id};Parent={transcript_id}",
                        ]
                    )
                    + "\n"
                )
                handle.write(
                    "\t".join(
                        [
                            model.seqid,
                            "synthetic",
                            "CDS",
                            str(start),
                            str(end),
                            ".",
                            model.strand,
                            "0",
                            f"ID={cds_id};Parent={transcript_id};protein_id={protein_id}",
                        ]
                    )
                    + "\n"
                )


def _with_seqid(model: GeneModel, seqid: str) -> GeneModel:
    return GeneModel(
        gene_id=model.gene_id,
        transcript_id=model.transcript_id,
        protein_id=model.protein_id,
        seqid=seqid,
        gene_start=model.gene_start,
        gene_end=model.gene_end,
        strand=model.strand,
        exons=model.exons,
    )


def _max_end(models: tuple[GeneModel, ...]) -> int:
    if not models:
        return 1000
    return max(model.gene_end for model in models)


def _seqids(models: tuple[GeneModel, ...]) -> tuple[str, ...]:
    seqids = tuple(sorted({model.seqid for model in models}))
    return seqids or ("chrSynthetic",)


def _shell_quote(value: str) -> str:
    if not value:
        return "''"
    safe = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_./:-")
    if all(char in safe for char in value):
        return value
    return "'" + value.replace("'", "'\"'\"'") + "'"

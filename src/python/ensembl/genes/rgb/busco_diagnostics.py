from __future__ import annotations

import argparse
import json
import os
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from ensembl.genes.rgb.io import ensure_dir


BUSCO_STATUS_RANK = {
    "Missing": 0,
    "Fragmented": 1,
    "Complete": 2,
    "Duplicated": 2,
}


BUSCO_DIAGNOSTIC_COLUMNS = [
    "busco_id",
    "orthodb_group",
    "taxid",
    "expected_length",
    "seq_region_name",
    "seq_region_start",
    "seq_region_end",
    "seq_region_strand",
    "genome_status",
    "protein_status",
    "problem_type",
    "diagnostic_class",
    "review_priority",
    "suggested_action",
    "core_gene_count",
    "core_protein_coding_count",
    "core_gene_ids",
    "core_biotypes",
    "layer_gene_count",
    "layer_biotypes",
    "layer_logic_names",
    "hmmer_status",
    "hmmer_evalue",
    "hmmer_hmm_coverage",
    "hmmer_protein_coverage",
    "hmmer_issues",
    "good_isoform_count",
    "good_isoforms",
    "summary",
]


ISOFORM_COLUMNS = [
    "busco_id",
    "seq_region_name",
    "seq_region_start",
    "seq_region_end",
    "gene_id",
    "transcript_id",
    "e_value",
    "hmm_coverage",
    "protein_coverage",
    "protein_length",
    "is_best_hit",
]


def parse_busco_id(busco_id: object) -> dict[str, str]:
    value = "" if pd.isna(busco_id) else str(busco_id)
    match = re.match(r"(\d+)at(\d+)", value)
    if match:
        return {"orthodb_group": match.group(1), "taxid": match.group(2)}
    return {"orthodb_group": value, "taxid": "unknown"}


def _safe_int(value: object, default: int = 0) -> int:
    try:
        if pd.isna(value):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _safe_float(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _string_value(value: object, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    return str(value)


def parse_busco_table(path: str, mode: str) -> pd.DataFrame:
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            busco_id = parts[0]
            status = parts[1]
            parsed = parse_busco_id(busco_id)
            row = {
                "busco_id": busco_id,
                "status": status,
                "orthodb_group": parsed["orthodb_group"],
                "taxid": parsed["taxid"],
                "seq_region_name": "",
                "seq_region_start": 0,
                "seq_region_end": 0,
                "seq_region_strand": "",
                "score": 0.0,
                "length": 0,
                "hit_id": "",
            }
            if mode == "genome" and status != "Missing" and len(parts) >= 8:
                row.update(
                    {
                        "seq_region_name": parts[2],
                        "seq_region_start": _safe_int(parts[3]),
                        "seq_region_end": _safe_int(parts[4]),
                        "seq_region_strand": parts[5],
                        "score": _safe_float(parts[6]),
                        "length": _safe_int(parts[7]),
                    }
                )
            elif mode == "protein" and status != "Missing" and len(parts) >= 3:
                row["hit_id"] = parts[2]
                if len(parts) >= 5:
                    row["score"] = _safe_float(parts[3])
                    row["length"] = _safe_int(parts[4])
            rows.append(row)
    return pd.DataFrame(rows)


def load_busco_metadata(
    busco_dataset_dir: Optional[str], busco_ids: Iterable[str]
) -> pd.DataFrame:
    rows = []
    lengths: dict[str, dict[str, float]] = {}
    if busco_dataset_dir:
        root = Path(busco_dataset_dir)
        candidates = [
            root / "lengths_cutoff",
            *root.glob("*/lengths_cutoff"),
            *root.glob("lineages/*/lengths_cutoff"),
        ]
        lengths_file = next((path for path in candidates if path.exists()), None)
        if lengths_file:
            with open(lengths_file, "r", encoding="utf-8") as handle:
                for raw in handle:
                    if raw.startswith("#") or not raw.strip():
                        continue
                    parts = raw.split()
                    if len(parts) >= 3:
                        lengths[parts[0]] = {
                            "length_sd": _safe_float(parts[1]),
                            "expected_length": _safe_float(parts[2]),
                        }
    for busco_id in busco_ids:
        parsed = parse_busco_id(busco_id)
        rows.append(
            {
                "busco_id": busco_id,
                "orthodb_group": parsed["orthodb_group"],
                "taxid": parsed["taxid"],
                "expected_length": lengths.get(busco_id, {}).get(
                    "expected_length", pd.NA
                ),
                "length_sd": lengths.get(busco_id, {}).get("length_sd", pd.NA),
            }
        )
    return pd.DataFrame(rows)


def identify_busco_problems(
    genome_busco: pd.DataFrame, protein_busco: pd.DataFrame
) -> pd.DataFrame:
    if genome_busco.empty or protein_busco.empty:
        return pd.DataFrame()
    genome = genome_busco.add_prefix("genome_")
    protein = protein_busco.add_prefix("protein_")
    merged = genome.merge(
        protein, left_on="genome_busco_id", right_on="protein_busco_id", how="inner"
    )
    rows = []
    for _, row in merged.iterrows():
        genome_status = _string_value(row.get("genome_status"))
        protein_status = _string_value(row.get("protein_status"))
        if BUSCO_STATUS_RANK.get(protein_status, 0) >= BUSCO_STATUS_RANK.get(
            genome_status, 0
        ):
            continue
        if genome_status == "Complete" and protein_status == "Missing":
            problem_type = "genome_complete_protein_missing"
        elif genome_status == "Complete" and protein_status == "Fragmented":
            problem_type = "genome_complete_protein_fragmented"
        elif genome_status == "Fragmented" and protein_status == "Missing":
            problem_type = "genome_fragmented_protein_missing"
        else:
            problem_type = "protein_mode_worse_than_genome_mode"
        rows.append(
            {
                "busco_id": row["genome_busco_id"],
                "orthodb_group": row["genome_orthodb_group"],
                "taxid": row["genome_taxid"],
                "seq_region_name": row["genome_seq_region_name"],
                "seq_region_start": row["genome_seq_region_start"],
                "seq_region_end": row["genome_seq_region_end"],
                "seq_region_strand": row["genome_seq_region_strand"],
                "genome_status": genome_status,
                "protein_status": protein_status,
                "problem_type": problem_type,
            }
        )
    return pd.DataFrame(rows)


def parse_hmmer_output(path: str, core_identifiers: Iterable[str]) -> dict:
    if not path or not os.path.exists(path):
        return {"status": "file_not_found"}
    allowed_ids = {
        str(identifier).split(".")[0]
        for identifier in core_identifiers
        if _string_value(identifier)
    }
    hits = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            parts = raw.split()
            if len(parts) < 19:
                continue
            transcript_id = parts[0]
            transcript_prefix = transcript_id.split(".")[0]
            gene_match = re.match(r"([A-Z]+G\d+)", transcript_prefix)
            gene_id = gene_match.group(1) if gene_match else transcript_prefix
            if (
                allowed_ids
                and gene_id not in allowed_ids
                and transcript_prefix not in allowed_ids
            ):
                continue
            hit = {
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "tlen": _safe_int(parts[2]),
                "qlen": _safe_int(parts[5]),
                "e_value": _safe_float(parts[6], float("inf")),
                "score": _safe_float(parts[7]),
                "bias": _safe_float(parts[8]),
                "domain_num": _safe_int(parts[9]),
                "domain_total": _safe_int(parts[10]),
                "hmm_from": _safe_int(parts[15]),
                "hmm_to": _safe_int(parts[16]),
                "ali_from": _safe_int(parts[17]),
                "ali_to": _safe_int(parts[18]),
            }
            hits.append(hit)
    if not hits:
        return {"status": "no_hits_for_core_genes"}
    best = min(hits, key=lambda item: item["e_value"])
    hmm_coverage = (
        (best["hmm_to"] - best["hmm_from"] + 1) / best["qlen"] * 100
        if best["qlen"]
        else 0.0
    )
    protein_coverage = (
        (best["ali_to"] - best["ali_from"] + 1) / best["tlen"] * 100
        if best["tlen"]
        else 0.0
    )
    issues = []
    if best["e_value"] > 1e-10:
        issues.append("poor_evalue")
    if hmm_coverage < 80:
        issues.append("low_hmm_coverage")
    if best["qlen"] and best["tlen"] < best["qlen"] * 0.5:
        issues.append("protein_too_short")
    if best["domain_total"] > 1:
        issues.append("fragmented_match")
    if protein_coverage < 50:
        issues.append("low_protein_coverage")

    good_isoforms = []
    for hit in hits:
        hit_hmm_cov = (
            (hit["hmm_to"] - hit["hmm_from"] + 1) / hit["qlen"] * 100
            if hit["qlen"]
            else 0.0
        )
        hit_prot_cov = (
            (hit["ali_to"] - hit["ali_from"] + 1) / hit["tlen"] * 100
            if hit["tlen"]
            else 0.0
        )
        if (
            hit["e_value"] <= 1e-10
            and hit_hmm_cov >= 80
            and hit["domain_total"] == 1
            and hit_prot_cov >= 50
        ):
            good_isoforms.append(
                {
                    "transcript_id": hit["transcript_id"],
                    "gene_id": hit["gene_id"],
                    "e_value": hit["e_value"],
                    "hmm_coverage": hit_hmm_cov,
                    "protein_coverage": hit_prot_cov,
                    "protein_length": hit["tlen"],
                    "is_best_hit": hit["transcript_id"] == best["transcript_id"],
                }
            )
    return {
        "status": "hit_found",
        "best_hit": best,
        "hmm_coverage": hmm_coverage,
        "protein_coverage": protein_coverage,
        "issues": issues,
        "good_isoforms": good_isoforms,
    }


def _overlap_rows(
    df: pd.DataFrame, seq_region_name: object, start: object, end: object
) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    required = {"seq_region_name", "seq_region_start", "seq_region_end"}
    if not required.issubset(df.columns):
        return pd.DataFrame(columns=df.columns)
    seq = _string_value(seq_region_name)
    qstart = _safe_int(start)
    qend = _safe_int(end)
    return df[
        (df["seq_region_name"].astype(str) == seq)
        & (df["seq_region_start"].map(_safe_int) <= qend)
        & (df["seq_region_end"].map(_safe_int) >= qstart)
    ].copy()


def _counts_text(counter: Counter) -> str:
    if not counter:
        return ""
    return ",".join(
        f"{key}({counter[key]})"
        for key in sorted(counter, key=lambda item: (-counter[item], item))
    )


def _has_translation(translations: pd.DataFrame, transcript_ids: set[int]) -> bool:
    if (
        translations.empty
        or not transcript_ids
        or "transcript_id" not in translations.columns
    ):
        return False
    return bool(
        set(translations["transcript_id"].dropna().map(_safe_int)) & transcript_ids
    )


def classify_busco_problem(
    problem: pd.Series,
    core_genes: pd.DataFrame,
    layer_genes: pd.DataFrame,
    core_transcripts: pd.DataFrame,
    core_translations: pd.DataFrame,
    hmmer_dir: Optional[str] = None,
) -> tuple[dict, list[dict]]:
    core_hits = _overlap_rows(
        core_genes,
        problem["seq_region_name"],
        problem["seq_region_start"],
        problem["seq_region_end"],
    )
    layer_hits = _overlap_rows(
        layer_genes,
        problem["seq_region_name"],
        problem["seq_region_start"],
        problem["seq_region_end"],
    )
    core_ids = (
        set(core_hits["gene_id"].dropna().map(_safe_int))
        if not core_hits.empty and "gene_id" in core_hits.columns
        else set()
    )
    tx_hits = (
        core_transcripts[
            core_transcripts["gene_id"].map(_safe_int).isin(core_ids)
        ].copy()
        if not core_transcripts.empty
        and core_ids
        and "gene_id" in core_transcripts.columns
        else pd.DataFrame()
    )
    transcript_ids = (
        set(tx_hits["transcript_id"].dropna().map(_safe_int))
        if not tx_hits.empty and "transcript_id" in tx_hits.columns
        else set()
    )
    core_biotypes = Counter(
        core_hits.get("biotype", pd.Series(dtype=str)).fillna("unknown").astype(str)
    )
    layer_biotypes = Counter(
        layer_hits.get("biotype", pd.Series(dtype=str)).fillna("unknown").astype(str)
    )
    layer_logic = Counter(
        layer_hits.get("logic_name", pd.Series(dtype=str)).fillna("unknown").astype(str)
    )
    protein_coding = int(core_biotypes.get("protein_coding", 0))
    translated = _has_translation(core_translations, transcript_ids)

    core_stable_ids = (
        [
            str(value)
            for value in core_hits.get("stable_id", pd.Series(dtype=str))
            .dropna()
            .tolist()
        ]
        if not core_hits.empty
        else []
    )
    tx_stable_ids = (
        [
            str(value)
            for value in tx_hits.get("stable_id", pd.Series(dtype=str))
            .dropna()
            .tolist()
        ]
        if not tx_hits.empty
        else []
    )
    hmmer = {"status": ""}
    if hmmer_dir and (core_stable_ids or tx_stable_ids):
        hmmer = parse_hmmer_output(
            os.path.join(hmmer_dir, f"{problem['busco_id']}.out"),
            core_stable_ids + tx_stable_ids,
        )

    if core_hits.empty:
        if layer_hits.empty:
            diagnostic_class = "no_core_or_layer_model"
            priority = "P1"
            action = "inspect_genome_busco_locus_and_assembly_context"
            summary = "Genome BUSCO finds this locus, but no overlapping core or layer gene model was extracted."
        else:
            diagnostic_class = "layer_evidence_not_selected"
            priority = "P1"
            action = "rescue_or_reprioritise_layer_candidate"
            summary = "Genome BUSCO finds this locus and layer models overlap it, but no core gene was built."
    elif protein_coding == 0:
        diagnostic_class = "non_coding_core_annotation"
        priority = "P1"
        action = "review_biotype_assignment_or_pseudogene_filtering"
        summary = "Core overlaps the BUSCO locus, but no overlapping core gene is protein_coding."
    elif not transcript_ids:
        diagnostic_class = "protein_coding_without_transcripts"
        priority = "P1"
        action = "inspect_core_transcript_export_or_gene_model_integrity"
        summary = "Overlapping core protein-coding genes have no extracted transcripts."
    elif not translated:
        diagnostic_class = "protein_coding_without_translation"
        priority = "P1"
        action = "inspect_cds_translation_and_frameshift_filters"
        summary = "Overlapping protein-coding transcripts lack extracted translations."
    elif hmmer.get("status") == "hit_found" and hmmer.get("good_isoforms"):
        diagnostic_class = "alternative_isoform_satisfies_busco"
        priority = "P1"
        action = "review_canonical_transcript_selection"
        summary = "HMMER found one or more non-canonical/core isoforms that satisfy BUSCO-like criteria."
    else:
        diagnostic_class = "protein_coding_not_busco_satisfying"
        priority = "P2"
        action = "inspect_cds_completeness_hmmer_and_model_structure"
        summary = "Core has protein-coding representation, but protein-mode BUSCO is worse than genome-mode BUSCO."

    good_isoforms = (
        hmmer.get("good_isoforms", []) if hmmer.get("status") == "hit_found" else []
    )
    isoform_rows = [
        {
            "busco_id": problem["busco_id"],
            "seq_region_name": problem["seq_region_name"],
            "seq_region_start": problem["seq_region_start"],
            "seq_region_end": problem["seq_region_end"],
            **isoform,
        }
        for isoform in good_isoforms
    ]
    hmmer_best = hmmer.get("best_hit", {}) if hmmer.get("status") == "hit_found" else {}
    row = {
        "busco_id": problem["busco_id"],
        "orthodb_group": problem.get("orthodb_group", ""),
        "taxid": problem.get("taxid", ""),
        "expected_length": problem.get("expected_length", pd.NA),
        "seq_region_name": problem["seq_region_name"],
        "seq_region_start": problem["seq_region_start"],
        "seq_region_end": problem["seq_region_end"],
        "seq_region_strand": problem["seq_region_strand"],
        "genome_status": problem["genome_status"],
        "protein_status": problem["protein_status"],
        "problem_type": problem["problem_type"],
        "diagnostic_class": diagnostic_class,
        "review_priority": priority,
        "suggested_action": action,
        "core_gene_count": len(core_hits),
        "core_protein_coding_count": protein_coding,
        "core_gene_ids": ",".join(core_stable_ids),
        "core_biotypes": _counts_text(core_biotypes),
        "layer_gene_count": len(layer_hits),
        "layer_biotypes": _counts_text(layer_biotypes),
        "layer_logic_names": _counts_text(layer_logic),
        "hmmer_status": hmmer.get("status", ""),
        "hmmer_evalue": hmmer_best.get("e_value", pd.NA),
        "hmmer_hmm_coverage": hmmer.get("hmm_coverage", pd.NA),
        "hmmer_protein_coverage": hmmer.get("protein_coverage", pd.NA),
        "hmmer_issues": ",".join(hmmer.get("issues", [])),
        "good_isoform_count": len(good_isoforms),
        "good_isoforms": ";".join(str(item["transcript_id"]) for item in good_isoforms),
        "summary": summary,
    }
    return row, isoform_rows


def build_busco_diagnostic_audit(
    genome_busco: pd.DataFrame,
    protein_busco: pd.DataFrame,
    core_genes: pd.DataFrame,
    layer_genes: pd.DataFrame,
    core_transcripts: pd.DataFrame,
    core_translations: pd.DataFrame,
    busco_metadata: Optional[pd.DataFrame] = None,
    hmmer_dir: Optional[str] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    problems = identify_busco_problems(genome_busco, protein_busco)
    if busco_metadata is not None and not busco_metadata.empty and not problems.empty:
        problems = problems.merge(
            busco_metadata, on=["busco_id", "orthodb_group", "taxid"], how="left"
        )
    rows = []
    isoforms = []
    for _, problem in problems.iterrows():
        row, isoform_rows = classify_busco_problem(
            problem,
            core_genes,
            layer_genes,
            core_transcripts,
            core_translations,
            hmmer_dir=hmmer_dir,
        )
        rows.append(row)
        isoforms.extend(isoform_rows)
    audit = pd.DataFrame(rows, columns=BUSCO_DIAGNOSTIC_COLUMNS)
    isoform_df = pd.DataFrame(isoforms, columns=ISOFORM_COLUMNS)
    summary = (
        audit.groupby(["diagnostic_class", "review_priority"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["review_priority", "count"], ascending=[True, False])
        if not audit.empty
        else pd.DataFrame(columns=["diagnostic_class", "review_priority", "count"])
    )
    return audit, isoform_df, summary


def write_busco_diagnostic_bed(audit: pd.DataFrame, path: str, pad_bp: int = 0) -> None:
    ensure_dir(os.path.dirname(path))
    columns = ["chrom", "start0", "end0", "name", "score", "strand"]
    if audit.empty:
        pd.DataFrame(columns=columns).to_csv(path, sep="\t", header=False, index=False)
        return
    rows = []
    for _, row in audit.iterrows():
        start = max(0, _safe_int(row["seq_region_start"]) - pad_bp - 1)
        end = _safe_int(row["seq_region_end"]) + pad_bp
        strand = {"+": "+", "-": "-", "1": "+", "-1": "-"}.get(
            _string_value(row["seq_region_strand"]), "."
        )
        rows.append(
            {
                "chrom": row["seq_region_name"],
                "start0": start,
                "end0": end,
                "name": f"BUSCO:{row['busco_id']}:{row['diagnostic_class']}",
                "score": 1000 if row["review_priority"] == "P1" else 700,
                "strand": strand,
            }
        )
    pd.DataFrame(rows, columns=columns).to_csv(
        path, sep="\t", header=False, index=False
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="BUSCO genome-vs-protein diagnostic audit over RGB extract tables"
    )
    parser.add_argument("--genome_busco", required=True)
    parser.add_argument("--protein_busco", required=True)
    parser.add_argument("--extract_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--format", choices=["tsv", "csv", "parquet"], default="tsv")
    parser.add_argument("--hmmer_dir")
    parser.add_argument("--busco_dataset_dir")
    parser.add_argument("--bed_pad_bp", type=int, default=0)
    args = parser.parse_args(argv)

    def read_named(base: str) -> pd.DataFrame:
        path = os.path.join(args.extract_dir, f"{base}.{args.format}")
        if args.format == "parquet":
            return pd.read_parquet(path)
        return pd.read_csv(path, sep="\t" if args.format == "tsv" else ",")

    genome = parse_busco_table(args.genome_busco, "genome")
    protein = parse_busco_table(args.protein_busco, "protein")
    metadata = load_busco_metadata(
        args.busco_dataset_dir, set(genome["busco_id"]) | set(protein["busco_id"])
    )
    audit, isoforms, summary = build_busco_diagnostic_audit(
        genome,
        protein,
        read_named("core_genes"),
        read_named("layer_genes"),
        read_named("core_transcripts"),
        read_named("core_translations"),
        metadata,
        args.hmmer_dir,
    )
    ensure_dir(args.output_dir)
    sep = "\t" if args.format == "tsv" else ","
    audit.to_csv(
        os.path.join(args.output_dir, f"busco_diagnostic_audit.{args.format}"),
        sep=sep,
        index=False,
    )
    isoforms.to_csv(
        os.path.join(args.output_dir, f"busco_isoform_suggestions.{args.format}"),
        sep=sep,
        index=False,
    )
    summary.to_csv(
        os.path.join(args.output_dir, f"busco_diagnostic_summary.{args.format}"),
        sep=sep,
        index=False,
    )
    write_busco_diagnostic_bed(
        audit,
        os.path.join(args.output_dir, "busco_diagnostic_issues.bed"),
        args.bed_pad_bp,
    )
    with open(
        os.path.join(args.output_dir, "busco_diagnostic_summary.json"),
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {"problem_count": len(audit), "isoform_suggestion_count": len(isoforms)},
            handle,
            indent=2,
        )
    print(f"[busco-diagnostics] wrote {len(audit)} problem rows -> {args.output_dir}")
    return 0

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
"""
Writers for the files added by output contract 1.0 of pairwise-compare.

    metrics.tsv                    canonical metric table (contract.py)
    evaluated_transcripts.tsv.gz   every parsed transcript, evaluated or excluded
    per_sequence.tsv               per-sequence counts (drill-down only)
    splice_sites_by_position.tsv   only when the reference selection keeps one
    cds_introns_by_length.tsv      transcript per gene (splice.py)

Coordinates are one-based inclusive. Undefined values are written as NA and
empty strings as empty fields. Every file written to the output directory is
listed with its SHA-256 and row count in comparison_manifest.json ("outputs").
"""

import gzip
import hashlib
import os
from pathlib import Path

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise import contract, splice
from ensembl.genes.annotation_qc.metrics.pairwise.evaluation_audit import (
    evaluated_transcripts_table,
)

METRICS_TSV = "metrics.tsv"
EVALUATED_TRANSCRIPTS = "evaluated_transcripts.tsv.gz"
PER_SEQUENCE_TSV = "per_sequence.tsv"
SPLICE_POSITION_TSV = "splice_sites_by_position.tsv"
INTRON_LENGTH_TSV = "cds_introns_by_length.tsv"
STOP_AUDIT_TSV = "stop_codon_harmonisation.tsv"

# Files whose layout the contract defines, with their schema versions.
SCHEMAS = {
    "comparison_summary.json": "summary/1.0",
    "comparison_summary.tsv": "summary_tsv/1.0",
    "comparison_details.tsv": "details/1.0",
    "consensus_transcript_labels.tsv": "transcript_labels/1.0",
    "gene_splits.tsv": "split_merge/1.0",
    "gene_merges.tsv": "split_merge/1.0",
    "cds_exact_one_to_one_pairs.tsv": "one_to_one_pairs/1.0",
    "reference_filter_audit.json": "filter_audit/1.0",
    "reference_filter_audit.tsv": "filter_audit/1.0",
    METRICS_TSV: "metrics/1.0",
    EVALUATED_TRANSCRIPTS: "evaluated_transcripts/1.0",
    PER_SEQUENCE_TSV: "per_sequence/1.0",
    SPLICE_POSITION_TSV: "splice_position/1.0",
    INTRON_LENGTH_TSV: "intron_length/1.0",
    "stop_codon_harmonisation.tsv": "stop_codon_audit/1.0",
    "query_evaluated.stop_harmonised.gtf": "gtf",
    "evidence_attribution_labeled.tsv": "evidence_attribution/1.0",
}


def _tsv(frame: pd.DataFrame, path: str, compress: bool = False) -> str:
    if compress:
        # mtime=0 keeps the gzip bytes, and so the checksum, reproducible.
        with (
            open(path, "wb") as raw,
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as handle,
        ):
            handle.write(frame.to_csv(sep="\t", index=False, na_rep="NA").encode())
    else:
        frame.to_csv(path, sep="\t", index=False, na_rep="NA")
    return path


def write_metrics(table: pd.DataFrame, outdir: str) -> str:
    """Write metrics.tsv; values unrounded, NA when the status is not ok."""
    return _tsv(table, os.path.join(outdir, METRICS_TSV))


def write_evaluated_transcripts(table: pd.DataFrame, outdir: str) -> str:
    """Write evaluated_transcripts.tsv.gz (both roles, every parsed transcript)."""
    return _tsv(table, os.path.join(outdir, EVALUATED_TRANSCRIPTS), compress=True)


def per_sequence_table(
    ref: pd.DataFrame, query: pd.DataFrame, pairs: pd.DataFrame
) -> pd.DataFrame:
    """Per sequence: reference/query genes, recovered and exact genes, one-to-one TP."""
    ref_counts = pd.DataFrame(
        {
            "reference_genes": ref.groupby("chrom").size(),
            "reference_cds_locus_recovered": ref.groupby("chrom")[
                "cds_locus_recovered"
            ].sum(),
            "reference_cds_coordinate_exact": ref.groupby("chrom")[
                "cds_coordinate_exact"
            ].sum(),
        }
    )
    query_counts = pd.DataFrame(
        {
            "query_genes": query.groupby("chrom").size(),
            "query_cds_coordinate_exact": query.groupby("chrom")[
                "cds_coordinate_exact"
            ].sum(),
        }
    )
    table = ref_counts.join(query_counts, how="outer")
    table["one_to_one_tp"] = pairs.groupby("chrom").size()
    table = table.fillna(0).astype("int64").reset_index(names="chrom")
    return table


def write_per_sequence(table: pd.DataFrame, outdir: str) -> str:
    """Write per_sequence.tsv."""
    return _tsv(table, os.path.join(outdir, PER_SEQUENCE_TSV))


def splice_tables(results: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flatten splice.analyse output into the position and length tables."""
    position, length = [], []
    for scope_id, result in results.items():
        position += [{"scope_id": scope_id, **row} for row in result["position"]]
        length += [{"scope_id": scope_id, **row} for row in result["length"]]
    length_table = pd.DataFrame(length)
    length_table.columns = [
        {"lo": "lo_bp", "hi": "hi_bp"}.get(c, c) for c in length_table.columns
    ]
    length_table["hi_bp"] = length_table["hi_bp"].astype("Int64")
    return pd.DataFrame(position), length_table


def write_splice_tables(results: dict, outdir: str) -> list[str]:
    """Write splice_sites_by_position.tsv and cds_introns_by_length.tsv."""
    position, length = splice_tables(results)
    return [
        _tsv(position, os.path.join(outdir, SPLICE_POSITION_TSV)),
        _tsv(length, os.path.join(outdir, INTRON_LENGTH_TSV)),
    ]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _rows(path: Path) -> int | None:
    if path.suffix == ".json":
        return None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
        lines = sum(1 for line in handle if line.strip() and not line.startswith("#"))
    return lines - (0 if path.name.endswith(".gtf") else 1)


def outputs_inventory(outdir: str, exclude: tuple = ()) -> dict:
    """SHA-256, size, data rows and schema of every file in the output directory."""
    inventory = {}
    for path in sorted(Path(outdir).iterdir()):
        if not path.is_file() or path.name in exclude:
            continue
        inventory[path.name] = {
            "sha256": _sha256(path),
            "bytes": path.stat().st_size,
            "rows": _rows(path),
            "schema": SCHEMAS.get(path.name),
        }
    return inventory


def write_contract_outputs(  # pylint: disable=too-many-arguments,too-many-locals
    options: dict,
    data: dict,
    models: dict,
    results: dict,
    pairs,
    scopes: dict,
    selection: dict,
    stop: dict,
    input_records: dict,
) -> dict:
    """
    Compute and write the files added by output contract 1.0.
    Args:
            options: outdir, export_evaluated_models, splice_available,
                    canonical_reference
            data: parsed/compared tables, exclusion trackers and seqname map
            models, results: GeneModels and classify_loci results per role
            pairs: One-to-one pair table
            scopes: analysis scope id -> sequence names
            selection: per role, the selection option and selected_by labels
            stop: stop-codon audit (None when not applied) and transformed ids
            input_records: manifest records of the reference and query files
    Returns:
            The manifest "analyses" block
    """
    ref, query = results["reference"], results["query"]
    tables = [
        evaluated_transcripts_table(
            role,
            data["parsed"][role],
            data["compared"][role],
            data["trackers"][role],
            selection[role]["mode"],
            selection[role]["labels"],
            data["mapping"],
            (input_records[role] or {}).get("sha256"),
            transformed=stop["transformed"] if role == "query" else None,
            transformation_audit=(
                STOP_AUDIT_TSV if role == "query" and stop["audit"] is not None else ""
            ),
            scopes=scopes,
            full=options["export_evaluated_models"] == "full",
        )
        for role in ("reference", "query")
    ]
    transcripts = pd.concat(tables, ignore_index=True)
    write_evaluated_transcripts(transcripts, options["outdir"])

    edges = contract.exact_edges(models["reference"], models["query"])
    rows = []
    for scope_id, sequences in {contract.ALL_SCOPE: None, **scopes}.items():
        rows += contract.scope_metric_rows(
            scope_id,
            sequences,
            ref,
            query,
            models["reference"],
            models["query"],
            pairs,
            edges,
            transcripts,
            stop["audit"],
            options["canonical_reference"],
        )
    write_metrics(contract.metric_table(rows), options["outdir"])
    write_per_sequence(per_sequence_table(ref, query, pairs), options["outdir"])

    analyses = {
        "core": {"status": "ok"},
        "one_to_one_exact_cds": {"status": "ok"},
        "evaluated_transcripts": {
            "status": "ok",
            "mode": options["export_evaluated_models"],
        },
        "stop_codon_harmonisation": (
            {"status": "ok"}
            if stop["audit"] is not None
            else {"status": "not_applicable", "reason": "transformation_not_applied"}
        ),
    }
    if options["splice_available"]:
        write_splice_tables(
            contract.splice_analysis(models["reference"], models["query"], ref, scopes),
            options["outdir"],
        )
        analyses["splice_position_length"] = {
            "status": "ok",
            "position_bins": splice.POSITION_BINS,
            "position_coordinates": list(splice.POSITION_COORDINATES),
            "position_min_intron_lengths": list(splice.POSITION_MIN_LENGTHS),
            "intron_length_edges": list(splice.INTRON_LENGTH_EDGES),
            "support_rule": splice.SUPPORT_RULE,
            "detection_rule": splice.DETECTION_RULE,
        }
    else:
        analyses["splice_position_length"] = {
            "status": "not_applicable",
            "reason": "requires_one_reference_transcript_per_gene",
        }
    return analyses

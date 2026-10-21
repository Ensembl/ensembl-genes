"""
Rebuild the exact gene models a comparison evaluated, for locus display.

Follows the comparator's documented procedure (output contract 1.0, "Rebuilding
the evaluated models"); transcripts are never selected again:

1. the input file must exist and its SHA-256 must equal ``inputs.<role>.sha256``
   in the manifest (and ``input_sha256`` of the export rows); otherwise the
   models are unavailable, never substituted by another file;
2. it is parsed with the comparator parser (``parse_annotation_for_comparison``,
   the format option the comparison used) and the recorded seqname map /
   assembly report are applied;
3. the transcripts listed in ``evaluated_transcripts.tsv.gz`` for that role are
   kept, with their per-policy status (evaluated or excluded, with the reason and
   ``selected_by``) taken from the export;
4. transcripts marked ``transformed = stop_codon_added`` get the new terminal CDS
   segment (and the inferred exon) from ``stop_codon_harmonisation.tsv``; the
   transformation is never inferred from a transformed file;
5. exon count, CDS segment count and CDS length of every evaluated transcript
   are checked against the export. A mismatch makes the models unavailable.

Coordinates in the returned rows are one-based inclusive.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ensembl.genes.annotation_qc.parsers.annotation import (
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.seqnames import (
    apply_seqname_mapping,
    build_seqname_mapping,
)

STOP_TRANSFORMATION = "stop_codon_added"


class ModelsUnavailable(RuntimeError):
    """The evaluated models cannot be rebuilt exactly; the message says why."""


def _input_record(manifest: dict, role: str) -> dict:
    record = (manifest.get("inputs") or {}).get(role)
    if not record or not record.get("path") or not record.get("sha256"):
        raise ModelsUnavailable(
            f"the manifest records no {role} input path and checksum"
        )
    return record


def parse_input(manifest: dict, role: str, checksums) -> pd.DataFrame:
    """Steps 1 and 2: checksum-verified parse with the recorded seqname mapping."""
    record = _input_record(manifest, role)
    path = Path(record["path"])
    if not path.exists():
        raise ModelsUnavailable(f"{role} input not found: {path}")
    actual = checksums.sha256(path)
    if actual != record["sha256"]:
        raise ModelsUnavailable(
            f"{role} input {path} has sha256 {actual[:12]}…, the comparison used "
            f"{record['sha256'][:12]}…; models not shown"
        )
    options = manifest.get("options") or {}
    fmt = options.get(f"{role}_format") or "auto"
    inputs = manifest.get("inputs") or {}
    seqmap = (inputs.get("seqname_map") or {}).get("path")
    report = (inputs.get("assembly_report") or {}).get("path")
    for label, path_text in (("seqname map", seqmap), ("assembly report", report)):
        if path_text and not Path(path_text).exists():
            raise ModelsUnavailable(f"recorded {label} not found: {path_text}")
    df = parse_annotation_for_comparison(str(path), fmt)
    mapping = build_seqname_mapping(seqmap, report)
    df, _ = apply_seqname_mapping(df, mapping)
    return df


def _intervals(rows: pd.DataFrame) -> pd.Series:
    """transcript_id -> 'start-end,...' (one-based inclusive, sorted)."""
    if rows.empty:
        return pd.Series(dtype=object)
    rows = rows.sort_values(["transcript_id", "Start", "End"], kind="stable")
    text = (rows.Start + 1).astype(str) + "-" + rows.End.astype(str)
    return text.groupby(rows.transcript_id, sort=False).agg(",".join)


def apply_stop_audit(features: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    """Step 4: move the terminal CDS (and an inferred exon) to the audited coordinates."""
    features = features.reset_index(drop=True)
    extended = audit[audit.status == "extended"]
    if extended.empty:
        return features
    moves = pd.DataFrame(
        {
            "transcript_id": extended.transcript_id.values,
            "Start": extended.terminal_cds_start.astype(int).values - 1,
            "End": extended.terminal_cds_end.astype(int).values,
            "new_start": extended.new_terminal_cds_start.astype(int).values - 1,
            "new_end": extended.new_terminal_cds_end.astype(int).values,
        }
    )
    eligible = (features.Feature == "CDS") | (
        (features.Feature == "exon")
        & features.source_feature.astype(str).str.startswith("inferred_")
    )
    candidates = features.loc[eligible, ["transcript_id", "Feature", "Start", "End"]]
    hits = candidates.reset_index().merge(moves, on=["transcript_id", "Start", "End"])
    found = set(hits.loc[hits.Feature == "CDS", "transcript_id"])
    absent = set(moves.transcript_id) - found
    if absent:
        first = moves[moves.transcript_id == sorted(absent)[0]].iloc[0]
        raise ModelsUnavailable(
            f"audited terminal CDS {first.Start + 1}-{first.End} of {first.transcript_id} "
            f"not found in the input ({len(absent)} transcripts)"
        )
    features = features.copy()
    features.loc[hits["index"].values, "Start"] = hits.new_start.values
    features.loc[hits["index"].values, "End"] = hits.new_end.values
    return features


def rebuild(
    parsed: pd.DataFrame,
    export: pd.DataFrame,
    audit: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Steps 3 to 5 for one role of one comparison.

    Args:
            parsed: parse_input() output for this role
            export: evaluated_transcripts rows of this role (one comparison)
            audit: stop_codon_harmonisation.tsv when the export marks transformed rows
    Returns:
            one row per exported transcript: gene_id, transcript_id,
            original_gene_id, original_transcript_id, chrom, start, end, strand,
            biotype, tags, exons, cds, transformed (coordinates one-based)
    Raises:
            ModelsUnavailable when a transformation lacks its audit or a check fails
    """
    wanted = set(export.transcript_id)
    features = parsed[parsed.transcript_id.isin(wanted)]
    absent = wanted - set(features.transcript_id)
    if absent:
        raise ModelsUnavailable(
            f"{len(absent)} exported transcripts not found in the input (e.g. {sorted(absent)[0]})"
        )
    transformed = export[export.transformed != ""]
    if not transformed.empty:
        kinds = set(transformed.transformed)
        if kinds != {STOP_TRANSFORMATION}:
            raise ModelsUnavailable(f"unsupported transformation(s) {sorted(kinds)}")
        if audit is None:
            raise ModelsUnavailable(
                "transcripts are marked stop_codon_added but the audit file is absent"
            )
        audit = audit[audit.transcript_id.isin(set(transformed.transcript_id))]
        features = apply_stop_audit(features, audit)
    exons = features[features.Feature == "exon"]
    cds = features[features.Feature == "CDS"]
    exon_iv, cds_iv = _intervals(exons), _intervals(cds)
    exon_n = exons.groupby("transcript_id").size()
    cds_n = cds.groupby("transcript_id").size()
    cds_len = (cds.End - cds.Start).groupby(cds.transcript_id).sum()
    rows = export.set_index("transcript_id")
    evaluated = rows[rows.status == "evaluated"]
    checks = {
        "exon_count": exon_n,
        "cds_segment_count": cds_n,
        "cds_length": cds_len,
    }
    for column, observed in checks.items():
        got = observed.reindex(evaluated.index).fillna(0).astype(int)
        bad = got != evaluated[column].astype(int)
        if bad.any():
            first = bad[bad].index[0]
            raise ModelsUnavailable(
                f"{int(bad.sum())} evaluated transcripts differ in {column} from the export "
                f"(e.g. {first}: rebuilt {got[first]}, exported {evaluated.at[first, column]})"
            )
    tx = features[features.Feature == "transcript"].drop_duplicates("transcript_id")
    tx = tx.set_index("transcript_id")
    out = pd.DataFrame(
        {
            "transcript_id": rows.index,
            "gene_id": rows.gene_id.values,
            "original_gene_id": rows.original_gene_id.values,
            "original_transcript_id": rows.original_transcript_id.values,
            "chrom": rows.chrom.values,
            "start": rows.start.astype(int).values,
            "end": rows.end.astype(int).values,
            "strand": rows.strand.values,
            "biotype": rows.transcript_biotype.values,
            "tags": tx.tags.reindex(rows.index).fillna("").values,
            "exons": exon_iv.reindex(rows.index).fillna("").values,
            "cds": cds_iv.reindex(rows.index).fillna("").values,
            "transformed": rows.transformed.values,
        }
    )
    return out


def evaluated_status(export: pd.DataFrame) -> pd.DataFrame:
    """Per-transcript status columns of one comparison's export (one role)."""
    return export[
        [
            "transcript_id",
            "gene_id",
            "status",
            "exclusion_reason",
            "selected_by",
            "exon_count",
            "cds_segment_count",
            "cds_length",
            "transformed",
        ]
    ].copy()

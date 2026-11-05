"""Feature-level annotation metrics.

These functions operate on the standardized PyRanges/DataFrame produced by
``parsers.annotation``. They return tables rather than printing or writing
files, so the same measurements can feed reports and later event analyses.
"""

from __future__ import annotations

import pandas as pd

from ensembl.genes.annotation_qc.metrics.events import derive_introns
from ensembl.genes.annotation_qc.parsers.annotation import standardize_annotation

FIVE_PRIME_UTR_TYPES = frozenset(
    {"five_prime_UTR", "five_prime_utr", "5_prime_utr", "5_utr", "5utr"}
)
THREE_PRIME_UTR_TYPES = frozenset(
    {"three_prime_UTR", "three_prime_utr", "3_prime_utr", "3_utr", "3utr"}
)
CODING_PHASES = frozenset({"0", "1", "2", 0, 1, 2})


def _interval_overlap(
    left_start: int, left_end: int, right_start: int, right_end: int
) -> int:
    return max(0, min(left_end, right_end) - max(left_start, right_start))


def _overlap_with_utrs(cds: pd.DataFrame, utr: pd.DataFrame) -> int:
    if cds.empty or utr.empty:
        return 0
    return sum(
        _interval_overlap(
            int(cds_row.Start),
            int(cds_row.End),
            int(utr_row.Start),
            int(utr_row.End),
        )
        for cds_row in cds.itertuples(index=False)
        for utr_row in utr.itertuples(index=False)
    )


def compute_feature_metrics(annotation) -> pd.DataFrame:
    """Compute one structural feature row per transcript."""
    df = standardize_annotation(annotation)
    transcripts = df[
        df["Feature"].isin({"transcript", "mRNA"}) & df["transcript_id"].notna()
    ].copy()
    if transcripts.empty:
        return pd.DataFrame(
            columns=[
                "feature_type",
                "feature_id",
                "gene_id",
                "exon_count",
                "exonic_length_bp",
                "intron_count",
                "five_prime_utr_length_bp",
                "three_prime_utr_length_bp",
                "cds_segment_count",
                "cds_length_nt",
                "phase_missing_count",
                "cds_frame_ok",
                "coordinate_valid",
                "child_order_valid",
                "cds_utr_overlap_bp",
            ]
        )

    exons = df[(df["Feature"] == "exon") & df["transcript_id"].notna()]
    cds = df[(df["Feature"] == "CDS") & df["transcript_id"].notna()]
    utr5 = df[df["Feature"].isin(FIVE_PRIME_UTR_TYPES) & df["transcript_id"].notna()]
    utr3 = df[df["Feature"].isin(THREE_PRIME_UTR_TYPES) & df["transcript_id"].notna()]
    introns = derive_introns(df)
    rows = []

    for transcript in transcripts.itertuples(index=False):
        transcript_id = str(transcript.transcript_id)
        transcript_exons = exons[exons["transcript_id"] == transcript_id]
        transcript_cds = cds[cds["transcript_id"] == transcript_id]
        transcript_utr5 = utr5[utr5["transcript_id"] == transcript_id]
        transcript_utr3 = utr3[utr3["transcript_id"] == transcript_id]
        transcript_introns = introns[introns["transcript_id"] == transcript_id]
        exon_order = transcript_exons.sort_values(
            "Start", ascending=str(transcript.Strand) != "-"
        )
        child_order_valid = True
        previous_end = None
        for exon in exon_order.itertuples(index=False):
            if int(exon.Start) >= int(exon.End):
                child_order_valid = False
            if previous_end is not None and int(exon.Start) < previous_end:
                child_order_valid = False
            previous_end = max(previous_end or int(exon.End), int(exon.End))

        phases = transcript_cds["phase"].dropna()
        phase_missing = int(
            len(transcript_cds) - len(phases[phases.isin(CODING_PHASES)])
        )
        cds_length = int(
            (
                transcript_cds["End"].astype(int) - transcript_cds["Start"].astype(int)
            ).sum()
        )
        rows.append(
            {
                "feature_type": "transcript",
                "feature_id": transcript_id,
                "gene_id": transcript.gene_id,
                "exon_count": len(transcript_exons),
                "exonic_length_bp": int(
                    (
                        transcript_exons["End"].astype(int)
                        - transcript_exons["Start"].astype(int)
                    ).sum()
                ),
                "intron_count": len(transcript_introns),
                "five_prime_utr_length_bp": int(
                    (
                        transcript_utr5["End"].astype(int)
                        - transcript_utr5["Start"].astype(int)
                    ).sum()
                ),
                "three_prime_utr_length_bp": int(
                    (
                        transcript_utr3["End"].astype(int)
                        - transcript_utr3["Start"].astype(int)
                    ).sum()
                ),
                "cds_segment_count": len(transcript_cds),
                "cds_length_nt": cds_length,
                "phase_missing_count": phase_missing,
                "cds_frame_ok": pd.NA if transcript_cds.empty else cds_length % 3 == 0,
                "coordinate_valid": bool(
                    int(transcript.Start) < int(transcript.End)
                    and str(transcript.Strand) in {"+", "-", "."}
                ),
                "child_order_valid": child_order_valid,
                "cds_utr_overlap_bp": _overlap_with_utrs(
                    transcript_cds, pd.concat([transcript_utr5, transcript_utr3])
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_feature_metrics(feature_metrics: pd.DataFrame) -> pd.DataFrame:
    """Return count, mean, and rate summaries from transcript feature rows."""
    if feature_metrics.empty:
        return pd.DataFrame(columns=["metric", "value"])

    total = len(feature_metrics)
    rows = [
        {"metric": "transcript_count", "value": total},
        {"metric": "mean_exon_count", "value": feature_metrics["exon_count"].mean()},
        {
            "metric": "mean_exonic_length_bp",
            "value": feature_metrics["exonic_length_bp"].mean(),
        },
        {
            "metric": "mean_intron_count",
            "value": feature_metrics["intron_count"].mean(),
        },
        {
            "metric": "mean_cds_length_nt",
            "value": feature_metrics["cds_length_nt"].mean(),
        },
        {
            "metric": "frame_consistent_transcripts",
            "value": int(feature_metrics["cds_frame_ok"].fillna(False).sum()),
        },
        {
            "metric": "coordinate_valid_transcripts",
            "value": int(feature_metrics["coordinate_valid"].sum()),
        },
        {
            "metric": "child_order_valid_transcripts",
            "value": int(feature_metrics["child_order_valid"].sum()),
        },
    ]
    return pd.DataFrame(rows)

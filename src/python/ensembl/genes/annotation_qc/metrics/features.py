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
    exon_by_transcript = {
        key: group for key, group in exons.groupby("transcript_id", sort=False)
    }
    cds_by_transcript = {
        key: group for key, group in cds.groupby("transcript_id", sort=False)
    }
    utr5_by_transcript = {
        key: group for key, group in utr5.groupby("transcript_id", sort=False)
    }
    utr3_by_transcript = {
        key: group for key, group in utr3.groupby("transcript_id", sort=False)
    }

    def grouped_count(frame: pd.DataFrame) -> pd.Series:
        return frame.groupby("transcript_id", sort=False).size()

    def grouped_length(frame: pd.DataFrame) -> pd.Series:
        lengths = frame["End"].astype(int) - frame["Start"].astype(int)
        return lengths.groupby(frame["transcript_id"], sort=False).sum()

    exon_counts = grouped_count(exons)
    exonic_lengths = grouped_length(exons)
    intron_counts = grouped_count(introns)
    utr5_lengths = grouped_length(utr5)
    utr3_lengths = grouped_length(utr3)
    cds_counts = grouped_count(cds)
    cds_lengths = grouped_length(cds)
    phase_missing_counts = (
        (~cds["phase"].isin(CODING_PHASES))
        .groupby(cds["transcript_id"], sort=False)
        .sum()
    )
    coordinate_valid = (
        transcripts["Start"].astype(int) < transcripts["End"].astype(int)
    ) & transcripts["Strand"].isin({"+", "-", "."})
    coordinate_valid_by_transcript = dict(
        zip(transcripts["transcript_id"], coordinate_valid, strict=False)
    )
    rows = []

    for transcript in transcripts.itertuples(index=False):
        transcript_id = str(transcript.transcript_id)
        transcript_exons = exon_by_transcript.get(transcript_id, exons.iloc[0:0])
        transcript_cds = cds_by_transcript.get(transcript_id, cds.iloc[0:0])
        transcript_utr5 = utr5_by_transcript.get(transcript_id, utr5.iloc[0:0])
        transcript_utr3 = utr3_by_transcript.get(transcript_id, utr3.iloc[0:0])
        exon_order = transcript_exons.sort_values(
            "Start", ascending=str(transcript.Strand) != "-"
        )
        child_order_valid = True
        previous_start = None
        previous_end = None
        for exon_start, exon_end in exon_order[["Start", "End"]].itertuples(
            index=False, name=None
        ):
            exon_start = int(exon_start)
            exon_end = int(exon_end)
            if exon_start >= exon_end:
                child_order_valid = False
            if previous_start is not None and str(transcript.Strand) == "-":
                child_order_valid &= exon_end <= previous_start
            elif previous_end is not None:
                child_order_valid &= exon_start >= previous_end
            previous_start = exon_start
            previous_end = exon_end

        cds_length = int(cds_lengths.get(transcript_id, 0))
        phase_missing = int(phase_missing_counts.get(transcript_id, 0))
        rows.append(
            {
                "feature_type": "transcript",
                "feature_id": transcript_id,
                "gene_id": transcript.gene_id,
                "exon_count": int(exon_counts.get(transcript_id, 0)),
                "exonic_length_bp": int(exonic_lengths.get(transcript_id, 0)),
                "intron_count": int(intron_counts.get(transcript_id, 0)),
                "five_prime_utr_length_bp": int(utr5_lengths.get(transcript_id, 0)),
                "three_prime_utr_length_bp": int(utr3_lengths.get(transcript_id, 0)),
                "cds_segment_count": int(cds_counts.get(transcript_id, 0)),
                "cds_length_nt": cds_length,
                "phase_missing_count": phase_missing,
                "cds_frame_ok": (
                    pd.NA
                    if transcript_id not in cds_counts.index
                    else cds_length % 3 == 0
                ),
                "coordinate_valid": bool(
                    coordinate_valid_by_transcript.get(transcript_id, False)
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

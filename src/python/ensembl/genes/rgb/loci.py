from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import pandas as pd


@dataclass
class Interval:
    start: int
    end: int

    def overlaps(self, other: "Interval") -> bool:
        return not (self.end < other.start or other.end < self.start)

    def gap_to(self, other: "Interval") -> int:
        # distance between inclusive intervals; if overlapping, negative/zero
        return other.start - self.end - 1


@dataclass
class Locus:
    seq: str
    strand: int
    start: int
    end: int
    index: int

    @property
    def locus_id(self) -> str:
        return f"{self.seq}:{self.strand}:{self.start}:{self.end}:{self.index}"


def _merge_overlaps(intervals: List[Interval]) -> List[Interval]:
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x.start, x.end))
    merged = [intervals[0]]
    for it in intervals[1:]:
        last = merged[-1]
        if last.overlaps(it) or it.start <= last.end:  # defensive
            last.end = max(last.end, it.end)
        else:
            merged.append(Interval(it.start, it.end))
    return merged


def _merge_with_gap(intervals: List[Interval], gap_bp: int) -> List[Interval]:
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x.start, x.end))
    merged = [Interval(intervals[0].start, intervals[0].end)]
    for it in intervals[1:]:
        last = merged[-1]
        gap = last.gap_to(it)
        if gap <= gap_bp:
            last.end = max(last.end, it.end)
        else:
            merged.append(Interval(it.start, it.end))
    return merged


def build_loci(
    core_df: pd.DataFrame,
    layer_df: pd.DataFrame,
    gap_bp: int = 5000,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build strict and expanded loci from union of core/layer genes.

    Returns (loci_strict_df, loci_expanded_df, gene_to_locus_df).
    """
    needed_cols = [
        "gene_id",
        "seq_region_name",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
    ]
    for name, df in ("core", core_df), ("layer", layer_df):
        missing = [c for c in needed_cols if c not in df.columns]
        if missing:
            raise ValueError(f"{name} genes missing columns: {missing}")

    # normalize dtypes
    def norm(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["seq_region_strand"] = out["seq_region_strand"].astype(int)
        out["seq_region_start"] = out["seq_region_start"].astype(int)
        out["seq_region_end"] = out["seq_region_end"].astype(int)
        return out

    core = norm(core_df)
    layer = norm(layer_df)

    # Prepare union intervals by (seq, strand)
    union_key_to_intervals: Dict[Tuple[str, int], List[Interval]] = {}
    for df in (core, layer):
        for (seq, strand), g in df.groupby(["seq_region_name", "seq_region_strand"]):
            lst = union_key_to_intervals.setdefault((seq, strand), [])
            lst += [Interval(int(s), int(e)) for s, e in zip(g.seq_region_start, g.seq_region_end)]

    strict_records: List[Tuple[str, int, int, int, int, int, str]] = []
    expanded_records: List[Tuple[str, int, int, int, int, int, str]] = []
    # Map genes → loci
    gene_map_records: List[Tuple[str, int, str, int, int, str, str]] = []

    for (seq, strand), intervals in union_key_to_intervals.items():
        strict_list = _merge_overlaps([Interval(i.start, i.end) for i in intervals])
        expanded_list = _merge_with_gap([Interval(i.start, i.end) for i in strict_list], gap_bp)

        # Build strict locus rows with counts
        # Pre-slice core/layer genes for this key
        core_key = core[(core.seq_region_name == seq) & (core.seq_region_strand == strand)]
        layer_key = layer[(layer.seq_region_name == seq) & (layer.seq_region_strand == strand)]

        # Strict loci
        strict_ids: List[Locus] = [Locus(seq, strand, it.start, it.end, i) for i, it in enumerate(strict_list)]
        for loc in strict_ids:
            c_count = ((core_key.seq_region_start <= loc.end) & (core_key.seq_region_end >= loc.start)).sum()
            l_count = ((layer_key.seq_region_start <= loc.end) & (layer_key.seq_region_end >= loc.start)).sum()
            strict_records.append(
                (
                    seq,
                    strand,
                    loc.start,
                    loc.end,
                    loc.end - loc.start + 1,
                    int(c_count),
                    int(l_count),
                    loc.locus_id,
                )
            )

        # Expanded loci
        expanded_ids: List[Locus] = [Locus(seq, strand, it.start, it.end, i) for i, it in enumerate(expanded_list)]
        for loc in expanded_ids:
            c_count = ((core_key.seq_region_start <= loc.end) & (core_key.seq_region_end >= loc.start)).sum()
            l_count = ((layer_key.seq_region_start <= loc.end) & (layer_key.seq_region_end >= loc.start)).sum()
            expanded_records.append(
                (
                    seq,
                    strand,
                    loc.start,
                    loc.end,
                    loc.end - loc.start + 1,
                    int(c_count),
                    int(l_count),
                    loc.locus_id,
                )
            )

        # Map each gene to its strict and expanded locus
        for _, row in pd.concat([core_key.assign(db_kind="core"), layer_key.assign(db_kind="layer")]).iterrows():
            iv = Interval(int(row.seq_region_start), int(row.seq_region_end))

            def find_locus(loci: List[Locus]) -> str:
                for loc in loci:
                    if not (iv.end < loc.start or iv.start > loc.end):
                        return loc.locus_id
                return ""

            strict_id = find_locus(strict_ids)
            expanded_id = find_locus(expanded_ids)
            gene_map_records.append(
                (
                    str(row.db_kind),
                    int(row.gene_id),
                    str(seq),
                    int(strand),
                    int(row.seq_region_start),
                    strict_id,
                    expanded_id,
                )
            )

    loci_cols = [
        "seq_region_name",
        "seq_region_strand",
        "locus_start",
        "locus_end",
        "locus_length",
        "core_gene_count",
        "layer_gene_count",
        "locus_id",
    ]
    strict_df = pd.DataFrame(strict_records, columns=loci_cols)
    expanded_df = pd.DataFrame(expanded_records, columns=loci_cols)

    map_df = pd.DataFrame(
        gene_map_records,
        columns=[
            "db_kind",
            "gene_id",
            "seq_region_name",
            "seq_region_strand",
            "gene_start",
            "locus_id_strict",
            "locus_id_expanded",
        ],
    )
    return strict_df, expanded_df, map_df

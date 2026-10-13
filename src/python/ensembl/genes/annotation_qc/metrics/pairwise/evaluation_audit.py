"""
Which transcripts a comparison evaluated, and why the others were left out.

Every transcript parsed from the reference or the query gets one row in
evaluated_transcripts.tsv.gz (output contract 1.0), with its status:

    evaluated   the transcript is in the gene models that were classified and
                matched
    excluded    exclusion_reason names the first step that removed it:
                    excluded_seqname            --genome-mismatch exclude
                    evaluation_mode:<mode>      reference evaluation preset
                    gene_biotype_filter         --reference-gene-biotypes
                    transcript_biotype_filter   --reference-transcript-biotypes
                    not_selected:<selection>    transcript selection policy
                    outside_regions             --region / --regions-file
                    no_exon_records             no exon (or CDS) row to compare

Evaluated transcripts also record how they were chosen (selected_by):

    all                              every transcript is kept
    longest_cds                      longest total CDS in its gene
    longest_exon_length              gene without CDS: longest total exon length
    canonical_tag                    the Ensembl_canonical-tagged transcript
    canonical_fallback_longest_cds   gene without a tagged transcript: longest
                                     CDS (or exon length) instead

The rules mirror selection.select_transcripts, which made the choice; they are
recomputed here from the same table only to label it. Together with the
comparator's own identifiers (gene_id/transcript_id, namespaced as
"<seqname>:<id>" when the file reuses an identifier on several sequences) and the
file identifiers, a consumer can rebuild exactly the evaluated models from the
input file without selecting transcripts again.
"""

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise.classify import GeneModels

CANONICAL_TAG = "Ensembl_canonical"

EXPORT_COLUMNS = [
    "role",
    "gene_id",
    "transcript_id",
    "original_gene_id",
    "original_transcript_id",
    "chrom",
    "source_seqname",
    "strand",
    "start",
    "end",
    "gene_biotype",
    "transcript_biotype",
    "source_feature",
    "record_origin",
    "canonical_tag",
    "selection",
    "status",
    "exclusion_reason",
    "selected_by",
    "exon_count",
    "cds_segment_count",
    "cds_length",
    "transformed",
    "transformation_audit",
    "input",
    "input_sha256",
    "scope_ids",
]
FULL_COLUMNS = ["exon_segments", "cds_segments"]


def transcript_ids(df: pd.DataFrame) -> set:
    """transcript_id of every transcript record."""
    return set(df.loc[df["Feature"] == "transcript", "transcript_id"])


class ExclusionTracker:
    """Record the first step that removes each transcript."""

    def __init__(self, df: pd.DataFrame):
        self.remaining = transcript_ids(df)
        self.reasons: dict[str, str] = {}

    def step(self, df: pd.DataFrame, reason: str) -> None:
        """Mark transcripts absent from df (but present before) with reason."""
        now = transcript_ids(df)
        for tid in self.remaining - now:
            self.reasons[tid] = reason
        self.remaining &= now

    def finish(self, models: GeneModels) -> None:
        """Transcripts that survived filtering but have no exon record."""
        modelled = {t for tids in models.genes["transcript_ids"] for t in tids}
        for tid in self.remaining - modelled:
            self.reasons[tid] = "no_exon_records"
        self.remaining &= modelled


def canonical_tag_count(df: pd.DataFrame) -> int:
    """Number of transcript records tagged Ensembl_canonical."""
    transcripts = df[df["Feature"] == "transcript"]
    return int(transcripts["tags"].str.contains(CANONICAL_TAG, regex=False).sum())


def selection_labels(df: pd.DataFrame, selected: set, mode: str) -> dict[str, str]:
    """
    selected_by for each selected transcript (see module docstring).
    Args:
            df: Comparison-schema table the selection was made from
            selected: transcript_ids kept by the selection
            mode: "all", "longest_cds" or "canonical"
    Returns:
            transcript_id -> selected_by
    """
    if mode == "all":
        return {tid: "all" for tid in selected}
    transcripts = df[df["Feature"] == "transcript"].drop_duplicates("transcript_id")
    lengths = (df["End"] - df["Start"]).rename("length")
    cds_length = lengths[df["Feature"] == "CDS"].groupby(df["transcript_id"]).sum()
    gene_cds = (
        transcripts["transcript_id"]
        .map(cds_length)
        .fillna(0)
        .groupby(transcripts["gene_id"])
        .max()
    )
    tagged = transcripts["tags"].str.contains(CANONICAL_TAG, regex=False)
    gene_tagged = tagged.groupby(transcripts["gene_id"]).any()
    gene_of = dict(zip(transcripts["transcript_id"], transcripts["gene_id"]))
    is_tagged = dict(zip(transcripts["transcript_id"], tagged))
    labels = {}
    for tid in selected:
        gene = gene_of.get(tid)
        longest = "longest_cds" if gene_cds.get(gene, 0) > 0 else "longest_exon_length"
        if mode == "canonical":
            if is_tagged.get(tid):
                labels[tid] = "canonical_tag"
            elif not gene_tagged.get(gene, False):
                labels[tid] = f"canonical_fallback_{longest}"
            else:  # pragma: no cover - select_transcripts never does this
                labels[tid] = "canonical_untagged"
        else:
            labels[tid] = longest
    return labels


def _segments(df: pd.DataFrame, feature: str, with_text: bool) -> pd.DataFrame:
    """Per transcript: count, total length and (optionally) one-based segments."""
    rows = df.loc[df["Feature"] == feature, ["transcript_id", "Start", "End"]]
    rows = rows.sort_values(["transcript_id", "Start", "End"], kind="stable")
    by_tid = rows.groupby("transcript_id", sort=False)
    out = pd.DataFrame(
        {
            "count": by_tid.size(),
            "length": (rows["End"] - rows["Start"])
            .groupby(rows["transcript_id"])
            .sum(),
        }
    )
    if with_text:
        text = (rows["Start"] + 1).astype(str) + "-" + rows["End"].astype(str)
        out["segments"] = text.groupby(rows["transcript_id"], sort=False).agg(",".join)
    return out


def evaluated_transcripts_table(  # pylint: disable=too-many-arguments,too-many-locals
    role: str,
    parsed: pd.DataFrame,
    compared: pd.DataFrame,
    tracker: ExclusionTracker,
    selection: str,
    selected_by: dict[str, str],
    mapping: dict[str, str],
    input_sha256: str | None,
    transformed: dict[str, str] | None = None,
    transformation_audit: str = "",
    scopes: dict[str, set] | None = None,
    full: bool = False,
) -> pd.DataFrame:
    """
    One row per parsed transcript of one annotation (EXPORT_COLUMNS).
    Args:
            role: "reference" or "query"
            parsed: Comparison-schema table as parsed (file sequence names)
            compared: The table after sequence-name mapping and any transformation,
                    before filtering; supplies comparison coordinates
            tracker: ExclusionTracker after finish()
            selection: The transcript selection option for this role
            selected_by: transcript_id -> selected_by label (evaluated transcripts)
            mapping: Seqname map that was applied
            input_sha256: SHA-256 of the input file (links to the manifest)
            transformed: transcript_id -> transformation id (e.g. stop_codon_added)
            transformation_audit: File name of the transformation audit, if any
            scopes: analysis scope id -> sequence names
            full: Also write exon_segments and cds_segments (one-based)
    Returns:
            DataFrame sorted by chrom, start, transcript_id
    """
    transformed = transformed or {}
    source = parsed[parsed["Feature"] == "transcript"].drop_duplicates("transcript_id")
    coords = compared[compared["Feature"] == "transcript"].drop_duplicates(
        "transcript_id"
    )
    coords = coords.set_index("transcript_id")[["Chromosome", "Start", "End"]]
    # Children from the compared table where present (transformed coordinates),
    # otherwise from the parsed table (e.g. sequences excluded before comparison).
    exons = _segments(parsed, "exon", full)
    cds = _segments(parsed, "CDS", full)
    exons.update(_segments(compared, "exon", full))
    cds.update(_segments(compared, "CDS", full))

    table = pd.DataFrame(
        {
            "role": role,
            "gene_id": source["gene_id"].to_numpy(),
            "transcript_id": source["transcript_id"].to_numpy(),
            "original_gene_id": source["original_gene_id"].to_numpy(),
            "original_transcript_id": source["original_transcript_id"].to_numpy(),
            "source_seqname": source["Chromosome"].to_numpy(),
            "strand": source["Strand"].to_numpy(),
            "gene_biotype": source["gene_biotype"].to_numpy(),
            "transcript_biotype": source["transcript_biotype"].to_numpy(),
            "source_feature": source["source_feature"].to_numpy(),
            "canonical_tag": source["tags"]
            .str.contains(CANONICAL_TAG, regex=False)
            .to_numpy(),
        }
    )
    tids = table["transcript_id"]
    table["chrom"] = [
        coords["Chromosome"].get(t, mapping.get(s, s))
        for t, s in zip(tids, table["source_seqname"])
    ]
    parsed_start = dict(zip(source["transcript_id"], source["Start"]))
    parsed_end = dict(zip(source["transcript_id"], source["End"]))
    table["start"] = [int(coords["Start"].get(t, parsed_start[t])) + 1 for t in tids]
    table["end"] = [int(coords["End"].get(t, parsed_end[t])) for t in tids]
    table["record_origin"] = [
        "inferred" if str(f).startswith("inferred_") else "explicit"
        for f in table["source_feature"]
    ]
    table["selection"] = selection
    table["status"] = [
        "excluded" if t in tracker.reasons else "evaluated" for t in tids
    ]
    table["exclusion_reason"] = [tracker.reasons.get(t, "") for t in tids]
    table["selected_by"] = [
        selected_by.get(t, "") if t not in tracker.reasons else "" for t in tids
    ]
    table["exon_count"] = tids.map(exons["count"]).fillna(0).astype("int64")
    table["cds_segment_count"] = tids.map(cds["count"]).fillna(0).astype("int64")
    table["cds_length"] = tids.map(cds["length"]).fillna(0).astype("int64")
    table["transformed"] = [transformed.get(t, "") for t in tids]
    table["transformation_audit"] = [
        transformation_audit if transformed.get(t) else "" for t in tids
    ]
    table["input"] = role
    table["input_sha256"] = input_sha256 or ""
    scopes = scopes or {}
    table["scope_ids"] = [
        ",".join(sid for sid, seqs in scopes.items() if chrom in seqs)
        for chrom in table["chrom"]
    ]
    columns = list(EXPORT_COLUMNS)
    if full:
        table["exon_segments"] = tids.map(exons["segments"]).fillna("")
        table["cds_segments"] = tids.map(cds["segments"]).fillna("")
        columns += FULL_COLUMNS
    return (
        table[columns]
        .sort_values(["chrom", "start", "transcript_id"], kind="stable")
        .reset_index(drop=True)
    )

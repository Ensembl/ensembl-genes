"""
Representative locus plots for a pairwise comparison (optional).

Requires matplotlib (pip install "ensembl-genes[plots]"). Writes qc/<category>_<n>_<chrom>_<start>.png
and qc/index.html. Each index row carries the gene_id, matched_id and one-based
locus that identify the matching row in comparison_details.tsv.
"""

import html
import os
from collections import defaultdict

import pandas as pd

TRACK_COLOURS = {"Reference": "#D9534F", "Query": "#337AB7"}
CDS_HEIGHT, UTR_HEIGHT = 0.4, 0.18


def import_pyplot():
    try:
        import matplotlib  # pylint: disable=import-outside-toplevel

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel
    except ImportError as error:
        raise ImportError(
            "--plots-per-category needs matplotlib: pip install 'ensembl-genes[plots]'"
        ) from error
    return plt


def _draw_track(ax, annotation, y, colour, chrom, start, end, max_transcripts=5):
    in_view = annotation[
        (annotation["Chromosome"] == chrom)
        & (annotation["End"] > start)
        & (annotation["Start"] < end)
    ]
    exons = in_view[in_view["Feature"] == "exon"]
    cds = in_view[in_view["Feature"] == "CDS"]
    for offset, (tid, tx_exons) in enumerate(
        exons.groupby("transcript_id", sort=False)
    ):
        if offset >= max_transcripts:
            break
        tx_cds = cds[cds["transcript_id"] == tid]
        row_y = y - offset * 0.12
        for exon in tx_exons.itertuples():
            height = UTR_HEIGHT if len(tx_cds) else CDS_HEIGHT
            ax.add_patch(
                _rect(
                    exon.Start,
                    exon.End,
                    row_y,
                    height,
                    colour,
                    0.3 if len(tx_cds) else 0.7,
                )
            )
        for part in tx_cds.itertuples():
            ax.add_patch(_rect(part.Start, part.End, row_y, CDS_HEIGHT, colour, 0.7))
        if len(tx_exons) > 1:
            ax.plot(
                [tx_exons["Start"].min(), tx_exons["End"].max()],
                [row_y, row_y],
                color=colour,
                linewidth=1,
                zorder=0,
            )


def _rect(start, end, y, height, colour, alpha):
    from matplotlib import patches  # pylint: disable=import-outside-toplevel

    return patches.Rectangle(
        (start, y - height / 2),
        end - start,
        height,
        linewidth=0.5,
        edgecolor=colour,
        facecolor=colour,
        alpha=alpha,
    )


def _plot_locus(plt, entry, category, tracks, path):
    pad = max(500, int((entry["end"] - entry["start"]) * 0.3))
    view_start, view_end = max(0, entry["start"] - pad), entry["end"] + pad
    fig, ax = plt.subplots(figsize=(16, max(4, 1.5 + len(tracks) * 0.8)))
    names = list(tracks)
    positions = {name: i for i, name in enumerate(reversed(names))}
    for name, annotation in tracks.items():
        _draw_track(
            ax,
            annotation,
            positions[name],
            TRACK_COLOURS.get(name, "#666666"),
            entry["chrom"],
            view_start,
            view_end,
        )
    ax.set_xlim(view_start, view_end)
    ax.set_ylim(-1, len(names))
    ax.set_yticks(list(positions.values()), list(positions))
    ax.axvspan(entry["start"], entry["end"], alpha=0.05, color="blue")
    ax.set_title(
        f"{category}: {entry['chrom']}:{entry['start'] + 1}-{entry['end']}\n"
        f"{entry['source']} {entry['gene_id']}  matched {entry['matched_id'] or '-'}",
        fontsize=10,
    )
    ax.set_xlabel("Genomic position (zero-based)")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_locus_plots(
    ref_results, query_results, reference, query, outdir, per_category=3
):
    """
    Plot the first loci (by position) of each reference class and of Novel query genes.
    Args:
            ref_results, query_results: Results from classify_loci
            reference, query: Comparison-schema DataFrames that were compared
            outdir: Output directory; plots go to outdir/qc
            per_category: Plots per classification
    Returns:
            Number of plots written
    """
    plt = import_pyplot()
    qc_dir = os.path.join(outdir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    by_class = defaultdict(list)
    for source, results in (("reference", ref_results), ("consensus", query_results)):
        for entry in results.to_dict("records"):
            if source == "reference" or entry["classification"] == "Novel":
                by_class[entry["classification"]].append({**entry, "source": source})

    tracks = {"Reference": reference, "Query": query}
    rows = []
    for category, entries in sorted(by_class.items()):
        entries.sort(key=lambda e: (e["chrom"], e["start"], e["gene_id"]))
        for number, entry in enumerate(entries[:per_category], start=1):
            name = (
                f"{category.lower()}_{number}_{entry['chrom']}_{entry['start'] + 1}.png"
            )
            _plot_locus(plt, entry, category, tracks, os.path.join(qc_dir, name))
            rows.append((category, name, entry))
    _write_index(qc_dir, rows)
    return len(rows)


def _write_index(qc_dir, rows):
    body = [
        "<!doctype html><html><head><meta charset='utf-8'><title>Pairwise comparison loci</title>",
        "<style>body{font-family:sans-serif;margin:2em}table{border-collapse:collapse}"
        "td,th{border:1px solid #ccc;padding:6px}img{max-width:480px}</style></head><body>",
        "<h1>Pairwise comparison loci</h1>",
        "<p>Match rows to comparison_details.tsv by source, gene_id and locus "
        "(one-based, inclusive). 'consensus' means the query annotation.</p>",
    ]
    grouped = pd.DataFrame(
        [
            (
                c,
                n,
                e["source"],
                e["gene_id"],
                e["matched_id"],
                e["chrom"],
                e["start"] + 1,
                e["end"],
            )
            for c, n, e in rows
        ],
        columns=[
            "category",
            "file",
            "source",
            "gene_id",
            "matched_id",
            "chrom",
            "start",
            "end",
        ],
    )
    for category, group in grouped.groupby("category", sort=True):
        body.append(
            f"<h2>{html.escape(category)} ({len(group)})</h2><table><tr><th>Plot</th>"
            "<th>source</th><th>gene_id</th><th>matched_id</th><th>locus</th></tr>"
        )
        for row in group.itertuples():
            body.append(
                f"<tr><td><a href='{row.file}'><img src='{row.file}'></a></td>"
                f"<td>{row.source}</td><td>{html.escape(str(row.gene_id))}</td>"
                f"<td>{html.escape(str(row.matched_id)) or '-'}</td>"
                f"<td>{html.escape(str(row.chrom))}:{row.start}-{row.end}</td></tr>"
            )
        body.append("</table>")
    body.append("</body></html>")
    with open(os.path.join(qc_dir, "index.html"), "w") as handle:
        handle.write("\n".join(body))

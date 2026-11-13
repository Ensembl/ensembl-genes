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
"""Format computed confidence summaries and optional QC figures."""
# Plot dependencies and report cross-links load only when figures are requested.
# pylint: disable=import-outside-toplevel
from __future__ import annotations
import json
from pathlib import Path
from ensembl.genes.annotation_qc.reports.high_confidence import write_records


def write_tsv(path, rows, fields):
    """Shared formatter for prepared summary/diagram records."""
    write_records(path, rows, fields)


def write_summary(result, outdir):
    """Write the same interoperable evidence-summary tables as the standalone tool."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        outdir / "transcripts.tsv", result["transcripts"], result["transcript_fields"]
    )
    for name, (rows, fields) in result["tables"].items():
        write_tsv(outdir / (name + ".tsv"), rows, fields)
    write_tsv(outdir / "overview.tsv", [result["overview"]], list(result["overview"]))
    (outdir / "provenance.json").write_text(
        json.dumps(result["provenance"], indent=2) + "\n"
    )


def write_comparison(results, outdir):
    """Combine prepared table rows; keep each model's original denominator."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for result in results:
        write_summary(result, outdir / "models" / result["model"])
    names = sorted(set().union(*(r["tables"] for r in results)))
    for name in names:
        contributors = [r["tables"][name] for r in results if name in r["tables"]]
        write_tsv(
            outdir / (name + ".tsv"),
            [row for data, _ in contributors for row in data],
            contributors[0][1],
        )
    overview = [r["overview"] for r in results]
    write_tsv(outdir / "overview.tsv", overview, list(overview[0]))
    (outdir / "comparison_provenance.json").write_text(
        json.dumps({"inputs": [r["provenance"] for r in results]}, indent=2) + "\n"
    )


# The orchestration/layout keeps related report operations together.
# pylint: disable-next=too-many-locals,too-many-statements
def write_figures(results, outdir, species=None):
    """Render already computed confidence distributions, pass rates and decision paths."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as error:
        raise ValueError(
            "Install ensembl-genes[annotation-qc-plots] to render confidence figures"
        ) from error
    from matplotlib.colors import LinearSegmentedColormap
    from ensembl.genes.annotation_qc.reports.high_confidence_sankey import draw_sankey

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    unit = (
        "Genes"
        if all(r["overview"]["genes"] == r["overview"]["transcripts"] for r in results)
        else "Transcripts"
    )
    suffix = " — " + species.replace("_", " ").capitalize() if species else ""
    classes = ["high", "medium", "low"]
    colors = ["#269585", "#E5AF38", "#CE5D63"]
    fig, ax = plt.subplots(figsize=(12, max(4, len(results) * 0.7 + 2)))
    left = np.zeros(len(results))
    for category, color in zip(classes, colors):
        values = [r["overview"][category + "_percent"] or 0 for r in results]
        bars = ax.barh(
            [r["model"] for r in results],
            values,
            left=left,
            color=color,
            label=category.capitalize(),
        )
        for segment, result in zip(bars, results):
            if segment.get_width() >= 5:
                ax.text(
                    segment.get_x() + segment.get_width() / 2,
                    segment.get_y() + segment.get_height() / 2,
                    f"{segment.get_width():.2f}%\n(n={result['overview'][category]:,})",
                    ha="center",
                    va="center",
                    fontsize=9,
                )
        left += values
    ax.set(
        xlim=(0, 100),
        xlabel=unit + " assessed (%)",
        title="Confidence distribution" + suffix,
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, frameon=False)
    fig.tight_layout()
    for extension in ["png", "svg"]:
        fig.savefig(
            outdir / f"confidence_distribution.{extension}",
            dpi=180,
            bbox_inches="tight",
        )
    plt.close(fig)
    names = {
        "transcript_cds_complete": "Complete CDS",
        "transcript_five_prime_utr": "5′ UTR",
        "transcript_three_prime_utr": "3′ UTR",
        "transcript_canonical_splicing": "Canonical splicing",
        "transcript_regular_intron_size": "Regular intron size",
        "transcript_intron_support": "Intron support",
    }
    flags = [f for f in names if any(f in r["plot_flags"] for r in results)]
    matrix = np.full((len(results), len(flags)), np.nan)
    for i, result in enumerate(results):
        lookup = {
            r["flag"]: r
            for r in result["tables"]["flag_distribution"][0]
            if r["group"] == "all" and r["state"] == "pass"
        }
        for j, flag in enumerate(flags):
            if flag in result["plot_flags"] and lookup[flag]["denominator"]:
                matrix[i, j] = lookup[flag]["percent"]
    cmap = LinearSegmentedColormap.from_list(
        "qc_teal", ["#EEF3F5", "#B7CCD4", "#648D9D", "#20566A"]
    )
    cmap.set_bad("#E4E7EA")
    fig, ax = plt.subplots(figsize=(14, 8))
    image = ax.imshow(
        np.ma.masked_invalid(matrix), vmin=0, vmax=100, cmap=cmap, aspect="auto"
    )
    for i, j in np.ndindex(matrix.shape):
        value = matrix[i, j]
        label = (
            "N/A"
            if np.isnan(value)
            else (
                ">99.99%"
                if 99.99 < value < 100
                else "<0.01%" if 0 < value < 0.01 else f"{value:.2f}%"
            )
        )
        ax.text(
            j,
            i,
            label,
            ha="center",
            va="center",
            color="white" if value >= 80 else "#19303E",
            fontsize=12,
        )
    ax.set_xticks(range(len(flags)), [names[f].replace(" ", "\n", 1) for f in flags])
    ax.set_yticks(range(len(results)), [r["model"] for r in results])
    ax.set_xticks(np.arange(-0.5, len(flags), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(results), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    ax.set_title("Evidence flags by method" + suffix, pad=20, fontsize=18)
    fig.colorbar(image, ax=ax, label=unit + " passing (%)")
    fig.tight_layout()
    for extension in ["png", "svg"]:
        fig.savefig(
            outdir / f"evidence_flags_heatmap.{extension}", dpi=180, bbox_inches="tight"
        )
    plt.close(fig)
    draw_sankey(results, outdir, species=species)

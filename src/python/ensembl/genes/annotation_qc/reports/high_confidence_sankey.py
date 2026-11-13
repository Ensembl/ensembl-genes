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
"""Render prepared confidence flow counts; biological rules are evaluated by metrics."""
# Plot dependencies and report cross-links load only when figures are requested.
# pylint: disable=import-outside-toplevel
from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path
import textwrap
from ensembl.genes.annotation_qc.config.high_confidence import RULES, CLASSES, COLORS


# The orchestration/layout keeps related report operations together.
# pylint: disable-next=too-many-branches,too-many-locals,too-many-statements
def draw_matched_rule_sankey(results, outdir, species=None):
    """Render model flows through matched rules to confidence classes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch, Rectangle
    from ensembl.genes.annotation_qc.reports.high_confidence_summary import write_tsv

    outdir = Path(outdir)
    counts = Counter()
    for result in results:
        counts.update(result["sankey_routes"])
    write_tsv(
        outdir / "confidence_sankey_matched_rules_flows.tsv",
        [
            {"model": m, "rule": r, "confidence": c, "count": n}
            for (m, r, c), n in sorted(counts.items())
        ],
        ["model", "rule", "confidence", "count"],
    )
    total = sum(counts.values())
    if not total:
        return
    methods = [r["model"] for r in results]
    rules = sorted(
        {(r, c) for m, r, c in counts}, key=lambda key: (CLASSES.index(key[1]), key[0])
    )
    stages = [
        methods,
        rules,
        [c for c in CLASSES if any(key[2] == c for key in counts)],
    ]
    totals = [Counter(), Counter(), Counter()]
    for (m, r, c), n in counts.items():
        for i, key in enumerate([m, (r, c), c]):
            totals[i][key] += n
    gap = min(0.043, 0.25 / max(1, max(map(len, stages)) - 1))
    scale = (0.90 - gap * (max(map(len, stages)) - 1)) / total
    positions = []
    for stage, tally in zip(stages, totals):
        height = total * scale + gap * (len(stage) - 1)
        y = (1 - height) / 2
        nodes = {}
        for key in stage:
            nodes[key] = (y, y + tally[key] * scale)
            y += tally[key] * scale + gap
        positions.append(nodes)
    fig, ax = plt.subplots(figsize=(18, 10))
    fig.subplots_adjust(left=0.16, right=0.83, bottom=0.12, top=0.85)
    ax.set(xlim=(-0.03, 1.03), ylim=(1, 0))
    ax.axis("off")
    xs = [0, 0.5, 1]
    node_width = 0.018
    source_offset = defaultdict(float)
    target_offset = defaultdict(float)

    def ribbon(stage, source, target, count, color):
        x0, x1 = xs[stage] + node_width / 2, xs[stage + 1] - node_width / 2
        y0 = positions[stage][source][0] + source_offset[stage, source]
        y1 = positions[stage + 1][target][0] + target_offset[stage + 1, target]
        h = count * scale
        mid = (x0 + x1) / 2
        verts = [
            (x0, y0),
            (mid, y0),
            (mid, y1),
            (x1, y1),
            (x1, y1 + h),
            (mid, y1 + h),
            (mid, y0 + h),
            (x0, y0 + h),
            (x0, y0),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.LINETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        ax.add_patch(
            PathPatch(
                MplPath(verts, codes), facecolor=color, alpha=0.35, edgecolor="none"
            )
        )
        source_offset[stage, source] += h
        target_offset[stage + 1, target] += h

    for m in methods:
        for r, c in rules:
            n = counts[m, r, c]
            if n:
                ribbon(0, m, (r, c), n, COLORS[c])
    for r, c in rules:
        ribbon(1, (r, c), c, totals[1][r, c], COLORS[c])
    # Keep rare-rule labels readable without changing count-proportional ribbons.
    rule_label_y = [sum(positions[1][key]) / 2 for key in stages[1]]
    minimum_gap = min(0.075, 0.9 / max(1, len(rule_label_y) - 1))
    for i in range(1, len(rule_label_y)):
        rule_label_y[i] = max(rule_label_y[i], rule_label_y[i - 1] + minimum_gap)
    if rule_label_y:
        rule_label_y[-1] = min(0.97, rule_label_y[-1])
        for i in range(len(rule_label_y) - 2, -1, -1):
            rule_label_y[i] = min(rule_label_y[i], rule_label_y[i + 1] - minimum_gap)
    for stage, nodes in enumerate(positions):
        for key, (lo, hi) in nodes.items():
            color = "#455D73" if stage == 0 else COLORS[key[1] if stage == 1 else key]
            ax.add_patch(
                Rectangle(
                    (xs[stage] - node_width / 2, lo),
                    node_width,
                    hi - lo,
                    facecolor=color,
                    edgecolor="white",
                    lw=0.5,
                )
            )
            if stage == 0:
                label = f"{key}\nn={totals[stage][key]:,}"
                x, align = -0.025, "right"
            elif stage == 1:
                label = (
                    textwrap.fill(RULES.get(key[0], key[0]), width=37)
                    + f"\nn={totals[stage][key]:,}"
                )
                x, align = 0.525, "left"
            else:
                label = (
                    f"{key.upper()}\n{totals[stage][key]:,}"
                    f" ({100*totals[stage][key]/total:.1f}%)"
                )
                x, align = 1.025, "left"
            label_y = (
                rule_label_y[stages[1].index(key)] if stage == 1 else (lo + hi) / 2
            )
            if stage == 1 and abs(label_y - (lo + hi) / 2) > 0.005:
                ax.plot(
                    [0.509, 0.522], [(lo + hi) / 2, label_y], color="#52677B", lw=0.6
                )
            ax.text(
                x,
                label_y,
                label,
                ha=align,
                va="center",
                fontsize=9 if stage == 1 else 12,
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.85,
                    "pad": 2,
                },
            )
    title = "Confidence assignment" + (
        " — " + species.replace("_", " ").capitalize() if species else ""
    )
    fig.suptitle(title, fontsize=22, y=0.96)
    for x, label in zip(
        xs, ["Prediction method", "Matched scoring rule", "Confidence"]
    ):
        ax.text(x, -0.055, label, ha="center", va="bottom", fontsize=14, weight="bold")
    fig.text(
        0.5,
        0.045,
        "Flow width = number of gene-model assessments; "
        "predictions are not matched between methods.\n"
        "High/medium routes require complete CDS and a homology hit. "
        "Intron checks = canonical splicing and regular intron size.",
        ha="center",
        fontsize=11,
        color="#52677B",
    )
    for suffix in ["png", "svg"]:
        fig.savefig(outdir / f"confidence_sankey_matched_rules.{suffix}", dpi=180)
    plt.close(fig)


# The orchestration/layout keeps related report operations together.
# pylint: disable-next=too-many-arguments,too-many-branches,too-many-locals,too-many-statements
def draw_workflow_sankey(
    results,
    outdir,
    species=None,
    source_field="model",
    output_stem="confidence_sankey",
    title=None,
):
    """Render sequential workflow decisions using precomputed paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch, Rectangle, Patch
    from ensembl.genes.annotation_qc.reports.high_confidence_summary import write_tsv

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = Counter()
    for result in results:
        paths.update(result["sankey_paths"])

    def exported(node):
        return (
            node.replace("method:", "species:", 1)
            if source_field == "species" and node.startswith("method:")
            else node
        )

    write_tsv(
        outdir / f"{output_stem}_paths.tsv",
        [
            {
                source_field: m,
                "path": " -> ".join(exported(node) for node in path),
                "confidence": c,
                "count": n,
            }
            for (m, path, c), n in paths.items()
        ],
        [source_field, "path", "confidence", "count"],
    )
    edges = Counter()
    nodes = Counter()
    for (m, path, c), n in paths.items():
        for node in path:
            nodes[node] += n
        for source, target in zip(path[:-1], path[1:]):
            edges[m, source, target] += n
    write_tsv(
        outdir / f"{output_stem}_flows.tsv",
        [
            {source_field: m, "source": exported(a), "target": exported(b), "count": n}
            for (m, a, b), n in edges.items()
        ],
        [source_field, "source", "target", "count"],
    )
    total = sum(n for (_, _, _), n in paths.items())
    if not total:
        return
    methods = [r["model"] for r in results]
    gap = 0.033
    scale = (0.9 - gap * max(4, len(methods) - 1)) / total
    x = {"method": 0, "required": 0.24, "high": 0.48, "medium": 0.74, "final": 1}
    position = {}
    cursor = 0.05
    for m in methods:
        node = "method:" + m
        position[node] = (0, cursor, cursor + nodes[node] * scale)
        cursor += nodes[node] * scale + gap
    high_height = nodes["final:high"] * scale
    medium_height = nodes["final:medium"] * scale
    # Empty regions let terminal ribbons bypass later tests without phantom nodes.
    starts = {
        "required:yes": 0.05,
        "required:no": 0.95 - nodes["required:no"] * scale,
        "high:yes": 0.05,
        "high:no": 0.05 + high_height + gap,
        "medium:yes": 0.05 + high_height + gap,
        "medium:no": 0.05 + high_height + gap + medium_height + gap,
        "final:high": 0.05,
        "final:medium": 0.05 + high_height + gap,
        "final:low": 0.95 - nodes["final:low"] * scale,
    }
    for node, start in starts.items():
        if nodes[node]:
            position[node] = (
                x[node.split(":", maxsplit=1)[0]],
                start,
                start + nodes[node] * scale,
            )
    fig, ax = plt.subplots(figsize=(24, 13.5))
    fig.subplots_adjust(
        left=0.17 if source_field == "species" else 0.10,
        right=0.87,
        bottom=0.16 if source_field == "species" else 0.14,
        top=0.79,
    )
    ax.set(xlim=(-0.02, 1.02), ylim=(1, 0))
    ax.axis("off")
    palette = [
        "#0072B2",
        "#E69F00",
        "#009E73",
        "#CC79A7",
        "#D55E00",
        "#56B4E9",
        "#6F5E9C",
    ]
    method_colors = {m: palette[i % len(palette)] for i, m in enumerate(methods)}
    source_offset = defaultdict(float)
    target_offset = defaultdict(float)
    width = 0.011
    for m in methods:
        model_edges = [(a, b, n) for (model, a, b), n in edges.items() if model == m]
        model_edges.sort(
            key=lambda edge: (
                position[edge[0]][0],
                position[edge[0]][1],
                position[edge[1]][1],
            )
        )
        for source, target, n in model_edges:
            sx, sy, _ = position[source]
            tx, ty, _ = position[target]
            y0 = sy + source_offset[source]
            y1 = ty + target_offset[target]
            h = n * scale
            x0, x1 = sx + width / 2, tx - width / 2
            mid = (x0 + x1) / 2
            points = [
                (x0, y0),
                (mid, y0),
                (mid, y1),
                (x1, y1),
                (x1, y1 + h),
                (mid, y1 + h),
                (mid, y0 + h),
                (x0, y0 + h),
                (x0, y0),
            ]
            codes = [
                MplPath.MOVETO,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.LINETO,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CLOSEPOLY,
            ]
            ax.add_patch(
                PathPatch(
                    MplPath(points, codes),
                    facecolor=method_colors[m],
                    alpha=0.40,
                    edgecolor="none",
                )
            )
            source_offset[source] += h
            target_offset[target] += h
    labels = {
        "required:yes": "Yes → test 2",
        "required:no": "No → LOW",
        "high:yes": "Yes → HIGH",
        "high:no": "No → test 3",
        "medium:yes": "Yes → MEDIUM",
        "medium:no": "No → LOW",
    }
    for node, (nx, lo, hi) in position.items():
        if node.startswith("method:"):
            color = method_colors[node.split(":", 1)[1]]
        elif node.startswith("final:"):
            color = COLORS[node.split(":")[1]]
        elif node in {"required:no", "medium:no"}:
            color = COLORS["low"]
        elif node == "high:yes":
            color = COLORS["high"]
        elif node == "medium:yes":
            color = COLORS["medium"]
        else:
            color = "#526A7C"
        ax.add_patch(
            Rectangle(
                (nx - width / 2, lo),
                width,
                hi - lo,
                facecolor=color,
                edgecolor="white",
                lw=0.6,
            )
        )
        if node.startswith("method:"):
            label = node.split(":", 1)[1] + f"\nn={nodes[node]:,}"
            tx, ha = nx - 0.012, "right"
        elif node.startswith("final:"):
            label = node.split(":")[1].upper() + f"\nn={nodes[node]:,}"
            tx, ha = nx + 0.012, "left"
        else:
            label = labels[node] + f"\nn={nodes[node]:,}"
            tx, ha = nx + 0.012, "left"
        ax.text(
            tx,
            (lo + hi) / 2,
            label,
            ha=ha,
            va="center",
            fontsize=11,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 2},
        )
    no_utr = all(r["overview"]["no_utr"] is True for r in results)
    titles = [
        (
            "Species\n(all methods combined)"
            if source_field == "species"
            else "Prediction\nmethod"
        ),
        "1  Required evidence\nComplete CDS +\nhomology hit",
        (
            "2  High-confidence test\nBoth intron checks +\nRNA support + error ≤10"
            if no_utr
            else "2  High-confidence test\nAny high-confidence rule\nunder the recorded UTR policy"
        ),
        (
            "3  Medium-confidence test\n≥1 intron check + either:\nRNA support &"
            " error ≤10\nOR error ≤5 & score ≥80"
        ),
        "Final\nconfidence",
    ]
    for nx, label in zip([0, 0.24, 0.48, 0.74, 1], titles):
        ax.text(
            nx,
            -0.045,
            label,
            ha="center",
            va="bottom",
            fontsize=13,
            weight="bold",
            linespacing=1.5,
        )
    title = title or "Confidence assignment" + (
        " — " + species.replace("_", " ").capitalize() if species else ""
    )
    fig.suptitle(title, fontsize=25, weight="bold", y=0.97)
    fig.legend(
        handles=[Patch(color=method_colors[m], label=m) for m in methods],
        loc="center",
        bbox_to_anchor=(0.5, 0.10 if source_field == "species" else 0.09),
        ncol=min(3 if source_field == "species" else 7, len(methods)),
        frameon=False,
        fontsize=12,
    )
    footer = (
        "Ribbon width = summed gene-model assessments across methods; "
        "colour = species. These are not unique genes.\n"
        "Later tests are skipped after assignment. "
        "Intron checks: canonical splicing and regular intron size."
        if source_field == "species"
        else "Ribbon width = gene-model assessments; colour = prediction method. "
        "A ribbon skips later tests once HIGH or LOW is assigned.\n"
        "Intron checks: canonical splicing and regular intron size. "
        "Predictions are not matched between methods."
    )
    fig.text(0.5, 0.035, footer, ha="center", fontsize=11, color="#52677B")
    for suffix in ["png", "svg"]:
        fig.savefig(outdir / f"{output_stem}.{suffix}", dpi=180)
    plt.close(fig)


def draw_sankey(results, outdir, species=None):
    """Keep the matched-rule view and make the sequential workflow the main Sankey."""
    draw_matched_rule_sankey(results, outdir, species=species)
    draw_workflow_sankey(results, outdir, species=species)

"""
One-to-one coordinate-exact CDS matching between reference and query genes.

Definitions:

    population   R = every reference gene and Q = every query gene in the
                 GeneModels passed in, i.e. after the evaluation mode, biotype
                 filters, transcript selection and region/sequence scope that
                 the runner applied. These are the comparator's
                 total_reference_genes and total_consensus_genes. Genes without
                 a comparable CDS (no transcript with exons and CDS) stay in R
                 and Q; they cannot be matched and count as unmatched.
    signature    a transcript's CDS signature is (sequence, gene strand, sorted
                 CDS intervals) with the coordinates exactly as parsed
                 (zero-based half-open, stop codon included or excluded as in the
                 file). Transcripts without CDS have no signature. The
                 transcripts used are those build_gene_models assigns to a gene
                 (at least one exon; CDS-only transcripts get exons at parse
                 time), the same ones the per-gene cds_coordinate_exact column
                 compares.
    edge         reference gene r and query gene q are joined when some
                 transcript of r and some transcript of q have the same
                 signature. Edges are not limited to the best pairs chosen by
                 classify_loci.
    matching     a maximum-cardinality bipartite matching of that graph: each
                 reference and each query gene is used at most once. Kuhn's
                 augmenting-path algorithm (breadth-first search per reference
                 gene) is exact for maximum cardinality.
    metrics      TP = matched pairs; precision = TP / Q; recall = TP / R;
                 F1 = 2 TP / (R + Q) (the harmonic mean of the two). A metric
                 whose denominator is 0 is None (undefined), not 0.

Tie resolution is deterministic: reference genes are processed in order of
(sequence, start, end, gene_id) and each gene's candidates are tried in the same
order, so the first augmenting path found in that order is used. Other maximum
matchings with the same TP may exist; they differ only inside connected
components with more than one edge (reported per pair as component sizes).
"""

from collections import defaultdict, deque

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise.classify import GeneModels

PAIR_COLUMNS = [
    "reference_gene_id",
    "reference_original_gene_id",
    "query_gene_id",
    "query_original_gene_id",
    "chrom",
    "strand",
    "reference_start",
    "reference_end",
    "query_start",
    "query_end",
    "cds_start",
    "cds_end",
    "cds_segments",
    "shared_cds_signatures",
    "supporting_reference_transcript_ids",
    "supporting_query_transcript_ids",
    "reference_exact_candidates",
    "query_exact_candidates",
    "component_reference_genes",
    "component_query_genes",
]

DEFINITION = {
    "edge": (
        "a reference and a query gene share at least one transcript pair with "
        "identical CDS segment coordinates on the same sequence and strand "
        "(coordinates as parsed; empty CDS excluded)"
    ),
    "matching": (
        "maximum-cardinality bipartite matching (Kuhn's augmenting paths); each "
        "reference and each query gene is used at most once"
    ),
    "tie_resolution": (
        "deterministic: reference genes in (sequence, start, end, gene_id) order, "
        "candidates in the same order; other maximum matchings with the same TP "
        "can exist within components that have more than one edge"
    ),
    "precision": "TP / total_consensus_genes (Q)",
    "recall": "TP / total_reference_genes (R)",
    "f1": "2 TP / (R + Q)",
    "undefined": "a value is null when its denominator is 0",
    "population": (
        "R and Q are all genes after filtering, transcript selection and scope; "
        "genes without a comparable CDS cannot match and count as unmatched"
    ),
}


def _signatures(models: GeneModels) -> dict:
    """signature -> [(gene row, transcript id)] for every transcript with CDS."""
    index = defaultdict(list)
    genes = models.genes
    for row, (chrom, strand, tids) in enumerate(
        zip(genes["Chromosome"], genes["Strand"], genes["transcript_ids"])
    ):
        for tid in tids:
            cds = models.cds.get(tid)
            if cds:
                index[(chrom, strand, cds)].append((row, tid))
    return index


def exact_cds_edges(reference: GeneModels, query: GeneModels) -> dict:
    """
    Gene pairs with at least one identical CDS signature.
    Args:
            reference, query: GeneModels after filtering, selection and scope
    Returns:
            {(reference row, query row): {"signatures": set of signatures,
             "reference_tids": set, "query_tids": set}}
    """
    query_index = _signatures(query)
    edges: dict = {}
    for signature, ref_members in _signatures(reference).items():
        query_members = query_index.get(signature)
        if not query_members:
            continue
        for ref_row, ref_tid in ref_members:
            for query_row, query_tid in query_members:
                edge = edges.setdefault(
                    (ref_row, query_row),
                    {"signatures": set(), "reference_tids": set(), "query_tids": set()},
                )
                edge["signatures"].add(signature)
                edge["reference_tids"].add(ref_tid)
                edge["query_tids"].add(query_tid)
    return edges


def maximum_matching(adjacency: dict, order: list) -> dict:
    """
    Maximum-cardinality bipartite matching (Kuhn's algorithm, BFS augmentation).
    Args:
            adjacency: left vertex -> right vertices, in the order to try them
            order: left vertices in processing order
    Returns:
            {left vertex: right vertex}
    """
    match_left: dict = {}
    match_right: dict = {}
    for root in order:
        parent: dict = {}
        queue = deque([root])
        free = None
        while queue and free is None:
            left = queue.popleft()
            for right in adjacency.get(left, ()):
                if right in parent:
                    continue
                parent[right] = left
                if right not in match_right:
                    free = right
                    break
                queue.append(match_right[right])
        # Flip the alternating path root -> ... -> free.
        right = free
        while right is not None:
            left = parent[right]
            previous = match_left.get(left)
            match_left[left], match_right[right] = right, left
            right = previous
    return match_left


def _components(edges: dict) -> tuple[dict, int, int]:
    """
    Connected components of the exact-CDS gene graph.
    Returns:
            ({edge: (reference genes, query genes) in its component},
             number of components, number of components with more than one edge)
    """
    parent: dict = {}

    def find(node):
        while parent.setdefault(node, node) != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for ref_row, query_row in edges:
        parent[find(("r", ref_row))] = find(("q", query_row))
    genes: dict = defaultdict(lambda: [0, 0])
    for node in list(parent):
        genes[find(node)][0 if node[0] == "r" else 1] += 1
    edge_counts: dict = defaultdict(int)
    for ref_row, _ in edges:
        edge_counts[find(("r", ref_row))] += 1
    sizes = {edge: tuple(genes[find(("r", edge[0]))]) for edge in edges}
    several = sum(1 for count in edge_counts.values() if count > 1)
    return sizes, len(edge_counts), several


def _gene_columns(genes: pd.DataFrame) -> dict:
    return {
        column: genes[column].tolist()
        for column in ("Chromosome", "Start", "End", "gene_id", "original_gene_id")
    }


def _ratio(count: int, total: int) -> float | None:
    return round(count / total, 4) if total else None


def _ordered(cols: dict, rows) -> list:
    """Gene rows sorted by (sequence, start, end, gene_id)."""
    return sorted(
        rows,
        key=lambda row: (
            cols["Chromosome"][row],
            cols["Start"][row],
            cols["End"][row],
            cols["gene_id"][row],
        ),
    )


def _adjacency(edges: dict, query_cols: dict) -> dict:
    """Reference row -> its exact query rows, in (sequence, start, end, gene_id) order."""
    adjacency: dict = defaultdict(list)
    for ref_row, query_row in edges:
        adjacency[ref_row].append(query_row)
    return {row: _ordered(query_cols, rows) for row, rows in adjacency.items()}


def _genes_with_cds(models: GeneModels) -> int:
    return int(
        sum(
            any(models.cds.get(tid) for tid in tids)
            for tids in models.genes["transcript_ids"]
        )
    )


def _pair_table(
    matching: dict, edges: dict, ref_cols: dict, query_cols: dict
) -> pd.DataFrame:
    """One row per matched pair (PAIR_COLUMNS), in reference gene order."""
    ref_degree: dict = defaultdict(int)
    query_degree: dict = defaultdict(int)
    for ref_row, query_row in edges:
        ref_degree[ref_row] += 1
        query_degree[query_row] += 1
    components = _components(edges)[0]
    records = []
    for ref_row in _ordered(ref_cols, matching):
        query_row = matching[ref_row]
        edge = edges[(ref_row, query_row)]
        chrom, strand, cds = min(edge["signatures"])
        records.append(
            (
                ref_cols["gene_id"][ref_row],
                ref_cols["original_gene_id"][ref_row],
                query_cols["gene_id"][query_row],
                query_cols["original_gene_id"][query_row],
                chrom,
                strand,
                ref_cols["Start"][ref_row],
                ref_cols["End"][ref_row],
                query_cols["Start"][query_row],
                query_cols["End"][query_row],
                cds[0][0],
                cds[-1][1],
                len(cds),
                len(edge["signatures"]),
                ",".join(sorted(edge["reference_tids"])),
                ",".join(sorted(edge["query_tids"])),
                ref_degree[ref_row],
                query_degree[query_row],
                *components[(ref_row, query_row)],
            )
        )
    return pd.DataFrame(records, columns=PAIR_COLUMNS)


def one_to_one_exact_cds(
    reference: GeneModels, query: GeneModels
) -> tuple[dict, pd.DataFrame]:
    """
    One-to-one coordinate-exact CDS precision, recall and F1.
    Args:
            reference, query: GeneModels after filtering, selection and scope
    Returns:
            (summary dict, DataFrame of matched pairs with PAIR_COLUMNS; gene and
            CDS coordinates zero-based half-open)
    """
    total_ref, total_query = len(reference.genes), len(query.genes)
    edges = exact_cds_edges(reference, query)
    # Gene rows are positions in the gene tables; read columns once.
    ref_cols = _gene_columns(reference.genes)
    query_cols = _gene_columns(query.genes)

    adjacency = _adjacency(edges, query_cols)
    matching = maximum_matching(adjacency, _ordered(ref_cols, adjacency))

    _, n_components, several_edges = _components(edges)
    with_candidate = len({query_row for _, query_row in edges})
    tp = len(matching)
    summary = {
        "reference_genes": total_ref,
        "query_genes": total_query,
        "reference_genes_with_cds": _genes_with_cds(reference),
        "query_genes_with_cds": _genes_with_cds(query),
        "exact_gene_pairs": len(edges),
        "reference_genes_with_exact_candidate": len(adjacency),
        "query_genes_with_exact_candidate": with_candidate,
        "true_positives": tp,
        "unmatched_reference_genes": total_ref - tp,
        "unmatched_query_genes": total_query - tp,
        "reference_genes_with_exact_candidate_unmatched": len(adjacency) - tp,
        "query_genes_with_exact_candidate_unmatched": with_candidate - tp,
        "exact_components": n_components,
        "exact_components_with_several_edges": several_edges,
        "precision": _ratio(tp, total_query),
        "recall": _ratio(tp, total_ref),
        "f1": _ratio(2 * tp, total_ref + total_query),
        "definition": DEFINITION,
    }
    return summary, _pair_table(matching, edges, ref_cols, query_cols)

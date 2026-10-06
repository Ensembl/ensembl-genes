"""
Independent CDS-signature analysis, generalised from the 2026-10-05 benchmark
(scripts/independent_metrics.py) to experiment configurations.

For one reference and its runs, per policy:

    reference eligibility   gene biotype filter (evaluation_mode protein_coding and
                            reference_gene_biotypes) and transcript biotype filter,
                            restricted to the scope sequences
    selection               all | longest_cds | canonical (Ensembl_canonical tag,
                            else longest CDS); ties -> first transcript in file
    per reference gene      cds_exact_any (some selected isoform's CDS chain equals a
                            query CDS chain on the same strand), multi_segment,
                            cds_intron_chain_any_pair (identical CDS intron chain with
                            ANY query transcript, not only the comparator's best pair)

The tables are written to ``<output_dir>/independent/<run_id>/<policy>.*`` and the
dashboard uses them only after they agree with the comparator per gene.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from ensembl.genes.annotation_qc.benchmark.independent import reader as ig


def cds_len(ivs):
    return sum(e - s + 1 for s, e, *_ in ivs)


def pick(tkeys, ann, mode):
    if mode == "canonical":
        for k in tkeys:
            if "Ensembl_canonical" in ann.tx[k]["tags"]:
                return [k]
    best = None
    for k in tkeys:
        if best is None or cds_len(ann.cds[k]) > cds_len(ann.cds[best]):
            best = k
    return [best]


def models(ann, scope, gene_ok, tx_ok, mode):
    by_gene = defaultdict(list)
    for tkey, t in ann.tx.items():
        if (
            t is None
            or (scope is not None and t["chrom"] not in scope)
            or not ann.cds.get(tkey)
        ):
            continue
        gene = ann.genes.get(t["gene"])
        if gene is None or not gene_ok(gene) or not tx_ok(t):
            continue
        by_gene[t["gene"]].append(tkey)
    out = {}
    for g, tkeys in by_gene.items():
        sel = tkeys if mode == "all" else pick(tkeys, ann, mode)
        out[g] = [
            (k, ann.tx[k]["chrom"], ann.tx[k]["strand"], ig.cds_chain(ann.cds[k]))
            for k in sel
        ]
    return out


def introns(chain):
    return tuple((chain[i][1] + 1, chain[i + 1][0] - 1) for i in range(len(chain) - 1))


def evaluate(ref_models, q_models):
    q_sig, q_ich = defaultdict(set), defaultdict(set)
    for gk, txs in q_models.items():
        g = gk[1]
        for _, c, s, ch in txs:
            q_sig[(c, s, ch)].add(g)
            if len(ch) > 1:
                q_ich[(c, s, introns(ch))].add(g)
    ref_sig = {(c, s, ch) for txs in ref_models.values() for _, c, s, ch in txs}
    ref_rows = []
    for g, txs in ref_models.items():
        multi = max(len(ch) for *_, ch in txs) > 1
        exact_q, chain_q = set(), set()
        for _, c, s, ch in txs:
            exact_q |= q_sig.get((c, s, ch), set())
            if len(ch) > 1:
                chain_q |= q_ich.get((c, s, introns(ch)), set())
        ref_rows.append(
            {
                "gene_id": g[1],
                "chrom": g[0],
                "multi_segment": multi,
                "cds_exact_any": bool(exact_q),
                "cds_exact_query_genes": ",".join(sorted(exact_q)),
                "cds_intron_chain_any_pair": (bool(chain_q) if multi else None),
            }
        )
    q_rows = [
        {
            "gene_id": g[1],
            "chrom": g[0],
            "cds_exact_any": any((c, s, ch) in ref_sig for _, c, s, ch in txs),
        }
        for g, txs in q_models.items()
    ]
    return ref_rows, q_rows


def _write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]) if rows else ["gene_id"], delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})


def run_independent(
    reference_annotation: Path,
    scope: set | None,
    evaluation: dict,
    runs: dict[str, Path],
    policies: dict[str, tuple[str, str]],
    out_root: Path,
    log=print,
) -> dict:
    """
    Args:
            reference_annotation: prepared reference annotation
            scope: sequence names in scope (None = all)
            evaluation: experiment evaluation block
            runs: run_id -> prepared query annotation
            policies: policy_id -> (reference selection, query selection)
            out_root: <output_dir>/independent
    """
    gene_biotypes = set(evaluation.get("reference_gene_biotypes") or [])
    if evaluation.get("evaluation_mode") == "protein_coding":
        gene_biotypes |= {"protein_coding"}
    tx_biotypes = set(evaluation.get("reference_transcript_biotypes") or [])
    gene_ok = (
        (lambda g: g["biotype"] in gene_biotypes) if gene_biotypes else (lambda g: True)
    )
    tx_ok = (lambda t: t["biotype"] in tx_biotypes) if tx_biotypes else (lambda t: True)
    keep = (
        (lambda a, f: (a.get("biotype") or a.get("gene_biotype", "")) in gene_biotypes)
        if gene_biotypes
        else None
    )
    ref = ig.read(str(reference_annotation), keep_gene=keep)
    summary = {}
    for run_id, path in runs.items():
        q = ig.read(str(path))
        for policy_id, (ref_sel, query_sel) in policies.items():
            rm = models(ref, scope, gene_ok, tx_ok, ref_sel)
            qm = models(q, scope, lambda g: True, lambda t: True, query_sel)
            rr, qr = evaluate(rm, qm)
            _write(out_root / run_id / f"{policy_id}.reference_genes.tsv", rr)
            _write(out_root / run_id / f"{policy_id}.query_genes.tsv", qr)
            summary[f"{run_id}/{policy_id}"] = {
                "reference_genes": len(rr),
                "query_genes": len(qr),
                "ref_cds_exact": sum(r["cds_exact_any"] for r in rr),
            }
            log(
                f"  independent {run_id}/{policy_id}: {summary[f'{run_id}/{policy_id}']}"
            )
    (out_root / "summary.json").parent.mkdir(parents=True, exist_ok=True)
    (out_root / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary

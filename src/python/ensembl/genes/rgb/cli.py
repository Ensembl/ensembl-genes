from __future__ import annotations

import argparse
import os
import sys
from typing import Optional

import pandas as pd

from .db import (
    DBParams,
    connect,
    extract_all_genes,
    extract_all_transcripts,
    extract_all_translations,
    extract_all_exon_chains,
    list_seq_regions,
)
from .io import append_tsv, write_df, write_manifest
from .loci import build_loci
from .utils import ensure_dir, make_run_id
from .summary import summarize_loci, load_mapping
from .layer_map import fetch_layer_static, parse_config_block, to_yaml
from .report import compute_retention_by_biotype
from .compete import analyze_competition


def _add_common_db_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--host", required=True)
    p.add_argument("--port", type=int, default=3306)
    p.add_argument("--user", required=True)
    p.add_argument("--password", default="")
    p.add_argument("--core_db", required=True)
    p.add_argument("--layer_db", required=True)
    p.add_argument("--coord_system_name", default=None, help="coord_system.name (chromosome|scaffold), 'toplevel', or 'auto' (default)")


def _add_common_out_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--output_dir", required=True, help="Absolute output directory")
    p.add_argument("--format", choices=["tsv", "csv", "parquet"], default="tsv")


def cmd_extract(args: argparse.Namespace) -> int:
    run_id = args.run_id or make_run_id(args.core_db, args.layer_db)
    out_root = os.path.join(args.output_dir, run_id)
    out_dir = os.path.join(out_root, "extract")
    ensure_dir(out_dir)

    core_params = DBParams(args.host, args.port, args.user, args.password, args.core_db)
    layer_params = DBParams(args.host, args.port, args.user, args.password, args.layer_db)

    # Determine seq_regions to process
    # Resolve coord_system selection
    requested = (args.coord_system_name or "auto").lower()
    resolved = None
    resolved_layer = None
    with connect(core_params) as core_conn, connect(layer_params) as layer_conn:
        from .db import detect_coord_system, has_toplevel
        if requested in (None, "", "auto"):
            cs = detect_coord_system(core_conn)
            if cs:
                resolved = cs
            elif has_toplevel(core_conn):
                resolved = "toplevel"
            else:
                resolved = None
        elif requested == "toplevel":
            resolved = "toplevel"
        else:
            resolved = requested

        # List regions under the resolved mode
        core_regions = set(list_seq_regions(core_conn, resolved))
        layer_regions = set(list_seq_regions(layer_conn, resolved))
        seq_regions = sorted(core_regions.union(layer_regions))

    manifest = {
        "phase": "extract",
        "run_id": run_id,
        "core_db": args.core_db,
        "layer_db": args.layer_db,
        "coord_system_requested": args.coord_system_name,
        "coord_system_resolved": resolved,
        "seq_regions": seq_regions,
        "format": args.format,
    }
    write_manifest(manifest, os.path.join(out_root, "manifest.extract.json"))

    # Extract genes in one go per DB (region‑chunked within helper)
    with connect(core_params) as core_conn:
        core_genes = extract_all_genes(core_conn, seq_regions, resolved)
        core_tx = extract_all_transcripts(core_conn, seq_regions, resolved)
        core_tr = extract_all_translations(core_conn)
        core_ex = extract_all_exon_chains(core_conn, seq_regions, resolved)
    with connect(layer_params) as layer_conn:
        layer_genes = extract_all_genes(layer_conn, seq_regions, resolved)
        layer_tx = extract_all_transcripts(layer_conn, seq_regions, resolved)
        layer_tr = extract_all_translations(layer_conn)
        layer_ex = extract_all_exon_chains(layer_conn, seq_regions, resolved)

    # Write outputs
    core_path = os.path.join(out_dir, f"core_genes.{args.format}")
    layer_path = os.path.join(out_dir, f"layer_genes.{args.format}")
    core_tx_path = os.path.join(out_dir, f"core_transcripts.{args.format}")
    layer_tx_path = os.path.join(out_dir, f"layer_transcripts.{args.format}")
    core_tr_path = os.path.join(out_dir, f"core_translations.{args.format}")
    layer_tr_path = os.path.join(out_dir, f"layer_translations.{args.format}")
    core_ex_path = os.path.join(out_dir, f"core_exons.{args.format}")
    layer_ex_path = os.path.join(out_dir, f"layer_exons.{args.format}")

    write_df(core_genes, core_path, fmt=args.format)
    write_df(layer_genes, layer_path, fmt=args.format)
    write_df(core_tx, core_tx_path, fmt=args.format)
    write_df(layer_tx, layer_tx_path, fmt=args.format)
    write_df(core_tr, core_tr_path, fmt=args.format)
    write_df(layer_tr, layer_tr_path, fmt=args.format)
    write_df(core_ex, core_ex_path, fmt=args.format)
    write_df(layer_ex, layer_ex_path, fmt=args.format)

    print(f"[extract] wrote {len(core_genes)} core genes → {core_path}")
    print(f"[extract] wrote {len(layer_genes)} layer genes → {layer_path}")
    print(f"[extract] wrote {len(core_tx)} core transcripts → {core_tx_path}")
    print(f"[extract] wrote {len(layer_tx)} layer transcripts → {layer_tx_path}")
    print(f"[extract] wrote {len(core_tr)} core translations → {core_tr_path}")
    print(f"[extract] wrote {len(layer_tr)} layer translations → {layer_tr_path}")
    print(f"[extract] wrote {len(core_ex)} core exons → {core_ex_path}")
    print(f"[extract] wrote {len(layer_ex)} layer exons → {layer_ex_path}")
    print(f"[extract] run_id={run_id}")
    return 0


def cmd_loci(args: argparse.Namespace) -> int:
    run_id = args.run_id
    if not run_id:
        print("--run_id is required (output from extract)", file=sys.stderr)
        return 2
    out_root = os.path.join(args.output_dir, run_id)
    extract_dir = os.path.join(out_root, "extract")
    loci_dir = os.path.join(out_root, "loci")
    ensure_dir(loci_dir)

    # Load gene tables
    def _load(path_base: str) -> pd.DataFrame:
        for fmt in (args.format, "tsv", "csv", "parquet"):
            path = os.path.join(extract_dir, f"{path_base}.{fmt}")
            if os.path.exists(path):
                if fmt == "parquet":
                    return pd.read_parquet(path)
                sep = "\t" if fmt == "tsv" else ","
                return pd.read_csv(path, sep=sep)
        raise FileNotFoundError(f"Could not find {path_base} with any supported extension in {extract_dir}")

    core_genes = _load("core_genes")
    layer_genes = _load("layer_genes")

    strict_df, expanded_df, map_df = build_loci(core_genes, layer_genes, args.locus_gap_bp)

    write_df(strict_df, os.path.join(loci_dir, f"loci.strict.{args.format}"), fmt=args.format)
    write_df(expanded_df, os.path.join(loci_dir, f"loci.expanded.{args.format}"), fmt=args.format)
    write_df(map_df, os.path.join(loci_dir, f"gene_to_locus.{args.format}"), fmt=args.format)

    manifest = {
        "phase": "loci",
        "run_id": run_id,
        "locus_gap_bp": args.locus_gap_bp,
        "inputs": {
            "core_genes": "core_genes",
            "layer_genes": "layer_genes",
        },
    }
    write_manifest(manifest, os.path.join(out_root, "manifest.loci.json"))

    print(
        f"[loci] strict:{len(strict_df)} expanded:{len(expanded_df)} mapped genes:{len(map_df)} → {loci_dir}"
    )
    return 0


def _load_table(extract_dir: str, base: str, prefer_fmt: str) -> pd.DataFrame:
    for fmt in (prefer_fmt, "tsv", "csv", "parquet"):
        path = os.path.join(extract_dir, f"{base}.{fmt}")
        if os.path.exists(path):
            if fmt == "parquet":
                return pd.read_parquet(path)
            sep = "\t" if fmt == "tsv" else ","
            return pd.read_csv(path, sep=sep)
    raise FileNotFoundError(f"Could not find {base} with any supported extension in {extract_dir}")


def cmd_summarize(args: argparse.Namespace) -> int:
    out_root = os.path.join(args.output_dir, args.run_id)
    extract_dir = os.path.join(out_root, "extract")
    loci_dir = os.path.join(out_root, "loci")
    summaries_dir = os.path.join(out_root, "summaries")
    ensure_dir(summaries_dir)

    # Load inputs
    core_genes = _load_table(extract_dir, "core_genes", args.format)
    layer_genes = _load_table(extract_dir, "layer_genes", args.format)
    core_tx = _load_table(extract_dir, "core_transcripts", args.format)
    layer_tx = _load_table(extract_dir, "layer_transcripts", args.format)
    loci_df = _load_table(loci_dir, "loci.strict", args.format)
    gene_map = _load_table(loci_dir, "gene_to_locus", args.format)

    evidence_map = load_mapping(args.evidence_class_map) if args.evidence_class_map else {}

    locus_summary = summarize_loci(
        loci_df=loci_df,
        gene_map_df=gene_map,
        core_genes=core_genes,
        layer_genes=layer_genes,
        core_tx=core_tx,
        layer_tx=layer_tx,
        evidence_map=evidence_map,
        locus_gap_bp=args.locus_gap_bp,
    )

    write_df(locus_summary, os.path.join(summaries_dir, f"locus_summary.{args.format}"), fmt=args.format)

    manifest = {
        "phase": "summarize",
        "run_id": args.run_id,
        "evidence_class_map": args.evidence_class_map,
        "locus_gap_bp": args.locus_gap_bp,
        "rows": len(locus_summary),
    }
    write_manifest(manifest, os.path.join(out_root, "manifest.summarize.json"))

    print(f"[summarize] wrote {len(locus_summary)} loci → {summaries_dir}")
    return 0


def cmd_export(args: argparse.Namespace) -> int:
    out_root = os.path.join(args.output_dir, args.run_id)
    summaries_dir = os.path.join(out_root, "summaries")
    review_dir = os.path.join(out_root, "review")
    ensure_dir(review_dir)

    # Load locus summary
    summ_path = os.path.join(summaries_dir, f"locus_summary.{args.format}")
    if not os.path.exists(summ_path):
        raise FileNotFoundError("Run summarize first; missing locus_summary")
    if args.format == "parquet":
        df = pd.read_parquet(summ_path)
    else:
        sep = "\t" if args.format == "tsv" else ","
        df = pd.read_csv(summ_path, sep=sep)

    def annotate_counts(df_in: pd.DataFrame) -> pd.DataFrame:
        # Load core genes to compute same/opp strand overlaps
        extract_dir = os.path.join(out_root, "extract")
        core_genes_path = os.path.join(extract_dir, f"core_genes.{args.format}")
        if args.format == "parquet" and os.path.exists(core_genes_path):
            cg = pd.read_parquet(core_genes_path)
        else:
            sep = "\t" if args.format == "tsv" else ","
            cg = pd.read_csv(core_genes_path, sep=sep)
        rows = []
        for r in df_in.itertuples():
            seq = str(r.seq_region_name)
            s = int(r.locus_start) - args.bed_pad_bp
            e = int(r.locus_end) + args.bed_pad_bp
            strand = int(r.seq_region_strand)
            sub = cg[cg.seq_region_name.astype(str) == seq]
            same = sub[sub.seq_region_strand == strand]
            opp = sub[sub.seq_region_strand == -strand]
            same_n = int(((same.seq_region_start <= e) & (same.seq_region_end >= s)).sum())
            opp_n = int(((opp.seq_region_start <= e) & (opp.seq_region_end >= s)).sum())
            rows.append((same_n, opp_n, s, e))
        out = df_in.copy()
        out["same_strand_core"] = [r[0] for r in rows]
        out["opp_strand_core"] = [r[1] for r in rows]
        out["bed_start0"] = [max(0, r[2]-1) for r in rows]
        out["bed_end0"] = [r[3] for r in rows]
        return out

    def to_bed(df_in: pd.DataFrame, path: str) -> None:
        tmp = df_in.copy()
        tmp = annotate_counts(tmp)
        tmp["chrom"] = tmp.seq_region_name
        tmp["name"] = (
            "LOCUS:" + tmp.seq_region_name.astype(str)
            + ":" + tmp.seq_region_strand.astype(str)
            + ":" + tmp.locus_start.astype(int).astype(str)
            + ":" + tmp.locus_end.astype(int).astype(str)
            + "|score:" + (tmp.get("no_core_score", 0.0).fillna(0.0)).map(lambda x: f"{x:.2f}")
            + "|same:" + tmp.same_strand_core.astype(int).astype(str)
            + "|opp:" + tmp.opp_strand_core.astype(int).astype(str)
        )
        tmp["score"] = (tmp.get("no_core_score", 0.0).fillna(0.0) * 1000).astype(int).clip(lower=0, upper=1000)
        tmp["strand"] = tmp.seq_region_strand.map({1: "+", -1: "-"}).fillna(".")
        cols = ["chrom", "bed_start0", "bed_end0", "name", "score", "strand"]
        tmp[cols].to_csv(path, sep="\t", header=False, index=False)

    # Evidence-rich no-core
    miss = df[(df.evidence_rich_no_core_flag == 1) & (df.no_core_score >= args.no_core_min_score)]
    miss = miss.sort_values(["no_core_score", "layer_span_bp"], ascending=[False, False]).head(args.top_n)
    # Apply strand-based filtering if requested
    if not miss.empty and args.strand_filter != "none":
        miss_ann = annotate_counts(miss)
        if args.strand_filter == "same_empty":
            miss = miss_ann[miss_ann.same_strand_core == 0]
        elif args.strand_filter == "both_empty":
            miss = miss_ann[(miss_ann.same_strand_core == 0) & (miss_ann.opp_strand_core == 0)]
    to_bed(miss, os.path.join(review_dir, "missing_gene_highconf.bed"))

    # Underbuilt loci
    under = df[df.underbuilt_locus_flag == 1]
    under = under.sort_values(["coverage_fraction", "layer_to_core_tx_ratio"], ascending=[True, False]).head(args.top_n)
    to_bed(under, os.path.join(review_dir, "underbuilt_loci.bed"))

    print(
        f"[export] wrote {len(miss)} missing_gene_highconf and {len(under)} underbuilt loci → {review_dir}"
    )
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    p = argparse.ArgumentParser(prog="rgb", description="Reverse Gene Builder CLI (Phase 1–2)")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_extract = sub.add_parser("extract", help="Extract gene tables from core and layer DBs")
    _add_common_db_args(p_extract)
    _add_common_out_args(p_extract)
    p_extract.add_argument("--run_id", default=None, help="Override auto run id")
    p_extract.set_defaults(func=cmd_extract)

    p_loci = sub.add_parser("loci", help="Build strict/expanded loci from extracted tables")
    _add_common_out_args(p_loci)
    p_loci.add_argument("--run_id", required=True, help="Run id from extract phase")
    p_loci.add_argument("--locus_gap_bp", type=int, default=5000)
    p_loci.set_defaults(func=cmd_loci)

    # summarize command placeholder wired below after function is defined in summary module
    p_sum = sub.add_parser("summarize", help="Compute locus metrics and flags from extracted tables and loci")
    _add_common_out_args(p_sum)
    p_sum.add_argument("--run_id", required=True)
    p_sum.add_argument("--evidence_class_map", default=None, help="YAML or 2-col TSV mapping logic_name to class")
    p_sum.add_argument("--locus_gap_bp", type=int, default=5000)
    p_sum.set_defaults(func=cmd_summarize)

    p_exp = sub.add_parser("export", help="Export BED review sets from locus summary")
    _add_common_out_args(p_exp)
    p_exp.add_argument("--run_id", required=True)
    p_exp.add_argument("--no_core_min_score", type=float, default=0.6)
    p_exp.add_argument("--top_n", type=int, default=500)
    p_exp.add_argument("--bed_pad_bp", type=int, default=0)
    p_exp.add_argument("--strand_filter", choices=["none", "same_empty", "both_empty"], default="same_empty",
                       help="Filter missing set by core presence on same/both strands")
    p_exp.set_defaults(func=cmd_export)

    p_lmap = sub.add_parser("layer-map", help="Parse LayerAnnotationStatic.pm into a YAML mapping (biotype→tier)")
    p_lmap.add_argument("--config_key", default="mammals_basic")
    p_lmap.add_argument("--layer_static_url", default="https://raw.githubusercontent.com/Ensembl/ensembl-analysis/refs/heads/main/modules/Bio/EnsEMBL/Analysis/Hive/Config/LayerAnnotationStatic.pm")
    p_lmap.add_argument("--out", required=True)
    def _do_lmap(a):
        text = fetch_layer_static(a.layer_static_url)
        lm = parse_config_block(text, a.config_key)
        mapping = lm.biotype_to_tier()
        yaml_text = to_yaml(mapping)
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write(yaml_text)
        print(f"[layer-map] wrote {len(mapping)} biotypes for '{a.config_key}' → {a.out}")
        return 0
    p_lmap.set_defaults(func=_do_lmap)

    p_rep = sub.add_parser("report", help="Produce retention/loss matrices by layer biotype/tier and a biotype crosswalk")
    _add_common_out_args(p_rep)
    p_rep.add_argument("--run_id", required=True)
    p_rep.add_argument("--layer_map", required=True, help="YAML from layer-map (biotype→tier)")
    p_rep.add_argument("--representation", choices=["overlap","intron_subset"], default="intron_subset",
                       help="Representation criterion: overlap (looser) or intron_subset (stricter, recommended)")
    def _do_report(a):
        out_root = os.path.join(a.output_dir, a.run_id)
        extract_dir = os.path.join(out_root, "extract")
        loci_dir = os.path.join(out_root, "loci")
        rep_dir = os.path.join(out_root, "report")
        ensure_dir(rep_dir)
        # Load data
        def _load(path):
            if a.format == "parquet" and os.path.exists(path + ".parquet"):
                return pd.read_parquet(path + ".parquet")
            sep = "\t" if a.format == "tsv" else ","
            df = pd.read_csv(path + f".{a.format}", sep=sep, low_memory=False)
            if "seq_region_name" in df.columns:
                df["seq_region_name"] = df["seq_region_name"].astype("string")
            return df
        gene_map = _load(os.path.join(loci_dir, "gene_to_locus"))
        core_tx = _load(os.path.join(extract_dir, "core_transcripts"))
        layer_tx = _load(os.path.join(extract_dir, "layer_transcripts"))
        layer_tr = _load(os.path.join(extract_dir, "layer_translations"))
        core_ex = _load(os.path.join(extract_dir, "core_exons")) if a.representation == "intron_subset" else None
        layer_ex = _load(os.path.join(extract_dir, "layer_exons")) if a.representation == "intron_subset" else None
        # Load layer map
        m = load_mapping(a.layer_map)
        by_bt, by_tier, cw, by_logic = compute_retention_by_biotype(
            gene_map, core_tx, layer_tx, layer_tr, m,
            representation=a.representation,
            core_exons=core_ex, layer_exons=layer_ex,
        )
        by_bt.to_csv(os.path.join(rep_dir, "source_retention_by_layer_biotype.tsv"), sep="\t", index=False)
        by_tier.to_csv(os.path.join(rep_dir, "source_retention_by_tier.tsv"), sep="\t", index=False)
        cw.to_csv(os.path.join(rep_dir, "biotype_crosswalk.tsv"), sep="\t", index=False)
        by_logic.to_csv(os.path.join(rep_dir, "source_retention_by_logic_name.tsv"), sep="\t", index=False)
        print(f"[report] wrote matrices → {rep_dir}")
        return 0
    p_rep.set_defaults(func=_do_report)

    p_comp = sub.add_parser("compete", help="Per-locus competitive evidence analysis (which tiers win; why no build)")
    _add_common_out_args(p_comp)
    p_comp.add_argument("--run_id", required=True)
    p_comp.add_argument("--layer_map", required=True, help="YAML from layer-map (biotype→tier)")
    def _do_compete(a):
        out_root = os.path.join(a.output_dir, a.run_id)
        extract_dir = os.path.join(out_root, "extract")
        loci_dir = os.path.join(out_root, "loci")
        comp_dir = os.path.join(out_root, "compete")
        ensure_dir(comp_dir)
        def _load(path):
            if a.format == "parquet" and os.path.exists(path + ".parquet"):
                return pd.read_parquet(path + ".parquet")
            sep = "\t" if a.format == "tsv" else ","
            df = pd.read_csv(path + f".{a.format}", sep=sep, low_memory=False)
            if "seq_region_name" in df.columns:
                df["seq_region_name"] = df["seq_region_name"].astype("string")
            return df
        loci = _load(os.path.join(loci_dir, "loci.strict"))
        mapping = _load(os.path.join(loci_dir, "gene_to_locus"))
        core_tx = _load(os.path.join(extract_dir, "core_transcripts"))
        layer_tx = _load(os.path.join(extract_dir, "layer_transcripts"))
        layer_tr = _load(os.path.join(extract_dir, "layer_translations"))
        b2t = load_mapping(a.layer_map)
        per_locus, agg = analyze_competition(loci, mapping, core_tx, layer_tx, layer_tr, b2t)
        per_locus.to_csv(os.path.join(comp_dir, "per_locus_competition.tsv"), sep="\t", index=False)
        agg.to_csv(os.path.join(comp_dir, "competitive_summary.tsv"), sep="\t", index=False)
        print(f"[compete] wrote per_locus_competition.tsv and competitive_summary.tsv → {comp_dir}")
        return 0
    p_comp.set_defaults(func=_do_compete)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

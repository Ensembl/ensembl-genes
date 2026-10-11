# main.py
from __future__ import annotations
import argparse
from mapper import map_ids

def cli():
    p = argparse.ArgumentParser("ensembl-idmap")
    p.add_argument("--ref-fasta", required=True)
    p.add_argument("--ref-gff", required=True)
    p.add_argument("--target-fasta", required=True)
    p.add_argument("--target-gff", required=True)
    p.add_argument("--output-gff", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--identity-min", type=float, default=0.80)
    args = p.parse_args()
    map_ids(
        ref_fasta=args.ref_fasta, ref_gff=args.ref_gff,
        tgt_fasta=args.target_fasta, tgt_gff=args.target_gff,
        out_gff=args.output_gff, report_txt=args.report,
        threads=args.threads, identity_min=args.identity_min
    )

if __name__ == "__main__":
    cli()

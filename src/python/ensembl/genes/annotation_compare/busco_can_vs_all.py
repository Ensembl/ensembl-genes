#!/usr/bin/env python3
"""
BUSCO Canonical vs All-Proteins Diagnostic (GFF-based)
Diagnose why canonical protein BUSCOs are lower than all-proteins.
Uses GFF to detect missing BUSCOs and potential extra introns in canonical transcripts.
"""

import argparse
from collections import defaultdict
import pandas as pd

# -----------------------------
# BUSCO parsing
# -----------------------------
def parse_busco_table(busco_file):
    """
    Parse BUSCO full_table.tsv
    Returns dict: BUSCO_ID -> {status, protein_id (version stripped)}
    """
    buscos = {}
    with open(busco_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            prot_id = parts[2] if len(parts) > 2 and parts[2] != "-" else None
            if prot_id:
                prot_id = prot_id.split(".")[0]  # strip version suffix
            buscos[parts[0]] = {
                "status": parts[1],
                "protein_id": prot_id
            }
    return buscos

# -----------------------------
# GFF parsing
# -----------------------------
def parse_gff(gff_file):
    """
    Parse GFF3 and build mappings:
        prot_to_tx: protein_id -> transcript_id
        tx_to_gene: transcript_id -> gene_id
        tx_cds_len: transcript_id -> total CDS length
        tx_cds_count: transcript_id -> number of CDS blocks
    """
    prot_to_tx = {}
    tx_to_gene = {}
    tx_cds_len = defaultdict(int)
    tx_cds_count = defaultdict(int)

    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 9:
                continue
            feature_type = cols[2]
            start, end = int(cols[3]), int(cols[4])
            attrs = dict(a.split("=", 1) for a in cols[8].split(";") if "=" in a)

            if feature_type in ("mRNA", "transcript"):
                tx = attrs.get("ID")
                gene = attrs.get("Parent")
                prot = attrs.get("protein_id")
                if prot:
                    prot = prot.split(".")[0]  # strip version
                if tx and gene:
                    tx_to_gene[tx] = gene
                if tx and prot:
                    prot_to_tx[prot] = tx

            elif feature_type == "CDS":
                parents = attrs.get("Parent", "").split(",")
                prot = attrs.get("protein_id")
                if prot:
                    prot = prot.split(".")[0]
                for parent in parents:
                    tx_cds_len[parent] += (end - start + 1)
                    tx_cds_count[parent] += 1
                    if prot:
                        prot_to_tx[prot] = parent

    return prot_to_tx, tx_to_gene, tx_cds_len, tx_cds_count

# -----------------------------
# Diagnosis logic
# -----------------------------
def diagnose(busco_can, busco_all, prot_to_tx, tx_to_gene, tx_cds_len, tx_cds_count):
    """
    Returns DataFrame with BUSCO diagnostic info.
    """
    complete_all = {b for b, v in busco_all.items() if v["status"] in ("Complete", "Duplicated")}
    complete_can = {b for b, v in busco_can.items() if v["status"] in ("Complete", "Duplicated")}

    lost_in_canonical = complete_all - complete_can
    missing_everywhere = set(busco_all.keys()) - complete_all

    rows = []

    def row_for_busco(bid, category):
        can = busco_can.get(bid, {})
        allp = busco_all.get(bid, {})

        prot_can = can.get("protein_id")
        prot_all = allp.get("protein_id")

        tx_can = prot_to_tx.get(prot_can)
        tx_all = prot_to_tx.get(prot_all)

        gene = tx_to_gene.get(tx_all) or tx_to_gene.get(tx_can)

        reasons = []

        # CDS info for intron detection
        cds_count_can = tx_cds_count.get(tx_can)
        cds_count_all = tx_cds_count.get(tx_all)
        cds_len_can = tx_cds_len.get(tx_can)
        cds_len_all = tx_cds_len.get(tx_all)

        if category == "Lost_in_canonical":
            reasons.append("Rescued by non-canonical isoform")
            if cds_count_can and cds_count_all and cds_count_can > cds_count_all:
                reasons.append(f"Canonical transcript has more CDS blocks ({cds_count_can} vs {cds_count_all}) → possible inserted intron")
            if cds_len_can and cds_len_all and cds_len_can < cds_len_all * 0.8:
                reasons.append(f"Canonical CDS shorter than rescued ({cds_len_can} vs {cds_len_all})")
            if not tx_can:
                reasons.append("Canonical transcript not linked to BUSCO protein")
        elif category == "Missing_everywhere":
            reasons.append("BUSCO missing in all isoforms; likely fragmented gene or assembly issue")

        return {
            "BUSCO_ID": bid,
            "Category": category,
            "Canonical_status": can.get("status"),
            "AllProteins_status": allp.get("status"),
            "Canonical_protein": prot_can,
            "Rescuing_protein": prot_all,
            "Gene_ID": gene,
            "Canonical_CDS_count": cds_count_can,
            "Rescued_CDS_count": cds_count_all,
            "Canonical_CDS_length": cds_len_can,
            "Rescued_CDS_length": cds_len_all,
            "Likely_reason": "; ".join(reasons)
        }

    for b in sorted(lost_in_canonical):
        rows.append(row_for_busco(b, "Lost_in_canonical"))
    for b in sorted(missing_everywhere):
        rows.append(row_for_busco(b, "Missing_everywhere"))

    return pd.DataFrame(rows)

# -----------------------------
# CLI
# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagnose canonical vs all-protein BUSCOs using GFF")
    parser.add_argument("--busco_canonical", required=True)
    parser.add_argument("--busco_all", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--out", default="busco_canonical_diagnostic.csv")
    args = parser.parse_args()

    print("Parsing BUSCO tables...")
    busco_can = parse_busco_table(args.busco_canonical)
    busco_all = parse_busco_table(args.busco_all)

    print("Parsing GFF...")
    prot_to_tx, tx_to_gene, tx_cds_len, tx_cds_count = parse_gff(args.gff)

    print("Running GFF-based diagnostic...")
    df = diagnose(busco_can, busco_all, prot_to_tx, tx_to_gene, tx_cds_len, tx_cds_count)

    df.to_csv(args.out, index=False)
    print(f"✓ Written {args.out} ({len(df)} BUSCOs diagnosed)")
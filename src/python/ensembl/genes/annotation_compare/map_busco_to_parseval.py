
import sys
import re

if len(sys.argv) != 3:
    sys.exit("Usage: map_busco_to_parseval.py busco_missing.tsv parseval.txt")

busco_map = {}

# BUSCO table: BUSCO_ID  Protein_ID ...
with open(sys.argv[1]) as fh:
    next(fh)
    for line in fh:
        cols = line.rstrip().split("\t")
        busco, protein = cols[0], cols[1]
        gene = protein.split("-")[0]
        busco_map[gene] = busco

current_locus = None
ref_genes = []

with open(sys.argv[2]) as fh:
    for line in fh:
        if line.startswith("|---- Locus:"):
            current_locus = line.strip().replace("|---- ", "")
            ref_genes = []

        elif line.strip().startswith("|    FUN_"):
            gene = line.strip().split()[0].replace("|", "")
            ref_genes.append(gene)

        elif line.startswith("|----------"):
            for g in ref_genes:
                if g in busco_map:
                    print(
                        busco_map[g],
                        g,
                        current_locus,
                        sep="\t"
                    )
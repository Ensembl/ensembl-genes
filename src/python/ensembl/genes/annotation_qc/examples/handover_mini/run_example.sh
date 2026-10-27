#!/usr/bin/env bash
# SYNTHETIC handover example: annotations -> comparator -> output bundles -> dashboard.
# Usage (from a copy of this folder, with annotation-qc on PATH):
#     bash run_example.sh            # comparisons, import, build
#     annotation-qc dashboard --config configs/mini_species.json --port 8765
# Mirrors the layout of the mammalian data: results/<run>/<variant>/<policy>/.
set -euo pipefail

policy_selections() {  # policy id -> "reference_selection query_selection"
  case "$1" in
    A) echo "all all" ;;
    B) echo "longest_cds longest_cds" ;;
    C) echo "canonical longest_cds" ;;
  esac
}

compare() {  # run_id variant policy query [extra pairwise-compare options...]
  local run=$1 variant=$2 policy=$3 query=$4; shift 4
  local sel; read -r -a sel <<< "$(policy_selections "$policy")"
  annotation-qc pairwise-compare \
    --reference inputs/reference.gff3 --query "$query" \
    --genome inputs/genome.fa --genome-mismatch error \
    --evaluation-mode protein_coding --reference-transcript-biotypes protein_coding \
    --reference-transcript-selection "${sel[0]}" --query-transcript-selection "${sel[1]}" \
    --analysis-scope chromosomes=1 \
    --export-evaluated-models ids \
    "$@" \
    --outdir "results/$run/$variant/$policy"
}

for policy in A B; do   # C needs Ensembl_canonical tags, which this example does not have
  compare mini_internal stop_harmonised "$policy" inputs/focus.gff3 --add-stop-codon --stop-codon-genetic-code 1
  compare mini_other original "$policy" inputs/other.gff3
done

annotation-qc benchmark validate --config configs/mini_species.json --prepare
annotation-qc benchmark import --config configs/mini_species.json
annotation-qc benchmark status --config configs/mini_species.json
annotation-qc benchmark build --config configs/mini_species.json

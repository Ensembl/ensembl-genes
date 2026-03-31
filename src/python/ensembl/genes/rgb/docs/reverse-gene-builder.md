# Reverse Gene Builder: engineering brief (Phase 1)

This document specifies a locus‑centric, genome‑wide audit that accounts for how upstream “layer” evidence is represented, compressed, or ignored in the “core” gene set. It implements M1–M2 first: full extraction and locus construction, with clear inputs, outputs, and acceptance criteria.

## Purpose
- Provide a comprehensive, reproducible accounting between layer evidence and core output to explain low core gene counts despite abundant evidence.

## Scope (Phase 1)
- In: two databases (core and layer), genome‑wide extraction; unified loci (strict and expanded); per‑locus counts and spans; review‑ready BEDs (scaffolding); aggregate toplines.
- Out: re‑alignment, de novo calling, per‑base browsers, automated corrections.

## Definitions
- Core entity: gene/transcript from core DB.
- Layer entity: gene/transcript from layer DB (provenance via `analysis.logic_name`).
- Evidence class: normalized category for `logic_name` (transcriptomic, homology, projection/synteny, ab initio, other).
- Locus
  - Strict: on (seq_region, strand), union of overlapping intervals (≥1 bp).
  - Expanded: merge strict loci if nearest edges are within `locus_gap_bp` (default 5,000 bp).

## Inputs
- DB: `host`, `port`, `user`, `password`, `core_db`, `layer_db`.
- Filters: `coord_system_name` (optional); `seq_region_name` allowlist (optional).
- Locus params: `locus_gap_bp` (default 5000), `locus_mode` (strict|expanded).
- Evidence classes: `evidence_class_map` YAML (logic_name → class).
- Output dir: absolute path.

## Extraction (M1)
- Tables: `gene`, `transcript`, `translation` (presence only), `analysis`, `seq_region`, `coord_system`.
- Never assume `seq_region_id` parity across DBs; join through `seq_region.name`.
- Restrict by `coord_system.name` when provided.
- Outputs (TSV by default; Parquet if `pyarrow` is available):
  - `extract/core_genes.*`: `gene_id,stable_id,seq_region_name,start,end,strand,biotype,canonical_transcript_id,logic_name`.
  - `extract/layer_genes.*`: same columns, layer DB.

## Locus construction (M2)
- Build from union of core+layer genes per (seq_region_name, strand).
- Strict: merge any overlapping intervals.
- Expanded: merge adjacent strict loci if gap ≤ `locus_gap_bp` (distance = `next.start - prev.end - 1`).
- Deterministic `locus_id`: `"{seq}:{strand}:{start}:{end}:{ordinal}"`.
- Outputs:
  - `loci/loci.strict.*` and `loci/loci.expanded.*` with columns `seq_region_name,strand,locus_start,locus_end,locus_length,core_gene_count,layer_gene_count`.
  - `loci/gene_to_locus.*`: `db_kind(core|layer),gene_id,seq_region_name,strand,start,end,locus_id_strict,locus_id_expanded`.

## Representation metrics (preview)
- Presence: `core_gene_count, layer_gene_count`.
- Span: `locus_length` (strict), reserved columns for coverage (populated in Phase 2).

## CLI
- `rgb extract --host ... --core_db ... --layer_db ... --coord_system_name chromosome --output_dir /abs/out --format tsv`
- `rgb loci --output_dir /abs/out --locus_gap_bp 5000 --run_id <from-extract>`
- `rgb summarize --output_dir /abs/out --run_id <from-extract> --evidence_class_map /abs/map.yaml`
 - `rgb layer-map --config_key mammals_basic --out /abs/out/<run_id>/layer_map.yaml`
- `rgb report --output_dir /abs/out --run_id <from-extract> --layer_map /abs/out/<run_id>/layer_map.yaml`
 - `rgb compete --output_dir /abs/out --run_id <from-extract> --layer_map /abs/out/<run_id>/layer_map.yaml`
 - `rgb export --output_dir /abs/out --run_id <from-extract> --no_core_min_score 0.6 --top_n 500`

## Acceptance (Phase 1–2)
- Deterministic outputs with manifests including parameters and input hashes.
- Human‑scale genome finishes ≤ 2h and ≤ 16 GB RAM (guideline; depends on DB and network).
- Loci counts and simple sanity checks logged per seq_region.
 - Locus summary contains presence, span, compression, coverage, diversity metrics, and coarse diagnostic flags.
- Report includes per-layer-biotype/tier retention matrices and crosswalks.
 - Compete includes per-locus winning tier/logic_name and no-build categorization.

## Risks and mitigations
- Logic name inflation: require mapping file in later phases; log top unmapped sources.
- Locus sensitivity: ship both strict and expanded; report deltas.
- Strand handling: loci never cross strands.

## Next
- Phase 2: per‑locus layer/core coverage, compression, diversity metrics; evidence‑rich no‑core scoring; underbuilt flags.
- Phase 3: merge/readthrough candidate signals and scoring.

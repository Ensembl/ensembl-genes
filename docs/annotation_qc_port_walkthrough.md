# Annotation QC Port: Script → Package

> **Goal of this walkthrough:** explain *why* the code moved the way it did, not just *what* changed.

---

## The One-Line Summary

> We moved from a monolithic script to a layered QC package: parser normalization, reusable metric functions, and runner/report separation — while preserving the legacy stat outputs where needed.

---

## Architecture Shift

### Old shape

```
parse GFF3 manually
  → build Gene / Transcript objects
    → compute stats inside one big function
      → write TSV / JSON / SQL
```

### New shape

```
parse annotation (GFF3 or GTF, gzip-transparent)
  → one standardized DataFrame
    → split by feature type (gene, transcript, exon, CDS, UTR…)
      → call isolated metric functions
        → runner collects results, writes outputs
```

**The key insight:** the data model changed from custom objects to a standardized annotation table. Code movement is secondary to that.

---

## Entry Points

| | Old | New |
|---|---|---|
| Script entry | [`scripts/ensembl_gff_metrics.py:1690`](../scripts/ensembl_gff_metrics.py) | [`annotation_qc/cli.py`](../src/python/ensembl/genes/annotation_qc/cli.py) |
| Orchestration | inline in `main()` | [`runners/genebuild_stats.py`](../src/python/ensembl/genes/annotation_qc/runners/genebuild_stats.py) |

---

## CLI Interaction

The package exposes a single top-level command: **`annotation-qc`**.

### How the CLI is wired

```
cli.py::main()
  ├── creates argparse parser with subparsers
  ├── iterates over RUNNERS = [assess_translation_validity, genebuild_stats]
  └── calls runner.register(subparsers) for each
```

Each runner registers itself as a subcommand and sets `func=_run`. When the user invokes the CLI, `args.func(args)` dispatches to the right runner's `_run`.

### Available subcommands

```bash
annotation-qc --help

# Subcommands:
#   feature-stats    Compute genebuild feature statistics
#   assess-translation-validity  (second runner)
```

### `feature-stats` — minimal invocation

```bash
annotation-qc feature-stats \
  --annotation annotation.gff3 \
  --outdir results/
```

**Always produces:**
- `results/feature_stats.tsv` — long-format summary (category | metric | value) for coding, noncoding, pseudogene
- `results/utr_stats.tsv` — per-transcript 5′/3′ UTR lengths
- `results/cds_utr5_overlap.tsv` — CDS / 5′ UTR exon overlaps

### `feature-stats` — with assembly metadata

```bash
annotation-qc feature-stats \
  --annotation annotation.gff3 \
  --assembly-accession GCA_000001405.29 \
  --scientific-name "Homo sapiens" \
  --taxon-id 9606 \
  --outdir results/
```

When `--assembly-accession` is provided the runner calls `ncbi_assembly_stats()` and writes the legacy compatibility TSVs:
- `coding_stats.tsv`, `noncoding_stats.tsv`, `pseudogene_stats.tsv`, `assembly_stats.tsv`

Add `--json` or `--sql` (with `--species-id`) to also write JSON and Ensembl meta-table INSERT IGNORE SQL.

### `feature-stats` — sequence-based QC

```bash
annotation-qc feature-stats \
  --annotation annotation.gff3 \
  --genome genome.fa \
  --genetic-code 1 \
  --outdir results/
```

Adding `--genome` unlocks two extra metric passes that require sequence but no alignment evidence:
- `results/translation_metrics.tsv` — per-CDS start/stop codon checks, internal stops
- `results/splice_junction_metrics.tsv` — per-intron GT-AG / GC-AG / AT-AC canonicality

`--genetic-code` takes an NCBI table ID (default 1 = standard). Useful for mitochondrial or non-standard genetic codes.

---

## Port Map

| Legacy responsibility | Old location | New location |
|---|---|---|
| Parse annotation file | `parse_ensembl_gff3` L385 | `parse_annotation` + `standardize_annotation` |
| gzip / raw attr handling | `open_maybe_gzip`, `parse_gff3_attrs` | delegated to PyRanges1 |
| Gene / Transcript model | `Gene` / `Transcript` classes L136 | replaced by one standardized DataFrame |
| Canonical transcript selection | `pick_canonical_transcript` L292 | `_canonical_info` (same ranking, over DataFrames) |
| Biotype grouping | `biotype_group` L343 | `_biotype_group` (logic preserved) |
| Intron derivation | `compute_introns` L116 | `_compute_introns` (vectorized) |
| Coding gene stats | `process_coding_genes` L549 | `compute_coding_stats` |
| Non-coding gene stats | `process_non_coding_genes` L732 | `compute_noncoding_stats` |
| Pseudogene stats | `process_pseudogenes` L874 | `compute_pseudogene_stats` |
| NCBI assembly metadata | `fetch_and_parse_ncbi_assembly_stats_report` L1311 | `ncbi_assembly_stats` |
| Output writers | `write_tsv`, `write_json_from_tsv`, `write_meta_sql` | `reports/meta_sql.py` + `_write_compat_reports` |

---

## Walkthrough: One Metric End-to-End

### Old flow for coding stats

1. Parse GFF3 line-by-line → populate `Gene` and `Transcript` objects
2. `pick_canonical_transcript` iterates over `gene.transcripts`
3. `process_coding_genes` loops over genes, accumulates totals in plain variables
4. Returns a flat dict; written directly to TSV from `main()`

### New flow for coding stats

1. `parse_annotation(path)` → one DataFrame via PyRanges1
2. Runner slices it: `gene_df`, `tx_df`, `exon_df`, `cds_df`
3. `_canonical_info(tx_df)` derives canonical transcript per gene using the same ranking rule, as a DataFrame join
4. `compute_coding_stats(gene_df, tx_df, exon_df, cds_df)` returns `(stats_dict, total_cds_len)`
5. Runner calls `_to_long(coding_stats, "coding")` and concatenates all categories into `feature_stats.tsv`

**Why this matters:** `compute_coding_stats` is now a pure function — testable in isolation, no file I/O, no object state.

---

## New Capabilities Beyond the Original Script

These metrics did not exist in the old script. They were only possible once parsing and metrics were separated:

| Metric | Function |
|---|---|
| Per-transcript 5′/3′ UTR lengths | `compute_transcript_utr_stats` |
| CDS / 5′ UTR exon overlap detection | `compute_cds_utr5_overlap` |
| Per-CDS translation QC (start/stop/internal stop) | `compute_translation_metrics` |
| Per-intron splice junction canonicality | `compute_splice_junction_metrics` |

---

## Parsing Layer in Detail

### `parse_annotation` (`parsers/annotation.py:129`)

```
file_path
  ├── .gff3 / .gff3.gz  →  pr.read_gff3()
  └── .gtf  / .gtf.gz   →  pr.read_gtf()
                               ↓
                    standardize_annotation()
                               ↓
                    normalized DataFrame
```

### `standardize_annotation` (`parsers/annotation.py:30`)

Fills in structural columns that downstream metrics rely on, regardless of input format:

| Column filled | Source |
|---|---|
| `gene_id` | from `ID` (gene rows) or `Parent` (transcript rows) |
| `transcript_id` | from `ID` (transcript rows) or `Parent` (child rows) |
| `biotype` | from `gene_type` if `biotype` absent (GTF compat) |
| `phase` | from `Frame` if `phase` absent |
| `feature_length` | `End - Start` |

This is the contract every metric function relies on — one place to fix format differences.

---

## Reporting Layer

`_write_compat_reports` in the runner writes the same TSV/JSON/SQL shape as the old script. This is deliberate: downstream processes that consume those files don't need to change.

The new `feature_stats.tsv` (long format) is additional output — easier to extend without changing column counts.

---

## Closing Point

The port makes parity testing straightforward: run both the old script and `annotation-qc feature-stats` on the same file and diff `coding_stats.tsv`. The metric logic is now independently testable because it is isolated from file I/O and CLI concerns.

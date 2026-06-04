# Annotation QC package design

This package is for annotation quality-control metrics that can be reused across
species, annotation sources, and reporting workflows.  Most of the original code
started as standalone scripts.  The package structure should make those scripts
easier to maintain without hiding the data model from contributors.

## Core idea

Use one standardized annotation table as the shared representation.

```text
annotation file
  -> parsers.annotation.parse_annotation()
  -> standardized pandas DataFrame
  -> runner filters the DataFrame once
  -> metric functions compute from plain DataFrames
  -> report functions write outputs
```

Do not introduce a second annotation object model unless there is a clear need.
The expected working style is DataFrame-based because it is easy to inspect,
debug, and teach.

## Package roles

`parsers/`
: Read input files and normalize their structure.  This includes GFF3/GTF
  annotations, FASTA files, and external tool outputs such as DIAMOND, STAR,
  BUSCO, or prediction TSVs.  Parser code is allowed to understand file-format
  details such as GFF3 `Parent`, GTF attributes, Ensembl prefixes, and optional
  metadata conventions.

`metrics/`
: Compute QC metrics from already parsed data.  Metric functions should accept
  DataFrames or parsed evidence objects and return dictionaries or DataFrames.
  They should avoid reading files, writing files, parsing attributes, or
  depending on CLI arguments.

`runners/`
: Wire parser, filtering, metric computation, and reporting together for a CLI
  workflow.  A runner can decide which feature slices it needs and should filter
  the standardized annotation table once before calling metrics.

`reports/`
: Convert metric outputs into TSV, JSON, SQL, plots, or text summaries.  Report
  code should not recompute metrics.

`metadata/`
: Fetch or derive external metadata, for example NCBI assembly statistics.

## Standard annotation table

`parse_annotation()` should return a pandas DataFrame that keeps parser output
columns and fills common structural columns where possible.

Expected core columns:

```text
Feature
Chromosome
Start
End
Strand
ID
Parent
gene_id
transcript_id
biotype
tag
phase
feature_length
```

Not every annotation source will provide every value.  Missing optional metadata
should remain missing rather than being guessed in metric code.

Parser normalization should handle common structural differences:

- GFF3 and GTF attribute syntax.
- Ensembl-style prefixes such as `gene:` and `transcript:`.
- `gene_type` / `transcript_type` as fallbacks for `biotype`.
- Deriving transcript IDs from child-feature parents where possible.
- Propagating `gene_id` from transcript rows to child features where possible.
- Multiple parent IDs in GFF3 child features.

For multiple-parent child features, prefer a table shape that keeps metric code
simple.  Duplicating one row per transcript parent is usually easier to work
with than requiring every metric to understand list-valued parent columns.

## Runner pattern

Runners should make the feature slices they pass to metrics explicit.

```python
annotation = parse_annotation(annotation_path)

gene_df = annotation[annotation["Feature"].isin(GENE_FEATURE_TYPES)]
tx_df = annotation[annotation["Feature"].isin(TRANSCRIPT_FEATURE_TYPES)]
exon_df = annotation[annotation["Feature"] == "exon"]
cds_df = annotation[annotation["Feature"] == "CDS"]
utr5_df = annotation[annotation["Feature"] == "five_prime_UTR"]
utr3_df = annotation[annotation["Feature"] == "three_prime_UTR"]

coding_stats, total_cds = compute_coding_stats(gene_df, tx_df, exon_df, cds_df)
utr_stats = compute_transcript_utr_stats(tx_df, gene_df, utr5_df, utr3_df)
translation = compute_translation_metrics(cds_df, fasta)
splice = compute_splice_junction_metrics(exon_df, fasta)
```

This makes dependencies visible.  For example, coding gene stats need CDS
features, but they also need genes, transcripts, exons, and introns.

## Metric function guidelines

Prefer metric functions that:

- Accept DataFrames or parsed evidence data.
- Return plain Python dictionaries or pandas DataFrames.
- Keep file IO and CLI concerns outside the metric layer.
- Treat missing optional metadata gracefully.
- Use parser-normalized columns instead of reparsing GFF/GTF attributes.
- Preserve behavior against the original script with focused parity tests.

Avoid metric functions that:

- Read annotation files directly.
- Write output files directly.
- Depend on global CLI state.
- Reimplement GFF/GTF parsing.
- Require contributors to understand a deep object model before they can add a
  simple metric.

## Porting script-style code

When moving an existing script into the package:

1. Identify its inputs: annotation, FASTA, external evidence, metadata, or report
   files.
2. Move input parsing into `parsers/` or reuse an existing parser.
3. Move calculations into `metrics/` as functions over DataFrames or parsed
   evidence.
4. Move output writing into `reports/` when the format is reusable.
5. Put CLI wiring in `runners/`.
6. Add tests for the metric functions and at least one parity check against the
   original script output where practical.

## Current script families

Genebuild feature stats
: Coding, non-coding, pseudogene, UTR, CDS/UTR overlap, translation, and splice
  summaries.  These belong in `metrics/features.py` with CLI wiring in
  `runners/genebuild_stats.py`.

Feature-wide metrics
: Exon/intron length distributions, phase summaries, UTR presence, feature GC,
  frame consistency, and flagged short/long features.  The script logic maps to
  `metrics/feature_wide.py`; CLI wiring should live in a runner.

Genome-wide metrics
: Genome size, soft-masked bases, annotation coverage, gene density, and
  stranded overlap counts.  The script logic maps to `metrics/genome_wide.py`;
  CLI wiring should live in a runner.

CDS/start-stop validation
: Start/stop codon checks are event-level metrics.  DIAMOND-backed validation
  should use a parser for DIAMOND output and a metric module for evidence-based
  classification.

UTR completeness
: STAR splice-junction parsing belongs in `parsers/`.  UTR completeness and
  junction-support calculations belong in a metric module.  Plotting belongs in
  `reports/`.

Gene boundary metrics
: Reference annotation parsing should use the standard annotation parser where
  possible.  Prediction TSV parsing should be separate.  Boundary matching and
  tolerance-based scoring belong in a metric module.

Data split utilities
: These are not annotation QC metrics unless they are needed by a QC/reporting
  workflow.  Keep them outside this package unless that connection becomes
  explicit.

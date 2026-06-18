# genome_wide_metrics.py

`genome_wide_metrics.py` calculates genome-wide QC metrics for selected
features in GFF3/GTF genome annotations.

The script reports feature density, genome coverage, soft-masked repeat
annotation, and stranded feature-overlap metrics to the terminal. It can also
write the same metric rows to a CSV file, write a table of soft-masked/feature
overlaps, and write a table of overlapping feature pairs for follow-up
inspection.

## Requirements

- Python 3.9 or newer
- No third-party Python packages
- A GFF3 or GTF annotation file, optionally gzip-compressed
- Either a genome size in base pairs or a genome FASTA file
- A soft-masked genome FASTA file, optionally gzip-compressed, for masking
  metrics

The FASTA and annotation seqids must refer to the same sequence regions for
soft-masked/feature overlap metrics. FASTA seqids are read from the first token
after `>`.

## Basic Usage

Calculate metrics with an explicit genome size:

```bash
python3 genome_wide_metrics.py annotation.gtf --genome-size 2500000000
```

Calculate genome size and soft-masking metrics from a FASTA file:

```bash
python3 genome_wide_metrics.py annotation.gtf --fasta genome.fa
```

Write metrics to a CSV file as well as the terminal:

```bash
python3 genome_wide_metrics.py annotation.gtf \
  --genome-size 2500000000 \
  --csv-output
```

Filter to protein-coding canonical genes and write CSV output:

```bash
python3 genome_wide_metrics.py annotation.gtf \
  --fasta genome.fa \
  --biotype protein_coding \
  --canonical-only \
  --csv-output
```

Write a table of soft-masked regions overlapping selected features:

```bash
python3 genome_wide_metrics.py annotation.gtf \
  --fasta genome.fa \
  --masked-overlaps-output masked_overlaps.tsv
```

Write a table of overlapping feature pairs:

```bash
python3 genome_wide_metrics.py annotation.gtf \
  --genome-size 2500000000 \
  --feature-overlaps-output feature_overlaps.tsv
```

The overlap outputs are tab-separated files, so `.tsv` is the recommended
extension.

## Options

| Option | Description |
| --- | --- |
| `gff` | Input GFF3/GTF annotation file. |
| `--genome-size N` | Genome size in base pairs. Required if `--fasta` is not supplied. |
| `--fasta FASTA` | Soft-masked genome FASTA used for genome size and masking metrics. |
| `--feature-type VALUE` | Annotation feature type to count. Default: `gene`. |
| `--canonical-only` | Include only features belonging to canonical transcripts. |
| `--biotype VALUE` | Include features with one of these biotypes. Can be repeated or comma-separated. |
| `--masked-overlaps-output PATH` | Write a TSV of soft-masked/feature overlap coordinates. Requires `--fasta`. |
| `--feature-overlaps-output PATH` | Write a TSV of overlapping feature pairs and their overlap coordinates. |
| `--csv-output` | Write terminal metrics to a generated CSV file. |

## CSV Output

When `--csv-output` is used, the output file is named:

```text
<annotation_base>.<biotype>.<canonical>.features.csv
```

Examples:

```text
annotation.all.all.features.csv
annotation.protein_coding.all.features.csv
annotation.protein_coding.canonical.features.csv
annotation.lncrna-protein_coding.canonical.features.csv
```

The CSV is written in the current working directory. The `annotation_base` is
the annotation filename with `.gtf`, `.gff`, `.gff3`, and optional `.gz` suffixes
removed.

If multiple `--biotype` values are supplied, they are sorted and joined with
hyphens. If no biotype filter is supplied, the biotype tag is `all`. If
`--canonical-only` is not used, the canonical tag is `all`.

Filename tags are sanitized so unusual characters do not break filenames.

The CSV has two columns:

```text
metrics_name,metrics_value
```

## Reported Metrics

The terminal report and optional CSV include:

- Genome size in base pairs
- Selected feature count
- Canonical-only filter status, when `--canonical-only` is requested
- Biotype filter, when `--biotype` is used
- Selected feature density per megabase
- Percentage of the genome covered by merged selected features
- Total soft-masked bases and percentage of genome soft-masked, when FASTA
  masking is available
- Percentage of soft-masked regions annotated by selected features
- Overlapping selected feature pairs on the same strand
- Antisense selected feature pairs on opposite strands

## Masked Overlaps Output

`--masked-overlaps-output` writes a TSV with these columns:

```text
chrom overlap_start overlap_end mask_start mask_end feature_id feature_start feature_end feature_strand overlap_bp
```

The rows describe each overlap between a soft-masked lowercase sequence run and
a selected annotation feature.

## Feature Overlaps Output

`--feature-overlaps-output` writes a TSV with these columns:

```text
chrom overlap_start overlap_end feature1_id feature1_start feature1_end feature1_strand feature2_id feature2_start feature2_end feature2_strand relation overlap_bp
```

The `relation` column is one of:

- `same_strand`: both overlapping features are on the same strand
- `antisense`: overlapping features are on opposite `+` and `-` strands
- `other_strand`: overlapping features have another strand relationship

## Notes

- GFF3/GTF coordinates are parsed as 1-based inclusive input coordinates.
- Internal calculations use zero-based half-open intervals.
- Overlap TSV coordinates are written as zero-based half-open intervals.
- If both `--genome-size` and `--fasta` are supplied, `--genome-size` is used
  for density and percentage denominators while the FASTA is still scanned for
  soft-masking metrics.
- Soft-masked lowercase runs are counted as annotated when any part of the run
  overlaps a selected feature.
- If no soft-masking is detected, the masked annotation metric is reported as
  `N/A`.
- Input files and optional TSV outputs ending in `.gz` are read or written as
  gzip-compressed files.


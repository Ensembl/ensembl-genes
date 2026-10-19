# feature_wide_metrics.py

`feature_wide_metrics.py` calculates feature-level QC metrics for GFF3/GTF
genome annotations using a matching genome FASTA file.

The script reports exon, intron, phase, GC, UTR, and CDS frame metrics to the
terminal. It can also write the same metric rows to a CSV file and write a
separate table of flagged features for follow-up inspection.

## Requirements

- Python 3.9 or newer
- No third-party Python packages
- A GFF3 or GTF annotation file, optionally gzip-compressed
- A genome FASTA file, optionally gzip-compressed

The FASTA and annotation seqids must refer to the same sequence regions. The
script handles standard first-token FASTA IDs and Ensembl-style FASTA headers
such as:

```text
>primary_assembly:Sscrofa_HxYL_paternal_Hampshire_1.0:1:1:288660608:1 primary_assembly 1
```

For this header, the script can match annotation seqid `1` as well as the full
first FASTA token.

## Basic Usage

```bash
python3 feature_wide_metrics.py annotation.gtf --fasta genome.fa
```

Write metrics to a CSV file as well as the terminal:

```bash
python3 feature_wide_metrics.py annotation.gtf --fasta genome.fa --csv-output
```

Filter to protein-coding canonical transcripts and write CSV output:

```bash
python3 feature_wide_metrics.py annotation.gtf \
  --fasta genome.fa \
  --biotype protein_coding \
  --canonical-only \
  --csv-output
```

Write a table of potentially problematic features:

```bash
python3 feature_wide_metrics.py annotation.gtf \
  --fasta genome.fa \
  --flagged-features flagged.tsv
```

`--flagged-features` writes a tab-separated file, so `.tsv` is the recommended
extension.

## Options

| Option | Description |
| --- | --- |
| `gff` | Input GFF3/GTF annotation file. |
| `--fasta FASTA` | Required genome FASTA file used for exon GC metrics. |
| `--format {auto,gff3,gff,gtf}` | Annotation format. Default: `auto`. |
| `--gc-bin-size N` | GC histogram bin size in percentage points. Default: `10`. |
| `--biotype VALUE` | Include transcripts with one of these biotypes. Can be repeated or comma-separated. |
| `--exclude-biotype VALUE` | Exclude transcripts with one of these biotypes. Can be repeated or comma-separated. |
| `--canonical-only` | Include only transcripts tagged as canonical, Ensembl_canonical, or MANE_Select. |
| `--csv-output` | Write terminal metrics to a generated CSV file. |
| `--flagged-features PATH` | Write a TSV of features/transcripts that may need review. |
| `--min-intron-length N` | Flag introns shorter than this threshold. Default: `20`. |
| `--min-exon-length N` | Flag exons shorter than this threshold. Default: `3`. |

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
```

The CSV is written in the same directory as the input annotation file. The
`annotation_base` is the annotation filename with `.gtf`, `.gff`, `.gff3`, and
optional `.gz` suffixes removed.

If multiple `--biotype` values are supplied, they are sorted and joined with
hyphens. If no biotype filter is supplied, the biotype tag is `all`. If
`--canonical-only` is not used, the canonical tag is `all`.

The filename biotype tag is based on `--biotype`; excluded biotypes are reported
inside the metrics table but are not added to the filename.

The CSV has two columns:

```text
metrics_name,metrics_value
```

## Reported Metrics

The terminal report and optional CSV include:

- Annotation format and active transcript filters
- Exon count and transcript counts
- Phase distribution, using exon phase when available and CDS phase otherwise
- Exon length summary and histogram
- Exon GC summary and histogram
- Number of exons skipped for GC calculation
- Annotation seqids missing from FASTA, when detected
- Intron length summary and histogram
- 5' and 3' UTR transcript presence, length, and exon-count summaries
- CDS frame-consistent and frame-inconsistent transcript counts
- Flagged-feature output path and thresholds, when `--flagged-features` is used

## Flagged Features Output

`--flagged-features` writes a TSV with these columns:

```text
flag_type transcript_id feature_type seqid start end strand length details
```

Flag types include:

- `frame_inconsistent_transcript`: CDS length is not divisible by 3
- `short_exon`: exon length is below `--min-exon-length`
- `short_intron`: intron length is below `--min-intron-length`
- `cds_without_usable_phase`: CDS phase is not `0`, `1`, or `2`

## Notes

- GFF3/GTF coordinates are parsed as 1-based inclusive input coordinates.
- Internal calculations use zero-based half-open intervals.
- Flagged feature coordinates are written back as 1-based inclusive positions.
- Only `A`, `C`, `G`, and `T` bases count toward exon GC percentages.
- Soft-masked lowercase bases are handled correctly because sequence chunks are
  uppercased before counting.

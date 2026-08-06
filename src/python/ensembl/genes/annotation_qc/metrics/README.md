# metrics

## Purpose

The `metrics` package contains all QC calculations.

A metric receives already parsed annotation data and returns calculated values.

Metrics should **never** read files or generate reports.

---

## Responsibilities

Typical responsibilities include:

* counting genes or transcripts
* feature-wide statistics
* exon distributions
* transcript statistics
* QC pass/fail calculations
* summary statistics

---

## What does NOT belong here?

The following belong elsewhere:

| Task             | Package    |
| ---------------- | ---------- |
| Reading GFF3/GTF | `parsers/` |
| Writing JSON/TSV | `reports/` |
| CLI handling     | `runners/` |
| Constants        | `config/`  |

---

## Input

Metric functions receive standardized annotation data produced by the parsers.

Example

```
annotation_df

Feature
Start
End
ID
Parent
gene_id
transcript_id
...
```

---

## Output

Metric functions should return plain Python objects such as

* dictionaries
* dataclasses
* pandas DataFrames

They should not write files.

---

## Example

```
annotation_df
        │
        ▼
compute_feature_wide_metrics()
        │
        ▼
{
    gene_count: ...
    transcript_count: ...
}
```

---

## Adding a new metric

Before creating a new module ask:

* Is this computing a statistic?
* Does it avoid file I/O? 
* Can it be unit tested independently?

If the answer is yes, it belongs here.

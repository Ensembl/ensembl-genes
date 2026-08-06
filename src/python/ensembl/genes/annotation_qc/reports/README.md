# reports

## Purpose

The `reports` package is responsible for presenting QC results.

Reports transform metrics into formats that can be consumed by users or downstream pipelines, such as TSV, JSON, HTML or summary tables.

This package should never calculate metrics itself.

---

# Responsibilities

Reports are responsible for:

* Formatting QC results
* Exporting JSON
* Exporting TSV
* Producing summary tables
* Creating human-readable output
* Generating visualisations (if required)

Reports consume already computed metrics.

---

# What does NOT belong here?

The following should never appear in this package:

* Annotation parsing
* QC calculations
* Feature validation
* Biological decision logic
* CLI argument parsing

Reports should only present existing information.

---

# Inputs

Reports receive outputs from the metrics package.

Examples include

* dictionaries
* pandas DataFrames
* dataclasses
* metric collections

Example

```python
{
    "gene_count": 24561,
    "transcript_count": 98214,
    "mean_exons": 8.7
}
```

---

# Outputs

Reports may generate

* TSV
* JSON
* CSV
* Markdown
* HTML
* plots
* summary text

Every report should be reproducible from the same metric object.

---

# Design principles

Reports should

* contain no biological logic
* avoid recomputing statistics
* separate formatting from calculations
* support multiple output formats

The workflow is

```
Metrics

↓

Formatting

↓

Output file
```

---

# Adding a new report

When implementing a new report

1. Accept metric objects as input.
2. Format them.
3. Export to the desired format.
4. Add unit tests.
5. Update documentation.

---

# Common mistakes

❌ Recalculating statistics

❌ Filtering annotation data

❌ Reading GFF files

❌ Calling parsers directly

---

# Testing

Reports should verify

* formatting
* output schema
* file generation
* expected values

without testing the metric calculations themselves.

---

# Related packages

* `metrics/`
* `runners/`

Reports should remain independent from parsing and QC implementation details.

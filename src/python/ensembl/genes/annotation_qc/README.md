# Annotation QC onboarding guideline

This guide explains how to structure code inside `annotation_qc/` 

## 1. Core rule

The package follows a single flow:

`input file -> parser -> standardized pandas DataFrame -> runner -> metric functions -> report output`

The important rule is:

* **parsers** normalize input;
* **metrics** compute QC logic;
* **runners** connect everything for CLI use;
* **reports** format results for export;
* **metadata** provides external reference data when needed.

Do not create a second annotation object model unless there is a strong reason.

## 2. Folder responsibilities

### `parsers/`

Use this folder for code that reads files and converts them into a standard internal representation.

Put here:

* GFF3/GTF parsing
* FASTA parsing
* DIAMOND/STAR/BUSCO/prediction TSV parsing
* attribute normalization
* ID cleanup
* structural normalization such as `gene_id`, `transcript_id`, `Parent`, `biotype`

Do not put here:

* QC decisions
* statistics calculation
* report writing
* CLI parsing

### `metrics/`

Use this folder for pure QC calculations.

Put here:

* feature-level metrics
* gene/transcript/exon summary calculations
* distribution checks
* classification logic
* tolerance-based scoring
* parity-safe refactors of old scripts

Rules for this folder:

* accept DataFrames or already parsed evidence objects
* return dicts or DataFrames
* avoid file I/O
* avoid CLI logic
* avoid parsing raw file formats again

### `runners/`

Use this folder for workflow orchestration.

Put here:

* CLI entry points
* input selection
* filtering slices from the parsed table
* calling metric functions in the right order
* passing results to reports

A runner should be the place where the workflow is assembled, not where the QC logic lives.

### `reports/`

Use this folder for turning metric results into output.

Put here:

* TSV export
* JSON export
* text summaries
* plots
* tables
* any reusable presentation layer

Do not recompute metrics here.

### `notebook/`

Use this folder for visualising and investigating result.

Put here any reusable notebook ideally one for each type of metrics.


### `config/`

Use this folder for configuration definitions and defaults.

Put here:

* thresholds
* constants
* runtime settings
* feature lists
* shared configuration values

## 3. Practical rule for deciding where code goes

Use this decision rule:

* If it **reads a file** → `parsers/`
* If it **calculates QC values** → `metrics/`
* If it **connects steps together** → `runners/`
* If it **formats output** → `reports/`
* If it **stores constants or thresholds** → `config/`
* If it **retrieves external reference data** → `metadata/`

## 4. Adding a new QC

* When implementing a new QC:
* Identify the required input.
* Reuse an existing parser where possible.
* Implement the calculation inside metrics/.
* Call the metric from a runner.
* Export the result through a report.
* Add tests.
* Update the relevant package documentation.

## 7. Template for a new metric module

Use this template for new files inside `metrics/`:

```python
"""
feature_wide.py

Purpose:
    Compute feature-wide annotation QC metrics from a standardized DataFrame.

Inputs:
    annotation_df: pandas.DataFrame
        Parsed and normalized annotation data.

Outputs:
    dict | pandas.DataFrame
        Metric summary suitable for reporting.
"""

from __future__ import annotations

import pandas as pd


def compute_feature_wide_metrics(annotation_df: pd.DataFrame) -> dict:
    """
    Compute feature-wide QC statistics.

    Parameters
    ----------
    annotation_df
        Standardized annotation table produced by parsers/.

    Returns
    -------
    dict
        Feature-wide metric values.

    Notes
    -----
    This function must not read files, write files, or perform CLI logic.
    """
    # 1. validate required columns
    # 2. derive slices if needed
    # 3. calculate metrics
    # 4. return plain Python data
    raise NotImplementedError
```

## 8. Template for a runner

```python
from ensembl.genes.annotation_qc.parsers.annotation import parse_annotation
from ensembl.genes.annotation_qc.metrics.feature_wide import compute_feature_wide_metrics
from ensembl.genes.annotation_qc.reports.feature_wide import write_feature_wide_report


def run_feature_wide(annotation_path: str, output_path: str) -> None:
    annotation_df = parse_annotation(annotation_path)
    metrics = compute_feature_wide_metrics(annotation_df)
    write_feature_wide_report(metrics, output_path)
```


## 12. Definition of done

The refactor is complete when:

* the metric module has no file I/O
* the runner owns workflow orchestration
* the parser owns normalization
* the report layer owns formatting
* tests cover the metric logic
* the documentation says exactly where new code should go
* a new starter can add a feature without guessing folder responsibilities

## 13. One-line summary for contributors

**Put parsing in `parsers/`, calculations in `metrics/`, orchestration in `runners/`, output in `reports/`, and configuration in `config/`.**

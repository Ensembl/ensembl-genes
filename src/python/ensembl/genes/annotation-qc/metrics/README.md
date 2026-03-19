# run_agat.py

## Overview

`run_agat.py` is a Python wrapper script used in the Ensembl gene annotation QC pipeline to compute annotation statistics using **AGAT (Another Gff Analysis Toolkit)**.

It processes genome annotation files (typically GFF3/GTF) and generates summary metrics describing gene models, transcripts, exons, CDS features, and other structural properties. These metrics are commonly used to assess annotation quality and completeness.

---

## Features

- Runs AGAT-based analysis on annotation files
- Computes structural metrics such as:
  - Gene counts
  - Transcript counts
  - Exon/CDS statistics
  - Feature lengths and distributions
- Supports configurable feature levels (e.g. gene, transcript, exon)
- Produces summary statistics in text format
- Designed for integration into Ensembl genebuild pipelines

---

## Requirements

### Software

- Python 3.x
- AGAT toolkit singularity

## Usage

```bash
python run_agat.py \
    --gff3 <annotation.gff3> \
    --outdir <output_dir> \
    --feature_levels <fetaure_levels.yaml> OPTIONAL \
    --agat_path <path_to_agat> OPTIONAL
```

 levels - path to feature levels config
 agat_path - path to agat singularity

## Output
Raw agat results in the output folder plus a csv.
 

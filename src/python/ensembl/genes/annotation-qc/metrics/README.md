# run_agat.py

## Overview

`features.py` is a Python wrapper script used in the Ensembl gene annotation QC pipeline to compute annotation statistics using **AGAT (Another Gff Analysis Toolkit)**.

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
- Produces summary statistics in text  and csv format

---

## Requirements

### Software

- Python 3.x
- AGAT toolkit singularity

## Usage

```bash
python run_agat.py \
    --input <annotation.gff3>       # Path to your input GFF3 file
    --output <output_dir>           # Directory where results will be saved
    --feature_levels <feature_levels.yaml>     # (Optional) Specify feature levels to extract
    --agat_path <path_to_agat>      # (Optional) Path to AGAT singularity, if not provided default container will be used
```

## Output
Raw agat results in the output folder plus a csv.

lll
## References
```
Dainat J. 2022. Another Gtf/Gff Analysis Toolkit (AGAT): Resolve interoperability issues and accomplish more with your annotations. Plant and Animal Genome XXIX Conference. https://github.com/NBISweden/AGAT.
```

 

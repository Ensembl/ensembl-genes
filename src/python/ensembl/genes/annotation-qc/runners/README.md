# This is the README for the temporary annotation comparison workflow

Annotations need to be on the same coordinate system aka same GCA!
## Genome Annotation Comparison Tool

A Python-based tool for comparing two genome annotations in **GFF3** or **GTF** format.  
It computes summary statistics, evaluates genomic overlap, and generates visual comparisons between annotations.


### Features

- Runs `gffcompare` for structural comparison
- Parses both GFF3 and GTF formats
- Computes annotation statistics:
  - Gene, transcript, exon, and CDS counts
  - Gene length distributions
  - Exons and transcripts per gene
- Identifies:
  - Shared genes (based on genomic coordinates)
  - Genes unique to each annotation
- Generates publication-ready plots:
  - Feature counts
  - Gene length distributions
  - Average metrics
  - Chromosome-level gene distribution


### Prerequisites

Before using this tool, ensure the following dependencies are installed:

#### 1. `gffcompare`

This tool relies on **gffcompare** for annotation comparison.

- Install from: https://github.com/gpertea/gffcompare
- Ensure the binary is accessible and provide its path via `--gffcompare_bin`

#### 2. Reference Annotation

You must provide a **reference annotation file** (`--gff_ref`).

- This does **not need to be a true biological reference**
- It can be **any annotation you want to compare against**
- The tool treats it as a baseline for comparison

#### 3. Python Dependencies

Install required Python packages from requirements.txt (repo base).


### Usage

```bash
python compare_annotations.py \
  --gff_ref reference.gff3 \
  --gff_test test.gff3 \
  --ref_name Reference \
  --test_name Test \
  --gffcompare_bin /path/to/gffcompare \
  --output_path results/
 ```

---
## BUSCO Missing Annotation Scanner (Advanced)

A tool for investigating **missing BUSCO genes** between two annotation sets.  
It identifies BUSCOs that are *complete in a reference* but *absent in a test annotation*, and then scans multiple annotation layers to determine whether those genes are actually present but missed, fragmented, or structurally different.


### Overview

This tool answers the question:

> *Are BUSCOs truly missing, or just not correctly annotated?*

For each BUSCO that is **complete in the reference but missing in the test**, the tool:

- Maps BUSCO IDs to reference transcripts and coordinates
- Scans one or more annotation files (GTF/GFF)
- Evaluates structural similarity to the reference
- Reports presence, coverage, exon structure, and multi-hit cases


### Features

- BUSCO-aware comparison between reference and test
- Supports **GFF3 and GTF** formats
- Multi-layer annotation analysis (e.g. transcriptomic, protein, orthology)
- Structural comparison based on:
  - Exon overlap
  - Intron count
  - CDS-only mode (optional)
- Classification of annotation presence:
  - `PRESENT_FULL`
  - `PRESENT_LONGER`
  - `PRESENT_SHORT`
  - `PRESENT_SHORT_INTRONS`
  - `ABSENT`
- Multi-hit detection (multiple candidate transcripts)
- Automated error categorisation (`status_check`)


### Prerequisites

#### 1. BUSCO Output Files

You need BUSCO result tables for:

- Reference annotation (`--busco_ref`)
- Test annotation (`--busco_test`)

These should be the standard BUSCO `full_table.tsv` files.


#### 2. Reference Annotation

A reference GFF/GTF file is required
- This provides the **ground truth coordinates and exon structure**
- It does **not need to be a perfect reference**
- It can be any annotation you want to compare against

The script is currently optimised for anno gtfs.

#### 3. Python Dependencies

Install required Python packages from requirements.txt (repo base).

### Usage

```bash
python busco_missing_scanner.py \
  --busco_ref ref_busco.tsv \
  --busco_test test_busco.tsv \
  --ref_annotation ref.gff3 \
  --test_annotation test.gff3
```

#### Optional annotation layers
You should include additional annotation sources to diagnose where genes are recovered

| Argument | Description                              |
|----------|------------------------------------------|
| `--transcriptomic` | Transciptomic annotation gtf             |
| `--protein` | Protein alignment (BUSCO) annotation gtf |
| `--protein_sel` | Selected protein (BUSCO) models gtf      |
| `--uniprot` | UniProt-based annotation gtf             |
| `--uniprot_sel` | Selected UniProt models gtf              |
| `--miniprot_ortho` | Orthology-based projection gtf           |
| `--miniprot_ortho_sel` | Selected orthology models gtf            |

#### CDS-only Mode

`--cds_only` Restricts comparison to CDS regions only. Removes UTR influence.

### Output

Each row corresponds to a BUSCO missing in the test annotation.

Core columns:

- `BUSCO_ID`
- `Reference_Transcript`
- `Chrom`, `Start`, `End`

Per annotation layer:

- `<layer>_status` — presence classification
- `<layer>_coverage_pct` — exon overlap with reference
- `<layer>_introns` — intron count
- `<layer>_hits` — number of candidate transcripts
- `<layer>_transcripts` — matching transcript IDs

#### Status Definitions

| Status | Meaning |
|--------|--------|
| `PRESENT_FULL` | Exact exon structure match |
| `PRESENT_LONGER` | Match exists but extends beyond reference |
| `PRESENT_SHORT` | Partial match (reduced coverage) |
| `PRESENT_SHORT_INTRONS` | Partial match with intron structure |
| `ABSENT` | No overlapping transcript |


#### Error Classification (`status_check`)

The script assigns high-level diagnostic labels. These mostly work but may need adjusting based on the layers used.

| Label | Interpretation |
|-------|----------------|
| `non_can_hit` | Found in test but not canonical |
| `miniprot_hit` | Recovered via orthology |
| `genebuild_error` | Present in evidence but missing in final build |
| `select_error` | Lost during model selection |
| `short_intron_error` | Fragmented or structurally incomplete |

### Notes
- BUSCO presence is evaluated using exon-level overlap, not just genomic span
- Multi-hit regions are explicitly tracked
- Coordinate consistency across annotations is assumed
- No chromosome name normalization (e.g. chr1 vs 1)
- GFF3 and GTF attribute parsing is handled automatically
- Structural comparison is exon-based; does not validate sequence correctness


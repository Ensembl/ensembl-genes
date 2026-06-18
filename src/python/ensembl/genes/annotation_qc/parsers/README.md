# Annotation Parser (GFF3 / GTF)

A small Python utility for parsing **GFF3** and **GTF** genome annotation files using the `pyranges1` library.  
The script reads the input file and loads it into a **PyRanges object**, allowing efficient genomic interval manipulation.

## Features

- Parses **GFF3** files
- Parses **GTF** files
- Automatically detects file type from extension
- Returns and prints a **PyRanges** object representation of the annotation

## Requirements

- Python 3.8+
- pyranges1

Install dependencies:

```bash
pip install pyranges1
```

## Usage

Run the script and provide the path to an annotation file.

```bash
python parse_annotation.py --file_path annotations.gff3
```

or

```bash
python parse_annotation.py --file_path annotations.gtf
```

## Arguments

| Argument | Description | Required |
|--------|-------------|--------|
| `--file_path` | Path to the input GFF3 or GTF file | Yes |

## Supported File Types

The script determines the parser based on the file extension:

- `.gff3` → parsed using `pr.read_gff3()`
- `.gtf` → parsed using `pr.read_gtf()`

If another file type is provided, the script raises an error.

## Example

```bash
python parse_annotation.py --file_path example_annotations.gff3
```

Example output:

```
Parsing GFF3... (example_annotations.gff3)
+--------------+-----------+-----------+------------+
| Chromosome   | Start     | End       | Feature    |
|--------------|-----------|-----------|------------|
| chr1         | 1000      | 2000      | gene       |
| chr1         | 1500      | 1800      | exon       |
...
```

## Script Overview

The script contains three main functions:

- `parse_gff3()` – reads a GFF3 file using `pyranges1.read_gff3`
- `parse_gtf()` – reads a GTF file using `pyranges1.read_gtf`
- `main()` – detects file type and calls the appropriate parser

## Notes

- The file extension must be `.gff3` or `.gtf`.
- The parsed data is returned as a **PyRanges object**, which can be used for downstream genomic interval analysis.

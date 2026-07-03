# parsers

## Purpose

The `parsers` package is responsible for reading external annotation files and converting them into a standardized internal representation that can be consumed by the rest of the Annotation QC framework.

Parsers isolate all file-format specific logic from the QC implementation. Once data has been parsed, the rest of the framework should no longer need to know whether it originated from a GFF3, GTF, TSV, or another supported format.

This package is the **only place where raw annotation files should be interpreted**.

---

# Responsibilities

Parsers are responsible for:

* Reading annotation files
* Parsing feature attributes
* Normalising identifiers
* Validating basic file structure
* Producing the framework's standard annotation representation
* Handling format-specific differences

Typical examples include:

* GFF3 parser
* GTF parser
* BED parser
* TSV evidence parser

---

# What does NOT belong here?

Do **not** place the following inside this package:

* QC metrics
* PASS / FAIL decisions
* Feature statistics
* Report generation
* CLI logic
* Plot generation

If code computes biological statistics, it belongs in `metrics/`.

---

# Inputs

Parsers operate on external files such as

* GFF3
* GTF
* BED
* FASTA
* TSV
* Other annotation formats

Example

```text
annotation.gff3
```

---

# Outputs

All parsers should return the same internal representation.

For example

```python
annotation_df
```

or another agreed internal object.

Every downstream package should be able to use this output without caring about the original file format.

---

# Design principles

A parser should

* perform parsing only
* avoid QC calculations
* avoid writing output files
* avoid plotting
* avoid command-line interaction

A parser should be deterministic:

```
Input file

↓

Same parsed output

↓

Independent of QC workflow
```

---

# Adding a new parser

When adding support for a new file format

1. Create a new parser module.
2. Read the external format.
3. Convert fields into the standard representation.
4. Validate required fields.
5. Add unit tests.
6. Update this documentation.

---

# Common mistakes

❌ Computing statistics while parsing

❌ Removing information that may be useful later

❌ Returning a different object type from other parsers

❌ Embedding report logic

---

# Testing

Every parser should have dedicated unit tests covering:

* valid input
* malformed files
* missing required fields
* edge cases
* expected output schema

---

# Related packages

* `metrics/`
* `runners/`
* `reports/`

The parser should be reusable by any workflow in the framework.


# Parsers

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
python annotation.py --file_path annotations.gff3
```

or

```bash
python annotation.py --file_path annotations.gtf
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
python annotation.py --file_path example_annotations.gff3
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

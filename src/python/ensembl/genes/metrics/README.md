# ensembl_gff_metrics.py

`ensembl_gff_metrics.py` analyses an Ensembl-style GFF3 annotation file and generates
summary metrics describing the genes and transcripts present in the annotation.

The script is designed to be run as a standalone command-line tool and is typically
used during Ensembl genebuild and annotation quality-control workflows.

It reads a GFF3 file, reconstructs gene and transcript structures, calculates metrics
for different gene types, retrieves assembly metadata from NCBI, and writes the
results to files in the requested format.

---

## What the script does

When run, the script will:

- Parse the input GFF3 file
- Build relationships between genes, transcripts, exons and CDS features
- Classify genes as coding, non-coding or pseudogenes
- Calculate metrics including:
  - numbers of genes and transcripts
  - exon counts
  - total coding sequence length
- Retrieve assembly metadata from NCBI using the provided assembly accession
- Write the calculated metrics to files in the output directory

---

## How to run the script

The script is run from the command line.

```text
usage: ensembl_gff_metrics.py [-h] --gff3 GFF3 --outdir OUTDIR
                              --assembly_accession ASSEMBLY_ACCESSION
                              [--scientific_name SCIENTIFIC_NAME]
                              [--taxon_id TAXON_ID]
                              [--strain STRAIN]
                              [--sex SEX]
                              [--keep_ncbi]
                              [--ncbi_debug_dir NCBI_DEBUG_DIR]
                              [--json]
                              [--sql]
                              [--dbname DBNAME]
                              [--species_id SPECIES_ID]
                              [--sql_filename SQL_FILENAME]
### Required arguments

- `--gff3`  
  Path to the input GFF3 annotation file.

- `--outdir`  
  Directory in which output files will be written.

- `--assembly_accession`  
  NCBI assembly accession used to retrieve assembly metadata.

### Optional arguments

- `--scientific_name`  
  Scientific name of the species. If not provided, it is retrieved from NCBI.

- `--taxon_id`  
  NCBI taxon ID. If not provided, it is retrieved from NCBI.

- `--strain`  
  Strain name associated with the assembly.

- `--sex`  
  Sex associated with the assembly.

- `--keep_ncbi`  
  Keep intermediate files downloaded from NCBI.

- `--ncbi_debug_dir`  
  Directory in which to store NCBI debug and intermediate files.

- `--json`  
  Write metrics in JSON format.

- `--sql`  
  Write metrics in SQL format.

- `--dbname`  
  Database name to associate with SQL output.

- `--species_id`  
  Ensembl species identifier used in SQL output.

- `--sql_filename`  
  Name of the SQL output file.

### Output

All output files are written to the directory specified by `--outdir`.

Depending on the command-line options used, output may include:

- Plain text summary metrics
- JSON-formatted metrics
- SQL statements suitable for loading into a database

Only the output formats explicitly requested are generated.

### Important notes

- The input GFF3 file must follow Ensembl conventions.
- Parent–child relationships are inferred from the `Parent` attribute.
- Features missing required relationships are skipped.
- Assembly metadata is retrieved once and reused during processing.

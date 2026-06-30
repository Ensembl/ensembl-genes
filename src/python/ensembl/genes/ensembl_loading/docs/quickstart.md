# Quickstart

## Running The CLI

Install the package in editable mode from the repository root:

```bash
pip install -e .
gff-loader --help
```

You can also run the CLI without installing the console script:

```bash
python src/python/ensembl/genes/ensembl_loading/gff_cli.py --help
```

Set logging verbosity with `--log-level`:

```bash
gff-loader --log-level DEBUG load-features annotations.gff3 ...
```

Write the same log output to a file with `--log-file`:

```bash
gff-loader --log-file gff-loader.log load-features annotations.gff3 ...
```

## General Use

There are two general loading modes.

`load-features`
: Load a GFF3 or supported GTF into an existing Ensembl core database. The
database, schema, coord_system rows, seq_region rows, and DNA must already
exist. This mode only inserts annotation feature rows.

`create-core`
: Create or reuse a core database, load schema SQL, bootstrap core metadata,
load FASTA sequence into `seq_region` and `dna`, then load GFF3/GTF features.
This mode is useful when the FASTA and annotation file already use matching
seq_region names.

### Loading Features Into An Existing Core

Use `load-features` when you already have an Ensembl core database with the
correct coordinate system and seq_region rows.

```bash
gff-loader load-features annotations.gff3 \
  --db-name mus_musculus_core_000001635_27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password>
```

All MySQL connection details are required for commands that connect to a core:
`--db-host`, `--db-port`, `--db-user`, and `--db-password`.
The documented CLI style uses hyphenated option names, but underscore aliases
such as `--db_port` are also accepted for DB connection options.

The non-DB defaults are:

```text
--source generic
--coord-system-name primary_assembly
```

Use `--coord-system-id` if you know the exact coordinate system ID:

```bash
gff-loader load-features annotations.gff3 \
  --db-name mus_musculus_core_000001635_27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-id 1 \
  --source generic
```

For Ensembl-style GFF3 exports, use the dedicated Ensembl source config:

```bash
gff-loader load-features Homo_sapiens.GRCh38.115.chromosome.20.gff3.gz \
  --db-name homo_sapiens_core_115_38 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name chromosome \
  --source ensembl
```

Use `--coord-system-name primary_assembly` instead if the target core stores the
GFF seqids under that coordinate system. The right value is the coordinate
system containing seq_region names such as `1`, `20`, `X`, `MT`, or scaffold
names from the GFF.

Use `--coord-system-name` and `--coord-system-version` when the target core has
multiple coordinate systems:

```bash
gff-loader load-features annotations.gff3 \
  --db-name example_core \
  --db-host mysql-host \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name scaffold \
  --coord-system-version GCA_000000000.1
```

For anno pipeline GTF files, use the dedicated `anno_gtf` and `ncrna_gtf`
source configs.

Important: `annotation_output/annotation.gtf` is only the main gene set. It
does not contain the Rfam and tRNAscan ncRNA genes. To reproduce the original
anno pipeline core load, load the main annotation GTF first with
`--source anno_gtf`, then load the separate ncRNA GTF outputs with
`--source ncrna_gtf`.

The simplest way to reproduce the full anno pipeline load is the combined
wrapper command. Point it at the anno run output directory that contains
subdirectories such as `annotation_output`, `rfam_output`, `dust_output`, and
`cpg_output`:

```bash
gff-loader load-anno-output /path/to/GCA_000000000.1 \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly
```

The wrapper requires a non-empty
`annotation_output/annotation.gtf`. Optional ncRNA, repeat, and simple-feature
outputs are loaded when their `annotation.gtf` files exist and are non-empty;
missing or empty optional files are skipped.

The usual anno pipeline load sequence is:

```bash
gff-loader load-features /path/to/annotation_output/annotation.gtf \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source anno_gtf

gff-loader load-features /path/to/rfam_output/annotation.gtf \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source ncrna_gtf

gff-loader load-features /path/to/trnascan_output/annotation.gtf \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source ncrna_gtf
```

The Rfam/tRNAscan loads insert genes with source `ensembl` and analysis logic
name `ncrna`, matching the original Perl loader. If either ncRNA output file is
empty or absent, skip that file.

Anno repeat/simple outputs are loaded separately with
`load-single-line-features`, matching the old Perl
`-load_type single_line_feature` mode:

```bash
gff-loader load-single-line-features /path/to/dust_output/annotation.gtf \
  --analysis-name dust \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly
```

Use the matching analysis name for each anno output:

```text
dust_output/annotation.gtf         -> dust
red_output/annotation.gtf          -> repeatdetector
trf_output/annotation.gtf          -> trf
repeatmasker_output/annotation.gtf -> repeatmask_repbase_human
cpg_output/annotation.gtf          -> cpg
eponine_output/annotation.gtf      -> eponine
```

The repeat analyses load into `repeat_consensus` and `repeat_feature`.
`cpg` and `eponine` load into `simple_feature`.

Set `--coord-system-name` to the coordinate system containing the GTF seqids.
For example, if the GTF first column contains names such as `17`, the target
core must already have matching `seq_region.name` values under the chosen
coordinate system.

`load-features` does not create a database, does not load schema SQL, does not
insert species or assembly metadata, and does not load FASTA sequence. It
parses the GFF3/GTF and inserts into these core feature tables:

```text
analysis
gene
transcript
exon
exon_transcript
translation
```

The `analysis` row is reused if `analysis.logic_name` already exists for the
selected source config. Otherwise it is created from the config.

### Creating A Core From FASTA And GFF3/GTF

Use `create-core` when you have a FASTA and GFF3/GTF whose sequence names
already match each other.

```bash
gff-loader create-core annotations.gff3 genome.fna \
  --species-name "Mus musculus" \
  --assembly-accession GCF_000001635.27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --source generic
```

By default, `create-core` loads the bundled schema:

```text
src/python/ensembl/genes/ensembl_loading/config/core_schema.sql
```

Use `--schema-sql-path /path/to/schema.sql` to override it. Use
`--schema-sql-path ""` to skip schema loading explicitly.

The database name is derived from species and assembly accession. Generic
sources keep the compact historical format:

```text
Mus musculus + GCF_000001635.27 -> mus_musculus_core_000001635_27
```

When `--source refseq` is used, the derived name follows the RefSeq production
style:

```text
Scientific name + GCF_037462849.1 -> scientific_name_gcf037462849v1_core_114_1
```

`create-core` performs these operations in one transaction:

1. Connect to the MySQL server without selecting a database.
2. `CREATE DATABASE IF NOT EXISTS <derived_core_db_name>`.
3. `USE <derived_core_db_name>`.
4. Load schema SQL from the bundled `config/core_schema.sql`, unless
   `--schema-sql-path` overrides or disables it.
5. Insert one `coord_system` row:
   `primary_assembly`, empty version, rank `1`,
   `default_version,sequence_level`.
6. Insert core `meta` rows for species name, assembly accession, assembly name,
   and genebuild/transcriptbuild/exonbuild levels.
7. Insert an `analysis` row from the selected source config.
8. Load FASTA sequences into `seq_region`, `dna`, and `seq_region_attrib`.
9. Parse and load GFF3 features.
10. Commit, or roll back the whole transaction if any step fails.

The optional positional `assembly_report` argument is retained for compatibility
with RefSeq-style calling code, but the core loader itself does not read the
assembly report. The FASTA and GFF3 must already contain the sequence names that
will be used as `seq_region.name`.

### Input Requirements

The loader expects GFF3-like rows with nine tab-separated columns. Comment rows
starting with `#` are skipped. Plain text and `.gz` GFF3 inputs are supported.

The parser expects standard `ID` and `Parent` attributes. Attributes are parsed
by splitting on semicolon and then on the first equals sign. This is simple and
matches the original loading assumptions; it does not URL-decode values and does
not split comma-separated multi-parent relationships.

Important input assumptions:

1. GFF seqids must match `seq_region.name` in the target core for
   `load-features`.
2. FASTA headers must match GFF seqids for `create-core`.
3. A feature with `strand` equal to `+` is loaded as `1`; any other strand value
   is currently loaded as `-1`.
4. Multi-parent exons or CDS rows should be pre-normalized before loading.
5. Features outside the configured gene/transcript/exon/CDS model are ignored.

### Why Coordinate System Options Are Needed

Ensembl stores a feature location as `seq_region_id`, start, end, and strand.
The loader therefore has to convert each GFF seqid into a `seq_region_id`.

For `load-features`, the `seq_region` rows already exist. The loader resolves
them using either:

```text
--coord-system-id
```

or:

```text
--coord-system-name
--coord-system-version
```

`--coord-system-id`
: Uses an exact `coord_system_id`. This is the safest option when the database
has several coordinate systems with similar names or repeated seq_region names.

`--coord-system-name`
: Looks up the first matching coordinate system by name, ordered by rank. The
default is `primary_assembly`.

`--coord-system-version`
: Adds a version filter to the name lookup.

If any parsed gene or transcript refers to a seqid that cannot be found in the
resolved coordinate system, loading fails before inserts are committed. The
error reports the first missing names and the total count if there are more
than ten.

`create-core` does not expose these options because it creates a new
`primary_assembly` coordinate system and inserts seq_regions from the FASTA.

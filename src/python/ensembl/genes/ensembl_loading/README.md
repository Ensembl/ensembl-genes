# Generic GFF/GTF Loader And RefSeq Import

This module provides `gff-loader`, a single CLI for loading GFF3/GTF annotation
into Ensembl core databases. It supports generic GFF/GTF loading, anno pipeline
GTF outputs, and RefSeq discovery/download/conversion/loading workflows.

The database loader is generic: once source-specific rules are captured in a
`GffSourceConfig`, the same core loading path is used for RefSeq, Ensembl-style
GFF3, anno GTF, ncRNA GTF, and other compatible sources.

## Start Here

| Goal | Command | Detailed doc |
| --- | --- | --- |
| Load GFF3/GTF into an existing core | `gff-loader load-features` | [docs/quickstart.md](docs/quickstart.md) |
| Create a core from FASTA plus GFF3/GTF | `gff-loader create-core` | [docs/quickstart.md](docs/quickstart.md) |
| Load a full anno pipeline output directory | `gff-loader load-anno-output` | [docs/quickstart.md](docs/quickstart.md) |
| Load repeat/simple feature GTF outputs | `gff-loader load-single-line-features` | [docs/source-configs.md](docs/source-configs.md) |
| List, download, convert, or load RefSeq | `gff-loader refseq ...` | [docs/refseq-import.md](docs/refseq-import.md) |
| Add support for another annotation source | edit `gff_source_config.py` | [docs/source-configs.md](docs/source-configs.md) |
| Understand parser, reconciliation, phases, or DB inserts | Python internals | [docs/loader-internals.md](docs/loader-internals.md) |
| Call the loader from Python | reusable functions | [docs/python-api-and-boundaries.md](docs/python-api-and-boundaries.md) |

## Install And Help

Install from the repository root:

```bash
pip install -e .
gff-loader --help
```

You can also run the CLI module directly:

```bash
python src/python/ensembl/genes/ensembl_loading/gff_cli.py --help
```

## Common Workflows

Load annotation into an existing core database:

```bash
gff-loader load-features annotations.gff3 \
  --db-name example_core \
  --db-host mysql-host \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source generic
```

Create a new core from matching FASTA and GFF3/GTF files:

```bash
gff-loader create-core annotations.gff3 genome.fna \
  --species-name "Mus musculus" \
  --assembly-accession GCF_000001635.27 \
  --db-host mysql-host \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --source generic
```

Download, convert, and load a RefSeq assembly into a core database:

```bash
gff-loader refseq run \
  --base-dir refseq_data \
  --assembly-acc GCF_000001635.27 \
  --load-core \
  --species-name "Mus musculus" \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 4532 \
  --db-user <write-user> \
  --db-password <password>
```

Load a complete anno pipeline output directory:

```bash
gff-loader load-anno-output /path/to/anno_run_output \
  --db-name example_core \
  --db-host mysql-host \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly
```

## Common Gotchas

Top-level options must appear before the command group. For example,
`--log-file` belongs to `gff-loader`, not to `refseq run`:

```bash
# Correct
gff-loader --log-file gff-loader.log refseq run --assembly-acc GCF_000001635.27 ...

# Incorrect
gff-loader refseq run --assembly-acc GCF_000001635.27 ... --log-file gff-loader.log
```

For `load-features`, the target core must already contain matching
`coord_system`, `seq_region`, and `dna` rows. Use `create-core` when the loader
should create the database, load schema SQL, load FASTA sequence, and then load
annotation features.

For `create-core`, FASTA headers and GFF seqids must already match. RefSeq files
usually need the RefSeq conversion workflow first because NCBI accessions such
as `NC_000001.11` need to be converted to Ensembl-style seq_region names such
as `1`.

Command-line passwords can appear in shell history and process listings. Prefer
a safer password mechanism when one is available in your environment.

## Module Layout

Generic loading modules:

```text
gff_cli.py
gff_core_loader.py
gff_models.py
gff_repeat_loader.py
gff_source_config.py
```

RefSeq-specific modules:

```text
refseq_ncbi.py
refseq_conversion.py
refseq_models.py
refseq_constants.py
```

See the focused docs for details:

- [docs/quickstart.md](docs/quickstart.md): install, command usage, generic loading, core creation, and coordinate system options.
- [docs/refseq-import.md](docs/refseq-import.md): RefSeq listing, download, conversion, end-to-end load, and RefSeq-specific feature rules.
- [docs/source-configs.md](docs/source-configs.md): source config fields, built-in configs, repeat/simple feature loading, and adding another source.
- [docs/loader-internals.md](docs/loader-internals.md): parser flow, reconciliation, biotype handling, exon phases, translations, inserts, and quality checks.
- [docs/python-api-and-boundaries.md](docs/python-api-and-boundaries.md): reusable Python functions and current boundaries.

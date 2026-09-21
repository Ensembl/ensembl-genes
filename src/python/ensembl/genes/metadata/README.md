# scripts/metadata

## About

This script was used to produce the metadata sql patch for the given database, <database_name>.sql.

It will also produce a log file with any CRITICAL all WARNING messages, <database_name>_metadata.log, you should check this before applying the patch.

## Running core_meta_data.py

**python core_meta_data.py -h**

usage: core_meta_data.py [-h] [-o OUTPUT_DIR] -d DB_NAME -s HOST -p PORT -t TEAM

Prepare SQL updates for core dbs

optional arguments:
  -h, --help            show this help message and exit

  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
     Path where the output and temp files will write to. Uses current dir by default

  -d DB_NAME, --db_name DB_NAME
     Database name

  -s HOST, --host HOST  Host server

  -p PORT, --port PORT  Host server port

  -v, --verbose         Enable verbose output (check that all required keys are not NULL/ empty)

  -t TEAM, --team TEAM  Team responsible for the database

**NOTE**

This script doesn't currently deal with collections, species_ids are hardcoded at **LN:207**. This should be updated, but for now, if you are dealing with a collection db, please hardcode the `species_id`.

## Running truth_checker.py

Compares the output from `core_meta_data.py` with the registry-backed `core_meta_reg.py` for the same core database. It runs both scripts, loads their generated JSON metadata, and fails if the metadata values differ. The `core_meta_data.py` script uses the registry to fetch the latest metadata, while `core_meta_reg.py` will use NCBI and ENA.

Use this when checking that changes in `core_meta_reg.py` still match the existing metadata output.

**Basic usage:**
```bash
python truth_checker.py \
  -d <database_name> \
  -s <core_host> \
  -p <core_port> \
  -t <team>
```

For collection databases, pass the production name:
```bash
python truth_checker.py \
  -d <database_name> \
  -s <core_host> \
  -p <core_port> \
  -t <team> \
  -n <production_name>
```

The registry connection defaults to:

- `--registry_host`: `$GBS1`
- `--registry_port`: `$GBP1`
- `--registry_user`: `ensro`
- `--registry_db`: `gb_assembly_metadata`

You can override these, and pass registry-specific selection options when needed:
```bash
python truth_checker.py \
  -d <database_name> \
  -s <core_host> \
  -p <core_port> \
  -t <team> \
  --assembly_accession <GCA_or_GCF_accession> \
  --genebuilder <genebuilder>
```

Useful debugging options:

- `--work_dir <path>` writes generated SQL and JSON files to a chosen directory.
- `--keep_outputs` keeps temporary output files when `--work_dir` is not provided.
- `--report_json <path>` writes a structured comparison report.
- `--verbose_scripts` passes `-v` through to both metadata scripts.

Successful runs print `TRUTH_CHECK_OK`. Differences print `TRUTH_CHECK_FAILED` with a JSON report showing missing, extra, and mismatched metadata keys.

## Running beta_patcher.py

Generates validation and patch SQL for beta metadata fixes in:

1. the production metadata database (`ensembl_genome_metadata`)
2. the matching beta core database `meta` table

The script takes a CSV of requested changes, resolves each `genome_uuid` to a
`production_name`, discovers the matching core database on ST6/ST5, and writes
separate `validate_*.sql` and `patch_*.sql` files.

**Setup:**

Install `ensembl-genes` in the active environment so the `beta_patcher`
console command is available:

```bash
# From a local checkout
cd ensembl-genes
pip install .
```

For development work, use an editable install instead:

```bash
cd ensembl-genes
pip install -e .
```

The metadata API package is optional, but preferred when available:

```bash
git clone https://github.com/Ensembl/ensembl-metadata-api.git
cd ensembl-metadata-api
pip install -e .
```

You also require connection details for the production metadata and taxonomy
databases. The metadata URI is used by the API path and by the direct SQLAlchemy
fallback.

```bash
export METADATA_URI="mysql+pymysql://user:pass@host:port/ensembl_genome_metadata"
export TAXONOMY_URI="mysql+pymysql://user:pass@host:port/ncbi_taxonomy"
```

If the metadata API is installed but incompatible with the current schema, force
direct SQLAlchemy access:

```bash
export BETA_PATCHER_FORCE_DIRECT=1
```

For metadata-only patch generation without metadata DB/core discovery, use
offline mode. This generates metadata SQL and skips core patches because the
core database name cannot be resolved.

```bash
export BETA_PATCHER_OFFLINE=1
```

**Usage:**

Assume `ensembl-genes` is installed in the active environment and use the
`beta_patcher` console command.

```bash
# Basic usage
beta_patcher \
  patches.csv \
  --jira-ticket EBD-1111 \
  --output-dir ./patches/

# Use explicit metadata/taxonomy URIs instead of environment variables
beta_patcher \
  patches.csv \
  --jira-ticket EBD-1111 \
  --metadata-uri "mysql+pymysql://user:pass@host:port/ensembl_genome_metadata" \
  --taxonomy-uri "mysql+pymysql://user:pass@host:port/ncbi_taxonomy" \
  --output-dir ./patches/

# With team filter (only applies patches where all affected genomes belong to specified team)
beta_patcher \
  patches.csv \
  --jira-ticket EBD-1111 \
  --team-filter Genebuild

# Override beta core server discovery
beta_patcher \
  patches.csv \
  --jira-ticket EBD-1111 \
  --st6-uri "mysql+pymysql://ensro:@mysql-ens-sta-6:4695/" \
  --st5-uri "mysql+pymysql://ensro:@mysql-ens-sta-5:4684/"
```

When working from a repo checkout without installing the package, use the module
entry point as a development fallback:

```bash
PYTHONPATH=src/python python3 -m ensembl.genes.metadata.beta_patcher \
  patches.csv \
  --jira-ticket EBD-1111 \
  --output-dir ./patches/
```

**CSV format:**

Required columns:

- `genome_uuid`
- `meta_key`
- `desired_meta_value`

Optional columns:

- `dataset_type` - defaults to `genebuild`
- `species_id` - defaults to `1`
- `table_location` - defaults to `dataset_attribute`

Valid `table_location` values are:

- `dataset_attribute` - updates metadata DB `dataset_attribute.value` and the
  core DB `meta.meta_value`
- `genome` - updates a column in the metadata DB `genome` table and the core DB
  `meta.meta_value`
- `organism` - updates a column in the metadata DB `organism` table and the core
  DB `meta.meta_value`
- `assembly` - updates a column in the metadata DB `assembly` table and the core
  DB `meta.meta_value`

Use `\N` in `desired_meta_value` to set SQL `NULL` in the metadata database.
Core DB patches store this as the string `'NULL'`, matching the current core
`meta` table convention.

Wrap CSV values containing commas in double quotes. If `species_id` is reported
as invalid, check for an unquoted comma in `desired_meta_value`.

Example:

```csv
genome_uuid,meta_key,desired_meta_value,dataset_type,species_id,table_location
a7335667-93e7-11ec-a39d-005056b38ce3,assembly.name,GRCh38.p14,genebuild,1,dataset_attribute
a7335667-93e7-11ec-a39d-005056b38ce3,genebuild_version,2024-01,genebuild,1,genome
b8446778-a4f8-22fd-b4ae-116167c49df4,strain,reference,genebuild,1,organism
d0668990-c6h0-44hf-d6cg-338389e6bfh6,genebuild.annotation_source,ensembl,genebuild,1,dataset_attribute
```

There is a template at `src/python/ensembl/genes/metadata/patches_template.csv`.

**Generated files:**

For a Jira ticket such as `EBD-1111`, the script writes:

- `validate_metadata_EBD-1111.sql`
- `patch_metadata_EBD-1111.sql`
- `validate_core_EBD-1111_<server>.sql`
- `patch_core_EBD-1111_<server>.sql`
- `patch_csv_<timestamp>.log`

Run the `validate_*.sql` files before applying the corresponding `patch_*.sql`
files.

Alternatively, for large patches it is advisable to make a copy of the metadata
DB and confirm the patches have the intended effect.

**Core database discovery:**

By default, the script searches ST6 first and ST5 second using:

- `mysql+pymysql://ensro:@mysql-ens-sta-6:4695/`
- `mysql+pymysql://ensro:@mysql-ens-sta-5:4684/`

The chosen core database must match `{production_name}%_core_%`. Exact matches
starting `{production_name}_core_` are preferred. If ST5/ST6 discovery is not
used, `--core-suffix` is appended to `production_name`; the default suffix is
`_core_114_1`.

**Safety checks:**

- Existing metadata values are checked. Matching values are written as
  `SKIPPED (no change)`.
- `dataset_attribute` patches update existing `dataset_attribute_id` rows when
  present; otherwise they insert the missing attribute.
- `organism` and `assembly` patches can affect multiple genomes that share the
  same organism or assembly. The output SQL includes affected genomes and teams.
- `--team-filter` skips shared `organism` or `assembly` patches unless all
  affected genomes belong to the requested team.
- The log reports whether any requested fields require THOAS updates.

### Finding genome_uuid for organism/assembly patches

When patching `organism` or `assembly` tables, you need to provide a genome_uuid. Use these queries to find genome UUIDs:

**Find all genomes for a specific assembly (by accession):**
```sql
SELECT DISTINCT
    genome.genome_uuid,
    genome.production_name,
    assembly.accession,
    assembly.name AS assembly_name,
    (SELECT da.value
     FROM genome_dataset gd
     JOIN dataset d ON gd.dataset_id = d.dataset_id AND d.name = 'genebuild'
     JOIN dataset_attribute da ON d.dataset_id = da.dataset_id
     JOIN attribute a ON da.attribute_id = a.attribute_id AND a.name = 'genebuild.team_responsible'
     WHERE gd.genome_id = genome.genome_id
     LIMIT 1) AS team_responsible
FROM genome
JOIN assembly ON genome.assembly_id = assembly.assembly_id
WHERE assembly.accession = 'GCA_000001405.14'
ORDER BY team_responsible, genome.production_name;
```

**Find all genomes for a specific organism (by biosample_id):**
```sql
SELECT DISTINCT
    genome.genome_uuid,
    genome.production_name,
    organism.biosample_id,
    organism.scientific_name,
    organism.strain,
    (SELECT da.value
     FROM genome_dataset gd
     JOIN dataset d ON gd.dataset_id = d.dataset_id AND d.name = 'genebuild'
     JOIN dataset_attribute da ON d.dataset_id = da.dataset_id
     JOIN attribute a ON da.attribute_id = a.attribute_id AND a.name = 'genebuild.team_responsible'
     WHERE gd.genome_id = genome.genome_id
     LIMIT 1) AS team_responsible
FROM genome
JOIN organism ON genome.organism_id = organism.organism_id
WHERE organism.biosample_id = 'SAMN04851098'
ORDER BY team_responsible, genome.production_name;
```

**Find genomes by organism strain:**
```sql
SELECT DISTINCT
    genome.genome_uuid,
    genome.production_name,
    organism.scientific_name,
    organism.strain,
    (SELECT da.value
     FROM genome_dataset gd
     JOIN dataset d ON gd.dataset_id = d.dataset_id AND d.name = 'genebuild'
     JOIN dataset_attribute da ON d.dataset_id = da.dataset_id
     JOIN attribute a ON da.attribute_id = a.attribute_id AND a.name = 'genebuild.team_responsible'
     WHERE gd.genome_id = genome.genome_id
     LIMIT 1) AS team_responsible
FROM genome
JOIN organism ON genome.organism_id = organism.organism_id
WHERE organism.scientific_name = 'Homo sapiens'
ORDER BY team_responsible, genome.production_name;
```

Pick any one of the returned genome_uuid values to use in your CSV. The script will automatically detect and warn about all other genomes sharing that organism/assembly.

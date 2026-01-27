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

## Running beta_patcher.py

Generates SQL patches for beta metadata fixes in both production metadata DB and core DBs.

**Setup:**
```bash
pip install -r requirements.txt
git clone https://github.com/Ensembl/ensembl-metadata-api.git
cd ensembl-metadata-api && pip install -r requirements.txt && pip install -e .
```

You also require two environemtal variables to enable interaction with the production metadata db. Substitute in the valid connection details for these database servers.
```bash
export METADATA_URI="mysql+pymysql://user:pass@host:port/ensembl_genome_metadata"
export TAXONOMY_URI="mysql+pymysql://user:pass@host:port/ncbi_taxonomy"
```

**Usage:**
```bash
# Basic usage
python beta_patcher.py patches.csv --jira-ticket EBD-1111 --output-dir ./patches/

# With team filter (only applies patches where all affected genomes belong to specified team)
python beta_patcher.py patches.csv --jira-ticket EBD-1111 --team-filter Genebuild
```

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
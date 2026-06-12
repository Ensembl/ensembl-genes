# Bioproject Tracking

## About

**Bioproject Tracking** is a Python script that fetches assembly accessions from the NCBI BioProject, a specified Taxonomy ID, or a named project portal (e.g. HPRC) and correlates them with Ensembl metadata. The script writes a summary report including annotation data, optional FTP link validation, and classification at a specified taxonomic rank.

This module sits upstream of the YAML generation pipeline in `ensembl.genes.projects`. Its output (accession lists and reports) can be used to prepare the input files consumed by `generate_project_yaml.py`.

## Modules

| Module | Purpose |
|---|---|
| `bioproject_tracking.py` | Main CLI script. Discovers assemblies via NCBI or project portals, looks them up in Ensembl metadata, optionally adds FTP and taxonomy info, and writes a TSV report. |
| `live_tracking.py` | Auxiliary tracking utilities (pre-existing on `main`). |

## Prerequisites

1. **Python 3.11+** (matching the repository-wide requirement in `pyproject.toml`).
2. **MySQL** connectivity to the Ensembl metadata database (read-only access via `ensro`).
3. A valid `bioproject_tracking_config.json` in the working directory, containing:
   - MySQL connection details for metadata and pre-release servers.
   - Base URLs for the NCBI Datasets API.
4. **NCBI Datasets CLI** (`datasets`) is optional but recommended — it retrieves all assembly versions, whereas the API fallback only returns the latest.
5. Dependencies are declared in the repository-level `pyproject.toml`. Install with:
   ```bash
   pip install -e .
   ```

## Inputs

- **`--bioproject_id`**: An NCBI BioProject ID (e.g. `PRJEB40665`).
- **`--taxon_id`**: A taxonomy ID (e.g. `9606`).
- **`--project_name`**: A named project portal for project-specific discovery. Currently supported: `hprc`.
- Exactly one of the above must be provided; they are mutually exclusive.

## Outputs

- **TSV report file** (default `./report_file.tsv`): One row per assembly/annotation match, with columns for accession, genome UUID, database name, and optionally taxonomy rank, FTP link, and FTP status.

## Configuration

The script reads `bioproject_tracking_config.json` from the current working directory at import time. This file defines:

- `server_details`: MySQL connection details for `meta` (beta metadata DB), `pre-release` (GB1 server), and others.
- `urls.datasets`: Base URLs for NCBI Datasets API endpoints (bioproject and taxon).

**Important**: Run the script from a directory containing `bioproject_tracking_config.json`, or update the config loading path. Database passwords are typically omitted when running from the EBI cluster with standard `.my.cnf` credentials.

## Usage

```bash
python -m ensembl.genes.tracking.bioproject_tracking -h
```

### Basic Examples

```bash
# Track assemblies from a BioProject
python -m ensembl.genes.tracking.bioproject_tracking \
  --bioproject_id PRJEB40665 \
  --report_file dtol_report.tsv

# Track assemblies for a taxonomy ID, haploid only, with classification
python -m ensembl.genes.tracking.bioproject_tracking \
  --bioproject_id PRJEB40665 \
  --haploid \
  --classification \
  --rank class \
  --report_file dtol_rank_report.tsv

# Include FTP link validation
python -m ensembl.genes.tracking.bioproject_tracking \
  --bioproject_id PRJEB40665 \
  --ftp \
  --report_file dtol_ftp_report.tsv

# Include pre-release databases for missing accessions
python -m ensembl.genes.tracking.bioproject_tracking \
  --bioproject_id PRJEB40665 \
  --pre_release \
  --report_file dtol_full_report.tsv
```

### HPRC Project-Specific Discovery

```bash
python -m ensembl.genes.tracking.bioproject_tracking \
  --project_name hprc \
  --report_file hprc_report.tsv
```

This will:
- Download the HPRC catalog JSON from the HPRC data explorer GitHub repository.
- Filter to release 2 assemblies with valid GCA accessions.
- De-duplicate and feed those GCAs into the Ensembl metadata lookup and report flow.

### CLI Options

| Flag | Description |
|---|---|
| `--bioproject_id` | NCBI BioProject ID |
| `--taxon_id` | Taxonomy ID |
| `--project_name` | Project-specific discovery (e.g. `hprc`). Mutually exclusive with the above. |
| `--haploid` | Restrict to haploid assemblies only |
| `--report_file` | Output TSV path (default: `./report_file.tsv`) |
| `--classification` | Include taxonomic classification breakdown |
| `--rank` | Taxonomic rank to retrieve (default: `order`) |
| `--ftp` | Validate and include FTP links in the report |
| `--pre_release` | Include pre-release database matches for missing accessions |

## Assumptions and Limitations

- The script loads `bioproject_tracking_config.json` from the current working directory at module import time. This means import will fail if the config file is not present.
- The `--classification` flag makes one NCBI API call per accession, which can be slow for large result sets.
- The NCBI Datasets API fallback (when the CLI is unavailable) only returns the latest assembly version per accession.
- FTP link validation (`--ftp`) makes HTTP HEAD/GET requests to check file existence, which can be slow for large datasets.
- Pre-release database lookup (`--pre_release`) uses `information_schema` pattern matching, which may return unexpected matches if naming conventions change.

## Relationship to `ensembl.genes.projects`

This tracking module and the `projects` YAML generation module serve different purposes:

- **Tracking** (`bioproject_tracking.py`): Discovers and reports on assemblies — answers "what assemblies exist and which have Ensembl annotations?"
- **Projects** (`generate_project_yaml.py`): Takes a list of genome identifiers and produces publishable YAML — answers "what should appear on the project page?"

A typical workflow is:
1. Run `bioproject_tracking.py` to discover assemblies and generate a tracking report.
2. Extract the genome UUIDs from the report to create an input file.
3. Run `generate_project_yaml.py` with that input file to produce the final YAML.

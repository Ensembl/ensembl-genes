# Projects Pages

Code for updating the Ensembl projects pages.

This repository handles the extraction of metadata from Ensembl tracking databases to generate publishable `species.yaml` files used by the various Ensembl project portals.

## Quick Start

The supported entrypoint for generating YAML is `generate_project_yaml.py`.

From the repo root:

```bash
export PYTHONPATH=$PWD/src/python

python -m ensembl.genes.projects.generate_project_yaml ./cbp_guuids.txt \
  --project cbp \
  --output cbp_species.yaml \
  --audit-file cbp_audit.tsv
```

## Environment Requirements

- **Python Version**: Python 3.9+ is recommended. (Note: Ensure type hint compatibility if using older versions like 3.8).
- **PYTHONPATH**: Must include the `src/python` directory of this repository so modules can be resolved (`export PYTHONPATH=$PWD/src/python`).
- **Required Packages**: 
  - `PyMySQL` (for DB connections)
  - `PyYAML` (for YAML generation)
  - `requests` (for FTP checks)
  - `xmltodict` (for NCBI taxonomy parsing)
- **Cluster Notes**: Database passwords may be omitted if you are running this from a cluster that has standard `.my.cnf` or generic Ensembl read-only credentials configured.

## Configuration

Configuration is managed via `config.py` and `server_config.json`:
- **`server_config.json`**: Lives in `src/python/ensembl/genes/projects/` and contains the DB connection details.
- **Metadata DB**: Configured under the `meta_beta` key in `server_config.json`. This tracks fully released genomes via their Genome UUIDs.
- **GB Registry DB**: Configured under the `gb1` key. This tracks pre-release genomes.
- **Updating Config**: Update the hosts or ports in `server_config.json` if servers are migrated.

## Input File Format

The script requires a positional input file (`input_guuids.txt`).
- **Format**: One identifier per line.
- **Identifiers**: Primarily, these should be Genome UUIDs (GUUIDs), which are standard UUID strings (e.g., `123e4567-e89b-12d3-a456-426614174000`).
- **Pre-release**: Core DB names may still be accepted for explicit pre-release handling, but this fallback mechanism is deprecated in favor of automated discovery.
- **Discovery**: Pre-release genomes that don't have GUUIDs yet are automatically discovered from the GB registry if the project is configured for scoping (e.g., `bioproject_scoping` in `config.py`).

## Supported Projects

The `--project` flag maps to configurations defined in `config.py`. Supported projects include:
- `cbp`
- `vgp`
- `dtol` (Darwin Tree of Life)
- `bge`
- `asg`
- `erga`
- `hprc`
- `mouse_genomes` (or `mouse`)
- `aquafaang`

## Output

The script generates two main outputs:
1. **`species.yaml`**: The rendered YAML file containing all genomes that successfully passed validation and FTP existence checks.
2. **Audit TSV (`--audit-file`)**: A highly recommended log file that records the decision for every candidate genome. 
   - Decisions include: `included_released`, `included_prerelease`, `excluded`, `excluded_duplicate`, or `kept_duplicate`.

## Publishability Rules

A genome is only included in the final YAML if it has valid FTP assets. A GUUID alone is not enough.
- **Released FTP**: The script checks the released FTP path.
- **Pre-release Fallback**: If the released FTP is not found, the script tries the pre-release FTP fallback.
- **Exclusion**: If neither exists, the genome is excluded from the YAML.
- **Dynamic Date Resolution**: The metadata database date is used as a hint, but the code dynamically lists the FTP directory to resolve the actual date if it differs.
- **Annotation Source**: The `annotation_source` from the metadata controls whether the FTP path points to `ensembl`, `braker`, or other sources.

## Images and Icons

The `icons.txt` file controls the image mappings for standard project pages.
- The script fetches the taxonomy lineage from NCBI.
- The most specific matching taxonomy classification in `icons.txt` wins.
- The default icon is `Metazoa.png` (or `Chordates.png` for chordates lacking a more specific mapping).
- **Important**: The species display name must always be the scientific name. It must not be polluted by strain, sample descriptions, or habitat text.

## HPRC and Mouse Differences

Different projects map to different underlying YAML schemas (`standard`, `hprc`, or `mouse`):
- **HPRC (`--project hprc`)**: 
  - Adds HPRC-specific columns such as `assembly_submitter` and `assembly_link`.
  - Determines `parent_of_origin` from the assembly report.
  - Automatically checks for and links to `variants_vep` if available on the FTP.
- **Mouse Genomes (`--project mouse_genomes`)**:
  - Automatically incorporates the `strain` metadata explicitly as a separate display field.
  - Includes `alternate` linkage if it is an alternate haplotype.

## File Roles

Understanding the codebase organization:
- **`generate_project_yaml.py`**: **The supported production CLI entrypoint.** Coordinates finding candidates and rendering them.
- **`yaml_renderer.py`**: Handles rendering `GenomeMetadata` models into validated dicts and performs actual FTP publishability checks.
- **`ftp_client.py`**: Core logic for querying the EBI/Ensembl FTP servers.
- **`gb_tracker.py`**: Client for GB registry pre-release discovery.
- **`metadata_db.py`**: Client for the metadata DB to look up released GUUIDs.
- **`bioproject_tracking.py`**: An auxiliary tool for candidate discovery (to generate the input file), but *not* for final YAML eligibility.
- **`write_yaml.py` / `hprc_write_yaml.py`**: Legacy/deprecated scripts kept for backwards compatibility. Do not use for new runs.

## Troubleshooting

- **`ModuleNotFoundError: No module named 'ensembl'`**: Ensure `PYTHONPATH` is set correctly (`export PYTHONPATH=$PWD/src/python`).
- **Python Type Hint Issues**: The code uses Python 3.9+ type hints. Upgrade your Python version if you see syntax errors around type annotations.
- **DB Connection Errors**: Check `server_config.json` hosts and your VPN/network connection. If on a cluster, ensure read-only MySQL permissions are loaded.
- **Missing FTP Assets**: If a genome is excluded, check the audit TSV. Usually, the files haven't been synchronized to the public EBI FTP yet.
- **Duplicate GUUID/Accession**: The script deduplicates automatically based on accession, species name, and date, keeping the most recent released copy. The excluded ones will show `excluded_duplicate` in the audit TSV.
- **Wrong Image Mapping**: Check that the NCBI taxonomy matches an entry in `icons.txt`. Add a new mapping if the lineage isn't covered.
- **Unexpected Species Names**: Check the core tracking database; the script strictly uses `species.scientific_name` to prevent formatting pollution.

## Example Commands

### CBP
```bash
python -m ensembl.genes.projects.generate_project_yaml ./cbp_guuids.txt \
  --project cbp --output cbp_species.yaml --audit-file cbp_audit.tsv
```

### HPRC
```bash
python -m ensembl.genes.projects.generate_project_yaml ./hprc_guuids.txt \
  --project hprc --output hprc_species.yaml --audit-file hprc_audit.tsv
```

### Mouse Genomes
```bash
python -m ensembl.genes.projects.generate_project_yaml ./mouse_guuids.txt \
  --project mouse_genomes --output mouse_species.yaml --audit-file mouse_audit.tsv
```

### Dry Run / Sanity Validation
To check what will be included without committing the YAML immediately, you can run the pipeline and inspect the audit file:
```bash
python -m ensembl.genes.projects.generate_project_yaml ./cbp_guuids.txt \
  --project cbp --output /dev/null --audit-file cbp_audit.tsv

# See what was dropped because FTP files are missing:
grep excluded cbp_audit.tsv

# See what genomes successfully fell back to pre-release FTP:
grep included_prerelease cbp_audit.tsv
```

### Changelog / Regression Detection
Compare a newly generated YAML against a previous version to detect added, removed, or modified entries:
```bash
python -m ensembl.genes.projects.generate_project_yaml ./cbp_guuids.txt \
  --project cbp \
  --output cbp_species.yaml \
  --changelog old_cbp_species.yaml
```

This prints a human-readable diff summary to stdout, for example:
```
=== CHANGELOG SUMMARY ===

Added (3):
  - GCA_XXXXX (Drosophila melanogaster)

Removed (2):
  - GCA_YYYYY (Bombus terrestris)

Modified (4):
  - GCA_ZZZZZ (Papilio machaon)
    - annotation_gtf:
        OLD: .../2025_09/genes.gtf.gz
        NEW: .../2025_10/genes.gtf.gz
    - image:
        OLD: Lepidoptera.png
        NEW: Arthropods.png

Total: 3 added, 2 removed, 4 modified.
```

To also write a machine-readable TSV changelog:
```bash
python -m ensembl.genes.projects.generate_project_yaml ./cbp_guuids.txt \
  --project cbp \
  --output cbp_species.yaml \
  --changelog old_cbp_species.yaml \
  --changelog-output cbp_changelog.tsv
```

The TSV contains columns: `status`, `accession`, `species`, `field`, `old_value`, `new_value`.

**Comparison details:**
- The comparison key is the `accession` field (or `assembly_accession` for HPRC, or `species`/`assembly` as fallback).
- URL trailing slashes and date format differences (`YYYY-MM` vs `YYYY_MM`) are normalised before comparison.
- Internal audit fields are excluded from the diff.
- The changelog is purely informational — it never alters the generated YAML or fails the run.
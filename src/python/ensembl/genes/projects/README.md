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

Icons are resolved from taxonomy lineage data via `icon_resolver.py`, using multiple data sources for robustness.

**Lineage data sources (tried in priority order):**
1. **Metadata DB lineage** — taxonomy classification pre-fetched from the Ensembl metadata database during data gathering (`GenomeMetadata.taxonomy_lineage`). This is the fastest and most reliable source.
2. **NCBI Entrez lineage** — live lookup by `taxon_id` via the NCBI Entrez API. Results are cached per `taxon_id` during a run.
3. **BUSCO lineage hint** — weak but reliable fallback parsing the `busco_lineage` string (e.g. `"insecta_odb10"` → Arthropods.png, `"mammalia_odb10"` → Mammals.png). Used when both metadata and NCBI sources fail.

The first source that produces a usable lineage is walked from most specific to broadest classification. The first taxonomy name that matches a configured rule determines the icon.

**Rule sources (in priority order):**
- `icons.txt` provides project-specific overrides or additions. Entries here win over built-in defaults for the same taxonomy name.
- Built-in default rules in `icon_resolver.py` cover all available icon files with broad taxonomy anchors (e.g. `Arthropoda -> Arthropods.png`, `Mammalia -> Mammals.png`, `Fungi -> Fungi.png`).

**Most-specific match wins:**
- Lepidoptera beats Arthropoda (→ `Lepidoptera.png`, not `Arthropods.png`)
- Mammalia beats Chordata (→ `Mammals.png`, not `Chordates.png`)
- Aves beats Chordata (→ `Birds.png`)
- Actinopterygii / Chondrichthyes beat Chordata (→ `Fish.png`)
- Amphibia beats Chordata (→ `Amphibians.png`)
- Testudines / Squamata / Serpentes / Crocodylia beat Chordata (→ `Reptiles.png`)

**Fallback behaviour:**
- If no specific rule matches but the lineage contains Chordata → `Chordates.png`
- If nothing matches at all → `Metazoa.png`
- If all lineage sources fail, a warning is logged and `Metazoa.png` is used.

**Audit/debugging:**
- The audit TSV includes an `image_source` column showing which lineage source was used and which taxonomy term matched (e.g. `metadata_lineage:Lepidoptera`, `busco_lineage:Insecta`).
- Debug-level logging shows the resolution path for each genome.

**Important constraints:**
- Icon selection must **never** use species-name parsing, sample descriptions, or any free-text metadata field.
- The species display name must always be the scientific name. It must not be polluted by strain, sample descriptions, or habitat text.

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
- **`config.py`**: Project configuration definitions (`ProjectConfig` dataclass) and factory function mapping project names to their data-fetch and schema rules.
- **`models.py`**: Shared data model (`GenomeMetadata` dataclass) used as the internal representation across all pipeline stages.
- **`yaml_renderer.py`**: Handles rendering `GenomeMetadata` models into validated dicts and performs actual FTP publishability checks.
- **`icon_resolver.py`**: Taxonomy-lineage-based icon resolver. Maps taxon IDs to icon filenames using NCBI lineage data, with per-run caching and `icons.txt` override support.
- **`ftp_client.py`**: Core logic for querying the EBI/Ensembl FTP servers.
- **`changelog.py`**: Compares old and new YAML outputs and generates human-readable and TSV changelogs.
- **`registry/gb_tracker.py`**: Client for GB registry pre-release discovery.
- **`registry/metadata_db.py`**: Client for the metadata DB to look up released GUUIDs.
- **`registry/ncbi_entrez.py`**: NCBI web scraping for assembly submitters, population data, and parent-of-origin.
- **`bioproject_tracking.py`** (in `ensembl.genes.tracking`): An auxiliary tool for candidate discovery (to generate the input file), but *not* for final YAML eligibility.
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
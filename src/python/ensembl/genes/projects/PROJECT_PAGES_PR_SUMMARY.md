# Project Pages PR Summary (`update/project_pages`)

A "shopping list" of everything introduced by this PR. For day-to-day usage
see `README.md`; for the field-by-field metadata sourcing see
`docs/metadata_source_map.md`.

## A. Overview

This PR refactors project-page `species.yaml` generation away from the legacy,
core-DB-scraping scripts (`write_yaml.py`, `hprc_write_yaml.py`) into a single
metadata-driven, auditable pipeline. Genome metadata is pulled from the Ensembl
metadata DB and the Genebuild (GB) pre-release registry, FTP publishability is
validated, icons and alternate-haplotype links are resolved, and every decision
is recorded in an audit TSV. The legacy scripts remain only for backwards
compatibility and emit a `DeprecationWarning`.

## B. Main entrypoint — `generate_project_yaml.py`

```bash
python -m ensembl.genes.projects.generate_project_yaml <input_guuids.txt> \
  --project <project> --output species.yaml --audit-file audit.tsv
```

- **Input semantics**: one identifier per line.
  - A 36-char UUID → looked up in the metadata DB (released genomes).
  - Anything else (e.g. a core DB name) → looked up in the GB registry
    (pre-release). Core-DB fallback is deprecated.
- **Discovery**: in addition to the input file, project-scoped pre-release
  genomes are auto-discovered from the GB registry (see D).
- **Output**: `--output` species.yaml (sorted by species/assembly).
- **Optional**: `--audit-file` (decision log TSV), `--changelog` /
  `--changelog-output` (diff vs a previous YAML).

## C. Metadata DB integration — `registry/metadata_db.py`

`MetadataDbClient` queries `ensembl_genome_metadata`. Fields sourced:
- species (scientific name), assembly accession, assembly name
- genome UUID (drives the beta link)
- annotation method / annotation source
- geneset last-update date
- BUSCO score + BUSCO lineage dataset
- taxonomy id and (where available) taxonomy lineage
- `get_genome_uuids_by_accessions()` — batch accession → UUID lookup used to
  turn alternate-haplotype accessions into beta links.

## D. GB registry pre-release integration — `registry/gb_tracker.py`

`GbTrackerClient` queries `gb_assembly_metadata`:
- `fetch_by_identifier()` — single core-DB/accession lookup (now also returns
  BUSCO score/lineage so icon resolution has a fallback signal).
- `fetch_project_pre_releases()` — discovers `gb_status = 'pre_released'`
  genomes scoped by `ProjectConfig` (`bioproject_scoping` /
  `custom_group_scoping`).
- Pre-release rows only become YAML rows when FTP files actually exist (see E);
  otherwise they are excluded or shown as `beta_link: Coming soon!`.

## E. FTP validation & dynamic date resolution — `yaml_renderer.py`, `ftp_client.py`

- Released FTP existence check on `ftp.ebi.ac.uk/pub/ensemblorganisms/...`.
- `annotation_source`-aware paths (e.g. `ensembl`, `braker`).
- Dynamic `YYYY_MM` date resolution: if the metadata date does not match the
  public FTP, the directory is listed and the newest matching date is used.
- Pre-release FTP fallback (`pub/databases/ensembl/pre-release/...`).
- A genome with no valid released or pre-release files is **excluded**
  (`audit_decision = excluded`).

## F. Beta link validation — `ftp_client.check_beta_species_status`

- Fetches the beta species page and inspects the **body**, because the beta
  "species not recognised" page still returns HTTP 200.
- Treated as unavailable on non-200 or on known markers
  ("We do not recognise the species identified",
  "Find available species in the Species selector").
- Only a confirmed page yields a real `beta_link`; otherwise
  `beta_link: Coming soon!`.
- **Never excludes** a genome that has valid FTP files.
- Pre-release genomes skip the check; results are cached per UUID.

## G. Changelog & audit — `changelog.py`, `--audit-file`

- `--audit-file`: one row per candidate with decision + reason, including
  `image_source`, `beta_link_status`, final `beta_link`, `taxon_id`,
  `busco_lineage`, `taxonomy_lineage`, and final `image`.
- `--changelog <old.yaml>`: human-readable added/removed/modified report;
  `--changelog-output` writes the machine-readable TSV.
- Decisions include: `included_released`, `included_prerelease`, `excluded`,
  `excluded_duplicate`, `kept_duplicate`.

## H. Duplicate handling — `generate_project_yaml.py`

- Candidates are grouped by accession (+ species + source for non-HPRC).
- Within a group the best record is kept based on date match and release
  status; the rest are marked `excluded_duplicate` / `kept_duplicate`.
- Guarantees no duplicate YAML rows for the same accession.

## I. Icon resolution — `icon_resolver.py`

- Resolves an icon from taxonomy lineage, strongest source first:
  1. metadata DB taxonomy lineage, 2. NCBI Entrez lineage (cached by taxon_id),
  3. BUSCO lineage hint (e.g. `eudicotyledons_odb12` → Plants.png), 4. default
  `Metazoa.png`.
- Most-specific match wins (Lepidoptera > Arthropoda; Mammalia/Aves/Fish/
  Reptiles/Amphibia > Chordata).
- `icons.txt` provides project overrides on top of built-in taxonomy anchors.
- Icons are **never** inferred from species names or free-text fields.

## J. Alternate haplotypes — `haplotype_resolver.py`

- `find_alternate_haplotypes()` pairs assemblies within the dataset using:
  1. shared NCBI BioSample accession (primary),
  2. shared sample/isolate name within the same taxon (secondary),
  3. assembly-name patterns (hap1/hap2, mat/pat, primary/alternate) — weak
     fallback, same-species only.
- The alternate accession is resolved to its genome UUID (local cache, then a
  batch metadata-DB lookup) and emitted as a beta link
  (`alternate: https://beta.ensembl.org/species/<uuid>`).
- If no UUID exists the `alternate` field is **omitted** (no raw accession, no
  broken link). Only assemblies present in the dataset are linked.

## K. Project-specific output support — `config.py`, `yaml_renderer.py`

- **Standard** projects (CBP, BGE, DToL, VGP, ASG, ERGA, AEGIS): image,
  accession, annotation method, FTP links, BUSCO, beta link, alternate.
- **HPRC** (`--project hprc`): `assembly_submitter`, `assembly_link`,
  `parent_of_origin` (from the NCBI assembly report), and `variants_vep`
  when present on FTP.
- **mouse_genomes** (`--project mouse_genomes` / `mouse`): explicit `strain`
  field and `alternate` haplotype linkage.

## L. Tracking script changes — `tracking/bioproject_tracking.py`

- Candidate-discovery tool (NOT the YAML-eligibility authority).
- `--project_name hprc` discovers HPRC release-2 GCAs from the HPRC catalog.
- Optional `--ftp` adds `ftp_status` columns; `--pre_release` adds GB matches.
- Config is now loaded lazily and resolved relative to the module, so the
  module imports and `--help` work from any working directory.

## M. Configuration — `config.py`, `server_config.json`

- `ProjectConfig` directs data sources, schema type, and pre-release scoping.
- `get_project_config()` maps project keys → config (hprc, mouse_genomes/mouse,
  aquafaang, vgp, darwin_tree_of_life, erga, cbp, bge, asg, and a sensible
  default for other standard projects such as AEGIS).
- `server_config.json` holds read-only DB connection details (`meta_beta`,
  `gb1`, etc.).

## N. Tests — `tests/ensembl/genes/projects/`

| File | Covers |
| --- | --- |
| `test_changelog.py` | YAML diff: added/removed/modified, date/slash normalisation, HPRC key handling |
| `test_icon_resolver.py` | metadata/NCBI/BUSCO source priority, most-specific match, plant BUSCO tokens, lineage normalisation, pre-release plant regression |
| `test_haplotype_resolver.py` | BioSample / sample-name / naming pairing, priority, species-safety, edge cases |
| `test_beta_link.py` | beta status (available/unavailable/error), caching, `_resolve_beta_link` branches, end-to-end render |

All tests mock network/DB, so they run offline. **Not covered locally** (require
cluster/internal access): live metadata DB / GB registry queries, live EBI FTP
existence checks, live NCBI Datasets/Entrez lookups, and live beta.ensembl.org
validation. Run the commands in section O on the cluster to exercise those.

## O. Validation commands (cluster)

```bash
export PYTHONPATH=$PWD/src/python

python -m ensembl.genes.projects.generate_project_yaml ./cbp_guuids.txt \
  --project cbp   --output cbp_species.yaml   --audit-file cbp_audit.tsv

python -m ensembl.genes.projects.generate_project_yaml ./bge_guuids.txt \
  --project bge   --output bge_species.yaml   --audit-file bge_audit.tsv

python -m ensembl.genes.projects.generate_project_yaml ./hprc_guuids.txt \
  --project hprc  --output hprc_species.yaml  --audit-file hprc_audit.tsv

python -m ensembl.genes.projects.generate_project_yaml ./mouse_guuids.txt \
  --project mouse_genomes --output mouse_species.yaml --audit-file mouse_audit.tsv
```

Spot-check examples:

```bash
grep -A8 "Geum rivale" bge_species.yaml          # plant icon
grep "GCA_964205265.1" bge_audit.tsv             # image_source / busco columns
grep -A15 "Achillea ptarmica" aegis_species.yaml # beta_link "Coming soon!" case
```

## P. Known limitations

- Live metadata DB, GB registry, EBI FTP, NCBI and beta.ensembl.org checks all
  require network / internal access; they cannot be exercised in offline CI.
- NCBI / Datasets lookups can fail or rate-limit; the pipeline falls back
  (e.g. BUSCO-based icon resolution) where possible and logs the path taken.
- Output YAML still depends on the `projects.ensembl.org` species.yaml table
  conventions (field names such as `image`, `beta_link`, `alternate`).

# Metadata Source Map

This document explicitly maps the metadata fields required by the project YAML generation pipeline to their exact sources in the Ensembl Metadata databases (`schema.sql` & `gb_schema.sql`), ensuring zero reliance on legacy core database parsing unless mathematically impossible.

| Output Field | Current Source (Legacy Code) | New Source DB & Table(s) | New Source Column(s) | Fallback Source | Notes / Assumptions |
| --- | --- | --- | --- | --- | --- |
| **Genome UUID** | N/A (Hardcoded or missing) | `schema.sql: genome` | `genome.genome_uuid` | None | Primary identifier for API routing. |
| **Database Name** | Script argument / file list | `schema.sql: dataset_source`, `dataset`, `genome_dataset` | `dataset_source.name` | Core DB arg | Derived by joining `genome_dataset` -> `dataset_source`. |
| **Assembly Accession** | Core DB `meta` table (`assembly.accession`) | `schema.sql: assembly` <br> *(Pre-release:* `gb_schema.sql: genebuild_status`*)* | `assembly.accession` <br> *(Pre-release: `gca_accession`)* | Core DB `meta` | Fetched dynamically instead of checking the `meta` table. |
| **Assembly Name** | Core DB `meta` table (`assembly.name`) | `schema.sql: assembly` <br> *(Pre-release:* `gb_schema.sql: assembly`*)* | `assembly.name` <br> *(Pre-release: `asm_name`)* | Core DB `meta` | |
| **Taxonomy ID** | Core DB `meta` table | `schema.sql: organism` | `organism.taxonomy_id` | Core DB `meta` | Required for phylogenetic tree grouping. |
| **Strain / Isolate** | Core DB `meta` table (`species.strain`) | `schema.sql: organism` | `organism.strain` | Core DB `meta` | Mapped precisely in schema instead of suffix matching species name. |
| **Pre-release Status** | `information_schema.schemata` database regex | `gb_schema.sql: genebuild_status` | `genebuild_status.gb_status` | Regex scanning | Formally eliminates fragile DB list scraping. We now look up `gb_status` = 'pre_released' or 'in_progress'. |
| **Alt Haplotype Link** | Assembly name suffix heuristics (`_mat`) | `schema.sql: organism` + `assembly` | Matching `organism.taxonomy_id` & `strain` w/ different `assembly.accession` | Name parsing | Relies on finding sibling entries in `schema.sql` with identical taxonomy & strain parameters. |
| **Project Membership**| Core DB list files (`hprc_cores.txt`) | `schema.sql: organism_group` | `organism_group.name` | List files | Mapped via `organism_group_member`. |

## Implementation Directives
1. **No Core DB Defaults:** The pipeline must query `ensembl_genome_metadata` (via `schema.sql`) first using the `genome_uuid` or `assembly_accession`.
2. **Pre-release Routing:** If the UUID/Accession is absent in `schema.sql`, query `gb_assembly_metadata` (via `gb_schema.sql`) tables: `assembly` JOIN `genebuild_status` JOIN `species`.
3. **Core DB Fallback:** The `core_db.py` module should only execute if `gb_schema` also yields no results, and it MUST emit an explicit warning log noting the fallback.

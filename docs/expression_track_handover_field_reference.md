# Expression-track handover field reference

This is one real run-level record copied from the generated pig-pilot handover
JSON. It is followed by an explanation of what each field is for and where the
value came from.

The example is a draft handover. Run-level `track_source` entries use the ENA
Browser URL convention from the earlier handover payloads. They are generated
from the run accession in draft mode and are not evidence that the alignment
submission has been verified. Aggregate records use a clearly marked project
placeholder until an overarching source URL is supplied. Final handover
generation rejects `--dummy-track-source`.

```json
{
  "assembly_uuid": "da38d82b-df50-419b-af74-886463b7bfa3",
  "dataset_attribute": [
    {
      "name": "expression.source",
      "value": "ENA/Genebuild pig pilot"
    }
  ],
  "dataset_source": {
    "name": "/hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/run_level/GCA_000003025.6/DRR169427/coverage_forward/DRR169427.forward.bw",
    "type": "bigwig"
  },
  "dataset_type": "expression_tracks",
  "description": "Sus scrofa RNA-seq rnaseq coverage track for run DRR169427 (forward strand) on GCA_000003025.6. Genebuild data-check pilot.",
  "label": "pig-pilot-20261104-r1",
  "name": "sus_scrofa_pilot_run_DRR169427_coverage_forward",
  "team_responsible": "genebuild",
  "track_files": {
    "summary": "/hps/nobackup/flicek/ensembl/genebuild/jackt/expression_tracks/pig-pilot/run_level/GCA_000003025.6/DRR169427/coverage_forward/DRR169427.forward.bw"
  },
  "track_metadata_payload": {
    "data_type": "rnaseq",
    "display_label": "DRR169427",
    "metadata": {
      "genotype": {
        "value": "Wild Type",
        "provenance": {
          "source": "ENA BioSample",
          "source_field": "genotype",
          "source_id": "SAMD00156646",
          "source_value": "wild type",
          "normalization": "explicit_display_label_rule"
        }
      },
      "organism": {
        "value": "Sus scrofa",
        "provenance": {
          "source": "ENA BioSample",
          "source_field": "organism",
          "source_id": "SAMD00156646",
          "source_value": "Sus scrofa"
        }
      },
      "sex": {
        "value": "Female",
        "provenance": {
          "source": "ENA BioSample",
          "source_field": "sex",
          "source_id": "SAMD00156646",
          "source_value": "female",
          "normalization": "explicit_display_label_rule"
        }
      },
      "study": {
        "value": "PRJDB7878",
        "provenance": {
          "source": "ENA run metadata",
          "source_field": "study_accession",
          "source_id": "DRR169427",
          "source_value": "PRJDB7878"
        }
      },
      "tissue": {
        "value": "MGC",
        "provenance": {
          "source": "ENA BioSample",
          "source_field": "tissue",
          "source_id": "SAMD00156646",
          "source_value": "MGC"
        }
      }
    },
    "provenance": {
      "biosample_accession": "SAMD00156646",
      "run_accessions": [
        "DRR169427"
      ],
      "study_accession": "PRJDB7878"
    },
    "schema_version": "transcriptomic-selection-0.3.1",
    "selection_id": "sus_scrofa_pilot_run_DRR169427_coverage_forward",
    "selection_hierarchy": [
      "tissue",
      "biosample",
      "run"
    ],
    "selection_level": "run",
    "target_genome": {
      "assembly": "GCA_000003025.6",
      "assembly_uuid": "da38d82b-df50-419b-af74-886463b7bfa3",
      "species": "Sus scrofa"
    }
  },
  "track_source": {
    "details": "Run DRR169427",
    "name": "ENA",
    "url": "https://www.ebi.ac.uk/ena/browser/view/DRR169427"
  },
  "user": "jackt@ebi.ac.uk"
}
```

## Outer handover fields

| Field | Purpose | Source in this record |
| --- | --- | --- |
| `team_responsible` | Identifies the submitting Ensembl team. | Command configuration: `--team`, default `genebuild`. |
| `user` | Identifies the submitting user. | Command configuration: `--user`. |
| `name` | Stable unique name for this submitted physical track. | Deterministically generated from species, run accession, and product slug. |
| `description` | Human-readable description for the track. | Generated from normalized species, run, strand, assembly, and product fields; paths are deliberately excluded. |
| `label` | Groups records belonging to one handover preparation/release. | Command configuration: `--label`. |
| `dataset_type` | Tells the receiver that this is an expression-track dataset. | Renderer constant: `expression_tracks`. |
| `assembly_uuid` | Identifies the Ensembl assembly/genome receiving the track. | Resolved from `GCA_000003025.6` through the metadata DB using `METADATA_URI`; not hard-coded. |
| `dataset_attribute` | Optional catalogue attributes attached to the dataset. | Command configuration: `--dataset-attribute`; the attribute name must be registered by the receiver. |
| `dataset_source` | Identifies the physical source track submitted by this record. | `name` comes from the processing manifest; `type` is derived from the product mapping. |
| `track_files` | Maps the physical file to the browser role. | `summary` is selected because a BigWig is a summary track. BigBeds use `details`. |
| `track_source` | Optional external provenance/landing location for the source data. | Run records link to the ENA Browser run page; aggregate records require an overarching project URL. |
| `track_metadata_payload` | Carries the project’s richer transcriptomic metadata convention. | Rendered by the handover builder from the joined manifests; its schema is currently project-defined because Confluence leaves it unconstrained. |

## `track_metadata_payload` fields

| Field | Purpose | Source in this record |
| --- | --- | --- |
| `schema_version` | Version of the internal transcriptomic payload convention. | Renderer/project convention: `transcriptomic-selection-0.3.1`. |
| `selection_id` | Internal identity of the selected track product. | Same deterministic identity as the outer `name`. |
| `selection_level` | States what biological entity the record represents. | Derived from the product manifest: `run`. Tissue aggregates use `tissue`. |
| `data_type` | Describes the assay/data family. | Project convention: `rnaseq`. |
| `display_label` | Presentable label suitable for a selector or browser UI. | Run accession for run-level records; standardized aggregate label for tissue-level records. Machine identifiers remain in `name` and `selection_id`. |
| `target_genome` | Makes the target genome explicit inside the project payload. | Species and assembly accession from the selected-run/processing manifests; `assembly_uuid` is copied from the required outer field and validated for equality. |
| `metadata` | Browser-facing biological/filter metadata with an extensible value shape. | Field names use canonical snake-case names and values use field-specific display rules where defined. Each field is `{ "value": ... , "provenance": ... }`; relationship identifiers remain in the payload-level `provenance`. |
| `provenance` | Keeps the source relationships needed to trace the record back. | ENA run/BioSample/study accessions and aggregate identifiers. Source-BAM lineage belongs in the outer `track_source`. |
| `selection_hierarchy` | Describes the intended navigation/relationship order. | Command default: `tissue,biosample,run`; this is separate from per-field metadata and does not imply `hierarchy_level` values. |
| `schema_version` | Versions this project-defined metadata payload. | Renderer/project convention: `transcriptomic-selection-0.3.1`. |

## `metadata` fields

| Field | Purpose | Source in this record |
| --- | --- | --- |
| `organism.value` | Biological species facet. | BioSample `organism` attribute. |
| `genotype.value` | Genotype facet when explicitly supplied. | BioSample `genotype` attribute. |
| `sex.value` | Sex facet when explicitly supplied. | BioSample `sex` attribute. |
| `tissue.value` | Tissue/anatomy facet. | BioSample tissue/organism-part/anatomy attribute; missing values remain missing. |
| `study.value` | Study/project facet when useful for display or filtering. | ENA run metadata `study_accession`. The source identifier is also retained in `provenance.study_accession`. |

Each metadata field also has a `provenance` object containing the source system,
source field, source identifier where available, and source value. When an
explicit display-label rule changes the value, `normalization` records that
fact while `source_value` preserves the original text.

## Draft Pydantic specification

The generator defines `TrackMetadataPayload` as a Pydantic model and validates
each payload before writing the handover JSON. It currently enforces the
`run`/`tissue` selection levels, `rnaseq` data type, UUID-shaped assembly
identifiers, canonical metadata field names, and the provenance requirement
for normalized values. File format and role are represented by the outer
handover fields rather than repeated inside this payload. Extra fields remain
allowed because the Confluence contract does not yet define this nested schema.

The JSON Schema can be emitted with:

```bash
uv run python -c 'import json; from ensembl.genes.transcriptomic_data.create_expression_track_handover import TrackMetadataPayload; print(json.dumps(TrackMetadataPayload.model_json_schema(), indent=2))'
```

# Genebuild observability user guide

This guide describes the implemented draft in
`ensembl.genes.rgb.observability`. It is intentionally table-based so it can be
used before deeper pipeline decision logging exists.

The audit answers three practical questions:

1. What layer evidence or candidate models were not represented in core?
2. Which expected genes from a reference-informed catalogue are missing,
   degraded, assembly-limited, or cleanly retained?
3. Which review actions are most likely to improve completeness?

## Current implementation status

Implemented:

- Single-command DB-to-report execution with `rgb_observability run`.
- Evidence fate audit from RGB `extract/` tables.
- Optional expected-gene presence audit from user-supplied TSV/CSV/Parquet
  catalogue and projection files.
- Automatic expected-gene catalogue generation from another Ensembl core DB.
- Automatic expected-gene catalogue generation from expected annotation GFF3.
- Same-coordinate expected-projection generation when the reference annotation
  DB is already on the audited assembly coordinate system.
- Same-assembly transcript structure audit from GffCompare `.tmap`.
- BUSCO-like reference protein-set audit from arbitrary protein mapping hits.
- Failure summaries.
- Source, expected-source, copy-number, BUSCO, biotype, locus, and feature
  profile tables.
- Beyond-BUSCO completeness panels comparing all expected genes, high-confidence
  expected genes, non-BUSCO high-confidence genes, BUSCO-linked genes, and
  same-species reference genes.
- Ranked `review_loci.tsv`.
- Browser-ready `review_loci.bed`.
- Separate browser-ready BEDs by important failure class.
- Human-readable `actionable_summary.md`.
- Recommendation table mapping failure classes to concrete review actions.

Not implemented yet:

- Full exon/intron/CDS structure scoring.
- Direct parsing of Hive decision logs.
- Parameter replay or automatic reruns.
- True cross-assembly projection/liftover. The command can copy coordinates
  only when the expected-reference DB already uses the audited assembly.
- Direct execution of GffCompare, Liftoff, LiftOn, miniprot, or MMseqs/BLAST.
  The audit consumes their normalized outputs.

The current implementation still produces useful review sets because it can
identify high-confidence span-level failures and separate likely build issues
from likely assembly/projection issues. However, completeness claims are only as
strong as the expected-gene catalogue and projection file supplied.

## One-command run

For a user who has core and layer databases, this is the main command:

```bash
rgb_observability run \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <core_db> \
  --layer_db <layer_db> \
  --coord_system_name chromosome \
  --output_dir /path/to/observability-out \
  --expected_genes /path/to/expected_genes.tsv \
  --expected_projections /path/to/expected_projections.tsv \
  --top_n 1000 \
  --bed_pad_bp 500
```

This command does all of the following:

1. Extracts core and layer genes/transcripts/translations.
2. Builds strict and expanded loci.
3. Runs evidence fate auditing.
4. Runs expected-gene completeness auditing if expected inputs are supplied.
5. Writes review tables, BED tracks, completeness profiles, and
   recommendations.

If you do not already have an expected-gene file, the command can create the
catalogue from another core DB:

```bash
rgb_observability run \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <new_core_db> \
  --layer_db <new_layer_db> \
  --coord_system_name chromosome \
  --output_dir /path/to/observability-out \
  --expected_core_db <previous_or_reference_core_db> \
  --expected_source_name prior_ensembl \
  --expected_confidence high
```

This writes:

```text
/path/to/observability-out/<run_id>/expected/expected_genes.tsv
```

With only `--expected_core_db`, the audit knows which genes are expected but
does not know where they should map on the audited assembly. Those rows are
classified as `projection_unmapped` until a projection file is supplied.

The command can also create expected genes from a GFF3 annotation:

```bash
rgb_observability run \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <new_core_db> \
  --layer_db <new_layer_db> \
  --coord_system_name chromosome \
  --output_dir /path/to/observability-out \
  --expected_gff3 /path/to/reference_annotation.gff3.gz \
  --expected_gff3_source_name refseq_gff3 \
  --expected_gff3_confidence high
```

This writes:

```text
/path/to/observability-out/<run_id>/expected/expected_genes.tsv
```

If the GFF3 is already on the audited assembly coordinates, add:

```bash
--expected_gff3_projection_mode same_coordinates
```

That writes `expected_projections.tsv` from the GFF3 gene/transcript
coordinates. Do not use this for a GFF3 on another assembly. For a different
assembly, project/liftover/map the GFF3 first, or provide an external
`--expected_projections` file.

If the reference core DB is already on the audited coordinate system, for
example a same-assembly previous annotation or a pre-projected annotation DB,
you can also auto-populate the projection file:

```bash
rgb_observability run \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <new_core_db> \
  --layer_db <new_layer_db> \
  --coord_system_name chromosome \
  --output_dir /path/to/observability-out \
  --expected_core_db <reference_core_db_on_target_coordinates> \
  --expected_projection_mode same_coordinates
```

This writes both:

```text
/path/to/observability-out/<run_id>/expected/expected_genes.tsv
/path/to/observability-out/<run_id>/expected/expected_projections.tsv
```

Do not use `--expected_projection_mode same_coordinates` for a normal old
assembly to new assembly comparison. Old assembly coordinates are not target
coordinates. In that case, create projections with an external
projection/liftover/alignment workflow and pass them with
`--expected_projections`.

For same-assembly annotation comparison, standardize both annotations to GFF3 or
GTF, run GffCompare externally, then pass the `.tmap`:

```bash
rgb_observability run \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <core_db> \
  --layer_db <layer_db> \
  --coord_system_name chromosome \
  --output_dir /path/to/observability-out \
  --gffcompare_tmap /path/to/gffcompare.annotated.tmap
```

For a BUSCO-like test using another annotation's protein set, map those
proteins to the audited assembly or audited proteins externally, normalize the
hits to the reference protein hit schema below, then run:

```bash
rgb_observability run \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <core_db> \
  --layer_db <layer_db> \
  --coord_system_name chromosome \
  --output_dir /path/to/observability-out \
  --expected_proteins /path/to/expected_proteins.tsv \
  --reference_protein_hits /path/to/reference_protein_hits.tsv \
  --busco_complete_percent 80.0 \
  --require_expected_genes \
  --protein_min_identity 0.30 \
  --protein_min_query_coverage 0.70
```

The command prints the generated `run_id`. The output layout is:

```text
/path/to/observability-out/<run_id>/
  extract/
  loci/
  observability/
    actionable_summary.md
    recommendations.tsv
    review_loci.tsv
    review_loci.bed
    completeness_profile.tsv
    ...
```

If you already have an RGB run directory, use:

```bash
rgb_observability audit \
  --run_dir /path/to/rgb-out/<run_id> \
  --expected_genes /path/to/expected_genes.tsv \
  --expected_projections /path/to/expected_projections.tsv
```

Validate before a large run:

```bash
rgb_observability validate-inputs \
  --run_dir /path/to/rgb-out/<run_id> \
  --expected_genes /path/to/expected_genes.tsv \
  --expected_projections /path/to/expected_projections.tsv
```

Create empty expected-gene and reference-protein templates:

```bash
rgb_observability init-expected-template --out_dir /path/to/templates
```

## Inputs

### Required RGB inputs when using `audit`

Run the existing RGB extract/loci stages first:

```bash
rgb extract \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <core_db> \
  --layer_db <layer_db> \
  --coord_system_name chromosome \
  --output_dir /path/to/rgb-out \
  --run_id <run_id>

rgb loci \
  --output_dir /path/to/rgb-out \
  --run_id <run_id> \
  --locus_gap_bp 5000
```

The observability audit reads:

- `/path/to/rgb-out/<run_id>/extract/core_genes.tsv`
- `/path/to/rgb-out/<run_id>/extract/layer_genes.tsv`
- `/path/to/rgb-out/<run_id>/extract/core_transcripts.tsv`
- `/path/to/rgb-out/<run_id>/extract/layer_transcripts.tsv`
- `/path/to/rgb-out/<run_id>/extract/core_translations.tsv`
- `/path/to/rgb-out/<run_id>/extract/layer_translations.tsv`

The transcript and translation tables are optional but improve prioritization.

### Optional expected gene catalogue

Use this when you want same-species regression or reference-informed expected
gene assessment.

You can supply it directly with `--expected_genes`, or let
`rgb_observability run` generate it from another core DB with
`--expected_core_db`. The automatic
catalogue uses core `gene.stable_id` as `expected_gene_id`, sets
`expected_source` from `--expected_source_name`, copies the gene biotype, sets
`expected_copy_number=1`, and leaves symbol, orthogroup, and BUSCO fields blank.
Those fields can be improved later by joining to xrefs, Compara, BUSCO, or
curated source tables.

You can also generate it from GFF3 with `--expected_gff3`. GFF3 `gene` features
are preferred. If no gene features exist, transcript-like features are grouped
by `Parent` or `ID`; if those are absent, CDS features are grouped similarly.
The parser uses attributes such as `ID`, `gene_id`, `Name`, `gene_name`,
`biotype`, `gene_biotype`, `orthogroup_id`, and `busco_id` when present.

Minimum columns:

| Column | Meaning |
| --- | --- |
| `expected_gene_id` | Stable audit id for the expected gene |
| `expected_source` | Example: `prior_ensembl`, `refseq`, `busco`, `compara` |
| `reference_stable_id` | Source gene id |
| `symbol` | Optional gene symbol |
| `biotype` | Expected biotype, commonly `protein_coding` |
| `orthogroup_id` | Optional orthology/family id |
| `busco_id` | Optional BUSCO id if applicable |
| `expected_copy_number` | Expected copy count, default `1` |
| `confidence` | `high`, `medium`, `low`, or `unknown` |

Example:

```text
expected_gene_id	expected_source	reference_stable_id	symbol	biotype	orthogroup_id	busco_id	expected_copy_number	confidence
EXP0001	prior_ensembl	ENSGPREV0001	GENE1	protein_coding	OG0001		1	high
EXP0002	busco	BUSCO:123at456		protein_coding		123at456	1	medium
```

### Optional expected gene projections

Minimum columns:

| Column | Meaning |
| --- | --- |
| `expected_gene_id` | Must match catalogue |
| `seq_region_name` | Target sequence on audited assembly |
| `seq_region_start` | 1-based start |
| `seq_region_end` | 1-based end |
| `seq_region_strand` | `1`, `-1`, `+`, or `-` |
| `projection_status` | `mapped`, `partial`, `split`, `unmapped` |
| `projection_identity` | Optional 0-1 score |
| `projection_coverage` | Optional 0-1 score |
| `assembly_gap_overlap_bp` | Optional gap overlap |
| `repeat_overlap_bp` | Optional repeat overlap |

Example:

```text
expected_gene_id	seq_region_name	seq_region_start	seq_region_end	seq_region_strand	projection_status	projection_identity	projection_coverage	assembly_gap_overlap_bp	repeat_overlap_bp
EXP0001	1	100000	105000	1	mapped	0.998	0.99	0	0
EXP0002	1	250000	252000	-1	mapped	0.950	0.80	500	0
```

Projection files can be supplied directly with `--expected_projections`.
Automatic projection generation is intentionally limited to
`--expected_projection_mode same_coordinates`, which copies coordinates from the
expected core DB only when the seq_region name exists in the target audit run.
Seq_regions not present in the target set are written as `unmapped`.

For another assembly of the same species, this tool still needs a real
projection file from outside the audit command. Good sources include an
annotation projection pipeline, whole-genome alignment plus chain/liftover, or
a reference-protein/transcript mapping workflow that emits target intervals and
coverage/identity scores. Without that file the audit can say "this gene is in
the expectation catalogue", but it cannot honestly say whether the new assembly
contains it.

### Optional same-assembly GffCompare input

Use this when two annotations are on the same assembly and the question is gene
structure consistency. The audit currently consumes GffCompare `.tmap` rows and
writes transcript-level structure classes.

Useful command shape outside the audit:

```bash
gffcompare -r reference.standardized.gff3 -o compare query.standardized.gff3
```

Then pass `compare.query.standardized.gff3.tmap` with `--gffcompare_tmap`.

The audit classifies exact intron-chain matches as `exact_intron_chain`,
compatible or contained matches separately, and partial/overlap-only class
codes as structure review candidates. This is the right track for same-assembly
annotation regression because coordinates are already comparable and transcript
structure is the signal.

### Optional reference protein-set audit

Use this when you want a BUSCO-like presence test from a richer protein set:
previous Ensembl proteins, RefSeq proteins, manually curated proteins, UniProt,
or a species/clade protein panel.

`expected_proteins.tsv` minimum columns:

| Column | Meaning |
| --- | --- |
| `expected_gene_id` | Stable expected gene/protein group id |
| `query_protein_id` | Protein id from the reference set |
| `expected_source` | Example: `prior_ensembl_protein`, `refseq_protein`, `uniprot` |
| `reference_stable_id` | Source gene or protein id |
| `symbol` | Optional symbol |
| `biotype` | Usually `protein_coding` |
| `orthogroup_id` | Optional family id |
| `busco_id` | Optional BUSCO id |
| `confidence` | `high`, `medium`, `low`, or `unknown` |

`reference_protein_hits.tsv` minimum useful columns:

| Column | Meaning |
| --- | --- |
| `expected_gene_id` | Must match expected protein catalogue when available |
| `query_protein_id` | Reference protein id |
| `target_gene_id` | Optional target gene id from mapper |
| `target_stable_id` | Optional target gene stable id |
| `target_transcript_id` | Optional target transcript id |
| `seq_region_name` | Target sequence on audited assembly |
| `seq_region_start` | Target hit start |
| `seq_region_end` | Target hit end |
| `seq_region_strand` | Target hit strand |
| `aligner` | Example: `miniprot`, `diamond`, `blastp`, `mmseqs` |
| `hit_rank` | Best hit should be `1` |
| `percent_identity` | 0-1 identity |
| `query_coverage` | 0-1 query protein coverage |
| `target_coverage` | Optional 0-1 target coverage |
| `alignment_score` | Optional aligner score |
| `frameshift_count` | Optional coding disruption count |
| `stop_codon_count` | Optional coding disruption count |

If `expected_proteins.tsv` is omitted, the audit can classify proteins that
have hits, but it cannot count no-hit proteins because it does not know the
full denominator. For a BUSCO-like completeness statement, always provide the
expected protein catalogue.

## Run the audit

Without expected genes:

```bash
rgb_observability audit \
  --output_dir /path/to/rgb-out \
  --run_id <run_id> \
  --format tsv \
  --top_n 500
```

With expected genes:

```bash
rgb_observability audit \
  --output_dir /path/to/rgb-out \
  --run_id <run_id> \
  --format tsv \
  --expected_genes /path/to/expected_genes.tsv \
  --expected_projections /path/to/expected_projections.tsv \
  --top_n 1000 \
  --bed_pad_bp 500
```

During local development before reinstalling the package, the same command can
be run as:

```bash
PYTHONPATH=src/python python -m ensembl.genes.rgb.observability audit ...
```

## Outputs

The audit writes to:

`/path/to/rgb-out/<run_id>/observability/`

| File | Use |
| --- | --- |
| `audit_loci.tsv` | Locus-level counts for core, layer, expected genes, and review hits |
| `gene_to_locus.tsv` | Gene-to-locus mapping used by the audit |
| `evidence_fate.tsv` | One row per layer gene/model with best core representation and failure class |
| `expected_gene_presence.tsv` | One row per expected gene projection with presence/degradation class |
| `source_profile.tsv` | Layer source/biotype profile: orphan counts, P1 counts, median core coverage |
| `expected_source_profile.tsv` | Expected source/confidence profile: retained/missing/assembly-limited counts |
| `busco_expected_crosswalk.tsv` | BUSCO-linked expected genes and their broader presence classifications |
| `completeness_profile.tsv` | Completeness by expected-gene panel, not just BUSCO |
| `non_busco_high_confidence_losses.tsv` | High-confidence expected-gene losses outside BUSCO |
| `busco_proxy_calibration.tsv` | BUSCO-linked vs non-BUSCO completeness comparison |
| `copy_number_audit.tsv` | Expected vs observed copy number by orthogroup, BUSCO id, or expected gene id |
| `biotype_transition.tsv` | Layer-to-core and expected-to-core biotype transitions |
| `same_assembly_structure.tsv` | GffCompare `.tmap` rows classified for same-assembly structure comparison |
| `same_assembly_structure_summary.tsv` | Counts by GffCompare class-derived structure class |
| `reference_protein_audit.tsv` | BUSCO-like audit for arbitrary expected protein sets |
| `reference_protein_summary.tsv` | Counts by reference protein hit class |
| `release_readiness.tsv` | PASS/WARN/FAIL gates for release sign-off |
| `feature_profile.tsv` | Long-form top-level feature metrics with values, denominators, and fractions |
| `failure_mode_summary.tsv` | Counts by audit track, class, and priority |
| `recommendations.tsv` | Aggregated action suggestions with rationale and target table/filter |
| `review_loci.tsv` | Ranked actionable loci from both audit tracks |
| `review_loci.bed` | BED6 browser track for the ranked loci |
| `actionable_summary.md` | Short human-readable summary and top review list |

Separate class-specific BED files are written under:

`/path/to/rgb-out/<run_id>/observability/review_beds/`

Current files include:

- `no_core_gene_built.bed`
- `missing_with_evidence.bed`
- `projection_only.bed`
- `assembly_limited.bed`
- `present_degraded.bed`
- `present_wrong_biotype.bed`
- `split.bed`
- `fused.bed`

Completeness-specific BED files are written under:

`/path/to/rgb-out/<run_id>/observability/completeness_beds/`

Current files include:

- `non_busco_high_confidence_losses.bed`
- `high_confidence_actionable_losses.bed`
- `reference_protein_issues.bed`

The run manifest is:

`/path/to/rgb-out/<run_id>/manifest.observability.json`

## How to interpret evidence fate

Important columns in `evidence_fate.tsv`:

| Column | Meaning |
| --- | --- |
| `layer_stable_id` | Layer model being audited |
| `layer_logic_name` | Evidence source/provenance |
| `best_core_stable_id` | Best overlapping same-strand core gene |
| `layer_span_coverage_by_core` | Fraction of layer span overlapped by best core model |
| `fate_class` | Broad fate of the layer model |
| `failure_class` | More actionable reason |
| `review_priority` | `P1`, `P2`, `P3`, `P4` |
| `suggested_action` | Next review action |

High-value classes:

| Class | Interpretation | Action |
| --- | --- | --- |
| `no_core_gene_built` | Layer evidence exists but no same-strand core gene overlaps | Inspect locus, consider candidate rescue |
| `coding_evidence_not_protein_coding` | Coding layer evidence maps to a non-protein-coding core model | Check biotype assignment and CDS |
| `core_span_underrepresents_evidence` | Core captures only a small part of the layer model | Inspect structure, split, truncation, or threshold loss |
| `opposite_strand_core_only` | Only opposite-strand core genes overlap | Check strand/evidence correctness |
| `possible_transcript_or_model_collapse` | Multiple layer models collapse into one core gene | Inspect model collapse and isoform handling |

## How to interpret expected gene presence

Important columns in `expected_gene_presence.tsv`:

| Column | Meaning |
| --- | --- |
| `expected_gene_id` | Expected gene under audit |
| `confidence` | Strength of expectation |
| `projection_status` | Whether expected gene maps to the assembly |
| `best_core_stable_id` | Best overlapping core gene |
| `best_layer_logic_name` | Best local layer support if any |
| `expected_span_coverage_by_core` | Fraction of projected expected span overlapped by core |
| `presence_class` | Presence/degradation class |
| `failure_class` | Actionable failure reason |
| `suggested_action` | Next review action |

High-value classes:

| Class | Interpretation | Action |
| --- | --- | --- |
| `present_clean` | Expected gene is represented at span level | Usually no immediate action |
| `present_degraded` | Expected gene has a core match, but coverage is weak | Inspect structure and protein length |
| `present_wrong_biotype` | Expected coding gene maps to incompatible core biotype | Check biotype/CDS assignment |
| `missing_with_evidence` | Expected gene and layer evidence exist, but no core model | Candidate rescue or layer/threshold review |
| `projection_only` | Expected gene projects to assembly, but no core/layer support | Inspect browser and evidence ingest |
| `assembly_limited` | Missing expected gene overlaps assembly gap | Likely assembly issue, not only genebuild |
| `split` | One expected gene overlaps multiple core genes | Inspect split gene candidate |
| `fused` | Multiple expected genes share one core gene | Inspect fusion/readthrough candidate |

## Review workflow

Start with:

1. Open `actionable_summary.md`.
2. Load `review_loci.bed` in a browser.
3. Sort `review_loci.tsv` by `review_priority`, then `audit_track`.
4. For each P1/P2 locus, inspect the matching row in either
   `evidence_fate.tsv` or `expected_gene_presence.tsv`.

Then use the profile tables to decide whether the issue is local or systematic:

| Question | Table |
| --- | --- |
| Which layer sources produce the most high-priority orphan evidence? | `source_profile.tsv` |
| Which expected sources regress most? | `expected_source_profile.tsv` |
| Are BUSCO failures isolated or part of broader expected-gene loss? | `busco_expected_crosswalk.tsv` and `feature_profile.tsv` |
| Is BUSCO a decent proxy for this build? | `completeness_profile.tsv` and `busco_proxy_calibration.tsv` |
| Which high-confidence losses would BUSCO miss? | `non_busco_high_confidence_losses.tsv` |
| Are gene families collapsed or expanded? | `copy_number_audit.tsv` |
| Are coding genes becoming non-coding or otherwise changing biotype? | `biotype_transition.tsv` |
| On the same assembly, are transcript structures consistent? | `same_assembly_structure.tsv` |
| Do reference proteins support genes that core missed? | `reference_protein_audit.tsv` |
| Which feature metric should be tracked across builds? | `feature_profile.tsv` |
| Should this build be blocked before release? | `release_readiness.tsv` |
| What should we try next? | `recommendations.tsv` |

Suggested triage:

| Suggested action | What to do |
| --- | --- |
| `try_candidate_rescue` | Check layer model and see whether it should be promoted or rescued |
| `check_biotype_assignment` | Inspect CDS/protein evidence and final biotype |
| `inspect_browser` | Review core, layer, projection, repeats, and gaps visually |
| `check_assembly_gap` | Confirm the gene is limited by assembly sequence |
| `inspect_split_candidate` | Check whether one expected gene was split into multiple core genes |
| `inspect_fusion_candidate` | Check whether multiple expected genes were fused |
| `inspect_collapse_group` | Check whether competing models/isoforms were collapsed too aggressively |
| `try_candidate_rescue_from_protein` | Check whether a confident reference protein alignment should seed or rescue a model |
| `inspect_isoform_structure` | Compare intron chain, terminal exons, and canonical transcript choice |

## Recommendation framework

`recommendations.tsv` translates observed classes into review actions. It is not
proof that a given action will rescue a gene; it is a prioritized hypothesis
list.

| Recommendation type | Trigger | What to try |
| --- | --- | --- |
| `candidate_rescue` | Coding layer models have no same-strand core representation | Inspect layer model and candidate rescue path |
| `expected_gene_rescue` | Expected gene has local evidence but no core model | Review source evidence, candidate generation, thresholds, and layer ordering |
| `biotype_or_cds_review` | Coding evidence maps to non-protein-coding core model | Check CDS, translation, pseudogene, and biotype logic |
| `structure_review` | Core only partially spans layer evidence | Inspect split/truncation/collapse; later use exon/CDS metrics |
| `evidence_ingest_or_projection_review` | Expected gene projects but lacks local layer/core support | Check projection and whether evidence entered layer DB |
| `assembly_quality_review` | Expected gene overlaps assembly gap signal | Treat as assembly-limited unless other evidence contradicts |
| `copy_number_review` | Expected and observed copy number differ | Review paralogue/family collapse or expansion |
| `beyond_busco_review` | High-confidence non-BUSCO expected genes are lost | Do not rely on BUSCO alone; inspect non-BUSCO losses |
| `busco_proxy_warning` | Non-BUSCO loss rate exceeds BUSCO-linked loss rate | State BUSCO is not representative for this build |
| `protein_supported_candidate_rescue` | Reference proteins map confidently but no core gene overlaps | Use protein alignments as rescue/evidence-ingest candidates |
| `protein_integrity_review` | Protein mapping has frameshift/stop signals | Separate assembly/projection disruption from final-model truncation |
| `transcript_structure_review` | Same-assembly transcript structures are not exact matches | Use GffCompare class codes to focus structural review |

Use the `target_table` and `target_filter` columns to find the rows behind each
recommendation.

## Release recovery workflow

When an automated genebuild has low BUSCO, missing expected multi-copy genes, or
user-reported missing families, do not start by rerunning the whole pipeline.
Start by freezing the current build and classifying the failure.

Run the audit with the strongest available expectations:

```bash
rgb_observability run \
  --host <host> \
  --port <port> \
  --user <user> \
  --password <password> \
  --core_db <core_db_under_review> \
  --layer_db <layer_db_under_review> \
  --coord_system_name chromosome \
  --output_dir /path/to/recovery-audit \
  --expected_genes /path/to/expected_genes.tsv \
  --expected_projections /path/to/expected_projections.tsv \
  --expected_proteins /path/to/expected_proteins.tsv \
  --reference_protein_hits /path/to/reference_protein_hits.tsv \
  --busco_complete_percent 80.0 \
  --require_expected_genes \
  --top_n 2000 \
  --bed_pad_bp 1000
```

Then triage in this order:

1. Open `release_readiness.tsv`. Any `FAIL` row is a release blocker unless it
   is explicitly reclassified as assembly-limited with evidence.
2. Open `recommendations.tsv`. Work through P1 rows first.
3. Open `review_loci.bed` and class-specific BEDs in the browser.
4. For user-reported missing multicopy genes, start with
   `copy_number_audit.tsv`, `reference_protein_audit.tsv`, and
   `expected_gene_presence.tsv`.
5. Classify each missing copy as one of:
   `assembly_limited`, `projection_or_mapping_failure`, `evidence_not_ingested`,
   `candidate_generated_but_rejected`, `collapsed_copy`, `wrong_biotype`,
   `structure_degraded`, or `true_absence_uncertain`.

Recovery actions should be targeted:

| Failure signal | Recovery action |
| --- | --- |
| `protein_supported_no_core_gene` | Rescue or seed a candidate from the protein alignment; check why protein evidence was not represented |
| `missing_with_evidence` | Promote/rescue layer candidate or adjust source/layer ordering |
| `collapsed_copy_number` | Inspect paralogue collapse, projection ambiguity, and competition thresholds |
| `projection_only` | Check projection quality and evidence ingest before calling true loss |
| `assembly_limited` | Separate from genebuild failure and report as assembly-limited |
| `coding_evidence_not_protein_coding` | Review CDS, pseudogene logic, and biotype assignment |
| `same_assembly_structure_regression` | Compare intron chains/canonical choice and repair structure regressions |

After rescue, rerun the audit on the patched build. The recovery is not complete
until the relevant `release_readiness.tsv` gates move from `FAIL` to `PASS` or
to a documented assembly-limited exception.

## Prevention gates

The release process should fail before public release when any configured
release gate fails:

| Gate | Default behavior |
| --- | --- |
| `busco_floor` | Fails if `--busco_complete_percent` is below `--busco_floor_percent` |
| `expected_catalogue_required` | Fails when `--require_expected_genes` is set and no expected genes were audited |
| `high_confidence_actionable_loss` | Fails on any high-confidence actionable expected-gene loss by default |
| `expected_genes_missing_with_evidence` | Fails when expected genes have evidence but no core model |
| `copy_number_regression` | Fails on any copy-number issue by default |
| `protein_supported_missing_core` | Fails when reference proteins map confidently but no core gene overlaps |
| `same_assembly_structure_regression` | Warns on same-assembly structural differences |

For early experimentation you can loosen thresholds, for example:

```bash
--busco_floor_percent 90 \
--max_high_confidence_actionable_loss_fraction 0.01 \
--max_copy_number_issue_fraction 0.02
```

For release sign-off, the stricter default is intentional: a known missing
high-confidence gene or collapsed key family should be an explicit exception,
not something hidden behind a passable aggregate BUSCO number.

## What makes an insight actionable?

An output row is actionable when it has:

- coordinates
- source/provenance
- best core match or explicit absence
- confidence or coding evidence
- failure class
- suggested action

Examples:

- A high-confidence previous same-species gene is projected cleanly, has layer
  evidence, and has no core model: candidate rescue or layer selection issue.
- A BUSCO marker fails protein mode, but the broader expected catalogue shows
  many non-BUSCO same-species genes lost too: broad build regression, not an
  isolated BUSCO event.
- A projected gene is missing but overlaps a large assembly gap: assembly-limited
  absence, not immediately a genebuild parameter issue.
- A coding layer model is represented by a non-coding core gene: biotype/CDS
  assignment issue.

## Feature profile metrics

`feature_profile.tsv` is a stable long-form table intended for comparisons
across builds. Each row has:

- `metric_group`
- `metric_name`
- `value`
- `denominator`
- `fraction`
- `description`

Important current metrics:

| Metric | Meaning |
| --- | --- |
| `layer_model_count` | Total layer models audited |
| `orphan_layer_model_count` | Layer models with no same-strand core representation |
| `coding_orphan_layer_model_count` | Coding layer models with no same-strand core representation |
| `p1_evidence_issue_count` | Layer evidence rows requiring high-priority review |
| `expected_gene_count` | Expected genes audited |
| `present_clean_count` | Expected genes represented cleanly at span level |
| `missing_with_evidence_count` | Expected genes with local evidence but no core model |
| `projection_only_count` | Expected genes projected to assembly but lacking core/layer support |
| `assembly_limited_count` | Expected genes blocked by assembly gap signal |
| `present_degraded_count` | Expected genes with weak span representation |
| `p1_expected_issue_count` | Expected-gene rows requiring high-priority review |
| `copy_number_issue_count` | Copy groups not matching expected copy number |
| `busco_expected_gene_count` | Expected genes carrying BUSCO ids |
| `busco_p1_or_p2_issue_count` | BUSCO-linked expected genes needing P1/P2 review |
| `high_confidence_present_clean_count` | Clean high-confidence expected genes |
| `high_confidence_actionable_loss_count` | Actionable high-confidence expected-gene losses |
| `high_confidence_non_busco_present_clean_count` | Clean high-confidence expected genes outside BUSCO |
| `high_confidence_non_busco_actionable_loss_count` | Actionable non-BUSCO high-confidence losses |
| `busco_linked_present_clean_count` | Clean BUSCO-linked expected genes |
| `busco_linked_actionable_loss_count` | Actionable BUSCO-linked expected-gene losses |
| `non_busco_high_confidence_loss_count` | High-confidence non-BUSCO expected genes with actionable loss/degradation |

These metrics are deliberately coarse. The next framework layer should add
exon/intron/CDS metrics without changing the existing metric names.

## Completeness beyond BUSCO

The completeness framework treats BUSCO as one panel in a broader expected-gene
catalogue. The key output is `completeness_profile.tsv`.

Current panels:

| Panel | Meaning |
| --- | --- |
| `all_expected` | Every expected gene supplied |
| `high_confidence` | Expected genes marked `confidence=high` |
| `high_confidence_non_busco` | High-confidence expected genes without a BUSCO id |
| `busco_linked` | Expected genes with a BUSCO id |
| `same_species_reference` | Expected genes from prior/same-species annotations |
| `source:<name>` | Per-source panels such as `source:prior_ensembl` or `source:refseq` |

Important completeness columns:

| Column | Meaning |
| --- | --- |
| `present_any_fraction` | Fraction with any core representation, including degraded/split/fused |
| `present_clean_fraction` | Fraction represented cleanly at current span-level resolution |
| `actionable_loss_fraction` | Fraction missing/degraded in a way that needs review |
| `assembly_limited_fraction` | Fraction likely limited by assembly gaps |

Use `busco_proxy_calibration.tsv` to compare BUSCO-linked genes with
high-confidence non-BUSCO expected genes. If non-BUSCO actionable loss is high
while BUSCO-linked loss is low, BUSCO is giving false reassurance for that
build. If both are high, the build likely has broad expected-gene completeness
problems.

## Caveats

- Current representation is span-level. Exon/intron/CDS precision is the next
  implementation step.
- Expected-gene quality depends on the quality of the catalogue and projection
  file.
- Same-assembly GffCompare input is transcript-structure evidence, not a
  cross-assembly projection.
- Reference protein-set audit is only BUSCO-like if you provide the full
  expected protein denominator, including proteins with no hits.
- `missing_no_evidence` is not proof of true absence.
- BUSCO should be one expected source, not the definition of completeness.
- Recommendations are hypotheses. They become stronger once decision logs,
  exon/CDS metrics, and rerun experiments are attached.

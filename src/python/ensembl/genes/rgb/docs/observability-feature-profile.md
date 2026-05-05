# Observability feature profile

Status: implemented draft

The feature profile is the stable metric layer for comparing genebuilds. The
full observability audit emits detailed row-level tables, but the feature
profile condenses those tables into metrics that can be tracked across builds,
assemblies, species, sources, and releases.

## Goals

- Provide a stable set of metrics before running on real data.
- Separate evidence-fate, expected-gene, BUSCO, copy-number, and review signals.
- Make it obvious whether a problem is local or systematic.
- Keep BUSCO in context rather than treating it as the full completeness score.
- Leave room for exon/intron/CDS metrics without changing the first metric
  contract.

## Generated framework outputs

The implemented draft writes these files under:

`<output_dir>/<run_id>/observability/`

| File | Status | Purpose |
| --- | --- | --- |
| `audit_loci.tsv` | implemented | Locus-level core/layer/expected/review counts |
| `gene_to_locus.tsv` | implemented | Mapping table copied/generated for downstream joins |
| `evidence_fate.tsv` | implemented | Layer model fate and failure classes |
| `expected_gene_presence.tsv` | implemented | Expected-gene presence/degradation classes |
| `source_profile.tsv` | implemented | Layer source and biotype aggregate profile |
| `expected_source_profile.tsv` | implemented | Expected-source and confidence aggregate profile |
| `busco_expected_crosswalk.tsv` | implemented | BUSCO-linked expected genes in broader context |
| `completeness_profile.tsv` | implemented | Completeness by expected-gene panel beyond BUSCO |
| `non_busco_high_confidence_losses.tsv` | implemented | High-confidence expected-gene losses BUSCO would miss |
| `busco_proxy_calibration.tsv` | implemented | BUSCO-linked vs non-BUSCO completeness comparison |
| `copy_number_audit.tsv` | implemented | Expected vs observed copy count by copy group |
| `biotype_transition.tsv` | implemented | Layer/expected biotype to core biotype transitions |
| `same_assembly_structure.tsv` | implemented | GffCompare `.tmap` structure comparison for same assembly |
| `same_assembly_structure_summary.tsv` | implemented | Same-assembly structure class counts |
| `reference_protein_audit.tsv` | implemented | BUSCO-like audit for arbitrary reference protein sets |
| `reference_protein_summary.tsv` | implemented | Reference protein hit class counts |
| `release_readiness.tsv` | implemented | PASS/WARN/FAIL release gates with required actions |
| `feature_profile.tsv` | implemented | Stable long-form metric table |
| `failure_mode_summary.tsv` | implemented | Counts by class and priority |
| `recommendations.tsv` | implemented | Aggregated action hypotheses and target filters |
| `review_loci.tsv` | implemented | Ranked actionable review loci |
| `review_loci.bed` | implemented | Browser track for all ranked review loci |
| `review_beds/*.bed` | implemented | Browser tracks split by failure class |
| `actionable_summary.md` | implemented | Human-readable first-pass report |

The one-command runner can also create `<output_dir>/<run_id>/expected/`
inputs from `--expected_core_db`. It writes an expected-gene catalogue from the
reference DB and, only with `--expected_projection_mode same_coordinates`, a
same-coordinate projection table. True cross-assembly projection remains an
external input because old assembly coordinates are not valid target
coordinates.

The same input directory can be generated from an expected annotation GFF3 with
`--expected_gff3`. As with core DB input, `--expected_gff3_projection_mode
same_coordinates` is valid only when the GFF3 coordinates already refer to the
audited assembly.

Same-assembly structure comparison and different-assembly protein evidence are
separate tracks. Same assembly should normally use standardized GFF/GTF plus
GffCompare. Different assembly should use projected expected genes and/or
reference protein mappings from tools such as Liftoff/LiftOn/miniprot before
the audit is run.

## `feature_profile.tsv` schema

| Column | Meaning |
| --- | --- |
| `metric_group` | High-level domain: evidence_fate, expected_presence, review, copy_number, busco |
| `metric_name` | Stable metric identifier |
| `value` | Metric numerator or count |
| `denominator` | Denominator where a fraction is meaningful |
| `fraction` | `value / denominator`, otherwise empty |
| `description` | Human-readable definition |

## Current metrics

### Evidence fate

| Metric | Meaning | Actionable when high |
| --- | --- | --- |
| `layer_model_count` | Number of layer models audited | Baseline denominator |
| `orphan_layer_model_count` | Layer models with no same-strand core representation | Evidence not making it into core |
| `coding_orphan_layer_model_count` | Coding layer models with no same-strand core representation | High-priority candidate rescue pool |
| `p1_evidence_issue_count` | Layer evidence rows marked P1 | Immediate review burden |

### Expected gene presence

| Metric | Meaning | Actionable when high |
| --- | --- | --- |
| `expected_gene_count` | Expected genes audited | Baseline denominator |
| `present_clean_count` | Expected genes represented cleanly at span level | Should be high |
| `missing_with_evidence_count` | Expected genes with local evidence but no core model | Build/layer/threshold issue candidates |
| `projection_only_count` | Expected genes projected to assembly with no core/layer support | Evidence ingest or true absence review |
| `assembly_limited_count` | Expected genes blocked by assembly gap signal | Assembly problem, not only genebuild |
| `present_degraded_count` | Expected genes weakly represented by core span | Structural review candidates |
| `p1_expected_issue_count` | Expected-gene rows marked P1 | Immediate expected-gene review burden |

### Review

| Metric | Meaning |
| --- | --- |
| `review_locus_count` | P1/P2 review loci emitted |
| `p1_review_locus_count` | P1 review loci emitted |

### Copy number

| Metric | Meaning |
| --- | --- |
| `copy_group_count` | Orthogroup/BUSCO/expected-gene copy groups audited |
| `copy_number_issue_count` | Copy groups with missing, collapsed, or expanded copy number |

### BUSCO context

| Metric | Meaning |
| --- | --- |
| `busco_expected_gene_count` | Expected genes with a BUSCO id |
| `busco_p1_or_p2_issue_count` | BUSCO-linked expected genes needing P1/P2 review |

### Completeness beyond BUSCO

| Metric | Meaning |
| --- | --- |
| `high_confidence_present_clean_count` | Clean high-confidence expected genes |
| `high_confidence_actionable_loss_count` | Actionable high-confidence expected-gene losses |
| `high_confidence_non_busco_present_clean_count` | Clean non-BUSCO high-confidence expected genes |
| `high_confidence_non_busco_actionable_loss_count` | Actionable non-BUSCO high-confidence losses |
| `busco_linked_present_clean_count` | Clean BUSCO-linked expected genes |
| `busco_linked_actionable_loss_count` | Actionable BUSCO-linked expected-gene losses |
| `non_busco_high_confidence_loss_count` | High-confidence non-BUSCO expected genes with actionable loss/degradation |

### Reference Protein Sets

| Metric | Meaning |
| --- | --- |
| `expected_protein_count` | Reference proteins audited with BUSCO-like protein-set logic |
| `protein_supported_built_count` | Reference proteins with confident hits overlapping a core gene |
| `protein_supported_no_core_gene_count` | Reference proteins with confident assembly hits but no overlapping core gene |
| `protein_hit_degraded_count` | Reference protein hits with frameshift or stop-codon signals |

### Same-Assembly Structure

| Metric | Meaning |
| --- | --- |
| `gffcompare_transcript_count` | GffCompare query transcript rows audited |
| `exact_intron_chain_count` | Query transcripts with exact intron-chain matches |
| `p1_or_p2_structure_issue_count` | Same-assembly transcript rows requiring structure review |

## Release readiness gates

`release_readiness.tsv` turns audit observations into release blockers. It is
designed for cases where a build could otherwise be released with a weak
aggregate score and hidden regressions.

| Gate | Blocks when |
| --- | --- |
| `busco_floor` | Supplied BUSCO complete percentage is below the configured floor |
| `expected_catalogue_required` | Expected-gene audit is required but missing |
| `high_confidence_actionable_loss` | High-confidence expected genes have actionable loss/degradation |
| `expected_genes_missing_with_evidence` | Expected genes have local evidence but no core gene |
| `copy_number_regression` | Expected copy-number groups are missing, collapsed, or expanded |
| `protein_supported_missing_core` | Reference proteins map confidently but no core gene overlaps |
| `same_assembly_structure_regression` | Same-assembly structures need review; warning by default |

## Aggregate profile tables

### `source_profile.tsv`

Use this to identify layer sources that are not making it into core.

| Column | Meaning |
| --- | --- |
| `layer_logic_name` | Layer evidence source |
| `layer_biotype` | Layer biotype |
| `n_layer_models` | Source model count |
| `n_coding_models` | Source models with translations |
| `n_orphan_models` | Source models with no same-strand core representation |
| `n_p1_models` | High-priority source models |
| `median_core_coverage` | Median span coverage by best core gene |
| `top_failure_class` | Most common failure class |

### `expected_source_profile.tsv`

Use this to compare expected sources and confidence tiers.

| Column | Meaning |
| --- | --- |
| `expected_source` | prior_ensembl, refseq, busco, compara, etc. |
| `confidence` | high, medium, low, unknown |
| `n_expected_genes` | Expected-gene rows |
| `n_present_clean` | Cleanly represented expected genes |
| `n_missing_with_evidence` | Expected genes with evidence but no core model |
| `n_projection_only` | Projected expected genes without local support |
| `n_assembly_limited` | Assembly-limited expected genes |
| `n_p1_genes` | High-priority expected-gene rows |
| `top_presence_class` | Most common presence class |

### `copy_number_audit.tsv`

Use this to find family-level collapses or expansions.

Copy groups are built from `orthogroup_id`, then `busco_id`, then
`expected_gene_id`.

| Class | Meaning |
| --- | --- |
| `copy_number_as_expected` | Observed core copy count matches expectation |
| `missing_all_copies` | No core copy observed |
| `collapsed_copy_number` | Fewer core copies than expected |
| `expanded_copy_number` | More core copies than expected |

### `busco_expected_crosswalk.tsv`

Use this to decide whether BUSCO is a useful proxy for broader loss.

The important comparison is not only BUSCO complete vs missing. It is whether
BUSCO-linked expected genes have the same failure distribution as non-BUSCO
high-confidence expected genes.

### `completeness_profile.tsv`

Use this as the primary completeness table. It treats BUSCO as one expected-gene
panel, not as the whole assessment.

Panels:

| Panel | Meaning |
| --- | --- |
| `all_expected` | Every expected gene supplied |
| `high_confidence` | Expected genes with `confidence=high` |
| `high_confidence_non_busco` | High-confidence expected genes without BUSCO ids |
| `busco_linked` | Expected genes carrying BUSCO ids |
| `same_species_reference` | Previous/same-species expected genes |
| `source:<name>` | Per-source expected-gene panels |

Core completeness columns:

| Column | Meaning |
| --- | --- |
| `present_any_fraction` | Any core representation, including degraded/split/fused |
| `present_clean_fraction` | Clean representation at current span-level resolution |
| `actionable_loss_fraction` | Missing/degraded/split/fused rows needing review |
| `assembly_limited_fraction` | Likely assembly-limited rows |

### `non_busco_high_confidence_losses.tsv`

This is the table to open when BUSCO looks good but the annotation may still be
losing important genes. It contains high-confidence expected genes without
BUSCO ids that have actionable loss/degradation classes.

### `busco_proxy_calibration.tsv`

This compares the BUSCO-linked panel to the high-confidence non-BUSCO panel.

Interpretation:

- BUSCO loss high and non-BUSCO loss high: likely broad completeness problem.
- BUSCO loss high and non-BUSCO loss low: sentinel-specific or BUSCO lineage
  issue is plausible.
- BUSCO loss low and non-BUSCO loss high: BUSCO is giving false reassurance.
- BUSCO loss low and non-BUSCO loss low: no broad presence/absence signal at
  current resolution.

### `biotype_transition.tsv`

Use this to detect systematic biotype changes.

Tracks:

- `layer_to_core`: layer biotype to best core biotype
- `expected_to_core`: expected biotype to best core biotype

High counts of `protein_coding -> no_core` or `protein_coding -> lncRNA` should
be reviewed before release.

## Review BED framework

The audit writes one combined BED and class-specific BED files.

Class-specific BEDs let reviewers load only the current question:

- `no_core_gene_built.bed`: layer evidence not represented in core
- `missing_with_evidence.bed`: expected gene plus evidence, no core
- `projection_only.bed`: expected projection but no core/layer support
- `assembly_limited.bed`: likely assembly-limited expected gene
- `present_degraded.bed`: weak expected-gene representation
- `present_wrong_biotype.bed`: expected coding/non-coding mismatch
- `split.bed`: split candidate
- `fused.bed`: fusion/readthrough candidate

## Next feature layer

## Recommendation layer

`recommendations.tsv` is the bridge from audit observation to action. Each row
contains:

- `recommendation_id`
- `priority`
- `scope`
- `recommendation_type`
- `trigger_count`
- `rationale`
- `next_action`
- `target_table`
- `target_filter`
- `confidence`

This table should be treated as a hypothesis queue. For example, a
`candidate_rescue` recommendation says that coding layer models lack core
representation, so rescue/layer/threshold review is plausible. It does not prove
that changing a threshold will improve the build. That stronger claim needs
decision logs or targeted reruns.

Current recommendation types:

- `candidate_rescue`
- `expected_gene_rescue`
- `biotype_or_cds_review`
- `structure_review`
- `evidence_ingest_or_projection_review`
- `assembly_quality_review`
- `copy_number_review`
- `beyond_busco_review`
- `busco_proxy_warning`

The next implementation layer should add structure-aware metrics while keeping
this output contract stable:

- exact intron-chain representation
- near intron-chain representation with splice tolerance
- intron subset/superset classification
- single-exon coverage
- CDS coverage
- protein length ratio
- stop-codon/premature-stop flags
- canonical transcript regression
- BUSCO genome/protein status cross-tab
- parameter-threshold proximity from decision logs

# Stable ID Mapper

Tools for generating Ensembl stable-ID mapping decisions between an old
reference annotation and a new target annotation.

The current recommended workflow is `run_stable_id_mapping.py`. It is a
single-species pipeline designed to be moved into Nextflow later. It uses
LiftOn projection as the primary source of old-feature evidence, structural
matching against the target annotation as the main acceptance signal, and
coordinate overlap only as a fallback.

This directory is intentionally still self-contained and separate from the
packaged `src/python/ensembl/genes` modules.

## Core Rules

- Stable-ID mapping is decided from reference/LiftOn/structure/locus evidence.
- Target stable IDs are treated as current row labels for reporting and SQL,
  not as mapping evidence. A target row carrying the same stable ID as an old
  reference row does not automatically map.
- Mapped genes, transcripts, and translations keep the old stable ID and get
  `old_version + 1`.
- New genes, transcripts, and translations receive IDs from the supplied stable
  ranges and always get version `1`.
- Exons are not stable-ID mapped by this pipeline. Exon/CDS rows are used only
  as transcript-structure evidence.
- One old feature can claim at most one target feature, and one target feature
  can be claimed by at most one old feature. Competition produces missing/new
  decisions that can be reviewed in the audit tables.

## Data Flow

```text
reference FASTA + reference GFF3
target FASTA    + target GFF3
        |
        v
LiftOn projects reference genes onto the target assembly
        |
        v
structural matcher compares LiftOn output to target GFF3
        |
        v
decision layer assigns mapped / missing / new stable-ID decisions
        |
        v
review TSV + score evidence + dry-run or executable SQL
        |
        v
optional audit tables for missing, coordinate-only, and new genes
```

The DB name is used only in generated SQL as a `USE` statement. **The Python
pipeline does not connect to MySQL.** This will be implemented in the nextflow pipeline.

## Requirements

- Python 3
- LiftOn available as `lifton`, or pass `--lifton-executable`

For local development and tests from the repository root:

```bash
python3 -m pytest pipelines/stable_id_mapper/tests
```

## Main Command

Run a full single-species mapping:

```bash
python3 pipelines/stable_id_mapper/run_stable_id_mapping.py \
  --ref-fasta ref.fa \
  --ref-gff ref.gff3 \
  --target-fasta tar.fa \
  --target-gff tar.gff3 \
  --db-name species_core_test \
  --mapping-session-id X \
  --gene-range ENSXXXXG:90000000000-90000099999 \
  --transcript-range ENSXXXXT:90000000000-90000099999 \
  --translation-range ENSXXXXP:90000000000-90000099999 \
  --output-dir out/stable-id-run \
  --lifton-threads 16
```

The default rules are in:

```text
pipelines/stable_id_mapper/stable_id_mapping_rules.json
```

`stable_id_mapping_rules_locus_calibration.json` is kept as the calibration
copy used during real-data tuning. It currently has the same intended defaults.

## Command-Line Parameters

All options accepted by `run_stable_id_mapping.py` are listed below.

### Required Inputs

`--ref-fasta PATH`
: Old/reference genome FASTA. Passed to LiftOn as the reference genome.

`--ref-gff PATH`
: Old/reference annotation GFF3. Stable IDs and old versions come from here.

`--target-fasta PATH`
: New/target genome FASTA. Passed to LiftOn as the target genome.

`--target-gff PATH`
: New/target annotation GFF3. Target rows to be claimed or marked new come from
  here.

`--db-name NAME`
: Core database name to write into the SQL `USE` statement.

`--output-dir PATH`
: Output directory for LiftOn output, matching tables, reports, decision TSV,
  score evidence, and SQL.

### Required for now, will be automatic in next stage

`--mapping-session-id INT`
: Stable-ID mapping session ID to write into SQL and decision TSV rows.

`--gene-range PREFIX:START-END`
: Stable-space range for newly assigned gene stable IDs.

`--transcript-range PREFIX:START-END`
: Stable-space range for newly assigned transcript stable IDs.

`--translation-range PREFIX:START-END`
: Stable-space range for newly assigned translation stable IDs.

### Rules And Thresholds

`--rules-config PATH`
: JSON rules file. Defaults to
  `pipelines/stable_id_mapper/stable_id_mapping_rules.json`.

`--min-overlap FLOAT`
: Override `coordinate_overlap.min_overlap` from the rules file for this run.
  Applies to coordinate fallback for genes, transcripts, and translations.

`--match-window INT`
: Override `structural_matching.window`.

`--match-topk INT`
: Override `structural_matching.topk`.

`--match-min-score FLOAT`
: Override `structural_matching.min_score`.

`--match-good FLOAT`
: Override `structural_matching.good_score`.

`--match-confident FLOAT`
: Override `structural_matching.confident_score`.

`--match-gene-fraction FLOAT`
: Override `structural_matching.gene_fraction`.

### Reusing Existing Outputs

`--existing-lifton-gff PATH`
: Reuse an existing LiftOn-projected GFF3 instead of running LiftOn.

`--existing-transcript-pairs PATH`
: Reuse an existing `lifton.transcript_pairs.tsv`. Must be supplied together
  with `--existing-gene-pairs`.

`--existing-gene-pairs PATH`
: Reuse an existing `lifton.gene_pairs.tsv`. Must be supplied together with
  `--existing-transcript-pairs`.

If only `--existing-lifton-gff` is supplied, the pipeline regenerates structural
matching and decisions. If all three existing files are supplied, it skips both
LiftOn and structural matching and regenerates decisions/SQL only.

### LiftOn Options

`--lifton-threads INT`
: Threads passed to LiftOn. Default: `8`.

`--lifton-executable CMD`
: LiftOn executable name or path. Default: `lifton`.

`--lifton-feature-types LIST`
: Comma-separated parent feature types for LiftOn `-f`. Default is `gene`.
  Child feature types such as `mRNA`, `transcript`, `exon`, `CDS`, and UTRs are
  rejected by validation because stable-ID mapping should project parent gene
  loci only.

`--lifton-feature-types-file PATH`
: Existing feature-type file to pass to LiftOn `-f`.

`--extra-lifton-arg ARG`
: Extra argument passed to LiftOn before the FASTA positional arguments. This
  option can be repeated.

`--dry-run-lifton-command`
: Validate inputs, write the LiftOn feature-type file if needed, print the
  LiftOn command, and exit without running the full pipeline.

### SQL And Output Control

`--no-translations`
: Skip translation decisions and SQL.

`--write-executable-sql`
: Write executable SQL ending in `COMMIT`. Without this flag, SQL is dry-run
  review SQL ending in `ROLLBACK`.

`--replace-events-for-session`
: In executable SQL only, delete existing `stable_id_event` rows for the same
  `mapping_session_id` before inserting new events.

`--backup-prefix TEXT`
: Prefix for backup tables in executable SQL. Defaults to a timestamped value.

`--batch-size INT`
: Number of rows per SQL `INSERT` batch. Default: `500`.

## Rules File

The default rules file is:

```text
pipelines/stable_id_mapper/stable_id_mapping_rules.json
```

Current calibrated defaults:

```json
{
  "coordinate_overlap": {
    "min_overlap": 0.75
  },
  "structural_matching": {
    "window": 100000,
    "topk": 5,
    "min_score": 0.3,
    "good_score": 0.45,
    "confident_score": 0.6,
    "gene_fraction": 0.6,
    "score_weights": {
      "span_containment": 0.55,
      "query_coverage": 0.25,
      "intron_sim": 0.05,
      "jacc_internal": 0.03,
      "jacc_all": 0.03,
      "exon_count_sim": 0.02,
      "boundary_sim": 0.05,
      "lifton_identity_prior": 0.02
    }
  }
}
```

Every rules-file parameter:

`coordinate_overlap.min_overlap`
: Minimum coordinate-overlap score for fallback mapping when structural
  evidence does not claim a target. The score is overlap length divided by the
  longer feature length, with contig and strand required to match. This value
  applies to gene, transcript, and translation coordinate fallback. The current
  calibrated default is `0.75`.

`structural_matching.window`
: Candidate search window, in base pairs. For each LiftOn transcript, target
  genes on the same contig and strand are considered if they overlap the
  projected transcript span expanded by this window.

`structural_matching.topk`
: Maximum number of transcript candidates retained per LiftOn transcript before
  one-to-one greedy assignment.

`structural_matching.min_score`
: Minimum transcript structural score needed for a transcript pair to enter the
  candidate set. Lower values increase sensitivity and review burden; higher
  values increase specificity.

`structural_matching.good_score`
: Score threshold used to label transcript pairs as `good` in
  `lifton.transcript_pairs.tsv`.

`structural_matching.confident_score`
: Score threshold used to label transcript pairs as `confident` in
  `lifton.transcript_pairs.tsv`.

`structural_matching.gene_fraction`
: Gene aggregation threshold. Transcript-pair scores are summed by
  LiftOn-gene/target-gene pair; the best target gene is accepted only if its
  fraction of the old gene's total transcript evidence is at least this value.
  This helps avoid accepting ambiguous split/merged evidence.

`structural_matching.score_weights`
: Weights for transcript structural scoring. The code normalizes the weights,
  so the values do not need to sum to `1`.

`score_weights.span_containment`
: Similarity of the projected transcript span and target transcript span,
  measured by overlap divided by the shorter span. This is the dominant current
  weight because it is robust to exon/CDS representation differences.

`score_weights.query_coverage`
: Fraction of the LiftOn transcript's exon/CDS interval covered by the target
  transcript. This is asymmetric: it asks how much of the projected old
  structure is represented in the target.

`score_weights.intron_sim`
: Exact intron-chain similarity. High value means splice sites agree.

`score_weights.jacc_internal`
: Jaccard similarity of internal exon intervals.

`score_weights.jacc_all`
: Jaccard similarity of all exon/CDS intervals.

`score_weights.exon_count_sim`
: Similarity based on exon count.

`score_weights.boundary_sim`
: Soft similarity of TSS/TES boundaries. Small boundary differences are
  tolerated.

`score_weights.lifton_identity_prior`
: Optional prior from LiftOn `protein_identity` and/or `dna_identity`
  attributes when present. If those attributes are absent, this component is
  `0`.

## Outputs

For `--output-dir out/stable-id-run`, the pipeline writes:

```text
out/stable-id-run/lifton/lifton_feature_types.txt
out/stable-id-run/lifton/projected_ref_on_target.gff3
out/stable-id-run/reports/missing_genes.txt
out/stable-id-run/matching/lifton.transcript_pairs.tsv
out/stable-id-run/matching/lifton.gene_pairs.tsv
out/stable-id-run/matching/lifton.gene_locus_comparison.tsv
out/stable-id-run/score_evidence.tsv
out/stable-id-run/stable_id_decisions.tsv
out/stable-id-run/sql/stable_id_updates.dry_run.sql
```

If `--write-executable-sql` is used, the SQL suffix is `.sql` instead of
`.dry_run.sql`.

### LiftOn Outputs

`lifton/lifton_feature_types.txt`
: Feature-type file passed to LiftOn `-f`.

`lifton/projected_ref_on_target.gff3`
: Reference genes projected onto the target assembly by LiftOn. This is not
  accepted directly as the stable-ID mapping; it is the old-feature query used
  for matching against the target annotation.

`reports/missing_genes.txt`
: Reference genes reported as not projected by LiftOn.

### Structural Matching Outputs

`matching/lifton.transcript_pairs.tsv`
: One-to-one transcript structural matches. Important columns:
  `lifton_tx`, `ref_tx`, `score`, `status`, component scores,
  `lifton_gene`, and `ref_gene`.

`matching/lifton.gene_pairs.tsv`
: Gene-level structural matches aggregated from transcript-pair evidence.
  Important columns: `lifton_gene`, `ref_gene`, `weighted_score`,
  `fraction_of_total`, and `n_transcripts`.

`matching/lifton.gene_locus_comparison.tsv`
: Audit/calibration table. For each LiftOn gene, it records the best target
  gene by whole-gene locus overlap and compares it with accepted structural
  evidence. Important columns include `target_gene_by_locus`, `locus_score`,
  `structure_accepted_target_gene`, `best_tx_structure_target_gene`, and
  `locus_vs_structure`.

### Decision Outputs

`score_evidence.tsv`
: Normalized evidence loaded from the structural pair tables. This is useful
  when checking whether a pair table row was actually usable by the decision
  layer.

`stable_id_decisions.tsv`
: Main review table. Columns:

```text
type                 gene, transcript, or translation
action               mapped, missing, or new
current_stable_id    target stable ID currently present in target GFF3
current_version      current target version
old_stable_id        old/reference stable ID, for mapped or missing rows
old_version          old/reference version
new_stable_id        stable ID to write after mapping
new_version          new version to write
mapping_session_id   supplied mapping session
score                confidence/evidence score
reason               short explanation of the decision path
```

Mapped rows have both `current_stable_id` and `old_stable_id`; the
`new_stable_id` is the old stable ID and `new_version` is `old_version + 1`.

Missing rows have an `old_stable_id` but no `current_stable_id` or
`new_stable_id`.

New rows have a `current_stable_id` and an allocated `new_stable_id` from the
supplied range with version `1`.

`sql/stable_id_updates.dry_run.sql`
: Review SQL. It creates temporary decision/update tables, emits count checks,
  and ends in `ROLLBACK`.

`sql/stable_id_updates.sql`
: Executable SQL, written only with `--write-executable-sql`. It creates backup
  tables, stages decisions, updates stable IDs with a temporary rename step to
  avoid collisions, inserts `stable_id_event`, and ends in `COMMIT`.

## Audit A Completed Run

After a run, write compact audit tables:

```bash
python3 pipelines/stable_id_mapper/audit_stable_id_run.py \
  --run-dir out/stable-id-run \
  --ref-gff ref.gff3 \
  --target-gff tar.gff3 \
  --limit 20
```

If the decision run reused matching outputs from another directory, pass the
locus-comparison file explicitly:

```bash
--locus-comparison previous-run/matching/lifton.gene_locus_comparison.tsv
```

Audit command options:

`--run-dir PATH`
: Completed stable-ID mapper output directory. Required.

`--ref-gff PATH`
: Reference GFF3. Used to add old gene biotypes and loci to the audit tables.

`--target-gff PATH`
: Target GFF3. Used to add target gene biotypes and loci to the audit tables.

`--locus-comparison PATH`
: Gene locus comparison TSV to use. Defaults to
  `RUN_DIR/matching/lifton.gene_locus_comparison.tsv`. Set this explicitly when
  a decision run reused matching outputs from another directory.

`--output-dir PATH`
: Audit output directory. Defaults to `RUN_DIR/reports/stable_id_audit`.

`--limit INT`
: Number of example rows printed per audit bucket.

Audit outputs:

```text
out/stable-id-run/reports/stable_id_audit/missing_genes.tsv
out/stable-id-run/reports/stable_id_audit/coordinate_mapped_genes.tsv
out/stable-id-run/reports/stable_id_audit/new_genes.tsv
```

How to interpret the audit:

- `structural_mapped`: best evidence bucket; usually expected to dominate.
- `coordinate_mapped`: fallback bucket. With the current `0.75` default, these
  should be strong locus-overlap mappings, but still deserve review if counts
  are large.
- `missing` with `reference gene was not projected by LiftOn`: old genes LiftOn
  did not project.
- `missing` with `no_locus_candidate`: projected old genes with no usable
  target gene nearby.
- `missing` with `same, claimed=yes`: the best target agreed by locus and
  structure but was already claimed by another old gene.
- `new` with `target_id_also_in_ref_gff=yes`: target current ID collides with a
  reference stable ID, but this is provenance only. Target IDs are not mapping
  evidence.

## Synthetic Calibration Cases

Generate deterministic small cases with known expected counts:

```bash
python3 pipelines/stable_id_mapper/generate_synthetic_stable_id_cases.py \
  --output-dir pipelines/stable_id_mapper/out/synthetic_cases
```

Current cases:

```text
high_identity              all old genes/transcripts/translations should map
isoform_shift_gene_only    gene locus maps, transcript/translation do not
unrelated_empty_lifton     old features missing; target features new
duplicated_competing_locus two old projections compete for one target locus
split_old_gene             one old gene overlaps two target genes
merged_old_genes           two old genes overlap one merged target gene
strand_mismatch            opposite strands should not map
contig_mismatch            same numeric spans on different contigs should not map
```

## Diagnostic Helpers

`diagnose_lifton_inputs.py`
: Checks reference FASTA/GFF3 and a failed LiftOn intermediate transcript FASTA
  for missing seqids, malformed FASTA, duplicate IDs, out-of-bounds features,
  and transcript child structure problems.

`diagnose_structural_matching_inputs.py`
: Checks whether target exon/CDS `Parent` attributes link to transcript rows by
  exact ID or core stable ID.

`diagnose_stable_id_feature_counts.py`
: Reports the feature universe seen by the stable-ID decision parser.

`pick_real_stable_id_examples.py`
: Selects representative examples from a completed run.

`pick_ftp_stable_id_truth_examples.py`
: Samples examples from old/new Ensembl-style GFF3 releases or direct GFF3
  paths. Use core DB `stable_id_event` as gold standard for complex official
  remaps.

## Compatibility And Prototype Scripts

`main_output_to_stable_id_event_sql.py`
: Compatibility wrapper around the modular SQL/decision code. It starts from an
  existing mapped/projected GFF3 and report rather than running LiftOn and
  structural matching.

`main.py` / `mapper.py`
: Older assembly-to-assembly prototype using gene-sequence alignment. It is
  useful as a prototype but is not the recommended production path.

`lifton_id_mapper.py`
: Standalone structural mapper used internally by `run_stable_id_mapping.py`.
  It can still be run directly for debugging or isolated matching experiments.

## Testing

From the repository root:

```bash
python3 -m pytest pipelines/stable_id_mapper/tests
python3 pipelines/stable_id_mapper/run_stable_id_mapping.py --help
python3 pipelines/stable_id_mapper/audit_stable_id_run.py --help
python3 pipelines/stable_id_mapper/run_lifton_projection.py --help
```

The compatibility SQL wrapper also has a built-in smoke test:

```bash
python3 pipelines/stable_id_mapper/main_output_to_stable_id_event_sql.py --test
```

## Known Limitations

- The current pipeline is single-species. Parallel processing across many
  species should be handled by a future Nextflow wrapper.
- Structural matching uses greedy one-to-one assignment, not a global optimum.
- Biotype is reported in audits but is not used as mapping evidence.
- Target stable IDs are intentionally ignored as evidence, even when they match
  reference IDs. This is to avoid confusion when the IDs assigned don't follow
  stable_id range selection.
- Coordinate fallback currently uses one shared `min_overlap` threshold for
  genes, transcripts, and translations.
- Large missing buckets can be correct when the target annotation knowingly
  contains fewer genes than the reference annotation.

# Stable ID Mapper

Prototype tools for mapping Ensembl stable IDs between gene annotations.

This pipeline lives in `pipelines/stable_id_mapper` and is currently separate
from the packaged `src/python/ensembl/genes` modules.

## Stable ID Mapping

Ensembl stable IDs identify genes, transcripts, exons, and proteins across
releases. When an annotation is updated, IDs for equivalent features should be
carried forward where possible, and new stable IDs should be created for new
features. Ensembl's public stable-ID documentation gives more background on the
ID format and versioning rules:

https://www.ensembl.org/info/genome/stable_ids/index.html

This branch includes two mapping workflows:

1. `main.py` / `mapper.py`: maps IDs from an old/reference FASTA+GFF3 to a
   new/target FASTA+GFF3 by aligning reference gene sequences to the target
   assembly.
2. `lifton_id_mapper.py`: maps a LiftOn GFF3 to a reference/Ensembl GFF3 that
   is already on the same target assembly, using transcript structure and local
   context.

This directory also includes small chicken chromosome 9 example inputs:

```text
ref.fa
ref.gff3
tar.fa
tar.gff3
```

## Requirements

The scripts are written for Python 3.

`lifton_id_mapper.py` uses only the Python standard library.

`main.py` can run with only the Python standard library, but realistic
assembly-to-assembly mapping requires `minimap2`. If `minimap2` is not
available, the code falls back to exact full-gene sequence search, which is
useful as a smoke test but usually too strict for real assembly updates.

Install `minimap2` with one of:

```bash
brew install minimap2
```

```bash
conda install -c bioconda minimap2
```

For repository development tooling, from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Single-Species Stable-ID Pipeline

The main one-species entrypoint is `run_stable_id_mapping.py`. It runs LiftOn,
writes structural matching score tables, loads those scores as mapping
evidence, generates a missing-gene report from the projected GFF3, and then
writes stable-ID decision TSV plus review SQL.

Required manual inputs:

- old/reference FASTA and GFF3
- new/target FASTA and GFF3
- core DB name
- mapping session ID
- gene, transcript, and translation stable-ID ranges
- output directory

Mapping thresholds and score weights live in `stable_id_mapping_rules.json`.
Use `--rules-config path/to/rules.json` to run with a copied/edited rules file.
Explicit CLI flags such as `--match-min-score` and `--min-overlap` override the
rules file for that run.

Example:

```bash
python3 run_stable_id_mapping.py \
  --ref-fasta ref.fa \
  --ref-gff ref.gff3 \
  --target-fasta tar.fa \
  --target-gff tar.gff3 \
  --db-name gallus_gallus_core_test \
  --mapping-session-id 1 \
  --gene-range ENSGALG:90000000000-90000099999 \
  --transcript-range ENSGALT:90000000000-90000099999 \
  --translation-range ENSGALP:90000000000-90000099999 \
  --output-dir out/stable-id-run \
  --lifton-threads 10
```

By default the SQL is dry-run/review SQL ending in `ROLLBACK`. Add
`--write-executable-sql` only after reviewing the TSV and SQL checks.

Expected outputs:

```text
out/stable-id-run/lifton/lifton_feature_types.txt
out/stable-id-run/lifton/projected_ref_on_target.gff3
out/stable-id-run/matching/lifton.transcript_pairs.tsv
out/stable-id-run/matching/lifton.gene_pairs.tsv
out/stable-id-run/matching/lifton.gene_locus_comparison.tsv
out/stable-id-run/reports/missing_genes.txt
out/stable-id-run/score_evidence.tsv
out/stable-id-run/stable_id_decisions.tsv
out/stable-id-run/sql/stable_id_updates.dry_run.sql
```

The default rules file currently controls:

```text
coordinate_overlap.min_overlap
structural_matching.window
structural_matching.topk
structural_matching.min_score
structural_matching.good_score
structural_matching.confident_score
structural_matching.gene_fraction
structural_matching.score_weights
```

The score weights tune how transcript pair confidence is calculated. The
pipeline still only maps genes, transcripts, and translations; exon/CDS rows are
used only as transcript-structure evidence.

Use `--dry-run-lifton-command` to print the LiftOn command without executing the
pipeline. This is useful on systems where LiftOn is provided by a module or
container. The dry run also writes `lifton/lifton_feature_types.txt`, because
LiftOn expects `-f` to point to a file rather than a comma-separated inline list.
The default feature file contains only `gene`, so LiftOn projects whole gene
models and does not try to lift child rows such as exons or UTRs as independent
loci.

After a successful LiftOn run, reuse existing outputs while debugging matching,
scoring, and SQL:

```bash
python3 run_stable_id_mapping.py \
  --ref-fasta ref.fa \
  --ref-gff ref.gff3 \
  --target-fasta tar.fa \
  --target-gff tar.gff3 \
  --db-name gallus_gallus_core_test \
  --mapping-session-id 1 \
  --gene-range ENSGALG:90000000000-90000099999 \
  --transcript-range ENSGALT:90000000000-90000099999 \
  --translation-range ENSGALP:90000000000-90000099999 \
  --output-dir out/stable-id-rerun \
  --existing-lifton-gff out/stable-id-run/lifton/projected_ref_on_target.gff3
```

To skip both LiftOn and structural matching, provide all three reuse files:

```bash
--existing-lifton-gff out/stable-id-run/lifton/projected_ref_on_target.gff3 \
--existing-transcript-pairs out/stable-id-run/matching/lifton.transcript_pairs.tsv \
--existing-gene-pairs out/stable-id-run/matching/lifton.gene_pairs.tsv
```

The LiftOn structural pair tables are the primary evidence for claiming a target
gene or transcript. If no usable pair-table match exists, the decision layer
falls back to coordinate overlap between the projected old feature and the
target GFF3 feature. Projected LiftOn IDs alone are not treated as mapped unless
they are tied to a concrete target feature. `score_evidence.tsv` lists the
normalized evidence loaded from the LiftOn gene and transcript pair tables,
including confidence labels and score source.

`matching/lifton.gene_locus_comparison.tsv` is a calibration report, not a
decision input yet. It compares each LiftOn-projected old gene with the best
nearby target gene by whole-gene span overlap and writes those locus metrics
next to the accepted transcript-structure evidence. This makes it easier to
spot cases where transcript structure is too strict because a different isoform
or CDS/exon representation was favored.

### Synthetic Calibration Cases

The bundled chicken chromosome 9 files are useful real-data examples, but they
are not a formal truth set unless the old and target annotations are known to be
near-equivalent. For behavior with known expectations, generate synthetic cases:

```bash
python3 pipelines/stable_id_mapper/generate_synthetic_stable_id_cases.py \
  --output-dir pipelines/stable_id_mapper/out/synthetic_cases
```

The generator writes small FASTA/GFF3/LiftOn-output triples and prints a command
for each case. The current cases are:

```text
high_identity              all old genes/transcripts/translations should map
isoform_shift_gene_only    gene locus should map, transcript structure should not
unrelated_empty_lifton     no old features should map; target features are new
duplicated_competing_locus two old projections compete for one target locus
split_old_gene             one old gene overlaps two target genes
merged_old_genes           two old genes overlap one merged target gene
strand_mismatch            overlapping spans on opposite strands should not map
contig_mismatch            matching numeric spans on different contigs should not map
```

Each case also writes `case.json` with the expected decision counts. These cases
are meant for regression checks and rule calibration before drawing conclusions
from messy real annotation examples.

For a completed real-data run, pick a small review set from the pipeline's own
outputs with:

```bash
python3 pipelines/stable_id_mapper/pick_real_stable_id_examples.py \
  --run-dir pipelines/stable_id_mapper/out/chr9-locus-calibration
```

This writes `reports/real_example_candidates.tsv` with representative examples
such as structural gene mappings, coordinate-fallback mappings, high locus
overlap with no accepted structure, locus/structure disagreements, and projected
genes with no target locus candidate.

That report is for inspecting a run; it is not an external truth set. To get
real examples from old/new annotation releases, compare GFF3 files from an
Ensembl-style FTP mirror:

```bash
python3 pipelines/stable_id_mapper/pick_ftp_stable_id_truth_examples.py \
  --ftp-root https://ftp.ensembl.org/pub \
  --old-release 114 \
  --new-release 115 \
  --species gallus_gallus \
  --feature-type gene \
  --seqid 9 \
  --output pipelines/stable_id_mapper/out/gallus_gallus_r114_r115_gene_truth.tsv \
  --manifest pipelines/stable_id_mapper/out/gallus_gallus_r114_r115_sources.json
```

The same tool also accepts direct local paths or URLs:

```bash
python3 pipelines/stable_id_mapper/pick_ftp_stable_id_truth_examples.py \
  --old-gff old.gff3.gz \
  --new-gff new.gff3.gz \
  --feature-type gene \
  --output real_gene_truth.tsv
```

It classifies examples as shared stable IDs at the same locus, shared stable IDs
with changed coordinates, version-changed stable IDs, old-only missing
candidates, and new-only new-feature candidates. For complex official remaps
where the old stable ID changes to a different new stable ID, use the core DB
`stable_id_event` table as the gold standard rather than GFF3 alone.

If LiftOn fails while indexing `lifton_output/intermediate_files/transcripts.fa`,
diagnose the reference inputs and the failed intermediate FASTA with:

```bash
python3 pipelines/stable_id_mapper/diagnose_lifton_inputs.py \
  --ref-fasta ref.fa \
  --ref-gff ref.gff3 \
  --lifton-transcripts-fa test_ii/lifton/lifton_output/intermediate_files/transcripts.fa \
  --output-tsv test_ii/lifton_input_diagnostics.tsv
```

This checks for GFF3 seqids missing from the reference FASTA, out-of-bounds
transcript/exon/CDS coordinates, transcript rows with no exon/CDS children,
duplicate feature IDs, and malformed or empty FASTA records in LiftOn's
intermediate transcript FASTA.

## Workflow 1: Assembly-to-Assembly Mapper

Use `main.py` when you have:

- old/reference FASTA
- old/reference GFF3
- new/target FASTA
- new/target GFF3

Show the command-line options:

```bash
python3 main.py --help
```

Run the bundled example:

```bash
mkdir -p out/root-smoke

python3 main.py \
  --ref-fasta ref.fa \
  --ref-gff ref.gff3 \
  --target-fasta tar.fa \
  --target-gff tar.gff3 \
  --output-gff out/root-smoke/mapped.gff3 \
  --report out/root-smoke/report.txt \
  --threads 8 \
  --identity-min 0.80
```

Outputs:

```text
out/root-smoke/mapped.gff3
out/root-smoke/report.txt
```

### Assembly Mapper Algorithm

At a high level, `main.py` calls `mapper.map_ids()` and performs these steps:

1. Load a minimal GFF3 hierarchy of `gene`, `mRNA`/`transcript`, `exon`, and
   `CDS` features.
2. Scan reference and target IDs to learn existing stable-ID prefixes, numeric
   widths, and namespaces such as `gene:` and `transcript:`.
3. Extract each reference gene sequence from the reference FASTA.
4. Align reference gene sequences to the target FASTA with `minimap2`, if it is
   installed.
5. Fall back to exact sequence matching if `minimap2` is unavailable.
6. Project transcript, exon, and CDS coordinates through the gene-level
   alignment CIGAR.
7. For mapped genes, keep the reference stable IDs and increment the `version`
   attribute.
8. For target genes that are not replaced by a mapped reference gene, create
   new stable IDs using target-derived prefixes.
9. Write a mapped GFF3 and a text report listing mapped, missing, and new genes.

Relevant files:

```text
main.py        CLI wrapper
mapper.py      mapping flow, prefix inference, version handling, report writing
aligner.py     minimap2 invocation and exact-match fallback
cigar_map.py   CIGAR interval projection
fasta_io.py    in-memory FASTA reader
gff_io.py      minimal GFF3 parser and writer
id_manager.py  stable-ID parsing and new-ID allocation
synteny.py     placeholder synteny scoring
```

### Assembly Mapper Smoke Test

If `minimap2` is not installed, the bundled example is expected to map very few
genes because it uses exact full-gene sequence matching. A fallback-path smoke
test produced:

```text
Genes in reference: 451
Mapped: 1 (0.1%)
Missing: 450 (53.9%)
New: 384
```

This proves the command runs, but it should not be interpreted as a useful
biological mapping result. Install `minimap2` before evaluating mapping quality.

### Core-DB SQL Generator

`main_output_to_stable_id_event_sql.py` is a compatibility wrapper around the
modular `stable_id_mapping` package. It converts a projected/mapped GFF3 and
mapping report to SQL that can be reviewed and then run against an Ensembl core
DB. It does not connect to MySQL itself; it only writes SQL and optional TSV
files.

Inputs needed:

- old/reference GFF3: the same file passed to `main.py --ref-gff`
- new/target GFF3: the same file passed to `main.py --target-gff`
- mapped GFF3: the file produced by `main.py --output-gff`
- mapping report: the file produced by `main.py --report`; if `--report` is
  not passed, the converter looks for `report.txt` next to `--mapped-gff`
- core DB name, used for the generated SQL `USE` statement
- `mapping_session_id`: the core-DB mapping session ID to write into
  `stable_id_event`
- gene, transcript, and translation stable-ID ranges, passed as
  `PREFIX:START-END`
- output SQL path
- optional output TSV path for review

The generated executable SQL:

1. backs up the core tables it will touch inside the same DB,
2. stages gene/transcript/translation stable-ID decisions,
3. resolves target rows to internal DB IDs before applying updates,
4. updates stable IDs with a two-phase temporary rename to avoid cascading
   `A -> B`, `B -> C` collisions,
5. inserts rows into `stable_id_event`.

It does not update exons.

#### Get the Stable-ID Range

The ID ranges come from `gb1` / `gb_assembly_metadata`. Query them with the GCA
accession for the assembly:

```sql
SELECT
    species_log.gca_accession,
    CONCAT(prefix.prefix, ':', stable.stable_space_start, '-', stable.stable_space_end) AS id_range
FROM stable_space_species_log AS species_log
INNER JOIN stable_space AS stable ON stable.stable_space_id = species_log.stable_space_id
INNER JOIN species_prefix AS prefix ON prefix.lowest_taxon_id = species_log.lowest_taxon_id
WHERE species_log.gca_accession = '<GCA_accession>';
```

Pass the returned ranges to the script as `--gene-range`,
`--transcript-range`, and `--translation-range`. `START` and `END` may include
leading zeroes; the widest width is preserved when generating IDs.

#### Run from Mapped/Projected Output

```bash
python3 main_output_to_stable_id_event_sql.py \
  --ref-gff ref.gff3 \
  --target-gff tar.gff3 \
  --mapped-gff out/root-smoke/mapped.gff3 \
  --report out/root-smoke/report.txt \
  --db-name gallus_gallus_core_test \
  --mapping-session-id 1 \
  --gene-range ENSGALG:90000000000-90000099999 \
  --transcript-range ENSGALT:90000000000-90000099999 \
  --translation-range ENSGALP:90000000000-90000099999 \
  --include-translations \
  --backup-prefix stable_id_mapper_backup_test \
  --output-sql out/root-smoke/core_updates.sql \
  --output-tsv out/root-smoke/core_updates.tsv
```

The script trusts the mapped/projected GFF3 as the source of successfully
mapped old IDs. It classifies the reference and target GFF3s by hierarchy:
`gene` records are genes, `mRNA`/`transcript` children of genes are transcripts,
and `CDS` records are translations. It does not update exons. Mapped stable IDs
always keep the old stable ID and increment the old version by one. New stable
IDs are allocated from the supplied stable-space range with version `1`.
Translation mappings are inferred only when one old and one target translation
sit under a mapped transcript; otherwise target translations receive new IDs
from `--translation-range`. When called through `run_stable_id_mapping.py`,
the generated decision and SQL scores are fed by LiftOn structural score
evidence where available.

Review the TSV before running the SQL against a core database. The SQL also
contains `SELECT` count checks for staged rows versus DB-matched rows.

Add `--dry-run` to generate validation-only SQL:

```bash
python3 main_output_to_stable_id_event_sql.py \
  --ref-gff ref.gff3 \
  --target-gff tar.gff3 \
  --mapped-gff out/root-smoke/mapped.gff3 \
  --report out/root-smoke/report.txt \
  --db-name gallus_gallus_core_test \
  --mapping-session-id 1 \
  --gene-range ENSGALG:90000000000-90000099999 \
  --transcript-range ENSGALT:90000000000-90000099999 \
  --translation-range ENSGALP:90000000000-90000099999 \
  --include-translations \
  --dry-run \
  --output-sql out/root-smoke/core_updates.dry_run.sql \
  --output-tsv out/root-smoke/core_updates.tsv
```

Dry-run SQL only creates temporary tables and emits `SELECT` checks. It does
not create backup tables, update core feature tables, delete previous events, or
insert into `stable_id_event`.

## Workflow 2: LiftOn-to-Reference Mapper

### Run LiftOn Projection

The stable-ID pipeline should run LiftOn itself. The wrapper keeps this as a
separate file-in/file-out step so it can be moved into Nextflow later without
changing the Python decision logic.

```bash
python3 run_lifton_projection.py \
  --ref-gff ref.gff3 \
  --ref-fasta ref.fa \
  --target-fasta tar.fa \
  --output-gff out/lifton/projected_ref_on_target.gff3 \
  --threads 10
```

The wrapper runs a command shaped like:

```bash
lifton -f out/lifton/lifton_feature_types.txt \
  -t 10 \
  -g ref.gff3 \
  -o out/lifton/projected_ref_on_target.gff3 \
  tar.fa \
  ref.fa
```

Use `--dry-run-command` to print the command without executing it. Use
`--feature-types` to override the default parent feature list (`gene`); the
wrapper writes the list to `lifton_feature_types.txt` before calling LiftOn. Use
`--feature-types-file` if you already have a LiftOn feature-type file.

Use `lifton_id_mapper.py` when you have:

- a LiftOn-projected GFF3 as the query
- a reference or Ensembl GFF3 on the same target assembly

Show the command-line options:

```bash
python3 lifton_id_mapper.py --help
```

Run the mapper:

```bash
mkdir -p out/lifton

python3 lifton_id_mapper.py \
  --lifton liftOn.gff3 \
  --reference ensembl_on_target.gff3 \
  --out-prefix out/lifton/mapping \
  --window 100000 \
  --topk 5 \
  --min-score 0.60 \
  --good 0.75 \
  --confident 0.85 \
  --gene-fraction 0.60 \
  --rewrite-reference \
  --rename-mode alias
```

Outputs:

```text
out/lifton/mapping.transcript_pairs.tsv
out/lifton/mapping.gene_pairs.tsv
out/lifton/mapping.mapped_reference.gff3
```

`--rename-mode alias` keeps reference IDs unchanged and adds
`lifton_gene_id` / `lifton_transcript_id` attributes.

`--rename-mode rename` replaces reference IDs with LiftOn IDs and updates
child `Parent` attributes. Use this mode only after reviewing the mapping
tables.

### LiftOn Mapper Algorithm

`lifton_id_mapper.py` performs structure-based transcript matching:

1. Parse both GFF3 files into genes and transcripts.
2. Build an interval index over reference genes by contig and strand.
3. For each LiftOn transcript, find candidate reference transcripts on the same
   contig and strand within `--window`.
4. Score candidate transcript pairs using:
   - intron-chain similarity, weight `0.50`
   - internal-exon Jaccard similarity, weight `0.25`
   - all-exon Jaccard similarity, weight `0.10`
   - exon-count similarity, weight `0.05`
   - TSS/TES boundary similarity, weight `0.05`
   - optional LiftOn `protein_identity` / `dna_identity` prior, weight `0.05`
5. Keep the top `--topk` candidates per LiftOn transcript.
6. Greedily select one-to-one transcript mappings by descending score.
7. Aggregate transcript mappings into gene mappings.
8. Optionally rewrite the reference GFF3 in alias or rename mode.

### LiftOn Mapper Smoke Test

A useful self-test is to map `tar.gff3` against itself:

```bash
mkdir -p out/lifton-self

python3 lifton_id_mapper.py \
  --lifton tar.gff3 \
  --reference tar.gff3 \
  --out-prefix out/lifton-self/self \
  --rewrite-reference \
  --rename-mode alias
```

Expected output for the bundled fixture:

```text
LiftOn transcripts with exons: 1977
Mapped transcripts (one-to-one): 1977
Top score: 0.950, median: 0.950
Gene pairs: 838 (fraction threshold 0.6)
```

The maximum score is `0.950` rather than `1.000` because the optional LiftOn
identity prior is not present in `tar.gff3`, and that prior contributes the
final `0.05` of the score.

## Testing

Basic checks:

```bash
python3 pipelines/stable_id_mapper/main_output_to_stable_id_event_sql.py --test
python3 -m pytest pipelines/stable_id_mapper/tests
python3 main.py --help
python3 run_stable_id_mapping.py --help
python3 run_lifton_projection.py --help
python3 lifton_id_mapper.py --help
python3 -m compileall .
```

If `pytest` is installed, run:

```bash
python3 -m pytest
```

## Known Limitations

The LiftOn mapper is currently the more self-contained and easier-to-test
workflow. The assembly-to-assembly mapper is useful as a prototype but needs
more validation before it should be used for release-quality GFF3 output.

Assembly mapper limitations:

- The assembly mapper GFF3 parser only handles exact feature types `gene`,
  `mRNA`, `transcript`, `exon`, and `CDS`. It does not currently preserve all
  Ensembl noncoding and pseudogene feature types.
- Exon and CDS attributes are not stored by `gff_io.py`, so exon/CDS metadata
  such as `exon_id`, `protein_id`, and phases may not be fully preserved.
- Generated transcript IDs are not always propagated to child exon/CDS `Parent`
  attributes in the output GFF3.
- Output rows are sorted by coordinate and feature type, so child rows can be
  emitted before their parent rows when they share the same start coordinate.
- `synteny_disambiguate()` currently computes neighbor scores but returns the
  input mappings unchanged.
- Without `minimap2`, exact-match fallback is too strict for most assembly
  updates.

LiftOn mapper limitations:

- Both inputs must already be on the same assembly coordinate system.
- Matching is greedy rather than a global optimum assignment.
- Biotype is recorded as metadata but ignored for scoring.
- `--rename-mode rename` is intentionally invasive and should be used only
  after reviewing `*.transcript_pairs.tsv` and `*.gene_pairs.tsv`.

## Suggested Next Work

Good next steps for making this branch production-ready:

1. Add unit tests for ID parsing, new-ID allocation, GFF3 parsing, scoring, and
   parent/child ID consistency.
2. Run real LiftOn smoke tests and tune score thresholds/reporting from those
   outputs.
3. Fix assembly mapper output so generated transcript IDs are consistently used
   by child features.
4. Use LiftOn and `minimap2` in smoke tests and record expected mapping
   summaries for the bundled fixtures.
5. Decide whether these scripts should stay as a standalone pipeline or move
   under `src/python/ensembl/genes/stable_id`.

# Annotation QC dashboard and locus browser

A local, reproducible way to compare gene predictions with a reference annotation using
`annotation-qc pairwise-compare`, see where predictors differ, and inspect individual loci.
Predictors, runs and reference genomes are added through a configuration file; no
application code changes.

```
experiment config ──► benchmark validate / prepare / import / run / independent ──► comparator outputs
                                                                                       │
                                     benchmark build (SQLite dataset) ◄────────────────┘
                                                   │
                                     annotation-qc dashboard (http://127.0.0.1:8765)
```

| layer | code | does |
|---|---|---|
| comparison | `metrics/pairwise`, `runners/pairwise_compare.py` (unchanged) | all classification and metric semantics |
| experiment preparation | `benchmark/` | configuration, recorded input preparation, validation, cached sequential runs, import of existing results, metric derivation, independent cross-check, dataset build |
| presentation | `dashboard/` | read-only local web server and single-page app |

**Why this design.** The dashboard uses only the Python standard library (`http.server`,
`sqlite3`) plus the package's existing dependencies; the page is plain HTML/JS/SVG with no
CDN or hosted service. Gene models are parsed once, with the comparator's own parser and
selection functions, into an indexed SQLite file, so the browser answers bounded queries
in milliseconds instead of re-reading 1.4 GB annotations. Nothing new is added to the core
dependencies; the optional `dashboard` extra only adds PyYAML for YAML configurations
(JSON works without it).

## 1. Install and launch

```bash
pip install -e "ensembl-genes"                # core: enough for JSON configs and the dashboard
pip install -e "ensembl-genes[dashboard]"     # optional: YAML configurations
```

```bash
annotation-qc dashboard --config experiment.json            # http://127.0.0.1:8765/
annotation-qc dashboard --config a.json --config b.json --port 8800
```

The server binds to `127.0.0.1` by default and has no authentication; do not bind it to a
public interface.

## 2. Open the existing human/mouse benchmark

The populated configuration is
`annotations_to_compare/qc_configs/ensembl2026_human_mouse.json`. It imports the
2026‑10‑05 benchmark (24 whole-genome comparisons; the 8 region pilots are recorded
separately) without rerunning anything.

```bash
cd annotations_to_compare
annotation-qc benchmark import   --config qc_configs/ensembl2026_human_mouse.json   # ~1 min: verifies prepared inputs, matches 24 results by checksum
annotation-qc benchmark validate --config qc_configs/ensembl2026_human_mouse.json   # ~1–2 min first time (cached afterwards)
annotation-qc benchmark build    --config qc_configs/ensembl2026_human_mouse.json   # ~6 min first time, ~10 GB peak; seconds when unchanged
annotation-qc dashboard          --config qc_configs/ensembl2026_human_mouse.json
annotation-qc benchmark verify-headline --config qc_configs/ensembl2026_human_mouse.json \
    --headline benchmark_runs/20261005_145546/headline_metrics.tsv                  # 1,368 values, 24 rows
```

Outputs go to `annotations_to_compare/qc_workspace/ensembl2026_human_mouse_ab_initio/`.
That dataset was built with store version 1, before the comparator's one-to-one
matching existed; the dashboard shows its one-to-one metrics as not available and its
earlier best-pair estimate under "Historical estimates" (a banner says so). Under the
current comparator code its results are stale (implementation changed), so rebuilding
it needs `benchmark run`, not only `build`.

The pig experiment (`qc_configs/pig_sscrofa11_ab_initio.json`) imports
`pig_comparison_runs/20261006_195011` (18 results, including supplementary
stop-harmonised runs) the same way: `import`, `validate`, `build`, `dashboard`.

## 3. Workflows

All commands take `--config EXPERIMENT.json`.

| command | purpose |
|---|---|
| `annotation-qc benchmark validate [--prepare]` | check configuration, inputs, formats, sequence names/bounds, scope and the known parser limitations; writes `validation/validation.json`; exit 1 on errors |
| `annotation-qc benchmark prepare [--force]` | recorded preparation (decompress, retype) and scope derivation |
| `annotation-qc benchmark import [--benchmark-dir DIR]` | register completed comparisons from an existing directory (no rerun, no copy) |
| `annotation-qc benchmark status` | status of every (run, policy): complete / stale / failed / not_run / unavailable, with the reason |
| `annotation-qc benchmark run [--runs …] [--policies …] [--dry-run] [--force]` | run missing, stale or failed comparisons one at a time |
| `annotation-qc benchmark independent [--runs …]` | optional independent cross-check (any-pair intron chain, exact-CDS agreement) |
| `annotation-qc benchmark build` | build or update the dashboard dataset (incremental) |
| `annotation-qc benchmark verify-headline --headline FILE` | compare dataset metrics with a `headline_metrics.tsv` (wide) table |
| `annotation-qc dashboard --config …` | serve the dashboard |

A new experiment: `validate --prepare` → `run` → (`independent`) → `build` → `dashboard`.

## 4. Configuration (format version 1)

JSON (or YAML with the extra). Relative paths resolve against the configuration file. A
commented template is at `benchmark/templates/experiment_template.json`; keys starting with
`_` are ignored.

| key | meaning |
|---|---|
| `config_version` | `1` |
| `experiment.id`, `title`, `description`, `synthetic` | identity; `synthetic: true` shows a warning banner everywhere |
| `output_dir` | workspace for every derived file (sources are never written) |
| `evaluation.evaluation_mode` | comparator preset: `all`, `protein_coding`, `cds_only`, `canonical` |
| `evaluation.reference_gene_biotypes`, `reference_transcript_biotypes` | extra filters (empty = none) |
| `evaluation.genome_mismatch` | `error` (default), `warn`, `exclude` |
| `policies.<ID>` | `label`, `description`, `reference_transcript_selection` and `query_transcript_selection` (`all` / `longest_cds` / `canonical`), `optional`, `plots_per_category` |
| `references[]` | `id` (unique), `species`, `assembly`, `assembly_accession`, `annotation_source`, `annotation_release`, `annotation`, `annotation_format`, `genome`, `annotation_preparation`, `genome_preparation`, `scope`, `seqname_map`, `notes` |
| `…preparation` | `steps`: `"decompress"`, `{"name": "retype_transcript_types", "types": [...]}`; `reuse_existing`: adopt an existing prepared file only if re-deriving it from the source gives the same sha256 |
| `references[].scope.type` | `all`, `reference_sequence_regions` (`##sequence-region` names), `sequences` (list), `regions_file` (path) |
| `runs[]` | `run_id` (unique identifier), `reference`, `tool` (grouping and colour), `label` (display name; default `tool`), `role` (`primary` / `supplementary`), `comparator_options`, `tool_version`, `model`, `annotation`, `format`, `preparation`, `provenance`, `caveats`, `aliases`, `color` |
| `runs[].comparator_options` | `add_stop_codon` (default false: original coordinates), `stop_codon_genetic_code` (NCBI table, default 1), `stop_codon_sequence_code` (`{sequence: table or "none"}`); passed to `pairwise-compare` and part of the run's parameters and cache key. Needs a reference genome |
| `policies.<ID>.role` | `primary` (default) or `supplementary` (shown on the policy selector) |
| `import` | `benchmark_dir`, `code_snapshot` (sha256sum-style lists of comparator files used then), `independent_tables` (pattern with `{alias}`, `{policy}`), `context_documents`, `paper_analogous` |

Unknown values stay unknown: leave `tool_version`/`model`/`provenance` null and the
dashboard shows "unknown".

### Adding a predictor

Append one object to `runs` with a new `run_id`, then:

```bash
annotation-qc benchmark run   --config experiment.json     # only the new run's comparisons execute
annotation-qc benchmark build --config experiment.json     # only the new components are built
```

Several runs of one tool (different models, versions, settings) are separate `run_id`s
with the same `tool`; they share the tool's colour, later runs in lighter tints.

### Comparing a new genome or reference

Add a `references[]` entry and its `runs[]`, then `validate --prepare`, `run`, `build`.
Validation reports what the parser cannot read yet (compressed annotations not named
`*.gz`, transcript types it does not list, verbatim duplicated rows, orphan exons, GTF
transcript IDs reused without transcript lines) and whether canonical tags exist. If the
annotation uses other sequence names than the FASTA, give a two-column `seqname_map`.

### Scope and transcript policies

- **Scope** is fixed per reference before comparing tools. `reference_sequence_regions`
  limits both annotations to the sequences the reference annotates; predictions elsewhere
  (alt haplotypes, patches) are counted as "outside scope" in the provenance view rather
  than inflating Novel.
- **Policy A** (`all`/`all`): a reference gene counts once if any eligible isoform is
  reproduced. **B** (`longest_cds` on both sides). **C** (`canonical`): the
  `Ensembl_canonical`-tagged transcript, falling back to the longest CDS per gene without a
  tag. C is *not* labelled MANE: for human Ensembl most, not all, canonical transcripts are
  MANE Select. A reference with no canonical tag at all makes C **unavailable** (the
  comparator would silently fall back to B for every gene); per-gene fallbacks are flagged.

### Supplementary sequence scopes

`references[].analysis_scopes` adds named sequence subsets, reported beside the whole
compared scope (`scope.label` names it, e.g. "Whole assembly"):

```json
"analysis_scopes": [{"id": "chr1_18_X", "label": "Chromosomes 1–18 and X",
                     "sequences": ["1", "2", "…", "18", "X"], "description": "…"}]
```

Define a subset from the assembly (sequence names), never from which sequences carry
predictions: a sequence without predictions may still have been processed. A subset
is applied to reference **and** predicted genes before any value is calculated
(`benchmark/scope.py`): R, Q, locus recovery, per-gene exact CDS, splits, merges and
intron chains from the per-gene outputs on those sequences; one-to-one TP from the
comparator's pairs on those sequences (exact-CDS edges never join sequences, so the
maximum matching decomposes by sequence); unique introns from the stored models. Each
build first repeats the derivation for the whole compared scope and compares it with
the comparator's own values; a difference marks the subset's rows *inconsistent*. A
subset naming a sequence that was not compared is reported unavailable. The whole
compared scope always shows the comparator's values.

### Display order and companion runs

`display.pinned_runs` lists runs shown first, in order; the other primary runs follow
in descending `display.order_metric` (default one-to-one exact-CDS recall) for the
selected policy and scope. A supplementary run with `companion_of: <run_id>` (e.g. the
stop-harmonised run of the same file) is listed directly after that run and is never
ranked. Without a `display` block the configured order is kept. Presentation only:
none of this enters a cache key.

### Original and stop-harmonised runs of one file

Configure them as two runs with the same `annotation` and different
`comparator_options` (and, for the harmonised one, a `label` such as "Tool — stop
harmonised" and `role: supplementary`). They have different parameters, commands and
cache keys; `import` assigns each manifest to the run whose stop-codon settings it
records (a manifest without the option compared original coordinates). For a
harmonised run the locus tracks show the comparator's own evaluated annotation
(`query_evaluated.stop_harmonised.gtf` of a complete result), with the per-transcript
audit in the tooltips, gene detail and a TSV download. Without that output the
tracks are reported unavailable; the original file is never shown in their place.

### Import versus run

`import` reuses results whose rebuilt cache key equals the key the configuration and the
current checkout would produce; it needs a `code_snapshot` because the benchmark ran on
uncommitted code (a Git commit alone does not identify it). Without a matching snapshot,
imported results are shown as **stale** and `run` recomputes them. Region-restricted runs
are kept as pilots and never mixed with whole-genome results.

## 5. Caching, provenance and resuming

A comparison is reused only when its **cache key** is unchanged. The key hashes:

- **inputs**: sha256 of the prepared query, reference annotation, genome, scope content and
  seqname map. Prepared files carry the source checksum and recipe.
- **parameters**: every effective `pairwise-compare` option.
- **implementation**: a content fingerprint of the comparator source files, plus the
  python, pandas and pyranges1 versions. The file list is derived, not maintained by
  hand: every in-package module that `runners/pairwise_compare.py` imports, directly or
  transitively (parsed from the sources; docstring-only package `__init__` files are
  left out). A new comparator module is covered as soon as it is imported. Uncommitted
  edits change the fingerprint. A code snapshot that lacks one of these files cannot
  verify a result, which is then imported as stale with the missing files named.

Consequences:

- Adding a run touches nothing else.
- Changing a policy invalidates only that policy.
- Changing a reference, the scope or the comparator code invalidates the affected runs.

The dataset build has its own component keys, which include a hash of `store.py`
and `derive.py`, so a change in a derivation rebuilds the affected rows.

Every attempt is appended to `comparisons/<run_id>/<policy>.attempts.jsonl`, recording the
argv, environment, exit status, wall time and per-process peak RSS. Output is written to
`comparisons/_partial/` and promoted only after exit 0 with every expected file present. A
failed run therefore never counts as a result and is retried by the next `run`.

```
<output_dir>/
  cache/                     checksum cache, annotation scans, FASTA lengths
  prepared/<reference_id>/   prepared annotation/genome (or verified reused files), scope.regions, *.preparation.json
  prepared/runs/<run_id>/    prepared queries (only if a run has preparation steps)
  validation/validation.json
  comparisons/<run_id>/<policy>/   comparator outputs + run_record.json
  imported/import_report.json      matched, pilots, unmatched, conflicts
  independent/<run_id>/            independent tables (computed) or *.import.json pointers
  dashboard/dashboard.sqlite
```

## 6. Using the dashboard

- **Overview** has three views. The analysis context (transcript policy, sequence
  scope, R and the stop-codon treatment) is shown above them; the scope selector is in
  the header.
  - **Coverage & exact CDS** (default): one compact row per run — loci recovered
    (same-strand CDS overlap), exact CDS recall and precision (one-to-one
    coordinate-exact pairs over R and Q) as percentages, and one-to-one exact-CDS F1 on
    a 0–1 scale. Counts and denominators are in the tooltips and the expandable
    *Counts and denominators* table. A built-in F1 that a dataset lacks is "not
    available" with its reason; the historical best-pair estimate is never shown there.
  - **Structural errors**: CDS intron-chain identity, unique reference / predicted CDS
    intron agreement, split and merge flags, each with its denominator in the header;
    then two exploratory plots (one reference transcript per gene: views B or C, chosen
    beside the plots):
    - *unsupported reference CDS splice sites by position* along the reference CDS
      (transcription direction; 0 = coding start, 1 = coding end; genomic CDS span or
      coding-sequence coordinate; donors and acceptors in separate panels; optionally
      leaving out introns < 60 bp);
    - *splice recovery by reference CDS intron length* (fixed bins, the same for every
      run; complete intron, donor or acceptor), with the number of reference introns
      per bin.

    Support is exact: the same sequence, strand and coordinate in any evaluated
    predicted transcript's CDS introns (the comparator's unique-intron rule). Sites of
    genes whose locus was not recovered are counted separately as missing-locus sites
    and excluded from the rates; empty bins are "not available". Binned counts download
    as TSV. These are associations, not causes.
  - **Detailed metrics**: every metric as before (whole compared scope). The CDS headline table shows every
  value as a percentage with *numerator / denominator* and its source. Failed, not-run,
  stale and unavailable results are labelled as such, and `n/a` (with the reason) is never
  zero. Further sections cover:
  - the policy-effect chart (A/B/C for the same reference);
  - CDS outcome composition, single vs multi-segment genes and partial buckets;
  - missed / Novel / split / merge and Novel-context categories;
  - UTR-inclusive metrics, labelled as not for ranking CDS-only predictors.

  Each table can be downloaded as TSV.
- **Genes & loci.** Search gene IDs, file IDs, names and transcript IDs, and filter by
  sequence, CDS structure, per-run outcome, splits, or a two-run contrast ("exact in run X,
  not exact in run Y"). Click a row to open the locus:
  - shared axis; UTR exons thin and outlined, CDS thick and filled, introns with strand chevrons;
  - evaluated vs non-evaluated isoforms; canonical ★ and compared-pair ◆ markers per run;
  - non-coding or pseudogene context hatched.

  The explanation table gives each run's comparator class, the compared isoform and query
  partner, the best-pair vs any-pair intron chain, and the CDS boundary or splice
  differences recomputed from the stored models. Windows are capped at 2 Mb and 60
  transcripts per track (with a note), and the gene table can be downloaded as TSV.
- **Predictions.** Per run: Novel context, merges, strand mismatch, exactness;
  namespaced IDs show the file ID.
- **Validation & provenance.** Validation issues, probed parser limitations,
  preparation records and checksums, scope, run provenance and caveats, outside-scope
  counts, per-policy run records (source, independent verification, runtime, memory),
  and the import report including pilots.
- **Metric dictionary** and **Context** (paper-analogous tables and documents, marked
  as supporting context, not a reproduction).

Deep links keep the state: `#exp=…&ref=…&policy=…&scope=…&view=coverage|structure|detail`
(plus `plot`, `pos`, `posmin`, `len` for the structural plots).

## 7. Metric interpretation (summary; the full dictionary is in the app)

- Lead with CDS metrics. Exon or UTR agreement is near zero by construction for CDS-only
  predictors.
- Historical `Exact_Match` is *structural*: same intron chain with ≥ 0.8 reciprocal
  overlap, and terminal coordinates may differ. Coordinate-exact CDS is the exactness
  metric.
- Gene-span locus detection pairs genes whose spans overlap, even without shared CDS.
  Same-strand CDS-overlap recovery is the coding measure.
- Intron chain (best pair) is scored on the comparator's single best pair. The independent
  any-pair value can be slightly higher.
- Intron-set metrics count unique introns, and their denominators grow with the number of
  evaluated isoforms. They are not per-gene metrics.
- Query intron support counts query introns outside the evaluated reference subset as
  unsupported.
- One-to-one exact-CDS TP, precision, recall and F1 are the comparator's own
  (`cds_exact_one_to_one`: maximum bipartite matching of genes sharing an identical
  CDS, each gene used once). Precision = TP/Q, recall = TP/R, F1 = 2TP/(R+Q); the table
  shows the numerator and denominator. The exported pairs
  (`cds_exact_one_to_one_pairs.tsv`) are stored, checked against TP and downloadable.
  The per-gene coordinate-exact counts can exceed TP when exact partners are shared.
- Results that predate the comparator's matching show these metrics as not available.
  The earlier best-pair estimate (min of unique genes among per-gene exact pairs) is
  kept as "historical" metrics and is never used in their place.
- Status is explicit: `not_available` (input or field absent), `undefined` (denominator
  0), `not_applicable` (e.g. stop-codon counts of an original-coordinate run),
  `inconsistent` (two outputs of one result disagree; shown with the reason). None is
  zero.
- Every value is for one reference, scope, policy and prediction run. The dashboard never
  pools references or policies and shows no combined "best tool" score.
- The benchmark's `headline_metrics.tsv` and `headline_metrics_long.tsv` hold the same
  data. The dataset is derived from comparator outputs; only the wide table is used, and
  only as a verification target.

## 8. What needs more than summary files

| capability | needs |
|---|---|
| overview metrics from the comparator summary | `comparison_summary.json` |
| one-to-one exact-CDS TP/P/R/F1, CDS-overlap locus recovery, stop-codon counts | `comparison_summary.json` (current comparator) |
| one-to-one pair table and download | `cds_exact_one_to_one_pairs.tsv` |
| per-gene table, partial buckets, historical estimate, Novel lists | `comparison_details.tsv` (+ `consensus_transcript_labels.tsv`) |
| stop-harmonised locus tracks and audit | `query_evaluated.stop_harmonised.gtf`, `stop_codon_harmonisation.tsv` |
| single/multi-segment split, evaluated isoforms, locus browser (reference side) | the prepared reference annotation |
| prediction tracks in the locus browser | the prediction files |
| any-pair intron chain | independent tables (`benchmark independent` or imported) |
| verified reuse of old results | a code snapshot of the comparator files used then |

Missing pieces are reported as `not_available` or "unavailable" with the reason; no gene-level
detail is invented.

## 9. Troubleshooting

- `UnicodeDecodeError` when parsing: a gzip/BGZF annotation not named `*.gz`. Add a
  `decompress` step.
- `Exon/CDS records whose Parent is not a transcript record`: a transcript type the parser
  does not list (add `retype_transcript_types`) or a genuinely orphaned exon (fix the file).
- Imported results show **stale**: run `benchmark status`, which lists the changed key part
  (inputs / parameters / implementation). Supply `import.code_snapshot` for old uncommitted
  code.
- Locus view says models unavailable: the reference was not prepared in this workspace, or
  the prediction file cannot be parsed (shown in Validation & provenance).
- Memory: a human whole-genome comparison peaks at about 7 GB and the human dataset build at
  about 10 GB. Run one at a time (the default).
- Port in use: `--port 8800`.

## 10. Synthetic demonstration

```bash
python -m ensembl.genes.annotation_qc.benchmark.examples.synthetic qc_synthetic_demo_NOT_REAL_DATA
cd qc_synthetic_demo_NOT_REAL_DATA
annotation-qc benchmark validate --config synthetic_demo.json --prepare   # reports the deliberately broken run (exit 1)
annotation-qc benchmark run --config synthetic_demo.json
annotation-qc benchmark independent --config synthetic_demo.json
annotation-qc benchmark build --config synthetic_demo.json
annotation-qc dashboard --config synthetic_demo.json
```

Every identifier starts with `SYNTH` and the experiment is flagged synthetic. It covers:

- two runs of one tool, and a second genome without canonical tags;
- a Tiberius-style file reusing gene IDs on two chromosomes;
- a broken run that fails and is recorded;
- an out-of-scope prediction and an Ensembl TEC transcript type;
- every outcome class.

`tests/annotation_qc/test_benchmark_dashboard.py` uses it to test configuration-only
extension, cache reuse and invalidation, import matching, collision handling, missing and
failed results, and gene-table ↔ locus linkage.

## 11. Known limitations

- The comparator parser limitations listed in the validation view are reported, not fixed.
  The dashboard never repairs biological structures. Preparation only decompresses or
  retypes transcript types, with an audit record.
- Validation's GTF transcript-reuse check is heuristic: strand or gene conflicts, and
  inferred transcripts longer than 2 Mb.
- The locus view shows at most 60 transcripts per track and windows of up to 2 Mb.
- Peak memory is the process maximum RSS for computed runs. It is not recorded for
  imported results.
- Paper-analogous tables are imported context; they are not recomputed for new
  experiments.

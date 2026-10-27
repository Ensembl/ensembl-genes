# Annotation QC dashboard and locus browser

A local, offline way to compare gene predictions with a reference annotation using the
results of `annotation-qc pairwise-compare`, see where prediction runs differ, and
inspect individual loci. Experiments, references, prediction runs (including several
runs of one tool) and supplementary variants are added through configuration only.

**One experiment per species.** Each configuration served by `annotation-qc dashboard`
has exactly one reference (one species, one assembly, one reference annotation) and any
number of prediction runs; the masthead selector lists the experiments by species and
every tab, table and plot follows it. Species may have different prediction tools: the
runs, metrics and available scopes come from each configuration and its comparator
outputs. A configuration with several references is refused with a message asking to
split it (one configuration per reference).

```
experiment config ──► benchmark validate ──► benchmark run (pairwise-compare)  ─┐
                                        └──► benchmark import (existing results) ┴─► comparator outputs (contract 1.0)
                                                                                          │
                                 benchmark build (SQLite index, exact models) ◄──────────┘
                                                   │
                                 annotation-qc dashboard (http://127.0.0.1:8765)
```

| layer | code | does |
|---|---|---|
| comparison | `metrics/pairwise`, `runners/pairwise_compare.py`, `reports/` | every classification, metric, scope, matching and splice table; writes the versioned output contract (package README, "Output contract") |
| orchestration | `benchmark/` | configuration, validation, cached sequential runs, import of existing results, implementation provenance, legacy import, dataset build |
| presentation | `dashboard/` | read-only local web server and single-page app |

**The dashboard computes no metric.** It reads `metrics.tsv` rows verbatim (value,
numerator, denominator, status, reason, populations, unit), the per-gene and pair tables,
the splice tables and the evaluated-transcript export; it indexes, filters and formats them.
Scope values are the comparator's own scope rows (never sums of `per_sequence.tsv`).
Earlier dashboard-side derivations (outcome categories, partial buckets, scoped metrics,
splice bins, CDS differences, transcript re-selection) were removed: the comparator
computes them.

Only the Python standard library (`http.server`, `sqlite3`) and the package's existing
dependencies are used; the page is plain HTML/JS/SVG with no CDN, web font or hosted
service. The optional `dashboard` extra only adds PyYAML for YAML configurations.

## 1. Install, start and restart

The package needs Python 3.12 or newer (`requires-python = ">=3.12"`; tested with 3.12
and 3.13). Name the interpreter explicitly: a bare `python`/`python3` may be an older
system Python (macOS ships 3.9), which pip rejects. A normal (non-editable) install is
enough, and JSON configurations need no extra:

```bash
python3.12 -m venv ~/venvs/annotation-qc          # or python3.13
~/venvs/annotation-qc/bin/pip install /path/to/ensembl-genes          # or a built wheel
~/venvs/annotation-qc/bin/pip install "/path/to/ensembl-genes[dashboard]"   # only for YAML configs
~/venvs/annotation-qc/bin/annotation-qc dashboard --help              # check the install
```

To run the test suite against that install, add `pytest` and run it from a copy of
`tests/` (outside the source tree, so the installed package is the one imported):
`~/venvs/annotation-qc/bin/pip install "pytest>=7,<9"`, then
`~/venvs/annotation-qc/bin/python -m pytest <copy>/annotation_qc`.

Serve one or several experiments (each needs a built dataset, section 3):

```bash
cd annotations_to_compare
~/venvs/annotation-qc/bin/annotation-qc dashboard --config qc_configs/pig_contract_v1.json
~/venvs/annotation-qc/bin/annotation-qc dashboard \
    --config qc_configs/mammal_dashboard/cat.json --config qc_configs/mammal_dashboard/dog.json \
    --config qc_configs/mammal_dashboard/rat.json --config qc_configs/mammal_dashboard/pig.json \
    --config qc_configs/mammal_dashboard/human.json --config qc_configs/mammal_dashboard/mouse.json \
    --port 8765
```

Open `http://127.0.0.1:8765/`. The server binds to `127.0.0.1` and has no
authentication; do not bind it to a public interface. **Restart from a normal terminal:**
stop it with Ctrl+C (or `kill <pid>`; `lsof -i :8765` shows the process) and run the same
command again. Nothing is lost: the datasets are files in each workspace. After a
`benchmark build`, restart the server so it reads the new dataset. Use `--port` when the
port is taken.

## 2. Configure an experiment (format version 2)

JSON (or YAML with the extra). Relative paths resolve against the configuration file.
Keep source files and folder names as they are (paths may use any spelling, e.g.
`fungi/z_triciti/`); the comparator reads gzip/BGZF annotations and FASTA directly.
A commented template is at `benchmark/templates/experiment_template.json`. Version-1
files still load unchanged.

| key | meaning |
|---|---|
| `config_version` | `2` (or `1`) |
| `experiment.id`, `title`, `description`, `synthetic` | identity; `synthetic: true` shows a warning banner |
| `output_dir` | workspace for every derived file (sources are never written) |
| `evaluation` | `evaluation_mode`, `reference_gene_biotypes`, `reference_transcript_biotypes`, `genome_mismatch` (comparator options) |
| `policies.<ID>` | `label`, `description`, `reference_transcript_selection`, `query_transcript_selection` (`all` / `longest_cds` / `canonical`) |
| `references[]` | `id`, `species`, `assembly`, `assembly_accession`, `annotation_source`, `annotation_release`, `annotation`, `genome`, `seqname_map` (applied to both files), `scope` (`type`, `label`), `analysis_scopes` (`id`, `label`, `sequences`, `description`; passed to the comparator as `--analysis-scope`), `caveats`, `notes` |
| `runs[]` | `run_id` (unique, stable), `reference`, `tool` (grouping and colour), `run_label` (e.g. "GMB run 2"), `label` (display text; default tool + run label), `tool_version`, `model`, `annotation`, `format`, `role` (`primary` / `supplementary`), `companion_of`, `comparator_options`, `display_role`, `source` (`original_path`, `sha256`, `recorded_in`), `legacy_ids`, `input_audit`, `provenance`, `caveats` |
| `runs[].display_role` | optional, presentation only: `focus` (the dataset of interest) or `primary_comparator`. Listed first (focus, then primary comparator), reserved colours (vermilion, blue), a text tag, their own marker shapes (filled circle, square; other runs open circles), heavier lines, drawn on top of the other runs. Never part of a comparison, cache key or metric; a species without such a run simply shows none |
| `runs[].comparator_options` | `add_stop_codon`, `stop_codon_genetic_code`, `stop_codon_sequence_code` (`{sequence: table}`); part of the run's parameters |
| `display` | `preferred_tools` (listed first, in order), `order_metric` (default `exact_cds_one_to_one_recall`), `pinned_runs` |
| `provenance.registries` | extra implementation registries (section 5) |
| `import` | `results_dirs` (headline results), `pilot_dirs` (always pilots), `independent_tables`, `context_documents`, `paper_analogous`; version-1 `benchmark_dir` and `code_snapshot` still work |

Validation of the configuration (`load_experiment`) reports every problem at once:
duplicate `run_id`s, two runs comparing the same file with the same options (their results
could not be told apart), a `companion_of` that names another file or the same options,
unknown tools in `preferred_tools`, malformed checksums, invalid policies, scopes or stop
options. `benchmark validate` checks the files: missing paths, `source.sha256` against the
annotation, `input_audit` against the annotation, formats, sequence names and bounds
against the FASTA (after the seqname map), canonical tags.

**A new prediction run** (another tool, or another run of the same tool): add one object
to `runs` with a new `run_id` (and a `run_label` for a further run of a tool). **A new
reference or experiment:** add a `references[]` entry and its runs, or write a new
configuration with its own `output_dir`. Only the new comparisons run; complete results
are reused (section 4).

**Supplementary variants.** A stop-harmonised variant is a second run on the same
`annotation` with `comparator_options.add_stop_codon`, `role: supplementary` and
`companion_of: <original run>`. It is listed beside its original run, labelled
"supplementary · stop-harmonised", never ranked, and can be hidden.

**A primary run with comparator preprocessing.** When only the corrected comparison of a
file should be shown (e.g. a prediction file whose CDS excludes the stop codon), configure
one primary run with `comparator_options.add_stop_codon` and no original-coordinate run.
It appears under its tool name like any other run; the option is listed in its tooltip,
in Provenance & validation, and the added stop codons are marked in the locus view.

## 3. Workflow

| step | command (all take `--config EXPERIMENT.json`) |
|---|---|
| 1. configure | write the configuration (section 2) |
| 2. validate inputs | `annotation-qc benchmark validate` (writes `validation/validation.json`; exit 1 on errors) |
| 3a. generate comparator outputs | `annotation-qc benchmark run [--runs …] [--policies …] [--dry-run]` (one comparison at a time) |
| 3b. or register existing outputs | `annotation-qc benchmark import` (no rerun, no copy) |
| check | `annotation-qc benchmark status`, `annotation-qc benchmark provenance` (read only: they never prepare inputs or write records, and stop with a message when a derived input is not prepared yet) |
| 4. build the dataset | `annotation-qc benchmark build` (incremental; indexes existing comparator outputs into `dashboard/dashboard.sqlite`, never reruns a comparison) |
| 5. serve | `annotation-qc dashboard --config A.json [--config B.json …] [--port 8765]` |
| 6. restart | Ctrl+C, then the same `dashboard` command (section 1) |

`run` passes the analysis scopes (`--analysis-scope`), the run's labels and source
provenance (`--provenance-json`, copied verbatim into the manifest) and
`--export-evaluated-models ids`, writes into `comparisons/_partial/` and promotes output
only after exit 0 with every expected contract file present.

## 4. Result status and reuse

A result is assessed the same way whether it was computed here or imported
(`benchmark/results.py`):

| status | meaning | shown |
|---|---|---|
| `complete` | contract result; inputs, parameters and definition version match the configuration; every output file matches the SHA-256 in its manifest; implementation verified or identical to this checkout | yes |
| `unverified` | as complete, but the implementation cannot be identified (section 5) | yes, flagged "unverified implementation"; `run` recomputes it |
| `legacy` | written before contract 1.0, matching inputs and parameters | yes, through the legacy importer (section 8) |
| `stale` | anything else (the reason names it: changed input, parameter or definition; a file that differs from its manifest checksum; a missing output) | no |
| `failed`, `not_run`, `unavailable` | last attempt failed; no result; policy cannot run (e.g. no canonical tags) | no |

The cache key covers inputs (checksums), comparator parameters (including stop options and
analysis scopes) and the **definition** (contract id, major version, definition version),
not the code fingerprint: code identity is assessed separately, so a verified historical
result stays valid when only dashboard or display code changes. Imports never overwrite a
computed result; `run` never overwrites an imported one.

Pilots: directories listed in `import.pilot_dirs`, and any result restricted with
`--region` or a different regions file, are listed in the import report and never shown as
headline results.

## 5. Implementation provenance

Each manifest records `implementation`: a fingerprint **scheme**, the fingerprint, the
SHA-256 of every covered file and the pandas/pyranges1/Python versions. A build is
identified by the pair (scheme, fingerprint); fingerprints of different schemes are never
compared.

| identity | rule |
|---|---|
| verified implementation `<build>` | (scheme, fingerprint) is in a registry, the fingerprint recomputes from the listed hashes, and the hashes equal the registered build's. The result is attributed to that historical build. |
| current code (this checkout) | scheme 2, not registered, every listed file identical to this checkout |
| unverified implementation | no fingerprint, unsupported scheme, fingerprint inconsistent with its file list, unknown pair, or files that differ from the registered build |

The packaged registry `benchmark/implementation_registry.json` holds:

- **d445054** (scheme `annotation_qc_loaded_modules_by_name/1`, the manifests without a
  `scheme` key; fingerprint `5c2c7089…`): the build behind every result of
  `contract_runs/20261012_151131_comparator_contract_v1`. Evidence: every file hash equals
  the git objects of d445054; the wheel used for the runs (sha256 `acf021d2…`) is
  file-identical to d445054; scheme 1 does not cover the runner, whose identity rests on
  that wheel check.
- **f0f838d** (scheme `annotation_qc_loaded_files_by_path/2`, fingerprint `c4de384a…`,
  29 files including the runner, as loaded by the `annotation-qc` entry point).

Add a build with
`annotation-qc benchmark register-build --manifest M.json … --repo REPO --commit SHA --build-id ID --output registry.json`
(every file hash is checked against `git show SHA:<file>`), then list the file under
`provenance.registries`. `benchmark provenance` prints the identity of every result.
Historical manifests are never rewritten.

## 6. Exact evaluated models (locus tracks)

`benchmark build` rebuilds the models each comparison evaluated, without selecting
transcripts again (`benchmark/models.py`, the comparator's documented procedure):

1. the input recorded in the manifest must exist and match its SHA-256 (and the export's
   `input_sha256`); otherwise the track is unavailable, never another file;
2. it is parsed with the comparator parser and the recorded seqname map;
3. the transcripts of `evaluated_transcripts.tsv.gz` are kept with their per-policy status
   (`evaluated`, or `excluded` with `exclusion_reason`), `selected_by`, original and
   comparison IDs;
4. transcripts marked `stop_codon_added` get the audited terminal CDS from
   `stop_codon_harmonisation.tsv` (missing audit = unavailable);
5. exon count, CDS segment count and CDS length of every evaluated transcript must equal
   the export, and every policy's export must reproduce the same coordinates.

The locus view draws evaluated isoforms solid and excluded ones faded and dashed (with the
reason), marks the compared isoform (◆), canonical tags (★) and audited stop extensions
(red tick). Pre-contract results have no export, so their tracks are unavailable.

## 7. Using the dashboard

The masthead selects the experiment (species; its reference is fixed); the toolbar selects
the transcript policy (default **A**, any-isoform recovery; B and C one transcript per
gene), the sequence scope (whole assembly, or the analysis scopes the comparator reported,
e.g. pig chromosomes 1–18 + X, human/mouse assembled chromosomes) and, when a
configuration has supplementary variants, whether they are shown. A context bar repeats
the reference (assembly, accession, annotation release), policy, scope (with R),
comparator definition version and implementation identity.

- **Headline CDS metrics** (default): one row per run — loci recovered (same-strand CDS
  overlap), one-to-one coordinate-exact CDS recall and precision, one-to-one exact-CDS F1;
  below it, *reference genes by CDS outcome*: one bar per run from the comparator's
  `cds_outcome_*` rows (exact, ends ≠, splice ≠, < 0.8, span, strand, missed), which
  partition R, so exact, partial and missed genes are compared side by side.
  Tooltips (also on keyboard focus) give numerator, denominator, populations, status and
  reason. *Counts and denominators* separates the one-to-one matching (R, Q, TP, edges)
  from the per-gene best-pair coordinate-exact counts, which are a different measure.
  Limitations relevant to the rows are listed under the table and marked ⓘ in the row.
- **Structural agreement**: CDS intron chain (best pair and any pair), unique intron
  recovery and support, splits and merges; splice-site position and intron-length plots
  from the comparator's tables. Under policy A they are not available (the manifest
  records `requires_one_reference_transcript_per_gene`); the page says so and offers B/C.
  Splice-site positions have two coordinates, both computed by the comparator in the
  direction of transcription (0 = coding start, 1 = coding end): **genomic CDS**
  (genomic bases from the coding start to the junction over the genomic CDS span, introns
  included) and **coding sequence** (CDS bases upstream of the junction over the CDS
  length; a donor and the acceptor of the same intron share a position).
- **All metrics**: every contract metric by group, with status and reason; policy-effect
  chart; downloads of the stored metric rows (verbatim), one-to-one pairs, per-sequence
  drill-down and stop-codon audits.
- **Genes & loci**: search, sequence, CDS structure, per-run outcome, split, one-to-one
  pair and two-run contrast filters; the locus view with exact models and a per-run
  explanation table from the comparator's per-gene columns (CDS start/stop differences,
  missing/extra introns, best-pair and any-pair chains).
- **Predictions**: per run, Novel context, merges, best-pair exactness and one-to-one pairs.
- **Provenance & validation**: registry builds and evidence, per-result implementation
  identity and output-checksum check, manifest provenance (prediction id, predictor, run
  label, original source path), input audits (processing scope), import report with
  pilots, validation issues.

Run order: runs with a `display_role` (focus, then primary comparator), then
`display.pinned_runs` and `display.preferred_tools` in order, then the other runs by
descending `display.order_metric` (default one-to-one exact-CDS recall) for the selected
policy and scope; a supplementary variant follows its original. Deep links keep the state
(`#exp=…&policy=…&scope=…&view=…&tab=…`, plus `pos`, `posmin`, `len` for the splice plots). Every control is reachable by keyboard
(Tab; arrow keys in tab lists and segmented controls; Enter opens a gene row).

## 8. Missing values and older results

- Comparator statuses are shown as written: `undefined` (`zero_denominator`),
  `not_applicable` (`transformation_not_applied`, `selection_not_canonical`,
  `requires_one_reference_transcript_per_gene`). Dashboard-side: `not available` (no row),
  and the result statuses of section 4. None is displayed as zero.
- **Legacy importer** (`benchmark/legacy.py`): pre-contract results (e.g. the human/mouse
  benchmark of 2026-10-05) map only summary fields whose definition is identical to a
  contract metric; every other contract metric is `not_available`
  (`pre_contract_result`). The earlier best-pair estimate of one-to-one TP is kept under
  `historical_*` ids and never fills the one-to-one fields. `benchmark verify-headline`
  checks a legacy dataset against a `headline_metrics.tsv` on the mapped columns only.
- Datasets built with store version 1–3 are refused with a message; build into a new
  workspace (they used dashboard-side derivations).

## 9. Limitations shown beside results

- Transcripts outside their own gene record (`gene_span_checks.query.transcripts_outside`
  in the summary) make span-paired per-gene exact counts potentially lower than an
  independent identical-CDS count; one-to-one values are unaffected (e.g. Z. tritici GMB
  run 2).
- Under policy A a stop-harmonised variant's best-pair intron-chain count can differ
  slightly from the original's; any-pair chains, unique introns, R and Q do not.
- Input audits (`runs[].input_audit`) report how many sequences carry predictions and
  whether the predictor recorded which sequences it processed; reference genes on
  sequences without prediction records count as missed, which does not establish that
  those sequences were not processed.
- Reference and run `caveats` from the configuration (e.g. F. graminearum: no
  mitochondrial reference genes, so mitochondrial predictions can only be Novel; GMB
  mitochondrial genetic-code and genome-FASTA provenance limitations).

## 10. Workspace layout

```
<output_dir>/
  cache/checksums.json             sha256 cache keyed by path, size, mtime
  prepared/<reference_id>/         preparation records (sources used as-is unless steps are configured)
  validation/validation.json
  comparisons/<run_id>/<policy>/   computed results + run_record.json (imported: run_record.json points elsewhere)
  comparisons/<run_id>/<policy>.attempts.jsonl
  imported/import_report.json      matched, pilots, unmatched, conflicts
  dashboard/dashboard.sqlite       the index the server reads (store version 4)
```

## 11. Troubleshooting

- A result shows **stale**: `benchmark status` names the changed part, or the output file
  whose checksum differs from its manifest.
- **unverified implementation**: `benchmark provenance` gives the reason; register the
  build if it is known (section 5) or rerun with `benchmark run`.
- Locus tracks unavailable: the reason is shown in the track and in Provenance (input
  moved or changed, audit missing, pre-contract result).
- `store version … is not supported`: build into a new `output_dir`.
- Memory: a human whole-genome comparison peaks at about 7 GB; run one at a time.

## 12. Synthetic demonstration

```bash
python -m ensembl.genes.annotation_qc.benchmark.examples.synthetic qc_synthetic_demo_NOT_REAL_DATA
cd qc_synthetic_demo_NOT_REAL_DATA
annotation-qc benchmark validate --config synthetic_species_A.json --prepare   # reports the deliberately broken run (exit 1)
annotation-qc benchmark run --config synthetic_species_A.json                  # ToolZ fails on purpose
annotation-qc benchmark build --config synthetic_species_A.json
annotation-qc dashboard --config synthetic_species_A.json
```

`synthetic_species_A.json` is the one-species experiment the dashboard serves (reference
SYNTH_A, display roles on ToolX model 1 and ToolY). `synthetic_demo.json` has two
references; the tests use it, and `annotation-qc dashboard` refuses it with a message.

Every identifier starts with `SYNTH` and the experiment is flagged synthetic. The tests
`tests/annotation_qc/test_benchmark_*.py` use it for configuration-only extension, cache
reuse and invalidation, import matching, provenance, missing and failed results and the
HTTP API.

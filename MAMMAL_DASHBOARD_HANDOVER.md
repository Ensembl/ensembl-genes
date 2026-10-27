# Mammalian Annotation QC dashboard: handover

For a developer who has not worked on the project. It covers installing the release,
running pairwise comparisons, launching the dashboard from existing comparison outputs,
adding data, validating, and the limits of what the dashboard shows. Every command here
was run on 2026-10-27 (macOS 15, Python 3.12.14 and 3.13.13).

More detail lives in three places:

| topic | file |
|---|---|
| comparator: inputs, options, outputs, output contract, metric definitions | `src/python/ensembl/genes/annotation_qc/README.md` |
| dashboard and benchmark workflow: configuration keys, result status, provenance | `src/python/ensembl/genes/annotation_qc/dashboard/README.md` |
| metric dictionary (also the dashboard's *Metric dictionary* tab) | `src/python/ensembl/genes/annotation_qc/benchmark/metric_dictionary.json` |

---

## A. Project overview

Two independent stages.

```
reference annotation ─┐
prediction annotation ┼─► annotation-qc pairwise-compare ─► output bundle (contract 1.0)
genome FASTA ─────────┘     (all scientific calculation)      metrics.tsv, per-gene tables,
                                                              splice tables, manifest
                                                                   │
species configuration (JSON) ─► annotation-qc benchmark import/build ─► dashboard.sqlite
                                   (indexes the bundles; computes no metric)  │
                                                                   annotation-qc dashboard
                                                                   (local web page, read-only)
```

- **Comparator** (`metrics/pairwise`, `runners/pairwise_compare.py`, `reports/`, `parsers/`).
  It compares one prediction (the *query*) with one trusted reference and writes a
  self-describing output folder. It needs no database and no dashboard.
- **Dashboard** (`benchmark/`, `dashboard/`). It reads the comparator's files, mainly
  `metrics.tsv` (every displayed number), the per-gene tables, the one-to-one pair table,
  the splice tables and the evaluated-transcript export. It indexes them into one SQLite
  file per experiment and serves a local page. **It never recalculates a comparison
  metric.** Values are shown as stored, and a missing value is labelled
  (`not available`, `not applicable`, `not run`), never shown as 0.
- **Species-level experiment.** One JSON configuration per species: one reference (species,
  assembly, annotation release), any number of prediction runs, the transcript policies
  (A/B/C) and the analysis scopes. Species can have different prediction tools. The
  masthead selector switches species; every tab follows.
- **Display roles** (`runs[].display_role`) are presentation only. `focus` (Internal
  predictions) and `primary_comparator` (Tiberius) are listed first, get reserved colours
  and marker shapes, and are drawn on top. They never change a value (tested).

## B. Repository and release

| item | value |
|---|---|
| repository | `https://github.com/Ensembl/ensembl-genes` |
| release branch | `feature/annotation-qc-mammal-dashboard` |
| release tag (proposed) | `mammal-dashboard-v1.0.0` (annotated; created after review) |
| comparator baseline | `feature/annotation-qc-pairwise-compare` @ `febc0ca`: output contract `1.0`, definition `2026-10.1` |
| relationship | the release branch is `febc0ca` + merge of the contract-v1 dashboard (`f81328f`, merge `e55a21c`) + dashboard commits. The comparator code (`metrics/`, `runners/`, `reports/`, `parsers/`) is byte-identical to `febc0ca` |
| comparator build used for the mammalian results | `febc0ca`, registered in `benchmark/implementation_registry.json` (fingerprint `69f2821f…`), so these results show as *verified implementation febc0ca* |
| package | `ensembl-genes` 0.1.0, console script `annotation-qc` |

Not part of this release: fungi/apicomplexa work (`feature/annotation-qc-dashboard-fungi-apicomplexa`),
older dashboards (`feature/annotation-qc-dashboard*`).

## C. Installation

Requirements: Python ≥ 3.12 (tested 3.12.14 and 3.13.13), macOS or Linux, `git`, `pip`.
There are no compiled build steps: pyranges1 pulls in `ruranges`, which has binary wheels
for macOS (arm64/x86_64) and Linux. Windows is untested. No database access is needed.

```bash
git clone https://github.com/Ensembl/ensembl-genes.git
cd ensembl-genes
git checkout mammal-dashboard-v1.0.0          # or: git checkout feature/annotation-qc-mammal-dashboard
python3.12 -m venv ~/venvs/annotation-qc       # any path; python3.13 also works
source ~/venvs/annotation-qc/bin/activate
pip install --upgrade pip
pip install .                                   # current dependency versions
pip install "pytest>=7,<9"                      # for the tests
annotation-qc --help                            # lists pairwise-compare, benchmark, dashboard
```

To **reproduce the mammalian comparisons bit-for-bit**, install the exact package versions
that produced them, using the constraints file in the package:

```bash
pip install -c src/python/ensembl/genes/annotation_qc/constraints-mammal-v1.txt .
```

The full test suite passes with both the pinned versions (pandas 3.0.6, pyranges1 1.4.5)
and the current ones (pandas 3.1.0, pyranges1 1.5.0). pandas 4 is untested. Optional
extras: `pip install ".[plots]"` (matplotlib, for `--plots-per-category`) and
`pip install ".[dashboard]"` (PyYAML, for YAML configurations).

Verify the installation without any private data:

```bash
python -m pytest -q tests/annotation_qc        # expect: all passed, 6 skipped (opt-in real-data tests)
```

## D. Running pairwise comparisons

### D.1 Inputs

| input | required | notes |
|---|---|---|
| `--reference` | yes | trusted annotation, GFF3 or GTF; plain, gzip or BGZF (detected from content) |
| `--query` | yes | prediction: GFF3, GTF or Tiberius GTF hybrid; same compression rules |
| `--genome` | recommended (required for `--add-stop-codon`) | FASTA, plain/gzip/BGZF; checks sequence names and bounds; `--genome-mismatch error` (default) stops on a mismatch |
| `--seqname-map` / `--assembly-report` | when names differ | renames annotation sequences to the FASTA names |

Ensembl GFF3 with `unconfirmed_transcript` or `gene_segment` transcript types must be
retyped first. The human and mouse references needed this; see the comparator README,
"Tested example". `annotation-qc benchmark prepare` applies the recorded recipe.

### D.2 A single comparison (tested example, no private data)

The release ships a tiny synthetic species: a genome with real ORFs, a reference and two
predictors. The focus predictor's CDS leaves out the stop codon.

```bash
cp -r src/python/ensembl/genes/annotation_qc/examples/handover_mini /tmp/handover_mini
cd /tmp/handover_mini

annotation-qc pairwise-compare \
  --reference inputs/reference.gff3 \
  --query inputs/other.gff3 \
  --genome inputs/genome.fa --genome-mismatch error \
  --evaluation-mode protein_coding --reference-transcript-biotypes protein_coding \
  --reference-transcript-selection longest_cds --query-transcript-selection longest_cds \
  --analysis-scope chromosomes=1 \
  --export-evaluated-models ids \
  --outdir results/mini_other/original/B
```

- **Transcript policy** (one output folder per policy). The mammalian results use:

  | policy | `--reference-transcript-selection` | `--query-transcript-selection` | meaning |
  |---|---|---|---|
  | A | `all` | `all` | any eligible reference isoform reproduced (default view) |
  | B | `longest_cds` | `longest_cds` | one transcript per gene on both sides |
  | C | `canonical` | `longest_cds` | the `Ensembl_canonical` reference transcript (needs the tag; the run stops without it) |

- **Scopes.** The whole compared set is always reported as scope `all`.
  `--analysis-scope ID=SEQ,SEQ…` adds named subsets, such as the chromosome scope, to
  `metrics.tsv` and the splice tables (repeatable, or `--analysis-scopes-file`).
  `--region` / `--regions-file` instead *restrict* the comparison: the human/mouse runs use
  `--regions-file` with the sequences the reference annotates, which excludes alt
  haplotypes and patches.
- **Stop codons** (optional, query only). For a prediction whose CDS excludes the stop
  codon, add `--add-stop-codon --stop-codon-genetic-code 1`, plus
  `--stop-codon-sequence-code MT=2` for a vertebrate mitochondrion. The comparator then
  writes `stop_codon_harmonisation.tsv` and `query_evaluated.stop_harmonised.gtf`. In the
  example, the focus predictor goes from 2 to 3 one-to-one exact genes. By convention the
  output goes under `<run>/stop_harmonised/<policy>` instead of `<run>/original/<policy>`.

### D.3 The whole example (comparisons → dashboard dataset)

```bash
cd /tmp/handover_mini
bash run_example.sh        # 2 predictors × policies A and B, then validate, import, status, build
annotation-qc dashboard --config configs/mini_species.json --port 8765
```

`run_example.sh` is short and shows the conventions: output layout
`results/<run_id>/<variant>/<policy>/`, the stop-codon options, and the four
`benchmark` commands.

### D.4 Did a comparison succeed?

1. Exit status 0, with the last log lines naming the output folder. On any input error
   the run exits non-zero and says why: unknown parent, sequence missing from the FASTA,
   or a canonical selection without tags.
2. `comparison_manifest.json` exists. It has `contract_version` `1.0`, `definition_version`
   `2026-10.1`, the SHA-256 of every input and output, all options, and the
   implementation fingerprint.
3. `annotation-qc benchmark import` then `annotation-qc benchmark status --config …` shows
   every (run, policy) as `complete`: inputs, options and definition match the
   configuration, and every output file matches its manifest checksum. Other statuses:
   - `stale`: something differs (the reason is named).
   - `unverified`: the code cannot be identified; register the build (dashboard README,
     section 5).
   - `not_run`: no result found.
   - `unavailable`: e.g. policy C without canonical tags.

### D.5 Output bundle

| file | used by the dashboard for |
|---|---|
| `metrics.tsv` | every headline, structural and all-metrics value, per scope |
| `comparison_details.tsv` | gene table and locus view (per-gene CDS outcome, best pair, CDS differences) |
| `cds_exact_one_to_one_pairs.tsv` | one-to-one pairs (download, pair check) |
| `consensus_transcript_labels.tsv`, `gene_splits.tsv`, `gene_merges.tsv` | Predictions tab, split/merge filters |
| `evaluated_transcripts.tsv.gz` | locus view (which isoforms were evaluated) |
| `splice_sites_by_position.tsv`, `cds_introns_by_length.tsv` | splice-site position and intron-length plots (policies B/C only) |
| `per_sequence.tsv` | per-sequence download |
| `stop_codon_harmonisation.tsv` | stop-codon audit, red ticks in the locus view |
| `comparison_summary.json` | limitations and diagnostics shown beside results |
| `comparison_manifest.json` | matching to configured runs (by checksum), status, provenance |

Column-level definitions are in the comparator README, "Outputs" and "Output contract".

## E. Running the dashboard

The dashboard needs only the species configurations and their built datasets
(`<output_dir>/dashboard/dashboard.sqlite`). All paths in a configuration are relative to
the configuration file, so a folder holding both can be moved or copied as a whole.
Layout of the mammalian data:

```
<DATA>/
  qc_configs/mammal_dashboard/{cat,dog,rat,pig,human,mouse}.json
  qc_workspace/mammal_<species>/dashboard/dashboard.sqlite
  contract_runs/20261021_*/04_runs/...      (comparator bundles; needed to rebuild or verify, not to view)
```

One species:

```bash
cd <DATA>
annotation-qc dashboard --config qc_configs/mammal_dashboard/pig.json --port 8765
```

All six species (the selector lists them in this order):

```bash
cd <DATA>
annotation-qc dashboard \
  --config qc_configs/mammal_dashboard/cat.json \
  --config qc_configs/mammal_dashboard/dog.json \
  --config qc_configs/mammal_dashboard/rat.json \
  --config qc_configs/mammal_dashboard/pig.json \
  --config qc_configs/mammal_dashboard/human.json \
  --config qc_configs/mammal_dashboard/mouse.json \
  --port 8765
```

Open `http://127.0.0.1:8765/`.

- **Another data folder:** pass the configurations by path, from any working directory,
  e.g. `--config /data/mammal_qc/qc_configs/mammal_dashboard/cat.json`. Datasets are found
  relative to each configuration. Tested by launching the six species from a different
  folder.
- **Port:** `--port N`. `lsof -i :8765` shows what holds a port.
- **Bind address:** `--host` (default `127.0.0.1`). The dashboard has no authentication;
  keep it local.
- **Restart** after every `benchmark build`.

Troubleshooting:

| message | cause and fix |
|---|---|
| `Dashboard dataset not built for {...}` | the configuration's `output_dir/dashboard/dashboard.sqlite` is missing: wrong folder layout, or run `annotation-qc benchmark build --config …` |
| `… has N references; the dashboard shows one species experiment per configuration` | split the configuration: one reference per file |
| `experiment id … is configured twice` | two configurations share `experiment.id` |
| a configuration error listing problems | JSON keys or values invalid; every problem is listed at once |
| `result directory not found: …` (from `benchmark import`) | an `import.results_dirs` entry does not exist at that path |
| a run shows `not run` | its result folders are not under `import.results_dirs`, or the checksums or options do not match the run: `annotation-qc benchmark status --config …` names the reason |
| a run shows `stale` | the result was made by another definition version, or with other inputs or options. Rerun it (section F) |
| `WARNING: … still has data in dashboard.sqlite-wal` (after a build) | a running dashboard was reading the dataset; stop it and run `benchmark build` again before copying the dataset |
| `Cannot listen on 127.0.0.1:8765 (Address already in use)` | another dashboard holds the port: choose another `--port` |

## F. Adding data

### F.1 A new prediction dataset for an existing species

1. Keep the prediction file unchanged, e.g. under `<DATA>/gene_prediction_outputs/`.
2. Run the comparator for policies A, B and C with the **same options as that species'
   other runs**. The exact command lines are recorded in
   `<DATA>/contract_runs/20261021_*/logs/command_ledger.tsv`: copy one, then change
   `--query` and `--outdir`. Human, for example (from `<DATA>`):

   ```bash
   for P in A B C; do
     case $P in A) R=all; Q=all;; B) R=longest_cds; Q=longest_cds;; C) R=canonical; Q=longest_cds;; esac
     annotation-qc pairwise-compare \
       --query gene_prediction_outputs/<new file> \
       --reference benchmark_runs/20261005_145546/07_derived_inputs/references/human_GRCh38.parser_compatible.gff3 \
       --genome benchmark_runs/20261005_145546/01_inputs/decompressed/human_GRCh38.fa --genome-mismatch error \
       --evaluation-mode protein_coding --reference-transcript-biotypes protein_coding \
       --reference-transcript-selection $R --query-transcript-selection $Q --export-evaluated-models ids \
       --regions-file benchmark_runs/20261005_145546/01_inputs/scope/human.reference_sequences.regions \
       --analysis-scope chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT \
       --outdir contract_runs/<new workspace>/04_runs/human/human_newtool/original/$P
   done
   ```

   For cat/dog/rat/pig use `references/<species>.gff3.bgz` and `references/<species>.fa.bgz`
   directly. Pig adds
   `--analysis-scope chr1_18_X=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X`.
3. **Register it** in the species configuration. Add a run with a new `run_id`,
   `reference`, `tool` (the name shown) and `annotation`. Add `comparator_options` if you
   used the stop-codon options. Add its folder (`…/<run>/<variant>`) to
   `import.results_dirs`.
4. Index and verify:

   ```bash
   annotation-qc benchmark validate --prepare --config qc_configs/mammal_dashboard/human.json
   annotation-qc benchmark import --config qc_configs/mammal_dashboard/human.json
   annotation-qc benchmark status --config qc_configs/mammal_dashboard/human.json
   annotation-qc benchmark build --config qc_configs/mammal_dashboard/human.json
   ```

   Then restart the dashboard. Results are matched to runs by content (query and reference
   SHA-256, stop-codon options, transcript selections), never by folder name.

**Display priority:** set `"display_role": "focus"` on the dataset of interest and
`"display_role": "primary_comparator"` on its main comparison, then rebuild. All other
runs are ordered by `display.order_metric` (default one-to-one exact-CDS recall).

**Human Internal predictions (open item):** add a run with
`"tool": "Internal predictions"` and `"display_role": "focus"` for the corrected file,
compare it as above (no stop-codon options when the file already includes stop codons),
then rebuild. No code change is needed.

### F.2 A new species experiment

Copy a configuration (e.g. `cat.json`). Give it a new `experiment.id` and `output_dir`,
then set `references[0]`: `id`, `species`, `assembly`, `assembly_accession`,
`annotation_release`, `annotation`, `genome`, `scope`, and optional `analysis_scopes`
(each also passed to the comparator as `--analysis-scope`). Then list the runs and
`import.results_dirs`. If the reference needs preparation (retyping, decompression),
declare it in `annotation_preparation`, as in `human.json`. Then run
validate / import / build and add `--config <new>.json` to the launch command. Template
with every key: `src/python/ensembl/genes/annotation_qc/benchmark/templates/experiment_template.json`.

When the comparator's `definition_version` changes, rerun **every** result of a species.
The dashboard shows results of another definition as `stale` and never mixes them.

## G. Validation

| what | command | needs |
|---|---|---|
| all automatic tests (comparator, dashboard, end-to-end smoke test) | `python -m pytest -q tests/annotation_qc` | nothing private |
| comparator only | `python -m pytest -q tests/annotation_qc/test_pairwise_*.py` | nothing private |
| dashboard and benchmark | `python -m pytest -q tests/annotation_qc/test_benchmark_*.py tests/annotation_qc/test_mammal_dashboard.py` | nothing private |
| end-to-end smoke test: annotations → comparator CLI → bundles → config → import/build → dashboard API; moved folder; display roles change no value; missing results | `python -m pytest -q tests/annotation_qc/test_release_smoke.py` | nothing private (uses `examples/handover_mini`) |
| real mammalian data: each configuration parses; one reference; runs and display roles as configured; every result `complete` with one definition and build; Internal predictions stop-corrected; every served value equals `metrics.tsv`; scopes offered = scopes reported; species differ in predictors | `ANNOTATION_QC_MAMMAL_CONFIGS=<DATA>/qc_configs/mammal_dashboard python -m pytest -q tests/annotation_qc/test_mammal_release_data.py` | the data package (skipped otherwise; the metric check is skipped when the comparator bundles are absent) |

Results on 2026-10-27:
- Source checkout: 227 passed, 6 skipped. The skips are the opt-in real-data tests.
- Clean installs from the wheel: 227 passed, 6 skipped, on Python 3.13 with current
  dependencies and on Python 3.12 with the constraints file.
- Opt-in real-data run: 31 passed.

## H. Scientific interpretation (summary)

Definitions are in the *Metric dictionary* tab and in the comparator README
("Output contract", "Interpreting the metrics"). The dashboard shows them as stored.

- **R and Q** are the evaluated reference and predicted genes in the selected policy and
  scope. Recall-type metrics divide by R (reference side). Precision-type metrics divide
  by Q (prediction side). F1 = 2TP / (R + Q).
- **Transcript policies** decide which transcripts are compared: A all, B longest CDS on
  both sides, C canonical reference and longest predicted. Each policy is a separate
  comparator run, so its values come from its own `metrics.tsv`.
- **Exact CDS:**
  - **One-to-one exact CDS** (headline recall, precision, F1) is a maximum one-to-one
    matching of genes with identical CDS (same start, every splice site, same stop).
  - **Coordinate-exact** counts (per gene, best CDS pair) are a different measure and can
    exceed TP.
- **Loci recovered** means a same-strand CDS overlap: the locus was found, not that the
  model is right.
- **CDS outcome bars** show the comparator's seven mutually exclusive per-gene outcomes
  (`cds_outcome_*`), which sum to R.
- **Intron chains:** a multi-segment gene whose CDS intron chain is identical in the best
  pair, or in any pair. Unique-intron recovery and support are sets of
  (sequence, strand, start, end).
- **Splice-site plots** are computed by the comparator for one reference transcript per
  gene, so they appear under policies B and C only. **Genomic CDS** = genomic distance
  from the coding start to the junction over the CDS span (introns included). **Coding
  sequence** = coding bases upstream of the junction over the CDS length. Both follow the
  direction of transcription. *Introns ≥ 60 bp* leaves out short gaps. The intron-length
  plot bins reference CDS introns by length.
- **Scopes** are the comparator's own rows for that subset (pairing and matching never
  cross sequences), never sums of per-sequence counts.

## I. Known limitations

1. **Human and mouse have no Internal predictions.** No such file existed when the
   results were produced. Human/mouse show Tiberius as primary comparison plus ANNEVO,
   Helixer and Vipsania.
2. **Policy A has no splice-site or intron-length plots.** The comparator computes them
   only with one reference transcript per gene; the page says so and offers B/C.
3. **Policy C needs `Ensembl_canonical` tags** in the reference. Without them, policy C
   is `unavailable`, as in the synthetic example.
4. **Display needs the built dataset.** Rebuilding a dataset (`benchmark build`) also
   needs the original annotation inputs, because locus models are rebuilt from them and
   checked against the evaluated-transcript export. The comparator bundles alone are not
   enough to rebuild.
5. **Provenance shows original paths.** The Provenance tab shows the absolute paths
   recorded when the data were built (e.g. `/Users/…`). They are records, not used to load
   anything.
6. **All results of a species must share one comparator definition.** Mixing
   `definition_version`s shows `stale` results.
7. **Comparator memory and time.** Human whole-genome runs need about 3.5–7 GB RAM and
   5–8 min per policy (about 40 min under memory pressure). Run human comparisons one at a
   time on a 16 GB machine.
8. **Prediction processing scope is not recorded.** For most predictors, reference genes
   on sequences without predictions count as missed. The input audits show these numbers
   where available.
9. **pandas 4 is untested.** pandas 3.1 already warns about two calls that pandas 4 will
   reject.
10. **Local only.** The server binds to 127.0.0.1 with no authentication. Not for shared
    hosting.
11. **Not covered by this release:** fungi/apicomplexa experiments, configurations with
    several references, YAML configurations without the `dashboard` extra, Windows.
12. **macOS-specific helper scripts.** The data-folder runner scripts
    (`contract_runs/*/scripts/run_comparisons.py`) call `/usr/bin/time -l` (macOS). On
    Linux use `-v`, or run `annotation-qc pairwise-compare` directly.

---

## Appendix 1. Data package for the mammalian dashboard

Sizes measured on 2026-10-27 in the original data folder. Nothing has been copied yet.

**A. Runtime (essential to view the dashboard exactly as now):** about 4.7 GB,
about 1 GB compressed (gzip -1 compresses `dashboard.sqlite` 4.9×).

| path | size |
|---|---|
| `qc_configs/mammal_dashboard/*.json` (6) | 92 KB |
| `qc_workspace/mammal_cat/dashboard/dashboard.sqlite` | 412 MB |
| `qc_workspace/mammal_dog/dashboard/dashboard.sqlite` | 439 MB |
| `qc_workspace/mammal_rat/dashboard/dashboard.sqlite` | 380 MB |
| `qc_workspace/mammal_pig/dashboard/dashboard.sqlite` | 499 MB |
| `qc_workspace/mammal_human/dashboard/dashboard.sqlite` (+ 866 MB `-wal`, see below) | 1,285 MB |
| `qc_workspace/mammal_mouse/dashboard/dashboard.sqlite` | 1,117 MB |

**Before copying,** stop every dashboard and run
`annotation-qc benchmark build --config qc_configs/mammal_dashboard/human.json`, so the
human dataset's write-ahead log is folded into `dashboard.sqlite`. The build skips
unchanged work, and the command prints a warning if a reader prevents it. Then copy only
`dashboard.sqlite`. The other five already have empty logs.

**A+. Comparator output bundles (recommended):** about 2.3 GB. They are the source of
truth, used by the real-data test's metric check and for re-indexing.

| path | size |
|---|---|
| `contract_runs/20261021_cat_dog_rat_pig_def2026-10.1/04_runs` (51 results) | 1.7 GB |
| `contract_runs/20261021_human_mouse_contract_v1/04_runs` (24 results) | 606 MB |
| `contract_runs/20261021_*/{00_env,scripts,logs,05_validation}` (wheel, pinned versions, command ledgers, validation) | ~1 MB |

Largest files: `comparison_details.tsv` (918 MB in total), `query_evaluated.stop_harmonised.gtf`
(801 MB), `evaluated_transcripts.tsv.gz` (305 MB). The stop-harmonised GTFs are only for
browsers and plots; the dashboard does not read them.

**B. Reproducibility inputs (to rerun the comparator):** about 5.3 GB minimum.

| path | size |
|---|---|
| `references/{cat,dog,rat,pig}.{gff3,fa}.bgz`, `references/{human_GRCh38,mouse_GRCm39}.{gff,fa}.bgz` | 4.1 GB |
| `gene_prediction_outputs/` (the 28 files named in the configurations) | 1.2 GB |
| `benchmark_runs/20261005_145546/01_inputs/scope/{human,mouse}.reference_sequences.regions` | 8 KB |
| `contract_runs/20261012_151131_comparator_contract_v1/02_audit/*.json` (input audits shown in Provenance) | 136 KB |

Derivable from B and not needed in the package:
- The human/mouse parser-compatible references (2.4 GB): `benchmark prepare`, recipe
  recorded, SHA-256 checked.
- The human/mouse decompressed FASTA (5.7 GB): `gzip -dc`.

**C. Not needed:**
- Superseded result sets: `contract_runs/20261012_*` (except `02_audit`),
  `contract_runs/20261014_*`, `pig_comparison_runs/`.
- `benchmark_runs/` (except the regions files and derived references above).
- Older dashboards and workspaces: `qc_workspace/` except `mammal_*/dashboard`,
  `dashboard_runs/`, `dashboard_fungi_apicomplexa/`.
- Investigations and demos: `investigations/`, `qc_planning/`,
  `qc_synthetic_demo_NOT_REAL_DATA/`.
- Small workspace bookkeeping: `qc_workspace/mammal_*/{cache,prepared,imported,validation}`
  (regenerated by import) and `verify_metrics.tsv`.

## Appendix 2. Release commands (run after review)

```bash
cd ensembl-genes                                  # the checkout on feature/annotation-qc-mammal-dashboard
git status --short --branch && git log --oneline -1          # 1. clean tree, expected HEAD
git add <files> && git commit -m "<message>"                 # 2. only if further changes were made
git tag -l mammal-dashboard-v1.0.0; git ls-remote --tags origin mammal-dashboard-v1.0.0   # must print nothing
git tag -a mammal-dashboard-v1.0.0 -m "Mammalian Annotation QC dashboard 1.0.0 (comparator febc0ca, contract 1.0, definition 2026-10.1)"   # 3.
git push -u origin feature/annotation-qc-mammal-dashboard    # 4. publish the branch
git push origin mammal-dashboard-v1.0.0                      # 5. publish the tag
# 6. on another computer
git clone https://github.com/Ensembl/ensembl-genes.git && cd ensembl-genes && git checkout mammal-dashboard-v1.0.0
```

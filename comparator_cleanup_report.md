# Comparator cleanup report: `annotation-qc pairwise-compare`

These fixes address the comparator issues found during the P. falciparum Vipsania vs Tiberius
evaluation (`ensembl_anno_apicomplexa/plasmodium_dl_compare_oct2/01_comparator_audit`).

- **Not changed:** locus pairing, the 0.8 reciprocal-overlap threshold, the exon-level class
  definitions and the CLI options.
- **Nothing committed or pushed.**

Environment: conda env `ensembl-genes-aqc` (Python 3.12, pandas 3.0.6, pyranges1 1.4.5), with this
checkout on the import path.

## Bugs fixed

### 1. `matched_consensus_count` double counting

**Old code (`summary.py`):**

```
matched_consensus_count = count(Matched) + len(query) − count(Novel) − count(Strand_Mismatch)
```

`len(query) − Novel − Strand_Mismatch` already includes the `Matched` query genes, so those genes
were counted twice. A `Matched` query gene is one with a same-strand reference partner but no
comparable transcripts. In a one-gene example the old code reported 2 matched query genes.

**New definition:** the number of query genes with a same-strand reference partner (classes
Exact_Match, Partial_Match, Structural_Mismatch or Matched), each counted once. Two invariants now
hold:

- `matched_consensus_count ≤ total_consensus_genes`
- matched + novel + `strand_mismatch_consensus_count` = total query genes

`matched_consensus_rate` is also reported, with the query genes as denominator.

On P. falciparum neither run had a `Matched` query gene, so the value is unchanged: Vipsania 5,262
and Tiberius 3,703. The new invariant sums hold: 5,262 + 156 + 29 = 5,447 and
3,703 + 95 + 29 = 3,827.

### 2. CDS class inherited from the best exon pair (FIX 4)

**Old behaviour:** `_best_match` ranked transcript pairs by exon similarity, which includes UTRs. It
then took the CDS class from that same pair. A query CDS identical to reference isoform B was
therefore scored against isoform A whenever A shared more UTR.

**New behaviour:** the best CDS pair is chosen in its own pass over all (own transcript, partner
transcript) pairs that both have CDS. Pairs are ranked by:

1. identical CDS
2. priority
3. reciprocal CDS overlap

The CDS pair is reported separately from the exon pair, using four columns:

- `cds_matched_id`
- `best_cds_match_transcript_id`
- `best_match_own_transcript_id`
- `best_cds_match_own_transcript_id`

Related change: `No_CDS` now means none of the gene's transcripts has CDS. Previously it meant the
best *exon* pair's own transcript had no CDS. The old rule could report a coding gene as `No_CDS`
when its best exon pair was a non-coding isoform.

### 3. Strand_Mismatch triggered by gene-span contact (FIX 5)

**Old behaviour:** a gene was `Strand_Mismatch` when its gene span overlapped only opposite-strand
genes, even if the only contact was through an intron or UTR.

**New behaviour:** a gene is `Strand_Mismatch` only when it shares CDS bases with an
opposite-strand gene. When either gene lacks CDS, shared exon bases are used instead. The new
`strand_mismatch_basis` column records `CDS` or `exon`.

Genes whose opposite-strand contact is span-only, intron-only or UTR-only are now `Missed`
(reference) or `Novel` (query), with `strand_mismatch_basis = no_feature_overlap`. They are counted
in `opposite_strand_no_feature_overlap_count`. The old gene-span number can still be recovered as
`strand_mismatch_count + opposite_strand_no_feature_overlap_count`.

## Semantic/reporting fixes

### Coordinate-exact metrics (FIX 2)

**What changed:**

- The `Exact_Match` class is kept so downstream consumers keep working, but it is now documented as
  *structural*: identical intron chain and ≥ 0.8 reciprocal overlap. Terminal coordinates may differ.
- New per-gene columns:
  - `cds_coordinate_exact`: the best CDS pair has identical CDS intervals, start and stop included.
  - `exon_coordinate_exact`: the same test on exons.
- New summary keys:
  - `cds_coordinate_exact_count` / `_rate`
  - `exon_coordinate_exact_count` / `_rate`
  - `cds_exact_match_not_coordinate_exact`
  - query-side `consensus_cds_coordinate_exact_*` and `consensus_exon_coordinate_exact_*`

### Single-exon intron chains (FIX 3)

**Old behaviour:** empty chain == empty chain counted as a match. Every overlapping single-exon
pair, even a partial overlap, was counted as "intron chain recovered".

**New behaviour:**

- `intron_chain_match` and `cds_intron_chain_match` are now nullable. They are `NA` when the gene's
  compared structure has no intron, and `comparison_details.tsv` writes them as `NA`.
- A multi-exon reference matched by a single-exon query is still `False`, not NA.
- Pairing and ranking are unchanged: single-exon pairs can still be Exact_Match.

**New denominators:**

| key | denominator |
|---|---|
| `multi_exon_reference_genes`, `multi_segment_cds_reference_genes` | (these are the denominators) |
| `*_intron_chain_sensitivity` | recovered / all multi-exon (multi-CDS-segment) reference genes |
| `*_intron_chain_rate` | recovered / detected reference genes whose compared structure has an intron |
| `*_intron_chain_matched` | same as `_evaluated`; old name kept |

`*_intron_chain_not_applicable_single_exon` reports how many detected genes were excluded as single
exon.

### Split and merge (FIX 6)

A reference gene and a query gene are **counterparts** when both of these hold:

- they are same-strand locus partners (the existing pairing)
- some transcript pair shares at least 10% of the shorter structure: CDS when both genes have CDS,
  exons otherwise

Gene-span, intron, UTR-only and opposite-strand contact do not qualify.

- **Split:** a reference gene with ≥ 2 query counterparts.
- **Merge:** a query gene with ≥ 2 reference counterparts.

Outputs:

- columns `counterpart_count` and `counterpart_ids`
- summary block `split_merge`: counts, genes involved, one-to-one count and the threshold
- new files `gene_splits.tsv` and `gene_merges.tsv`, with gene ID, seqid, one-based coordinates,
  strand, class, number of counterparts, counterpart IDs and counterpart loci

The 10% threshold is a new parameter (`SPLIT_MERGE_MIN_SHARED_FRACTION`). It is the definition the
independent analysis already used, and it was not tuned on these results.

### Novel context (FIX 7)

**Old behaviour:** in `protein_coding` mode, a prediction on a reference pseudogene was reported only
as `Novel`.

**New behaviour:**

- Novel query genes keep the `Novel` class and are never counted as matches.
- Each Novel gene gets a `novel_category` from exon overlap, on either strand, with the reference
  **before** mode/biotype filtering:
  - `Overlaps_reference_pseudogene`
  - `Overlaps_other_non_protein_coding_reference`
  - `Overlaps_unevaluated_protein_coding_reference`
  - `Overlaps_evaluated_reference_opposite_strand`: UTR/exon contact only, new because of FIX 5
  - `Novel_no_reference_overlap`
- The runner passes the unfiltered reference as `reference_context`.
- The categories appear in `comparison_details.tsv`, in `consensus_transcript_labels.tsv`, in the
  `novel_categories` summary block and in the console summary.

### Denominators

- `comparison_summary.json` has a new `rate_denominators` block that names the denominator of every
  rate.
- **Reference-based (sensitivity-like) rates** divide by `total_reference_genes`, or by the
  multi-exon counts for intron-chain sensitivity. They live in `sensitivity`, `sensitivity_cds`,
  `intron_chain` and `locus_detection_exonic`.
- **Query-based (precision-like) rates** divide by `total_consensus_genes` and live in
  `specificity`. Every rate there is new, and each has a matching count.
- **Intron support (new):** the `intron_support` block compares unique introns. It reports
  reference introns recovered / reference introns and query introns supported / query introns,
  separately for exon-level and CDS-level introns.
- **Console summary:** split into a reference-based section and a query-based section, and each
  intron-chain line now names its denominator.
- **TSV summary:** now also writes `intron_chain`, `novel_categories`, `split_merge` and
  `intron_support`. The old TSV left out `intron_chain`.

## Tests added

The new file `tests/annotation_qc/test_pairwise_regressions.py` has 24 tests. Two existing tests in
`test_pairwise_metrics.py` were updated to the corrected semantics:

- single-exon intron chain is now NA
- the intron-chain and specificity blocks have new keys

| # | requirement | test |
|---|---|---|
| 1 | matched_consensus_count arithmetic | `TestMatchedConsensusCount` (Matched class counted once; all five query classes, sum invariant, ≤ total) |
| 2 | true coordinate-exact CDS | `TestCoordinateExact::test_identical_coordinates` |
| 3 | terminal CDS difference not coordinate-exact | `test_same_intron_chain_different_start`, `test_different_terminal_exon_end`, `test_partial_overlap` |
| 4 | single-exon intron chain = N/A | `TestSingleExonIntronChain::test_single_exon_cases_are_not_applicable` (exact, partial and no overlap) |
| 5 | multi-exon exact intron chain | `test_multi_exon_exact_and_single_exon_query_on_multi_exon_reference` |
| 6 | independent CDS vs exon best pair | `test_cds_pair_independent_of_exon_pair` (isoform A best exon pair, isoform B exact CDS → B reported), `test_no_cds_uses_any_coding_isoform` |
| 7 | real CDS strand mismatch | `TestStrandMismatch::test_opposite_strand_cds_overlap`, `test_exon_overlap_when_one_side_has_no_cds` |
| 8 | UTR-only opposite strand not misclassified | `test_opposite_strand_utr_only_overlap`, `test_neighbouring_opposite_strand_gene_span_overlap`, `test_same_strand_cds_overlap` |
| 9 | gene split | `TestSplitMerge::test_split`, `test_one_to_one` |
| 10 | gene merge | `TestSplitMerge::test_merge` |
| – | non-qualifying neighbours | `test_non_qualifying_neighbours` (UTR-only span contact; 5 bp CDS contact), `test_opposite_strand_is_not_a_counterpart` |
| 11–13 | true novel, pseudogene overlap, non-PC overlap, not counted as matches | `TestNovelContext` |
| 14 | denominator sanity | `test_denominators` (reference vs query denominators, multi-exon intron-chain denominators, intron support, all rates in [0, 1]) |
| – | CLI end-to-end | `test_runner_writes_refined_outputs` (split/merge TSVs, new detail/label columns, NA formatting, summary TSV rows) |

### Results

| run | result |
|---|---|
| Baseline, before any change (`tests/annotation_qc`) | 55 passed, 1 skipped (entry-point test skipped: `annotation-qc` not on PATH) |
| After (`tests/annotation_qc`) | 79 passed, 1 skipped |
| After, with the env's `bin/` on PATH (`tests/annotation_qc`) | **80 passed, 0 skipped**, including the installed CLI entry point |
| After, whole `tests/` directory, env on PATH | **81 passed** |
| `--plots-per-category 2` on the real Tiberius data | ran cleanly; plots and `qc/index.html` written |
| `black --check` | clean. Only files changed here were reformatted; the originals were already clean. |

The bugs were also reproduced against a snapshot of the original code (scratch script, outside the
repo). The same synthetic fixtures gave:

- `matched_consensus_count` 2 for 1 query gene
- single-exon partial overlap reported as intron chain recovered (1/1)
- CDS identical to isoform B classified `Partial_Match` (0.65)
- 3′-UTR-only antisense contact reported as `Strand_Mismatch` on both sides
- no split reported
- pseudogene-overlapping prediction reported as plain `Novel`

The fixed code gives the expected result in every case.

## Before/after behaviour (P. falciparum, protein_coding, all isoforms)

**Inputs:** identical to `scripts/p5_run_comparator.sh`, without plots.

- "Before": the existing `02_vipsania_compare` / `03_tiberius_compare` outputs. Their source
  checksums match the pre-change code exactly (`comparator_source_checksums.tsv`).
- "Independent": `04_independent_metrics` and `comparator_vs_independent.md`.
- The rerun outputs were written to the session scratchpad and are not in the repo. Rerun the
  command from `p5_run_comparator.sh` to regenerate them.

| metric | Vipsania before → after | independent | Tiberius before → after | independent |
|---|---|---|---|---|
| matched_consensus_count (of Q) | 5,262 → 5,262 (of 5,447) | – | 3,703 → 3,703 (of 3,827) | – |
| CDS coordinate-exact (of R = 5,285) | – → **4,741** | 4,741 | – → **1,521** | 1,521 |
| CDS Exact_Match (structural) | 5,003 → 5,032 | – | 1,741 → 1,743 | – |
| …structural but not coordinate-exact | – → 291 | 289 (vs old 5,003) | – → 222 | 222 |
| coordinate-exact CDS on a different ref isoform than the best exon pair | – → 28 | – | – → 3 | – |
| …of which previously not CDS Exact_Match | 27 | 27 | 2 | 2 |
| CDS intron chain recovered | 5,003 / 5,259 detected → **2,646 / 2,859** multi-segment (92.5%) | 2,646 / 2,859 | 1,914 / 3,413 → **213 / 2,859** (7.5%) | 214 / 2,859 |
| exon intron chain recovered | 4,827 / 5,259 (91.8%) → 2,527 / 2,948 multi-exon (85.7%) | – | 1,853 / 3,413 (54.3%) → 204 / 2,948 (6.9%) | – |
| Strand_Mismatch (reference) | 1 → **0** | 0 | 45 → **8** (all CDS) | 8 |
| Missed (reference) | 25 → 26 | – | 1,827 → 1,864 (+37 former span/UTR "strand mismatches") | – |
| Strand_Mismatch (query) | 34 → 29 | – | 32 → 29 | – |
| gene splits | proxy 18 → **7** | 7 | proxy 265 → **263** | 263 |
| gene merges | proxy 15 → **16** | 16 | proxy 1 → **0** | 0 |
| Novel | 151 → 156 | – | 92 → 95 | – |
| …overlaps pseudogene | (plain Novel) → **135** | 135 | (plain Novel) → **84** | 83 + 1 opposite strand |
| …overlaps other non-PC | → 1 | 1 ncRNA | → 0 | – |
| …evaluated PC, opposite strand, UTR/exon only | → 4 | – | → 3 | – |
| …no reference overlap | → **16** | 16 | → **8** | 8 |
| CDS introns: ref recovered / query supported | – → 96.6% / 95.8% | 96.6% recovered | – → 4.8% / 54.4% | 4.8% recovered |

**Per-gene agreement with the independent tables:**

| comparison | agreement |
|---|---|
| `cds_coordinate_exact` vs `exact_cds` | identical for all 5,285 reference genes, both tools |
| split flag vs `is_split` | identical, both tools |
| `cds_intron_chain_match` vs `cds_intron_chain` | identical for Vipsania; one Tiberius gene differs (see limitation 1) |

**Unchanged by design:** locus detection (5,259 / 3,413), exon-level classes apart from the
strand-mismatch reclassification, `best_match_transcript_id` and `exon_overlap`. Zero genes
changed in those columns.

## Remaining known limitations

1. **Intron chains are still evaluated on the single best pair.** Ranking puts a ≥ 0.8-overlap pair
   ahead of a lower-overlap pair whose intron chain is identical. Examples:
   - Tiberius PF3D7_1461800 is split into g1991 + g1992. The better-overlapping piece is reported,
     and the piece with the identical CDS chain is not.
   - At exon level, 26 (Vipsania) and 3 (Tiberius) reference genes have an identical intron chain in
     some non-best pair. The independent "any pair" metric counts them.

   Changing this would change the ranking, which was out of scope.
2. **Locus pairing is still gene-span based.** Partial_Match can rest on span-only contact; this is
   already reported as `span_only_detected_count`, which is 0 here.
3. **Class names were kept for compatibility.**
   - `Exact_Match` remains the structural class. Use `cds_coordinate_exact` / `exon_coordinate_exact`
     for true coordinate identity.
   - Novel query genes keep `classification = Novel`. The refinement is in `novel_category`.
   - `exon_intron_chain_matched` / `cds_intron_chain_matched` keep their old names but now count
     only multi-exon (multi-CDS-segment) genes.
4. **Novel context depends on `reference_context`.** It uses exon overlap on either strand, so a
   prediction antisense to a pseudogene is labelled `Overlaps_reference_pseudogene`. When
   `classify_loci` is called from Python without `reference_context`, the categories are computed
   against the evaluated reference only; the CLI always passes the unfiltered reference.
5. **The split/merge threshold is fixed in code.** The 10% shared fraction is not a CLI option.
6. **Same-strand partners still take precedence.** A gene with any same-strand partner is never
   Strand_Mismatch, even if it also has opposite-strand CDS overlap. This is unchanged behaviour.
7. **Isoform choice affects the intron-chain denominator.** A multi-exon reference gene whose best
   exon pair uses a single-exon isoform gets `NA`, but it is still counted in the multi-exon
   denominator as not recovered.
8. **Intron support covers the whole query.** `query_introns_supported_rate` counts query introns at
   loci outside the evaluated reference set (e.g. pseudogene loci in `protein_coding` mode) as
   unsupported.

## Files changed

All changes are uncommitted, on branch `feature/annotation-qc-pairwise-compare`.

| file | change |
|---|---|
| `src/python/ensembl/genes/annotation_qc/metrics/pairwise/classify.py` | independent CDS pair; coordinate-exact flags; nullable intron-chain status; feature-aware strand mismatch; counterparts (split/merge); Novel categories; `max_cds_count`/biotype in gene models; docstring definitions |
| `src/python/ensembl/genes/annotation_qc/metrics/pairwise/summary.py` | `matched_consensus_count` fix; coordinate-exact, intron-chain denominators, query-based rates, `novel_categories`, `split_merge`, `rate_denominators`; `summarise_intron_support`; extra label columns |
| `src/python/ensembl/genes/annotation_qc/reports/pairwise_compare.py` | new detail columns appended after `original_gene_id`; `NA` formatting; `gene_splits.tsv` / `gene_merges.tsv`; extra summary TSV rows; console summary with explicit denominators |
| `src/python/ensembl/genes/annotation_qc/runners/pairwise_compare.py` | passes the unfiltered (region-subset) reference as Novel context; adds `intron_support`; writes split/merge tables; help text |
| `src/python/ensembl/genes/annotation_qc/README.md` | metric definitions, outputs, denominators, change log entry 7 |
| `tests/annotation_qc/test_pairwise_metrics.py` | two tests updated to the corrected semantics |
| `tests/annotation_qc/test_pairwise_regressions.py` | new: 24 regression tests |
| `comparator_cleanup_report.md` | this report |

The parsers, selection policies, plotting, CLI options, source annotations and GMB were not changed.

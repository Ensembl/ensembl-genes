# Annotation QC

## Compare an annotation to a reference

`annotation-qc pairwise-compare` compares any **query** annotation (gene predictions from
Tiberius, ANNEVO, Helixer, Vipsania, AUGUSTUS, a Gene Model Builder build, …) with a
**trusted reference** annotation of the *same assembly*. It reports how many reference
genes the query reproduces (reference-based, recall-like) and how many query genes have
no reference counterpart (query-based, precision-like). The reference is used for
evaluation only. The command needs no GMB installation, build directory or metadata.

### Environment

Install `ensembl-genes` into a Python ≥ 3.12 environment. The command uses only the
package's declared dependencies (pandas, pyranges1 ≥ 1.3, pyfaidx). Plots need the
optional `plots` extra (matplotlib).

```bash
pip install -e "ensembl-genes"            # core comparison
pip install -e "ensembl-genes[plots]"     # optional: --plots-per-category
```

Resources: a whole human Ensembl reference (≈ 1.4 GB GFF3, ≈ 650,000 transcripts) takes
about 3–4 minutes and 3.5–4.5 GB of RAM per run on a laptop, almost all of it spent parsing
the reference. Run comparisons one at a time on machines with ≤ 16 GB RAM.

### Quick start: any prediction against any reference

```bash
annotation-qc pairwise-compare \
  --query predictions.gtf \
  --reference reference.gff3 \
  --genome genome.fa \
  --evaluation-mode protein_coding \
  --reference-transcript-biotypes protein_coding \
  --outdir comparison
```

Then read, in this order:

1. `comparison_summary.json` → `sensitivity_cds.cds_coordinate_exact_count` /
   `total_reference_genes`: reference coding genes whose CDS the query reproduces exactly
   (start, every splice site and stop identical) for at least one eligible isoform.
2. `specificity.consensus_cds_coordinate_exact_count` / `total_consensus_genes`: the same
   from the query side ("consensus" always means the query).
3. `cds_exact_one_to_one`: one-to-one coordinate-exact CDS precision, recall and F1
   (each reference and query gene used once), and `cds_overlap_locus_recovery`: reference
   loci with any same-strand CDS overlap. See
   [One-to-one exact-CDS F1 and locus recovery](#one-to-one-exact-cds-f1-and-locus-recovery).
4. `sensitivity_cds.cds_intron_chain_recovered` / `multi_segment_cds_reference_genes`,
   `intron_support.cds`, `split_merge`, `novel_categories`.

Most ab initio predictors annotate CDS only (no UTRs). For them, use the CDS metrics as the
headline. Exon-level `Exact_Match` and `exon_coordinate_exact` also require the UTRs to
agree, so for these predictors they stay near zero by construction.

### Tested example: human and mouse Ensembl references

This exact sequence was run on Ensembl GRCh38.p14 / GRCm39 references (genebuild
2026-04) with Tiberius, ANNEVO, Helixer and Vipsania predictions. Shown for mouse; human
is identical with the GRCh38 files.

```bash
# 1. Annotations: decompress or name gzip files *.gz. A gzip/BGZF annotation named *.bgz fails
#    (see "Input formats"). The FASTA may stay BGZF/gzip compressed.
# 2. Current Ensembl GFF3 uses transcript types the parser does not list yet
#    (unconfirmed_transcript for TEC, gene_segment). Retype them; IDs and Parent links are kept.
gzip -dc mouse_GRCm39.gff.bgz \
  | awk 'BEGIN{FS=OFS="\t"} !/^#/ && ($3=="unconfirmed_transcript" || $3=="gene_segment") {$9=$9";original_type="$3; $3="transcript"} {print}' \
  > mouse_GRCm39.gff3
# 3. Evaluation scope: the sequences the reference annotates (excludes e.g. alt haplotypes
#    and patches that a predictor may have annotated but the reference does not).
grep '^##sequence-region' mouse_GRCm39.gff3 | awk '{print $2}' > reference_sequences.txt
# 4. Compare. Repeat with identical options for every predictor.
annotation-qc pairwise-compare \
  --query mouse_GRCm39_tiberius.gtf \
  --reference mouse_GRCm39.gff3 \
  --genome mouse_GRCm39.fa.bgz \
  --evaluation-mode protein_coding \
  --reference-transcript-biotypes protein_coding \
  --regions-file reference_sequences.txt \
  --outdir mouse_tiberius_all_isoforms
# 5. Optional second view: one representative isoform per gene on both sides.
#    Add  --reference-transcript-selection longest_cds --query-transcript-selection longest_cds
```

Comparing several predictors: keep the reference, `--genome`, map, mode, biotype filters,
selection options and regions identical. Numbers are only comparable between runs with
identical options.

### GMB example (optional use case)

```bash
# Run after GMB preflight → build → finalise.
annotation-qc pairwise-compare \
  --query "$OUT/finalise/consensus.gff3" \
  --reference "$REFERENCE_GFF3" \
  --genome "$GENOME_FASTA" \
  --evaluation-mode protein_coding \
  --outdir "$OUT/reference_comparison"
```

- GMB's own hard QC (FASTA/UTR/boundary checks) asks whether an annotation is *internally
  consistent*. This command asks whether it *agrees with a reference*. Passing one says
  nothing about the other, and nothing computed here feeds back into GMB.
- **All isoforms (default):** `finalise/consensus.gff3` holds every GMB transcript.
- **Canonical output:** use `finalise/canonical/consensus.canonical_annotated.gff3`
  **with `--query-transcript-selection canonical`**. That file still contains every
  isoform and only *tags* the canonical one. Without the option, the result is the same as
  the all-isoform run.
- **Backbone only:** run again with the ab initio backbone (e.g. the Tiberius GTF) as
  `--query` and every other option unchanged.
- When the query sits in a GMB `finalise/` directory, the manifest also records the GMB
  `handover_manifest.json` and checks the query's SHA-256 against it. For other inputs,
  `gmb_finalise` is `null`.

### Input formats, compression and identifiers

- **Formats** are detected from content, not from the file name: GFF3 (`ID`/`Parent`), GTF
  (`gene_id`/`transcript_id`; gene and transcript lines may be missing and are then
  inferred from exons) and the Tiberius GTF hybrid (bare IDs on gene/transcript lines).
  Override with `--query-format` / `--reference-format`.
- **Compression:** annotations may be plain or gzip/BGZF **named `*.gz`**. A compressed
  annotation with another suffix (e.g. `*.gff.bgz`) currently fails with a
  `UnicodeDecodeError`. The FASTA given to `--genome` may be plain, gzip or BGZF under any
  name. No index is written next to it.
- **Feature types used:** genes (`gene`, `ncRNA_gene`, `pseudogene`), transcripts
  (`mRNA`, `transcript`, the RNA types including `Y_RNA`, `*_gene_segment`,
  `pseudogenic_transcript`, `processed_transcript`), `exon` and `CDS`. Other types (UTRs,
  introns, codons) are ignored. A transcript type outside this list (Ensembl's
  `unconfirmed_transcript`, `gene_segment`) leaves its exons without a parent, and the run
  stops: retype it as in the tested example.
- **CDS-only files.** A transcript with CDS lines but no exon lines (e.g. a GTF holding
  only `gene` and `CDS` lines) gets one exon per CDS line, with the same coordinates
  (`source_feature` `inferred_exon_from_CDS`). The run log and the manifest's
  `parse_diagnostics` report `cds_only_transcripts` and `inferred_exon_rows`. Transcripts
  with at least one exon line are used as written. Exon metrics then describe coding exons
  only, as for any predictor without UTRs.
- **Exon lines that omit UTRs.** Some predictors write `exon` lines for the coding part
  only and describe UTRs solely with `five_prime_UTR`/`three_prime_UTR` lines (transcript
  and gene lines then extend beyond the exons). UTR lines are not read, so such files are
  compared as written: exon metrics see no UTR, while gene spans (used for locus pairing)
  still include it.
- **CDS convention:** CDS coordinates are compared as written. Ensembl GFF3, Tiberius,
  ANNEVO, Helixer, Vipsania and OrionGene include the stop codon in the CDS. If one file
  excludes it (e.g. GENCODE/GTF with separate `stop_codon` lines, or a CDS-only GTF whose
  gene lines end 3 bp after the last CDS), every complete gene differs by 3 bp and is never
  coordinate-exact, while CDS `Exact_Match` (same intron chain) is unaffected. Check the
  three genomic bases after each CDS end before interpreting terminal differences, and
  keep any stop-harmonised re-run separate from the comparison of original coordinates.
  `--add-stop-codon` performs that re-run for the query (see
  [Optional stop-codon harmonisation](#optional-stop-codon-harmonisation---add-stop-codon)).
- **Identifiers.** Identity fields are `ID`/`Parent` (GFF3) and `gene_id`/`transcript_id`
  (GTF). Display names (`Name`, `gene_name`) are never used for linking, so repeated
  names are harmless. An identity value reused on **different sequences** (e.g. Tiberius
  numbering `g1, g2, …` per chromosome) is namespaced as `seqname:id` in every occurrence,
  with the file value kept in `original_gene_id`. Unique IDs are left unchanged. An
  identity value reused on the **same sequence** stops the run, as do children without a
  resolvable parent. Exons/CDS with several parents (`Parent=t1,t2`) are counted for each
  parent. Not detected: verbatim duplicated exon/CDS lines (counted twice, which breaks
  coordinate-exact and intron-chain results for that transcript) and, in a GTF *without*
  transcript lines, one `transcript_id` reused at two loci on one sequence (silently fused
  into one transcript). Deduplicate such files first.
- **Sequence names** must match the FASTA. Use `--seqname-map` (or `--assembly-report`)
  otherwise. Matching names and in-bounds coordinates do not prove the same assembly.
  Check the assembly accession of every input.

### Optional stop-codon harmonisation (`--add-stop-codon`)

Some annotations end the CDS before the stop codon (a *stop-excluded* convention);
Ensembl and most predictors include it. Compared as written, every complete model then
differs by 3 bp: coordinate-exact CDS, one-to-one TP, precision, recall and F1 drop to
about 0, while locus recovery and CDS intron chains are unaffected. `--add-stop-codon`
re-runs the comparison with the **query** converted to the stop-included convention.
The reference is never changed. Without the option (default) the original coordinates
are compared; report harmonised runs separately and label them as such. Harmonisation
aligns a convention: it does not show that an ORF is valid, and it repairs nothing.

```bash
annotation-qc pairwise-compare --query cds_only.gtf --reference reference.gff3 \
  --genome genome.fa --evaluation-mode protein_coding \
  --add-stop-codon --outdir comparison_stop_harmonised
```

**What changes.** Only the terminal CDS segment of a transcript, extended by 3 bp in
transcript orientation (+ strand: end + 3; − strand: start − 3), when the next three
genomic bases, read in transcript orientation, are a stop codon of the sequence's genetic
code. Exons inferred from CDS (CDS-only files) and transcript/gene records inferred from
their children are extended with it. Identifiers, parents, other CDS segments and all
explicit records stay as written. Rules depend only on coordinates and sequence, never on
tool names, file names or identifiers.

**Frame evidence (eligibility).** CDS phase is not used, so missing phase is fine.
The frame is taken to start at the first CDS base and must be supported by the sequence.
Checks, in order, with the first failure reported as the reason:
1. sequence present in the genome (`sequence_not_available`);
2. a genetic code assigned (`genetic_code_not_specified`);
3. `+` or `−` strand (`strand_not_defined`);
4. no overlapping CDS segments (`overlapping_cds_segments`);
5. CDS inside the sequence (`cds_outside_sequence`);
6. only A/C/G/T in the CDS (`ambiguous_bases_in_cds`);
7. CDS length a multiple of 3 (`cds_length_not_multiple_of_3`);
8. first codon ATG (`no_start_codon`): 5′-partial models and non-ATG starts are not
   guessed;
9. no in-frame stop before the last codon (`internal_stop_codon`).

An adjacent stop triplet after a disrupted frame is therefore never used.

**Outcomes for eligible models.**

| status | reason | meaning |
|---|---|---|
| unchanged | `stop_already_in_cds` | last codon is a stop. A second application is therefore a no-op. |
| unchanged | `no_adjacent_stop` | the next three bases are not a stop codon |
| skipped | `split_stop_codon_unsupported` | the CDS ends at (or within 2 bp of) its exon end and another exon follows. A stop split by an intron is not added. |
| skipped | `extension_beyond_explicit_exons`, `cds_outside_explicit_exons` | explicit exon records would not contain the extended CDS. Exons are never rewritten. |
| skipped | `extension_out_of_bounds`, `ambiguous_bases_after_cds` | the three bases are past the sequence end or not A/C/G/T |
| skipped | `extension_outside_transcript_record`, `extension_outside_gene_record` | an explicit transcript/gene record does not contain the extended CDS |
| extended | `stop_codon_added` | terminal CDS (and inferred records) extended by 3 bp |

**Genetic code.** `--stop-codon-genetic-code` sets the NCBI table for all sequences
(default 1, standard; supported: 1–6, 9–14, 16, 21–26, 29, 30, 33). Tables whose stops
are context-dependent (27, 28, 31) are not supported. `--stop-codon-sequence-code
MT=2` sets one sequence's table, and `NAME=none` skips it. Sequences named like organelle
genomes (`MT`, `M`, `chrM`, `chrMT`, `mito`, `mitochondrion`, `Pt`, `chrC`, `chrPt`,
`chloroplast`, `plastid`; case-insensitive) never receive the default table: without an
override their models are skipped (`genetic_code_not_specified`).

**Outputs and provenance.**
- `stop_codon_harmonisation.tsv`: every query transcript with CDS (status, reason,
  genetic code and its source, last and next codon, original and new terminal CDS
  coordinates, records changed).
- `query_evaluated.stop_harmonised.gtf`: the query as compared. Locus plots
  (`--plots-per-category`) also draw these coordinates.
- `comparison_summary.json` → `query_stop_codon_harmonisation`: `applied` (`false` in
  default runs), counts by status and reason, and the policy.
- `comparison_manifest.json` → `stop_codon_harmonisation`: the same, plus the SHA-256 of
  both output files. The query file itself is recorded, unchanged, under `inputs.query`.

Harmonisation runs before query transcript selection, so `longest_cds` sees the
harmonised lengths.

### Reference scope and isoform policy

- **Scope.** Pairing is restricted to whatever both annotations contain. A query gene on a
  sequence the reference does not annotate (alternative haplotype, patch, unplaced
  scaffold) is `Novel`, and a reference gene on a sequence the predictor skipped is
  `Missed`. Choose the scope before comparing tools, independently of their results. A
  whole-sequence regions file (one seqname per line) is the simplest way.
- **Reference genes and isoforms.** `--evaluation-mode protein_coding` keeps
  protein-coding genes with *all* their transcripts (including NMD and retained-intron
  isoforms). Add `--reference-transcript-biotypes protein_coding` to keep only
  protein-coding isoforms. Genes that then have no isoform are dropped, which changes
  `total_reference_genes`.
- **Isoforms per gene.** With `--reference-transcript-selection all` (default), a reference
  gene counts as CDS-exact when *any* kept isoform is reproduced. Each gene is counted
  once, whatever the number of matching isoforms or query transcripts. With modern
  references (often > 10 coding isoforms per gene) this is much more lenient than scoring
  one representative. `longest_cds` and `canonical` select one isoform per gene (longest
  total CDS, ties → first in file; `Ensembl_canonical` tag, falling back to longest CDS)
  independently of the query. A correct alternative isoform then counts as a miss.
  `--query-transcript-selection` applies the same rules to the query and is a no-op for
  predictors that emit one transcript per gene.

### Arguments

| argument | meaning |
|---|---|
| `--query` (required) | annotation to evaluate: GFF3, GTF or Tiberius hybrid; plain or gzip named `*.gz` |
| `--reference` (required) | trusted reference: GFF3 or GTF; plain or gzip named `*.gz` |
| `--outdir` (required) | output directory |
| `--genome` | genome FASTA (plain, gzip or BGZF). Checks that every sequence exists and every feature lies within it. Recommended. |
| `--seqname-map` | two-column TSV/CSV map, see below |
| `--assembly-report` | NCBI assembly report (GenBank accession → assigned molecule). `--seqname-map` entries win on conflict. |
| `--genome-mismatch {error,warn,exclude}` | what to do when, with `--genome`, an annotation sequence is missing from the FASTA or has features past its end. `error` (default) stops. `exclude` removes those sequences from **both** annotations and records them in the audit and manifest. `warn` keeps them, which is the old `gmb-compare` behaviour; their reference genes then count as Missed. |
| `--evaluation-mode {all,protein_coding,cds_only,canonical}` | reference preset (default `all`). `protein_coding` keeps genes with `gene_biotype`/`biotype` `protein_coding` and all their transcripts. |
| `--reference-transcript-selection {all,longest_cds,canonical}` | reference transcripts per gene (default `all`). `canonical` takes the `Ensembl_canonical`-tagged transcript, falling back to longest CDS. |
| `--reference-gene-biotypes`, `--reference-transcript-biotypes` | extra comma-separated biotype filters |
| `--query-transcript-selection {all,longest_cds,canonical}` | query transcripts per gene (default `all`) |
| `--query-format`, `--reference-format` | `auto` (default: detected from content), `gff3`, `gtf`, `tiberius` |
| `--region`, `--regions-file` | restrict both annotations (and the Novel context) to whole genes overlapping `seqname[:start-end]` (one-based, inclusive, comparison seqnames). A line with only a seqname keeps the whole sequence. `--region` can be repeated. |
| `--evidence-attribution` | GMB `build/evidence_attribution.tsv`; writes it back with each transcript's label |
| `--plots-per-category N` | write N locus plots per class to `qc/` with `qc/index.html` (needs matplotlib) |
| `--add-stop-codon` | query only, needs `--genome`: add an adjacent in-frame stop codon to CDS that end before it (below). Default: original coordinates |
| `--stop-codon-genetic-code TABLE` | NCBI translation table for `--add-stop-codon` (default 1) |
| `--stop-codon-sequence-code SEQNAME=TABLE` | per-sequence table (comparison seqname), `none` to skip a sequence; repeatable |

**Seqname map direction.** Column 1 is the name used in an annotation file and column
2 is the name to compare under (normally the FASTA header), e.g. `Pf3D7_01_v3	1`. The
same map is applied to both annotations. Unlisted names are left unchanged. A header
row (`from_seqname	to_seqname`) and `#` comments are allowed. A map that merges two
sequences into one name is rejected. If the reference and query share no sequence
name after mapping, the command stops rather than reporting every gene as Missed.

### Outputs

| file | content |
|---|---|
| `comparison_summary.json` / `.tsv` | headline metrics (below) |
| `comparison_details.tsv` | one row per reference gene and per query gene. `source` is `reference` or `consensus`. Coordinates are one-based inclusive. |
| `consensus_transcript_labels.tsv` | one row per query transcript: Matched / Strand_Mismatch / Novel, plus the best CDS reference transcript and `novel_category` |
| `gene_splits.tsv` / `gene_merges.tsv` | reference genes with ≥ 2 query counterparts / query genes with ≥ 2 reference counterparts, with counterpart IDs and one-based loci |
| `cds_exact_one_to_one_pairs.tsv` | the gene pairs of the one-to-one exact-CDS matching (columns below) |
| `stop_codon_harmonisation.tsv` | only with `--add-stop-codon`: one row per query transcript with CDS (status, reason, genetic code, codons, original/new terminal CDS, one-based) |
| `query_evaluated.stop_harmonised.gtf` | only with `--add-stop-codon`: the query exactly as compared (after harmonisation, transcript selection and regions), for plots and browsers |
| `reference_filter_audit.json` / `.tsv` | reference counts before and after filtering, biotypes, excluded sequences |
| `evidence_attribution_labeled.tsv` | only with `--evidence-attribution` |
| `comparison_manifest.json` | query/reference/genome/map paths **with SHA-256**, all options, parse diagnostics, seqname checks, code version, and the GMB finalise manifest and build directory when the query is a GMB output (`query_matches_manifest` confirms the checksum) |

File and key names are kept from `gmb-compare` so existing consumers keep working. In
them, **"consensus" always means the query annotation** (`total_consensus_genes`,
`consensus_classification`, `novel_consensus_count`,
`consensus_transcript_labels.tsv`).

### Interpreting the metrics

| question | field(s) | denominator |
|---|---|---|
| reference coding genes with an exactly reproduced CDS (start, splice sites, stop) | `sensitivity_cds.cds_coordinate_exact_count` | `total_reference_genes` |
| one-to-one exact-CDS precision, recall, F1 | `cds_exact_one_to_one.precision`, `.recall`, `.f1` (TP = `.true_positives`) | Q, R, R + Q |
| reference loci recovered by same-strand CDS overlap | `cds_overlap_locus_recovery.recovered_count`, `.rate` | `total_reference_genes` |
| query genes whose CDS exactly equals an eligible reference isoform | `specificity.consensus_cds_coordinate_exact_count` | `total_consensus_genes` |
| same CDS intron chain, terminal (start/stop) coordinates may differ | `sensitivity_cds.cds_exact_match_count`; `cds_exact_match_not_coordinate_exact` | `total_reference_genes` |
| CDS intron chain recovered (best pair) | `sensitivity_cds.cds_intron_chain_recovered` | `multi_segment_cds_reference_genes` (`_sensitivity`) or detected multi-segment genes (`_rate`) |
| unique CDS introns recovered / supported | `intron_support.cds` | `reference_introns` / `query_introns` |
| locus detected (same-strand **gene-span** overlap) | `sensitivity.locus_detected_count`; `locus_detection_exonic` | `total_reference_genes` |
| missed / unmatched predictions | `sensitivity.missed_count`; `specificity.novel_consensus_count`, `novel_categories` | R / Q |
| splits, merges, strand mismatches | `split_merge`, `sensitivity.strand_mismatch_count`, `specificity.strand_mismatch_consensus_count` | counts |
| UTR-inclusive agreement | `sensitivity.exact_match_count` (structural), `exon_coordinate_exact_count` | `total_reference_genes` |

Every reference gene gets one class, and every query gene gets one class:

| class | meaning |
|---|---|
| Exact_Match | *structural* exact: best transcript pair has ≥ 0.8 reciprocal exonic overlap **and** an identical intron chain. Terminal coordinates may differ. |
| Structural_Mismatch | ≥ 0.8 reciprocal overlap, different intron chain |
| Partial_Match | same-strand gene spans overlap, but the pair is below 0.8 overlap (can be 0 exonic overlap, see below) |
| Strand_Mismatch | no same-strand partner, and an opposite-strand gene shares **CDS** bases (exon bases when either gene lacks CDS). `strand_mismatch_basis` says which. |
| Missed / Novel | no same-strand partner and no opposite-strand feature overlap (reference / query). Opposite-strand contact through spans, introns or UTRs only is recorded as `strand_mismatch_basis = no_feature_overlap`. |

The CDS class applies the same rules to the **best CDS pair**, which is chosen
independently of the best exon pair (`cds_matched_id`, `best_cds_match_transcript_id`).
`No_CDS` means none of the gene's transcripts has CDS. For a reference gene, CDS
`Missed` means no query partner transcript has CDS.

Coordinate-exact columns: `cds_coordinate_exact` is True only when the best CDS pair
has identical CDS intervals, start and stop included; `exon_coordinate_exact` is the
same for exons. CDS `Exact_Match` without `cds_coordinate_exact` usually means a
different start codon.

Novel query genes keep the `Novel` class (they never count as matches) and get a
`novel_category` from exon overlap, on either strand, with the reference **before**
mode/biotype filtering: `Overlaps_reference_pseudogene`,
`Overlaps_other_non_protein_coding_reference`,
`Overlaps_unevaluated_protein_coding_reference`,
`Overlaps_evaluated_reference_opposite_strand` or `Novel_no_reference_overlap`.

Splits and merges: a reference and a query gene are *counterparts* when they are
same-strand locus partners and some transcript pair shares ≥ 10% of the shorter CDS
(exons when either gene lacks CDS). A split is a reference gene with ≥ 2 query
counterparts; a merge is a query gene with ≥ 2 reference counterparts.

- **Denominators.** `rate_denominators` in `comparison_summary.json` lists them.
  Reference-based (sensitivity-like) rates divide by the reference genes after
  filtering (`total_reference_genes`). Query-based (precision-like) rates in
  `specificity` divide by `total_consensus_genes`. `intron_support` gives reference
  introns recovered / reference introns and query introns supported / query introns.
- **Intron chains** apply only to structures with an intron. `intron_chain_match` and
  `cds_intron_chain_match` are `NA` for single-exon / single-CDS-segment structures.
  `exon_intron_chain_sensitivity` divides by all multi-exon reference genes;
  `exon_intron_chain_rate` (and `exon_intron_chain_matched`) uses detected reference
  genes whose compared transcript has an intron. The CDS keys work the same way.
- `matched_consensus_count` is the number of query genes with a same-strand reference
  partner (Exact/Partial/Structural/Matched), each counted once.
- **CDS coordinate-exact** is the strictest measure of coding accuracy. **Exact_Match**
  also requires UTR extent to agree, so it drops when the reference has longer UTRs
  than the query.
- **Locus detection** counts Exact + Structural + Partial. Pairing uses gene spans, so a
  query gene wider than its transcripts can "detect" a reference gene without sharing
  any exon with it. `locus_detection_exonic` gives the count with real exonic overlap
  and the span-only difference. `gene_span_checks` lists genes whose span extends
  beyond their exons.
- **Reference completeness and selection change the numbers.** Missing or partial
  reference genes inflate Novel and deflate precision. Genes the reference doesn't have
  count against the query. Restricting the reference (`protein_coding`, `canonical`,
  `longest_cds`) changes the denominator and which transcripts can match exactly. With
  all reference isoforms, the query only needs to reproduce one of them per gene.

### One-to-one exact-CDS F1 and locus recovery

**Populations.** R = `total_reference_genes` and Q = `total_consensus_genes`: every gene
left after the evaluation mode, biotype filters, `--*-transcript-selection` and
`--region`/`--regions-file` scope (and `--genome-mismatch exclude`). Genes without a
comparable CDS (no transcript with exons and CDS, e.g. a query lncRNA) stay in R and Q.
They cannot be matched and count as unmatched, so a predictor's non-coding genes lower
its precision.

**Exact-CDS graph.** A transcript's CDS signature is (sequence, strand, all CDS segment
coordinates), compared exactly as parsed. Reference gene *r* and query gene *q* are joined
when some selected transcript of *r* and some transcript of *q* have the same signature.
Transcripts without CDS have no signature. Every transcript is considered, not only the
best pairs reported in `comparison_details.tsv`.

**Matching and metrics.** A maximum-cardinality bipartite matching (Kuhn's augmenting
paths) uses each reference and each query gene at most once:

```
TP        = number of matched gene pairs
precision = TP / Q
recall    = TP / R
F1        = 2·TP / (R + Q)            (= harmonic mean of precision and recall)
unmatched_reference_genes = R − TP,   unmatched_query_genes = Q − TP
```

A value whose denominator is 0 is `null` in JSON and `NA` in the TSV, not 0. The block also
records how many genes have any exact candidate, the number of edges, and the number of
connected components with more than one edge.

**Example with competing exact matches.** Reference gene R1 has isoforms with CDS X and Y,
reference gene R2 has CDS X, query gene Q1 has CDS X and query gene Q2 has CDS Y. The
edges are R1–Q1, R1–Q2 and R2–Q1.

- Per gene, all four genes are coordinate-exact (`cds_coordinate_exact_count` = 2,
  `consensus_cds_coordinate_exact_count` = 2).
- Both reference genes' best CDS pair names Q1, so deduplicating best pairs gives TP = 1,
  as does a greedy first-come pairing.
- The maximum matching is R1–Q2 and R2–Q1: TP = 2, precision = recall = F1 = 1.

If instead R1 and R2 had the same CDS and the query had a single gene Q with that CDS,
both reference genes would be coordinate-exact, but TP = 1: recall 0.5, precision 1.0,
F1 = 2·1/(2+1) = 0.667. This is why the historical per-gene counts can exceed TP: identical
models annotated in several genes (duplicated reference genes, or a predictor emitting
the same model twice).

**Pairs file** (`cds_exact_one_to_one_pairs.tsv`, one row per TP, coordinates one-based):

| column | meaning |
|---|---|
| `reference_gene_id`, `query_gene_id` | comparison IDs (`seqname:id` when namespaced) |
| `reference_original_gene_id`, `query_original_gene_id` | IDs as written in the files |
| `chrom`, `strand` | shared sequence and strand |
| `reference_start`/`_end`, `query_start`/`_end` | gene spans |
| `cds_start`, `cds_end`, `cds_segments` | span and segment count of the shared CDS (the lowest one if several) |
| `shared_cds_signatures` | number of distinct identical CDS between the two genes |
| `supporting_reference_transcript_ids`, `supporting_query_transcript_ids` | every transcript of each gene with one of those CDS |
| `reference_exact_candidates`, `query_exact_candidates` | number of exact partners each gene has |
| `component_reference_genes`, `component_query_genes` | size of the connected component |

**Ties.** The matching is deterministic. Reference genes are processed in
(sequence, start, end, gene_id) order, and candidates are tried in the same order. Other
maximum matchings with the same TP can exist wherever a component has more than one
edge (`component_*_genes` > 1). They can pair different genes but never change TP,
precision, recall or F1.

**CDS-overlap locus recovery.** `cds_overlap_locus_recovery.recovered_count` counts
reference genes with a same-strand query locus partner whose transcripts share at least
one CDS base. It equals the number of reference rows with `cds_overlap > 0` (the best CDS
pair's reciprocal overlap is > 0 exactly when some pair shares a base); `rate` divides by
R. `gene_span_detected_without_cds_overlap` counts genes that gene-span detection
(`sensitivity.locus_detected_count`) finds without any shared CDS. Limitations:
- Partners are found through gene-span overlap, so a gene line that does not cover its own
  CDS can miss a CDS overlap.
- `comparison_details.tsv` rounds `cds_overlap` to 4 decimals. A one-base overlap between
  very long CDS can print as `0.0`. The summary count uses the unrounded value.

**Three different questions.**

| metric | asks | strictness |
|---|---|---|
| `cds_overlap_locus_recovery` | is the coding locus found at all? | any shared CDS base, same strand |
| `cds_intron_chain_recovered` | are the CDS splice sites right? | identical CDS intron chain in the best pair; start/stop may differ; multi-segment genes only |
| `cds_exact_one_to_one` | is the whole coding model right, one prediction per gene? | identical CDS including start and stop; one-to-one |

**Sensitivity of these numbers.**
- **Reference isoform selection:** with all isoforms, a query gene matches if it reproduces
  any eligible isoform. `longest_cds` or `canonical` leave one target per gene and lower
  TP. A correct alternative isoform then counts as a miss.
- **Reference completeness:** missing or partial reference genes turn correct predictions
  into unmatched query genes (lower precision). Duplicated reference models cap TP below
  the per-gene exact counts.
- **Stop-codon convention:** an annotation that excludes the stop codon from the CDS never
  matches one that includes it. TP, precision, recall and F1 are then near 0, while CDS
  intron chains and locus recovery are unaffected. By default the comparator does not
  harmonise this. Check it (see "CDS convention" above); `--add-stop-codon` re-runs with
  the query converted, and that run must be reported separately from the
  original-coordinate run.

### Limitations to keep in mind

- **Locus pairing uses gene spans.** A same-strand span overlap pairs two genes even without
  shared exons or CDS. `locus_detection_exonic.span_only_detected_count` counts genes
  detected without exonic overlap. For a coding view, use `cds_overlap_locus_recovery`.
- **Intron chains are scored on one best pair.** The best pair ranks ≥ 0.8 reciprocal
  overlap above an identical intron chain with lower overlap. A gene split into pieces, or
  matched by a longer and a shorter prediction, can report `cds_intron_chain_match = False`
  although some other pair has the identical chain. An "any pair" count can therefore be
  slightly higher.
- **Query intron support covers the whole query in scope.** Query introns at loci outside
  the evaluated reference (pseudogenes, non-coding or unevaluated genes) count as
  unsupported.
- **Exon/intron sets of the reference grow with its isoforms.** With all isoforms,
  `intron_support.*.reference_introns` includes every isoform's introns, so the
  recovered fraction falls as the reference gets richer. It is not a per-gene measure.
- **F1 is gene-level and exact-CDS only.** `cds_exact_one_to_one` is the only F1. Do not
  combine `cds_coordinate_exact_count` and `consensus_cds_coordinate_exact_count` into an
  F1: they count genes independently and can both exceed the one-to-one TP. No exon-level
  (CDS segment) precision/recall is computed.
- **Splits/merges** use a fixed 10 % shared-CDS threshold (not a CLI option).
- **Novel is relative to the reference.** Predictions at loci the reference lacks or
  annotates as pseudogenes are Novel. `novel_category` says which, but they still lower
  query-based rates.
- **Same-strand partners take precedence.** A gene with any same-strand partner is never
  Strand_Mismatch, even with opposite-strand CDS overlap.

### Coordinates and parsing

Annotations are parsed with pyranges1 and converted once, at the parser boundary, to
zero-based half-open intervals (`Start = GFF start − 1`, `End = GFF end`). Overlap
fractions use true lengths, and features sharing a single base count as overlapping.
Reports convert back to one-based inclusive. Gene → transcript → exon/CDS links come
from explicit `Parent`/`gene_id`/`transcript_id` values. Children with several parents
are counted once per parent. Unresolvable links (orphan exons, transcripts without a
gene, duplicate gene IDs on one sequence) stop the run. IDs reused on different
sequences (Tiberius `g1`, `g2`, … per chromosome) are namespaced as `seqname:id`, and
the original IDs are kept in `original_gene_id`.

### Differences from `gmb-compare`

Validated on P. falciparum and T. gondii GMB builds (every gene row compared):

1. `gmb-compare` kept GFF one-based starts in half-open arithmetic, so each exon was
   1 bp short and genes sharing only their boundary base did not pair. Correcting this
   moves overlap values by at most ±0.004, reclassifies a handful of genes near the 0.8
   threshold or at touching boundaries, and never changed an Exact_Match or CDS-exact
   count in validation.
2. When a reference gene ties between query genes with byte-identical spans, the
   earlier file record wins. `gmb-compare`'s choice there depended on an unstable sort.
3. A sequence absent from the genome, or a feature past its end, stops the run by
   default (`--genome-mismatch`). `gmb-compare` only warned.
4. Regions keep whole genes. `gmb-compare` cut exons at the region edge.
5. Multi-parent exons/CDS are counted for every parent. `gmb-compare` dropped them.
   GTF files without gene/transcript lines get them inferred from their exons, and
   CDS-only transcripts get exons from their CDS lines.
6. Audit biotype counts are integers, not strings. Plots show the reference and query
   only; evidence-track plots, `tool_performance_analysis` and `validate_annotation` are
   not ported yet.
7. Correctness and reporting fixes from the P. falciparum Vipsania/Tiberius
   evaluation (pairing thresholds unchanged): `matched_consensus_count` no longer
   counts `Matched` query genes twice; the CDS pair is chosen independently of the exon
   pair; Strand_Mismatch needs CDS (or exon) overlap, not gene-span contact; single-exon
   intron chains are `NA` rather than a match; coordinate-exact, split/merge,
   Novel-context and query-based rates are new. The `comparison_details.tsv` columns
   up to `original_gene_id` are unchanged, and new columns are appended after it.

---

# Annotation QC onboarding guideline

This guide explains how to structure code inside `annotation_qc/` 

## 1. Core rule

The package follows a single flow:

`input file -> parser -> standardized pandas DataFrame -> runner -> metric functions -> report output`

The important rule is:

* **parsers** normalize input;
* **metrics** compute QC logic;
* **runners** connect everything for CLI use;
* **reports** format results for export;
* **metadata** provides external reference data when needed.

Do not create a second annotation object model unless there is a strong reason.

## 2. Folder responsibilities

### `parsers/`

Use this folder for code that reads files and converts them into a standard internal representation.

Put here:

* GFF3/GTF parsing
* FASTA parsing
* DIAMOND/STAR/BUSCO/prediction TSV parsing
* attribute normalization
* ID cleanup
* structural normalization such as `gene_id`, `transcript_id`, `Parent`, `biotype`

Do not put here:

* QC decisions
* statistics calculation
* report writing
* CLI parsing

### `metrics/`

Use this folder for pure QC calculations.

Put here:

* feature-level metrics
* gene/transcript/exon summary calculations
* distribution checks
* classification logic
* tolerance-based scoring
* parity-safe refactors of old scripts

Rules for this folder:

* accept DataFrames or already parsed evidence objects
* return dicts or DataFrames
* avoid file I/O
* avoid CLI logic
* avoid parsing raw file formats again

### `runners/`

Use this folder for workflow orchestration.

Put here:

* CLI entry points
* input selection
* filtering slices from the parsed table
* calling metric functions in the right order
* passing results to reports

A runner should be the place where the workflow is assembled, not where the QC logic lives.

### `reports/`

Use this folder for turning metric results into output.

Put here:

* TSV export
* JSON export
* text summaries
* plots
* tables
* any reusable presentation layer

Do not recompute metrics here.

### `notebook/`

Use this folder for visualising and investigating result.

Put here any reusable notebook ideally one for each type of metrics.


### `config/`

Use this folder for configuration definitions and defaults.

Put here:

* thresholds
* constants
* runtime settings
* feature lists
* shared configuration values

## 3. Practical rule for deciding where code goes

Use this decision rule:

* If it **reads a file** → `parsers/`
* If it **calculates QC values** → `metrics/`
* If it **connects steps together** → `runners/`
* If it **formats output** → `reports/`
* If it **stores constants or thresholds** → `config/`
* If it **retrieves external reference data** → `metadata/`

## 4. Adding a new QC

* When implementing a new QC:
* Identify the required input.
* Reuse an existing parser where possible.
* Implement the calculation inside metrics/.
* Call the metric from a runner.
* Export the result through a report.
* Add tests.
* Update the relevant package documentation.

## 7. Template for a new metric module

Use this template for new files inside `metrics/`:

```python
"""
feature_wide.py

Purpose:
    Compute feature-wide annotation QC metrics from a standardized DataFrame.

Inputs:
    annotation_df: pandas.DataFrame
        Parsed and normalized annotation data.

Outputs:
    dict | pandas.DataFrame
        Metric summary suitable for reporting.
"""

from __future__ import annotations

import pandas as pd


def compute_feature_wide_metrics(annotation_df: pd.DataFrame) -> dict:
    """
    Compute feature-wide QC statistics.

    Parameters
    ----------
    annotation_df
        Standardized annotation table produced by parsers/.

    Returns
    -------
    dict
        Feature-wide metric values.

    Notes
    -----
    This function must not read files, write files, or perform CLI logic.
    """
    # 1. validate required columns
    # 2. derive slices if needed
    # 3. calculate metrics
    # 4. return plain Python data
    raise NotImplementedError
```

## 8. Template for a runner

```python
from ensembl.genes.annotation_qc.parsers.annotation import parse_annotation
from ensembl.genes.annotation_qc.metrics.feature_wide import compute_feature_wide_metrics
from ensembl.genes.annotation_qc.reports.feature_wide import write_feature_wide_report


def run_feature_wide(annotation_path: str, output_path: str) -> None:
    annotation_df = parse_annotation(annotation_path)
    metrics = compute_feature_wide_metrics(annotation_df)
    write_feature_wide_report(metrics, output_path)
```


## 12. Definition of done

The refactor is complete when:

* the metric module has no file I/O
* the runner owns workflow orchestration
* the parser owns normalization
* the report layer owns formatting
* tests cover the metric logic
* the documentation says exactly where new code should go
* a new starter can add a feature without guessing folder responsibilities

## 13. One-line summary for contributors

**Put parsing in `parsers/`, calculations in `metrics/`, orchestration in `runners/`, output in `reports/`, and configuration in `config/`.**

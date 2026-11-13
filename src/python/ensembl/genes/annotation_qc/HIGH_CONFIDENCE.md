# High-confidence gene assessment

High-confidence QC follows the package's parser → DataFrame → metric → runner → report flow. It assesses one representative transcript per protein-coding gene and attaches the resulting confidence class to the gene. Evidence tables retain transcript IDs so each gene-level decision can be traced to the assessed model.

## Code layout

| Layer | Module | Responsibility |
| --- | --- | --- |
| Configuration | `config/high_confidence.py` | Flags, thresholds, output schemas, confidence policy version and diagram labels |
| Parsing | `parsers/high_confidence.py` | GTF/GFF3 normalization, FASTA access, junction/helper/alignment tables, stored evidence and model configurations |
| Metrics | `metrics/high_confidence.py` | Representative selection, CDS/protein reconstruction, intron/splice checks, homology arithmetic, classification and rescoring |
| Metrics | `metrics/high_confidence_summary.py` | Flag distributions, numeric summaries, evidence patterns, policy comparisons and validated Sankey paths |
| Reports | `reports/high_confidence.py` | TSV, FASTA, GFF3, manifests and rescored output formatting |
| Reports | `reports/high_confidence_summary.py` | Per-model/comparative tables, confidence plots and annotated evidence heatmaps |
| Reports | `reports/high_confidence_sankey.py` | Layout/rendering of precomputed method → decisions → confidence flows |
| Runners | `runners/high_confidence.py` | Coordinate the complete assessment and optional DIAMOND execution |
| Runners | `runners/rescore_high_confidence.py` | Reuse saved evidence without extracting sequences or rerunning DIAMOND |
| Runners | `runners/high_confidence_summary.py` | Summarize one model or compare several prediction methods |

Every runner registers a subcommand through `register(subparsers)` in the existing `annotation-qc` CLI. Modules use absolute package imports and support `python -m ensembl.genes.annotation_qc.runners.<module>` as well.

## Annotation boundary and coordinate conventions

`parse_high_confidence_annotation()` returns a pandas DataFrame with the framework's standard `Chromosome`, `Feature`, `Start`, `End`, `Strand`, `Frame`, `gene_id` and `transcript_id` columns. `Start` is zero-based and `End` is exclusive. Lossless repeated attributes and normalized hierarchy metadata are retained in `Attributes` and DataFrame metadata. Exported GFF3, junction TSVs and evidence coordinates remain one-based and inclusive.

The generic PyRanges annotation parser remains unchanged. This specialized adapter is needed for predictor bare-ID rows, order-independent GFF3 `ID`/`Parent` links, shared children, percent-escaped identifiers, repeated tags, CDS phase continuity and separately annotated/split stop codons. A small **private** reconstruction-record bridge (`parsers/_high_confidence_records.py`) preserves the tested frame/UTR/translation behavior. It is not a second public annotation API: parser and runner boundaries use DataFrames. The bridge is scoped to reconstruction and must not become a competing shared annotation representation.

FASTA I/O happens in the parser before sequence assessment. The metric receives an in-memory region for the current transcript, including three flanking bases for the existing bounded stop-codon normalization. It does not open a genome, annotation, junction or DIAMOND file. The entire genome is not loaded into memory.

The existing generic DIAMOND parser expects a different 13-column format. High-confidence QC therefore provides an explicitly named nine-column legacy adapter; it does not reinterpret the generic parser's output.

## Installation

The repository requires Python 3.12 or newer. Pandas and Biopython are already base dependencies. Install the optional plotting extra when rendering figures:

```bash
python -m pip install -e '.[annotation-qc-plots]'
```

DIAMOND is an external executable and is required only for `--run-diamond`. Scoring from existing helper results or alignment tables does not launch it. Plotting imports Matplotlib lazily; table-only workflows do not require it.

## Assess an annotation

```bash
annotation-qc high-confidence \
  --annotation annotation.gff3 \
  --genome genome.fa \
  --junctions junctions.tsv \
  --diamond-results species.prots.res \
  --no-utr \
  --output-prefix results/species
```

GTF, GFF3 and compressed annotations are accepted. The genome must be plain, regularly wrapped FASTA. The old `--gtf` and `--no_utr` spellings remain aliases. Defaults for canonical selection, automatic normalization, translation exceptions, mitochondrial contigs and coordinate-error policy are retained.

Other homology modes are mutually exclusive:

* `--diamond-alignments alignments.tsv`: nine fields, no header: `qseqid sseqid qstart qend slen pident nident length gaps`.
* `--run-diamond --protein-db reference.dmnd`: checked native DIAMOND execution using the same fields and masking option.
* `--diamond-script get_diamond_coverage.py --protein-db reference.dmnd`: explicit legacy-helper compatibility mode. Native execution is preferable because the original external helper builds an unquoted shell command.

The seven-column helper result format remains headerless: `query_id hit_id error perc_id score query_coverage target_coverage`. IDs must exactly match selected transcript IDs. Missing homology results in low confidence. Supplied values must be finite; duplicate/unknown query IDs are rejected.

Junction TSVs have a header: `seqname start end strand read_count`. Coordinates describe the intron itself, one-based and inclusive. Every CDS junction must meet `--min-junction-reads` for the support flag to pass. Without a junction file, spliced models have unknown support; absent rows in a supplied file mean no observed support. Intronless models pass support vacuously, preserving the existing convention.

## Confidence policy

`confidence_policy_version` is **2**. The default UTR-aware policy preserves the original ordered Perl branches. With `--no-utr`, the UTR-only high route is disabled and **every high-confidence assignment requires intron support**. High confidence then requires a complete CDS, a homology hit, both intron checks, support and error ≤10. Strong homology without support may qualify as medium if at least one intron check passes. Unknown support cannot satisfy support-dependent rules.

The CDS-complete flag preserves the original ATG start and TAA/TAG/TGA endpoint heuristic. Frame diagnostics, internal stops and ambiguous residues are reported but do not silently introduce new confidence thresholds. Standard and configured vertebrate mitochondrial translation, including coordinate-annotated selenocysteine, retain the existing reconstruction rules.

Homology arithmetic is deliberately unchanged for parity: target coverage uses the **query alignment span / target protein length**, and can exceed 100%. Error is `abs(query_coverage - 100) + abs(target_coverage - 100)`. This is not target-coordinate coverage. Correcting the approximation would require a separate, explicit scoring-policy change. Best-hit selection, first-hit ties and two-decimal rounding are preserved.

## Outputs and provenance

The existing suffixes remain compatible with downstream workflows: `.evidence.tsv`, `.introns.tsv`, `.cds.fa`, `.prots.fa`, `.annotations.gff3`, `.normalization.tsv`, `.input_validation.tsv`, `.translation_exceptions.tsv` and `.manifest.json`. Native homology execution also writes `.diamond.tsv` and `.prots.res`.

The manifest is written last as a completion marker. It includes input hashes, options, policy version, support requirement, canonical-selection policy and hashes of all implementation layers. Coordinate and translation failures write diagnostics without a success manifest. Rerunning an assessment overwrites its generated outputs; rescoring requires a fresh destination.

## Rescore stored evidence

```bash
annotation-qc high-confidence-rescore \
  --input-root original_results \
  --output-dir results_rna_required_v2 \
  --no-utr
```

The batch layout is `ROOT/SPECIES/METHOD/*.evidence.tsv`; actual file prefixes are retained. A single run uses `--evidence FILE` instead of `--input-root`. Rescoring updates evidence/rule columns, versioned metadata and generated GFF3 `gene_quality` attributes when a sibling exported GFF is present. Unchanged sequence/junction companions are linked to the original files, so retain the source directory. Original files are preserved and output collisions are rejected. Neither sequence extraction nor DIAMOND is rerun.

## Summaries and optional reports

```bash
annotation-qc high-confidence-summary \
  --evidence results/species.evidence.tsv \
  --label Tiberius \
  --species canis_lupus_familiaris \
  --output-dir reports/tiberius \
  --plots
```

For comparisons, replace `--evidence` with `--models models.tsv`. The configuration has `model` and `evidence` columns and optional `manifest`/`introns` columns. Paths are relative to the TSV. Models must have unique safe labels. Summaries retain unknown evidence separately from failure and use each model's own denominator.

Machine-readable outputs preserve the established `transcripts.tsv`, `overview.tsv`, flag/category/numeric distributions, evidence patterns, flag pairs, confidence strata, policy transitions and provenance schemas. They can be consumed by the existing multi-species/presentation analysis scripts. The repository report layer additionally provides confidence distributions, annotated evidence heatmaps and both sequential/matched-rule Sankeys. Sankey routes are calculated and validated in metrics before rendering; outdated stored labels must be rescored rather than relabelled by a plotting function.

## Tests and checks

```bash
python -m pip install -e '.[cicd,annotation-qc-plots]'
python -m pytest tests/ensembl/genes/annotation_qc -q
```

Tests cover the public DataFrame boundary, GFF3 hierarchy, both strands, split codons, stop-codon inclusion, RNA-required high confidence, missing support, DIAMOND parity/ties, coordinate failures, output consistency, rescoring, summary denominators and shared-CLI registration. Plot tests run when the optional dependency is installed. An architecture check prevents file/CLI operations from entering metric modules.

## Final EAM gene confidence CSV

After scoring or rescoring, create the final two-column CSV:

```bash
annotation-qc high-confidence-gene-summary \
  --evidence results/species.evidence.tsv \
  --output results/species.gene_confidence_summary.csv
```

The input must contain `gene_id` and `confidence` columns. The command counts
**unique genes present in the scored evidence**. This excludes annotation genes
that were not selected or were skipped during scoring. It does not count every
transcript or unassessed gene in the original GTF. Identical gene/class rows count
once; conflicting classes for one gene, missing gene IDs and invalid classes are
rejected. No scoring or DIAMOND alignment is rerun.

The CSV contains these seven metrics in this order:

```csv
metric,value
eam.total_genes,100
eam.total_high_confidence_genes,60
eam.total_high_confidence_genes_pct,60.00
eam.total_medium_confidence_genes,30
eam.total_medium_confidence_genes_pct,30.00
eam.total_low_confidence_genes,10
eam.total_low_confidence_genes_pct,10.00
```

Percentages use `100 × category_gene_count / eam.total_genes`, rounded to two
decimal places. A header-only input produces zero counts and `0.00` percentages;
individually rounded percentages may sum to 99.99 or 100.01. Calculation lives in
`metrics/gene_confidence_summary.py`, CSV formatting in
`reports/gene_confidence_summary.py` and CLI orchestration in
`runners/gene_confidence_summary.py`.

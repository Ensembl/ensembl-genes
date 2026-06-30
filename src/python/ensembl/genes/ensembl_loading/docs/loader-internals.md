# Loader Internals

## Generic Feature Flow

The loader converts GFF rows into a normalized in-memory annotation before any
database inserts happen.

```text
GFF3/GTF
  -> parse_converted_gff3()
  -> ParsedAnnotation
  -> reconcile_annotation()
  -> compute_exon_phases()
  -> apply_biotype_overrides()
  -> database inserts
```

### Attribute Parsing And ID Normalization

`parse_gff3_attributes()` converts the ninth GFF column into a dictionary:

```text
ID=gene1;Name=GeneOne;gene_biotype=protein_coding
```

becomes:

```python
{
    "ID": "gene1",
    "Name": "GeneOne",
    "gene_biotype": "protein_coding",
}
```

`normalize_id()` then strips each configured prefix from IDs. The generic config
does not strip anything. The RefSeq config strips `gene-`, `rna-`, `cds-`, and
`exon-`.

`parent_id()` applies the same normalization to `Parent`.

### Gene Rows

A GFF row with `type == gene` becomes a `GeneRecord`.

The loader sets:

`seq_name`
: GFF column 1 after any upstream conversion.

`start`, `end`
: GFF columns 4 and 5.

`strand`
: `1` for `+`, `-1` otherwise.

`stable_id`
: Normalized `ID`.

`name`
: First configured gene name attribute that exists. By default this is `Name`
then `gene`. If neither exists, the normalized gene ID is used.

`xref_geneid`
: The value after the configured `gene_xref_prefix` inside `Dbxref`, if present.
The current insert path stores the value in memory only; it does not populate
`xref` or `display_xref_id`.

`biotype`
: Resolved by `resolve_biotype()`, then possibly rewritten by
`apply_biotype_overrides()`.

### Transcript Rows

A GFF row becomes a `TranscriptRecord` if its feature type is present in
`source_config.parsed_transcript_feature_types`.

The loader sets:

`gene_id`
: Normalized `Parent`.

`stable_id`
: First configured transcript stable ID attribute that exists. By default this
is `Name`. If it is missing, the normalized transcript ID is used.

`biotype`
: Resolved from attributes and source config.

The transcript is linked to its parent gene by normalized IDs. If the parent
gene was not present in the GFF, reconciliation creates a synthetic gene later.

### Exon Rows

A GFF row with `type == exon` is attached to the transcript identified by the
normalized `Parent`.

If the transcript already exists, the exon is appended to that transcript.

If the transcript does not exist, the loader creates a synthetic relationship:

1. It treats the normalized exon parent as a gene ID.
2. If that gene does not exist, it creates a synthetic `GeneRecord` at the exon
   coordinates.
3. It creates or reuses a dummy transcript named `<parent>_dTx`.
4. It attaches the exon to that dummy transcript.

This preserves exon-only or gene-exon structures well enough for insertion into
the Ensembl core schema, which requires transcript and exon_transcript rows.

### CDS Rows

A GFF row with `type == CDS` becomes a `CdsSegment`. CDS rows are grouped by
normalized `Parent`.

CDS rows are used for:

1. Computing exon `phase` and `end_phase`.
2. Creating `translation` rows.
3. Updating `transcript.canonical_translation_id`.

CDS rows do not create genes or transcripts by themselves. If a CDS parent does
not match a parsed transcript, that CDS group is ignored during phase and
translation insertion.

### Ignored Features

The parser ignores feature types that are not:

```text
gene
exon
CDS
<configured transcript feature type>
```

For example, UTR, start_codon, stop_codon, region, match, and repeat features
are ignored unless support is added to the parser or represented through the
existing model.

## Generic Biotype Handling

Biotype resolution happens in two steps:

1. `resolve_biotype()` chooses an initial biotype while parsing each feature.
2. `apply_biotype_overrides()` performs final gene and transcript rewrites.

### Initial Biotype Resolution

The initial rule order is:

1. Read `gbkey`, `pseudo`, `transcript_biotype`, and `gene_biotype` from the
   configured attribute names.
2. If `gbkey` ends with the configured segment suffix, create a segment biotype
   from the first character of `gbkey`. With the RefSeq defaults, `V_segment`
   becomes `IG_V_gene`.
3. If `gbkey` contains the configured transcribed-pseudogene token, return
   `transcribed_pseudogene`.
4. If `pseudo=true`, return `pseudogene`.
5. If `transcript_biotype` exists, return that value.
6. If the feature type is in `biotype_transcript_feature_types`, return an
   explicit mapped biotype if configured; otherwise return the feature type.
7. If the feature is an exon and `gbkey` is in `exon_gbkey_biotype_map`, return
   that mapped biotype.
8. If `gene_biotype` exists, return that value.
9. Fall back to `default_biotype`.

This order matters. For example, `pseudo=true` is handled before
`transcript_biotype`, and configured transcript feature maps are handled before
`gene_biotype`.

### Final Biotype Overrides

After parsing, reconciliation, and phase calculation, the loader applies final
overrides:

1. If a gene biotype is in `gene_biotype_overrides`, the gene biotype is
   rewritten.
2. If a transcript's parent gene biotype is in
   `transcript_biotype_overrides`, the transcript biotype is rewritten to the
   configured value.

This is how a source can parse one raw value and then normalize it to the final
Ensembl biotype expected by downstream core consumers.

## Reconciliation And Edge Cases

`reconcile_annotation()` makes sure every parsed object can be inserted into an
Ensembl core schema.

Transcript without exons:
: A synthetic exon is added covering the full transcript coordinates.

Transcript whose parent gene is missing:
: A synthetic gene is created using the transcript coordinates and biotype.

Gene without transcripts:
: A dummy transcript named `<gene_id>_dTx` is created. It contains one synthetic
exon covering the full gene coordinates.

Exon whose parent transcript is missing:
: The parser creates a dummy transcript named `<parent>_dTx`. If needed, it also
creates a synthetic gene named `<parent>`.

CDS whose parent transcript is missing:
: The CDS group remains in memory but is skipped for phase and translation
creation because there is no transcript to attach it to.

Missing seq_regions:
: `load-features` fails before committing if any parsed gene or transcript
seqid is not present in the target core database for the selected coord_system.

Duplicate IDs:
: The in-memory annotation stores genes and transcripts in dictionaries keyed by
normalized ID. Later rows with the same normalized ID replace earlier gene or
transcript records.

Multiple Parent values:
: The current parser does not split comma-separated `Parent` values. A value
such as `Parent=tx1,tx2` is treated as one literal parent ID. Pre-normalize such
GFFs if multi-parent features need to be loaded separately.

## Exon Phase And Translation Handling

`compute_exon_phases()` calculates exon phase information from CDS segments.

For each transcript with CDS:

1. Exons are sorted in transcript order. Positive-strand transcripts sort by
   start ascending. Negative-strand transcripts sort by start descending.
2. CDS segments are sorted in coding order. Positive-strand CDS segments sort
   by start ascending. Negative-strand CDS segments sort by end descending.
3. The first valid CDS phase is read from the GFF phase column.
4. All transcript exons are initialized to `phase=-1` and `end_phase=-1`.
5. Each exon is checked for overlap with the ordered CDS segments.
6. Non-coding exons keep `-1` internally.
7. The first coding exon uses the first valid CDS phase if one was present.
8. Later coding exon phases are calculated from the number of coding bases seen
   so far.
9. `end_phase` is calculated as `coding_bases % 3`.

When exons are inserted into the database, `None` and `-1` phases are currently
stored as `0`. Coding phases calculated as `0`, `1`, or `2` are stored as-is.

`insert_translations()` creates translations from CDS groups:

1. For positive strand transcripts, translation start is the minimum CDS start
   and translation end is the maximum CDS end.
2. For negative strand transcripts, translation start is the maximum CDS end
   and translation end is the minimum CDS start.
3. The loader finds the exon containing the translation start and the exon
   containing the translation end.
4. It calculates `seq_start` and `seq_end` offsets inside those exons.
5. It inserts a `translation` row with stable ID `<transcript_id>_prot`.
6. It updates `transcript.canonical_translation_id`.

If the translation start or end cannot be placed inside an exon, no translation
row is inserted for that transcript.

## Database Insert Details

The database insert flow is shared by `load-features` and `create-core` after
input preparation.

`analysis`
: The loader uses `source_config.analysis_logic_name`. Existing analysis rows
are reused in `load-features`; `create-core` inserts one during bootstrap.

`gene`
: Inserts one row per `GeneRecord`. `stable_id` is the normalized gene ID.
`source` is `source_config.source_label`. `canonical_transcript_id` is the first
transcript found for that gene in parsed insertion order. `display_xref_id` is
currently `NULL`.

`transcript`
: Inserts one row per `TranscriptRecord`. `stable_id` comes from configured
attributes, usually `Name`, then falls back to normalized transcript ID.
`source` is `source_config.source_label`. `display_xref_id` is currently
`NULL`.

`exon`
: Inserts exon rows per transcript. Exons are not deduplicated across
transcripts. Positive-strand exons are ranked by start ascending; negative
strand exons are ranked by start descending. `stable_id` is currently `NULL`.

`exon_transcript`
: Links each inserted exon to its transcript with the calculated rank.

`translation`
: Inserted only when a matching transcript and CDS group allow start and end
placement inside exons.

Numeric ID allocation differs by mode:

`create-core`
: Gene and transcript numeric IDs start at `1` and are allocated
deterministically from sorted stable IDs.

`load-features`
: Gene, transcript, and exon numeric IDs start at current table max plus one in
the existing core database.

All database writes run with `autocommit=False`. Any exception rolls back the
transaction and is logged.

## Post-Load Quality Check

Every load path runs a post-load quality check after feature rows have been
inserted and before the transaction is committed. This applies to:

```text
gff-loader load-features
gff-loader create-core
gff-loader refseq run --load-core
```

The check compares the parsed GFF model against the rows inserted into the core
database during that load. It uses the exact numeric IDs allocated by the loader,
so existing rows in the target core do not affect the result.

Checked counts:

```text
gene
transcript
exon
exon_transcript
translation
```

Checked integrity:

```text
gene stable_id and biotype
transcript stable_id and biotype
exon stable_id when the selected source config captures one
exon_transcript links
translation stable_id
gene biotype distribution
transcript biotype distribution
```

The report is printed at the end of the load and also sent through the logger.
Use `--log-file` if the report should be written to a file:

```bash
gff-loader --log-file import.log load-features annotations.gff3 ...
```

Successful report shape:

```text
GFF core quality check: PASS
source: ensembl
expected vs core rows:
  genes: expected=48958 observed=48958 missing=0 unexpected=0 mismatched=0
  transcripts: expected=76019 observed=76019 missing=0 unexpected=0 mismatched=0
  exons: expected=365661 observed=365661 missing=0 unexpected=0 mismatched=0
  exon_transcript_links: expected=365661 observed=365661 missing=0 unexpected=0 mismatched=0
  translations: expected=39824 observed=39824 missing=0 unexpected=0 mismatched=0
```

If anything is missing, unexpected, or mismatched, the report is printed and
logged, then the loader raises an error and rolls back the transaction.

# Generic GFF/GTF Loader And RefSeq Import

This directory contains a modular GFF/GTF loader. The entry point is a single
CLI, `gff-loader`, which can load generic GFF3 features and supported GTF
annotations into Ensembl core databases and can also run RefSeq-specific
discovery, download, conversion, and loading workflows.

The important design point is that the database loader is generic. It does not
need to know whether the annotation came from RefSeq, GENCODE, MAKER, an anno
GTF pipeline, or another producer once the source-specific rules have been
captured in a `GffSourceConfig`.

## Module Layout

Generic modules:

```text
gff_cli.py
gff_core_loader.py
gff_models.py
gff_repeat_loader.py
gff_source_config.py
```

`gff_cli.py`
: The unified CLI. It exposes generic loading commands directly and RefSeq
commands under the `refseq` subcommand group.

`gff_core_loader.py`
: The parser and Ensembl core database loader. It contains the generic feature
flow: parse GFF3/GTF, reconcile missing relationships, compute exon phases,
apply biotype overrides, resolve seq_regions, allocate IDs, and insert rows.

`gff_models.py`
: Typed in-memory records used by the generic loader: `GeneRecord`,
`TranscriptRecord`, `ExonRecord`, `CdsSegment`, and `ParsedAnnotation`.

`gff_repeat_loader.py`
: Loader for anno pipeline `single_line_feature` GTF outputs. It mirrors the
old Perl repeat/simple-feature path by loading repeat analyses into
`repeat_consensus` and `repeat_feature`, and `cpg`/`eponine` into
`simple_feature`.

`gff_source_config.py`
: Source configuration. This is where feature-type sets, ID normalization,
biotype rules, source labels, and analysis metadata are defined. The current
registered configs are `anno_gtf`, `ncrna_gtf`, `generic`, `ensembl`, and
`refseq`.

RefSeq-specific modules:

```text
refseq_ncbi.py
refseq_conversion.py
refseq_models.py
refseq_constants.py
```

`refseq_ncbi.py`
: NCBI RefSeq discovery and download helpers. It reads NCBI
`assembly_summary.txt` files and downloads genomic GFF3, genomic FASTA, and
assembly report files.

`refseq_conversion.py`
: RefSeq conversion helpers. It uses the NCBI assembly report to convert RefSeq
accessions in FASTA headers and GFF3 seqids to Ensembl-style seq_region names.

`refseq_models.py`
: RefSeq/NCBI assembly metadata and downloaded-path records.

`refseq_constants.py`
: NCBI group names and metadata field constants.

## Running The CLI

Install the package in editable mode from the repository root:

```bash
pip install -e .
gff-loader --help
```

You can also run the CLI without installing the console script:

```bash
python src/python/ensembl/genes/ensembl_loading/gff_cli.py --help
```

Set logging verbosity with `--log-level`:

```bash
gff-loader --log-level DEBUG load-features annotations.gff3 ...
```

Write the same log output to a file with `--log-file`:

```bash
gff-loader --log-file gff-loader.log load-features annotations.gff3 ...
```

## General Use

There are two general loading modes.

`load-features`
: Load a GFF3 or supported GTF into an existing Ensembl core database. The
database, schema, coord_system rows, seq_region rows, and DNA must already
exist. This mode only inserts annotation feature rows.

`create-core`
: Create or reuse a core database, load schema SQL, bootstrap core metadata,
load FASTA sequence into `seq_region` and `dna`, then load GFF3/GTF features.
This mode is useful when the FASTA and annotation file already use matching
seq_region names.

### Loading Features Into An Existing Core

Use `load-features` when you already have an Ensembl core database with the
correct coordinate system and seq_region rows.

```bash
gff-loader load-features annotations.gff3 \
  --db-name mus_musculus_core_000001635_27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password>
```

All MySQL connection details are required for commands that connect to a core:
`--db-host`, `--db-port`, `--db-user`, and `--db-password`.
The documented CLI style uses hyphenated option names, but underscore aliases
such as `--db_port` are also accepted for DB connection options.

The non-DB defaults are:

```text
--source generic
--coord-system-name primary_assembly
```

Use `--coord-system-id` if you know the exact coordinate system ID:

```bash
gff-loader load-features annotations.gff3 \
  --db-name mus_musculus_core_000001635_27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-id 1 \
  --source generic
```

For Ensembl-style GFF3 exports, use the dedicated Ensembl source config:

```bash
gff-loader load-features Homo_sapiens.GRCh38.115.chromosome.20.gff3.gz \
  --db-name homo_sapiens_core_115_38 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name chromosome \
  --source ensembl
```

Use `--coord-system-name primary_assembly` instead if the target core stores the
GFF seqids under that coordinate system. The right value is the coordinate
system containing seq_region names such as `1`, `20`, `X`, `MT`, or scaffold
names from the GFF.

Use `--coord-system-name` and `--coord-system-version` when the target core has
multiple coordinate systems:

```bash
gff-loader load-features annotations.gff3 \
  --db-name example_core \
  --db-host mysql-host \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name scaffold \
  --coord-system-version GCA_000000000.1
```

For anno pipeline GTF files, use the dedicated `anno_gtf` and `ncrna_gtf`
source configs.

Important: `annotation_output/annotation.gtf` is only the main gene set. It
does not contain the Rfam and tRNAscan ncRNA genes. To reproduce the original
anno pipeline core load, load the main annotation GTF first with
`--source anno_gtf`, then load the separate ncRNA GTF outputs with
`--source ncrna_gtf`.

The simplest way to reproduce the full anno pipeline load is the combined
wrapper command. Point it at the anno run output directory that contains
subdirectories such as `annotation_output`, `rfam_output`, `dust_output`, and
`cpg_output`:

```bash
gff-loader load-anno-output /path/to/GCA_000000000.1 \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly
```

The wrapper requires a non-empty
`annotation_output/annotation.gtf`. Optional ncRNA, repeat, and simple-feature
outputs are loaded when their `annotation.gtf` files exist and are non-empty;
missing or empty optional files are skipped.

The usual anno pipeline load sequence is:

```bash
gff-loader load-features /path/to/annotation_output/annotation.gtf \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source anno_gtf

gff-loader load-features /path/to/rfam_output/annotation.gtf \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source ncrna_gtf

gff-loader load-features /path/to/trnascan_output/annotation.gtf \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source ncrna_gtf
```

The Rfam/tRNAscan loads insert genes with source `ensembl` and analysis logic
name `ncrna`, matching the original Perl loader. If either ncRNA output file is
empty or absent, skip that file.

Anno repeat/simple outputs are loaded separately with
`load-single-line-features`, matching the old Perl
`-load_type single_line_feature` mode:

```bash
gff-loader load-single-line-features /path/to/dust_output/annotation.gtf \
  --analysis-name dust \
  --db-name lepidoptera_example_core \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly
```

Use the matching analysis name for each anno output:

```text
dust_output/annotation.gtf         -> dust
red_output/annotation.gtf          -> repeatdetector
trf_output/annotation.gtf          -> trf
repeatmasker_output/annotation.gtf -> repeatmask_repbase_human
cpg_output/annotation.gtf          -> cpg
eponine_output/annotation.gtf      -> eponine
```

The repeat analyses load into `repeat_consensus` and `repeat_feature`.
`cpg` and `eponine` load into `simple_feature`.

Set `--coord-system-name` to the coordinate system containing the GTF seqids.
For example, if the GTF first column contains names such as `17`, the target
core must already have matching `seq_region.name` values under the chosen
coordinate system.

`load-features` does not create a database, does not load schema SQL, does not
insert species or assembly metadata, and does not load FASTA sequence. It
parses the GFF3/GTF and inserts into these core feature tables:

```text
analysis
gene
transcript
exon
exon_transcript
translation
```

The `analysis` row is reused if `analysis.logic_name` already exists for the
selected source config. Otherwise it is created from the config.

### Creating A Core From FASTA And GFF3/GTF

Use `create-core` when you have a FASTA and GFF3/GTF whose sequence names
already match each other.

```bash
gff-loader create-core annotations.gff3 genome.fna \
  --species-name "Mus musculus" \
  --assembly-accession GCF_000001635.27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --source generic
```

By default, `create-core` loads the bundled schema:

```text
src/python/ensembl/genes/ensembl_loading/config/core_schema.sql
```

Use `--schema-sql-path /path/to/schema.sql` to override it. Use
`--schema-sql-path ""` to skip schema loading explicitly.

The database name is derived from species and assembly accession. Generic
sources keep the compact historical format:

```text
Mus musculus + GCF_000001635.27 -> mus_musculus_core_000001635_27
```

When `--source refseq` is used, the derived name follows the RefSeq production
style:

```text
Scientific name + GCF_037462849.1 -> scientific_name_gcf037462849v1_core_114_1
```

`create-core` performs these operations in one transaction:

1. Connect to the MySQL server without selecting a database.
2. `CREATE DATABASE IF NOT EXISTS <derived_core_db_name>`.
3. `USE <derived_core_db_name>`.
4. Load schema SQL from the bundled `config/core_schema.sql`, unless
   `--schema-sql-path` overrides or disables it.
5. Insert one `coord_system` row:
   `primary_assembly`, empty version, rank `1`,
   `default_version,sequence_level`.
6. Insert core `meta` rows for species name, assembly accession, assembly name,
   and genebuild/transcriptbuild/exonbuild levels.
7. Insert an `analysis` row from the selected source config.
8. Load FASTA sequences into `seq_region`, `dna`, and `seq_region_attrib`.
9. Parse and load GFF3 features.
10. Commit, or roll back the whole transaction if any step fails.

The optional positional `assembly_report` argument is retained for compatibility
with RefSeq-style calling code, but the core loader itself does not read the
assembly report. The FASTA and GFF3 must already contain the sequence names that
will be used as `seq_region.name`.

### Input Requirements

The loader expects GFF3-like rows with nine tab-separated columns. Comment rows
starting with `#` are skipped. Plain text and `.gz` GFF3 inputs are supported.

The parser expects standard `ID` and `Parent` attributes. Attributes are parsed
by splitting on semicolon and then on the first equals sign. This is simple and
matches the original loading assumptions; it does not URL-decode values and does
not split comma-separated multi-parent relationships.

Important input assumptions:

1. GFF seqids must match `seq_region.name` in the target core for
   `load-features`.
2. FASTA headers must match GFF seqids for `create-core`.
3. A feature with `strand` equal to `+` is loaded as `1`; any other strand value
   is currently loaded as `-1`.
4. Multi-parent exons or CDS rows should be pre-normalized before loading.
5. Features outside the configured gene/transcript/exon/CDS model are ignored.

### Why Coordinate System Options Are Needed

Ensembl stores a feature location as `seq_region_id`, start, end, and strand.
The loader therefore has to convert each GFF seqid into a `seq_region_id`.

For `load-features`, the `seq_region` rows already exist. The loader resolves
them using either:

```text
--coord-system-id
```

or:

```text
--coord-system-name
--coord-system-version
```

`--coord-system-id`
: Uses an exact `coord_system_id`. This is the safest option when the database
has several coordinate systems with similar names or repeated seq_region names.

`--coord-system-name`
: Looks up the first matching coordinate system by name, ordered by rank. The
default is `primary_assembly`.

`--coord-system-version`
: Adds a version filter to the name lookup.

If any parsed gene or transcript refers to a seqid that cannot be found in the
resolved coordinate system, loading fails before inserts are committed. The
error reports the first missing names and the total count if there are more
than ten.

`create-core` does not expose these options because it creates a new
`primary_assembly` coordinate system and inserts seq_regions from the FASTA.

## Source Configs

Feature interpretation is controlled by `GffSourceConfig` in
`gff_source_config.py`.

The config defines:

`name`
: CLI selector used by `--source`.

`source_label`
: Value inserted into `gene.source` and `transcript.source`.

`analysis_logic_name`
: Value used to find or create the `analysis` row.

`analysis_program`
: Value inserted into `analysis.program` when a new analysis row is created.

`parsed_transcript_feature_types`
: GFF feature types that become `TranscriptRecord` objects.

`parsed_gene_feature_types`
: GFF feature types that become `GeneRecord` objects. This is needed because
some sources use Sequence Ontology terms such as `ncRNA_gene` or `pseudogene`
instead of only the literal type `gene`.

`biotype_transcript_feature_types`
: GFF feature types that can directly contribute a transcript biotype.

`transcript_feature_biotype_map`
: Explicit feature-type to biotype mappings. For example, `mRNA` can become
`protein_coding` instead of the literal feature type.

`exon_gbkey_biotype_map`
: Optional `gbkey` to biotype mappings used when resolving exon-derived
synthetic features.

`gene_biotype_overrides`
: Final gene biotype rewrites applied after parsing.

`transcript_biotype_overrides`
: Final transcript biotype rewrites based on the parent gene biotype.

`id_prefixes_to_strip`
: Prefixes removed from GFF IDs and Parent values before internal IDs are used.

`attribute_format`
: Attribute syntax parser. `gff3` reads `key=value` pairs; `gtf` reads
`key "value";` pairs.

`gene_id_attribute`, `transcript_id_attribute`, `parent_gene_attribute`,
`exon_parent_attribute`, `cds_parent_attribute`
: Attribute names used to identify feature IDs and parent-child relationships.
GFF3 configs usually use `ID` and `Parent`; GTF configs usually use `gene_id`
and `transcript_id`.

`gene_name_attributes`
: Attribute priority used to choose `GeneRecord.name`.

`transcript_stable_id_attributes`
: Attribute priority used to choose `TranscriptRecord.stable_id`.

`exon_stable_id_attributes`
: Attribute priority used to choose `ExonRecord.stable_id`. If no configured
attribute is found, exon stable IDs are inserted as `NULL`.

`translation_stable_id_attributes`
: Attribute priority used to choose `translation.stable_id` from CDS rows. If
no configured attribute is found, the fallback stable ID is
`<transcript_id>_prot`.

`translation_coords_attribute`
: Optional transcript attribute used to synthesize CDS segments when the input
has transcript/exon rows but no CDS rows. The anno GTF config uses
`translation_coords`.

`transcript_rows_define_genes`
: Whether transcript rows should create and expand parent gene records. This is
used for anno GTF files that do not contain explicit gene rows.

`gene_xref_prefix`
: Prefix used to detect a GeneID inside `Dbxref`. The value is captured in the
in-memory gene record, although the current DB insertion path does not create
xref or display_xref rows.

`default_biotype`
: Fallback biotype when no source-specific rule applies.

`toplevel_attrib_type_id`
: Attribute type ID inserted into `seq_region_attrib` when `create-core` loads
FASTA seq_regions.

### Generic Source Config

The default generic source is intentionally conservative:

```text
name: generic
source_label: gff
analysis_logic_name: gff_import
analysis_program: GFF3
id_prefixes_to_strip: ()
default_biotype: protein_coding
```

Generic transcript-like feature types:

```text
mRNA
transcript
lnc_RNA
snRNA
rRNA
snoRNA
ncRNA
antisense_RNA
scRNA
piRNA
siRNA
tRNA
```

Generic biotype rules:

1. `pseudo=true` becomes `pseudogene`.
2. `transcript_biotype=<value>` wins when present.
3. `mRNA` becomes `protein_coding`.
4. Other configured transcript feature types keep their feature type as the
   biotype unless mapped.
5. `gene_biotype=<value>` is used after transcript/exon-specific checks.
6. If nothing matches, the biotype falls back to `protein_coding`.

The generic config does not strip ID prefixes. A GFF row with `ID=gene1` remains
`gene1`; a row with `ID=gene-gene1` also remains `gene-gene1`.

### Ensembl Source Config

Use `--source ensembl` for GFF3 files exported from Ensembl FTP or equivalent
Ensembl-style dumps.

The Ensembl config is different from the generic config in three important
ways:

1. It treats `gene`, `ncRNA_gene`, and `pseudogene` rows as gene rows.
2. It strips Ensembl GFF3 ID prefixes such as `gene:`, `transcript:`, `exon:`,
   `CDS:`, `protein:`, and `chromosome:`.
3. It reads the Ensembl `biotype` attribute for both genes and transcripts.
4. It reads exon stable IDs from `exon_id`, then `Name`, then `ID`.
5. It reads translation stable IDs from CDS `protein_id`, then `Name`, then
   `ID`.

Core metadata:

```text
name: ensembl
source_label: ensembl
analysis_logic_name: ensembl_gff_import
analysis_program: Ensembl_GFF3
default_biotype: protein_coding
```

Recognized Ensembl gene feature types:

```text
gene
ncRNA_gene
pseudogene
```

Recognized Ensembl transcript-like feature types:

```text
mRNA
transcript
pseudogenic_transcript
lnc_RNA
snRNA
rRNA
snoRNA
ncRNA
antisense_RNA
scRNA
piRNA
siRNA
miRNA
tRNA
vault_RNA
Y_RNA
RNase_MRP_RNA
RNase_P_RNA
telomerase_RNA
D_gene_segment
V_gene_segment
J_gene_segment
C_gene_segment
C_region
```

For a typical Ensembl GFF3 gene block:

```text
20 havana gene 100 900 . + . ID=gene:ENSG000001;Name=ABC1;biotype=protein_coding;gene_id=ENSG000001
20 havana mRNA 100 900 . + . ID=transcript:ENST000001;Parent=gene:ENSG000001;Name=ABC1-201;biotype=protein_coding;transcript_id=ENST000001
20 havana exon 100 200 . + . Parent=transcript:ENST000001;Name=ENSE000001;exon_id=ENSE000001;rank=1
20 havana CDS 120 200 . + 0 Parent=transcript:ENST000001;protein_id=ENSP000001
```

the loader stores:

```text
gene stable_id:       ENSG000001
transcript stable_id: ENST000001
exon stable_id:       ENSE000001
translation stable_id: ENSP000001
gene biotype:         protein_coding
transcript biotype:   protein_coding
```

Rows such as `chromosome`, `biological_region`, UTRs, and other non-core
annotation features are ignored by the current feature loader.

### Anno GTF Source Config

Use `--source anno_gtf` for GTF files produced by the anno annotation pipeline.
These files contain `transcript` and `exon` rows, with relationships expressed
through GTF attributes:

```text
17 ensembl transcript 169704 183277 . + . gene_id "annotation_1"; transcript_id "annotation_1"; biotype "protein_coding"; translation_coords "169704:169787:11:183093:183277:109";
17 ensembl exon       169704 169787 . + . gene_id "annotation_1"; transcript_id "annotation_1"; exon_number "1";
```

Core metadata:

```text
name: anno_gtf
source_label: ensembl
analysis_logic_name: ensembl
analysis_program: Anno_GTF
attribute_format: gtf
default_biotype: not_set
```

The loader handles anno GTF files as follows:

1. `gene_id` becomes the gene stable ID.
2. `transcript_id` becomes the transcript stable ID.
3. Parent gene records are synthesized from transcript rows and expanded to the
   min/max span of all transcripts with the same `gene_id`.
4. Transcript `biotype` is used directly when present.
5. If `biotype` is absent, transcripts with `translation_coords` become
   `protein_coding`; transcripts without `translation_coords` become `not_set`.
6. Transcript `translation_coords` is converted into CDS segments, exon phases,
   and translation rows. If no `translation_coords` value is present, no
   translation is created.

This config only covers the main `annotation_output/annotation.gtf` file.
Rfam and tRNAscan ncRNA genes are produced as separate pipeline outputs and
must be loaded separately with `--source ncrna_gtf`.

### ncRNA GTF Source Config

Use `--source ncrna_gtf` for the anno pipeline's Rfam and tRNAscan GTF files,
for example `rfam_output/annotation.gtf` and
`trnascan_output/annotation.gtf`.

Core metadata:

```text
name: ncrna_gtf
source_label: ensembl
analysis_logic_name: ncrna
analysis_program: Anno_ncRNA_GTF
attribute_format: gtf
default_biotype: misc_RNA
```

This config uses the same GTF relationship rules as `anno_gtf`: `gene_id`
defines the parent gene, `transcript_id` defines the transcript, and gene rows
are synthesized from transcript rows when needed. The transcript `biotype`
attribute is used directly, so Rfam values such as `misc_RNA`, `rRNA`,
`ribozyme`, `snRNA`, and `snoRNA` are preserved. If a transcript has no
`biotype`, it defaults to `misc_RNA`.

### Single-Line Repeat/Simple Feature Loader

Use `load-single-line-features` for the anno pipeline outputs that the old Perl
loader handled with `-load_type single_line_feature`.

Supported repeat analyses:

```text
dust
repeatdetector
repeatmask_repbase_human
trf
```

For repeats, the loader reads GTF attributes `repeat_name`, `repeat_class`,
`repeat_type`, `repeat_consensus`, and `score` when present. Missing
`repeat_name` and `repeat_class` default to the analysis name, missing
`repeat_consensus` defaults to `N`, and repeat coordinates are loaded with
`repeat_start = 1` and `repeat_end = feature length`, matching the Perl loader.

Supported simple-feature analyses:

```text
cpg
eponine
```

For simple features, rows are inserted into `simple_feature` with an empty
display label and score `0`, matching the old Perl behavior.

### Adding Another GFF Source

To add another source, create a new `GffSourceConfig` in
`gff_source_config.py` and add it to `SOURCE_CONFIGS`.

Example shape:

```python
CUSTOM_CONFIG = GffSourceConfig(
    name="custom",
    source_label="custom",
    analysis_logic_name="custom_gff_import",
    analysis_program="CustomGFF",
    parsed_gene_feature_types=frozenset({"gene"}),
    parsed_transcript_feature_types=frozenset({"mRNA", "transcript"}),
    biotype_transcript_feature_types=frozenset({"mRNA", "transcript"}),
    transcript_feature_biotype_map={"mRNA": "protein_coding"},
    id_prefixes_to_strip=("gene:", "transcript:"),
)

SOURCE_CONFIGS[CUSTOM_CONFIG.name] = CUSTOM_CONFIG
```

Use it from the CLI:

```bash
gff-loader load-features custom.gff3 \
  --db-name custom_core \
  --db-host mysql-host \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --source custom
```

You should only need parser changes if the source does not express features as
gene, transcript-like feature, exon, and CDS rows with configured relationship
attributes, or as anno-style GTF transcript/exon rows with `translation_coords`.

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

## RefSeq Import

RefSeq support is built as a source-specific layer around the generic GFF
loader. It adds:

1. Discovery of available NCBI RefSeq assemblies.
2. Download of genomic GFF3, genomic FASTA, and assembly report files.
3. Conversion of RefSeq sequence accessions to Ensembl-style seq_region names.
4. RefSeq-specific feature and biotype rules through `REFSEQ_CONFIG`.
5. A `refseq run` command for download plus conversion, with optional core DB
   loading.

The generic database insertion code is the same once the converted RefSeq GFF3
has been parsed.

### RefSeq Source Config

`REFSEQ_CONFIG` is registered under `--source refseq`.

Core metadata:

```text
name: refseq
source_label: refseq
analysis_logic_name: refseq_import
analysis_program: NCBI_RefSeq
id_prefixes_to_strip: gene-, rna-, cds-, exon-
default_biotype: protein_coding
```

RefSeq transcript-like feature types:

```text
mRNA
transcript
lnc_RNA
snRNA
rRNA
snoRNA
ncRNA
antisense_RNA
scRNA
telomerase_RNA
RNase_P_RNA
SRP_RNA
RNase_MRP_RNA
piRNA
siRNA
tRNA
pseudogenic_tRNA
D_gene_segment
V_gene_segment
J_gene_segment
C_gene_segment
C_region
```

RefSeq feature types used for direct transcript-biotype resolution:

```text
mRNA
transcript
lnc_RNA
snRNA
rRNA
snoRNA
ncRNA
antisense_RNA
scRNA
telomerase_RNA
RNase_P_RNA
SRP_RNA
RNase_MRP_RNA
piRNA
siRNA
tRNA
pseudogenic_tRNA
C_region
precursor_RNA
```

Explicit RefSeq feature-type mappings:

```text
mRNA    -> protein_coding
C_region -> IG_C_gene
```

RefSeq exon `gbkey` mappings:

```text
ncRNA         -> ncRNA
precursor_RNA -> precursor_RNA
C_region      -> IG_C_gene
```

RefSeq final gene biotype overrides:

```text
V_segment -> IG_V_gene
D_segment -> IG_D_gene
J_segment -> IG_J_gene
C_segment -> IG_C_gene
C_region  -> IG_C_gene
```

RefSeq final transcript biotype overrides based on parent gene biotype:

```text
lncRNA                 -> lncRNA
antisense_RNA          -> antisense_RNA
pseudogene             -> pseudogene
transcribed_pseudogene -> transcribed_pseudogene
rRNA                   -> rRNA
snRNA                  -> snRNA
snoRNA                 -> snoRNA
tRNA                   -> tRNA
miRNA                  -> miRNA
ncRNA                  -> ncRNA
misc_RNA               -> misc_RNA
telomerase_RNA         -> telomerase_RNA
RNase_P_RNA            -> RNase_P_RNA
SRP_RNA                -> SRP_RNA
RNase_MRP_RNA          -> RNase_MRP_RNA
IG_V_gene              -> IG_V_gene
IG_D_gene              -> IG_D_gene
IG_J_gene              -> IG_J_gene
IG_C_gene              -> IG_C_gene
```

### RefSeq Feature Handling

RefSeq IDs often contain source prefixes:

```text
ID=gene-GeneA
ID=rna-NM_001
ID=cds-XP_001
ID=exon-123
Parent=gene-GeneA
Parent=rna-NM_001
```

With `--source refseq`, the loader strips `gene-`, `rna-`, `cds-`, and `exon-`
from IDs and Parent values before linking records. This lets a transcript
parent of `gene-GeneA` link to the parsed gene ID `GeneA`, and a CDS parent of
`rna-NM_001` link to transcript `NM_001`.

RefSeq genes:

1. `type=gene` becomes a `GeneRecord`.
2. `type=pseudogene` also becomes a `GeneRecord`.
3. Stable ID is the normalized `ID`.
4. Gene name is taken from `Name`, then `gene`, then normalized ID.
5. `Dbxref=GeneID:<value>` is captured in memory.
6. Biotype is resolved from `gbkey`, `pseudo`, `transcript_biotype`,
   `gene_biotype`, and RefSeq override maps.

RefSeq transcripts:

1. Any configured transcript-like type becomes a `TranscriptRecord`.
2. Stable ID is usually `Name`, which often preserves accessions such as
   `NM_...`, `NR_...`, or `XM_...`.
3. Parent is normalized and used as the gene ID.
4. Biotype is resolved using the RefSeq rule order.

RefSeq exons:

1. `type=exon` rows attach to the normalized transcript parent.
2. If the exon parent does not match a transcript, dummy transcript handling is
   used as described in the generic reconciliation section.
3. Exon `gbkey` can drive synthetic biotypes for exon-only structures.

RefSeq CDS:

1. `type=CDS` rows are grouped by normalized transcript parent.
2. CDS phase values drive exon phase calculation.
3. CDS coordinates drive `translation` insertion.
4. CDS rows whose parent does not match a parsed transcript do not create a
   translation.

RefSeq immunoglobulin and segment features:

1. `gbkey` values ending in `_segment` are converted to IG segment biotypes
   using the first character of `gbkey`.
2. `V_segment` becomes `IG_V_gene`.
3. `D_segment` becomes `IG_D_gene`.
4. `J_segment` becomes `IG_J_gene`.
5. `C_segment` and `C_region` become `IG_C_gene`.
6. Transcript biotypes are then normalized from the parent gene biotype by
   `transcript_biotype_overrides`.

RefSeq pseudogene handling:

1. If `gbkey` contains `Transcribed_Pseudogene`, the biotype becomes
   `transcribed_pseudogene`.
2. If `pseudo=true`, the biotype becomes `pseudogene`.
3. These checks happen before `transcript_biotype` and feature-type mappings.

### RefSeq Sequence Name Conversion

RefSeq GFF3 and FASTA files usually use RefSeq accessions such as `NC_000001.11`
or `NT_...` as sequence IDs. Ensembl-style cores normally use chromosome or
scaffold names. The conversion module uses the NCBI assembly report.

`load_refseq_name_map()` reads each non-comment assembly report row:

```text
column 0: sequence name
column 2: assigned molecule
column 6: RefSeq accession
```

The mapping rule is:

1. If assigned molecule is not `na`, use assigned molecule.
2. Otherwise use sequence name.
3. Key the mapping by RefSeq accession.

Example:

```text
NC_000001.11 -> 1
NT_187361.1  -> HSCHR1_CTG1_UNLOCALIZED
```

`convert_gff_to_ensembl()` applies the mapping to:

1. GFF feature row column 1.
2. `##sequence-region` directives.

If a feature seqid is not in the assembly report, the original seqid is kept
and a warning count is logged. If `--chrom-filter` is supplied, rows are kept
when either the original RefSeq accession or the converted name matches the
filter. Malformed feature rows with fewer than nine columns are skipped.

`convert_fna_headers()` applies the same mapping to FASTA headers. It follows
the original script's plain-text FASTA handling, so pass a decompressed FASTA
file when calling it directly.

### Listing RefSeq Assemblies

List assemblies across all configured NCBI groups:

```bash
gff-loader refseq list \
  --base-dir refseq_data \
  --max-print 20 \
  --output-tsv refseq_annotation_metadata.tsv
```

Restrict to one group:

```bash
gff-loader refseq list \
  --base-dir refseq_data \
  --group vertebrate_mammalian \
  --max-print 20
```

Use an empty TSV path to skip metadata output:

```bash
gff-loader refseq list --output-tsv ""
```

Implementation flow:

```text
list_available_annotations()
  -> fetch_groups()
  -> fetch_assembly_summary()
  -> parse_assembly_summary()
  -> summarize_annotations()
  -> write_annotation_metadata_tsv()
```

`parse_assembly_summary()` skips comments, blank rows, rows shorter than the
expected NCBI column count, and rows whose `ftp_path` is `na`. The species name
is derived from the first two words of NCBI's organism name.

### Downloading RefSeq Assemblies

Download one assembly by accession:

```bash
gff-loader refseq download \
  --base-dir refseq_data \
  --assembly-acc GCF_000001635.27 \
  --max-workers 2
```

Download the first matching species:

```bash
gff-loader refseq download \
  --base-dir refseq_data \
  --species-name "Mus musculus"
```

Download all latest assemblies from a group:

```bash
gff-loader refseq download \
  --base-dir refseq_data \
  --group vertebrate_mammalian \
  --max-workers 2
```

Implementation flow:

```text
download_annotations()
  -> find_annotation_targets()
  -> download_assembly()
  -> download_file()
```

Target selection:

`--assembly-acc`
: Searches configured RefSeq groups until the exact accession is found.

`--species-name`
: Returns the first assembly whose NCBI organism name starts with the requested
species string, case-insensitively.

`--group`
: Downloads records in that group whose `version_status` is `latest`.

Downloaded files are stored under an NCBI-style accession path:

```text
<base_dir>/GCF/000/001/635/GCF_000001635.27/
```

For one assembly, the downloader expects these files:

```text
<ftp_base>_genomic.gff.gz
<ftp_base>_genomic.fna.gz
<ftp_base>_assembly_report.txt
```

Existing local files are reused. Empty downloads raise an error. Multiple
downloads use `ThreadPoolExecutor` with `--max-workers`, capped to the number of
targets.

### Converting RefSeq Files

Convert RefSeq GFF3 seqids:

```bash
gff-loader refseq convert-gff \
  refseq_data/GCF/000/001/635/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.gff.gz \
  refseq_data/GCF/000/001/635/GCF_000001635.27/GCF_000001635.27_GRCm39_assembly_report.txt \
  --output refseq_data/GCF/000/001/635/GCF_000001635.27/GCF_000001635.27_GRCm39_ensembl.gff3
```

Filter to a converted chromosome name or original RefSeq accession:

```bash
gff-loader refseq convert-gff input.gff.gz assembly_report.txt \
  --output chr10.gff3 \
  --chrom-filter 10 \
  --chrom-filter NC_000076.7
```

Convert RefSeq FASTA headers:

```bash
gff-loader refseq convert-fna \
  refseq_data/GCF/000/001/635/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna.gz \
  refseq_data/GCF/000/001/635/GCF_000001635.27/GCF_000001635.27_GRCm39_assembly_report.txt \
  refseq_data/GCF/000/001/635/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic_ensembl.fna
```

The GFF and FASTA converter accepts plain text or `.gz` input.

### End-To-End RefSeq Run

`refseq run` combines download and conversion:

```bash
gff-loader refseq run \
  --base-dir refseq_data \
  --assembly-acc GCF_000146045.2
```

Default converted outputs are written beside the downloaded files:

```text
<ftp_base>_genomic_ensembl.fna
<ftp_base>_ensembl.gff3
```

You can override converted output paths for a single assembly:

```bash
gff-loader refseq run \
  --base-dir refseq_data \
  --assembly-acc GCF_000001635.27 \
  --converted-fna /path/to/genome_ensembl.fna \
  --converted-gff /path/to/annotation_ensembl.gff3
```

Custom converted output paths are rejected when more than one assembly is being
processed, because a single path cannot represent multiple assemblies.

`refseq run --load-core` downloads, converts, and then calls the generic
`create-core` path with `--source refseq` by default:

```bash
gff-loader refseq run \
  --base-dir refseq_data \
  --assembly-acc GCF_000001635.27 \
  --load-core \
  --species-name "Mus musculus" \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 4532 \
  --db-user <write-user> \
  --db-password <password>
```

`--load-core` requires:

```text
--species-name
--db-host
--db-port
--db-user
--db-password
```

`--assembly-acc` selects the RefSeq assembly to download. When `--load-core` is
used, that same downloaded assembly accession is also used for core metadata and
the derived core DB name. RefSeq core names use this shape:

```text
<scientific_name>_<lowercase_accession_without_underscore_and_dot_as_v>_core_<schema_version>_<assembly_version>
```

For example, `Scientific name` with `GCF_037462849.1` becomes
`scientific_name_gcf037462849v1_core_114_1`. `--assembly-accession` remains
available only as an override for unusual cases where the loaded core metadata
should use a different assembly identifier.

`--load-core` with `--group` is not supported, because group downloads can
produce many assemblies and the current loader creates one core database per
assembly load.

### Loading Converted RefSeq Files Manually

After converting RefSeq files, you can load them explicitly through the generic
core creation command:

```bash
gff-loader create-core \
  GCF_000001635.27_GRCm39_ensembl.gff3 \
  GCF_000001635.27_GRCm39_genomic_ensembl.fna \
  GCF_000001635.27_GRCm39_assembly_report.txt \
  --species-name "Mus musculus" \
  --assembly-accession GCF_000001635.27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user ensadmin \
  --db-password <password> \
  --source refseq
```

You can also load a converted RefSeq GFF into an existing core:

```bash
gff-loader load-features GCF_000001635.27_GRCm39_ensembl.gff3 \
  --db-name mus_musculus_gcf000001635v27_core_114_27 \
  --db-host mysql-ens-genebuild-prod-6 \
  --db-port 3306 \
  --db-user <write-user> \
  --db-password <password> \
  --coord-system-name primary_assembly \
  --source refseq
```

Use this when the target database already contains matching seq_regions and DNA.

## Python Functions

The CLI is a wrapper around reusable functions.
The DB-loading functions require all MySQL connection fields explicitly,
including `db_port`.

Generic loading functions:

```python
from ensembl.genes.ensembl_loading.gff_core_loader import (
    load_gff_features_to_core,
    load_to_ensembl_core,
    prepare_annotation_for_load,
)
from ensembl.genes.ensembl_loading.gff_source_config import get_source_config

source_config = get_source_config("generic")

summary = load_gff_features_to_core(
    gff_path="annotations.gff3",
    db_name="example_core",
    db_host="mysql-host",
    db_user="write-user",
    db_password="password",
    db_port=3306,
    source_config=source_config,
)
```

RefSeq functions:

```python
from ensembl.genes.ensembl_loading.refseq_ncbi import (
    list_available_annotations,
    download_annotations,
)
from ensembl.genes.ensembl_loading.refseq_conversion import (
    convert_fna_headers,
    convert_gff_to_ensembl,
)
```

These functions use standard `logging` rather than printing directly, so callers
can configure verbosity, capture logs in tests, or integrate them into a larger
pipeline.

## Current Boundaries

The current modular loader is generic at the GFF feature-loading layer. RefSeq
download and assembly-report conversion remain RefSeq-specific.

The generic loader is suitable for another source when:

1. The source can be represented as gene, transcript-like, exon, and CDS rows.
2. IDs and Parent values can be normalized with prefix stripping.
3. Biotypes can be resolved from feature type and attributes using
   `GffSourceConfig`.
4. GFF seqids already match target core seq_regions, or a matching FASTA is
   available for `create-core`.

Parser or loader changes may be needed for:

1. Multiple Parent relationships.
2. URL-decoded or escaped GFF attributes.
3. Unstranded features that should not be represented as `-1`.
4. Feature types that need dedicated Ensembl tables outside gene, transcript,
   exon, exon_transcript, and translation.
5. Xref/display label insertion.
6. Stable ID version tables or external DB metadata.

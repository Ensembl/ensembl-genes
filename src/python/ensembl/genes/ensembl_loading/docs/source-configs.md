# Source Configs

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

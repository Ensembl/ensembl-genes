# RefSeq Import

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
<scientific_name>_<lowercase_accession_without_underscore_and_dot_as_v>_rs_core_<schema_version>_1
```

For example, `Scientific name` with `GCF_037462849.1` becomes
`scientific_name_gcf037462849v1_rs_core_114_1`. `--assembly-accession` remains
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

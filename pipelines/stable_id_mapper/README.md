# Stable ID Mapper

Prototype tools for mapping Ensembl stable IDs between gene annotations.

This pipeline lives in `pipelines/stable_id_mapper` and is currently separate
from the packaged `src/python/ensembl/genes` modules.

## Stable ID Mapping

Ensembl stable IDs identify genes, transcripts, exons, and proteins across
releases. When an annotation is updated, IDs for equivalent features should be
carried forward where possible, and new stable IDs should be created for new
features. Ensembl's public stable-ID documentation gives more background on the
ID format and versioning rules:

https://www.ensembl.org/info/genome/stable_ids/index.html

This branch includes two mapping workflows:

1. `main.py` / `mapper.py`: maps IDs from an old/reference FASTA+GFF3 to a
   new/target FASTA+GFF3 by aligning reference gene sequences to the target
   assembly.
2. `lifton_id_mapper.py`: maps a LiftOn GFF3 to a reference/Ensembl GFF3 that
   is already on the same target assembly, using transcript structure and local
   context.

This directory also includes small chicken chromosome 9 example inputs:

```text
ref.fa
ref.gff3
tar.fa
tar.gff3
```

## Requirements

The scripts are written for Python 3.

`lifton_id_mapper.py` uses only the Python standard library.

`main.py` can run with only the Python standard library, but realistic
assembly-to-assembly mapping requires `minimap2`. If `minimap2` is not
available, the code falls back to exact full-gene sequence search, which is
useful as a smoke test but usually too strict for real assembly updates.

Install `minimap2` with one of:

```bash
brew install minimap2
```

```bash
conda install -c bioconda minimap2
```

For repository development tooling, from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Workflow 1: Assembly-to-Assembly Mapper

Use `main.py` when you have:

- old/reference FASTA
- old/reference GFF3
- new/target FASTA
- new/target GFF3

Show the command-line options:

```bash
python3 main.py --help
```

Run the bundled example:

```bash
mkdir -p out/root-smoke

python3 main.py \
  --ref-fasta ref.fa \
  --ref-gff ref.gff3 \
  --target-fasta tar.fa \
  --target-gff tar.gff3 \
  --output-gff out/root-smoke/mapped.gff3 \
  --report out/root-smoke/report.txt \
  --threads 8 \
  --identity-min 0.80
```

Outputs:

```text
out/root-smoke/mapped.gff3
out/root-smoke/report.txt
```

### Assembly Mapper Algorithm

At a high level, `main.py` calls `mapper.map_ids()` and performs these steps:

1. Load a minimal GFF3 hierarchy of `gene`, `mRNA`/`transcript`, `exon`, and
   `CDS` features.
2. Scan reference and target IDs to learn existing stable-ID prefixes, numeric
   widths, and namespaces such as `gene:` and `transcript:`.
3. Extract each reference gene sequence from the reference FASTA.
4. Align reference gene sequences to the target FASTA with `minimap2`, if it is
   installed.
5. Fall back to exact sequence matching if `minimap2` is unavailable.
6. Project transcript, exon, and CDS coordinates through the gene-level
   alignment CIGAR.
7. For mapped genes, keep the reference stable IDs and increment the `version`
   attribute.
8. For target genes that are not replaced by a mapped reference gene, create
   new stable IDs using target-derived prefixes.
9. Write a mapped GFF3 and a text report listing mapped, missing, and new genes.

Relevant files:

```text
main.py        CLI wrapper
mapper.py      mapping flow, prefix inference, version handling, report writing
aligner.py     minimap2 invocation and exact-match fallback
cigar_map.py   CIGAR interval projection
fasta_io.py    in-memory FASTA reader
gff_io.py      minimal GFF3 parser and writer
id_manager.py  stable-ID parsing and new-ID allocation
synteny.py     placeholder synteny scoring
```

### Assembly Mapper Smoke Test

If `minimap2` is not installed, the bundled example is expected to map very few
genes because it uses exact full-gene sequence matching. A fallback-path smoke
test produced:

```text
Genes in reference: 451
Mapped: 1 (0.1%)
Missing: 450 (53.9%)
New: 384
```

This proves the command runs, but it should not be interpreted as a useful
biological mapping result. Install `minimap2` before evaluating mapping quality.

### Core-DB SQL Generator

`main_output_to_stable_id_event_sql.py` is a quick converter from `main.py`
output to SQL that can be reviewed and then run against an Ensembl core DB.
The script does not connect to MySQL itself; it only writes SQL and optional TSV
files.

Inputs needed:

- old/reference GFF3: the same file passed to `main.py --ref-gff`
- new/target GFF3: the same file passed to `main.py --target-gff`
- mapped GFF3: the file produced by `main.py --output-gff`
- mapping report: the file produced by `main.py --report`; if `--report` is
  not passed, the converter looks for `report.txt` next to `--mapped-gff`
- `mapping_session_id`: the core-DB mapping session ID to write into
  `stable_id_event`
- gene, transcript, and translation stable-ID ranges, passed as
  `PREFIX:START-END`
- output SQL path
- optional output TSV path for review

The generated executable SQL:

1. backs up the core tables it will touch inside the same DB,
2. stages gene/transcript/translation stable-ID decisions,
3. resolves target rows to internal DB IDs before applying updates,
4. updates stable IDs with a two-phase temporary rename to avoid cascading
   `A -> B`, `B -> C` collisions,
5. inserts rows into `stable_id_event`.

It does not update exons.

#### Get the Stable-ID Range

The ID ranges come from `gb1` / `gb_assembly_metadata`. Query them with the GCA
accession for the assembly:

```sql
SELECT
    species_log.gca_accession,
    CONCAT(prefix.prefix, ':', stable.stable_space_start, '-', stable.stable_space_end) AS id_range
FROM stable_space_species_log AS species_log
INNER JOIN stable_space AS stable ON stable.stable_space_id = species_log.stable_space_id
INNER JOIN species_prefix AS prefix ON prefix.lowest_taxon_id = species_log.lowest_taxon_id
WHERE species_log.gca_accession = '<GCA_accession>';
```

Pass the returned ranges to the script as `--gene-range`,
`--transcript-range`, and `--translation-range`. `START` and `END` may include
leading zeroes; the widest width is preserved when generating IDs.

#### Run from `main.py` Output

```bash
python3 main_output_to_stable_id_event_sql.py \
  --ref-gff ref.gff3 \
  --target-gff tar.gff3 \
  --mapped-gff out/root-smoke/mapped.gff3 \
  --report out/root-smoke/report.txt \
  --mapping-session-id 1 \
  --gene-range ENSGALG:90000000000-90000099999 \
  --transcript-range ENSGALT:90000000000-90000099999 \
  --translation-range ENSGALP:90000000000-90000099999 \
  --include-translations \
  --backup-prefix stable_id_mapper_backup_test \
  --output-sql out/root-smoke/core_updates.sql \
  --output-tsv out/root-smoke/core_updates.tsv
```

The script trusts `main.py`'s mapped GFF3 as the source of successfully mapped
old IDs. It classifies the reference and target GFF3s by hierarchy: top-level
features are genes, children of genes are transcripts, and `CDS` records are
translations. It does not update exons. Translation mappings are inferred only
when one old and one target translation sit under a mapped transcript;
otherwise target translations receive new IDs from `--translation-range`.

Review the TSV before running the SQL against a core database. The SQL also
contains `SELECT` count checks for staged rows versus DB-matched rows.

Add `--dry-run` to generate validation-only SQL:

```bash
python3 main_output_to_stable_id_event_sql.py \
  --ref-gff ref.gff3 \
  --target-gff tar.gff3 \
  --mapped-gff out/root-smoke/mapped.gff3 \
  --report out/root-smoke/report.txt \
  --mapping-session-id 1 \
  --gene-range ENSGALG:90000000000-90000099999 \
  --transcript-range ENSGALT:90000000000-90000099999 \
  --translation-range ENSGALP:90000000000-90000099999 \
  --include-translations \
  --dry-run \
  --output-sql out/root-smoke/core_updates.dry_run.sql \
  --output-tsv out/root-smoke/core_updates.tsv
```

Dry-run SQL only creates temporary tables and emits `SELECT` checks. It does
not create backup tables, update core feature tables, delete previous events, or
insert into `stable_id_event`.

## Workflow 2: LiftOn-to-Reference Mapper

Use `lifton_id_mapper.py` when you have:

- a LiftOn-projected GFF3 as the query
- a reference or Ensembl GFF3 on the same target assembly

Show the command-line options:

```bash
python3 lifton_id_mapper.py --help
```

Run the mapper:

```bash
mkdir -p out/lifton

python3 lifton_id_mapper.py \
  --lifton liftOn.gff3 \
  --reference ensembl_on_target.gff3 \
  --out-prefix out/lifton/mapping \
  --window 100000 \
  --topk 5 \
  --min-score 0.60 \
  --good 0.75 \
  --confident 0.85 \
  --gene-fraction 0.60 \
  --rewrite-reference \
  --rename-mode alias
```

Outputs:

```text
out/lifton/mapping.transcript_pairs.tsv
out/lifton/mapping.gene_pairs.tsv
out/lifton/mapping.mapped_reference.gff3
```

`--rename-mode alias` keeps reference IDs unchanged and adds
`lifton_gene_id` / `lifton_transcript_id` attributes.

`--rename-mode rename` replaces reference IDs with LiftOn IDs and updates
child `Parent` attributes. Use this mode only after reviewing the mapping
tables.

### LiftOn Mapper Algorithm

`lifton_id_mapper.py` performs structure-based transcript matching:

1. Parse both GFF3 files into genes and transcripts.
2. Build an interval index over reference genes by contig and strand.
3. For each LiftOn transcript, find candidate reference transcripts on the same
   contig and strand within `--window`.
4. Score candidate transcript pairs using:
   - intron-chain similarity, weight `0.50`
   - internal-exon Jaccard similarity, weight `0.25`
   - all-exon Jaccard similarity, weight `0.10`
   - exon-count similarity, weight `0.05`
   - TSS/TES boundary similarity, weight `0.05`
   - optional LiftOn `protein_identity` / `dna_identity` prior, weight `0.05`
5. Keep the top `--topk` candidates per LiftOn transcript.
6. Greedily select one-to-one transcript mappings by descending score.
7. Aggregate transcript mappings into gene mappings.
8. Optionally rewrite the reference GFF3 in alias or rename mode.

### LiftOn Mapper Smoke Test

A useful self-test is to map `tar.gff3` against itself:

```bash
mkdir -p out/lifton-self

python3 lifton_id_mapper.py \
  --lifton tar.gff3 \
  --reference tar.gff3 \
  --out-prefix out/lifton-self/self \
  --rewrite-reference \
  --rename-mode alias
```

Expected output for the bundled fixture:

```text
LiftOn transcripts with exons: 1977
Mapped transcripts (one-to-one): 1977
Top score: 0.950, median: 0.950
Gene pairs: 838 (fraction threshold 0.6)
```

The maximum score is `0.950` rather than `1.000` because the optional LiftOn
identity prior is not present in `tar.gff3`, and that prior contributes the
final `0.05` of the score.

## Testing

There are currently no committed unit tests for these stable-ID mapping
scripts.

Basic checks:

```bash
python3 main.py --help
python3 lifton_id_mapper.py --help
python3 -m compileall .
```

If `pytest` is installed, run:

```bash
python3 -m pytest
```

At the time this README was written, the branch contains `pytest.ini` but no
test files, so `pytest` is expected to collect no tests unless tests are added.

## Known Limitations

The LiftOn mapper is currently the more self-contained and easier-to-test
workflow. The assembly-to-assembly mapper is useful as a prototype but needs
more validation before it should be used for release-quality GFF3 output.

Assembly mapper limitations:

- The assembly mapper GFF3 parser only handles exact feature types `gene`,
  `mRNA`, `transcript`, `exon`, and `CDS`. It does not currently preserve all
  Ensembl noncoding and pseudogene feature types.
- Exon and CDS attributes are not stored by `gff_io.py`, so exon/CDS metadata
  such as `exon_id`, `protein_id`, and phases may not be fully preserved.
- Generated transcript IDs are not always propagated to child exon/CDS `Parent`
  attributes in the output GFF3.
- Output rows are sorted by coordinate and feature type, so child rows can be
  emitted before their parent rows when they share the same start coordinate.
- `synteny_disambiguate()` currently computes neighbor scores but returns the
  input mappings unchanged.
- Without `minimap2`, exact-match fallback is too strict for most assembly
  updates.

LiftOn mapper limitations:

- Both inputs must already be on the same assembly coordinate system.
- Matching is greedy rather than a global optimum assignment.
- Biotype is recorded as metadata but ignored for scoring.
- `--rename-mode rename` is intentionally invasive and should be used only
  after reviewing `*.transcript_pairs.tsv` and `*.gene_pairs.tsv`.

## Suggested Next Work

Good next steps for making this branch production-ready:

1. Add unit tests for ID parsing, new-ID allocation, GFF3 parsing, scoring, and
   parent/child ID consistency.
2. Broaden `gff_io.py` support for Ensembl feature types and preserve exon/CDS
   attributes.
3. Fix assembly mapper output so generated transcript IDs are consistently used
   by child features.
4. Use `minimap2` in smoke tests and record expected mapping summaries for the
   bundled fixtures.
5. Decide whether these scripts should stay as a standalone pipeline or move
   under `src/python/ensembl/genes/stable_id`.

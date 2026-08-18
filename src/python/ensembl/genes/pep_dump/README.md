# Ensembl protein FASTA dump reproducer

Pure-Python reproduction of the Ensembl production protein FASTA dump
(`pep.fa`) straight from an Ensembl Core MySQL/MariaDB database.

The reference implementation is the Perl pipeline
`Bio::EnsEMBL::Production::Pipeline::FileDump::Geneset_FASTA`. This project
reimplements its `pep` output closely enough to be byte-identical, without
Perl, BioPerl, the Ensembl Perl API, or Biopython.

## The two tools

| | |
|---|---|
| **`ensembl_pep_dump.py`** | the dump tool — reads a Core DB, writes `pep.fa` |
| **`compare_pep_fastas.py`** | the validation tool — compares a generated FASTA against an official one |

They are strictly separate. The dump tool does not compare, download, or
touch the network beyond its database connection. The comparison tool does
not touch a database.

## What it does

* Selects top-level slices exactly as the Perl API does, splits them into
  reference chromosomes / reference non-chromosomal / non-reference, and
  writes the first two in descending length order.
* Rebuilds each transcript's spliced sequence from exons, applies `_rna_edit`
  transcript attributes, extracts the CDS, applies start-exon phase padding,
  trims a trailing partial codon, translates with the seq_region's NCBI codon
  table, drops the terminal stop, applies the initial-Met rewrite, and applies
  `initial_met` / `_selenocysteine` / `amino_acid_sub` / `_stop_codon_rt`
  translation SeqEdits.
* Emits Ensembl `pep` headers and wraps sequence lines at exactly 60
  characters.
* Optionally also writes `pep-including_alt.fa`, the same records followed by
  the non-reference (ALT and PATCH) ones.

Details and the Perl sources each rule comes from are in
[docs/perl_behaviour.md](docs/perl_behaviour.md).

## What it does not do

* No `cdna.fa`, `cds.fa`, GTF, GFF3 or EMBL output — protein only.
* No BLAST index generation.
* No FTP downloads, no comparison logic, no Perl shell-out.
* No support for genebuilds stored below the top-level coord system. Exons
  must live on the same seq_region as their slice, which is true for any
  database with `genebuild.level = toplevel`. Anything else exits with code 4
  and a clear message rather than emitting wrong sequence.

## Validated databases

| Database | Assembly | Result |
|---|---|---|
| `homo_sapiens_gca009914755v4_core_110_1` | T2T-CHM13v2.0 | exact parity, 104,957 records, byte-identical |
| `homo_sapiens_core_116_38` | GRCh38 | exact parity, 382,428 records, byte-identical (against `--alt-output`) |

Species other than human are untested. The slice and translation logic is
species-agnostic and reads the codon table per seq_region, so other
vertebrates are expected to work, but nothing beyond the two databases above
has been verified.

## Install

```bash
python3 -m pip install -r requirements.txt
```

Python 3.8+ and PyMySQL are the only requirements. 


## Running a dump

```bash
python3 ensembl_pep_dump.py \
  --host mysql-host \
  --port 3306 \
  --user ensro \
  --password-env MYSQL_PWD \
  --database homo_sapiens_gca009914755v4_core_110_1 \
  --output pep_generated.fa
```

Other options:

| Option | Purpose |
|---|---|
| `--alt-output PATH` | also write `pep-including_alt.fa`: the same records followed by the non-reference ones |
| `--start-codon-set {ncbi,atg-only}` | which codons initiate translation under NCBI table 1 (see below) |
| `--species-id N` | `coord_system.species_id` for multi-species collection DBs |
| `--connect-timeout N` | connection timeout in seconds (default 30) |
| `--debug-id STABLE_ID` | print diagnostics for one transcript or translation to stderr while dumping |
| `--no-progress` | suppress the stderr progress line, for batch logs |

`--debug-id` only prints; it never changes the output file and never compares
anything. It accepts a transcript or translation stable ID, with or without
the version:

```bash
python3 ensembl_pep_dump.py ... --output /dev/null --debug-id ENSP05220066032
```

### `--start-codon-set`: read this before dumping a new database

NCBI table 1 nominates `TTG`, `CTG` and `ATG` as initiators, and Ensembl
rewrites the first residue to `M` for any of them. Whether the production
pipeline honoured `TTG` and `CTG` has changed over time, and the two official
files this tool is validated against disagree:

| Start set | T2T geneset 2022_07 | GRCh38 geneset 2026_04 |
|---|---|---|
| `ncbi` (default) | 246 spurious initial `M` | exact parity |
| `atg-only` | exact parity | 230 records short an initial `M` |

No single setting reproduces both. `ncbi` matches current Ensembl production
and is the default; `atg-only` is needed for the 2022_07 T2T file. Both
wrapper scripts set the flag explicitly so neither depends on the default.
Only NCBI table 1 is affected — table 2's alternative starts are honoured
either way. Details in [docs/perl_behaviour.md](docs/perl_behaviour.md) §9.

### Exit codes

| Code | Meaning |
|---|---|
| 0 | dump completed |
| 1 | usage error |
| 2 | database connection failure |
| 3 | missing tables or required metadata |
| 4 | unrecoverable error during the dump |


## Running parity checks

```bash
python3 compare_pep_fastas.py \
  --expected t2t_chm13.pep.fa.gz \
  --observed t2t_chm13.generated.pep.fa \
  --report t2t_chm13.parity_report.tsv \
  --summary t2t_chm13.parity_summary.json
```


### Interpreting a parity report

Exit code:

| Code | Meaning |
|---|---|
| 0 | exact parity — same records, same order, byte-identical decompressed content |
| 1 | every record matches, but the bytes differ — a formatting problem such as line width or trailing newlines |
| 2 | a biological, header or ordering mismatch |
| 3 | a file could not be read or parsed, or the arguments were invalid |

The **JSON summary** carries every count, both input paths, a UTC timestamp,
and `passed`. The **TSV report** contains a header row plus one row per
mismatch and nothing else, so an empty-but-for-the-header report means clean.

TSV categories:

| Category | Meaning |
|---|---|
| `missing` | in the expected file, absent from the generated one |
| `extra` | generated but not expected |
| `header_mismatch` | same record, different header; `detail` gives the first differing character offset |
| `sequence_mismatch` | same record, different peptide; `detail` gives the 1-based first differing residue and both lengths, and the columns give ±10 residues of context |
| `order_mismatch` | the first position where the shared records appear in a different order |
| `duplicate_in_*` | a record ID appears more than once in that file |

Read them in that order. A pile of `sequence_mismatch` rows sharing a first
differing position of 1 points at initial-Met handling; a `first_diff_char`
clustered in the coordinate field of headers points at coding-region bounds;
a single `order_mismatch` at index 0 of a slice points at slice ordering.


## Known limitations

* Protein only; no cDNA or CDS output.
* Validated on two human databases. Other species are plausible but unproven.
* Exons must be on the same seq_region as their top-level slice.
* Slice ordering among equal-length slices depends on the order the server
  returns `seq_region` rows in. This matches Perl's stable sort over the same
  query, but it is not guaranteed across servers with different query plans.
* The NCBI table 1 start set has to be chosen per target file; there is no
  setting that reproduces both validated dumps. See `--start-codon-set` above.
* The published GRCh38 protein FASTA carries the non-reference records
  appended, so parity for it is checked against `--alt-output`, not the
  reference-only `pep.fa`.



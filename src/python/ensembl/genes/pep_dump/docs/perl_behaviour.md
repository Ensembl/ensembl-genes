# Ensembl Perl behaviour reproduced by `ensembl_pep_dump.py`

Every rule below was derived from the Perl production code and confirmed
against official dumps. The relevant sources are vendored in the repository
root for reference: `Base.pm`, `Base_Filetype.pm`, `Geneset_FASTA.pm`,
`Slice.pm`, `SliceAdaptor.pm`, `SequenceAdaptor.pm`, `Transcript.pm`,
`TranscriptAdaptor.pm`, `Translation.pm`, `FASTASerializer.pm`.

The pipeline entry point is
`Bio::EnsEMBL::Production::Pipeline::FileDump::Geneset_FASTA`.

---

## 1. Which file this reproduces

`Geneset_FASTA.pm` writes `pep.fa` and, when non-reference regions exist, a
second `pep-including_alt.fa`. This tool reproduces `pep.fa` only.

The `pep.fa` published under
`ftp.ebi.ac.uk/pub/ensemblorganisms/<species>/<accession>/ensembl/geneset/<date>/`
is this pipeline's output.

`ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz`
is **not**. It comes from the older `FASTA::DumpFile` dumper and differs in
three visible ways:

| | Geneset_FASTA `pep.fa` | legacy `pep.all.fa` |
|---|---|---|
| location field | `GRCh38:1:13893666:14304658:1` | `chromosome:GRCh38:1:...` |
| gene description | quoted, `description:"..."` | unquoted |
| record set | reference regions only | reference plus ALT/patch |

Use the geneset file as the parity target.

One further wrinkle, established empirically. `Geneset_FASTA::run` writes
`pep.fa` from `[@$chr, @$non_chr]` and, when non-reference slices exist,
copies it to `pep-including_alt.fa` and appends the non-reference records.
The geneset directories publish only **one** protein FASTA, and for GRCh38 its
content is the `-including_alt` variant: 382,428 records, of which the last
11,513 are on ALT and PATCH regions (`HSCHR6_MHC_*`, `HG*_PATCH`, and so on),
appended after all 370,915 reference records. This tool's `--alt-output`
option reproduces that file.

T2T-CHM13v2.0 has no non-reference top-level slices at all, so its `pep.fa`
and `pep-including_alt.fa` would be identical and the distinction is
invisible there.

## 2. Slice selection

`Geneset_FASTA::run` calls `Base::get_slices($dba)`, which does:

```perl
my @slices = @{ $sa->fetch_all('toplevel', undef, 1, undef, undef) };
```

The positional arguments are `($cs_name, $cs_version, $include_non_reference,
$include_duplicates, $include_lrg)`. Two consequences follow from the
`SliceAdaptor::fetch_all` implementation:

* `include_non_reference = 1`, so `non_ref` regions are returned. They are
  separated out later and written to the `-including_alt` file, not `pep.fa`.
* `include_lrg` is `undef`, so **LRG seq_regions are excluded**. This matters:
  GRCh38 release 116 has 2,030 top-level seq_regions of which 1,324 are LRG.
  Missing this filter inflates the record set with duplicated LRG genes.
  T2T-CHM13v2.0 has no LRG regions, so the omission is invisible there.

`fetch_all` also restricts to `coord_system.species_id = $dba->species_id()`,
which matters only for multi-species collection databases.

Every slice returned spans its whole seq_region: `start = 1`,
`end = seq_region.length`, `strand = 1`.

## 3. GRCh38 chromosome Y / PAR filtering

```perl
if ($dba->species eq 'homo_sapiens') {
  if ($assembly eq 'GRCh38') {
    @slices = grep { !($_->seq_region_name() eq 'Y' &&
                       ($_->end() < 2781480 || $_->start() > 56887902)) } @slices;
  } else {
    $self->throw("Cannot filter PAR region for Human assembly $assembly");
  }
}
```

Two points that are easy to get wrong:

* The guard is on `$dba->species`, the production name. The T2T database's
  production name is `homo_sapiens_gca009914755v4`, not `homo_sapiens`, so the
  whole block — including the `throw` — is skipped for it. The T2T dump is not
  "a human dump that happens to skip PAR filtering"; the code never enters the
  branch at all.
* Because top-level slices always have `start = 1`, the `start > 56887902` arm
  can never fire, and `end < 2781480` reduces to "length shorter than the PAR1
  boundary". On GRCh38 release 116 the only top-level `Y` region is the full
  57,227,415 bp chromosome, so the filter removes nothing. It is implemented
  anyway, in `filter_grch38_y_par()`, for fidelity.

Release 116 GRCh38 also has an **empty `assembly_exception` table**, so no
PAR or haplotype sequence redirection is needed when building chromosome
sequence.

## 4. Slice ordering and classification

```perl
foreach my $slice (sort { $b->length <=> $a->length } @slices) {
  if ($slice->is_reference) {
    if ($slice->is_chromosome) { push @chr, $slice }
    else                       { push @non_chr, $slice }
  } else                       { push @non_ref, $slice }
}
```

* reference = **no** `non_ref` seq_region attribute
* chromosome = **has** a `karyotype_rank` seq_region attribute
* `pep.fa` receives `[@$chr, @$non_chr]`, in that order

Sorting is by descending length. Perl's `sort` has used a stable mergesort
since 5.8, and Python's `sorted(..., reverse=True)` is also stable, so slices
of equal length keep the order the `seq_region` query returned them in. Two
pairs of equal-length reference slices exist in GRCh38 release 116 (1,428 bp
and 1,048 bp); neither carries a coding gene.

## 5. Transcript ordering

`Geneset_FASTA::print_to_file` calls `$slice->get_all_Transcripts()`, which
reaches `TranscriptAdaptor::fetch_all_by_Slice` with the constraint
`t.is_current = 1` and **no `ORDER BY`**. The row order is therefore whatever
the server returns, which in practice is primary-key order.

This implementation makes that explicit with `ORDER BY t.transcript_id ASC`.

Evidence:

* **T2T-CHM13v2.0** — 104,957 records, 0 ordering mismatches, decompressed
  output byte-identical to the official `pep.fa`.
* **GRCh38 release 116** — 382,428 records, 0 ordering mismatches and 0
  header mismatches across the full record set. See
  `docs/handover_notes.md` for the recorded run.

Note that transcript_id order is *not* coordinate order. The first GRCh38
chromosome-1 record is `ENSP00000489835` at 13,893,666, not the
lowest-coordinate coding transcript, which is what an implementation ordering
by `seq_region_start` would emit.

## 6. Sequence retrieval

`SequenceAdaptor::fetch_by_Slice_start_end_strand` projects a slice down to
the sequence-level coord system and concatenates the pieces.

* **Sequence-level top-level regions** (T2T chromosomes) have their own row in
  `dna`.
* **Assembled regions** (GRCh38 chromosomes and scaffolds) are rebuilt from
  `assembly`, joining each component's `dna` row, reverse-complementing
  components with `ori = -1`. Positions no component covers stay `N`.

This tool materialises the whole active slice once per slice rather than
issuing a query per exon; the result is identical and the query count drops
from one-per-exon to one-per-slice.

Reverse-strand exons are reverse-complemented with an IUPAC-aware mapping that
preserves case.

## 7. Splicing and `_rna_edit`

`Transcript::spliced_seq` concatenates exon sequences in `exon_transcript.rank`
order, then applies `_rna_edit` transcript attributes.

Attribute values are `start end alt_seq` in 1-based inclusive cDNA
coordinates; a two-field value is a deletion. Edits are applied
**highest-coordinate first** so earlier edits do not shift the coordinates of
later ones.

## 8. CDS extraction and phase padding

`Transcript::translateable_seq`:

```perl
my $mrna  = $self->spliced_seq();
my $start = $self->cdna_coding_start();
my $end   = $self->cdna_coding_end();
$mrna = substr($mrna, $start - 1, $end - $start + 1);
my $start_phase = $self->translation->start_Exon->phase();
if ($start_phase > 0) { $mrna = "N" x $start_phase . $mrna }
```

`cdna_coding_start` walks exons in rank order accumulating lengths until it
reaches `translation.start_exon_id`, then adds `translation.seq_start`;
`cdna_coding_end` does the same with `end_exon_id` and `seq_end`. A
`_transl_start` / `_transl_end` transcript attribute overrides the computed
value outright.

`_rna_edit` length changes then shift the coordinates:

* coding **start** shifts for edits starting strictly before it (`< start`)
* coding **end** shifts for edits starting at or before one base past it
  (`<= end + 1`), so an edit abutting the 3' end still extends the CDS

Phase padding prepends 1 or 2 `N` characters, which translate to `X`.

## 9. Translation

`Transcript::translate` in order:

1. Read the codon table from the slice's `codon_table` seq_region attribute,
   defaulting to NCBI table 1.
2. Trim a trailing partial codon. `complete_codon` is false in the production
   pipeline, so the remainder is chopped rather than padded with `N`.
3. Return `undef` if nothing is left.
4. Translate with `Bio::Tools::CodonTable->translate($mrna, 0)`. The `0`
   marks the CDS as incomplete, so BioPerl performs no start-codon or
   terminal-stop handling of its own and does not reject internal stops —
   internal `*` must survive until translation SeqEdits are applied.
5. Drop the last residue if the final codon is a terminator.
6. Rewrite the first residue as `M` if it is not already `M` **and** the first
   codon is a start codon for the table in use.
7. Apply translation-level SeqEdits.

Ambiguous codons follow BioPerl: expand over IUPAC possibilities, collapse to
a single residue only when every expansion agrees, otherwise `X`.

### Start codons: the one behaviour that is not fixed

Start codons come from the table's NCBI `starts` string. NCBI table 1
nominates three initiators — `TTG`, `CTG` and `ATG` — and table 2 nominates
`ATT`, `ATC`, `ATA`, `ATG` and `GTG`.

Table 2's alternatives are honoured in both published dumps. Table 1's
`TTG` / `CTG` are **not consistent between them**, and this was established by
running the dump both ways against both files:

| Start set | T2T geneset 2022_07 | GRCh38 geneset 2026_04 |
|---|---|---|
| `ATG` only | exact parity | 230 records short an initial `M` |
| NCBI (`TTG`/`CTG`/`ATG`) | 246 records with a spurious initial `M` | exact parity |

The same gene shows it directly: translation `ENSP05220032019` in the T2T file
begins `LGFLYYSAWKL`, while its GRCh38 counterpart `ENSP00000432965` begins
`MGFLYYSAWKL`. Both CDSs start with the same non-`ATG` initiator; only the
dumping pipeline's BioPerl differs.

No single behaviour reproduces both files, so the table 1 start set is a
runtime choice: `--start-codon-set ncbi` (the default, matching current
Ensembl production) or `--start-codon-set atg-only` (matching the 2022_07 T2T
dump). Only table 1 is affected.

### `cds_start_NF` does not suppress the rewrite

`Transcript::translate` has **no** `cds_start_NF` check, and the published
data agrees: of the 230 GRCh38 release-116 records that depend on the
initial-Met rewrite, **224 are flagged `cds_start_NF`** and the official dump
rewrites every one of them.

An earlier iteration of this tool suppressed the rewrite for `cds_start_NF`
transcripts. That guard was provably dead under the `atg-only` start set — a
rewrite needs `peptide[0] != 'M'` together with a start-codon first codon, and
`ATG` always translates to `M` — so it never affected the T2T result, and it
would have broken GRCh38. It has been removed.

## 10. Translation SeqEdits

`Translation::modify_translation` applies these translation attributes:

* `initial_met`
* `_selenocysteine`
* `amino_acid_sub`
* `_stop_codon_rt`

Same `start end alt_seq` encoding as `_rna_edit`, same highest-coordinate-first
application order, but in amino-acid coordinates.

## 11. Header format

`Geneset_FASTA::header($transcript, 'pep')` produces, space separated:

```
<translation_stable_id>.<version>
pep
<coord_system_version>:<seq_region_name>:<coding_start>:<coding_end>:<strand>
gene:<gene_stable_id>.<version>
transcript:<transcript_stable_id>.<version>
gene_biotype:<gene_biotype>
transcript_biotype:<transcript_biotype>
gene_symbol:<display_xref_label>        # only if the gene has a display xref
description:"<gene_description>"        # only if the gene has a description
```

* A version of `0`, `NULL` or empty is omitted along with the dot.
* Coding coordinates are genomic, low value first, from
  `Transcript::coding_region_start` / `coding_region_end`. The strand is the
  transcript's, so a reverse-strand record still lists the lower coordinate
  first and ends `:-1`.
* Only the coord system **version** appears (`GRCh38`, `T2T-CHM13v2.0`), not
  its name.
* The description is wrapped in double quotes and emitted verbatim, including
  any embedded `[Source:...]` and `;parent_gene_display_xref=` fragments.

## 12. Line wrapping

`FASTASerializer` wraps at 60 characters. Each record is a `>` header line
followed by the sequence in 60-character lines, each newline-terminated, with
no blank separator between records.

## 13. Record inclusion

`print_to_file` writes a `pep` record only when both hold:

* `$transcript->translateable_seq ne ''`
* `$transcript->translate` returns a defined value

Transcripts with no translation, no exons, or an empty CDS are skipped.

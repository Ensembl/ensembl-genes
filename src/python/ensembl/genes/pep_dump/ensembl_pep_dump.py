#!/usr/bin/env python3
"""
ensembl_pep_dump.py - dump protein FASTA from an Ensembl Core database.

Pure-Python reimplementation of the Ensembl production protein FASTA dump
(Bio::EnsEMBL::Production::Pipeline::FileDump::Geneset_FASTA, 'pep' output).
Reads an Ensembl Core MySQL/MariaDB schema and writes a pep.fa-style file.

This tool only dumps. Validation against a reference FASTA lives in
compare_pep_fastas.py.

Exit codes:
    0  dump completed
    1  usage error (argparse)
    2  database connection failure
    3  missing tables or required metadata
    4  unrecoverable error during the dump
"""

import argparse
import os
import sys
import time

import pymysql
import pymysql.cursors

EXIT_OK = 0
EXIT_USAGE = 1
EXIT_CONNECTION = 2
EXIT_METADATA = 3
EXIT_DUMP = 4

FASTA_LINE_WIDTH = 60

# Tables the dump reads. Checked up-front so a partially imported database
# fails with EXIT_METADATA rather than mid-dump.
REQUIRED_TABLES = (
    'meta', 'coord_system', 'seq_region', 'seq_region_attrib', 'attrib_type',
    'assembly', 'dna', 'gene', 'transcript', 'transcript_attrib',
    'translation', 'translation_attrib', 'exon', 'exon_transcript', 'xref',
)


class CoreDbError(Exception):
    """Schema or metadata problem in the Core database."""


class DumpError(Exception):
    """Unrecoverable error while producing the dump."""


# ---------------------------------------------------------------------------
# Codon tables
# ---------------------------------------------------------------------------
#
# NCBI genetic codes in the layout used by BioPerl's Bio::Tools::CodonTable.
# The 'starts' string marks which codons count as translation initiators; the
# Ensembl API consults it when deciding whether to rewrite the first residue
# as M. Table 1 is given the official NCBI start set (TTG, CTG, ATG); see
# set_table1_start_codons() for why that is selectable.

CODON_TABLES = {}

_BASE1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_BASE2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
_BASE3 = "TCAG" * 16

assert len(_BASE1) == 64 and len(_BASE2) == 64 and len(_BASE3) == 64


def _build_codon_table(table_id, aa_string, starts_string):
    table = {}
    start_codons = set()
    stop_codons = set()
    for i in range(64):
        codon = _BASE1[i] + _BASE2[i] + _BASE3[i]
        aa = aa_string[i]
        table[codon] = aa
        if starts_string[i] == 'M':
            start_codons.add(codon)
        if aa == '*':
            stop_codons.add(codon)
    CODON_TABLES[table_id] = {
        'table': table,
        'starts': start_codons,
        'stops': stop_codons,
    }


# 1: Standard
_build_codon_table(1,
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M---------------M---------------M----------------------------")
# 2: Vertebrate Mitochondrial
_build_codon_table(2,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    "----------**--------------------MMMM----------**---M------------")
# 3: Yeast Mitochondrial
_build_codon_table(3,
    "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**----------------------MM----------------------------")
# 4: Mold, Protozoan, Coelenterate Mitochondrial
_build_codon_table(4,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--MM------**-------M------------MMMM---------------M------------")
# 5: Invertebrate Mitochondrial
_build_codon_table(5,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    "---M------**--------------------MMMM---------------M------------")
# 6: Ciliate, Dasycladacean and Hexamita Nuclear
_build_codon_table(6,
    "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--------------*----M--------------------------------------------")
# 9: Echinoderm and Flatworm Mitochondrial
_build_codon_table(9,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "----------**-----------------------M---------------M------------")
# 10: Euplotid Nuclear
_build_codon_table(10,
    "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**-----------------------M----------------------------")
# 11: Bacterial, Archaeal and Plant Plastid
_build_codon_table(11,
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**--*----M------------MMMM---------------M------------")
# 12: Alternative Yeast Nuclear
_build_codon_table(12,
    "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------**--*----M---------------M----------------------------")
# 13: Ascidian Mitochondrial
_build_codon_table(13,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    "---M------**----------------------MM---------------M------------")
# 14: Alternative Flatworm Mitochondrial
_build_codon_table(14,
    "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "-----------*-----------------------M----------------------------")
# 16: Chlorophycean Mitochondrial
_build_codon_table(16,
    "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "----------*---*----M--------------------------------------------")
# 21: Trematode Mitochondrial
_build_codon_table(21,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "----------**-----------------------M---------------M------------")
# 22: Scenedesmus obliquus Mitochondrial
_build_codon_table(22,
    "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "------*---*---*----M--------------------------------------------")
# 23: Thraustochytrium Mitochondrial
_build_codon_table(23,
    "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "--*-------**--*-----------------M--M---------------M------------")
# 24: Rhabdopleuridae Mitochondrial
_build_codon_table(24,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    "---M------**-------M---------------M---------------M------------")
# 25: Candidate Division SR1 and Gracilibacteria
_build_codon_table(25,
    "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M------**-----------------------M----------------------------")


_CODON_CACHE = {}

START_CODON_SETS = ('ncbi', 'atg-only')

# NCBI table 1 nominates TTG, CTG and ATG as initiators. Whether the Perl
# pipeline honoured the two alternatives has changed over time, and it is
# directly visible in the published dumps:
#
#   * GRCh38 geneset 2026_04 rewrites a leading L to M for TTG/CTG starts
#     (230 records), including transcripts flagged cds_start_NF.
#   * T2T-CHM13v2.0 geneset 2022_07 leaves those records as L (246 records).
#
# No single behaviour reproduces both files, so the start set is selectable.
# 'ncbi' matches current Ensembl production; 'atg-only' matches the 2022_07
# T2T dump. Only table 1 is affected: table 2's alternative starts are
# honoured in both published dumps.

_TABLE1_NCBI_STARTS = frozenset(CODON_TABLES[1]['starts'])


def set_table1_start_codons(mode):
    """Select which codons initiate translation under NCBI table 1."""
    if mode not in START_CODON_SETS:
        raise ValueError("unknown start codon set: {0}".format(mode))
    CODON_TABLES[1]['starts'] = (
        set(_TABLE1_NCBI_STARTS) if mode == 'ncbi' else {'ATG'}
    )
    _CODON_CACHE.clear()


IUPAC = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
    'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
}

def translate_codon(codon, codon_table_id=1):
    """
    Translate one codon. Ambiguous codons are expanded over their IUPAC
    possibilities and collapse to a single residue only when every expansion
    agrees, matching Bio::Tools::CodonTable->translate; otherwise 'X'.
    """
    key = (codon, codon_table_id)
    cached = _CODON_CACHE.get(key)
    if cached is not None:
        return cached

    norm = codon.upper().replace('U', 'T')
    ct = CODON_TABLES[codon_table_id]['table']

    aa = ct.get(norm)
    if aa is None:
        aas = set()
        bases = [IUPAC.get(b, [b]) for b in norm]
        if len(bases) == 3:
            for b1 in bases[0]:
                for b2 in bases[1]:
                    for b3 in bases[2]:
                        resolved = ct.get(b1 + b2 + b3)
                        if resolved is not None:
                            aas.add(resolved)
        aa = aas.pop() if len(aas) == 1 else 'X'

    _CODON_CACHE[key] = aa
    return aa


def translate_sequence(mrna, codon_table_id=1):
    """
    Translate an mRNA string codon by codon.

    Equivalent to Bio::Tools::CodonTable->translate($mrna, 0): the CDS is
    treated as incomplete, so no start-codon or terminal-stop special casing
    happens here. Trailing bases that do not form a full codon are ignored.
    """
    return ''.join(
        translate_codon(mrna[i:i + 3], codon_table_id)
        for i in range(0, len(mrna) - 2, 3)
    )


def is_start_codon(codon, codon_table_id=1):
    codon = codon.upper().replace('U', 'T')
    return codon in CODON_TABLES[codon_table_id]['starts']


def is_ter_codon(codon, codon_table_id=1):
    codon = codon.upper().replace('U', 'T')
    return codon in CODON_TABLES[codon_table_id]['stops']


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans(
    'ACGTacgtRYSWKMBDHVNryswkmbdhvn',
    'TGCAtgcaYRSWMKVHDBNyrswmkvhdbn'
)

_COMPLEMENT_BYTES = bytes.maketrans(
    b'ACGTacgtRYSWKMBDHVNryswkmbdhvn',
    b'TGCAtgcaYRSWMKVHDBNyrswmkvhdbn'
)


def reverse_complement(seq):
    """Reverse complement, preserving case and IUPAC ambiguity codes."""
    if isinstance(seq, (bytes, bytearray)):
        return bytes(seq).translate(_COMPLEMENT_BYTES)[::-1]
    return seq.translate(_COMPLEMENT)[::-1]


def _as_bytes(value):
    """DNA columns come back as str or bytes depending on driver settings."""
    if isinstance(value, (bytes, bytearray)):
        return bytes(value)
    return value.encode('latin1')


def apply_seq_edit(seq, start, end, alt_seq):
    """Replace the 1-based inclusive span [start, end] with alt_seq."""
    return seq[:start - 1] + alt_seq + seq[end:]


def parse_seq_edit(value):
    """
    Parse a SeqEdit attribute value ('start end alt_seq').

    A two-field value is a deletion. Returns (start, end, alt_seq) or None if
    the value is not a usable edit.
    """
    parts = value.split()
    if len(parts) >= 3:
        return int(parts[0]), int(parts[1]), parts[2]
    if len(parts) == 2:
        return int(parts[0]), int(parts[1]), ''
    return None


def collect_seq_edits(values):
    """
    Parse SeqEdit attribute values and order them for safe application.

    Edits are applied highest-coordinate first so that earlier edits do not
    shift the coordinates of later ones.
    """
    edits = [e for e in (parse_seq_edit(v) for v in values) if e is not None]
    edits.sort(key=lambda e: e[0], reverse=True)
    return edits


def stable_id_version(stable_id, version):
    """Format 'stable_id.version', dropping an absent or zero version."""
    if version:
        return "{0}.{1}".format(stable_id, version)
    return stable_id


def write_fasta_record(fh, header, sequence, line_width=FASTA_LINE_WIDTH):
    """
    Write one FASTA record in Bio::EnsEMBL::Utils::IO::FASTASerializer layout:
    a '>' header line then the sequence hard-wrapped at line_width.
    """
    fh.write('>')
    fh.write(header)
    fh.write('\n')
    for i in range(0, len(sequence), line_width):
        fh.write(sequence[i:i + line_width])
        fh.write('\n')


# ---------------------------------------------------------------------------
# Slice classification
# ---------------------------------------------------------------------------

# GRCh38 chromosome Y pseudoautosomal boundaries used by the Perl pipeline
# (Base.pm get_slices).
GRCH38_Y_PAR1_END = 2781480
GRCH38_Y_PAR2_START = 56887902


def filter_grch38_y_par(slices, species, assembly):
    """
    Drop chromosome Y regions that lie wholly outside the GRCh38 PAR bounds.

    Mirrors the grep in Base.pm get_slices(), which runs only when the Core DB
    species is literally 'homo_sapiens'. Assembly-specific databases such as
    homo_sapiens_gca009914755v4 (T2T-CHM13v2.0) have a different production
    name, so the Perl block - including its throw() on a non-GRCh38 human
    assembly - never fires for them.

    Top-level slices always span the whole seq_region (start == 1), so in
    practice this only removes Y seq_regions shorter than the PAR1 boundary.
    On GRCh38 release 116 the only top-level Y region is the full chromosome,
    which is retained; the filter is kept for fidelity with the Perl source.
    """
    if species != 'homo_sapiens':
        return list(slices)
    if assembly != 'GRCh38':
        raise CoreDbError(
            "Cannot filter PAR region for Human assembly {0}".format(assembly)
        )
    return [
        s for s in slices
        if not (
            s['name'] == 'Y'
            and (s['end'] < GRCH38_Y_PAR1_END or s['start'] > GRCH38_Y_PAR2_START)
        )
    ]


def classify_slices(slices, is_reference, is_chromosome):
    """
    Sort top-level slices by descending length and split them into
    (reference chromosomes, reference non-chromosomal, non-reference).

    Mirrors Base.pm get_slices(). The sort is stable, so slices of equal
    length keep the order the seq_region query returned them in, matching
    Perl's stable mergesort.
    """
    ordered = sorted(slices, key=lambda s: s['length'], reverse=True)

    chromosomes = []
    non_chromosomal = []
    non_reference = []
    for s in ordered:
        if is_reference(s):
            if is_chromosome(s):
                chromosomes.append(s)
            else:
                non_chromosomal.append(s)
        else:
            non_reference.append(s)
    return chromosomes, non_chromosomal, non_reference


# ---------------------------------------------------------------------------
# Database access
# ---------------------------------------------------------------------------

class EnsemblCoreDB:
    """Read-only access to an Ensembl Core schema."""

    def __init__(self, host, port, user, password, database,
                 defaults_file=None, species_id=None, connect_timeout=30):
        connect_kwargs = dict(
            database=database,
            charset='latin1',
            cursorclass=pymysql.cursors.DictCursor,
            connect_timeout=connect_timeout,
            read_timeout=3600,
            autocommit=True,
        )
        if defaults_file:
            connect_kwargs['read_default_file'] = os.path.expanduser(defaults_file)
        # pymysql lets an explicit argument win over the defaults file, so
        # host/port/user are only sent when they were actually resolved. That
        # lets a my.cnf supply them when the caller did not.
        for key, value in (('host', host), ('port', port), ('user', user),
                           ('password', password)):
            if value is not None:
                connect_kwargs[key] = value

        try:
            self.conn = pymysql.connect(**connect_kwargs)
        except pymysql.Error as exc:
            raise ConnectionError(str(exc))

        self.database = database
        self._attrib_types = {}
        self._coord_systems = {}
        self._seq_region_attribs = {}
        self._genes = {}
        self._xrefs = {}
        self._species_id = species_id

        # Per-slice state, reset by set_active_slice().
        self._active_seq_region_id = None
        self._active_length = 0
        self._active_sequence = None
        self._translations = {}
        self._exons = {}
        self._transcript_attribs = {}
        self._translation_attribs = {}

    def close(self):
        try:
            self.conn.close()
        except pymysql.Error:
            pass

    def _query(self, sql, params=None):
        with self.conn.cursor() as cur:
            cur.execute(sql, params or ())
            return cur.fetchall()

    def _query_one(self, sql, params=None):
        rows = self._query(sql, params)
        return rows[0] if rows else None

    # --- schema and metadata ---

    def check_schema(self):
        """Raise CoreDbError if any table the dump reads is absent."""
        rows = self._query(
            "SELECT table_name AS t FROM information_schema.tables "
            "WHERE table_schema = %s", (self.database,)
        )
        present = {r['t'].lower() for r in rows}
        missing = [t for t in REQUIRED_TABLES if t not in present]
        if missing:
            raise CoreDbError(
                "database {0} is missing required tables: {1}".format(
                    self.database, ', '.join(missing))
            )

    def meta_value(self, key):
        row = self._query_one(
            "SELECT meta_value FROM meta WHERE meta_key = %s "
            "ORDER BY meta_id LIMIT 1", (key,)
        )
        return row['meta_value'] if row else None

    @property
    def species_id(self):
        """
        The coord_system species_id to restrict slice queries to.

        Single-species Core DBs use 1; collection databases hold several, in
        which case one must be named explicitly.
        """
        if self._species_id is None:
            rows = self._query("SELECT DISTINCT species_id FROM coord_system")
            ids = sorted(r['species_id'] for r in rows if r['species_id'] is not None)
            if len(ids) == 1:
                self._species_id = ids[0]
            elif not ids:
                raise CoreDbError("coord_system table has no species_id values")
            else:
                raise CoreDbError(
                    "multi-species database: pass --species-id (found {0})".format(
                        ', '.join(str(i) for i in ids))
                )
        return self._species_id

    def _load_attrib_types(self):
        if not self._attrib_types:
            for row in self._query("SELECT attrib_type_id, code FROM attrib_type"):
                self._attrib_types[row['attrib_type_id']] = row['code']
        return self._attrib_types

    def coord_systems(self):
        if not self._coord_systems:
            for row in self._query(
                "SELECT coord_system_id, name, version, `rank`, attrib "
                "FROM coord_system"
            ):
                self._coord_systems[row['coord_system_id']] = row
        return self._coord_systems

    # --- slices ---

    def _seq_region_ids_with_attrib(self, code):
        return {
            r['seq_region_id'] for r in self._query(
                "SELECT sr.seq_region_id FROM seq_region sr "
                "JOIN seq_region_attrib sra USING (seq_region_id) "
                "JOIN attrib_type at USING (attrib_type_id) "
                "JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id "
                "WHERE at.code = %s AND cs.species_id = %s",
                (code, self.species_id)
            )
        }

    def fetch_toplevel_slices(self):
        """
        Top-level slices, equivalent to
        $sa->fetch_all('toplevel', undef, 1, undef, undef).

        include_non_reference is on, so non_ref regions are kept, but
        include_lrg is off, so LRG regions are excluded. Row order is the
        seq_region query's natural order, which classify_slices() relies on to
        break length ties the same way the Perl sort does.
        """
        rows = self._query(
            "SELECT sr.seq_region_id, sr.name, sr.length, sr.coord_system_id "
            "FROM seq_region sr "
            "JOIN seq_region_attrib sra USING (seq_region_id) "
            "JOIN attrib_type at USING (attrib_type_id) "
            "JOIN coord_system cs ON cs.coord_system_id = sr.coord_system_id "
            "WHERE at.code = 'toplevel' AND cs.species_id = %s",
            (self.species_id,)
        )
        lrg_ids = self._seq_region_ids_with_attrib('LRG')
        coord_systems = self.coord_systems()

        slices = []
        for row in rows:
            if row['seq_region_id'] in lrg_ids:
                continue
            cs = coord_systems.get(row['coord_system_id'])
            if cs is None:
                raise CoreDbError(
                    "seq_region {0} references unknown coord_system {1}".format(
                        row['name'], row['coord_system_id'])
                )
            slices.append({
                'seq_region_id': row['seq_region_id'],
                'name': row['name'],
                'length': row['length'],
                'start': 1,
                'end': row['length'],
                'coord_system_id': row['coord_system_id'],
                'cs_name': cs['name'],
                'cs_version': cs['version'] or '',
            })
        return slices

    def get_slices_for_dump(self):
        """Return (chromosomes, non_chromosomal, non_reference) for this DB."""
        slices = self.fetch_toplevel_slices()
        slices = filter_grch38_y_par(
            slices,
            self.meta_value('species.production_name'),
            self.meta_value('assembly.default'),
        )

        non_ref_ids = self._seq_region_ids_with_attrib('non_ref')
        chromosome_ids = self._seq_region_ids_with_attrib('karyotype_rank')

        return classify_slices(
            slices,
            is_reference=lambda s: s['seq_region_id'] not in non_ref_ids,
            is_chromosome=lambda s: s['seq_region_id'] in chromosome_ids,
        )

    def seq_region_attribs(self, seq_region_id):
        """All attributes of a seq_region as {code: [values]}."""
        cached = self._seq_region_attribs.get(seq_region_id)
        if cached is not None:
            return cached

        self._load_attrib_types()
        attribs = {}
        for row in self._query(
            "SELECT attrib_type_id, value FROM seq_region_attrib "
            "WHERE seq_region_id = %s", (seq_region_id,)
        ):
            code = self._attrib_types.get(row['attrib_type_id'])
            attribs.setdefault(code, []).append(row['value'])

        self._seq_region_attribs[seq_region_id] = attribs
        return attribs

    def codon_table_id(self, seq_region_id):
        """Codon table for a seq_region; NCBI table 1 unless overridden."""
        values = self.seq_region_attribs(seq_region_id).get('codon_table', [])
        return int(values[0]) if values else 1

    # --- per-slice preloading ---

    def set_active_slice(self, seq_region_id, length):
        """
        Point the per-slice caches at one seq_region and bulk-load its
        transcript metadata. Slice DNA is loaded lazily on first use.
        """
        self.release_active_slice()
        self._active_seq_region_id = seq_region_id
        self._active_length = length
        self._load_attrib_types()

        for row in self._query(
            "SELECT tl.translation_id, tl.transcript_id, tl.seq_start, tl.seq_end, "
            "       tl.start_exon_id, tl.end_exon_id, tl.stable_id, tl.version "
            "FROM translation tl "
            "JOIN transcript t ON tl.transcript_id = t.transcript_id "
            "WHERE t.seq_region_id = %s AND t.is_current = 1",
            (seq_region_id,)
        ):
            self._translations[row['transcript_id']] = row

        for row in self._query(
            "SELECT e.exon_id, e.seq_region_id, e.seq_region_start, e.seq_region_end, "
            "       e.seq_region_strand, e.phase, e.end_phase, e.stable_id, e.version, "
            "       et.transcript_id, et.rank "
            "FROM exon e "
            "JOIN exon_transcript et ON e.exon_id = et.exon_id "
            "JOIN transcript t ON et.transcript_id = t.transcript_id "
            "WHERE t.seq_region_id = %s AND t.is_current = 1 "
            "ORDER BY et.transcript_id, et.rank ASC",
            (seq_region_id,)
        ):
            self._exons.setdefault(row['transcript_id'], []).append(row)

        for row in self._query(
            "SELECT ta.transcript_id, ta.attrib_type_id, ta.value "
            "FROM transcript_attrib ta "
            "JOIN transcript t ON ta.transcript_id = t.transcript_id "
            "WHERE t.seq_region_id = %s AND t.is_current = 1",
            (seq_region_id,)
        ):
            code = self._attrib_types.get(row['attrib_type_id'])
            self._transcript_attribs.setdefault(
                row['transcript_id'], {}).setdefault(code, []).append(row['value'])

        for row in self._query(
            "SELECT tla.translation_id, tla.attrib_type_id, tla.value "
            "FROM translation_attrib tla "
            "JOIN translation tl ON tla.translation_id = tl.translation_id "
            "JOIN transcript t ON tl.transcript_id = t.transcript_id "
            "WHERE t.seq_region_id = %s AND t.is_current = 1",
            (seq_region_id,)
        ):
            code = self._attrib_types.get(row['attrib_type_id'])
            self._translation_attribs.setdefault(
                row['translation_id'], {}).setdefault(code, []).append(row['value'])

    def release_active_slice(self):
        """Drop per-slice caches, including the slice DNA."""
        self._active_seq_region_id = None
        self._active_length = 0
        self._active_sequence = None
        self._translations = {}
        self._exons = {}
        self._transcript_attribs = {}
        self._translation_attribs = {}

    def has_translations(self):
        return bool(self._translations)

    def transcripts_for_active_slice(self):
        """
        Current transcripts on the active slice in transcript_id order.

        The Perl API issues no ORDER BY here; transcript_id ascending is the
        order the official dumps are in. See docs/perl_behaviour.md.
        """
        return self._query(
            "SELECT t.transcript_id, t.gene_id, t.seq_region_id, "
            "       t.seq_region_start, t.seq_region_end, t.seq_region_strand, "
            "       t.stable_id, t.version, t.biotype, t.source, t.description "
            "FROM transcript t "
            "WHERE t.seq_region_id = %s AND t.is_current = 1 "
            "ORDER BY t.transcript_id ASC",
            (self._active_seq_region_id,)
        )

    def translation_for(self, transcript_id):
        return self._translations.get(transcript_id)

    def exons_for(self, transcript_id):
        return self._exons.get(transcript_id, [])

    def transcript_attribs(self, transcript_id, code):
        return self._transcript_attribs.get(transcript_id, {}).get(code, [])

    def translation_attribs(self, translation_id, code):
        return self._translation_attribs.get(translation_id, {}).get(code, [])

    # --- genes ---

    def gene(self, gene_id):
        cached = self._genes.get(gene_id)
        if cached is not None:
            return cached
        row = self._query_one(
            "SELECT gene_id, stable_id, version, biotype, description, "
            "       display_xref_id "
            "FROM gene WHERE gene_id = %s", (gene_id,)
        )
        self._genes[gene_id] = row
        return row

    def display_xref_label(self, xref_id):
        if xref_id is None:
            return None
        if xref_id not in self._xrefs:
            row = self._query_one(
                "SELECT display_label FROM xref WHERE xref_id = %s", (xref_id,)
            )
            self._xrefs[xref_id] = row['display_label'] if row else None
        return self._xrefs[xref_id]

    # --- DNA ---

    def _load_active_sequence(self):
        """
        Materialise the active slice's forward-strand DNA as bytes.

        Sequence-level top-level regions (for example the T2T chromosomes)
        have a row in dna. Assembled regions such as GRCh38 chromosomes are
        rebuilt from their components via the assembly table; gaps stay as N.
        """
        seq_region_id = self._active_seq_region_id
        row = self._query_one(
            "SELECT sequence FROM dna WHERE seq_region_id = %s", (seq_region_id,)
        )
        if row is not None:
            self._active_sequence = _as_bytes(row['sequence'])
            return

        buf = bytearray(b'N' * self._active_length)
        with self.conn.cursor(pymysql.cursors.SSDictCursor) as cur:
            cur.execute(
                "SELECT a.asm_start, a.asm_end, a.cmp_start, a.cmp_end, a.ori, "
                "       d.sequence "
                "FROM assembly a "
                "JOIN dna d ON d.seq_region_id = a.cmp_seq_region_id "
                "WHERE a.asm_seq_region_id = %s",
                (seq_region_id,)
            )
            for row in cur:
                component = _as_bytes(row['sequence'])[
                    row['cmp_start'] - 1:row['cmp_end']
                ]
                if row['ori'] == -1:
                    component = reverse_complement(component)
                span = row['asm_end'] - row['asm_start'] + 1
                if len(component) != span:
                    raise DumpError(
                        "assembly component for seq_region {0} at {1}-{2} has "
                        "length {3}, expected {4}".format(
                            seq_region_id, row['asm_start'], row['asm_end'],
                            len(component), span)
                    )
                buf[row['asm_start'] - 1:row['asm_end']] = component

        self._active_sequence = bytes(buf)

    def exon_sequence(self, exon):
        """Strand-corrected genomic sequence of one exon."""
        seq_region_id = exon['seq_region_id']
        start = exon['seq_region_start']
        end = exon['seq_region_end']

        if seq_region_id != self._active_seq_region_id:
            raise DumpError(
                "exon {0} is on seq_region {1}, outside the active slice {2}".format(
                    exon['exon_id'], seq_region_id, self._active_seq_region_id)
            )

        if self._active_sequence is None:
            self._load_active_sequence()

        seq = self._active_sequence[start - 1:end].decode('latin1')
        if len(seq) < (end - start + 1):
            seq += 'N' * ((end - start + 1) - len(seq))
        if exon['seq_region_strand'] == -1:
            seq = reverse_complement(seq)
        return seq


# ---------------------------------------------------------------------------
# Transcript processing
# ---------------------------------------------------------------------------

def build_spliced_seq(db, transcript_id, exons):
    """
    Exon sequences concatenated in rank order, with _rna_edit SeqEdits applied.
    Equivalent to Transcript::spliced_seq with edits enabled.
    """
    seq = ''.join(db.exon_sequence(exon) for exon in exons)
    for start, end, alt in collect_seq_edits(
        db.transcript_attribs(transcript_id, '_rna_edit')
    ):
        seq = apply_seq_edit(seq, start, end, alt)
    return seq


def _rna_edit_length_deltas(db, transcript_id):
    """(edit_start, length_change) per _rna_edit, highest coordinate first."""
    deltas = []
    for start, end, alt in collect_seq_edits(
        db.transcript_attribs(transcript_id, '_rna_edit')
    ):
        deltas.append((start, len(alt) - (end - start + 1)))
    return deltas


def compute_cdna_coding_start(db, transcript_id, exons, translation):
    """
    Position of the first coding base in cDNA coordinates.
    Matches Transcript::cdna_coding_start.
    """
    override = db.transcript_attribs(transcript_id, '_transl_start')
    if override:
        return int(override[0])

    start = 0
    for exon in exons:
        if exon['exon_id'] == translation['start_exon_id']:
            start += translation['seq_start']
            break
        start += exon['seq_region_end'] - exon['seq_region_start'] + 1

    for edit_start, delta in _rna_edit_length_deltas(db, transcript_id):
        if edit_start < start:
            start += delta
    return start


def compute_cdna_coding_end(db, transcript_id, exons, translation):
    """
    Position of the last coding base in cDNA coordinates.
    Matches Transcript::cdna_coding_end.
    """
    override = db.transcript_attribs(transcript_id, '_transl_end')
    if override:
        return int(override[0])

    end = 0
    for exon in exons:
        if exon['exon_id'] == translation['end_exon_id']:
            end += translation['seq_end']
            break
        end += exon['seq_region_end'] - exon['seq_region_start'] + 1

    # An edit starting exactly one base past the CDS still extends it, hence
    # the <= end + 1 test rather than < end.
    for edit_start, delta in _rna_edit_length_deltas(db, transcript_id):
        if edit_start <= end + 1:
            end += delta
    return end


def build_translateable_seq(db, transcript_id, exons, translation):
    """
    The CDS as Transcript::translateable_seq returns it: the coding span of
    the spliced transcript, prefixed with N padding for a non-zero start-exon
    phase so the reading frame starts on a codon boundary.
    """
    if translation is None:
        return ''

    cdna_start = compute_cdna_coding_start(db, transcript_id, exons, translation)
    cdna_end = compute_cdna_coding_end(db, transcript_id, exons, translation)
    if not cdna_start or not cdna_end:
        return ''

    mrna = build_spliced_seq(db, transcript_id, exons)
    cds = mrna[cdna_start - 1:cdna_end]

    start_exon = next(
        (e for e in exons if e['exon_id'] == translation['start_exon_id']), None
    )
    if start_exon is not None and start_exon['phase'] > 0:
        cds = 'N' * start_exon['phase'] + cds

    return cds


TRANSLATION_SEQ_EDIT_CODES = (
    'initial_met', '_selenocysteine', 'amino_acid_sub', '_stop_codon_rt',
)


def translate_transcript(db, transcript_id, exons, translation, seq_region_id,
                         translateable=None):
    """
    Peptide for one transcript, following Transcript::translate.

    Returns None where the Perl API returns undef.
    """
    if translation is None:
        return None

    mrna = (build_translateable_seq(db, transcript_id, exons, translation)
            if translateable is None else translateable)
    if not mrna:
        return None

    # complete_codon is off in the production pipeline, so a trailing partial
    # codon is trimmed rather than padded with N.
    delta = len(mrna) % 3
    if delta:
        mrna = mrna[:-delta]
    if not mrna:
        return None

    codon_table_id = db.codon_table_id(seq_region_id)
    first_codon = mrna[:3]
    last_codon = mrna[-3:]

    peptide = translate_sequence(mrna, codon_table_id)

    # The peptide never carries the terminal stop.
    if is_ter_codon(last_codon, codon_table_id):
        peptide = peptide[:-1]

    # Rewrite the first residue as M for a genuine initiator codon.
    # Transcript::translate applies this unconditionally: a cds_start_NF flag
    # does NOT suppress it. 224 of the 230 GRCh38 release-116 records that
    # depend on this branch are cds_start_NF, and the official dump rewrites
    # all of them.
    if peptide and peptide[0] != 'M' and is_start_codon(first_codon, codon_table_id):
        peptide = 'M' + peptide[1:]

    edit_values = []
    for code in TRANSLATION_SEQ_EDIT_CODES:
        edit_values.extend(db.translation_attribs(translation['translation_id'], code))
    for start, end, alt in collect_seq_edits(edit_values):
        peptide = apply_seq_edit(peptide, start, end, alt)

    return peptide


def coding_region_bounds(exons, translation):
    """
    Genomic (start, end) of the CDS, low coordinate first, as
    Transcript::coding_region_start / coding_region_end report them.
    """
    if translation is None:
        return None, None

    start_exon = next(
        (e for e in exons if e['exon_id'] == translation['start_exon_id']), None)
    end_exon = next(
        (e for e in exons if e['exon_id'] == translation['end_exon_id']), None)
    if start_exon is None or end_exon is None:
        return None, None

    if start_exon['seq_region_strand'] == 1:
        return (start_exon['seq_region_start'] + translation['seq_start'] - 1,
                end_exon['seq_region_start'] + translation['seq_end'] - 1)
    return (end_exon['seq_region_end'] - translation['seq_end'] + 1,
            start_exon['seq_region_end'] - translation['seq_start'] + 1)


def build_pep_header(db, transcript, gene, translation, exons,
                     cs_version, seq_region_name):
    """
    The pep FASTA header, as Geneset_FASTA::header($transcript, 'pep') builds it.
    """
    coding_start, coding_end = coding_region_bounds(exons, translation)

    parts = [
        stable_id_version(translation['stable_id'], translation['version']),
        'pep',
        "{0}:{1}:{2}:{3}:{4}".format(
            cs_version, seq_region_name, coding_start, coding_end,
            transcript['seq_region_strand']),
        "gene:" + stable_id_version(gene['stable_id'], gene['version']),
        "transcript:" + stable_id_version(transcript['stable_id'], transcript['version']),
        "gene_biotype:" + gene['biotype'],
        "transcript_biotype:" + transcript['biotype'],
    ]

    symbol = db.display_xref_label(gene['display_xref_id'])
    if symbol:
        parts.append("gene_symbol:" + symbol)
    if gene['description'] is not None:
        parts.append('description:"{0}"'.format(gene['description']))

    return ' '.join(parts)


# ---------------------------------------------------------------------------
# Dump
# ---------------------------------------------------------------------------

def dump_pep(db, output_path, alt_output_path=None, progress=True,
             debug_id=None):
    """
    Write the pep FASTA for this database.

    Follows Geneset_FASTA::run/print_to_file for the 'pep' output.
    output_path receives reference chromosomes then reference non-chromosomal
    regions, each sorted by descending length.

    alt_output_path, when given, receives the same records followed by the
    non-reference ones, reproducing pep-including_alt.fa. The Perl pipeline
    copies the finished pep.fa and appends to the copy; writing both handles
    in one pass is equivalent and avoids re-reading the sequence.
    """
    chromosomes, non_chromosomal, non_reference = db.get_slices_for_dump()
    reference = chromosomes + non_chromosomal
    if not reference:
        raise CoreDbError("no reference top-level slices found")

    started = time.time()
    alt_fh = None
    try:
        with open(output_path, 'w') as fh:
            if alt_output_path:
                alt_fh = open(alt_output_path, 'w')

            handles = [fh] + ([alt_fh] if alt_fh else [])
            reference_records = _dump_slices(
                db, reference, handles, progress, debug_id, 'reference')

            non_reference_records = 0
            if alt_fh and non_reference:
                non_reference_records = _dump_slices(
                    db, non_reference, [alt_fh], progress, debug_id,
                    'non-reference')
    finally:
        if alt_fh is not None:
            alt_fh.close()

    if progress:
        print("\r".ljust(78), end='', file=sys.stderr)
        print("\rWrote {0:,} records to {1} in {2:.1f}s".format(
            reference_records, output_path, time.time() - started),
            file=sys.stderr)
        if alt_output_path:
            print("Wrote {0:,} records to {1}".format(
                reference_records + non_reference_records, alt_output_path),
                file=sys.stderr)

    return reference_records


def _dump_slices(db, slices, handles, progress, debug_id, label):
    """Write pep records for a list of slices to every open handle."""
    total = 0
    for index, sl in enumerate(slices, start=1):
        seq_region_id = sl['seq_region_id']
        db.set_active_slice(seq_region_id, sl['length'])

        if progress:
            print("\r[{0} {1}/{2}] {3} ({4:,} bp), {5:,} records".format(
                label, index, len(slices), sl['name'], sl['length'], total
            ).ljust(78), end='', file=sys.stderr, flush=True)

        if db.has_translations():
            for transcript in db.transcripts_for_active_slice():
                transcript_id = transcript['transcript_id']
                translation = db.translation_for(transcript_id)
                if translation is None:
                    continue

                exons = db.exons_for(transcript_id)
                if not exons:
                    continue

                translateable = build_translateable_seq(
                    db, transcript_id, exons, translation)
                if not translateable:
                    continue

                peptide = translate_transcript(
                    db, transcript_id, exons, translation, seq_region_id,
                    translateable)
                if peptide is None:
                    continue

                gene = db.gene(transcript['gene_id'])
                if gene is None:
                    raise DumpError(
                        "transcript {0} references missing gene {1}".format(
                            transcript['stable_id'], transcript['gene_id'])
                    )

                if debug_id and _matches_debug_id(
                        debug_id, transcript, translation):
                    print_debug_record(
                        db, transcript, translation, exons, translateable,
                        peptide, sl)

                header = build_pep_header(
                    db, transcript, gene, translation, exons,
                    sl['cs_version'], sl['name'])
                for handle in handles:
                    write_fasta_record(handle, header, peptide)
                total += 1

        db.release_active_slice()

    return total


def _matches_debug_id(debug_id, transcript, translation):
    return debug_id in (
        transcript['stable_id'],
        stable_id_version(transcript['stable_id'], transcript['version']),
        translation['stable_id'],
        stable_id_version(translation['stable_id'], translation['version']),
    )


def print_debug_record(db, transcript, translation, exons, translateable,
                       peptide, sl):
    """Print diagnostics for one transcript to stderr. Used by --debug-id."""
    out = sys.stderr
    coding_start, coding_end = coding_region_bounds(exons, translation)
    transcript_id = transcript['transcript_id']

    print("\n" + "=" * 78, file=out)
    print("transcript      {0}".format(
        stable_id_version(transcript['stable_id'], transcript['version'])), file=out)
    print("translation     {0}".format(
        stable_id_version(translation['stable_id'], translation['version'])), file=out)
    print("slice           {0} ({1}), strand {2}".format(
        sl['name'], sl['cs_version'], transcript['seq_region_strand']), file=out)
    print("codon table     {0}".format(db.codon_table_id(sl['seq_region_id'])), file=out)
    print("coding region   {0}-{1}".format(coding_start, coding_end), file=out)
    print("cdna coding     {0}-{1}".format(
        compute_cdna_coding_start(db, transcript_id, exons, translation),
        compute_cdna_coding_end(db, transcript_id, exons, translation)), file=out)

    print("exons           {0}".format(len(exons)), file=out)
    for exon in exons:
        flags = ''
        if exon['exon_id'] == translation['start_exon_id']:
            flags += ' [start_exon offset={0}]'.format(translation['seq_start'])
        if exon['exon_id'] == translation['end_exon_id']:
            flags += ' [end_exon offset={0}]'.format(translation['seq_end'])
        print("  rank {0:<3} {1} {2}-{3}:{4} phase={5} end_phase={6}{7}".format(
            exon['rank'], stable_id_version(exon['stable_id'], exon['version']),
            exon['seq_region_start'], exon['seq_region_end'],
            exon['seq_region_strand'], exon['phase'], exon['end_phase'], flags),
            file=out)

    rna_edits = db.transcript_attribs(transcript_id, '_rna_edit')
    if rna_edits:
        print("_rna_edit       {0}".format(rna_edits), file=out)
    if db.transcript_attribs(transcript_id, 'cds_start_NF'):
        print("cds_start_NF    set", file=out)
    for code in TRANSLATION_SEQ_EDIT_CODES:
        values = db.translation_attribs(translation['translation_id'], code)
        if values:
            print("{0:<15} {1}".format(code, values), file=out)

    print("cds length      {0}".format(len(translateable)), file=out)
    print("cds 5'          {0}".format(translateable[:60]), file=out)
    print("cds 3'          {0}".format(translateable[-60:]), file=out)
    print("peptide length  {0}".format(len(peptide)), file=out)
    print("peptide N-term  {0}".format(peptide[:60]), file=out)
    print("peptide C-term  {0}".format(peptide[-60:]), file=out)
    print("=" * 78, file=out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def resolve_password(args):
    """
    Work out the MySQL password from the chosen authentication option.

    Returns None when the password should come from a defaults file or when
    the server needs no password, so pymysql is left to its own default.
    """
    if args.password is not None:
        return args.password
    if args.password_env:
        if args.password_env not in os.environ:
            raise CoreDbError(
                "environment variable {0} is not set".format(args.password_env)
            )
        return os.environ[args.password_env]
    if args.mysql_defaults_file:
        return None
    return ''


class _ArgumentParser(argparse.ArgumentParser):
    """argparse exits 2 on a usage error; 2 means connection failure here."""

    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(EXIT_USAGE, "{0}: error: {1}\n".format(self.prog, message))


def build_parser():
    parser = _ArgumentParser(
        prog='ensembl_pep_dump.py',
        description='Dump protein FASTA (pep.fa) from an Ensembl Core database.',
    )
    parser.add_argument('--host', default=None,
                        help='MySQL host (default 127.0.0.1, or from '
                             '--mysql-defaults-file)')
    parser.add_argument('--port', type=int, default=None,
                        help='MySQL port (default 3306, or from '
                             '--mysql-defaults-file)')
    parser.add_argument('--user', default=None,
                        help='MySQL user (default ensro, or from '
                             '--mysql-defaults-file)')
    parser.add_argument('--database', required=True, help='Core database name')
    parser.add_argument('--output', required=True, help='Output FASTA path')

    auth = parser.add_mutually_exclusive_group()
    auth.add_argument('--password', default=None,
                      help='MySQL password (visible in the process list; '
                           'prefer --password-env or --mysql-defaults-file)')
    auth.add_argument('--password-env', metavar='VAR',
                      help='name of an environment variable holding the password')
    auth.add_argument('--mysql-defaults-file', metavar='PATH',
                      help='my.cnf-style file to read credentials from')

    parser.add_argument('--alt-output', metavar='PATH',
                        help='also write pep-including_alt.fa here: the same '
                             'records followed by the non-reference ones')
    parser.add_argument('--start-codon-set', choices=START_CODON_SETS,
                        default='ncbi',
                        help="which codons initiate translation under NCBI "
                             "table 1: 'ncbi' honours TTG and CTG as well as "
                             "ATG (default, matches current Ensembl "
                             "production); 'atg-only' matches the 2022_07 "
                             "T2T-CHM13v2.0 dump")
    parser.add_argument('--species-id', type=int, default=None,
                        help='coord_system species_id (needed for collection DBs)')
    parser.add_argument('--connect-timeout', type=int, default=30,
                        help='connection timeout in seconds')
    parser.add_argument('--debug-id', metavar='STABLE_ID',
                        help='print diagnostics for one transcript or '
                             'translation stable ID while dumping')
    parser.add_argument('--no-progress', action='store_true',
                        help='suppress progress output on stderr')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    try:
        password = resolve_password(args)
    except CoreDbError as exc:
        print("error: {0}".format(exc), file=sys.stderr)
        return EXIT_METADATA

    # Without a defaults file there is nothing else to fall back on, so the
    # built-in defaults apply; with one, unset options are left to the file.
    host, port, user = args.host, args.port, args.user
    if not args.mysql_defaults_file:
        host = host or '127.0.0.1'
        port = port or 3306
        user = user or 'ensro'

    try:
        db = EnsemblCoreDB(
            host=host, port=port, user=user,
            password=password, database=args.database,
            defaults_file=args.mysql_defaults_file,
            species_id=args.species_id,
            connect_timeout=args.connect_timeout,
        )
    except ConnectionError as exc:
        print("error: cannot connect to {0}:{1}/{2}: {3}".format(
            host or 'localhost', port or 3306, args.database, exc),
            file=sys.stderr)
        return EXIT_CONNECTION

    try:
        db.check_schema()
        species = db.meta_value('species.production_name')
        assembly = db.meta_value('assembly.default')
        if species is None or assembly is None:
            raise CoreDbError(
                "meta table lacks species.production_name or assembly.default"
            )
        if not args.no_progress:
            print("{0}: {1} / {2}".format(args.database, species, assembly),
                  file=sys.stderr)

        set_table1_start_codons(args.start_codon_set)
        dump_pep(db, args.output,
                 alt_output_path=args.alt_output,
                 progress=not args.no_progress,
                 debug_id=args.debug_id)
    except CoreDbError as exc:
        print("error: {0}".format(exc), file=sys.stderr)
        return EXIT_METADATA
    except DumpError as exc:
        print("error: {0}".format(exc), file=sys.stderr)
        return EXIT_DUMP
    except pymysql.Error as exc:
        print("error: database error during dump: {0}".format(exc), file=sys.stderr)
        return EXIT_DUMP
    except OSError as exc:
        print("error: cannot write {0}: {1}".format(args.output, exc), file=sys.stderr)
        return EXIT_DUMP
    finally:
        db.close()

    return EXIT_OK


if __name__ == '__main__':
    sys.exit(main())

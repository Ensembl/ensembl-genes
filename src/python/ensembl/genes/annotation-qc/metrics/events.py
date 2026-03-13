
from dataclasses import dataclass
from typing import FrozenSet

from Bio.Data import CodonTable


@dataclass(frozen=True)
class GeneticCode:
	"""
	NCBI genetic code table definition.

	Attributes
	----------
	table_id : int
		NCBI translation table number.
	name : str
		Primary name of the table (e.g. "Standard", "Vertebrate Mitochondrial").
	stop_codons : frozenset[str]
		Codons that terminate translation under this code.
	"""
	table_id: int
	name: str
	stop_codons: FrozenSet[str]


def get_genetic_code(table_id: int) -> GeneticCode:
	"""
	Return a GeneticCode for the given NCBI translation table ID.

	Parameters
	----------
	table_id : int
		NCBI genetic code table number (e.g. 1 = standard, 2 = vertebrate mt).

	Raises
	------
	ValueError
		If the table ID is not recognised by Biopython.
	"""
	try:
		table = CodonTable.unambiguous_dna_by_id[table_id]
	except KeyError:
		raise ValueError(f"Unknown NCBI genetic code table: {table_id}")
	return GeneticCode(
		table_id=table_id,
		name=table.names[0],
		stop_codons=frozenset(table.stop_codons),
	)


STANDARD_CODE = get_genetic_code(1)

# Documented alternative initiation codons (near cognates).
# Position-3 variants (ATA, ATC, ATT) are NOT included: wobble-position
# substitutions are not recognised by the initiator tRNA.
NEAR_COGNATE_STARTS: FrozenSet[str] = frozenset({
    "CTG", "GTG", "TTG",  # position 1 variants (A→C/G/T)
    "ACG",                 # position 2 variant  (T→C)
})


@dataclass
class CDSMetrics:
	"""
	Metrics about a CDS sequence computed from nucleotide FASTA.

	Attributes
	----------
	cds_len_nt : int
		Length in nucleotides.
	start_codon : str
		First 3 nucleotides if available, else "".
	stop_codon : str
		Last 3 nucleotides if available, else "".
	start_status : str
		"canonical" (ATG), "near_cognate" (1-nt substitution from ATG),
		"non_cognate" (all other non-ATG codons), or "missing".
	stop_status : str
		"canonical" (valid stop for the genetic code), "noncanonical", or "missing".
	frame_ok : bool
		True if length divisible by 3.
	frame_error : bool
		True if length not divisible by 3.
	has_internal_stop : bool
		True if an in-frame stop codon exists before the terminal codon.
	"""
	cds_len_nt: int
	start_codon: str
	stop_codon: str
	start_status: str
	stop_status: str
	frame_ok: bool
	frame_error: bool
	has_internal_stop: bool


def compute_cds_metrics(seq_nt: str, genetic_code: GeneticCode = STANDARD_CODE) -> CDSMetrics:
	"""
	Compute QC metrics for a putative coding sequence (CDS) from nucleotides.

	Evaluates:
	  - CDS length (nt)
	  - start codon and status (ATG => canonical, 1-nt variant => near_cognate, other => non_cognate)
	  - stop codon and status (valid stop for genetic code => canonical)
	  - frame correctness (length % 3 == 0)
	  - internal in-frame stops (excluding terminal codon)

	Parameters
	----------
	seq_nt : str
		Nucleotide sequence (None-like treated as "").
	genetic_code : GeneticCode
		NCBI genetic code table to use. Defaults to the standard code (table 1).

	Returns
	-------
	CDSMetrics
		Computed metrics.
	"""
	seq = (seq_nt or "").upper()
	length_of_cds = len(seq)

	start_codon = seq[:3] if length_of_cds >= 3 else ""
	stop_codon = seq[-3:] if length_of_cds >= 3 else ""

	if length_of_cds < 3:
		start_status = "missing"
		stop_status = "missing"
	else:
		if start_codon == "ATG":
			start_status = "canonical"
		elif start_codon in NEAR_COGNATE_STARTS:
			start_status = "near_cognate"
		else:
			start_status = "non_cognate"
		stop_status = "canonical" if stop_codon in genetic_code.stop_codons else "noncanonical"

	frame_error = (length_of_cds % 3) != 0
	frame_ok = not frame_error

	has_internal_stop = False
	if length_of_cds >= 6:
		last_full_codon_start = (length_of_cds // 3) * 3 - 3
		for i in range(0, max(0, last_full_codon_start), 3):
			codon = seq[i:i + 3]
			if codon in genetic_code.stop_codons:
				has_internal_stop = True
				break

	return CDSMetrics(
		cds_len_nt=length_of_cds,
		start_codon=start_codon,
		stop_codon=stop_codon,
		start_status=start_status,
		stop_status=stop_status,
		frame_ok=frame_ok,
		frame_error=frame_error,
		has_internal_stop=has_internal_stop
	)


# Canonical splice-site dinucleotide pairs (donor–acceptor).
# GT-AG: U2 spliceosome, major class (~99 % of junctions)
# GC-AG: U2 spliceosome, minor class (~1 %)
# AT-AC: U12 spliceosome, minor class (<0.1 %)
CANONICAL_SPLICE_PAIRS: FrozenSet[str] = frozenset({"GT-AG", "GC-AG", "AT-AC"})


@dataclass
class SpliceJunctionMetrics:
	"""
	Metrics for a single splice junction derived from the intron sequence.

	Attributes
	----------
	donor_dinucleotide : str
		First two nucleotides of the intron (5' splice site), e.g. "GT".
		Empty string if the intron sequence is too short.
	acceptor_dinucleotide : str
		Last two nucleotides of the intron (3' splice site), e.g. "AG".
		Empty string if the intron sequence is too short.
	junction_type : str
		Donor and acceptor joined by a hyphen, e.g. "GT-AG", or "unknown"
		if either dinucleotide could not be extracted.
	is_canonical : bool
		True if the junction type is one of GT-AG, GC-AG, or AT-AC.
	"""
	donor_dinucleotide: str
	acceptor_dinucleotide: str
	junction_type: str
	is_canonical: bool


def assess_splice_junction(intron_seq: str) -> SpliceJunctionMetrics:
	"""
	Assess the canonicality of a splice junction from the intron sequence.

	Extracts the donor dinucleotide (first 2 nt) and acceptor dinucleotide
	(last 2 nt) of the intron, classifies the junction type, and flags
	whether it is canonical.

	Parameters
	----------
	intron_seq : str
		Nucleotide sequence of the intron (None-like treated as "").
		Must be at least 4 nt for non-overlapping donor/acceptor sites.

	Returns
	-------
	SpliceJunctionMetrics
		Computed splice junction metrics.
	"""
	seq = (intron_seq or "").upper()

	if len(seq) < 4:
		return SpliceJunctionMetrics(
			donor_dinucleotide="",
			acceptor_dinucleotide="",
			junction_type="unknown",
			is_canonical=False,
		)

	donor = seq[:2]
	acceptor = seq[-2:]
	junction_type = f"{donor}-{acceptor}"

	return SpliceJunctionMetrics(
		donor_dinucleotide=donor,
		acceptor_dinucleotide=acceptor,
		junction_type=junction_type,
		is_canonical=junction_type in CANONICAL_SPLICE_PAIRS,
	)

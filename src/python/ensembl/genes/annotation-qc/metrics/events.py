
from dataclasses import dataclass

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
        "canonical" (ATG), "non_ATG", or "missing".
    stop_status : str
        "canonical" (TAA/TAG/TGA), "noncanonical", or "missing".
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


def compute_cds_metrics(seq_nt: str) -> CDSMetrics:
    """
    Compute QC metrics for a putative coding sequence (CDS) from nucleotides.

    Evaluates:
      - CDS length (nt)
      - start codon and status (ATG => canonical)
      - stop codon and status (TAA/TAG/TGA => canonical)
      - frame correctness (length % 3 == 0)
      - internal in-frame stops (excluding terminal codon)

    Parameters
    ----------
    seq_nt : str
        Nucleotide sequence (None-like treated as "").

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
        start_status = "canonical" if start_codon == "ATG" else "non_ATG"
        stop_status = "canonical" if stop_codon in {"TAA", "TAG", "TGA"} else "noncanonical"

    frame_error = (length_of_cds % 3) != 0
    frame_ok = not frame_error

    has_internal_stop = False
    if length_of_cds >= 6:
        last_full_codon_start = (length_of_cds // 3) * 3 - 3
        for i in range(0, max(0, last_full_codon_start), 3):
            codon = seq[i:i + 3]
            if codon in {"TAA", "TAG", "TGA"}:
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
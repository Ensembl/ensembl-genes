# fasta_io.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

class FastaReader:
    """
    Lightweight FASTA reader with in-memory contigs and O(1) slicing.
    """
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self._seqs: Dict[str, str] = {}

    def _ensure_loaded(self) -> None:
        if self._seqs:
            return
        cur = None
        buf: list[str] = []
        with self.path.open() as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln:
                    continue
                if ln.startswith(">"):
                    if cur is not None:
                        self._seqs[cur] = "".join(buf).upper()
                    cur = ln[1:].split()[0]
                    buf = []
                else:
                    buf.append(ln)
            if cur is not None:
                self._seqs[cur] = "".join(buf).upper()

    @property
    def chromosomes(self) -> Iterable[str]:
        self._ensure_loaded()
        return self._seqs.keys()

    def get(self, chrom: str) -> str:
        self._ensure_loaded()
        return self._seqs[chrom]

    def slice(self, chrom: str, start_1: int, end_1: int) -> str:
        """Inclusive, 1-based slice."""
        self._ensure_loaded()
        s = self._seqs[chrom]
        if not (1 <= start_1 <= end_1 <= len(s)):
            raise ValueError(f"Invalid slice {chrom}:{start_1}-{end_1} (len={len(s)})")
        return s[start_1-1:end_1]

def revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                       "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
    return seq.translate(tr)[::-1]

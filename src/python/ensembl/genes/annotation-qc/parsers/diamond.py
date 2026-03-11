import pandas as pd


def parse_diamond_hits(path: str) -> pd.DataFrame:
    """Parse a DIAMOND tabular hits file into a DataFrame.

    Expects a headerless TSV with columns:
    qseqid, sseqid, pident, length, qstart, qend, qlen, qcovhsp,
    sstart, send, slen, evalue, bitscore
    """
    diamond_cols = [
        "qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "qcovhsp",
        "sstart", "send", "slen", "evalue", "bitscore"
    ]
    return pd.read_csv(path, sep="\t", header=None, names=diamond_cols)
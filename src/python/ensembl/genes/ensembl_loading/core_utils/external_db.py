"""Helpers for resolving production-controlled external databases."""

from __future__ import annotations

from typing import Any

REFSEQ_GENOMIC_DB_NAME = "RefSeq_genomic"
INSDC_DB_NAME = "INSDC"
REFSEQ_GENE_DB_NAME = "RefSeq_gene_name"
REFSEQ_TRANSCRIPT_DB_BY_PREFIX = {
    "NM_": "RefSeq_mRNA",
    "XM_": "RefSeq_mRNA_predicted",
    "NR_": "RefSeq_ncRNA",
    "XR_": "RefSeq_ncRNA_predicted",
}
REFSEQ_TRANSLATION_DB_BY_PREFIX = {
    "NP_": "RefSeq_peptide",
    "XP_": "RefSeq_peptide_predicted",
}


def get_external_db_id(cursor: Any, db_name: str) -> int:
    """Resolve an external database ID from the synchronized core table."""
    cursor.execute(
        "SELECT external_db_id FROM external_db WHERE db_name = %s LIMIT 1",
        (db_name,),
    )
    row = cursor.fetchone()
    if row is None:
        raise ValueError(
            f"Required external database {db_name!r} was not found in core.external_db"
        )
    value = row["external_db_id"] if isinstance(row, dict) else row[0]
    return int(value)


def get_refseq_external_db_name(accession: str, object_type: str) -> str:
    """Resolve the RefSeq external database name from an accession prefix."""
    mapping = (
        REFSEQ_TRANSCRIPT_DB_BY_PREFIX
        if object_type == "Transcript"
        else REFSEQ_TRANSLATION_DB_BY_PREFIX
    )
    for prefix, db_name in mapping.items():
        if accession.startswith(prefix):
            return db_name
    raise ValueError(f"Unsupported RefSeq {object_type} accession: {accession!r}")


def is_refseq_accession(accession: str, object_type: str) -> bool:
    """Return whether an accession has a supported RefSeq feature prefix."""
    mapping = (
        REFSEQ_TRANSCRIPT_DB_BY_PREFIX
        if object_type == "Transcript"
        else REFSEQ_TRANSLATION_DB_BY_PREFIX
    )
    return any(accession.startswith(prefix) for prefix in mapping)

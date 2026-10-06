"""Helpers for writing small GFF3/GTF fixtures."""

import gzip

import pytest


def gff3_text(records: list[tuple]) -> str:
    """records: (seqid, type, start, end, strand, attributes), one-based inclusive."""
    lines = ["##gff-version 3"]
    for seqid, feature, start, end, strand, attributes in records:
        lines.append(
            f"{seqid}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}"
        )
    return "\n".join(lines) + "\n"


def gene_records(
    seqid, gene, strand, transcripts, biotype="protein_coding", with_cds=True
):
    """
    Build gene/mRNA/exon/CDS records.
    transcripts: {transcript_id: [(start, end), ...]} exon coordinates; CDS equals exons
    unless with_cds is False or a dict {transcript_id: [(start, end), ...]}.
    """
    starts = [s for exons in transcripts.values() for s, _ in exons]
    ends = [e for exons in transcripts.values() for _, e in exons]
    records = [
        (seqid, "gene", min(starts), max(ends), strand, f"ID={gene};biotype={biotype}")
    ]
    for tid, exons in transcripts.items():
        records.append(
            (
                seqid,
                "mRNA",
                min(s for s, _ in exons),
                max(e for _, e in exons),
                strand,
                f"ID={tid};Parent={gene}",
            )
        )
        for start, end in exons:
            records.append((seqid, "exon", start, end, strand, f"Parent={tid}"))
        cds = (
            with_cds.get(tid, [])
            if isinstance(with_cds, dict)
            else (exons if with_cds else [])
        )
        for start, end in cds:
            records.append((seqid, "CDS", start, end, strand, f"Parent={tid}"))
    return records


@pytest.fixture
def write_gff3(tmp_path):
    def _write(name, records, compress=False):
        path = tmp_path / name
        text = gff3_text(records)
        if compress:
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return str(path)

    return _write


@pytest.fixture
def write_text(tmp_path):
    def _write(name, text, compress=False):
        path = tmp_path / name
        if compress:
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return str(path)

    return _write

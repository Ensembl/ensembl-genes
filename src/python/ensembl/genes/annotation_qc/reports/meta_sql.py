"""Report writers for annotation QC genebuild statistics."""

import json
import re
from typing import Any, Dict, List


def write_tsv(path: str, headers: List[str], row: Dict[str, object]) -> None:
    """
    Write a single-row TSV file with a specified header order.

    This function writes a tab-separated values (TSV) file containing a header
    row followed by exactly one data row. Output values are written in the
    order specified by ``headers``. Missing values (None) are written as empty
    fields.

    Parameters
    ----------
    path : str
        Path to the output TSV file.
    headers : List[str]
        Ordered list of column names to write as the header row.
    row : Dict[str, object]
        Mapping from column names to values. Values are converted to strings;
        keys missing from the mapping or explicitly set to None result in empty
        fields in the output.

    Returns
    -------
    None
    """
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        fh.write(
            "\t".join("" if row.get(h) is None else str(row.get(h)) for h in headers)
            + "\n"
        )


def write_json_from_tsv(tsv_path: str, json_path: str) -> None:
    """
    Convert a TSV file into a JSON array of row objects.

    The input TSV file is expected to contain a single header row followed by
    one or more data rows. Each data row is converted into a dictionary keyed
    by the header fields, and all rows are written as a JSON array.

    Empty fields in the TSV are preserved as empty strings in the resulting
    JSON, matching TSV semantics and avoiding implicit type coercion.

    Parameters
    ----------
    tsv_path : str
        Path to the input TSV file.
    json_path : str
        Path to the output JSON file.

    Returns
    -------
    None
    """

    rows: List[Dict[str, Any]] = []

    with open(tsv_path, "rt", encoding="utf-8") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if header is None:
                header = fields
                continue

            row = {}
            for k, v in zip(header, fields):
                # Keep empty strings as empty strings (matches TSV semantics)
                row[k] = v
            rows.append(row)

    with open(json_path, "wt", encoding="utf-8") as out:
        json.dump(rows, out, indent=2, sort_keys=False)


def write_meta_sql(
    out_path: str,
    dbname: str,
    species_id: str,
    coding_stats: Dict[str, object],
    noncoding_stats: Dict[str, object],
    pseudogene_stats: Dict[str, object],
    assembly_row: Dict[str, object],
) -> None:
    """
    Write SQL statements to populate the Ensembl-style `meta` table with
    gene build and assembly statistics.

    This function generates a SQL script containing `INSERT IGNORE` statements
    for numeric metadata values derived from coding, non-coding, pseudogene,
    and assembly statistics. The output is intended to be executed against a
    database that follows the Ensembl schema conventions.

    Only values that appear to be numeric (integers or floats, including
    negatives) are written. Empty values, missing keys, and non-numeric values
    are silently skipped.

    The resulting file:
      * Optionally begins with a `USE <dbname>;` statement
      * Contains one `INSERT IGNORE INTO meta (...)` statement per valid value
      * Never overwrites existing rows in the `meta` table

    Parameters
    ----------
    out_path : str
        Path to the output `.sql` file to be written. Any existing file at this
        path will be overwritten.

    dbname : str
        Name of the target database. If provided (non-empty), a `USE <dbname>;`
        statement is written at the top of the file. If empty, no `USE`
        statement is emitted.

    species_id : str
        Numeric species identifier used in the `meta.species_id` column.
        This value is written verbatim into the SQL statements and must be
        non-empty. A `ValueError` is raised if it is missing.

    coding_stats : Dict[str, object]
        Mapping of human-readable column names to values for coding gene
        statistics (e.g. "Coding genes", "Average exon length"). Only the keys
        referenced internally by this function are consulted.

    noncoding_stats : Dict[str, object]
        Mapping of human-readable column names to values for non-coding gene
        statistics.

    pseudogene_stats : Dict[str, object]
        Mapping of human-readable column names to values for pseudogene
        statistics.

    assembly_row : Dict[str, object]
        Mapping of human-readable column names to values for genome assembly
        statistics (e.g. N50, total genome length, GC percentage).

    Returns
    -------
    None
        This function does not return a value. Its sole side effect is writing
        a SQL script to `out_path`.

    Raises
    ------
    ValueError
        If `species_id` is empty or evaluates to False.

    Notes
    -----
    * Values are normalized by converting them to strings and stripping
      surrounding whitespace.
    * A value is considered valid only if it matches the regular expression
      for a signed integer or floating-point number.
    * SQL values are emitted without quotes, assuming numeric semantics.
    * The function intentionally uses `INSERT IGNORE` to avoid conflicts with
      existing meta keys already present in the database.

    Example
    -------
    >>> write_meta_sql(
    ...     out_path="meta.sql",
    ...     dbname="homo_sapiens_core",
    ...     species_id="1",
    ...     coding_stats={"Coding genes": 20000},
    ...     noncoding_stats={},
    ...     pseudogene_stats={},
    ...     assembly_row={"Contig N50": 50000000},
    ... )
    """

    if not species_id:
        raise ValueError("species_id is required to write meta SQL")

    coding_map = {
        "genebuild.stats.coding_genes": "Coding genes",
        "genebuild.stats.average_genomic_span": "Average genomic span",
        "genebuild.stats.average_sequence_length": "Average sequence length",
        "genebuild.stats.average_cds_length": "Average CDS length",
        "genebuild.stats.shortest_gene_length": "Shortest gene",
        "genebuild.stats.longest_gene_length": "Longest gene",
        "genebuild.stats.total_transcripts": "Total transcripts",
        "genebuild.stats.coding_transcripts": "Coding transcripts",
        "genebuild.stats.transcripts_per_gene": "Transcripts per gene",
        "genebuild.stats.coding_transcripts_per_gene": "Coding transcripts per gene",
        "genebuild.stats.total_exons": "Total exons",
        "genebuild.stats.total_coding_exons": "Total coding exons",
        "genebuild.stats.average_exon_length": "Average exon length",
        "genebuild.stats.average_coding_exon_length": "Average coding exon length",
        "genebuild.stats.average_coding_exons_per_transcript": "Average exons per transcript",
        "genebuild.stats.average_coding_exons_per_coding_transcript": "Average coding exons per coding transcript",
        "genebuild.stats.total_introns": "Total introns",
        "genebuild.stats.average_intron_length": "Average intron length",
    }

    noncoding_map = {
        "genebuild.stats.nc_non_coding_genes": "Non-coding genes",
        "genebuild.stats.nc_small_non_coding_genes": "Small non-coding genes",
        "genebuild.stats.nc_long_non_coding_genes": "Long non-coding genes",
        "genebuild.stats.nc_misc_non_coding_genes": "Misc non-coding genes",
        "genebuild.stats.nc_average_genomic_span": "Average genomic span",
        "genebuild.stats.nc_average_sequence_length": "Average sequence length",
        "genebuild.stats.nc_shortest_gene_length": "Shortest gene",
        "genebuild.stats.nc_longest_gene_length": "Longest gene",
        "genebuild.stats.nc_total_transcripts": "Total transcripts",
        "genebuild.stats.nc_transcripts_per_gene": "Transcripts per gene",
        "genebuild.stats.nc_total_exons": "Total exons",
        "genebuild.stats.nc_average_exon_length": "Average exon length",
        "genebuild.stats.nc_average_exons_per_transcript": "Average exons per transcript",
        "genebuild.stats.nc_total_introns": "Total introns",
        "genebuild.stats.nc_average_intron_length": "Average intron length",
    }

    pseudogene_map = {
        "genebuild.stats.ps_pseudogenes": "Pseudogenes",
        "genebuild.stats.ps_average_genomic_span": "Average genomic span",
        "genebuild.stats.ps_average_sequence_length": "Average sequence length",
        "genebuild.stats.ps_shortest_gene_length": "Shortest gene",
        "genebuild.stats.ps_longest_gene_length": "Longest gene",
        "genebuild.stats.ps_total_transcripts": "Total transcripts",
        "genebuild.stats.ps_transcripts_per_gene": "Transcripts per gene",
        "genebuild.stats.ps_total_exons": "Total exons",
        "genebuild.stats.ps_average_exon_length": "Average exon length",
        "genebuild.stats.ps_average_exons_per_transcript": "Average exons per transcript",
        "genebuild.stats.ps_total_introns": "Total introns",
        "genebuild.stats.ps_average_intron_length": "Average intron length",
    }

    assembly_map = {
        "assembly.stats.contig_n50": "Contig N50",
        "assembly.stats.total_genome_length": "Total genome length",
        "assembly.stats.total_coding_sequence_length": "Total coding sequence length",
        "assembly.stats.total_gap_length": "Total gap length",
        "assembly.stats.spanned_gaps": "Spanned gaps",
        "assembly.stats.chromosomes": "Chromosomes",
        "assembly.stats.toplevel_sequences": "Toplevel sequences",
        "assembly.stats.component_sequences": "Component sequences",
        "assembly.stats.gc_percentage": "% GC",
    }

    def normalize_value(v) -> str:
        if v is None:
            return ""
        s = str(v).strip()
        return s

    def is_empty(v: str) -> bool:
        return v == ""

    # Keep values if they look numeric, else we skip
    num_re = re.compile(r"^-?\d+(\.\d+)?$")

    def emit(fh, meta_key: str, raw_val: str) -> None:
        if is_empty(raw_val):
            return
        if not num_re.match(raw_val):
            return
        fh.write(
            "INSERT IGNORE INTO meta (species_id, meta_key, meta_value) "
            f"VALUES({species_id}, '{meta_key}', {raw_val});\n"
        )

    with open(out_path, "w", encoding="utf-8") as fh:
        if dbname:
            fh.write(f"USE {dbname};\n")

        for mk, col in coding_map.items():
            emit(fh, mk, normalize_value(coding_stats.get(col, "")))

        for mk, col in noncoding_map.items():
            emit(fh, mk, normalize_value(noncoding_stats.get(col, "")))

        for mk, col in pseudogene_map.items():
            emit(fh, mk, normalize_value(pseudogene_stats.get(col, "")))

        for mk, col in assembly_map.items():
            emit(fh, mk, normalize_value(assembly_row.get(col, "")))

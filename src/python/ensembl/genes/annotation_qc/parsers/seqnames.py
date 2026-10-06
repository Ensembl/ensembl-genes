"""
Read sequence-name maps and check annotation sequence names against a genome.

Map direction is always annotation name -> comparison name: the first column is
the name used in an annotation file, the second the name to use for the
comparison (normally the genome FASTA header). The same map is applied to the
reference and the query; names absent from the map are left unchanged.

Precedence when both sources are given: an explicit --seqname-map entry
overrides an NCBI assembly-report entry for the same name.
"""

import pandas as pd

_HEADER_FROM = {"from_seqname", "from", "seqname", "old", "old_name"}


def load_seqname_map(file_path: str) -> dict[str, str]:
    """
    Read a two-column seqname map (TSV, CSV or whitespace separated).
    Blank lines, '#' comments and a from/to header row are ignored.
    Args:
            file_path: Path to the map
    Returns:
            dict mapping annotation seqname to comparison seqname
    """
    mapping: dict[str, str] = {}
    with open(file_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            for separator in ("\t", ",", None):
                parts = [p.strip() for p in line.split(separator)]
                if len(parts) >= 2:
                    break
            if len(parts) < 2:
                raise ValueError(
                    f"{file_path}:{line_number}: expected two columns, got '{line}'"
                )
            source, target = parts[0], parts[1]
            if source.lower() in _HEADER_FROM and "to" in target.lower():
                continue
            if source in mapping and mapping[source] != target:
                raise ValueError(
                    f"{file_path}:{line_number}: '{source}' is mapped to both "
                    f"'{mapping[source]}' and '{target}'"
                )
            mapping[source] = target
    return mapping


def load_assembly_report(file_path: str) -> dict[str, str]:
    """
    Read an NCBI assembly report into GenBank accession -> assigned molecule.
    Rows whose assigned molecule is 'na' are skipped.
    Args:
            file_path: Path to *_assembly_report.txt
    Returns:
            dict mapping GenBank accession to chromosome name
    """
    mapping: dict[str, str] = {}
    with open(file_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[2] != "na":
                mapping[parts[4]] = parts[2]
    return mapping


def build_seqname_mapping(
    seqname_map: str | None = None, assembly_report: str | None = None
) -> dict[str, str]:
    """
    Merge an assembly report and a custom map; the custom map wins on collisions.
    Args:
            seqname_map: Optional custom map path
            assembly_report: Optional NCBI assembly report path
    Returns:
            dict mapping annotation seqname to comparison seqname (may be empty)
    """
    mapping: dict[str, str] = {}
    if assembly_report:
        mapping.update(load_assembly_report(assembly_report))
    if seqname_map:
        mapping.update(load_seqname_map(seqname_map))
    return mapping


def apply_seqname_mapping(
    annotation: pd.DataFrame, mapping: dict[str, str]
) -> tuple[pd.DataFrame, dict]:
    """
    Rename sequences in a comparison-schema DataFrame.

    Raises ValueError when two different annotation seqnames would be merged
    into one comparison seqname.
    Args:
            annotation: DataFrame with a Chromosome column
            mapping: annotation seqname -> comparison seqname
    Returns:
            (renamed DataFrame, report with mapped/unmapped seqnames)
    """
    names = sorted(annotation["Chromosome"].unique())
    mapped = {name: mapping[name] for name in names if name in mapping}
    targets = pd.Series({name: mapping.get(name, name) for name in names}, dtype=object)
    collisions = targets[targets.duplicated(keep=False)]
    if not collisions.empty:
        pairs = ", ".join(f"{s}->{t}" for s, t in collisions.items())
        raise ValueError(f"Seqname mapping merges distinct sequences: {pairs}")

    renamed = annotation.copy()
    if mapped:
        renamed["Chromosome"] = renamed["Chromosome"].map(
            lambda name: mapping.get(name, name)
        )
    report = {
        "seqnames": len(names),
        "mapped": mapped,
        "unmapped": [name for name in names if name not in mapping],
    }
    return renamed, report


def check_seqnames_against_genome(
    annotation: pd.DataFrame, sequence_lengths: dict[str, int]
) -> dict:
    """
    Compare annotation sequence names and coordinates with genome sequences.
    Args:
            annotation: Comparison-schema DataFrame (after seqname mapping)
            sequence_lengths: FASTA sequence name -> length
    Returns:
            dict with 'missing_seqnames' and 'out_of_bounds_seqnames' (each
            name -> gene count on that sequence), 'out_of_bounds_features' and
            'out_of_bounds_examples'
    """
    genes = annotation[annotation["Feature"] == "gene"]
    gene_counts = genes["Chromosome"].value_counts()
    missing = sorted(set(annotation["Chromosome"]) - set(sequence_lengths))
    lengths = annotation["Chromosome"].map(sequence_lengths)
    beyond = annotation[lengths.notna() & (annotation["End"] > lengths)]
    return {
        "missing_seqnames": {name: int(gene_counts.get(name, 0)) for name in missing},
        "out_of_bounds_seqnames": {
            name: int(gene_counts.get(name, 0))
            for name in sorted(set(beyond["Chromosome"]))
        },
        "out_of_bounds_features": len(beyond),
        "out_of_bounds_examples": [
            f"{r.Chromosome}:{r.Start + 1}-{r.End} ({r.Feature}; sequence length "
            f"{sequence_lengths[r.Chromosome]})"
            for r in beyond.head(3).itertuples()
        ],
    }

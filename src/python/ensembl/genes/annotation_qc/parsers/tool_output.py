"""
Read small tool-output files used alongside annotation comparisons.
"""

import json
from pathlib import Path

import pandas as pd


def read_evidence_attribution(file_path: str) -> pd.DataFrame:
    """
    Read a GMB evidence_attribution.tsv (one row per transcript).
    Args:
            file_path: Path to the TSV
    Returns:
            DataFrame; must contain a transcript_id column
    """
    table = pd.read_csv(file_path, sep="\t", dtype={"transcript_id": str})
    if "transcript_id" not in table.columns:
        raise ValueError(f"{file_path} has no transcript_id column")
    return table


def find_gmb_handover_manifest(annotation_path: str) -> tuple[str, dict] | None:
    """
    Find the GMB finalise handover_manifest.json that lists an annotation.

    GMB writes finalise/consensus.gff3 and
    finalise/canonical/consensus.canonical_annotated.gff3, with the manifest in
    finalise/. The annotation's directory and its parent are searched.
    Args:
            annotation_path: Path to a GMB finalise GFF3
    Returns:
            (manifest path, parsed manifest) or None when not found
    """
    annotation = Path(annotation_path).resolve()
    for directory in (annotation.parent, annotation.parent.parent):
        manifest = directory / "handover_manifest.json"
        if manifest.is_file():
            with open(manifest, "rt") as handle:
                return str(manifest), json.load(handle)
    return None

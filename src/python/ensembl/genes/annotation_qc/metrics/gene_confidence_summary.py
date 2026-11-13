# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Compute final confidence counts and percentages using unique scored genes."""
import pandas as pd

from ensembl.genes.annotation_qc.config.high_confidence import CLASSES


def summarize_gene_confidence(evidence: pd.DataFrame) -> dict[str, int | float]:
    """Count unique genes; reject missing labels and conflicting gene assignments.

    The denominator is genes present in the scored evidence, not all genes in the
    original annotation. Repeated identical gene/class assignments count once.
    Percentages are rounded to two decimal places; an empty input gives zeros.
    """
    required = {"gene_id", "confidence"}
    missing = required - set(evidence.columns)
    if missing:
        raise ValueError(f"Missing required evidence columns: {sorted(missing)}")
    assignments = evidence[["gene_id", "confidence"]].copy()
    valid_ids = assignments["gene_id"].map(
        lambda value: isinstance(value, str)
        and bool(value.strip())
        and value not in {".", "NA", "N/A"}
    )
    if not valid_ids.all():
        raise ValueError("Evidence contains missing or invalid gene_id values")
    if not assignments["confidence"].isin(CLASSES).all():
        raise ValueError("Confidence must be high, medium or low for every gene")
    assignments = assignments.drop_duplicates()
    conflicts = assignments.loc[
        assignments["gene_id"].duplicated(keep=False), "gene_id"
    ].unique()
    if len(conflicts):
        raise ValueError(
            f"Conflicting confidence assignments for genes: {', '.join(conflicts[:5])}"
        )
    total = len(assignments)
    counts = assignments["confidence"].value_counts()
    result: dict[str, int | float] = {"eam.total_genes": total}
    for category in CLASSES:
        count = int(counts.get(category, 0))
        key = f"eam.total_{category}_confidence_genes"
        result[key] = count
        result[key + "_pct"] = round(100 * count / total, 2) if total else 0.0
    return result

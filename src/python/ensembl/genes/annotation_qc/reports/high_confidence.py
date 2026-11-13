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
"""Export previously computed high-confidence evidence and annotations."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from urllib.parse import quote
from ensembl.genes.annotation_qc.config.high_confidence import HOMOLOGY_FIELDS
from ensembl.genes.annotation_qc.parsers.high_confidence import single_attr


from ensembl.genes.annotation_qc.config.high_confidence import VALIDATION_FIELDS
from ensembl.genes.annotation_qc.config.high_confidence import (
    NORMALIZATION_FIELDS,
    EXCEPTION_FIELDS,
)


def write_homology(path: Path, results: dict) -> None:
    """Write the legacy seven-column helper result format."""
    with path.open("w", encoding="utf-8") as handle:
        for query, hit in results.items():
            handle.write(
                "\t".join(
                    [query, hit["hit_id"]]
                    + [f"{hit[key]:.2f}" for key in HOMOLOGY_FIELDS[1:]]
                )
                + "\n"
            )


def cell(value):
    """Format missing and boolean values for evidence exports."""
    return "." if value is None else int(value) if isinstance(value, bool) else value


def write_fasta(handle, name: str, sequence: str) -> None:
    """Write a FASTA entry with lines wrapped at 60 bases."""
    handle.write(f">{name}\n")
    for index in range(0, len(sequence), 60):
        handle.write(sequence[index : index + 60] + "\n")


# The orchestration/layout keeps related report operations together.
# pylint: disable-next=too-many-arguments
def gff_line(handle, seqname, kind, start, end, strand, attrs, phase="."):
    """Write one escaped GFF3 row with inclusive coordinates."""

    def escape(value):
        return quote(str(value), safe=".:^*$@!+_?-|")

    attributes = ";".join(f"{key}={escape(value)}" for key, value in attrs.items())
    handle.write(
        "\t".join(
            map(
                str,
                [
                    escape(seqname),
                    "gtf_confidence",
                    kind,
                    start,
                    end,
                    ".",
                    strand,
                    phase,
                    attributes,
                ],
            )
        )
        + "\n"
    )


# The orchestration/layout keeps related report operations together.
# pylint: disable-next=too-many-locals
def write_gff(handle, gene, tx, result, confidence):
    # Type prefixes prevent collisions between gene, transcript and feature IDs.
    """Format precomputed gene, transcript and coding feature records."""
    gene_id, tx_id = f"gene:{gene.id}", f"transcript:{tx.id}"
    gff_line(
        handle,
        gene.seqname,
        "gene",
        gene.start,
        gene.end,
        gene.strand,
        {"ID": gene_id, "gene_id": gene.id, "gene_quality": confidence},
    )
    exons = result["exons"]
    transcript_attrs = {
        "ID": tx_id,
        "Parent": gene_id,
        "transcript_id": tx.id,
        "transl_table": result["evidence"]["translation_table"],
        "exon_structure_source": result["evidence"]["exon_structure_source"],
        "cds_phase_source": result["evidence"]["cds_phase_source"],
        "stop_codon_source": result["evidence"]["stop_codon_source"],
    }
    if result["evidence"]["selenocysteine_count"]:
        transcript_attrs.update(
            translation_status="selenocysteine_from_gtf",
            selenocysteine_count=result["evidence"]["selenocysteine_count"],
        )
    gff_line(
        handle,
        tx.seqname,
        "mRNA",
        min(f.start for f in exons),
        max(f.end for f in exons),
        tx.strand,
        transcript_attrs,
    )
    for index, exon in enumerate(exons, 1):
        attrs = {"ID": f"exon:{tx.id}:{index}", "Parent": tx_id}
        exon_id = single_attr(exon.attrs, "exon_id")
        if exon_id:
            attrs["exon_id"] = exon_id
        gff_line(handle, tx.seqname, "exon", exon.start, exon.end, tx.strand, attrs)
    for kind in ("five_prime_UTR", "three_prime_UTR"):
        for index, (start, end) in enumerate(result[kind], 1):
            gff_line(
                handle,
                tx.seqname,
                kind,
                start,
                end,
                tx.strand,
                {"ID": f"{kind}:{tx.id}:{index}", "Parent": tx_id},
            )
    for index, ((start, end), phase) in enumerate(
        zip(result["cds_spans"], result["gff_cds_phases"]), 1
    ):
        gff_line(
            handle,
            tx.seqname,
            "CDS",
            start,
            end,
            tx.strand,
            {"ID": f"CDS:{tx.id}:{index}", "Parent": tx_id},
            phase,
        )
    sec_features = sorted(
        (f for f in tx.features if f.kind.lower() == "selenocysteine"),
        key=lambda f: f.start,
        reverse=tx.strand == "-",
    )
    for index, feature in enumerate(sec_features, 1):
        gff_line(
            handle,
            tx.seqname,
            "Selenocysteine",
            feature.start,
            feature.end,
            tx.strand,
            {"ID": f"selenocysteine:{tx.id}:{index}", "Parent": tx_id},
        )
    for index, intron in enumerate(result["introns"], 1):
        gff_line(
            handle,
            tx.seqname,
            "intron",
            intron["start"],
            intron["end"],
            tx.strand,
            {"ID": f"intron:{tx.id}:{index}", "Parent": tx_id},
        )


def write_records(path, rows, fields):
    """Export records using the established 1/0/. evidence representation."""
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: cell(row.get(key)) for key in fields})


class HighConfidenceReports:
    """Own file creation/formatting; metric evaluation remains in the runner."""

    def __init__(self, paths):
        self.paths = paths
        paths["manifest"].parent.mkdir(parents=True, exist_ok=True)
        paths["manifest"].unlink(missing_ok=True)

    def write_validation(self, issues, policy):
        """Write coordinate diagnostics, including on failed runs."""

        write_records(
            self.paths["input_validation"],
            [dict(row, policy=policy) for row in issues],
            VALIDATION_FIELDS,
        )

    def write_extraction_diagnostics(self, normalization, exceptions):
        """Write already computed normalization and translation-exception records."""

        write_records(self.paths["normalization"], normalization, NORMALIZATION_FIELDS)
        write_records(
            self.paths["translation_exceptions"], exceptions, EXCEPTION_FIELDS
        )

    def write_sequences(self, results):
        """Export reconstructed sequences for homology search and downstream use."""
        with (
            self.paths["cds"].open("w", encoding="utf-8") as cds,
            self.paths["proteins"].open("w", encoding="utf-8") as proteins,
        ):
            for result in results:
                write_fasta(cds, result["tx"].id, result["cds"])
                write_fasta(proteins, result["tx"].id, result["protein"])

    def write_evidence(self, tables):
        """Export metric tables and formatted GFF3 without applying scoring rules."""
        for name in ["evidence", "introns"]:
            frame = tables[name]
            write_records(
                self.paths[name], frame.to_dict("records"), list(frame.columns)
            )
        with self.paths["gff"].open("w", encoding="utf-8") as handle:
            handle.write("##gff-version 3\n")
            for result in tables["gff_records"]:
                write_gff(
                    handle, result["gene"], result["tx"], result, result["confidence"]
                )

    def write_manifest(self, manifest):
        """Write the successful-run marker last."""
        self.paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")


def write_rescored_outputs(paths, result, manifest, companions=None, gff=None):
    """Write already classified evidence and update only GFF gene quality attributes."""
    paths["manifest"].parent.mkdir(parents=True, exist_ok=True)
    frame = result["evidence"]
    write_records(paths["evidence"], frame.to_dict("records"), list(frame.columns))
    write_records(
        paths["changes"],
        result["changes"],
        [
            "gene_id",
            "transcript_id",
            "old_confidence",
            "new_confidence",
            "old_rule",
            "new_rule",
            "class_changed",
        ],
    )
    if gff is not None:
        confidence = dict(zip(frame["gene_id"], frame["confidence"]))
        with paths["gff"].open("w", encoding="utf-8") as handle:
            for row in gff.to_dict("records"):
                if row["gene_id"] is None:
                    handle.write(row["raw_line"])
                else:
                    columns = list(row["columns"])
                    attrs = [
                        part
                        for part in columns[8].rstrip(";").split(";")
                        if not part.startswith("gene_quality=")
                    ]
                    columns[8] = ";".join(
                        attrs + [f"gene_quality={confidence[row['gene_id']]}"]
                    )
                    handle.write("\t".join(columns) + "\n")
    for target, source in (companions or {}).items():
        target.symlink_to(source.resolve())
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")

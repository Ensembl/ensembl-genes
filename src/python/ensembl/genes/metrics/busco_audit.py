#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
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

"""
BUSCO Evidence Audit Tool

This tool compares BUSCO genome mode and protein mode results to identify
problematic loci where gene models are missing or fragmented, and audits
the evidence available in the layer database to diagnose the cause.
"""

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import pymysql


def parse_busco_table(file_path: str, mode: str) -> Dict[str, Dict]:
    """
    Parse a BUSCO full_table.tsv file.

    Args:
        file_path (str): Path to BUSCO full_table.tsv
        mode (str): Either "genome" or "protein" - determines parsing format

    Returns:
        Dict[str, Dict]: Dictionary of BUSCO ID to entry data
    """
    entries = {}
    with open(file_path, "r") as f:  # pylint: disable=unspecified-encoding
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            busco_id = parts[0]
            status = parts[1]  # Complete, Duplicated, Fragmented, Missing

            if mode == "genome":
                # Genome mode format: BUSCO_ID Status Sequence Start End Strand Score Length
                if status == "Missing":
                    entries[busco_id] = {
                        "busco_id": busco_id,
                        "status": status,
                        "sequence": "",
                        "start": 0,
                        "end": 0,
                        "strand": "",
                        "score": 0.0,
                        "length": 0,
                    }
                elif len(parts) >= 8:
                    entries[busco_id] = {
                        "busco_id": busco_id,
                        "status": status,
                        "sequence": parts[2],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[5],
                        "score": float(parts[6]),
                        "length": int(parts[7]),
                    }
            elif mode == "protein":
                # Protein mode format: BUSCO_ID Status Gene/Transcript_ID Score Length [URL] [Description]
                # We only care about status for protein mode
                gene_id = parts[2] if len(parts) > 2 and status != "Missing" else ""
                entries[busco_id] = {
                    "busco_id": busco_id,
                    "status": status,
                    "gene_id": gene_id,  # Store gene/transcript ID from protein mode
                    "sequence": "",  # No genomic coordinates in protein mode
                    "start": 0,
                    "end": 0,
                    "strand": "",
                    "score": 0.0,
                    "length": 0,
                }

    return entries


def identify_problems(
    genome_entries: Dict[str, Dict], protein_entries: Dict[str, Dict]
) -> List[Dict]:
    """
    Identify problematic loci where protein BUSCO is worse than genome BUSCO.
    Goal: protein completeness should be >= genome completeness.

    Args:
        genome_entries (Dict): BUSCO entries from genome mode
        protein_entries (Dict): BUSCO entries from protein mode

    Returns:
        List[Dict]: List of problem loci where protein annotation is failing
    """
    problems = []
    all_busco_ids = set(genome_entries.keys()) | set(protein_entries.keys())

    for busco_id in all_busco_ids:
        genome_entry = genome_entries.get(busco_id)
        protein_entry = protein_entries.get(busco_id)

        if not genome_entry or not protein_entry:
            continue

        genome_status = genome_entry["status"]
        protein_status = protein_entry["status"]

        # Identify problematic cases where protein is WORSE than genome
        problem_type = None

        # Case 1: Genome is Complete, but Protein is Missing
        # This means gene exists but no valid protein
        if genome_status == "Complete" and protein_status == "Missing":
            problem_type = "GenomeComplete→ProteinMissing"
        # Case 2: Genome is Complete, but Protein is Fragmented
        # Gene exists but protein is incomplete
        elif genome_status == "Complete" and protein_status == "Fragmented":
            problem_type = "GenomeComplete→ProteinFragmented"
        # Case 3: Genome is Fragmented, but Protein is Missing
        # Partial gene exists but no protein at all
        elif genome_status == "Fragmented" and protein_status == "Missing":
            problem_type = "GenomeFragmented→ProteinMissing"

        if problem_type:
            # Use GENOME coordinates since genome BUSCO found it
            problems.append(
                {
                    "busco_id": busco_id,
                    "genome_status": genome_status,
                    "protein_status": protein_status,
                    "sequence": genome_entry["sequence"],
                    "start": genome_entry["start"],
                    "end": genome_entry["end"],
                    "strand": genome_entry["strand"],
                    "problem_type": problem_type,
                    "classification": "",
                    "core_genes": [],
                    "layer_genes": [],
                    "notes": "",
                }
            )

    return problems


def get_seq_region_id(sequence: str, host: str, port: int, user: str, password: str, core_db: str) -> Optional[int]:
    """
    Get seq_region_id from seq_region name.

    Args:
        sequence (str): Sequence/contig name
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password (can be empty string for read-only)
        core_db (str): Core database name

    Returns:
        Optional[int]: seq_region_id or None
    """
    try:
        if password:
            conn = pymysql.connect(
                host=host, user=user, password=password, port=port, database=core_db
            )
        else:
            conn = pymysql.connect(
                host=host, user=user, port=port, database=core_db
            )
        cursor = conn.cursor()
        cursor.execute("SELECT seq_region_id FROM seq_region WHERE name = %s", (sequence,))
        result = cursor.fetchone()
        cursor.close()
        conn.close()
        return result[0] if result else None
    except pymysql.Error as err:
        print(f"Error querying seq_region: {err}")
        return None


def get_genes_at_locus(
    seq_region_id: int, start: int, end: int, host: str, port: int, user: str, password: str, core_db: str
) -> List[Dict]:
    """
    Get all genes overlapping a locus.

    Args:
        seq_region_id (int): seq_region_id
        start (int): Start coordinate
        end (int): End coordinate
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password (can be empty string for read-only)
        core_db (str): Core database name

    Returns:
        List[Dict]: List of genes at locus
    """
    try:
        if password:
            conn = pymysql.connect(
                host=host, user=user, password=password, port=port, database=core_db
            )
        else:
            conn = pymysql.connect(
                host=host, user=user, port=port, database=core_db
            )
        cursor = conn.cursor(pymysql.cursors.DictCursor)

        query = """
            SELECT
                g.gene_id,
                g.stable_id,
                sr.name as seq_region_name,
                g.seq_region_start,
                g.seq_region_end,
                g.seq_region_strand,
                g.biotype
            FROM gene g
            JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
            WHERE g.seq_region_id = %s
            AND g.seq_region_start <= %s
            AND g.seq_region_end >= %s
        """
        cursor.execute(query, (seq_region_id, end, start))
        genes = cursor.fetchall()

        # Get transcripts for each gene
        for gene in genes:
            cursor.execute(
                """
                SELECT transcript_id, stable_id, seq_region_start, seq_region_end, biotype
                FROM transcript
                WHERE gene_id = %s
                """,
                (gene["gene_id"],),
            )
            gene["transcripts"] = cursor.fetchall()

        cursor.close()
        conn.close()
        return genes

    except pymysql.Error as err:
        print(f"Error querying genes: {err}")
        return []


def get_layer_genes_at_locus(
    seq_region_id: int, start: int, end: int, host: str, port: int, user: str, password: str, layer_db: str
) -> List[Dict]:
    """
    Get all genes at a locus in the layer database.

    Args:
        seq_region_id (int): seq_region_id
        start (int): Start coordinate
        end (int): End coordinate
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password (can be empty string for read-only)
        layer_db (str): Layer database name

    Returns:
        List[Dict]: List of genes in layer database at locus
    """
    try:
        if password:
            conn = pymysql.connect(
                host=host, user=user, password=password, port=port, database=layer_db
            )
        else:
            conn = pymysql.connect(
                host=host, user=user, port=port, database=layer_db
            )
        cursor = conn.cursor(pymysql.cursors.DictCursor)

        query = """
            SELECT
                g.gene_id,
                g.stable_id,
                sr.name as seq_region_name,
                g.seq_region_start,
                g.seq_region_end,
                g.seq_region_strand,
                g.biotype,
                a.logic_name
            FROM gene g
            JOIN seq_region sr ON g.seq_region_id = sr.seq_region_id
            LEFT JOIN analysis a ON g.analysis_id = a.analysis_id
            WHERE g.seq_region_id = %s
            AND g.seq_region_start <= %s
            AND g.seq_region_end >= %s
        """
        cursor.execute(query, (seq_region_id, end, start))
        genes = cursor.fetchall()

        # Get transcripts for each gene
        for gene in genes:
            cursor.execute(
                """
                SELECT transcript_id, stable_id, seq_region_start, seq_region_end, biotype
                FROM transcript
                WHERE gene_id = %s
                """,
                (gene["gene_id"],),
            )
            gene["transcripts"] = cursor.fetchall()

        cursor.close()
        conn.close()
        return genes

    except pymysql.Error as err:
        print(f"Error querying layer genes: {err}")
        return []


def audit_evidence(
    problems: List[Dict], host: str, port: int, user: str, password: str, core_db: str, layer_db: str
) -> List[Dict]:
    """
    Audit genes at each problematic locus by comparing core vs layer databases.

    Args:
        problems (List[Dict]): List of problem loci
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password
        core_db (str): Core database name
        layer_db (str): Layer database name

    Returns:
        List[Dict]: Updated problem loci with gene comparison data
    """
    total = len(problems)
    for idx, problem in enumerate(problems, 1):
        print(f"Auditing {idx}/{total}: {problem['busco_id']}", end="\r")

        # Get seq_region_id
        seq_region_id = get_seq_region_id(
            problem["sequence"], host, port, user, password, core_db
        )
        if not seq_region_id:
            problem["classification"] = "UNKNOWN"
            problem["notes"] = f"Could not map sequence '{problem['sequence']}' to seq_region_id"
            continue

        # Get genes at locus in core database
        problem["core_genes"] = get_genes_at_locus(
            seq_region_id, problem["start"], problem["end"], host, port, user, password, core_db
        )

        # Get genes at locus in layer database
        problem["layer_genes"] = get_layer_genes_at_locus(
            seq_region_id, problem["start"], problem["end"], host, port, user, password, layer_db
        )

        # Classify the problem
        classify_problem(problem)

    print(f"\nCompleted evidence audit for {total} loci")
    return problems


def classify_problem(problem: Dict):
    """
    Classify why protein annotation is failing at this locus.
    Focus: Core biotypes vs Layer evidence sources.

    Args:
        problem (Dict): Problem locus data (modified in place)
    """
    core_genes = problem["core_genes"]
    layer_genes = problem["layer_genes"]

    # Get biotype counts for core (these are real biotypes)
    core_biotype_counts = {}
    for g in core_genes:
        biotype = g["biotype"] or "unknown"
        core_biotype_counts[biotype] = core_biotype_counts.get(biotype, 0) + 1

    # Get biotype counts for layer (what's stored in layer biotype field - often logic_names)
    layer_biotype_counts = {}
    layer_logic_name_counts = {}
    for g in layer_genes:
        biotype = g["biotype"] or "unknown"
        logic_name = g.get("logic_name", "unknown")
        layer_biotype_counts[biotype] = layer_biotype_counts.get(biotype, 0) + 1
        layer_logic_name_counts[logic_name] = layer_logic_name_counts.get(logic_name, 0) + 1

    # Store detailed breakdowns
    problem["core_biotype_counts"] = core_biotype_counts
    problem["layer_biotype_counts"] = layer_biotype_counts
    problem["layer_logic_name_counts"] = layer_logic_name_counts
    problem["layer_evidence_count"] = len(layer_genes)  # Total evidence alignments

    # Count protein_coding specifically in core
    core_protein_coding = core_biotype_counts.get("protein_coding", 0)

    # Classification logic
    if not core_genes:
        if layer_genes:
            problem["classification"] = "NO_GENE_BUILT"
            problem["notes"] = (
                f"No genes in core. Layer has {len(layer_genes)} evidence alignments from sources: "
                f"{dict(layer_logic_name_counts)}. Evidence not merged/selected into core gene."
            )
        else:
            problem["classification"] = "NO_GENES_FOUND"
            problem["notes"] = "Genome BUSCO found locus but no genes in core or layer"
        return

    # Check if core only has non-coding
    if core_protein_coding == 0:
        problem["classification"] = "NON_CODING_BIOTYPE"
        problem["notes"] = (
            f"Core has {len(core_genes)} gene(s) with non-coding biotypes: {dict(core_biotype_counts)}. "
            f"Layer has {len(layer_genes)} evidence alignments from: {dict(layer_logic_name_counts)}. "
            f"Gene classified as non-coding (pseudogene/lncRNA), so no valid protein for BUSCO."
        )
        return

    # Core has protein_coding but BUSCO protein mode fails
    genes_without_transcripts = [g for g in core_genes if not g.get("transcripts")]
    if genes_without_transcripts:
        problem["classification"] = "NO_TRANSCRIPTS"
        problem["notes"] = (
            f"Core has {core_protein_coding} protein_coding gene(s) but {len(genes_without_transcripts)} lack transcripts. "
            f"Core biotypes: {dict(core_biotype_counts)}. "
            f"Layer evidence: {len(layer_genes)} alignments from {dict(layer_logic_name_counts)}"
        )
        return

    # Core has protein_coding genes but they don't satisfy BUSCO
    problem["classification"] = "PROTEIN_CODING_BUT_NOT_BUSCO_SATISFYING"

    # Get evidence source summary
    evidence_summary = ", ".join(f"{k}({v})" for k, v in sorted(layer_logic_name_counts.items(), key=lambda x: x[1], reverse=True)[:5])
    if len(layer_logic_name_counts) > 5:
        evidence_summary += f" and {len(layer_logic_name_counts) - 5} more"

    problem["notes"] = (
        f"Core has {core_protein_coding} protein_coding gene(s) but BUSCO protein mode fails. "
        f"Possible issues: incomplete CDS, frameshifts, premature stops, or poor canonical selection. "
        f"Core biotypes: {dict(core_biotype_counts)}. "
        f"Layer: {len(layer_genes)} evidence alignments from {evidence_summary}"
    )


def generate_report(problems: List[Dict], output_tsv: str, output_json: Optional[str] = None):
    """
    Generate TSV report and optional JSON output.

    Args:
        problems (List[Dict]): List of problem loci with evidence
        output_tsv (str): Output TSV file path
        output_json (Optional[str]): Optional JSON output file path
    """
    # Generate TSV report
    with open(output_tsv, "w", newline="") as f:  # pylint: disable=unspecified-encoding
        writer = csv.writer(f, delimiter="\t")

        # Header
        writer.writerow(
            [
                "BUSCO_ID",
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "Genome_BUSCO",
                "Protein_BUSCO",
                "Issue",
                "Core_Genes",
                "Core_Protein_Coding",
                "Core_Biotypes",
                "Core_Gene_IDs",
                "Layer_Evidence_Total",
                "Layer_Top_Biotypes",
                "Layer_Top_Analyses",
                "Summary",
            ]
        )

        # Data rows
        for problem in problems:
            core_gene_ids = ",".join(g["stable_id"] or "None" for g in problem["core_genes"]) if problem["core_genes"] else "-"

            # Format biotype breakdowns
            core_biotype_counts = problem.get("core_biotype_counts", {})
            layer_biotype_counts = problem.get("layer_biotype_counts", {})
            layer_logic_name_counts = problem.get("layer_logic_name_counts", {})
            layer_evidence_count = problem.get("layer_evidence_count", 0)

            # Core biotypes - just list unique ones
            core_biotypes_str = ",".join(sorted(core_biotype_counts.keys())) if core_biotype_counts else "-"

            # Layer biotypes - show top 3 (what's in layer biotype field)
            top_layer_biotypes = sorted(layer_biotype_counts.items(), key=lambda x: x[1], reverse=True)[:3]
            layer_biotype_str = ", ".join(f"{k}({v})" for k, v in top_layer_biotypes) if top_layer_biotypes else "-"

            # Layer analyses - show top 3 logic names
            top_analyses = sorted(layer_logic_name_counts.items(), key=lambda x: x[1], reverse=True)[:3]
            layer_analysis_str = ", ".join(f"{k}({v})" for k, v in top_analyses) if top_analyses else "-"

            core_protein_coding = core_biotype_counts.get("protein_coding", 0)

            # Simplified summary based on classification
            classification = problem["classification"]
            if classification == "PROTEIN_CODING_BUT_NOT_BUSCO_SATISFYING":
                summary = f"{core_protein_coding} protein_coding gene(s), likely incomplete CDS or poor canonical"
            elif classification == "NON_CODING_BIOTYPE":
                summary = f"Gene classified as {core_biotypes_str}, no valid protein"
            elif classification == "NO_GENE_BUILT":
                summary = f"{layer_evidence_count} evidence alignments not merged into core gene"
            elif classification == "NO_GENES_FOUND":
                summary = "BUSCO found locus but no genes/evidence in databases"
            elif classification == "NO_TRANSCRIPTS":
                summary = f"{core_protein_coding} protein_coding gene(s) but missing transcripts"
            else:
                summary = problem["notes"][:100]

            writer.writerow(
                [
                    problem["busco_id"],
                    problem["sequence"],
                    problem["start"],
                    problem["end"],
                    problem["strand"],
                    problem["genome_status"],
                    problem["protein_status"],
                    classification,
                    len(problem["core_genes"]),
                    core_protein_coding,
                    core_biotypes_str,
                    core_gene_ids,
                    layer_evidence_count,
                    layer_biotype_str,
                    layer_analysis_str,
                    summary,
                ]
            )

    print(f"TSV report written to: {output_tsv}")

    # Generate summary statistics
    print_summary(problems)

    # Generate JSON output if requested
    if output_json:
        generate_json_report(problems, output_json)


def print_summary(problems: List[Dict]):
    """
    Print summary statistics.

    Args:
        problems (List[Dict]): List of problem loci
    """
    print("\n" + "=" * 80)
    print("BUSCO EVIDENCE AUDIT SUMMARY")
    print("=" * 80)

    # Classification breakdown
    classification_counts = defaultdict(int)
    problem_type_counts = defaultdict(int)

    for problem in problems:
        classification_counts[problem["classification"]] += 1
        problem_type_counts[problem["problem_type"]] += 1

    print("\nProblem Types:")
    for ptype, count in sorted(problem_type_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {ptype}: {count}")

    print("\nClassifications:")
    for classification, count in sorted(classification_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {classification}: {count}")

    # Gene statistics
    total_core_genes = sum(len(p["core_genes"]) for p in problems)
    total_layer_evidence = sum(p.get("layer_evidence_count", 0) for p in problems)

    # Protein coding statistics
    total_core_protein_coding = sum(p.get("core_biotype_counts", {}).get("protein_coding", 0) for p in problems)

    # Aggregate all biotype counts and evidence sources
    core_biotype_totals = defaultdict(int)
    layer_biotype_totals = defaultdict(int)
    layer_logic_name_totals = defaultdict(int)

    for problem in problems:
        for biotype, count in problem.get("core_biotype_counts", {}).items():
            core_biotype_totals[biotype] += count
        for biotype, count in problem.get("layer_biotype_counts", {}).items():
            layer_biotype_totals[biotype] += count
        for logic_name, count in problem.get("layer_logic_name_counts", {}).items():
            layer_logic_name_totals[logic_name] += count

    print(f"\nCore Gene Statistics (at problem loci):")
    print(f"  Total genes in core: {total_core_genes}")
    print(f"  Protein_coding genes: {total_core_protein_coding} ({100*total_core_protein_coding/total_core_genes if total_core_genes > 0 else 0:.1f}%)")
    print(f"  Non-coding genes: {total_core_genes - total_core_protein_coding} ({100*(total_core_genes - total_core_protein_coding)/total_core_genes if total_core_genes > 0 else 0:.1f}%)")

    print(f"\nLayer Evidence Statistics (at problem loci):")
    print(f"  Total evidence alignments: {total_layer_evidence}")
    print(f"  Evidence alignments per core gene: {total_layer_evidence/total_core_genes if total_core_genes > 0 else 0:.1f}")
    print(f"  Unique analysis sources: {len(layer_logic_name_totals)}")

    print(f"\nCore Biotype Breakdown:")
    for biotype, count in sorted(core_biotype_totals.items(), key=lambda x: x[1], reverse=True):
        pct = 100 * count / sum(core_biotype_totals.values()) if sum(core_biotype_totals.values()) > 0 else 0
        print(f"  {biotype}: {count} ({pct:.1f}%)")

    print(f"\nLayer Biotype Breakdown (top 20):")
    for idx, (biotype, count) in enumerate(sorted(layer_biotype_totals.items(), key=lambda x: x[1], reverse=True)[:20], 1):
        pct = 100 * count / sum(layer_biotype_totals.values()) if sum(layer_biotype_totals.values()) > 0 else 0
        print(f"  {idx}. {biotype}: {count} ({pct:.1f}%)")

    print(f"\nLayer Analysis Source Breakdown (top 20):")
    for idx, (logic_name, count) in enumerate(sorted(layer_logic_name_totals.items(), key=lambda x: x[1], reverse=True)[:20], 1):
        pct = 100 * count / sum(layer_logic_name_totals.values()) if sum(layer_logic_name_totals.values()) > 0 else 0
        print(f"  {idx}. {logic_name}: {count} ({pct:.1f}%)")

    print("=" * 80)


def generate_json_report(problems: List[Dict], output_json: str):
    """
    Generate detailed JSON report.

    Args:
        problems (List[Dict]): List of problem loci
        output_json (str): Output JSON file path
    """
    report = {
        "summary": {
            "total_problems": len(problems),
            "classifications": {},
            "problem_types": {},
        },
        "problems": [],
    }

    # Summary statistics
    for problem in problems:
        report["summary"]["classifications"][problem["classification"]] = (
            report["summary"]["classifications"].get(problem["classification"], 0) + 1
        )
        report["summary"]["problem_types"][problem["problem_type"]] = (
            report["summary"]["problem_types"].get(problem["problem_type"], 0) + 1
        )

    # Detailed problem data
    for problem in problems:
        problem_data = {
            "busco_id": problem["busco_id"],
            "problem_type": problem["problem_type"],
            "classification": problem["classification"],
            "coordinates": {
                "sequence": problem["sequence"],
                "start": problem["start"],
                "end": problem["end"],
                "strand": problem["strand"],
            },
            "statuses": {
                "genome": problem["genome_status"],
                "protein": problem["protein_status"],
            },
            "core_genes": [
                {
                    "stable_id": g["stable_id"] or "None",
                    "biotype": g["biotype"] or "unknown",
                    "coordinates": {
                        "start": g["seq_region_start"],
                        "end": g["seq_region_end"],
                        "strand": g["seq_region_strand"],
                    },
                }
                for g in problem["core_genes"]
            ],
            "layer_genes": [
                {
                    "stable_id": g["stable_id"] or "None",
                    "biotype": g["biotype"] or "unknown",
                    "logic_name": g.get("logic_name", "unknown"),
                    "coordinates": {
                        "start": g["seq_region_start"],
                        "end": g["seq_region_end"],
                        "strand": g["seq_region_strand"],
                    },
                }
                for g in problem["layer_genes"]
            ],
            "notes": problem["notes"],
        }
        report["problems"].append(problem_data)

    with open(output_json, "w") as f:  # pylint: disable=unspecified-encoding
        json.dump(report, f, indent=2)

    print(f"JSON report written to: {output_json}")


def parse_db_name(db_name: str) -> Tuple[str, str, str]:
    """
    Parse database name to extract species, version, and assembly.

    Args:
        db_name (str): Database name (e.g., homo_sapiens_core_113_38 or jackt_gca049306965v1_core_114)

    Returns:
        Tuple[str, str, str]: (species, version, assembly)

    Raises:
        ValueError: If database name format is invalid
    """
    # Try format with assembly: species_core_version_assembly
    pattern_with_assembly = r"^(.+)_core_(\d+)_(\d+)$"
    match = re.match(pattern_with_assembly, db_name)
    if match:
        return match.groups()

    # Try format without assembly: species_core_version
    pattern_without_assembly = r"^(.+)_core_(\d+)$"
    match = re.match(pattern_without_assembly, db_name)
    if match:
        species, version = match.groups()
        return (species, version, "")

    raise ValueError(
        f"Invalid database name format: {db_name}. Expected format: species_core_version or species_core_version_assembly"
    )


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="BUSCO Evidence Audit - Compare genome and protein mode results and audit evidence"
    )

    # Input files
    parser.add_argument(
        "--genome_busco",
        required=True,
        type=str,
        help="Path to BUSCO genome mode full_table.tsv",
    )
    parser.add_argument(
        "--protein_busco",
        required=True,
        type=str,
        help="Path to BUSCO protein mode full_table.tsv",
    )

    # Database parameters
    parser.add_argument(
        "--core_db",
        required=True,
        type=str,
        help="Core database name (e.g., homo_sapiens_core_113_38)",
    )
    parser.add_argument(
        "--layer_db",
        type=str,
        help="Layer database name (default: auto-derived from core_db)",
    )
    parser.add_argument("--host", required=True, type=str, help="Database host")
    parser.add_argument("--port", required=True, type=int, help="Database port (e.g., 3306)")
    parser.add_argument("--user", required=True, type=str, help="Database user")
    parser.add_argument("--password", type=str, default="", help="Database password (optional for read-only users)")

    # Output
    parser.add_argument("--output_tsv", required=True, type=str, help="Output TSV report file")
    parser.add_argument("--output_json", type=str, help="Optional detailed JSON output file")

    args = parser.parse_args()

    # Parse database names
    try:
        species, version, assembly = parse_db_name(args.core_db)
        core_db = args.core_db
        if args.layer_db:
            layer_db = args.layer_db
        elif assembly:
            layer_db = f"{species}_layer_{version}_{assembly}"
        else:
            layer_db = f"{species}_layer_{version}"
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    print("=" * 80)
    print("BUSCO EVIDENCE AUDIT")
    print("=" * 80)
    print(f"Core database: {core_db}")
    print(f"Layer database: {layer_db}")
    print(f"Database host: {args.host}:{args.port}")
    print("=" * 80)

    # Parse BUSCO tables
    print("\nParsing BUSCO results...")
    genome_entries = parse_busco_table(args.genome_busco, "genome")
    print(f"Parsed {len(genome_entries)} BUSCO entries from genome mode")
    protein_entries = parse_busco_table(args.protein_busco, "protein")
    print(f"Parsed {len(protein_entries)} BUSCO entries from protein mode")

    # Identify problems
    print("\nIdentifying problematic loci...")
    problems = identify_problems(genome_entries, protein_entries)
    print(f"Identified {len(problems)} problematic loci")

    # Audit evidence
    print("\nAuditing evidence at problematic loci...")
    problems = audit_evidence(problems, args.host, args.port, args.user, args.password, core_db, layer_db)

    # Generate reports
    print("\nGenerating reports...")
    generate_report(problems, args.output_tsv, args.output_json)

    print("\nAudit complete!")


if __name__ == "__main__":
    main()

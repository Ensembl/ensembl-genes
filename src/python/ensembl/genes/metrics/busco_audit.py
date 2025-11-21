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
from typing import Dict, List, Optional, Set, Tuple

import pymysql


def parse_busco_table(file_path: str) -> Dict[str, Dict]:
    """
    Parse a BUSCO full_table.tsv file.

    Args:
        file_path (str): Path to BUSCO full_table.tsv

    Returns:
        Dict[str, Dict]: Dictionary of BUSCO ID to entry data
    """
    entries = {}
    with open(file_path, "r") as f:  # pylint: disable=unspecified-encoding
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue

            busco_id = parts[0]
            status = parts[1]  # Complete, Duplicated, Fragmented, Missing

            # Handle Missing entries (no coordinates)
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

    return entries


def identify_problems(
    genome_entries: Dict[str, Dict], protein_entries: Dict[str, Dict]
) -> List[Dict]:
    """
    Identify problematic loci where genome and protein results differ.

    Args:
        genome_entries (Dict): BUSCO entries from genome mode
        protein_entries (Dict): BUSCO entries from protein mode

    Returns:
        List[Dict]: List of problem loci
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

        # Identify problematic cases
        problem_type = None

        # Case 1: Protein is Complete, but Genome is Missing
        if protein_status == "Complete" and genome_status == "Missing":
            problem_type = "Complete→Missing"
        # Case 2: Protein is Complete, but Genome is Fragmented
        elif protein_status == "Complete" and genome_status == "Fragmented":
            problem_type = "Complete→Fragmented"
        # Case 3: Protein is Fragmented, but Genome is Missing
        elif protein_status == "Fragmented" and genome_status == "Missing":
            problem_type = "Fragmented→Missing"

        if problem_type:
            # Use protein coordinates for investigation
            problems.append(
                {
                    "busco_id": busco_id,
                    "genome_status": genome_status,
                    "protein_status": protein_status,
                    "sequence": protein_entry["sequence"],
                    "start": protein_entry["start"],
                    "end": protein_entry["end"],
                    "strand": protein_entry["strand"],
                    "problem_type": problem_type,
                    "classification": "",
                    "genes_at_locus": [],
                    "protein_evidence": [],
                    "dna_evidence": [],
                    "supporting_evidence_used": set(),
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
        password (str): Database password
        core_db (str): Core database name

    Returns:
        Optional[int]: seq_region_id or None
    """
    try:
        conn = pymysql.connect(
            host=host, user=user, password=password, port=port, database=core_db
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
        password (str): Database password
        core_db (str): Core database name

    Returns:
        List[Dict]: List of genes at locus
    """
    try:
        conn = pymysql.connect(
            host=host, user=user, password=password, port=port, database=core_db
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


def get_protein_evidence_at_locus(
    seq_region_id: int, start: int, end: int, host: str, port: int, user: str, password: str, layer_db: str
) -> List[Dict]:
    """
    Get protein alignment evidence from layer database.

    Args:
        seq_region_id (int): seq_region_id
        start (int): Start coordinate
        end (int): End coordinate
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password
        layer_db (str): Layer database name

    Returns:
        List[Dict]: List of protein evidence alignments
    """
    try:
        conn = pymysql.connect(
            host=host, user=user, password=password, port=port, database=layer_db
        )
        cursor = conn.cursor(pymysql.cursors.DictCursor)

        query = """
            SELECT
                paf.hit_name,
                sr.name as seq_region_name,
                paf.seq_region_start,
                paf.seq_region_end,
                paf.seq_region_strand,
                paf.score,
                paf.perc_ident,
                a.logic_name,
                ed.db_name as external_db_name
            FROM protein_align_feature paf
            JOIN seq_region sr ON paf.seq_region_id = sr.seq_region_id
            JOIN analysis a ON paf.analysis_id = a.analysis_id
            LEFT JOIN external_db ed ON paf.external_db_id = ed.external_db_id
            WHERE paf.seq_region_id = %s
            AND paf.seq_region_start <= %s
            AND paf.seq_region_end >= %s
            ORDER BY paf.score DESC
        """
        cursor.execute(query, (seq_region_id, end, start))
        results = cursor.fetchall()
        cursor.close()
        conn.close()
        return results

    except pymysql.Error as err:
        print(f"Error querying protein evidence: {err}")
        return []


def get_dna_evidence_at_locus(
    seq_region_id: int, start: int, end: int, host: str, port: int, user: str, password: str, layer_db: str
) -> List[Dict]:
    """
    Get DNA alignment evidence from layer database.

    Args:
        seq_region_id (int): seq_region_id
        start (int): Start coordinate
        end (int): End coordinate
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password
        layer_db (str): Layer database name

    Returns:
        List[Dict]: List of DNA evidence alignments
    """
    try:
        conn = pymysql.connect(
            host=host, user=user, password=password, port=port, database=layer_db
        )
        cursor = conn.cursor(pymysql.cursors.DictCursor)

        query = """
            SELECT
                daf.hit_name,
                sr.name as seq_region_name,
                daf.seq_region_start,
                daf.seq_region_end,
                daf.seq_region_strand,
                daf.score,
                daf.perc_ident,
                a.logic_name,
                ed.db_name as external_db_name
            FROM dna_align_feature daf
            JOIN seq_region sr ON daf.seq_region_id = sr.seq_region_id
            JOIN analysis a ON daf.analysis_id = a.analysis_id
            LEFT JOIN external_db ed ON daf.external_db_id = ed.external_db_id
            WHERE daf.seq_region_id = %s
            AND daf.seq_region_start <= %s
            AND daf.seq_region_end >= %s
            ORDER BY daf.score DESC
        """
        cursor.execute(query, (seq_region_id, end, start))
        results = cursor.fetchall()
        cursor.close()
        conn.close()
        return results

    except pymysql.Error as err:
        print(f"Error querying DNA evidence: {err}")
        return []


def get_supporting_evidence_for_genes(
    gene_ids: List[int], host: str, port: int, user: str, password: str, core_db: str
) -> Set[str]:
    """
    Get supporting evidence used in final gene models.

    Args:
        gene_ids (List[int]): List of gene_ids
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password
        core_db (str): Core database name

    Returns:
        Set[str]: Set of feature IDs used as supporting evidence
    """
    if not gene_ids:
        return set()

    try:
        conn = pymysql.connect(
            host=host, user=user, password=password, port=port, database=core_db
        )
        cursor = conn.cursor()

        placeholders = ",".join(["%s"] * len(gene_ids))
        query = f"""
            SELECT DISTINCT tse.feature_id
            FROM transcript_supporting_evidence tse
            JOIN transcript t ON tse.transcript_id = t.transcript_id
            WHERE t.gene_id IN ({placeholders})
        """
        cursor.execute(query, tuple(gene_ids))
        results = cursor.fetchall()
        cursor.close()
        conn.close()
        return {str(row[0]) for row in results}

    except pymysql.Error as err:
        print(f"Error querying supporting evidence: {err}")
        return set()


def audit_evidence(
    problems: List[Dict], host: str, port: int, user: str, password: str, core_db: str, layer_db: str
) -> List[Dict]:
    """
    Audit evidence for each problematic locus.

    Args:
        problems (List[Dict]): List of problem loci
        host (str): Database host
        port (int): Database port
        user (str): Database user
        password (str): Database password
        core_db (str): Core database name
        layer_db (str): Layer database name

    Returns:
        List[Dict]: Updated problem loci with evidence data
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

        # Get genes at locus
        problem["genes_at_locus"] = get_genes_at_locus(
            seq_region_id, problem["start"], problem["end"], host, port, user, password, core_db
        )

        # Get evidence from layer database
        problem["protein_evidence"] = get_protein_evidence_at_locus(
            seq_region_id, problem["start"], problem["end"], host, port, user, password, layer_db
        )
        problem["dna_evidence"] = get_dna_evidence_at_locus(
            seq_region_id, problem["start"], problem["end"], host, port, user, password, layer_db
        )

        # Get supporting evidence used in final models
        if problem["genes_at_locus"]:
            gene_ids = [g["gene_id"] for g in problem["genes_at_locus"]]
            problem["supporting_evidence_used"] = get_supporting_evidence_for_genes(
                gene_ids, host, port, user, password, core_db
            )

        # Classify the problem
        classify_problem(problem)

    print(f"\nCompleted evidence audit for {total} loci")
    return problems


def classify_problem(problem: Dict):
    """
    Classify the problem based on evidence.

    Args:
        problem (Dict): Problem locus data (modified in place)
    """
    # Case 1: No evidence in layer database
    if not problem["protein_evidence"] and not problem["dna_evidence"]:
        problem["classification"] = "NO_EVIDENCE"
        problem["notes"] = "No protein or DNA evidence found in layer database at this locus"
        return

    # Case 2: Evidence exists but no genes built
    if not problem["genes_at_locus"]:
        if problem["protein_evidence"]:
            problem["classification"] = "EVIDENCE_NOT_USED"
            problem["notes"] = f"Found {len(problem['protein_evidence'])} protein evidence alignments but no genes built"
        else:
            problem["classification"] = "EVIDENCE_NOT_USED"
            problem["notes"] = f"Found {len(problem['dna_evidence'])} DNA evidence alignments but no genes built"
        return

    # Case 3: Genes exist but wrong biotype/filtering
    wrong_biotypes = [
        g for g in problem["genes_at_locus"] if g["biotype"] in ("pseudogene", "artifact", "TEC")
    ]
    if wrong_biotypes:
        problem["classification"] = "WRONG_BIOTYPE"
        problem["notes"] = f"Genes built but classified as: {', '.join(g['biotype'] for g in wrong_biotypes)}"
        return

    # Case 4: Genes exist but evidence not used
    if problem["protein_evidence"]:
        evidence_ids = {e["hit_name"] for e in problem["protein_evidence"][:10]}
        used_evidence = evidence_ids & problem["supporting_evidence_used"]
        if not used_evidence:
            problem["classification"] = "EVIDENCE_IGNORED"
            problem["notes"] = (
                f"Found {len(problem['protein_evidence'])} protein alignments but none used in final models"
            )
            return

    # Case 5: Evidence was used but model is incomplete
    problem["classification"] = "MODEL_INCOMPLETE"
    problem["notes"] = (
        f"Evidence used ({len(problem['supporting_evidence_used'])} features) but resulting model may be fragmented"
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
                "Problem_Type",
                "Classification",
                "Sequence",
                "Start",
                "End",
                "Strand",
                "Genome_Status",
                "Protein_Status",
                "Genes_Count",
                "Gene_Stable_IDs",
                "Gene_Biotypes",
                "Protein_Evidence_Count",
                "DNA_Evidence_Count",
                "Evidence_Used_Count",
                "Top_Evidence_Logic_Names",
                "Notes",
            ]
        )

        # Data rows
        for problem in problems:
            gene_ids = (
                ",".join(g["stable_id"] for g in problem["genes_at_locus"])
                if problem["genes_at_locus"]
                else ""
            )
            gene_biotypes = (
                ",".join(set(g["biotype"] for g in problem["genes_at_locus"]))
                if problem["genes_at_locus"]
                else ""
            )
            top_logic_names = (
                ",".join(
                    set(
                        e["logic_name"]
                        for e in (problem["protein_evidence"][:5] + problem["dna_evidence"][:5])
                    )
                )
                if (problem["protein_evidence"] or problem["dna_evidence"])
                else ""
            )

            writer.writerow(
                [
                    problem["busco_id"],
                    problem["problem_type"],
                    problem["classification"],
                    problem["sequence"],
                    problem["start"],
                    problem["end"],
                    problem["strand"],
                    problem["genome_status"],
                    problem["protein_status"],
                    len(problem["genes_at_locus"]),
                    gene_ids,
                    gene_biotypes,
                    len(problem["protein_evidence"]),
                    len(problem["dna_evidence"]),
                    len(problem["supporting_evidence_used"]),
                    top_logic_names,
                    problem["notes"],
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

    # Evidence statistics
    total_protein_evidence = sum(len(p["protein_evidence"]) for p in problems)
    total_dna_evidence = sum(len(p["dna_evidence"]) for p in problems)
    total_genes = sum(len(p["genes_at_locus"]) for p in problems)

    print(f"\nEvidence Statistics:")
    print(f"  Total protein evidence alignments: {total_protein_evidence}")
    print(f"  Total DNA evidence alignments: {total_dna_evidence}")
    print(f"  Total genes at problem loci: {total_genes}")

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

    # Detailed problem data (convert sets to lists for JSON serialization)
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
            "genes": [
                {
                    "stable_id": g["stable_id"],
                    "biotype": g["biotype"],
                    "coordinates": {
                        "start": g["seq_region_start"],
                        "end": g["seq_region_end"],
                        "strand": g["seq_region_strand"],
                    },
                }
                for g in problem["genes_at_locus"]
            ],
            "evidence_counts": {
                "protein": len(problem["protein_evidence"]),
                "dna": len(problem["dna_evidence"]),
                "used": len(problem["supporting_evidence_used"]),
            },
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
        db_name (str): Database name (e.g., homo_sapiens_core_113_38)

    Returns:
        Tuple[str, str, str]: (species, version, assembly)

    Raises:
        ValueError: If database name format is invalid
    """
    pattern = r"^(.+)_core_(\d+)_(\d+)$"
    match = re.match(pattern, db_name)
    if not match:
        raise ValueError(
            f"Invalid database name format: {db_name}. Expected format: species_core_version_assembly"
        )
    return match.groups()


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
    parser.add_argument("--password", required=True, type=str, help="Database password")

    # Output
    parser.add_argument("--output_tsv", required=True, type=str, help="Output TSV report file")
    parser.add_argument("--output_json", type=str, help="Optional detailed JSON output file")

    args = parser.parse_args()

    # Parse database names
    try:
        species, version, assembly = parse_db_name(args.core_db)
        core_db = args.core_db
        layer_db = args.layer_db or f"{species}_layer_{version}_{assembly}"
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
    genome_entries = parse_busco_table(args.genome_busco)
    print(f"Parsed {len(genome_entries)} BUSCO entries from genome mode")
    protein_entries = parse_busco_table(args.protein_busco)
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

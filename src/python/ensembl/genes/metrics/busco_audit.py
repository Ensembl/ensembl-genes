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
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pymysql


def parse_busco_id(busco_id: str) -> Dict[str, str]:
    """
    Parse BUSCO ID to extract OrthoDB group and taxonomic info.

    Args:
        busco_id (str): BUSCO ID like "135504at7898"

    Returns:
        Dict with orthodb_group and taxid
    """
    match = re.match(r"(\d+)at(\d+)", busco_id)
    if match:
        return {
            "orthodb_group": match.group(1),
            "taxid": match.group(2),
        }
    return {"orthodb_group": busco_id, "taxid": "unknown"}


def get_busco_metadata(busco_dataset_dir: str, busco_ids: List[str]) -> Dict[str, Dict]:
    """
    Extract BUSCO metadata from the lineage dataset directory.

    Args:
        busco_dataset_dir (str): Path to BUSCO lineage directory (e.g., actinopterygii_odb10)
        busco_ids (List[str]): List of BUSCO IDs to get metadata for

    Returns:
        Dict mapping BUSCO ID to metadata (name, length, etc.)
    """
    metadata = {}

    if not busco_dataset_dir:
        print("Warning: No BUSCO dataset directory provided")
        return metadata

    busco_path = Path(busco_dataset_dir)
    if not busco_path.exists():
        print(f"Warning: BUSCO dataset directory does not exist: {busco_dataset_dir}")
        return metadata

    # Read lengths_cutoff file which has expected protein lengths
    lengths_file = busco_path / "lengths_cutoff"

    print(f"Looking for lengths_cutoff at: {lengths_file}")
    print(f"File exists: {lengths_file.exists()}")

    busco_lengths = {}

    # If direct path doesn't exist, try common BUSCO dataset structures
    if not lengths_file.exists():
        print(f"Checking for common BUSCO dataset structures...")
        possible_locations = [
            busco_path / "lineages" / "actinopterygii_odb10" / "lengths_cutoff",
            busco_path / "actinopterygii_odb10" / "lengths_cutoff",
        ]
        for loc in possible_locations:
            print(f"  Trying: {loc}")
            if loc.exists():
                print(f"  Found it!")
                lengths_file = loc
                break

    if lengths_file.exists():
        try:
            with open(lengths_file, "r") as f:  # pylint: disable=unspecified-encoding
                line_count = 0
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        busco_id = parts[0]
                        busco_lengths[busco_id] = {
                            "length_sd": float(parts[1]),
                            "length": float(parts[2]),
                        }
                        line_count += 1
                print(f"Loaded {line_count} BUSCO length entries from lengths_cutoff")
        except Exception as e:  # pylint: disable=broad-except
            print(f"Warning: Could not parse lengths_cutoff: {e}")
    else:
        print(f"Warning: Could not find lengths_cutoff file in any expected location")

    for busco_id in busco_ids:
        parsed = parse_busco_id(busco_id)
        metadata[busco_id] = {
            "orthodb_group": parsed["orthodb_group"],
            "taxid": parsed["taxid"],
            "length": busco_lengths.get(busco_id, {}).get("length"),
            "length_sd": busco_lengths.get(busco_id, {}).get("length_sd"),
        }

    return metadata


def parse_hmmer_output(hmmer_file: str, core_gene_ids: List[str]) -> Dict:
    """
    Parse HMMER output file to understand why BUSCO failed.
    Also identifies alternative transcript isoforms that would satisfy BUSCO.

    The HMMER output contains results for ALL transcripts that were tested by BUSCO,
    not just the canonical. This function identifies which transcripts would satisfy
    BUSCO criteria, allowing you to see if a better isoform exists.

    Args:
        hmmer_file (str): Path to HMMER .out file
        core_gene_ids (List[str]): List of gene stable IDs to look for

    Returns:
        Dict: HMMER results with best hit info and alternative isoforms
    """
    if not Path(hmmer_file).exists():
        return {"status": "file_not_found"}

    try:
        with open(hmmer_file, "r") as f:  # pylint: disable=unspecified-encoding
            lines = f.readlines()

        # Find hits for our genes (matches transcript IDs, but we have gene IDs)
        hits = []
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()
            if len(parts) < 19:
                continue

            transcript_id = parts[0]
            # Extract gene ID from transcript (e.g., ENSDARG00160000965.1 from ENSDART00160002563.1)
            gene_id_match = re.match(r"(ENSDARG\d+)", transcript_id)
            if not gene_id_match:
                continue

            gene_id = gene_id_match.group(1)
            if gene_id not in core_gene_ids:
                continue

            hits.append({
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "tlen": int(parts[2]),  # Target (protein) length
                "qlen": int(parts[5]),  # Query (BUSCO HMM) length
                "e_value": float(parts[6]),
                "score": float(parts[7]),
                "bias": float(parts[8]),
                "domain_num": int(parts[9]),
                "domain_total": int(parts[10]),
                "c_evalue": float(parts[11]) if parts[11] != "-" else None,
                "hmm_from": int(parts[15]),
                "hmm_to": int(parts[16]),
                "ali_from": int(parts[17]),
                "ali_to": int(parts[18]),
            })

        if not hits:
            return {"status": "no_hits_for_genes"}

        # Get best hit (lowest e-value)
        best_hit = min(hits, key=lambda x: x["e_value"])

        # Calculate coverage for best hit
        hmm_coverage = (best_hit["hmm_to"] - best_hit["hmm_from"] + 1) / best_hit["qlen"] * 100
        protein_coverage = (best_hit["ali_to"] - best_hit["ali_from"] + 1) / best_hit["tlen"] * 100

        # Diagnose issues with best hit
        issues = []
        if best_hit["e_value"] > 1e-10:
            issues.append("poor_evalue")
        if hmm_coverage < 80:
            issues.append("low_hmm_coverage")
        if best_hit["tlen"] < best_hit["qlen"] * 0.5:
            issues.append("protein_too_short")
        if best_hit["domain_total"] > 1:
            issues.append("fragmented_match")
        if protein_coverage < 50:
            issues.append("low_protein_coverage")

        # Find isoforms that would satisfy BUSCO criteria
        # Criteria: e-value <= 1e-10, HMM coverage >= 80%, single domain, protein coverage >= 50%
        # These are transcripts that exist in the HMMER results (meaning they were tested)
        # and would pass BUSCO if they were selected as canonical
        good_isoforms = []
        for hit in hits:
            hit_hmm_cov = (hit["hmm_to"] - hit["hmm_from"] + 1) / hit["qlen"] * 100
            hit_prot_cov = (hit["ali_to"] - hit["ali_from"] + 1) / hit["tlen"] * 100

            if (hit["e_value"] <= 1e-10 and
                hit_hmm_cov >= 80 and
                hit["domain_total"] == 1 and
                hit_prot_cov >= 50):

                # Mark if this is the best hit (likely the current canonical or best transcript)
                is_best = (hit["transcript_id"] == best_hit["transcript_id"])

                good_isoforms.append({
                    "transcript_id": hit["transcript_id"],
                    "e_value": hit["e_value"],
                    "hmm_coverage": hit_hmm_cov,
                    "protein_coverage": hit_prot_cov,
                    "protein_length": hit["tlen"],
                    "is_best_hit": is_best,
                })

        return {
            "status": "hit_found",
            "best_hit": best_hit,
            "total_hits": len(hits),
            "hmm_coverage": hmm_coverage,
            "protein_coverage": protein_coverage,
            "issues": issues,
            "good_isoforms": good_isoforms,  # Isoforms that would satisfy BUSCO
        }

    except Exception as e:  # pylint: disable=broad-except
        return {"status": "parse_error", "error": str(e)}


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
    Get all genes overlapping a locus with transcript and translation details.

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
        List[Dict]: List of genes at locus with transcript details including protein lengths
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
                g.canonical_transcript_id,
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

        # Get transcripts with translation details for each gene
        # Note: We calculate protein length from exon sequences
        for gene in genes:
            cursor.execute(
                """
                SELECT
                    t.transcript_id,
                    t.stable_id,
                    t.seq_region_start,
                    t.seq_region_end,
                    t.biotype,
                    tl.translation_id,
                    tl.seq_start as cds_start,
                    tl.start_exon_id,
                    tl.seq_end as cds_end,
                    tl.end_exon_id
                FROM transcript t
                LEFT JOIN translation tl ON t.transcript_id = tl.transcript_id
                WHERE t.gene_id = %s
                ORDER BY t.transcript_id
                """,
                (gene["gene_id"],),
            )
            transcripts = cursor.fetchall()

            # Calculate protein length from coding exon sequences
            for transcript in transcripts:
                # Get all exons for this transcript in order
                cursor.execute(
                    """
                    SELECT
                        e.exon_id,
                        e.seq_region_start,
                        e.seq_region_end,
                        et.rank
                    FROM exon_transcript et
                    JOIN exon e ON et.exon_id = e.exon_id
                    WHERE et.transcript_id = %s
                    ORDER BY et.rank
                    """,
                    (transcript["transcript_id"],),
                )
                exons = cursor.fetchall()
                # Store exons in transcript for later analysis
                transcript["exons"] = exons

                if transcript.get("translation_id"):
                    # Calculate CDS length accounting for translation boundaries
                    cds_length = 0
                    start_exon_id = transcript.get("start_exon_id")
                    end_exon_id = transcript.get("end_exon_id")
                    cds_start = transcript.get("cds_start")
                    cds_end = transcript.get("cds_end")

                    in_cds = False
                    for exon in exons:
                        exon_length = exon["seq_region_end"] - exon["seq_region_start"] + 1

                        if exon["exon_id"] == start_exon_id:
                            # First coding exon - adjust for cds_start
                            in_cds = True
                            if exon["exon_id"] == end_exon_id:
                                # Single exon CDS
                                cds_length += cds_end - cds_start + 1
                            else:
                                cds_length += exon_length - cds_start + 1
                        elif exon["exon_id"] == end_exon_id:
                            # Last coding exon - adjust for cds_end
                            cds_length += cds_end
                            in_cds = False
                        elif in_cds:
                            # Middle coding exon - use full length
                            cds_length += exon_length

                    if cds_length > 0:
                        # Protein length = CDS length / 3 (codons)
                        transcript["protein_length"] = cds_length // 3
                    else:
                        transcript["protein_length"] = None
                else:
                    transcript["protein_length"] = None

            gene["transcripts"] = transcripts

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
    problems: List[Dict], host: str, port: int, user: str, password: str, core_db: str, layer_db: str, hmmer_dir: Optional[str] = None
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
        hmmer_dir (Optional[str]): Path to HMMER output directory

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

        # Parse HMMER output if available and genes exist in core
        problem["hmmer_results"] = {}
        if hmmer_dir and problem["core_genes"]:
            hmmer_file = Path(hmmer_dir) / f"{problem['busco_id']}.out"
            core_gene_ids = [g["stable_id"] for g in problem["core_genes"] if g.get("stable_id")]
            if core_gene_ids:
                problem["hmmer_results"] = parse_hmmer_output(str(hmmer_file), core_gene_ids)

        # Classify the problem
        classify_problem(problem)

        # Identify alternative isoforms worth testing
        problem["alternative_isoform_suggestions"] = identify_alternative_isoforms(problem)

    print(f"\nCompleted evidence audit for {total} loci")
    return problems


def identify_alternative_isoforms(problem: Dict) -> List[Dict]:
    """
    Identify alternative isoforms that might be worth testing for BUSCO.

    Suggests non-canonical isoforms that have substantially different protein sequences
    based on length differences and unique exon coverage.

    Args:
        problem (Dict): Problem locus data

    Returns:
        List[Dict]: Suggested alternative isoforms with rationale
    """
    suggestions = []

    for gene in problem.get("core_genes", []):
        canonical_id = gene.get("canonical_transcript_id")
        transcripts = gene.get("transcripts", [])

        if not transcripts or not canonical_id:
            continue

        # Find canonical transcript
        canonical = None
        for t in transcripts:
            if t["transcript_id"] == canonical_id:
                canonical = t
                break

        if not canonical or not canonical.get("protein_length"):
            continue

        canonical_length = canonical["protein_length"]

        # Get canonical exon positions
        canonical_exons = canonical.get("exons", [])
        canonical_positions = set()
        for exon in canonical_exons:
            for pos in range(exon["seq_region_start"], exon["seq_region_end"] + 1):
                canonical_positions.add(pos)

        # Find alternative isoforms with significantly different protein lengths
        for transcript in transcripts:
            if transcript["transcript_id"] == canonical_id:
                continue  # Skip canonical

            if not transcript.get("protein_length"):
                continue  # Skip non-coding

            alt_length = transcript["protein_length"]
            length_diff_pct = abs(alt_length - canonical_length) / canonical_length * 100

            # Suggest if >20% length difference or substantially longer
            if length_diff_pct > 20 or (alt_length > canonical_length * 1.5):
                reason = []
                if alt_length > canonical_length:
                    reason.append(f"{length_diff_pct:.0f}% longer ({alt_length} vs {canonical_length} aa)")
                else:
                    reason.append(f"{length_diff_pct:.0f}% shorter ({alt_length} vs {canonical_length} aa)")

                # Calculate unique genomic positions covered by alternative isoform
                alt_exons = transcript.get("exons", [])
                alt_positions = set()
                for exon in alt_exons:
                    for pos in range(exon["seq_region_start"], exon["seq_region_end"] + 1):
                        alt_positions.add(pos)

                # Find unique positions in alternative that canonical doesn't have
                unique_to_alt = alt_positions - canonical_positions
                unique_to_canonical = canonical_positions - alt_positions
                shared_positions = alt_positions & canonical_positions

                if unique_to_alt:
                    reason.append(f"{len(unique_to_alt)}bp unique to alt")
                if unique_to_canonical:
                    reason.append(f"{len(unique_to_canonical)}bp unique to canonical")

                # Calculate overlap percentage
                if alt_positions:
                    overlap_pct = len(shared_positions) / len(alt_positions) * 100
                    reason.append(f"{overlap_pct:.0f}% overlap")

                suggestions.append({
                    "gene_id": gene["stable_id"],
                    "transcript_id": transcript["stable_id"],
                    "protein_length": alt_length,
                    "canonical_length": canonical_length,
                    "length_diff_pct": length_diff_pct,
                    "unique_bp_to_alt": len(unique_to_alt),
                    "unique_bp_to_canonical": len(unique_to_canonical),
                    "overlap_pct": overlap_pct if alt_positions else 0,
                    "reason": "; ".join(reason),
                })

    return suggestions


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


def generate_report(problems: List[Dict], output_tsv: str, output_json: Optional[str] = None, busco_metadata: Optional[Dict] = None):
    """
    Generate TSV report and optional JSON output.

    Args:
        problems (List[Dict]): List of problem loci with evidence
        output_tsv (str): Output TSV file path
        output_json (Optional[str]): Optional JSON output file path
        busco_metadata (Optional[Dict]): BUSCO metadata (expected lengths, etc.)
    """
    if busco_metadata is None:
        busco_metadata = {}
    # Generate TSV report
    with open(output_tsv, "w", newline="") as f:  # pylint: disable=unspecified-encoding
        writer = csv.writer(f, delimiter="\t")

        # Header
        writer.writerow(
            [
                "BUSCO_ID",
                "OrthoDB_Group",
                "TaxID",
                "Expected_Length",
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
                "HMMER_Status",
                "HMMER_E-value",
                "HMMER_Coverage",
                "HMMER_Issues",
                "Alternative_Isoforms",
                "Summary",
            ]
        )

        # Data rows
        for problem in problems:
            # Get BUSCO metadata
            busco_id = problem["busco_id"]
            meta = busco_metadata.get(busco_id, {})
            orthodb_group = meta.get("orthodb_group", "-")
            taxid = meta.get("taxid", "-")
            expected_length = meta.get("length")
            expected_length_str = f"{expected_length:.0f}" if expected_length else "-"

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

            # Parse HMMER results
            hmmer = problem.get("hmmer_results", {})
            hmmer_status = hmmer.get("status", "-")
            alternative_isoforms_str = "-"

            if hmmer_status == "hit_found":
                best_hit = hmmer["best_hit"]
                hmmer_evalue = f"{best_hit['e_value']:.2e}"
                hmmer_coverage = f"HMM:{hmmer['hmm_coverage']:.0f}%,Prot:{hmmer['protein_coverage']:.0f}%"
                hmmer_issues = ",".join(hmmer["issues"]) if hmmer["issues"] else "none"

                # Format alternative isoforms
                good_isoforms = hmmer.get("good_isoforms", [])
                if good_isoforms:
                    # Show up to 3 good isoforms with their stats
                    # Mark the best hit with an asterisk to show it's likely the canonical
                    isoform_strs = []
                    for iso in good_isoforms[:3]:
                        marker = "*" if iso.get("is_best_hit") else ""
                        isoform_strs.append(
                            f"{iso['transcript_id']}{marker}(E:{iso['e_value']:.1e},HMM:{iso['hmm_coverage']:.0f}%)"
                        )
                    alternative_isoforms_str = "; ".join(isoform_strs)
                    if len(good_isoforms) > 3:
                        alternative_isoforms_str += f" +{len(good_isoforms)-3} more"
            else:
                hmmer_evalue = "-"
                hmmer_coverage = "-"
                hmmer_issues = hmmer_status

            # Simplified summary based on classification
            classification = problem["classification"]
            if classification == "PROTEIN_CODING_BUT_NOT_BUSCO_SATISFYING":
                if hmmer_status == "hit_found" and hmmer["issues"]:
                    good_isoforms = hmmer.get("good_isoforms", [])
                    if good_isoforms:
                        summary = f"{core_protein_coding} protein_coding gene(s), HMMER: {', '.join(hmmer['issues'])} BUT {len(good_isoforms)} good isoform(s) exist!"
                    else:
                        summary = f"{core_protein_coding} protein_coding gene(s), HMMER: {', '.join(hmmer['issues'])}"
                else:
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
                    orthodb_group,
                    taxid,
                    expected_length_str,
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
                    hmmer_status,
                    hmmer_evalue,
                    hmmer_coverage,
                    hmmer_issues,
                    alternative_isoforms_str,
                    summary,
                ]
            )

    print(f"TSV report written to: {output_tsv}")

    # Generate alternative isoforms suggestions file
    isoform_file = output_tsv.replace(".tsv", "_isoform_suggestions.tsv")
    write_isoform_suggestions(problems, isoform_file)

    # Generate summary statistics
    print_summary(problems)

    # Generate JSON output if requested
    if output_json:
        generate_json_report(problems, output_json)


def write_isoform_suggestions(problems: List[Dict], output_file: str):
    """
    Write alternative isoform suggestions to a separate TSV file.

    Args:
        problems (List[Dict]): List of problem loci
        output_file (str): Output file path for isoform suggestions
    """
    all_suggestions = []

    for problem in problems:
        suggestions = problem.get("alternative_isoform_suggestions", [])
        if suggestions:
            for sug in suggestions:
                all_suggestions.append({
                    "busco_id": problem["busco_id"],
                    "chromosome": problem["sequence"],
                    "locus_start": problem["start"],
                    "locus_end": problem["end"],
                    **sug
                })

    if not all_suggestions:
        print("No alternative isoform suggestions to write")
        return

    with open(output_file, "w", newline="") as f:  # pylint: disable=unspecified-encoding
        writer = csv.writer(f, delimiter="\t")

        # Header
        writer.writerow([
            "BUSCO_ID",
            "Chromosome",
            "Locus_Start",
            "Locus_End",
            "Gene_ID",
            "Alternative_Transcript_ID",
            "Alt_Protein_Length",
            "Canonical_Protein_Length",
            "Length_Diff_Percent",
            "Unique_bp_to_Alt",
            "Unique_bp_to_Canonical",
            "Overlap_Percent",
            "Reason",
        ])

        # Data rows
        for sug in all_suggestions:
            writer.writerow([
                sug["busco_id"],
                sug["chromosome"],
                sug["locus_start"],
                sug["locus_end"],
                sug["gene_id"],
                sug["transcript_id"],
                sug["protein_length"],
                sug["canonical_length"],
                f"{sug['length_diff_pct']:.1f}%",
                sug["unique_bp_to_alt"],
                sug["unique_bp_to_canonical"],
                f"{sug['overlap_pct']:.1f}%",
                sug["reason"],
            ])

    print(f"Alternative isoform suggestions written to: {output_file}")
    print(f"  Total suggestions: {len(all_suggestions)} isoforms across {len([p for p in problems if p.get('alternative_isoform_suggestions')])} genes")


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

    # HMMER isoform statistics
    genes_with_hmmer = sum(1 for p in problems if p.get("hmmer_results", {}).get("status") == "hit_found")
    genes_with_good_isoforms = sum(1 for p in problems if p.get("hmmer_results", {}).get("good_isoforms"))
    total_good_isoforms = sum(len(p.get("hmmer_results", {}).get("good_isoforms", [])) for p in problems)

    if genes_with_hmmer > 0:
        print(f"\nHMMER Analysis Results:")
        print(f"  Genes with HMMER hits: {genes_with_hmmer}")
        print(f"  Genes with alternative good isoforms: {genes_with_good_isoforms} ({100*genes_with_good_isoforms/genes_with_hmmer:.1f}%)")
        print(f"  Total good alternative isoforms found: {total_good_isoforms}")
        if genes_with_good_isoforms > 0:
            print(f"  Average good isoforms per gene: {total_good_isoforms/genes_with_good_isoforms:.1f}")

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

    # HMMER analysis
    parser.add_argument(
        "--hmmer_dir",
        type=str,
        help="Path to HMMER output directory (e.g., busco_protein/run_actinopterygii_odb10/hmmer_output/initial_run_results/)",
    )

    # BUSCO dataset metadata
    parser.add_argument(
        "--busco_dataset_dir",
        type=str,
        help="Path to BUSCO lineage dataset directory (e.g., actinopterygii_odb10) for OrthoDB metadata",
    )

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

    # Extract BUSCO metadata if dataset directory provided
    busco_metadata = {}
    if args.busco_dataset_dir:
        print("\nExtracting BUSCO metadata...")
        busco_ids = [p["busco_id"] for p in problems]
        busco_metadata = get_busco_metadata(args.busco_dataset_dir, busco_ids)
        print(f"Retrieved metadata for {len(busco_metadata)} BUSCO IDs")

    # Audit evidence
    print("\nAuditing evidence at problematic loci...")
    problems = audit_evidence(problems, args.host, args.port, args.user, args.password, core_db, layer_db, args.hmmer_dir)

    # Generate reports
    print("\nGenerating reports...")
    generate_report(problems, args.output_tsv, args.output_json, busco_metadata)

    print("\nAudit complete!")


if __name__ == "__main__":
    main()

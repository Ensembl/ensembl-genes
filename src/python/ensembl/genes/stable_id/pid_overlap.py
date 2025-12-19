"""
Script to compare protein IDs between two GFF3 files and generate a report
"""

import csv
import re
import pandas as pd


# Function to extract protein information from a GFF3 file
# Merges all CDS records for the same protein into a single entry
def extract_proteins_from_gff3(gff3_file: str) -> pd.DataFrame:
    """
    Extract protein information from GFF3 file.

    Args:
        gff3_file (str): Path to the GFF3 file.

    Returns:
        pd.DataFrame: DataFrame containing protein information.
    """
    protein_data = {}
    with open(gff3_file, "r") as file:  # pylint: disable=unspecified-encoding
        for line in file:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            if columns[2] == "CDS":
                attributes = columns[8]
                match = re.search(r"ID=CDS:([^;]+)", attributes)
                if match:
                    protein_id = match.group(1)
                    if protein_id not in protein_data:
                        protein_data[protein_id] = {
                            "ProteinID": protein_id,
                            "Chromosome": columns[0],
                            "Start": int(columns[3]),
                            "End": int(columns[4]),
                            "Strand": columns[6],
                            "Attributes": [attributes],
                        }
                    else:
                        # Update the start and end to encompass all CDS regions for the protein
                        protein_data[protein_id]["Start"] = min(
                            protein_data[protein_id]["Start"], int(columns[3])
                        )
                        protein_data[protein_id]["End"] = max(
                            protein_data[protein_id]["End"], int(columns[4])
                        )
                        protein_data[protein_id]["Attributes"].append(attributes)  # type: ignore[union-attr] # pylint:disable=line-too-long
    # Flatten attributes into a single string for each protein
    for protein in protein_data.values():
        protein["Attributes"] = ";".join(protein["Attributes"])  # type: ignore[arg-type]
    return pd.DataFrame(protein_data.values())


# Function to compute shared protein statistics and generate the report
def compare_protein_sets(
    gff3_file_old: str,  # pylint:disable=redefined-outer-name
    gff3_file_new: str,  # pylint:disable=redefined-outer-name
    output_report: str,  # pylint:disable=redefined-outer-name
) -> None:  # pylint:disable=redefined-outer-name
    """
    Compare protein sets from two GFF3 files and generate a report.

    Args:
        gff3_file_old (str): Path to the older GFF3 file.
        gff3_file_new (str): Path to the newer GFF3 file.
        output_report (str): Path to save the output report.

    Returns:
        None
    """
    # Extract protein information from both files
    proteins_old = extract_proteins_from_gff3(gff3_file_old)
    proteins_new = extract_proteins_from_gff3(gff3_file_new)

    # Find shared protein IDs
    shared_proteins = pd.merge(
        proteins_old, proteins_new, on="ProteinID", suffixes=("_old", "_new")
    )

    # Compute proportions
    proportion_old = (
        len(shared_proteins) / len(proteins_old) if len(proteins_old) > 0 else 0
    )
    proportion_new = (
        len(shared_proteins) / len(proteins_new) if len(proteins_new) > 0 else 0
    )

    print(
        f"Proportion of proteins in the old GFF3 that have shared IDs: {proportion_old:.2%}"
    )
    print(
        f"Proportion of proteins in the new GFF3 that have shared IDs: {proportion_new:.2%}"
    )

    # Add columns for genomic locations in the report
    shared_proteins = shared_proteins[
        [
            "ProteinID",
            "Chromosome_old",
            "Start_old",
            "End_old",
            "Strand_old",
            "Chromosome_new",
            "Start_new",
            "End_new",
            "Strand_new",
        ]
    ]

    # Save report
    shared_proteins.to_csv(
        output_report, index=False, sep="\t", quoting=csv.QUOTE_MINIMAL
    )
    print(f"Report saved to {output_report}")


# Example usage
if __name__ == "__main__":
    gff3_file_old = input("Enter the path to the older GFF3 file: ").strip()
    gff3_file_new = input("Enter the path to the newer GFF3 file: ").strip()
    output_report = input("Enter the path for the output report file: ").strip()

    compare_protein_sets(gff3_file_old, gff3_file_new, output_report)

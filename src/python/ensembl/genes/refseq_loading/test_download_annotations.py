#!/usr/bin/env python3
"""
Test script to download annotation data for specific RefSeq assemblies (GCF IDs).
"""

from refseq2ensembl import download_annotations


def main():
    base_dir = "refseq_data"

    print("Testing annotation downloads using GCF accessions...\n")

    assemblies_to_test = [
        "GCF_000001635.27",  # Mus musculus (GRCm39)
        "GCF_000001405.40"   # Homo sapiens (GRCh38)
    ]

    for asm in assemblies_to_test:
        print(f"→ Downloading annotation for: {asm}")
        try:
            download_annotations(
                base_dir=base_dir,
                assembly_acc=asm
            )
            print(f"[✓] Finished downloading for {asm}\n")
        except Exception as e:
            print(f"[✗] Failed to download {asm}: {e}\n")


if __name__ == "__main__":
    main()

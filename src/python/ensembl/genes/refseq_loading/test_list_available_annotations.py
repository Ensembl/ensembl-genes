#!/usr/bin/env python
"""
Smoke-test for list_available_annotations().

Place this beside the module that defines the function (e.g. refseq2ensembl.py)
and run:  python test_list_available_annotations.py
"""

import sys
from refseq2ensembl import list_available_annotations

def main():
    # Go through all groups and print every species (no truncation)
    listings = list_available_annotations(
        max_print=None,
        return_dict=True
    )

    print("\nSUMMARY")
    print("-------")
    for group, species_list in listings.items():
        print(f"{group}: {len(species_list)} species listed")
        if not species_list:
            print(f"[WARNING] No species found for group {group}")

    print("\nALL GROUPS PROCESSED")


if __name__ == "__main__":
    main()

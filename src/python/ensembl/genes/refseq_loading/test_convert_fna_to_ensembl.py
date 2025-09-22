#!/usr/bin/env python3

import os
import sys

from refseq2ensembl import convert_fna_headers  # Assumes function is in your module

def main():
    if len(sys.argv) != 4:
        print("Usage: test_convert_fna_headers.py <input.fna> <assembly_report.txt> <output.fna>")
        sys.exit(1)

    fna_in = sys.argv[1]
    rpt = sys.argv[2]
    fna_out = sys.argv[3]

    print("→ Converting FNA headers:")
    print(f"   Input      : {fna_in}")
    print(f"   Report     : {rpt}")
    print(f"   Output     : {fna_out}")

    if not os.path.exists(fna_in):
        print(f"[Error] Input FNA file not found: {fna_in}")
        sys.exit(1)
    if not os.path.exists(rpt):
        print(f"[Error] Assembly report file not found: {rpt}")
        sys.exit(1)

    convert_fna_headers(fna_in, rpt, fna_out)

    print("✓ Conversion completed.")

if __name__ == "__main__":
    main()

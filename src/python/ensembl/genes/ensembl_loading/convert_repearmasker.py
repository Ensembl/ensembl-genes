"""Command-line wrapper for RefSeq RepeatMasker-to-GTF conversion."""

from __future__ import annotations

import argparse

from .refseq_conversion import convert_repeatmasker_to_gtf


def main() -> int:
    """Convert one RepeatMasker output file to GTF."""

    parser = argparse.ArgumentParser()
    parser.add_argument("repeatmasker", help="Input *_rm.out or *_rm.out.gz")
    parser.add_argument("assembly_report", help="NCBI assembly report")
    parser.add_argument("output", help="Output RepeatMasker GTF")
    args = parser.parse_args()
    convert_repeatmasker_to_gtf(
        args.repeatmasker,
        args.assembly_report,
        args.output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

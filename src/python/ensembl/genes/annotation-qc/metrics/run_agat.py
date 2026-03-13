"""
GFF3 annotation statistics pipeline using AGAT.

This script processes a GFF3 genome annotation file by:
  1. Splitting it into a canonical-transcript-only file and a full (all-transcripts) file.
  2. Running AGAT (Another Gff Analysis Toolkit) statistics on each file via Singularity.
  3. Parsing the resulting AGAT statistics reports into structured CSV files.

Typical usage:
    python script.py --gff3 annotation.gff3 --outdir ./output

Dependencies:
    - Singularity (with AGAT SIF image at the path defined by AGAT_SIF)
    - AGAT >= 1.4.2

Output files (written to --outdir):
    - <stem>_canonical.gff3         : GFF3 with only Ensembl canonical transcripts
    - <stem>_all.gff3               : Full GFF3 (copy of input)
    - <stem>_canonical_agat_stats.txt / .csv
    - <stem>_all_agat_stats.txt     / .csv
"""

import argparse
import os
import subprocess
from pathlib import Path
import csv
import re

AGAT_SIF = "/hps/nobackup/flicek/ensembl/genebuild/lazar/metrics/agat_1.4.2--pl5321hdfd78af_0.sif"


def run_agat(
    gff_file: Path, output_txt: Path, workdir: Path, my_level: Path, agat_path: Path
) -> None:
    """Run AGAT statistics on a GFF3 file via Singularity.

    Executes agat_sp_statistics.pl inside the AGAT Singularity container,
    writing a human-readable statistics report to output_txt.
    The step is skipped if output_txt already exists.

    Args:
        gff_file: Path to the GFF3 file to analyse.
        output_txt: Destination path for the AGAT statistics report.
        workdir: Working direcory for AGAT (log storing).
        my_level: Path to custom feture_levels.yaml.
        agat_path: Path to AGAT singularity.

    Returns:
        None

    Raises:
        subprocess.CalledProcessError: If the Singularity/AGAT command exits
            with a non-zero return code.
    """

    if output_txt.exists():
        print(f"{output_txt} exists. Skipping AGAT.")
        return

    print(f"Running AGAT for {gff_file.name}")

    # Base command
    cmd = [
        "singularity",
        "exec",
        "--cleanenv",
    ]

    # Add bind if a custom levels file is provided
    if my_level:
        bind_spec = f"{my_level.resolve()}:/usr/local/lib/perl5/site_perl/auto/share/dist/AGAT/feature_levels.yaml"
        cmd += ["-B", bind_spec]

    # Add SIF and AGAT command
    cmd += [
        agat_path,
        "agat_sp_statistics.pl",
        "--gff",
        str(gff_file),
        "-o",
        output_txt.name,
    ]

    subprocess.run(cmd, check=True, cwd=workdir)


def parse_agat_stats(stats_file: Path) -> list[tuple[str, str]]:
    """Parse an AGAT statistics report into (metric, value) pairs.

    Processes the plain-text report produced by agat_sp_statistics.pl,
    handling:
      - Section headers (delimited by dashes) that are prepended to metric names.
      - "Shortest isoforms excluded" blocks, which append the suffix
        "_minus_shortest_isoform" to affected metric names.

    Metric names are normalised to lowercase with spaces replaced by
    underscores, and section prefixes are joined with an underscore.

    Args:
        stats_file: Path to the AGAT statistics text file to parse.

    Returns:
        A list of (metric, value) tuples where *metric* is the
        normalised metric name (str) and *value* is the raw numeric
        string as it appears in the report.
    """

    rows = []
    current_section = None
    minus_isoform = False

    section_pattern = re.compile(r"-+\s*(.+?)\s*-+")
    value_pattern = re.compile(r"(.+?)\s+([0-9.]+)$")

    with open(stats_file) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # Detect isoform message
            if "shortest isoforms excluded" in line.lower():
                minus_isoform = True
                continue

            # Detect section header
            section_match = section_pattern.match(line)
            if section_match:
                current_section = section_match.group(1).strip()
                minus_isoform = False
                continue

            match = value_pattern.match(line)
            if match:
                metric = match.group(1).strip()
                value = match.group(2)

                if current_section:
                    metric = f"{current_section}_{metric}"

                metric = metric.replace(" ", "_").lower()

                if minus_isoform:
                    metric = f"{metric}_minus_shortest_isoform"

                rows.append((metric, value))

    return rows


def write_csv(rows: list[tuple[str, str]], output_csv: Path) -> None:
    """Write parsed AGAT statistics to a two-column CSV file.

    Creates a CSV with a ``metric`` and ``value`` header row followed by one
    row per entry in *rows*. The file is skipped if it already exists.

    Args:
        rows: Sequence of ``(metric, value)`` tuples as returned by
            :func:`parse_agat_stats`.
        output_csv: Destination path for the CSV file.

    Returns:
        None
    """

    if output_csv.exists():
        print(f"{output_csv} exists. Skipping CSV creation.")
        return

    print(f"Writing {output_csv}")

    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])

        for metric, value in rows:
            writer.writerow([metric, value])


def main() -> None:
    """Entry point for the GFF3 annotation statistics pipeline.

    Parses command-line arguments, orchestrates the GFF3 split, AGAT runs,
    and CSV export steps, skipping any step whose output already exists.

    Command-line arguments:
        --gff3 (str, required): Path to the input GFF3 annotation file.
        --outdir (str, optional): Directory for all output files. Defaults to
            the current working directory.

    Returns:
        None
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", required=True, help="Path to GFF3/GTF file")
    parser.add_argument("--outdir", default=".", help="Path to output directory")
    parser.add_argument(
        "--feature_levels", help="Path to custom fetaure_levels.yaml file"
    )
    parser.add_argument("--agat_path", help="Path to AGAT singularity")

    args = parser.parse_args()

    input_gff3 = Path(args.gff3)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    my_levels = Path(args.feature_levels)
    agat_path = Path(args.agat_path) if args.agat_path else Path(AGAT_SIF)

    full_stats_txt = outdir / f"{input_gff3.stem}_agat_stats.txt"

    run_agat(input_gff3, full_stats_txt, outdir, my_levels, agat_path)
    full_csv = outdir / f"{input_gff3.stem}_agat_stats.csv"
    if not full_csv.exists():
        full_rows = parse_agat_stats(full_stats_txt)
        write_csv(full_rows, full_csv)


if __name__ == "__main__":
    main()

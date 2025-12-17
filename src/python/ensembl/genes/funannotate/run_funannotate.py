import subprocess
from pathlib import Path
import argparse
import gzip
import shutil

def unzip_keep_original(file_path):
    """
    Unzip a .gz file but keep the original.
    Returns the path to the uncompressed file.
    """
    file_path = Path(file_path)
    unzipped_file = file_path.with_suffix('')  # remove .gz
    if not unzipped_file.exists():
        print(f"Unzipping {file_path} -> {unzipped_file}")
        with gzip.open(file_path, 'rb') as f_in, open(unzipped_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return unzipped_file

def run_funannotate(gcas, base_dir, output_dir, cpus=20):
    """
    Run funannotate predict for a list of GCAs.

    Parameters:
        gcas (list of str): List of GCA IDs.
        base_dir (str or Path): Base directory containing GCA folders.
        output_dir (str or Path): Base output directory where results will be stored.
        cpus (int): Number of CPUs to use.
    """
    base_dir = Path(base_dir)
    output_dir = Path(output_dir)
    funannotate_cmd = "/hps/nobackup/flicek/ensembl/genebuild/lazar/bbr/funannotate/funannotate-singularity" 

    for gca in gcas:
        gca_dir = base_dir / gca
        output_dir_gca = output_dir / gca
        output_dir_gca.mkdir(parents=True, exist_ok=True)

        # Find the genome FASTA file
        genome_files = list(gca_dir.glob("*_softmasked_toplevel.fa.gz"))
        if not genome_files:
            print(f"No genome FASTA found for {gca}, skipping...")
            continue
        genome_fa = unzip_keep_original(genome_files[0])

        # Extract species name from filename
        species_parts = genome_fa.stem.split("_")[:2]
        species = " ".join(part.capitalize() for part in species_parts)

        # Assume other files are in subfolders
        rna_bam_files = list((gca_dir / "star_output").glob("*.bam"))
        stringtie_gtf_files = list((gca_dir / "stringtie_output").glob("*.stringtie.gtf"))

        if not rna_bam_files or not stringtie_gtf_files:
            print(f"No RNA BAM or StringTie GTF found for {gca}, skipping...")
            continue

        rna_bam = rna_bam_files[0]
        stringtie_gtf = stringtie_gtf_files[0]
        protein_gtf = gca_dir / "genblast_output" / "annotation.gtf"
        trna_gtf = gca_dir / "trnascan_output" / "annotation.gtf"

        cmd = [
            funannotate_cmd, "predict",
            "-i", str(genome_fa),
            "-o", str(output_dir_gca),
            "--species", species,
            "--rna_bam", str(rna_bam),
            "--stringtie", str(stringtie_gtf),
            "--protein_alignments", str(protein_gtf),
            "--trnascan", str(trna_gtf),
            "--cpus", str(cpus)
        ]

        print(f"Running funannotate for {gca} ({species})...")
        subprocess.run(cmd, check=True)
        print(f"Finished {gca}\n")


def main():
    parser = argparse.ArgumentParser(description="Run funannotate for multiple GCAs")
    parser.add_argument("--gca_file", required=True, help="Text file with one GCA per line")
    parser.add_argument("--base_dir", required=True, help="Base directory containing GCA folders")
    parser.add_argument("--output_dir", required=True, help="Base output directory for funannotate results")
    parser.add_argument("--cpus", type=int, default=20, help="Number of CPUs to use")

    args = parser.parse_args()

    # Read GCAs from the text file
    with open(args.gca_file) as f:
        gcas = [line.strip() for line in f if line.strip()]

    run_funannotate(gcas, args.base_dir, args.output_dir, args.cpus)


if __name__ == "__main__":
    main()
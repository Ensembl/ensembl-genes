#!/bin/bash
#SBATCH --job-name=intron_prediction    # Job name
#SBATCH --output=lift_logs/output_%j.log     # Standard output (stdout) log file
#SBATCH --error=lift_logs/error_%j.log       # Standard error (stderr) log file
#SBATCH --ntasks=1                      # Number of tasks (1 job instance)
#SBATCH --cpus-per-task=10              # Request 10 CPU cores
#SBATCH --mem=30G                       # Request 20GB RAM
#SBATCH --time=12:00:00                 # Max runtime (hh:mm:ss) (adjust as needed)
#SBATCH --partition=production          # Partition name (adjust if needed)

# Capture arguments passed to the script
REF_GFF3=$1
REF_FA=$2
TAR_FA=$3
OUTPUT=$4
#GFF3_FILE=$1
#FASTA_FILE=$2
#CORES=${3:-${SLURM_CPUS_PER_TASK:-10}}
#CANONICAL_ONLY=${4:-1}
#USE_MOUSE=${5:-0}                 # Default: 0 (i.e. use Homo sapiens by default)
#REFSEQ=${6:-0}
#QUALITY_RATINGS=$7
#VERBOSE=${8:-"--verbose"}         # Default: Verbose mode enabled


# Load Python (if necessary, depending on your system)
#module load python/3.10.9  # Adjust Python version if needed

# Run the Python script with SLURM-allocated resources
lifton -f feature_types_small.txt -t 10 -g $REF_GFF3 -o $OUTPUT $REF_FA $TAR_FA

#lifton -f CDS,exon,five_prime_UTR,gene,lnc_RNA,miRNA,mRNA,ncRNA_gene,pseudogene,pseudogenic_transcript,rRNA,scRNA,snoRNA,three_prime_UTR -t 10 -g bref.gff3 cref.fa bref.fa

#!/usr/bin/env bash
# Download test data for pangenome mapping (chr20 subset)
#
# Downloads:
# - GRCh38 chr20 sequence
# - CHM13v2 chr20 sequence  
# - GENCODE v44 annotation (chr20 only)
#
# Usage: ./download_test_data.sh [output_dir]

set -euo pipefail

OUTPUT_DIR="${1:-data}"
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "Downloading test data for chr20"
echo "Output directory: ${OUTPUT_DIR}"
echo "=========================================="

# URLs
GRCH38_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz"
CHM13_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gff3.gz"

# Download GRCh38 chr20
echo ""
echo "1. Downloading GRCh38 chr20..."
if [[ ! -f "${OUTPUT_DIR}/GRCh38.chr20.fa" ]]; then
    curl -L "${GRCH38_URL}" | gunzip > "${OUTPUT_DIR}/GRCh38.chr20.fa"
    echo "   Done: ${OUTPUT_DIR}/GRCh38.chr20.fa"
else
    echo "   Already exists: ${OUTPUT_DIR}/GRCh38.chr20.fa"
fi

# Download CHM13 full genome and extract chr20
echo ""
echo "2. Downloading CHM13v2 (full genome, then extracting chr20)..."
if [[ ! -f "${OUTPUT_DIR}/CHM13v2.chr20.fa" ]]; then
    if [[ ! -f "${OUTPUT_DIR}/chm13v2.0.fa" ]]; then
        echo "   Downloading full CHM13v2 genome..."
        curl -L "${CHM13_URL}" | gunzip > "${OUTPUT_DIR}/chm13v2.0.fa"
    fi
    
    # Create index and extract chr20
    echo "   Indexing and extracting chr20..."
    samtools faidx "${OUTPUT_DIR}/chm13v2.0.fa"
    samtools faidx "${OUTPUT_DIR}/chm13v2.0.fa" chr20 > "${OUTPUT_DIR}/CHM13v2.chr20.fa"
    echo "   Done: ${OUTPUT_DIR}/CHM13v2.chr20.fa"
else
    echo "   Already exists: ${OUTPUT_DIR}/CHM13v2.chr20.fa"
fi

# Download GENCODE annotation and extract chr20
echo ""
echo "3. Downloading GENCODE v44 annotation (chr20 only)..."
if [[ ! -f "${OUTPUT_DIR}/gencode.v44.chr20.gff3" ]]; then
    if [[ ! -f "${OUTPUT_DIR}/gencode.v44.primary_assembly.annotation.gff3" ]]; then
        echo "   Downloading full annotation..."
        curl -L "${GENCODE_URL}" | gunzip > "${OUTPUT_DIR}/gencode.v44.primary_assembly.annotation.gff3"
    fi
    
    # Extract chr20 features
    echo "   Extracting chr20 features..."
    grep -E "^(#|chr20\t)" "${OUTPUT_DIR}/gencode.v44.primary_assembly.annotation.gff3" > "${OUTPUT_DIR}/gencode.v44.chr20.gff3"
    echo "   Done: ${OUTPUT_DIR}/gencode.v44.chr20.gff3"
else
    echo "   Already exists: ${OUTPUT_DIR}/gencode.v44.chr20.gff3"
fi

# Index FASTA files
echo ""
echo "4. Indexing FASTA files..."
samtools faidx "${OUTPUT_DIR}/GRCh38.chr20.fa" 2>/dev/null || true
samtools faidx "${OUTPUT_DIR}/CHM13v2.chr20.fa" 2>/dev/null || true

# Summary
echo ""
echo "=========================================="
echo "Download complete!"
echo "=========================================="
echo ""
echo "Files created:"
ls -lh "${OUTPUT_DIR}"/*.fa "${OUTPUT_DIR}"/*.gff3 2>/dev/null || true

echo ""
echo "To run the mapping pipeline:"
echo ""
echo "  python cli.py map \\"
echo "    --ref-fasta ${OUTPUT_DIR}/GRCh38.chr20.fa \\"
echo "    --ref-gff ${OUTPUT_DIR}/gencode.v44.chr20.gff3 \\"
echo "    --target-fasta ${OUTPUT_DIR}/CHM13v2.chr20.fa \\"
echo "    --output-gff ${OUTPUT_DIR}/output/chm13.chr20.mapped.gff3 \\"
echo "    --output-stats ${OUTPUT_DIR}/output/chr20.stats.json \\"
echo "    --chromosomes chr20 \\"
echo "    --threads 4"
echo ""

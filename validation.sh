#!/bin/bash
set -e

echo "🚀 Starting validation of the new configuration-driven YAML pipeline"

export PYTHONPATH=src/python
INPUT_FILE="dummy_input.txt"

echo "Creating test input file $INPUT_FILE..."
cat << EOF > $INPUT_FILE
homo_sapiens_core_110_38
caenorhabditis_elegans_core_57_110_282
GCA_009914755.4
EOF

echo ""
echo "🧪 Testing VGP/Standard configuration..."
python src/python/ensembl/genes/projects/generate_project_yaml.py --project vgp $INPUT_FILE --output vgp_test_final.yaml
head -n 20 vgp_test_final.yaml

echo ""
echo "🧪 Testing HPRC configuration (with BioSample scraping and submitter overrides)..."
python src/python/ensembl/genes/projects/generate_project_yaml.py --project hprc $INPUT_FILE --output hprc_test_final.yaml
head -n 20 hprc_test_final.yaml

echo ""
echo "🧪 Testing Mouse Pangenome configuration..."
python src/python/ensembl/genes/projects/generate_project_yaml.py --project mouse_genomes $INPUT_FILE --output mouse_test_final.yaml
head -n 20 mouse_test_final.yaml

echo ""
echo "✅ Validation complete!"

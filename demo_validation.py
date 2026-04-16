import sys
import yaml
import logging

from ensembl.genes.projects.config import get_project_config
from ensembl.genes.projects.generate_project_yaml import main
from unittest.mock import patch
from ensembl.genes.projects.registry.gb_tracker import GbTrackerClient

# We will create test_inputs_val.txt with just ONE released UUID (CBP)
with open('test_inputs_val.txt', 'w') as f:
    f.write('63990222-c62d-4918-8d10-dc4e806a1909\n') # GCA_007570785.1

original_fetch = GbTrackerClient.fetch_project_pre_releases

def mock_fetch(self, config):
    entries = original_fetch(self, config)
    # create a mock entry for GCA_007570785.1 (the released genome) appearing as if it was found in GB
    import copy
    if entries:
        mock_dup = copy.deepcopy(entries[0])
        mock_dup.accession = 'GCA_007570785.1'
        mock_dup.species_name = 'Mock Colliding Genome'
        entries.append(mock_dup)
    return entries

with patch.object(GbTrackerClient, 'fetch_project_pre_releases', new=mock_fetch):
    sys.argv = ['generate_project_yaml.py', '--project', 'cbp', 'test_inputs_val.txt', '--output', 'val_out.yaml']
    logging.basicConfig(level=logging.INFO)
    main()

with open('val_out.yaml') as f:
    data = yaml.safe_load(f)
    print("\n--- YAML OUTPUT HEAD ---")
    for doc in data[:3]:
        print(doc)
    print(f"Total entries in YAML: {len(data)}")

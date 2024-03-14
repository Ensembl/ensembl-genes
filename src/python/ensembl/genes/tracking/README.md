#

## About


## Running bioproject_tracking.py

**python bioproject_tracking.py  -h**

usage: bioproject_tracking.py [-h] [--haploid] [--report_file REPORT_FILE] [--rank RANK] bioproject_id

Fetch assembly accessions from NCBI BioProject and report the number of corresponding annotations in rapid.ensembl.org.

positional arguments:

  bioproject_id         NCBI BioProject ID

optional arguments:

  -h, --help                     show this help message and exit

      --haploid                  Fetch only haploid assemblies

      --report_file REPORT_FILE  Where to write report to

      --rank RANK                Taxonomic rank to classify

## Example Usage

I want to report the number of DToL assemblies for which I can find annotation on rapid.ensembl.org. I only want to report for primary haplotypes (use --haploid flag) and I want a breakdown by taxonomic class:

**python bioproject_tracking.py --haploid PRJEB40665 --rank class --report_file ./dtol_rank_report_file.csv**

Found 1637 assemblies under BioProject ID PRJEB40665

Found 631 annotations in rapid.ensembl.org for 612 unique species

Breakdown:

Counter({'Insecta': 591, 'Mammalia': 10, 'Magnoliopsida': 6, 'Actinopteri': 5, 'Gastropoda': 4, 'Bivalvia': 4, 'Amphibia': 2, 'Anthozoa': 2, 'Clitellata': 2, 'Aves': 2, 'Asteroidea': 1, 'Staurozoa': 1, 'Pilidiophora': 1})
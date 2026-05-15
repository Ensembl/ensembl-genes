# Bioproject Tracking

## About

**Bioproject Tracking** is a Python script that fetches assembly accessions from the NCBI BioProject, a specified Taxonomy ID, or a named project portal (e.g. HPRC) and correlates them with metadata for Ensembl. The script then writes out a summary of findings, including annotation data, an optional FTP link, and classification at a specified taxonomic rank.

The script relies on:

- The NCBI Datasets API to retrieve assembly accessions.
- A MySQL database (as configured in `bioproject_tracking_config.json`) to check how many of those assemblies have corresponding annotation in Ensembl.
- A second query to the NCBI Datasets API to retrieve taxonomy information at your specified taxonomic rank (e.g., `order`, `class`, `phylum`).

## Prerequisites

1. **Python 3.7+** (recommended).
2. **MySQL** server or a connection to a MySQL database containing the Ensembl metadata tables.
3. A valid `bioproject_tracking_config.json` configuration file, which should include:
   - MySQL connection details for your Ensembl metadata database (host, user, port, db_name).
   - The base URLs for the NCBI Datasets API queries (e.g., for BioProject and Taxon).
4. Install the Python dependencies listed in `requirements.txt`:
   ```bash
   pip install -r requirements.txt
   ```

## Running `bioproject_tracking.py`

```bash
python -m ensembl.genes.tracking.bioproject_tracking -h
```

### Usage

```
usage: bioproject_tracking.py [-h] [--bioproject_id BIOPROJECT_ID]
                              [--taxon_id TAXON_ID]
                              [--project_name PROJECT_NAME] [--haploid]
	                      [--report_file REPORT_FILE] [--classification]
			      [--rank RANK] [--ftp] [--pre_release]

This script fetches assembly accessions from NCBI BioProject (or a named
project portal) and reports the number of corresponding annotations in
the beta.ensembl.org metadata.

optional arguments:
  -h, --help            show this help message and exit
  --bioproject_id BIOPROJECT_ID NCBI BioProject ID
  --taxon_id TAXON_ID   Taxonomy ID
  --project_name PROJECT_NAME
                        Project-specific discovery (e.g. 'hprc').
                        Mutually exclusive with --bioproject_id / --taxon_id.
                        Currently supported: hprc
  --haploid             Fetch only haploid assemblies
  --report_file REPORT_FILE
  --classification      Provide breakdown of taxonomic classification [note: this will significantly slow down the script, only use if needed]
  --rank RANK
  --ftp                 Include FTP links in the report
  --pre_release         Include list of pre-release databases in the report
```

## Example Usage

You can combine flags to tailor your query. For example, to report the number of DToL assemblies for which you can find annotation on `beta.ensembl.org`, restricting to primary haplotypes (`--haploid`), provding taxonmic classification (`--classification`) and summarising by taxonomic rank=class (`--rank class`):

```bash
python -m ensembl.genes.tracking.bioproject_tracking --bioproject_id PRJEB40665 --haploid --classification --rank class --report_file ./dtol_rank_report_file.csv
```

Example output might look like:

```
Found 1637 assemblies under BioProject ID PRJEB40665
Found 631 annotations in rapid.ensembl.org for 612 unique species

Breakdown:
Counter({'Insecta': 591, 'Mammalia': 10, 'Magnoliopsida': 6, 'Actinopteri': 5,
         'Gastropoda': 4, 'Bivalvia': 4, 'Amphibia': 2, 'Anthozoa': 2,
         'Clitellata': 2, 'Aves': 2, 'Asteroidea': 1, 'Staurozoa': 1,
         'Pilidiophora': 1})
```

1. **`--haploid`** ensures only haploid assemblies are reported.
2. **`--classification`** tells the script to report taoxnomic classifications 
3. **`--rank`** (e.g., `class`) tells the script which taxonomic rank to retrieve from NCBI.
4. **`--report_file`** specifies where to save the tab-separated output.

### HPRC Project-Specific Discovery

To discover all HPRC release 2 assemblies dynamically from the HPRC data portal:

```bash
python -m ensembl.genes.tracking.bioproject_tracking --project_name hprc --report_file hprc_report.tsv
```

This will:
- Download the HPRC catalog JSON from the [HPRC data explorer](https://data.humanpangenome.org/assemblies) GitHub repo.
- Filter to release 2 assemblies with valid GCA accessions.
- De-duplicate and log how many GCAs were found.
- Feed those GCAs into the same Ensembl metadata lookup and report generation flow.

The output report has the same columns as the BioProject path.

---

**Note**: Make sure you update the `bioproject_tracking_config.json` file with your specific MySQL credentials and NCBI API endpoint details before running the script.

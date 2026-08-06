# How to install this repository

This project requires Python 3.12+.

## Basic installation

To install the package from GitHub using pip:

```bash
pip install git+https://github.com/Ensembl/ensembl-genes.git
```

## Development installation

Create and activate a virtual environment, then install in editable mode:

```bash
python3 -m venv .venv
source .venv/bin/activate
git clone https://github.com/Ensembl/ensembl-genes.git
pip install -e ensembl-genes/.[cicd,docs]
```

Documentation is generated with Sphinx. The source lives in `sphinx`.

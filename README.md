# ensembl-genes
Python Ensembl Gene Annotation code source repository

This repository contains Python code supporting Ensembl Genes workflows, including data processing, validation, analysis, and pipeline utilities used across gene build and annotation tasks.

The repository is structured as an installable Python package, enabling code reuse across pipelines, notebooks, and command-line tools.

## Installation & Dependency Management

All Python dependencies should be declared in pyproject.toml.
This ensures the repository can be installed consistently and reproducibly.

For development, install the package in editable mode:

pip install -e .


This allows:

Immediate reflection of code changes

Use of the package across scripts and workflows

Cleaner dependency management via pyproject.toml

Please do not manage dependencies via ad-hoc requirements.txt files unless explicitly required for external tooling.

## Repository Structure

Each top-level directory corresponds to a specific functional area or workflow

Each directory contains its own README describing:

Purpose

Inputs / outputs

Usage details

Any workflow-specific assumptions

This top-level README provides context only — detailed documentation can be found in each folder.



## Contributing Guidelines 

When adding or modifying code:

Add dependencies to pyproject.toml

Keep functionality scoped and reusable

Update or add a README in the relevant directory if behavior changes
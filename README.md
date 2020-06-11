# Ensembl Python Template

[![Documentation Status](https://readthedocs.org/projects/template-python/badge/?version=latest)](http://template-python.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/Ensembl/templat-python.svg?branch=master)](https://travis-ci.org/Ensembl/template-python)

Example Ensembl python module.

This repo structure can be forked and used as the base template for new tools and workflows. It should have all of the base functionality and is set up for unit testing and with pylint to ensure code clarity.

Creating a New Repository
-------------------------

```
git clone --depth 1 -b master https://github.com/Ensembl/template-python.git

rm -rf template-python/.git
mv template-python <NEW_PROJECT_NAME>
cd <NEW_PROJECT_NAME>

mv ensembl_template_py <NEW_PROJECT_NAME>

git init
git add .
git commit -m 'Initial commit'

git remote add origin https://github.com/Ensembl/<NEW_PROJECT_NAME>.git
git remote -v

git push origin master
```

Once cloned the following files will need to be customised:
- README.md
- NOTICE
- setup.py
- `__init__.py`
- docs/conf.py

The files in docs contain boilerplate for the installation instructions. New files will need to be added for all the modules so that the function documentation can be imported correctly.

# Requirements
- pyenv and pyenv-virtualenv
- Python 2.7.12+
- Python 3.6+
- Python Modules:
  - pylint
  - pytest

Installation
------------

Directly from GitHub:

```
cd ${HOME}/code

git clone https://github.com/Ensembl/<NEW_PROJECT_NAME>.git

cd <NEW_PROJECT_NAME>
```

Create the Python environment

```
pyenv-virtualenv 2.7.12 test_area
pyenv activate test_area
pip install -e .
pip install -r requirements.txt
```

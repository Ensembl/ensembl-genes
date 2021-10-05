# Gene Symbol Classifier eHive pipeline
# https://github.com/Ensembl/ensembl-genes/pipelines/gene_symbol_classifier

FROM python:3.9

LABEL maintainer="William Stark <william@ebi.ac.uk>"

# disable caching
ENV \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

# Poetry installation environment variables
ENV \
    POETRY_HOME="/opt/poetry" \
    POETRY_VERSION=1.1.11 \
    POETRY_VIRTUALENVS_CREATE=false

# add Poetry bin directory to PATH
ENV PATH="${POETRY_HOME}/bin:${PATH}"

# install Poetry
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python -

# specify working directory
WORKDIR /app

# copy project dependencies files
COPY \
    pyproject.toml \
    poetry.lock \
    /app/

# install project dependencies
RUN poetry install --no-dev

# create /app/data directory
RUN mkdir --verbose /app/data

# copy pipeline program files
COPY \
    filter_assignments.py \
    /app/

VOLUME /app/data

ENTRYPOINT ["python", "/app/filter_assignments.py"]

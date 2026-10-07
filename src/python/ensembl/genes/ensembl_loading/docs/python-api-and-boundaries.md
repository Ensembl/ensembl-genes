# Python API And Boundaries

## Python Functions

The CLI is a wrapper around reusable functions.
The DB-loading functions require all MySQL connection fields explicitly,
including `db_port`.

Generic loading functions:

```python
from ensembl.genes.ensembl_loading.gff_core_loader import (
    load_gff_features_to_core,
    load_to_ensembl_core,
    prepare_annotation_for_load,
)
from ensembl.genes.ensembl_loading.gff_source_config import get_source_config

source_config = get_source_config("generic")

summary = load_gff_features_to_core(
    gff_path="annotations.gff3",
    db_name="example_core",
    db_host="mysql-host",
    db_user="write-user",
    db_password="password",
    db_port=3306,
    source_config=source_config,
)
```

RefSeq functions:

```python
from ensembl.genes.ensembl_loading.refseq_ncbi import (
    list_available_annotations,
    download_annotations,
)
from ensembl.genes.ensembl_loading.refseq_conversion import (
    convert_fna_headers,
    convert_gff_to_ensembl,
)
```

These functions use standard `logging` rather than printing directly, so callers
can configure verbosity, capture logs in tests, or integrate them into a larger
pipeline.

The public loading facade is `gff_core_loader.py`. Its implementation delegates
assembly-report metadata to `gff_metadata.py`, GFF feature parsing and
normalization to `gff_annotation.py`, and schema/database row operations to
`gff_core_database.py`. Existing callers can continue importing the documented
loader functions from `gff_core_loader`.

For RefSeq GCF assemblies, resolve the official annotation URL first:

```bash
python -m ensembl.genes.ensembl_loading.ncbi_annotation_report \
  GCF_000001635.27
```

The loader calls this resolver automatically for GCF assemblies and stores the
returned URL as `genebuild.provider_url`. The resolver is isolated so it can be
removed later without changing the core loader.

## Current Boundaries

The current modular loader is generic at the GFF feature-loading layer. RefSeq
download and assembly-report conversion remain RefSeq-specific.

The generic loader is suitable for another source when:

1. The source can be represented as gene, transcript-like, exon, and CDS rows.
2. IDs and Parent values can be normalized with prefix stripping.
3. Biotypes can be resolved from feature type and attributes using
   `GffSourceConfig`.
4. GFF seqids already match target core seq_regions, or a matching FASTA is
   available for `create-core`.

Parser or loader changes may be needed for:

1. Multiple Parent relationships.
2. URL-decoded or escaped GFF attributes.
3. Unstranded features that should not be represented as `-1`.
4. Feature types that need dedicated Ensembl tables outside gene, transcript,
   exon, exon_transcript, and translation.
5. Xref/display label insertion.
6. Stable ID version tables or external DB metadata.

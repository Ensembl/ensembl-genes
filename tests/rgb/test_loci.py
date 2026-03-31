from __future__ import annotations

import pandas as pd

from ensembl.genes.rgb.loci import build_loci


def test_build_loci_basic():
    core = pd.DataFrame(
        [
            {
                "gene_id": 1,
                "seq_region_name": "chr1",
                "seq_region_start": 100,
                "seq_region_end": 200,
                "seq_region_strand": 1,
            },
            {
                "gene_id": 2,
                "seq_region_name": "chr1",
                "seq_region_start": 500,
                "seq_region_end": 600,
                "seq_region_strand": 1,
            },
        ]
    )

    layer = pd.DataFrame(
        [
            {
                "gene_id": 10,
                "seq_region_name": "chr1",
                "seq_region_start": 180,
                "seq_region_end": 250,
                "seq_region_strand": 1,
            },
            {
                "gene_id": 11,
                "seq_region_name": "chr1",
                "seq_region_start": 800,
                "seq_region_end": 900,
                "seq_region_strand": 1,
            },
        ]
    )

    strict, expanded, mapping = build_loci(core, layer, gap_bp=100)

    # chr1:+ should yield two strict loci: [100,250] and [500,600] and [800,900] => actually 3
    assert len(strict) == 3
    # expanded with gap 100 keeps the 3 since 600->800 distance is 199
    assert len(expanded) == 3

    # counts
    first = strict.sort_values(["locus_start"]).iloc[0]
    assert first.core_gene_count == 1
    assert first.layer_gene_count == 1

    # mapping includes all genes
    assert set(mapping.db_kind.unique()) == {"core", "layer"}
    assert mapping.shape[0] == 4

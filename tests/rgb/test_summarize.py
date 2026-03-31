from __future__ import annotations

import pandas as pd

from ensembl.genes.rgb.summary import summarize_loci


def test_summarize_minimal():
    loci = pd.DataFrame(
        [
            {
                "seq_region_name": "chr1",
                "seq_region_strand": 1,
                "locus_start": 100,
                "locus_end": 300,
                "locus_length": 201,
                "core_gene_count": 1,
                "layer_gene_count": 1,
            }
        ]
    )
    loci["locus_id"] = "chr1:1:100:300:0"

    gene_map = pd.DataFrame(
        [
            {
                "db_kind": "core",
                "gene_id": 1,
                "seq_region_name": "chr1",
                "seq_region_strand": 1,
                "gene_start": 120,
                "locus_id_strict": "chr1:1:100:300:0",
                "locus_id_expanded": "chr1:1:100:300:0",
            },
            {
                "db_kind": "layer",
                "gene_id": 10,
                "seq_region_name": "chr1",
                "seq_region_strand": 1,
                "gene_start": 180,
                "locus_id_strict": "chr1:1:100:300:0",
                "locus_id_expanded": "chr1:1:100:300:0",
            },
        ]
    )

    core_genes = pd.DataFrame(
        [
            {
                "gene_id": 1,
                "seq_region_name": "chr1",
                "seq_region_start": 120,
                "seq_region_end": 200,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "canonical_transcript_id": 1001,
                "logic_name": "core_pipeline",
            }
        ]
    )
    layer_genes = pd.DataFrame(
        [
            {
                "gene_id": 10,
                "seq_region_name": "chr1",
                "seq_region_start": 180,
                "seq_region_end": 280,
                "seq_region_strand": 1,
                "biotype": "",
                "canonical_transcript_id": None,
                "logic_name": "rnaseq_cufflinks",
            }
        ]
    )
    core_tx = pd.DataFrame(
        [
            {
                "transcript_id": 1001,
                "gene_id": 1,
                "seq_region_name": "chr1",
                "seq_region_start": 120,
                "seq_region_end": 200,
                "seq_region_strand": 1,
                "biotype": "protein_coding",
                "logic_name": "core_pipeline",
            }
        ]
    )
    layer_tx = pd.DataFrame(
        [
            {
                "transcript_id": 2001,
                "gene_id": 10,
                "seq_region_name": "chr1",
                "seq_region_start": 180,
                "seq_region_end": 280,
                "seq_region_strand": 1,
                "biotype": "",
                "logic_name": "rnaseq_cufflinks",
            }
        ]
    )

    summary = summarize_loci(
        loci_df=loci,
        gene_map_df=gene_map,
        core_genes=core_genes,
        layer_genes=layer_genes,
        core_tx=core_tx,
        layer_tx=layer_tx,
        evidence_map={"rnaseq_cufflinks": "transcriptomic"},
        locus_gap_bp=5000,
    )
    assert summary.shape[0] == 1
    row = summary.iloc[0]
    assert row.layer_span_bp >= 100
    assert row.core_span_bp > 0
    assert row.layer_bp_covered_by_core > 0
    assert row.layer_logic_name_count == 1
    assert row.layer_evidence_class_count == 1

"""
Pairwise comparison of a query annotation against a reference annotation.

selection   evaluation modes, biotype filters, transcript selection and region
            subsetting applied to comparison-schema DataFrames
classify    locus pairing and exon/CDS/intron-chain classification
summary     headline metrics, stratified counts and per-transcript labels

All functions are pure: they take DataFrames produced by
parsers.annotation.parse_annotation_for_comparison and return new objects.
"""

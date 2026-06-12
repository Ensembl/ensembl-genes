#!/usr/bin/env python3
"""Smoke-test for the changelog comparison module.

Creates two synthetic YAML files (old and new) with known differences,
runs the comparison, and checks the output.

Usage:
    python test_changelog.py
"""
import os
import sys
import tempfile
import importlib.util

# Load changelog.py directly by path to avoid namespace package issues
_CHANGELOG_PATH = os.path.join(os.path.dirname(__file__), "changelog.py")
_spec = importlib.util.spec_from_file_location("changelog", _CHANGELOG_PATH)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

load_yaml_as_keyed_dict = _mod.load_yaml_as_keyed_dict
compare_yamls = _mod.compare_yamls
format_changelog = _mod.format_changelog
write_changelog_tsv = _mod.write_changelog_tsv
_normalise_value = _mod._normalise_value
_extract_key = _mod._extract_key

OLD_YAML = """\
- species: Drosophila melanogaster
  accession: GCA_000001
  annotation_gtf: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Drosophila_melanogaster/GCA_000001/ensembl/geneset/2025_09/genes.gtf.gz
  annotation_gff3: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Drosophila_melanogaster/GCA_000001/ensembl/geneset/2025_09/genes.gff3.gz
  image: Arthropods.png
  ftp_dumps: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Drosophila_melanogaster/GCA_000001/
  beta_link: Coming soon!

- species: Bombus terrestris
  accession: GCA_000002
  annotation_gtf: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Bombus_terrestris/GCA_000002/ensembl/geneset/2025_08/genes.gtf.gz
  image: Hymenoptera.png

- species: Papilio machaon
  accession: GCA_000003
  annotation_gtf: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Papilio_machaon/GCA_000003/ensembl/geneset/2025_07/genes.gtf.gz
  image: Lepidoptera.png
  busco_score: "92.3"
"""

NEW_YAML = """\
- species: Drosophila melanogaster
  accession: GCA_000001
  annotation_gtf: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Drosophila_melanogaster/GCA_000001/ensembl/geneset/2025_10/genes.gtf.gz
  annotation_gff3: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Drosophila_melanogaster/GCA_000001/ensembl/geneset/2025_10/genes.gff3.gz
  image: Arthropods.png
  ftp_dumps: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Drosophila_melanogaster/GCA_000001/
  beta_link: https://beta.ensembl.org/species/some-uuid

- species: Papilio machaon
  accession: GCA_000003
  annotation_gtf: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Papilio_machaon/GCA_000003/ensembl/geneset/2025_07/genes.gtf.gz
  image: Arthropods.png
  busco_score: "93.1"

- species: Apis mellifera
  accession: GCA_000004
  annotation_gtf: https://ftp.ebi.ac.uk/pub/ensemblorganisms/Apis_mellifera/GCA_000004/ensembl/geneset/2025_10/genes.gtf.gz
  image: Hymenoptera.png
"""


def test_normalise_value():
    assert _normalise_value(None) == ""
    assert _normalise_value("  foo/  ") == "foo"
    assert _normalise_value("2025-09") == "2025_09"
    assert _normalise_value("2025_09") == "2025_09"
    print("  [PASS] _normalise_value")


def test_extract_key():
    assert _extract_key({"accession": "GCA_1"}) == "GCA_1"
    assert _extract_key({"assembly_accession": "GCA_2"}) == "GCA_2"
    assert _extract_key({"species": "Homo sapiens"}) == "Homo sapiens"
    assert _extract_key({}) is None
    print("  [PASS] _extract_key")


def test_compare():
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(OLD_YAML)
        old_path = f.name
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(NEW_YAML)
        new_path = f.name

    try:
        old_docs = load_yaml_as_keyed_dict(old_path)
        new_docs = load_yaml_as_keyed_dict(new_path)

        assert "GCA_000001" in old_docs, f"Missing GCA_000001 in old: {old_docs.keys()}"
        assert "GCA_000002" in old_docs
        assert "GCA_000003" in old_docs
        assert len(old_docs) == 3

        assert "GCA_000001" in new_docs
        assert "GCA_000004" in new_docs
        assert "GCA_000002" not in new_docs
        assert len(new_docs) == 3

        added, removed, modified = compare_yamls(old_docs, new_docs)

        # Added: GCA_000004
        assert len(added) == 1, f"Expected 1 added, got {len(added)}: {added}"
        assert added[0][0] == "GCA_000004"
        print("  [PASS] Added detection")

        # Removed: GCA_000002
        assert len(removed) == 1, f"Expected 1 removed, got {len(removed)}: {removed}"
        assert removed[0][0] == "GCA_000002"
        print("  [PASS] Removed detection")

        # Modified: GCA_000001 and GCA_000003
        mod_keys = [m[0] for m in modified]
        assert "GCA_000001" in mod_keys, f"GCA_000001 not in modified: {mod_keys}"
        assert "GCA_000003" in mod_keys, f"GCA_000003 not in modified: {mod_keys}"
        print("  [PASS] Modified detection")

        # Check field-level diffs for GCA_000001
        gca1_mod = [m for m in modified if m[0] == "GCA_000001"][0]
        diff_fields = {d[0] for d in gca1_mod[2]}
        assert (
            "annotation_gtf" in diff_fields
        ), f"annotation_gtf not in diffs: {diff_fields}"
        assert "beta_link" in diff_fields, f"beta_link not in diffs: {diff_fields}"
        # ftp_dumps has trailing slash normalisation so should NOT differ
        assert (
            "ftp_dumps" not in diff_fields
        ), f"ftp_dumps should not differ: {diff_fields}"
        print(
            "  [PASS] Field-level diff (GCA_000001: gtf date change, beta_link transition)"
        )

        # GCA_000003: image change + busco_score change
        gca3_mod = [m for m in modified if m[0] == "GCA_000003"][0]
        diff_fields_3 = {d[0] for d in gca3_mod[2]}
        assert "image" in diff_fields_3, f"image not in diffs: {diff_fields_3}"
        assert (
            "busco_score" in diff_fields_3
        ), f"busco_score not in diffs: {diff_fields_3}"
        print(
            "  [PASS] Field-level diff (GCA_000003: image change, busco_score change)"
        )

        # Format the report and check it's non-empty
        report = format_changelog(added, removed, modified)
        assert "=== CHANGELOG SUMMARY ===" in report
        assert "Added (1)" in report
        assert "Removed (1)" in report
        assert "Modified (2)" in report
        print("  [PASS] Report formatting")

        # TSV output
        tsv_path = new_path + ".changelog.tsv"
        write_changelog_tsv(tsv_path, added, removed, modified)
        with open(tsv_path) as tf:
            lines = tf.readlines()
        assert lines[0].startswith("status\t"), f"Bad TSV header: {lines[0]}"
        assert any("added" in l for l in lines)
        assert any("removed" in l for l in lines)
        assert any("modified" in l for l in lines)
        print("  [PASS] TSV output")
        os.unlink(tsv_path)

    finally:
        os.unlink(old_path)
        os.unlink(new_path)


def test_no_differences():
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(OLD_YAML)
        path = f.name
    try:
        docs = load_yaml_as_keyed_dict(path)
        added, removed, modified = compare_yamls(docs, docs)
        assert len(added) == 0
        assert len(removed) == 0
        assert len(modified) == 0
        report = format_changelog(added, removed, modified)
        assert "No differences detected" in report
        print("  [PASS] No differences (identical inputs)")
    finally:
        os.unlink(path)


def test_date_normalisation_ignored():
    """Verify that YYYY-MM vs YYYY_MM date-only differences in URLs are not flagged."""
    yaml_a = """\
- accession: GCA_TEST
  species: Test species
  annotation_gtf: https://example.com/2025-09/genes.gtf.gz
"""
    yaml_b = """\
- accession: GCA_TEST
  species: Test species
  annotation_gtf: https://example.com/2025_09/genes.gtf.gz
"""
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(yaml_a)
        path_a = f.name
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(yaml_b)
        path_b = f.name
    try:
        docs_a = load_yaml_as_keyed_dict(path_a)
        docs_b = load_yaml_as_keyed_dict(path_b)
        _, _, modified = compare_yamls(docs_a, docs_b)
        assert len(modified) == 0, f"Date normalisation failed, got diffs: {modified}"
        print("  [PASS] Date normalisation (YYYY-MM vs YYYY_MM)")
    finally:
        os.unlink(path_a)
        os.unlink(path_b)


def test_hprc_schema():
    """Verify that HPRC docs (assembly_accession instead of accession) are handled."""
    yaml_old = """\
- assembly: HG002_hap1_hprc_r2
  assembly_accession: GCA_018852605.3
  annotation_gtf: old_url
"""
    yaml_new = """\
- assembly: HG002_hap1_hprc_r2
  assembly_accession: GCA_018852605.3
  annotation_gtf: new_url
"""
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(yaml_old)
        path_a = f.name
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(yaml_new)
        path_b = f.name
    try:
        docs_a = load_yaml_as_keyed_dict(path_a)
        docs_b = load_yaml_as_keyed_dict(path_b)
        assert "GCA_018852605.3" in docs_a
        _, _, modified = compare_yamls(docs_a, docs_b)
        assert len(modified) == 1
        assert modified[0][0] == "GCA_018852605.3"
        print("  [PASS] HPRC schema (assembly_accession key)")
    finally:
        os.unlink(path_a)
        os.unlink(path_b)


if __name__ == "__main__":
    print("Running changelog tests...\n")
    test_normalise_value()
    test_extract_key()
    test_compare()
    test_no_differences()
    test_date_normalisation_ignored()
    test_hprc_schema()
    print("\nAll tests passed.")

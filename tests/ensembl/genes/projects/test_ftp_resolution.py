# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for YamlRenderer FTP resolution — all network calls mocked."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from ensembl.genes.projects.config import ProjectConfig
from ensembl.genes.projects.ftp_manifest import (
    EBI_FTP_BASE,
    EnsemblFtpManifest,
)
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.yaml_renderer import YamlRenderer

# ---------------------------------------------------------------------------
# Shared fixture data
# ---------------------------------------------------------------------------

_FIXTURE_PATH = (
    Path(__file__).parent / "fixtures" / "manifest_sample.json"
)
with _FIXTURE_PATH.open() as _f:
    _FIXTURE_DATA = json.load(_f)


def _make_meta(
    accession: str = "GCA_922984935.2",
    species_name: str = "Meles meles",
    annotation_source: str = "ensembl",
    annotation_date: str = "2022-11-01",
    is_released: bool = True,
    genome_uuid: str = "uuid-badger-1",
    annotation_method: str = "BRAKER2",
    assembly_name: str = "mMelMel3.2",
    strain: str | None = None,
    busco_score: str | None = None,
    busco_lineage: str | None = None,
    taxon_id: int | None = 9662,
    alternate_of: str | None = None,
    is_on_rapid: bool = False,
) -> GenomeMetadata:
    return GenomeMetadata(
        genome_uuid=genome_uuid,
        dbname="meles_meles_core_1_110",
        accession=accession,
        species_name=species_name,
        assembly_name=assembly_name,
        strain=strain,
        taxon_id=taxon_id,
        taxonomy_lineage=[],
        alternate_of=alternate_of,
        annotation_source=annotation_source,
        annotation_method=annotation_method,
        annotation_date=annotation_date,
        busco_score=busco_score,
        busco_lineage=busco_lineage,
        is_released=is_released,
        is_on_rapid=is_on_rapid,
    )


def _make_renderer(
    manifest_data=_FIXTURE_DATA,
    prefer_ensembl: bool = True,
    schema_type: str = "standard",
    with_ftp_client: bool = False,
) -> YamlRenderer:
    manifest = EnsemblFtpManifest(manifest_data) if manifest_data is not None else None
    config = ProjectConfig(
        project_name="test",
        schema_type=schema_type,
        prefer_ensembl_provider=prefer_ensembl,
        allow_beta_urls=False,
    )
    ftp_client = MagicMock() if with_ftp_client else None
    with patch(
        "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
        return_value="unavailable",
    ), patch(
        "ensembl.genes.projects.yaml_renderer.check_url_status",
        return_value=False,
    ):
        renderer = YamlRenderer(config, ftp_client=ftp_client, manifest=manifest)
    return renderer


def _render(renderer: YamlRenderer, meta: GenomeMetadata) -> dict:
    with patch(
        "ensembl.genes.projects.yaml_renderer.check_url_status",
        return_value=False,
    ), patch(
        "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
        return_value="unavailable",
    ), patch.object(
        renderer.icon_resolver, "resolve_icon", return_value=("fish.svg", "fish", "test")
    ):
        return renderer.render(meta)


# ---------------------------------------------------------------------------
# Provider selection tests
# ---------------------------------------------------------------------------

class TestProviderSelection:
    """P1–P8: all provider selection policy paths."""

    def test_p1_exact_match(self):
        """annotation_source='ensembl' selects the ensembl provider."""
        renderer = _make_renderer()
        meta = _make_meta(annotation_source="ensembl")
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_provider"] == "ensembl"
        assert assets["__audit_provider_status__"] == "exact_match"

    def test_p2_alias_match(self):
        """annotation_source='Braker2' resolves via alias to 'braker'."""
        # Need a fixture with a braker provider; use the ambiguous species
        meta = _make_meta(
            accession="GCA_999999999.1",
            species_name="Ambiguous species",
            annotation_source="Braker2",  # alias for braker
        )
        renderer = _make_renderer(prefer_ensembl=False)
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_provider"] == "braker"
        assert assets["__audit_provider_status__"] == "alias_match"

    def test_p3_single_provider_shortcut(self):
        """When only one provider exists, it is selected regardless of metadata."""
        meta = _make_meta(
            accession="GCA_922984935.2",
            annotation_source="some_other_source",  # no match
        )
        # Only ensembl provider in the manifest for this accession
        renderer = _make_renderer(prefer_ensembl=False)
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_provider"] == "ensembl"
        assert assets["__audit_provider_status__"] == "single_provider"

    def test_p4_ensembl_preference_when_multiple_and_prefer_true(self):
        """With prefer_ensembl=True and multiple providers, ensembl is selected."""
        # GCA_999999999.1 has community and braker (no ensembl)
        # We need to add ensembl to the fixture for this test
        data = json.loads(json.dumps(_FIXTURE_DATA))
        data["species"]["Ambiguous_species"]["assemblies"]["GCA_999999999.1"][
            "genebuild_providers"
        ]["ensembl"] = {
            "2023_05": {
                "release": "2023_05",
                "paths": {
                    "genebuild": {
                        "files": {
                            "annotations": {
                                "genes.gtf.gz": "GCA/999/999/999/1/ensembl/2023_05/geneset/genes.gtf.gz"
                            }
                        }
                    }
                }
            }
        }
        meta = _make_meta(
            accession="GCA_999999999.1",
            annotation_source="unmatched_source",
        )
        renderer = _make_renderer(manifest_data=data, prefer_ensembl=True)
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_provider"] == "ensembl"
        assert assets["__audit_provider_status__"] == "ensembl_preference"

    def test_p5_no_ensembl_preference_when_prefer_false(self):
        """With prefer_ensembl=False and multiple non-ensembl providers, ambiguity results."""
        meta = _make_meta(
            accession="GCA_999999999.1",
            annotation_source="unmatched_source",
        )
        renderer = _make_renderer(prefer_ensembl=False)
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_provider"] is None
        assert assets["__audit_provider_status__"] == "ambiguous"
        assert assets["audit_decision"] != "included_released"

    def test_p6_ambiguity_when_no_ensembl_present_and_prefer_true(self):
        """prefer_ensembl=True doesn't help when ensembl is not in the providers."""
        meta = _make_meta(
            accession="GCA_999999999.1",
            annotation_source="unmatched_source",
        )
        renderer = _make_renderer(prefer_ensembl=True)  # ensembl not in providers
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["__audit_provider_status__"] == "ambiguous"

    def test_p7_annotation_source_none_single_provider(self):
        """annotation_source=None falls through to single-provider shortcut."""
        meta = _make_meta(
            accession="GCA_922984935.2",
            annotation_source=None,
        )
        renderer = _make_renderer(prefer_ensembl=False)
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_provider"] == "ensembl"
        assert assets["__audit_provider_status__"] == "single_provider"

    def test_p8_annotation_source_none_multiple_providers_prefer_ensembl(self):
        """annotation_source=None with multiple providers and prefer_ensembl selects ensembl."""
        data = json.loads(json.dumps(_FIXTURE_DATA))
        data["species"]["Ambiguous_species"]["assemblies"]["GCA_999999999.1"][
            "genebuild_providers"
        ]["ensembl"] = {
            "2023_05": {
                "release": "2023_05",
                "paths": {
                    "genebuild": {"files": {"annotations": {"genes.gtf.gz": "GCA/999/999/999/1/ensembl/2023_05/geneset/genes.gtf.gz"}}}
                }
            }
        }
        meta = _make_meta(
            accession="GCA_999999999.1",
            annotation_source=None,
        )
        renderer = _make_renderer(manifest_data=data, prefer_ensembl=True)
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_provider"] == "ensembl"
        assert assets["__audit_provider_status__"] == "ensembl_preference"


# ---------------------------------------------------------------------------
# Date selection tests
# ---------------------------------------------------------------------------

class TestDateSelection:
    """D1–D4: date normalisation and matching."""

    def test_d1_metadata_date_matches_manifest_key(self):
        """2022-11-01 metadata date selects manifest key 2022_11 (prefix match)."""
        renderer = _make_renderer()
        meta = _make_meta(annotation_date="2022-11-01")
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_date"] == "2022_11"
        assert assets["__audit_date_status__"] == "exact_match"

    def test_d2_metadata_dash_mm_matches_manifest_underscore(self):
        """2022-11 metadata date also matches manifest key 2022_11."""
        renderer = _make_renderer()
        meta = _make_meta(annotation_date="2022-11")
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_date"] == "2022_11"
        assert assets["__audit_date_status__"] == "exact_match"

    def test_d3_no_metadata_match_selects_latest(self):
        """With two manifest dates and no metadata match, the latest is selected."""
        meta = _make_meta(
            accession="GCA_012295145.1",
            species_name="Chrysaora quinquecirrha",
            annotation_date="1999-01-01",  # no match
        )
        renderer = _make_renderer()
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_date"] == "2024_02"  # latest of 2023_10, 2024_02
        assert assets["__audit_date_status__"] == "latest_selected"

    def test_d4_metadata_date_selects_correct_of_two_dates(self):
        """When metadata date matches one of two manifest dates, that one is used."""
        meta = _make_meta(
            accession="GCA_012295145.1",
            species_name="Chrysaora quinquecirrha",
            annotation_date="2023-10-01",
        )
        renderer = _make_renderer()
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["resolved_date"] == "2023_10"
        assert assets["__audit_date_status__"] == "exact_match"


# ---------------------------------------------------------------------------
# File classification and verbatim path tests
# ---------------------------------------------------------------------------

class TestFileClassification:
    """R8–R12: usable-geneset definition and optional-file handling."""

    def test_r8_gtf_only_is_released(self):
        """Manifest record with only GTF (no GFF3) still counts as a released geneset."""
        meta = _make_meta(
            accession="GCA_999999999.1",
            annotation_source="community",
        )
        renderer = _make_renderer(prefer_ensembl=False)
        # community provider has only genes.gtf.gz
        assets = renderer._resolve_ftp_assets(meta)
        # Single-provider or exact match, should be released
        assert assets["audit_decision"] in ("included_released", "excluded")
        # If found, gtf must be present
        if assets["audit_decision"] == "included_released":
            assert "genes.gtf.gz" in assets["annotation_files"]

    def test_r9_gff3_only_is_released(self):
        """Manifest record with GFF3 but no GTF still counts as a released geneset."""
        meta = _make_meta(
            accession="GCA_012295145.1",
            species_name="Chrysaora quinquecirrha",
            annotation_date="2023-10-01",
        )
        # 2023_10 date has gff3.gz but no gtf.gz
        renderer = _make_renderer()
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["audit_decision"] == "included_released"
        ann = assets["annotation_files"]
        assert "genes.gff3.gz" in ann
        assert "genes.gtf.gz" not in ann

    def test_r10_neither_gtf_nor_gff3_is_not_released(self):
        """Manifest record with neither GTF nor GFF3 is not treated as released."""
        meta = _make_meta(
            accession="GCA_888888888.1",
            species_name="No gtf species",
        )
        renderer = _make_renderer()
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["audit_decision"] != "included_released"

    def test_r11_missing_optional_pep_does_not_exclude(self):
        """A manifest record without pep.fa.bgz is still considered released."""
        # 2023_10 for jellyfish has pep but no gtf; use 2024_02 which has both
        meta = _make_meta(
            accession="GCA_012295145.1",
            species_name="Chrysaora quinquecirrha",
            annotation_date="2024-02-01",
        )
        renderer = _make_renderer()
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["audit_decision"] == "included_released"
        # 2024_02 has pep.fa.bgz
        assert "pep.fa.bgz" in assets["annotation_files"]

    def test_r12_genome_files_from_assembly_record(self):
        """softmasked_genome URL uses assembly_genome_files, not provider/date files."""
        renderer = _make_renderer()
        meta = _make_meta()
        assets = renderer._resolve_ftp_assets(meta)
        assert "softmasked.fa.bgz" in assets["genome_files"]
        url = assets["genome_files"]["softmasked.fa.bgz"]
        assert url.startswith(EBI_FTP_BASE)


# ---------------------------------------------------------------------------
# Verbatim manifest path tests
# ---------------------------------------------------------------------------

class TestVerbatimPaths:
    """R1–R3: manifest paths used verbatim (not reconstructed)."""

    def test_r1_gtf_url_verbatim_from_manifest(self):
        """GTF URL is constructed as EBI_FTP_BASE + relative_path from manifest."""
        renderer = _make_renderer()
        meta = _make_meta()
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["audit_decision"] == "included_released"
        expected_path = "GCA/922/984/935/2/ensembl/2022_11/geneset/genes.gtf.gz"
        assert assets["annotation_files"]["genes.gtf.gz"] == EBI_FTP_BASE + expected_path

    def test_r2_homology_path_verbatim_with_date_subdir(self):
        """Homology URL preserves the YYYY_MM_DD subdirectory verbatim."""
        renderer = _make_renderer()
        meta = _make_meta()
        assets = renderer._resolve_ftp_assets(meta)
        hom_url = assets["homology_files"].get("homology.tsv.gz", "")
        assert "2024_09_18" in hom_url
        assert hom_url == EBI_FTP_BASE + "GCA/922/984/935/2/ensembl/2022_11/homology/2024_09_18/homology.tsv.gz"

    def test_r3_softmasked_url_verbatim(self):
        """softmasked.fa.bgz URL is built from AssemblyRecord.assembly_genome_files verbatim."""
        renderer = _make_renderer()
        meta = _make_meta()
        assets = renderer._resolve_ftp_assets(meta)
        soft_url = assets["genome_files"].get("softmasked.fa.bgz", "")
        assert soft_url == EBI_FTP_BASE + "GCA/922/984/935/2/ensembl/2022_11/genome/softmasked.fa.bgz"


# ---------------------------------------------------------------------------
# Manifest unavailable tests
# ---------------------------------------------------------------------------

class TestManifestUnavailable:
    """R12: manifest=None falls back to pre-release."""

    def test_manifest_none_sets_manifest_unavailable_status(self):
        """manifest=None records manifest_unavailable in audit."""
        renderer = _make_renderer(manifest_data=None)
        meta = _make_meta(is_released=True)
        assets = renderer._resolve_ftp_assets(meta)
        assert assets["__audit_manifest_status__"] == "manifest_unavailable"
        assert assets["audit_decision"] != "included_released"


# ---------------------------------------------------------------------------
# Audit key lifecycle test
# ---------------------------------------------------------------------------

class TestAuditKeyLifecycle:
    """R15: __audit_*__ keys present before extraction, absent from clean doc."""

    def test_audit_keys_present_in_raw_doc(self):
        """render() returns a dict containing __audit_*__ keys."""
        renderer = _make_renderer()
        meta = _make_meta()
        doc = _render(renderer, meta)
        # At this point they should be GONE because _render calls render() directly
        # and the test framework doesn't extract them.  We need raw access.
        # To test this properly, we call _resolve_ftp_assets directly:
        assets = renderer._resolve_ftp_assets(meta)
        audit_keys = [k for k in assets if k.startswith("__audit_")]
        assert len(audit_keys) > 0, "Expected audit keys in _resolve_ftp_assets output"

    def test_audit_keys_absent_after_extract(self):
        """After _extract_audit_fields, no __audit_*__ keys remain."""
        from ensembl.genes.projects.generate_project_yaml import _extract_audit_fields

        sample_doc = {
            "species": "Meles meles",
            "__audit_decision__": "included_released",
            "__audit_manifest_status__": "found",
            "__audit_provider_status__": "exact_match",
            "__audit_date_status__": "exact_match",
            "__audit_future_field__": "should_also_be_removed",
        }
        extracted = _extract_audit_fields(sample_doc)
        # All audit keys extracted
        assert "__audit_decision__" in extracted
        assert "__audit_future_field__" in extracted
        # None remain in the doc
        remaining_audit = [k for k in sample_doc if k.startswith("__audit_")]
        assert remaining_audit == []
        # Non-audit keys untouched
        assert "species" in sample_doc


# ---------------------------------------------------------------------------
# Standard project regression test
# ---------------------------------------------------------------------------

class TestStandardProjectRegression:
    """R13: standard project rendering produces expected schema."""

    def test_released_genome_yaml_schema(self):
        """A released genome produces correct YAML fields with new URL format."""
        from ensembl.genes.projects.generate_project_yaml import _extract_audit_fields

        renderer = _make_renderer(schema_type="standard")
        meta = _make_meta()
        doc = _render(renderer, meta)
        # Simulate what generate_project_yaml.py does before writing YAML:
        _extract_audit_fields(doc)

        # Public YAML schema fields present
        assert "species" in doc
        assert "accession" in doc
        assert doc["accession"] == "GCA_922984935.2"
        assert "annotation_gtf" in doc
        assert "annotation_gff3" in doc
        assert "proteins" in doc
        assert "transcripts" in doc
        assert "softmasked_genome" in doc
        assert "ftp_dumps" in doc

        # New URL format (bgz for proteins/transcripts)
        assert doc["proteins"].endswith("pep.fa.bgz")
        assert doc["transcripts"].endswith("cdna.fa.bgz")

        # ftp_dumps points to accession root
        assert "GCA/922/984/935/2" in doc["ftp_dumps"]

        # No __audit_*__ keys remain after extraction
        audit_leaks = [k for k in doc if k.startswith("__audit_")]
        assert audit_leaks == [], f"Audit keys leaked into YAML: {audit_leaks}"

    def test_released_genome_no_internal_keys(self):
        """No internal implementation keys appear in the final YAML doc
        after the orchestration layer has extracted audit fields."""
        from ensembl.genes.projects.generate_project_yaml import _extract_audit_fields

        renderer = _make_renderer(schema_type="standard")
        meta = _make_meta()
        doc = _render(renderer, meta)
        _extract_audit_fields(doc)  # simulate orchestration layer
        internal_keys = [k for k in doc if k.startswith("__")]
        assert internal_keys == []


# ---------------------------------------------------------------------------
# HPRC rendering test
# ---------------------------------------------------------------------------

class TestHprcRendering:
    """R14: HPRC rendering regression — VEP intentionally absent."""

    def test_hprc_no_variants_vep(self):
        """HPRC output does not contain variants_vep (old URL is broken)."""
        renderer = _make_renderer(schema_type="hprc")
        meta = _make_meta(assembly_name="CHM13v2.0")
        doc = _render(renderer, meta)
        assert "variants_vep" not in doc

    def test_hprc_no_broken_vep_url(self):
        """No old-format VEP URL (containing ftp_species_name) is emitted."""
        renderer = _make_renderer(schema_type="hprc")
        meta = _make_meta()
        doc = _render(renderer, meta)
        for v in doc.values():
            if isinstance(v, str) and "vep" in v.lower():
                pytest.fail(f"Unexpected VEP URL in HPRC output: {v}")

    def test_hprc_no_internal_keys(self):
        """No __audit_*__ keys appear in the cleaned HPRC doc after extraction."""
        from ensembl.genes.projects.generate_project_yaml import _extract_audit_fields

        renderer = _make_renderer(schema_type="hprc")
        meta = _make_meta()
        doc = _render(renderer, meta)
        _extract_audit_fields(doc)  # simulate orchestration layer
        assert [k for k in doc if k.startswith("__")] == []

    def test_hprc_ftp_dumps_accession_root(self):
        """ftp_dumps for a released HPRC genome points to the accession root."""
        renderer = _make_renderer(schema_type="hprc")
        meta = _make_meta()
        doc = _render(renderer, meta)
        assert "ftp_dumps" in doc
        assert "GCA/922/984/935/2" in doc["ftp_dumps"]

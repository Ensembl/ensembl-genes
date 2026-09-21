"""
Tests for HPRC variation VCF and VEP link resolution.

Verifies:
- Standard and irregular/duplicated-date variation VCF paths
- Deterministic selection across multiple variation dates
- Fallback to older date when newest candidate has no VCF
- Network/URL validation for variation_vcf and variants_vep
- Strict independence of variants_vep and variation_vcf (VEP only, VCF only, both, neither)
- Audit status codes and date tracking
- Preservation of containing directory URLs with trailing slashes

All network access is mocked.
"""

# pylint: disable=missing-class-docstring,missing-function-docstring
# pylint: disable=protected-access,too-many-arguments,too-many-positional-arguments
# pylint: disable=too-few-public-methods,line-too-long

import json
from pathlib import Path
from unittest.mock import patch

from ensembl.genes.projects.config import ProjectConfig, get_project_config
from ensembl.genes.projects.ftp_manifest import EnsemblFtpManifest
from ensembl.genes.projects.generate_project_yaml import _extract_audit_fields
from ensembl.genes.projects.legacy_vep_manifest import LegacyVepManifest
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.yaml_renderer import YamlRenderer

# ---------------------------------------------------------------------------
# Fixture loading
# ---------------------------------------------------------------------------

_FIXTURES_DIR = Path(__file__).parent / "fixtures"

with (_FIXTURES_DIR / "manifest_sample.json").open() as _f:
    _MANIFEST_DATA = json.load(_f)

with (_FIXTURES_DIR / "legacy_species_vep_sample.json").open() as _f:
    _LEGACY_VEP_DATA = json.load(_f)


def _make_meta(
    accession: str,
    species_name: str = "Homo sapiens",
    assembly_name: str = "Test_Assembly",
    annotation_source: str = "ensembl",
    annotation_date: str = "2025-08",
    is_released: bool = True,
) -> GenomeMetadata:
    return GenomeMetadata(
        genome_uuid=f"uuid-{accession.lower().replace('.', '_')}",
        dbname="homo_sapiens_core",
        accession=accession,
        species_name=species_name,
        assembly_name=assembly_name,
        annotation_source=annotation_source,
        annotation_date=annotation_date,
        is_released=is_released,
        taxon_id=9606,
        taxonomy_lineage=[],
    )


def _make_hprc_renderer(
    manifest_data: dict | None = _MANIFEST_DATA,
    legacy_vep_data: dict | None = _LEGACY_VEP_DATA,
    *,
    include_variation_vcf: bool = True,
    use_legacy_vep_fallback: bool = True,
) -> YamlRenderer:
    config = get_project_config("hprc")
    config.include_variation_vcf = include_variation_vcf
    config.use_legacy_vep_fallback = use_legacy_vep_fallback

    manifest = EnsemblFtpManifest(manifest_data) if manifest_data else None
    legacy = LegacyVepManifest(legacy_vep_data) if legacy_vep_data else None

    return YamlRenderer(
        config=config,
        ftp_client=None,
        manifest=manifest,
        legacy_vep_manifest=legacy,
    )


# ---------------------------------------------------------------------------
# Scenarios 1 & 2: Exact nested path and parent-only exclusion
# ---------------------------------------------------------------------------


class TestScenario1And2NestedVcfPath:
    """Scenarios 1 & 2: exact nested VCF path and parent-only path rejection."""

    def test_scenario_1_exact_nested_vcf_path(self):
        """Scenario 1: Duplicated-date components in manifest path are preserved exactly."""
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018469665.2")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        expected_dir = (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "GCA/018/469/665/2/ensembl/2025_08/variation/2025_12_15/2025_12_15/"
        )
        assert doc["variation_vcf"] == expected_dir
        assert doc["variation_vcf"].endswith("/")
        assert not doc["variation_vcf"].endswith("variation.vcf.gz")
        assert "/2025_12_15/2025_12_15/" in doc["variation_vcf"]

    def test_scenario_1_second_nested_accession(self):
        """Scenario 1 (variant): GCA_018506975.2 also preserves duplicated date path."""
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018506975.2")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        expected_dir = (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "GCA/018/506/975/2/ensembl/2025_08/variation/2025_12_15/2025_12_15/"
        )
        assert doc["variation_vcf"] == expected_dir
        assert doc["variation_vcf"].endswith("/")

    def test_scenario_2_parent_only_path_not_emitted(self):
        """Scenario 2: Explicitly assert YAML does NOT contain parent-only URL when file is deeper."""
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018469665.2")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        truncated_parent = (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "GCA/018/469/665/2/ensembl/2025_08/variation/2025_12_15/"
        )
        assert doc["variation_vcf"] != truncated_parent
        assert not doc["variation_vcf"].endswith("/variation/2025_12_15/")


# ---------------------------------------------------------------------------
# Scenarios 3 & 4: VCF file validation vs directory validation
# ---------------------------------------------------------------------------


class TestScenario3And4VcfFileValidation:
    """Scenarios 3 & 4: validation of actual variation.vcf.gz FILE URL."""

    def test_scenario_3_vcf_file_validation_failure_omits_field(self):
        """Scenario 3: Parent directory returns HTTP 200 but variation.vcf.gz fails.

        Expected: variation_vcf omitted, genome retained, audit url_unavailable.
        """
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018469665.2")

        def _mock_check(url: str) -> bool:
            # File URL fails (404), but parent directory would succeed (200)
            if url.endswith(".vcf.gz"):
                return False
            return True

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            side_effect=_mock_check,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variation_vcf" not in doc
        assert doc["assembly_accession"] == "GCA_018469665.2"
        assert audit.get("__audit_decision__") == "included_released"
        assert audit.get("__audit_variation_status__") == "url_unavailable"

    def test_scenario_4_valid_vcf_file_emits_directory(self):
        """Scenario 4: variation.vcf.gz validation succeeds → directory URL emitted."""
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018469665.2")

        captured_urls: list[str] = []

        def _mock_check(url: str) -> bool:
            captured_urls.append(url)
            return True

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            side_effect=_mock_check,
        ):
            doc = renderer.render(meta)

        _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        assert doc["variation_vcf"].endswith("/")
        # Verify check_url_status was passed the actual file URL
        assert any(u.endswith("variation.vcf.gz") for u in captured_urls)


# ---------------------------------------------------------------------------
# Scenario 5: Variation release date modelling
# ---------------------------------------------------------------------------


class TestScenario5VariationDate:
    """Scenario 5: variation release date distinct from annotation date."""

    def test_scenario_5_variation_date_separate_from_annotation_date(self):
        """Scenario 5: Annotation date 2025_08 vs Variation date 2025_12_15.

        audit_variation_date must report 2025_12_15, never 2025_08.
        """
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018469665.2", annotation_date="2025-08")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        assert audit.get("__audit_variation_status__") == "available"
        assert audit.get("__audit_variation_date__") == "2025_12_15"
        assert audit.get("__audit_variation_date__") != "2025_08"


# ---------------------------------------------------------------------------
# Scenarios 6 & 7: Multiple variation dates and fallback
# ---------------------------------------------------------------------------


class TestScenario6And7MultipleVariationDates:
    """Scenarios 6 & 7: deterministic newest selection and fallback on failure."""

    def test_scenario_6_multiple_variation_dates_selects_newest(self):
        """Scenario 6: GRCh38 has 2025_08_13 and 2026_07_13 → selects 2026_07_13."""
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_000001405.29", annotation_date="2026-04")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        assert "2026_07_13/" in doc["variation_vcf"]
        assert audit.get("__audit_variation_status__") == "available"
        assert audit.get("__audit_variation_date__") == "2026_07_13"

    def test_scenario_7_newest_invalid_selects_older_valid(self):
        """Scenario 7: Newest date 2026_07_13 file fails check, older 2025_08_13 succeeds.

        Expected: older valid VCF selected deterministically.
        """
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_000001405.29", annotation_date="2026-04")

        def _mock_check(url: str) -> bool:
            # Newest candidate fails, older candidate succeeds
            if "2026_07_13" in url:
                return False
            return True

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            side_effect=_mock_check,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        assert "2025_08_13/" in doc["variation_vcf"]
        assert "2026_07_13/" not in doc["variation_vcf"]
        assert audit.get("__audit_variation_status__") == "available"
        assert audit.get("__audit_variation_date__") == "2025_08_13"


# ---------------------------------------------------------------------------
# Scenario 8: No variation data
# ---------------------------------------------------------------------------


class TestScenario8NoVcf:
    """Scenario 8: assembly with no variation data in manifest."""

    def test_scenario_8_no_vcf_omits_field_and_retains_genome(self):
        """Scenario 8: GCA_018852605.1 has variation_data: {} in manifest."""
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018852605.1", annotation_date="2022-07")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variation_vcf" not in doc
        assert doc["assembly_accession"] == "GCA_018852605.1"
        assert audit.get("__audit_decision__") == "included_released"
        assert audit.get("__audit_variation_status__") == "no_vcf"


# ---------------------------------------------------------------------------
# Scenarios 9 & 10: VEP available and file validation
# ---------------------------------------------------------------------------


class TestScenario9And10VepResolutionAndValidation:
    """Scenarios 9 & 10: VEP resolution for known HPRC assemblies and file validation."""

    def test_scenario_9_vep_available(self):
        """Scenario 9: Known HPRC example GCA_046332015.1 with genes.gff3.bgz valid (authoritative path rule)."""
        renderer = _make_hprc_renderer(legacy_vep_data=None)
        meta = _make_meta(
            "GCA_046332015.1",
            assembly_name="HG002.alt.mat.v1",
            annotation_date="2025-12",
        )

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variants_vep" in doc
        expected_vep_dir = (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "Homo_sapiens/GCA_046332015.1/vep/ensembl/geneset/2025_12/"
        )
        assert doc["variants_vep"] == expected_vep_dir
        assert doc["variants_vep"].endswith("/")
        assert "genes.gff3.bgz" not in doc["variants_vep"]
        assert audit.get("__audit_vep_status__") == "available_legacy_manifest"

    def test_scenario_9_vep_available_with_manifest_source(self):
        """Scenario 9 (variant): With authoritative LegacyVepManifest source containing GCA_046332015.1."""
        vep_data = {
            "species": {
                "Homo_sapiens": {
                    "assemblies": {
                        "GCA_046332015.1": {
                            "genebuild_providers": {
                                "ensembl": {
                                    "2025_12": {
                                        "paths": {
                                            "genebuild": {
                                                "files": {
                                                    "vep": {
                                                        "genes.gff3.bgz": "Homo_sapiens/GCA_046332015.1/vep/ensembl/geneset/2025_12/genes.gff3.bgz"
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        renderer = _make_hprc_renderer(legacy_vep_data=vep_data)
        meta = _make_meta(
            "GCA_046332015.1",
            assembly_name="HG002.alt.mat.v1",
            annotation_date="2025-12",
        )

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variants_vep" in doc
        expected_vep_dir = (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "Homo_sapiens/GCA_046332015.1/vep/ensembl/geneset/2025_12/"
        )
        assert doc["variants_vep"] == expected_vep_dir
        assert audit.get("__audit_vep_status__") == "available_legacy_manifest"

    def test_scenario_10_vep_parent_exists_but_file_does_not(self):
        """Scenario 10: Parent directory exists but genes.gff3.bgz fails URL check.

        Expected: variants_vep omitted, genome retained.
        """
        renderer = _make_hprc_renderer(legacy_vep_data=None)
        meta = _make_meta(
            "GCA_046332015.1",
            assembly_name="HG002.alt.mat.v1",
            annotation_date="2025-12",
        )

        def _mock_check(url: str) -> bool:
            # genes.gff3.bgz fails, but parent directory would succeed
            if url.endswith(".gff3.bgz"):
                return False
            return True

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            side_effect=_mock_check,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variants_vep" not in doc
        assert doc["assembly_accession"] == "GCA_046332015.1"
        assert audit.get("__audit_decision__") == "included_released"


# ---------------------------------------------------------------------------
# Scenarios 11 to 14: Independence of VEP and variation_vcf
# ---------------------------------------------------------------------------


class TestScenario11To14Independence:
    """Scenarios 11-14: strict independence of variants_vep and variation_vcf."""

    def test_scenario_11_vep_only(self):
        """Scenario 11: Genome has VEP, but no variation VCF."""
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_018852605.1", annotation_date="2022-07")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variants_vep" in doc
        assert "variation_vcf" not in doc
        assert audit.get("__audit_vep_status__") == "available_legacy_manifest"
        assert audit.get("__audit_variation_status__") == "no_vcf"

    def test_scenario_12_vcf_only(self):
        """Scenario 12: Genome has variation VCF, but no VEP."""
        # GCA_000001405.29 has variation in manifest, but is non-HPRC so no VEP probe
        renderer = _make_hprc_renderer()
        meta = _make_meta("GCA_000001405.29", annotation_date="2026-04")

        def _mock_check(url: str) -> bool:
            if "vep" in url:
                return False
            return True

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            side_effect=_mock_check,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        assert "variants_vep" not in doc
        assert audit.get("__audit_variation_status__") == "available"

    def test_scenario_13_both_present(self):
        """Scenario 13: Genome has both VEP and variation VCF present independently."""
        renderer = _make_hprc_renderer(legacy_vep_data=None)
        meta = _make_meta("GCA_018469665.2", annotation_date="2025-08")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variants_vep" in doc
        assert "variation_vcf" in doc
        assert (
            doc["variation_vcf"]
            == "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/018/469/665/2/ensembl/2025_08/variation/2025_12_15/2025_12_15/"
        )
        assert (
            doc["variants_vep"]
            == "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/GCA_018469665.2/vep/ensembl/geneset/2025_08/"
        )
        assert audit.get("__audit_vep_status__") == "available_legacy_manifest"
        assert audit.get("__audit_variation_status__") == "available"

    def test_scenario_14_neither_present(self):
        """Scenario 14: Genome has neither VEP nor VCF; genome is still included."""
        renderer = _make_hprc_renderer()
        meta = _make_meta(
            "GCA_922984935.2",
            species_name="Meles meles",
            annotation_date="2022-11",
        )

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=False,
        ):
            doc = renderer.render(meta)

        audit = _extract_audit_fields(doc)
        assert "variants_vep" not in doc
        assert "variation_vcf" not in doc
        assert audit.get("__audit_decision__") == "included_released"
        assert doc["assembly_accession"] == "GCA_922984935.2"


# ---------------------------------------------------------------------------
# Scenario 15: T2T-style nested variation path regression
# ---------------------------------------------------------------------------


class TestScenario15T2TNestedVariation:
    """Scenario 15: T2T-style extra subdirectories in manifest path are respected."""

    def test_scenario_15_t2t_nested_variation_path(self):
        """Scenario 15: GCA_009914755.4 manifest path has /variation/2025_01_16/nested_extra/.

        Expected: containing directory preserves /nested_extra/ verbatim.
        """
        renderer = _make_hprc_renderer()
        meta = _make_meta(
            "GCA_009914755.4",
            assembly_name="T2T-CHM13v2.0",
            annotation_date="2022-07",
        )

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        _extract_audit_fields(doc)
        assert "variation_vcf" in doc
        expected_dir = (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "GCA/009/914/755/4/ensembl/2022_07/variation/2025_01_16/nested_extra/"
        )
        assert doc["variation_vcf"] == expected_dir
        assert doc["variation_vcf"].endswith("/")
        assert "/nested_extra/" in doc["variation_vcf"]


# ---------------------------------------------------------------------------
# Scenario 16: Non-HPRC project behavior
# ---------------------------------------------------------------------------


class TestScenario16NonHprcVariationConfig:
    """Ensure variation_vcf is omitted for non-HPRC projects by default."""

    def test_variation_vcf_disabled_for_standard_project(self):
        config = ProjectConfig(
            project_name="vgp",
            schema_type="standard",
            include_variation_vcf=False,
        )
        manifest = EnsemblFtpManifest(_MANIFEST_DATA)
        renderer = YamlRenderer(config=config, manifest=manifest)
        meta = _make_meta("GCA_018469665.2")

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            doc = renderer.render(meta)

        _extract_audit_fields(doc)
        assert "variation_vcf" not in doc

"""
Tests for HPRC VEP resolution via the legacy manifest fallback.

All network calls (URL checks, manifest downloads) are mocked.
Run with::

    pytest tests/ensembl/genes/projects/test_hprc_vep_resolution.py -v
"""

# pylint: disable=missing-class-docstring,missing-function-docstring
# pylint: disable=protected-access

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from ensembl.genes.projects.config import ProjectConfig, get_project_config
from ensembl.genes.projects.ftp_manifest import EnsemblFtpManifest
from ensembl.genes.projects.generate_project_yaml import _extract_audit_fields
from ensembl.genes.projects.legacy_vep_manifest import LegacyVepManifest
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.yaml_renderer import YamlRenderer

# ---------------------------------------------------------------------------
# Fixture data
# ---------------------------------------------------------------------------

_LEGACY_FIXTURE = Path(__file__).parent / "fixtures" / "legacy_species_vep_sample.json"

with _LEGACY_FIXTURE.open() as _f:
    _LEGACY_MANIFEST_DATA = json.load(_f)

# A real HPRC accession that appears in the legacy fixture.
_HPRC_ACCESSION = "GCA_018852605.1"


def _make_hprc_meta(
    accession: str = _HPRC_ACCESSION,
    annotation_source: str = "ensembl",
    annotation_date: str = "2022-07",
    is_released: bool = True,
) -> GenomeMetadata:
    return GenomeMetadata(
        genome_uuid="uuid-hg002-1",
        dbname="homo_sapiens_core",
        accession=accession,
        species_name="Homo sapiens",
        assembly_name="HG002.alt.pat.f1_v2",
        annotation_source=annotation_source,
        annotation_date=annotation_date,
        is_released=is_released,
        taxon_id=9606,
        taxonomy_lineage=[],
    )


def _make_hprc_renderer(
    *,
    new_manifest_data=None,
    legacy_manifest_data=None,
) -> YamlRenderer:
    """Build a YamlRenderer for HPRC with controlled manifests."""
    config = get_project_config("hprc")

    new_manifest = (
        EnsemblFtpManifest(new_manifest_data) if new_manifest_data is not None else None
    )
    legacy_manifest = (
        LegacyVepManifest(legacy_manifest_data)
        if legacy_manifest_data is not None
        else None
    )

    return YamlRenderer(
        config,
        ftp_client=None,
        manifest=new_manifest,
        legacy_vep_manifest=legacy_manifest,
    )


# ---------------------------------------------------------------------------
# TestVepResolutionOrder
# ---------------------------------------------------------------------------


class TestVepResolutionOrder:

    def test_vep_not_applicable_for_non_hprc_project(self):
        """Standard project has use_legacy_vep_fallback=False → not_applicable."""
        config = ProjectConfig(
            project_name="dtol",
            schema_type="standard",
            use_legacy_vep_fallback=False,
        )
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        meta = _make_hprc_meta()

        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2022_07"}

        _url, status = renderer._resolve_vep_url(meta, ftp_assets)
        assert status == "not_applicable"
        assert _url is None

    def test_vep_legacy_manifest_unavailable_when_not_loaded(self):
        """HPRC config, but no legacy manifest loaded → legacy_manifest_unavailable."""
        config = get_project_config("hprc")
        renderer = YamlRenderer(config, legacy_vep_manifest=None)
        meta = _make_hprc_meta()
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2022_07"}

        _url, status = renderer._resolve_vep_url(meta, ftp_assets)
        assert status == "legacy_manifest_unavailable"
        assert _url is None

    def test_vep_not_found_when_accession_absent_from_legacy(self):
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        meta = _make_hprc_meta(accession="GCA_000000000.1")
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2022_07"}

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            _url, status = renderer._resolve_vep_url(meta, ftp_assets)
        assert status == "not_found"
        assert _url is None

    def test_vep_not_found_when_accession_has_no_vep(self):
        """GCA_999000001.1 is in fixture but has no vep section."""
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        meta = _make_hprc_meta(accession="GCA_999000001.1")
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2023_01"}

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            _url, status = renderer._resolve_vep_url(meta, ftp_assets)
        assert status == "not_found"
        assert _url is None

    def test_vep_available_legacy_manifest_when_url_resolves(self):
        """GCA_018852605.1 has VEP in legacy fixture; URL check passes.
        The emitted URL must point to the dated release *directory*."""
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        meta = _make_hprc_meta()
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2022_07"}

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            vep_url, status = renderer._resolve_vep_url(meta, ftp_assets)

        assert status == "available_legacy_manifest"
        assert vep_url is not None
        assert vep_url.startswith("https://ftp.ebi.ac.uk/pub/ensemblorganisms/")
        # Must link to the directory, not the file.
        assert vep_url.endswith("/")
        assert "genes.gff3.bgz" not in vep_url
        assert "genes.gff3.bgz.csi" not in vep_url
        assert _HPRC_ACCESSION in vep_url
        # Must contain the exact date component from the manifest.
        assert "2022_07" in vep_url

    def test_vep_url_is_directory_derived_from_manifest_path(self):
        """The VEP URL is the parent directory of the manifest's primary VEP file.

        The manifest stores the relative file path; the emitted URL must be
        the containing dated directory (trailing /) with no filename appended.
        """
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        meta = _make_hprc_meta()
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2022_07"}

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            vep_url, _ = renderer._resolve_vep_url(meta, ftp_assets)

        expected_directory = (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/"
        )
        assert vep_url == expected_directory
        assert vep_url.endswith("/")
        assert "genes.gff3.bgz" not in vep_url

    def test_vep_url_unavailable_when_check_fails(self):
        """URL check returns False → legacy_url_unavailable, no field emitted."""
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        meta = _make_hprc_meta()
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2022_07"}

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=False,
        ):
            vep_url, status = renderer._resolve_vep_url(meta, ftp_assets)

        assert status == "legacy_url_unavailable"
        assert vep_url is None

    def test_ambiguous_legacy_record_produces_no_url(self):
        """GCA_888000001.1 has two providers; no provider hint → ambiguous."""
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        meta = _make_hprc_meta(accession="GCA_888000001.1")
        ftp_assets = {
            "resolved_provider": None,
            "resolved_date": "2022_07",
        }

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            return_value=True,
        ):
            vep_url, status = renderer._resolve_vep_url(meta, ftp_assets)

        assert status == "ambiguous_legacy_record"
        assert vep_url is None

    def test_directory_url_passed_to_check_url_status(self):
        """check_url_status must receive the directory URL (trailing /), not the
        genes.gff3.bgz file URL."""
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)
        meta = _make_hprc_meta()
        ftp_assets = {"resolved_provider": "ensembl", "resolved_date": "2022_07"}

        captured_urls: list[str] = []

        def _capture(url: str) -> bool:
            captured_urls.append(url)
            return True

        with patch(
            "ensembl.genes.projects.yaml_renderer.check_url_status",
            side_effect=_capture,
        ):
            renderer._resolve_vep_url(meta, ftp_assets)

        assert len(captured_urls) == 1
        checked = captured_urls[0]
        assert checked.endswith(
            "/"
        ), f"check_url_status received a non-directory URL: {checked!r}"
        assert (
            "genes.gff3.bgz" not in checked
        ), f"check_url_status received a file URL instead of a directory: {checked!r}"
        assert "2022_07" in checked


# ---------------------------------------------------------------------------
# TestHprcRendererOutput
# ---------------------------------------------------------------------------


# A pre-built ftp_resolution dict representing a successfully resolved HPRC
# assembly.  Used in full-render tests that patch _resolve_ftp_assets so the
# renderer does not short-circuit at 'excluded' before reaching VEP resolution.
_HPRC_FTP_RESOLUTION = {
    "is_released": True,
    "ftp_species_name": "Homo_sapiens",
    "resolved_date": "2022_07",
    "resolved_provider": "ensembl",
    "acc_ftp_path": "GCA/018/852/605/1",
    "audit_decision": "included_released",
    "audit_reason": "Resolved from new manifest.",
    "annotation_files": {
        "genes.gtf.gz": (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "GCA/018/852/605/1/ensembl/2022_07/geneset/genes.gtf.gz"
        ),
        "genes.gff3.gz": (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "GCA/018/852/605/1/ensembl/2022_07/geneset/genes.gff3.gz"
        ),
    },
    "genome_files": {},
    "homology_files": {},
    "variation_files": {},
}


class TestHprcRendererOutput:

    def test_variants_vep_present_when_legacy_resolved(self):
        """Full render: variants_vep emitted when legacy manifest resolves."""
        renderer = _make_hprc_renderer(
            legacy_manifest_data=_LEGACY_MANIFEST_DATA,
        )
        meta = _make_hprc_meta()

        with (
            patch.object(
                renderer, "_resolve_ftp_assets", return_value=_HPRC_FTP_RESOLUTION
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
        ):
            doc = renderer.render(meta)

        assert "variants_vep" in doc
        vep = doc["variants_vep"]
        assert vep.startswith("https://")
        assert vep.endswith(
            "/"
        ), "variants_vep must link to a directory (trailing slash)"
        assert "genes.gff3.bgz" not in vep, "variants_vep must not contain a filename"
        assert "genes.gff3.bgz.csi" not in vep
        assert _HPRC_ACCESSION in vep
        assert "2022_07" in vep
        assert doc["__audit_vep_status__"] == "available_legacy_manifest"

    def test_variants_vep_absent_when_url_check_fails(self):
        """Full render: variants_vep omitted when URL check fails."""
        renderer = _make_hprc_renderer(
            legacy_manifest_data=_LEGACY_MANIFEST_DATA,
        )
        meta = _make_hprc_meta()

        with (
            patch.object(
                renderer, "_resolve_ftp_assets", return_value=_HPRC_FTP_RESOLUTION
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=False,
            ),
        ):
            doc = renderer.render(meta)

        assert "variants_vep" not in doc
        assert doc["__audit_vep_status__"] == "legacy_url_unavailable"

    def test_variants_vep_absent_when_legacy_manifest_none(self):
        """Full render: no legacy manifest → variants_vep absent."""
        renderer = _make_hprc_renderer(legacy_manifest_data=None)
        meta = _make_hprc_meta()

        with (
            patch.object(
                renderer, "_resolve_ftp_assets", return_value=_HPRC_FTP_RESOLUTION
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
        ):
            doc = renderer.render(meta)

        assert "variants_vep" not in doc
        assert doc["__audit_vep_status__"] == "legacy_manifest_unavailable"

    def test_variants_vep_absent_when_accession_not_in_legacy(self):
        """Full render: accession absent from legacy manifest → not_found."""
        renderer = _make_hprc_renderer(
            legacy_manifest_data=_LEGACY_MANIFEST_DATA,
        )
        meta = _make_hprc_meta(accession="GCA_999000001.1")
        no_vep_resolution = {
            **_HPRC_FTP_RESOLUTION,
            "resolved_provider": None,
            "resolved_date": "2023_01",
        }

        with (
            patch.object(
                renderer, "_resolve_ftp_assets", return_value=no_vep_resolution
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
        ):
            doc = renderer.render(meta)

        assert "variants_vep" not in doc
        assert doc["__audit_vep_status__"] == "not_found"

    def test_hprc_genome_not_excluded_when_vep_missing(self):
        """The genome must still appear even if VEP cannot be resolved."""
        renderer = _make_hprc_renderer(legacy_manifest_data=None)
        meta = _make_hprc_meta()

        with (
            patch.object(
                renderer, "_resolve_ftp_assets", return_value=_HPRC_FTP_RESOLUTION
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
        ):
            doc = renderer.render(meta)

        # Core HPRC fields must still be present
        assert "assembly_accession" in doc
        assert doc["assembly_accession"] == _HPRC_ACCESSION
        assert doc["__audit_decision__"] != "excluded"

    def test_no_ftp_species_name_used_for_vep(self):
        """VEP URL must come from the manifest, not be built from ftp_species_name."""
        renderer = _make_hprc_renderer(
            legacy_manifest_data=_LEGACY_MANIFEST_DATA,
        )
        meta = _make_hprc_meta()

        with (
            patch.object(
                renderer, "_resolve_ftp_assets", return_value=_HPRC_FTP_RESOLUTION
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
        ):
            doc = renderer.render(meta)

        vep = doc.get("variants_vep", "")
        # The real URL is derived from the manifest path, which includes the
        # accession.  It must be a directory (trailing slash), with no filename.
        assert _HPRC_ACCESSION in vep, (
            "VEP URL must use the manifest path (which includes the accession), "
            "not a guessed species-name path."
        )
        assert vep.endswith("/"), "variants_vep must link to a directory"
        assert "genes.gff3.bgz" not in vep, "variants_vep must not contain a filename"

    def test_audit_keys_stripped_after_extract(self):
        """All __audit_*__ keys disappear after _extract_audit_fields."""
        renderer = _make_hprc_renderer(
            legacy_manifest_data=_LEGACY_MANIFEST_DATA,
        )
        meta = _make_hprc_meta()

        with (
            patch.object(
                renderer, "_resolve_ftp_assets", return_value=_HPRC_FTP_RESOLUTION
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
        ):
            doc = renderer.render(meta)

        _extract_audit_fields(doc)
        assert not any(k.startswith("__") for k in doc)


# ---------------------------------------------------------------------------
# TestNonHprcUnchanged
# ---------------------------------------------------------------------------


class TestNonHprcUnchanged:

    def test_standard_project_does_not_require_legacy_manifest(self):
        """Standard project config must not have use_legacy_vep_fallback."""
        config = get_project_config("dtol")
        assert config.use_legacy_vep_fallback is False

    def test_mouse_project_does_not_require_legacy_manifest(self):
        config = get_project_config("mouse_genomes")
        assert config.use_legacy_vep_fallback is False

    def test_hprc_project_has_legacy_vep_fallback_enabled(self):
        config = get_project_config("hprc")
        assert config.use_legacy_vep_fallback is True

    def test_standard_render_unaffected_by_legacy_manifest(self):
        """Passing a legacy_vep_manifest to a standard renderer has no effect."""
        import json
        from pathlib import Path

        from ensembl.genes.projects.ftp_manifest import EnsemblFtpManifest

        fixture_path = Path(__file__).parent / "fixtures" / "manifest_sample.json"
        with fixture_path.open() as fh:
            new_data = json.load(fh)

        config = ProjectConfig(
            project_name="dtol",
            schema_type="standard",
            allow_beta_urls=False,
        )
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(
            config,
            manifest=EnsemblFtpManifest(new_data),
            legacy_vep_manifest=legacy,
        )

        # A badger genome from the new manifest fixture
        meta = GenomeMetadata(
            genome_uuid="uuid-badger-1",
            dbname="meles_meles_core_1_110",
            accession="GCA_922984935.2",
            species_name="Meles meles",
            assembly_name="mMelMel3.2",
            annotation_source="ensembl",
            annotation_date="2022-11-01",
            is_released=True,
            taxon_id=9662,
            taxonomy_lineage=[],
        )

        with (
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=False,
            ),
        ):
            doc = renderer.render(meta)

        # Standard schema must not emit variants_vep
        assert "variants_vep" not in doc
        # Standard schema must not have a vep_status audit key
        assert "__audit_vep_status__" not in doc

    def test_renderer_reuses_injected_legacy_manifest(self):
        """The renderer holds and reuses one LegacyVepManifest instance."""
        config = get_project_config("hprc")
        legacy = LegacyVepManifest(_LEGACY_MANIFEST_DATA)
        renderer = YamlRenderer(config, legacy_vep_manifest=legacy)

        assert renderer.legacy_vep_manifest is legacy

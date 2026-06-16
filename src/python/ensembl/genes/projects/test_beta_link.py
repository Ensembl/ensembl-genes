"""Tests for beta.ensembl.org species-page availability checking.

All network calls are mocked. Run with:
    pytest src/python/ensembl/genes/projects/test_beta_link.py -v
"""

from unittest.mock import MagicMock, patch

from ensembl.genes.projects.config import ProjectConfig
from ensembl.genes.projects.ftp_client import check_beta_species_status
from ensembl.genes.projects.models import GenomeMetadata
from ensembl.genes.projects.yaml_renderer import YamlRenderer


def _resp(status_code=200, body=""):
    m = MagicMock()
    m.status_code = status_code
    m.text = body
    return m


def _meta(uuid="525cc33d-82e2-4df5-90f2-70ef460cb418"):
    return GenomeMetadata(
        genome_uuid=uuid,
        dbname="db",
        accession="GCA_982375345.1",
        species_name="Achillea ptarmica",
        assembly_name="asm",
    )


# ---------------------------------------------------------------------------
# check_beta_species_status
# ---------------------------------------------------------------------------


class TestCheckBetaSpeciesStatus:

    def test_available_normal_page(self):
        with patch(
            "ensembl.genes.projects.ftp_client.requests.get",
            return_value=_resp(200, "<html>Achillea ptarmica genome browser</html>"),
        ):
            assert check_beta_species_status("uuid-ok") == "available"

    def test_unavailable_not_recognised_marker(self):
        body = "<html>We do not recognise the species identified as ...</html>"
        with patch(
            "ensembl.genes.projects.ftp_client.requests.get",
            return_value=_resp(200, body),
        ):
            assert check_beta_species_status("uuid-x") == "unavailable"

    def test_unavailable_species_selector_marker(self):
        body = "<html>Find available species in the Species selector</html>"
        with patch(
            "ensembl.genes.projects.ftp_client.requests.get",
            return_value=_resp(200, body),
        ):
            assert check_beta_species_status("uuid-y") == "unavailable"

    def test_unavailable_non_200(self):
        with patch(
            "ensembl.genes.projects.ftp_client.requests.get",
            return_value=_resp(404, "not found"),
        ):
            assert check_beta_species_status("uuid-z") == "unavailable"

    def test_error_on_exception(self):
        import requests as real_requests

        with patch(
            "ensembl.genes.projects.ftp_client.requests.get",
            side_effect=real_requests.Timeout("timeout"),
        ):
            assert check_beta_species_status("uuid-timeout") == "error"

    def test_error_on_missing_uuid(self):
        assert check_beta_species_status("") == "error"
        assert check_beta_species_status("unknown") == "error"


# ---------------------------------------------------------------------------
# Renderer caching
# ---------------------------------------------------------------------------


class TestRendererCaching:

    def test_beta_status_cached_per_uuid(self):
        config = ProjectConfig(project_name="dtol", schema_type="standard")
        renderer = YamlRenderer(config)
        with patch(
            "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
            return_value="available",
        ) as mock_check:
            assert renderer._check_beta_status("uuid-1") == "available"
            assert renderer._check_beta_status("uuid-1") == "available"
            assert renderer._check_beta_status("uuid-1") == "available"
            # Only one underlying network check despite three calls
            assert mock_check.call_count == 1


# ---------------------------------------------------------------------------
# _resolve_beta_link
# ---------------------------------------------------------------------------


class TestResolveBetaLink:

    def _renderer(self):
        return YamlRenderer(ProjectConfig(project_name="dtol"))

    def test_available_emits_real_url(self):
        renderer = self._renderer()
        with patch(
            "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
            return_value="available",
        ):
            link, status = renderer._resolve_beta_link(_meta(), target_released=True)
        assert link == (
            "https://beta.ensembl.org/species/" "525cc33d-82e2-4df5-90f2-70ef460cb418"
        )
        assert status == "available"

    def test_unavailable_emits_coming_soon(self):
        """The Achillea ptarmica acceptance case."""
        renderer = self._renderer()
        with patch(
            "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
            return_value="unavailable",
        ):
            link, status = renderer._resolve_beta_link(_meta(), target_released=True)
        assert link == "Coming soon!"
        assert status == "unavailable"

    def test_error_emits_coming_soon(self):
        renderer = self._renderer()
        with patch(
            "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
            return_value="error",
        ):
            link, status = renderer._resolve_beta_link(_meta(), target_released=True)
        assert link == "Coming soon!"
        assert status == "error"

    def test_prerelease_skips_check(self):
        renderer = self._renderer()
        with patch(
            "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
            return_value="available",
        ) as mock_check:
            link, status = renderer._resolve_beta_link(_meta(), target_released=False)
        assert link == "Coming soon!"
        assert status == "skipped_prerelease"
        # No network check performed for pre-release genomes
        assert mock_check.call_count == 0

    def test_missing_uuid_skips_check(self):
        renderer = self._renderer()
        meta = _meta(uuid="unknown")
        with patch(
            "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
            return_value="available",
        ) as mock_check:
            link, status = renderer._resolve_beta_link(meta, target_released=True)
        assert link == "Coming soon!"
        assert status == "skipped_no_uuid"
        assert mock_check.call_count == 0


# ---------------------------------------------------------------------------
# End-to-end: standard render emits Coming soon! for unavailable beta
# ---------------------------------------------------------------------------


class TestStandardRenderIntegration:

    def test_unavailable_beta_keeps_genome_with_coming_soon(self):
        """Achillea ptarmica: FTP valid, beta unavailable -> Coming soon!,
        genome NOT excluded."""
        config = ProjectConfig(project_name="dtol", schema_type="standard")
        renderer = YamlRenderer(config)
        meta = _meta()
        meta.is_released = True

        # Force the FTP resolution to look released with valid assets, and the
        # beta page to be unavailable.
        ftp_resolution = {
            "is_released": True,
            "ftp_species_name": "Achillea_ptarmica",
            "resolved_date": "2024_01",
            "audit_decision": "included_released",
            "audit_reason": "Found released FTP assets.",
        }
        with (
            patch.object(renderer, "_resolve_ftp_assets", return_value=ftp_resolution),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="unavailable",
            ),
        ):
            doc = renderer.render(meta)

        # Genome retained (has FTP links) but beta_link downgraded
        assert doc["beta_link"] == "Coming soon!"
        assert doc["accession"] == "GCA_982375345.1"
        assert "annotation_gtf" in doc
        assert doc["__audit_beta_status__"] == "unavailable"

    def test_available_beta_emits_real_url(self):
        config = ProjectConfig(project_name="dtol", schema_type="standard")
        renderer = YamlRenderer(config)
        meta = _meta()
        meta.is_released = True

        ftp_resolution = {
            "is_released": True,
            "ftp_species_name": "Achillea_ptarmica",
            "resolved_date": "2024_01",
            "audit_decision": "included_released",
            "audit_reason": "Found released FTP assets.",
        }
        with (
            patch.object(renderer, "_resolve_ftp_assets", return_value=ftp_resolution),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_url_status",
                return_value=True,
            ),
            patch(
                "ensembl.genes.projects.yaml_renderer.check_beta_species_status",
                return_value="available",
            ),
        ):
            doc = renderer.render(meta)

        assert doc["beta_link"] == (
            "https://beta.ensembl.org/species/" "525cc33d-82e2-4df5-90f2-70ef460cb418"
        )
        assert doc["__audit_beta_status__"] == "available"

"""
Tests for legacy_vep_manifest.py — the reader for species.json VEP fallback data.

All network calls are mocked.  Run with::

    pytest tests/ensembl/genes/projects/test_legacy_vep_manifest.py -v
"""

# pylint: disable=missing-class-docstring,missing-function-docstring
# pylint: disable=protected-access

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from ensembl.genes.projects.legacy_vep_manifest import (
    EBI_FTP_BASE,
    LegacyVepManifest,
    LegacyVepManifestError,
    LegacyVepRecord,
    VEP_PRIMARY_FILE,
    _manifest_url,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_FIXTURE_PATH = (
    Path(__file__).parent / "fixtures" / "legacy_species_vep_sample.json"
)


def _load_fixture() -> dict:
    with open(_FIXTURE_PATH, encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# TestManifestUrl
# ---------------------------------------------------------------------------


class TestManifestUrl:
    def test_trailing_slash_on_base(self):
        url = _manifest_url("https://ftp.ebi.ac.uk/pub/ensemblorganisms/", "a/b.bgz")
        assert url == "https://ftp.ebi.ac.uk/pub/ensemblorganisms/a/b.bgz"

    def test_no_trailing_slash_on_base(self):
        url = _manifest_url("https://ftp.ebi.ac.uk/pub/ensemblorganisms", "a/b.bgz")
        assert url == "https://ftp.ebi.ac.uk/pub/ensemblorganisms/a/b.bgz"

    def test_leading_slash_on_path(self):
        url = _manifest_url("https://ftp.ebi.ac.uk/pub/ensemblorganisms", "/a/b.bgz")
        assert url == "https://ftp.ebi.ac.uk/pub/ensemblorganisms/a/b.bgz"

    def test_both_slashes_normalised(self):
        url = _manifest_url(
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/", "/a/b.bgz"
        )
        assert url == "https://ftp.ebi.ac.uk/pub/ensemblorganisms/a/b.bgz"

    def test_no_double_slash_in_middle(self):
        url = _manifest_url("https://base.com/pub/", "/relative/path")
        assert "//" not in url.replace("https://", "")


# ---------------------------------------------------------------------------
# TestLegacyVepManifestInit — valid fixture
# ---------------------------------------------------------------------------


class TestLegacyVepManifestInit:
    def _make(self) -> LegacyVepManifest:
        return LegacyVepManifest(_load_fixture())

    def test_missing_species_key_raises(self):
        with pytest.raises(LegacyVepManifestError, match="'species'"):
            LegacyVepManifest({"other": {}})

    def test_empty_dict_raises(self):
        with pytest.raises(LegacyVepManifestError, match="'species'"):
            LegacyVepManifest({})

    def test_non_dict_raises(self):
        with pytest.raises(LegacyVepManifestError):
            LegacyVepManifest(["not", "a", "dict"])  # type: ignore[arg-type]

    def test_len_counts_vep_records(self):
        """Fixture has VEP records for GCA_018852605.1, GCA_018472705.1,
        and two-provider GCA_888000001.1 (2 records) — total 4."""
        m = self._make()
        assert len(m) == 4

    def test_index_built_for_accession_with_vep(self):
        m = self._make()
        assert "GCA_018852605.1" in m._index

    def test_accession_without_vep_not_in_index(self):
        """GCA_999000001.1 has no vep section in the fixture."""
        m = self._make()
        assert "GCA_999000001.1" not in m._index


# ---------------------------------------------------------------------------
# TestLegacyVepManifestLookup
# ---------------------------------------------------------------------------


class TestLegacyVepManifestLookup:
    def _make(self) -> LegacyVepManifest:
        return LegacyVepManifest(_load_fixture())

    def test_lookup_accession_with_vep_returns_record(self):
        m = self._make()
        rec = m.lookup_vep("GCA_018852605.1")
        assert isinstance(rec, LegacyVepRecord)
        assert rec.accession == "GCA_018852605.1"
        assert rec.provider == "ensembl"
        assert rec.date_key == "2022_07"

    def test_lookup_returns_correct_primary_vep_file(self):
        m = self._make()
        rec = m.lookup_vep("GCA_018852605.1")
        assert rec is not None
        assert VEP_PRIMARY_FILE in rec.all_vep_files
        assert rec.vep_relative_path == rec.all_vep_files[VEP_PRIMARY_FILE]

    def test_vep_relative_path_exact_from_manifest(self):
        """The path must come verbatim from the manifest, not constructed."""
        m = self._make()
        rec = m.lookup_vep("GCA_018852605.1")
        assert rec is not None
        expected = (
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/genes.gff3.bgz"
        )
        assert rec.vep_relative_path == expected

    def test_lookup_accession_absent_returns_none(self):
        m = self._make()
        assert m.lookup_vep("GCA_000000000.1") is None

    def test_lookup_accession_no_vep_returns_none(self):
        """GCA_999000001.1 is in the fixture but has no vep section."""
        m = self._make()
        assert m.lookup_vep("GCA_999000001.1") is None

    def test_lookup_second_hprc_accession(self):
        """GCA_018472705.1 has ensembl/2022_08 VEP."""
        m = self._make()
        rec = m.lookup_vep("GCA_018472705.1")
        assert rec is not None
        assert rec.provider == "ensembl"
        assert rec.date_key == "2022_08"

    def test_provider_hint_resolves_ambiguous_record(self):
        """GCA_888000001.1 has ensembl + community; provider hint must pick correctly."""
        m = self._make()
        rec = m.lookup_vep("GCA_888000001.1", provider="ensembl")
        assert rec is not None
        assert rec.provider == "ensembl"

    def test_community_provider_hint_resolves_ambiguous_record(self):
        m = self._make()
        rec = m.lookup_vep("GCA_888000001.1", provider="community")
        assert rec is not None
        assert rec.provider == "community"

    def test_no_hint_on_ambiguous_returns_none(self):
        """Without a disambiguating hint, ambiguous records return None."""
        m = self._make()
        assert m.lookup_vep("GCA_888000001.1") is None

    def test_date_hint_resolves_single_provider(self):
        m = self._make()
        rec = m.lookup_vep("GCA_018852605.1", annotation_date="2022-07")
        assert rec is not None
        assert rec.date_key == "2022_07"

    def test_date_hint_dash_normalised(self):
        """2022-07 matches date_key 2022_07."""
        m = self._make()
        rec = m.lookup_vep("GCA_018472705.1", annotation_date="2022-08")
        assert rec is not None
        assert rec.date_key == "2022_08"

    def test_is_ambiguous_true_for_two_providers(self):
        m = self._make()
        assert m.is_ambiguous("GCA_888000001.1") is True

    def test_is_ambiguous_false_for_single_record(self):
        m = self._make()
        assert m.is_ambiguous("GCA_018852605.1") is False

    def test_is_ambiguous_false_for_absent(self):
        m = self._make()
        assert m.is_ambiguous("GCA_000000000.1") is False

    def test_is_ambiguous_resolved_by_provider(self):
        """After provider filtering, only one record remains → not ambiguous."""
        m = self._make()
        assert m.is_ambiguous("GCA_888000001.1", provider="ensembl") is False


# ---------------------------------------------------------------------------
# TestFromUrl — success and all failure modes
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# TestDirectoryPath — directory derivation from the primary VEP file path
# ---------------------------------------------------------------------------


class TestDirectoryPath:
    """Tests for LegacyVepRecord.directory_path."""

    def _record(
        self,
        vep_relative_path: str = (
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/genes.gff3.bgz"
        ),
    ) -> LegacyVepRecord:
        return LegacyVepRecord(
            accession="GCA_018852605.1",
            provider="ensembl",
            date_key="2022_07",
            vep_relative_path=vep_relative_path,
        )

    def test_directory_path_strips_filename(self):
        rec = self._record()
        assert rec.directory_path == (
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/"
        )

    def test_directory_path_ends_with_slash(self):
        assert self._record().directory_path.endswith("/")

    def test_directory_path_contains_no_filename(self):
        dp = self._record().directory_path
        assert "genes.gff3.bgz" not in dp

    def test_directory_path_preserves_date(self):
        dp = self._record().directory_path
        assert "2022_07" in dp

    def test_directory_path_preserves_accession(self):
        dp = self._record().directory_path
        assert "GCA_018852605.1" in dp

    def test_directory_path_preserves_provider(self):
        dp = self._record().directory_path
        assert "ensembl" in dp

    def test_directory_path_different_date(self):
        rec = self._record(
            "Homo_sapiens/GCA_018472705.1/vep/ensembl/geneset/2022_08/genes.gff3.bgz"
        )
        assert rec.directory_path == (
            "Homo_sapiens/GCA_018472705.1/vep/ensembl/geneset/2022_08/"
        )

    def test_directory_path_bare_filename_raises(self):
        """A path with no directory component must raise LegacyVepManifestError."""
        rec = self._record(vep_relative_path="genes.gff3.bgz")
        with pytest.raises(LegacyVepManifestError, match="Cannot derive VEP directory"):
            _ = rec.directory_path

    def test_full_url_from_directory_path_ends_with_slash(self):
        """Joining directory_path with EBI_FTP_BASE must produce a trailing-slash URL."""
        rec = self._record()
        full_url = _manifest_url(EBI_FTP_BASE, rec.directory_path)
        assert full_url == (
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/"
        )
        assert full_url.endswith("/")
        assert "genes.gff3.bgz" not in full_url

    def test_full_url_no_double_slash(self):
        rec = self._record()
        full_url = _manifest_url(EBI_FTP_BASE, rec.directory_path)
        assert "//" not in full_url.replace("https://", "")


# ---------------------------------------------------------------------------
# TestFromUrl — success and all failure modes
# ---------------------------------------------------------------------------


class TestFromUrl:
    def _make_response(self, data: dict, status_code: int = 200) -> MagicMock:
        mock = MagicMock()
        mock.status_code = status_code
        mock.json.return_value = data
        mock.raise_for_status = MagicMock()
        if status_code >= 400:
            mock.raise_for_status.side_effect = requests.exceptions.HTTPError(
                f"{status_code} Error"
            )
        return mock

    def test_success(self):
        with patch(
            "ensembl.genes.projects.legacy_vep_manifest.requests.get",
            return_value=self._make_response(_load_fixture()),
        ):
            m = LegacyVepManifest.from_url()
        assert isinstance(m, LegacyVepManifest)
        assert len(m) > 0

    def test_timeout_raises_legacy_vep_manifest_error(self):
        with patch(
            "ensembl.genes.projects.legacy_vep_manifest.requests.get",
            side_effect=requests.exceptions.Timeout("timed out"),
        ):
            with pytest.raises(LegacyVepManifestError, match="Timed out"):
                LegacyVepManifest.from_url()

    def test_http_error_raises_legacy_vep_manifest_error(self):
        with patch(
            "ensembl.genes.projects.legacy_vep_manifest.requests.get",
            return_value=self._make_response({}, status_code=503),
        ):
            with pytest.raises(LegacyVepManifestError, match="Could not download"):
                LegacyVepManifest.from_url()

    def test_malformed_json_raises_legacy_vep_manifest_error(self):
        mock = MagicMock()
        mock.raise_for_status = MagicMock()
        mock.json.side_effect = ValueError("No JSON object")
        with patch(
            "ensembl.genes.projects.legacy_vep_manifest.requests.get",
            return_value=mock,
        ):
            with pytest.raises(LegacyVepManifestError, match="not valid JSON"):
                LegacyVepManifest.from_url()

    def test_missing_species_key_raises_legacy_vep_manifest_error(self):
        with patch(
            "ensembl.genes.projects.legacy_vep_manifest.requests.get",
            return_value=self._make_response({"wrong_key": {}}),
        ):
            with pytest.raises(LegacyVepManifestError, match="'species'"):
                LegacyVepManifest.from_url()

    def test_exception_is_chained(self):
        with patch(
            "ensembl.genes.projects.legacy_vep_manifest.requests.get",
            side_effect=requests.exceptions.Timeout("timed out"),
        ):
            with pytest.raises(LegacyVepManifestError) as exc_info:
                LegacyVepManifest.from_url()
        assert exc_info.value.__cause__ is not None

    def test_relative_path_vs_full_url_handling(self):
        """_manifest_url correctly joins base with relative path from manifest."""
        relative = (
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/genes.gff3.bgz"
        )
        full = _manifest_url(EBI_FTP_BASE, relative)
        assert full.startswith("https://")
        assert "ensemblorganisms" in full
        assert "//" not in full.replace("https://", "")

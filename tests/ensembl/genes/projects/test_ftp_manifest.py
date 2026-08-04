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
"""Tests for ftp_manifest — all network access is mocked."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from ensembl.genes.projects.ftp_manifest import (
    EBI_FTP_BASE,
    EnsemblFtpManifest,
    ManifestError,
    _parse_manifest_date,
)

# Load the shared fixture once at module level.
_FIXTURE_PATH = Path(__file__).parent / "fixtures" / "manifest_sample.json"
with _FIXTURE_PATH.open() as _f:
    _FIXTURE_DATA = json.load(_f)


class TestParseManifestedDate:
    """Unit tests for _parse_manifest_date."""

    def test_yyyy_mm(self):
        assert _parse_manifest_date("2022_11") == (2022, 11)

    def test_yyyy_mm_dd(self):
        assert _parse_manifest_date("2023_10_18") == (2023, 10, 18)

    def test_yyyy_dash_mm(self):
        assert _parse_manifest_date("2022-11") == (2022, 11)

    def test_yyyy_dash_mm_dash_dd(self):
        assert _parse_manifest_date("2022-11-01") == (2022, 11, 1)

    def test_invalid_month_returns_none(self):
        assert _parse_manifest_date("2024_99") is None

    def test_invalid_day_returns_none(self):
        assert _parse_manifest_date("2024_02_45") is None

    def test_garbage_returns_none(self):
        assert _parse_manifest_date("bad") is None

    def test_empty_returns_none(self):
        assert _parse_manifest_date("") is None


class TestEnsemblFtpManifestInit:
    """Tests for EnsemblFtpManifest construction from dict."""

    def test_valid_fixture_loads(self):
        manifest = EnsemblFtpManifest(_FIXTURE_DATA)
        assert len(manifest) > 0

    def test_empty_dict_raises(self):
        with pytest.raises(ManifestError, match="'species'"):
            EnsemblFtpManifest({})

    def test_missing_species_key_raises(self):
        with pytest.raises(ManifestError, match="'species'"):
            EnsemblFtpManifest({"other_key": {}})

    def test_duplicate_accession_raises(self):
        """Duplicate accessions across species entries must raise ManifestError."""
        dup_data = {
            "species": {
                "Species_a": {
                    "assemblies": {
                        "GCA_000000001.1": {
                            "genebuild_providers": {},
                            "assembly": {"files": {"genome_sequences": {}}},
                        }
                    }
                },
                "Species_b": {
                    "assemblies": {
                        "GCA_000000001.1": {
                            "genebuild_providers": {},
                            "assembly": {"files": {"genome_sequences": {}}},
                        }
                    }
                },
            }
        }
        with pytest.raises(ManifestError, match="Duplicate accession"):
            EnsemblFtpManifest(dup_data)


class TestEnsemblFtpManifestLookup:
    """Tests for EnsemblFtpManifest.lookup()."""

    def setup_method(self):
        self.manifest = EnsemblFtpManifest(_FIXTURE_DATA)

    def test_lookup_gca_ensembl(self):
        """GCA ensembl-provider assembly is found and has expected files."""
        record = self.manifest.lookup("GCA_922984935.2")
        assert record is not None
        assert "ensembl" in record.providers
        pdr = record.providers["ensembl"]["2022_11"]
        assert "genes.gtf.gz" in pdr.annotation_files
        assert pdr.annotation_files["genes.gtf.gz"].startswith("GCA/922/984/935/2/")

    def test_lookup_gcf_refseq(self):
        """GCF refseq-provider assembly is found and has expected files."""
        record = self.manifest.lookup("GCF_021464435.1")
        assert record is not None
        assert "refseq" in record.providers

    def test_lookup_community_via_ambiguous_species(self):
        """Assembly with two non-ensembl providers (community, braker) is found."""
        record = self.manifest.lookup("GCA_999999999.1")
        assert record is not None
        assert "community" in record.providers
        assert "braker" in record.providers

    def test_lookup_absent_accession_returns_none(self):
        """Accession not in manifest returns None, not an exception."""
        assert self.manifest.lookup("GCA_000000000.0") is None

    def test_assembly_genome_files_on_record(self):
        """Genome sequence files are stored on AssemblyRecord, not ProviderDateRecord."""
        record = self.manifest.lookup("GCA_922984935.2")
        assert "softmasked.fa.bgz" in record.assembly_genome_files
        # Genome files must NOT appear inside ProviderDateRecord
        pdr = record.providers["ensembl"]["2022_11"]
        assert not hasattr(pdr, "genome_files")

    def test_homology_path_verbatim(self):
        """Homology relative path (including YYYY_MM_DD subdir) is preserved verbatim."""
        record = self.manifest.lookup("GCA_922984935.2")
        pdr = record.providers["ensembl"]["2022_11"]
        homology_path = pdr.homology_files.get("homology.tsv.gz", "")
        assert "2024_09_18" in homology_path

    def test_multiple_dates_for_one_provider(self):
        """Assembly with two annotation dates exposes both."""
        record = self.manifest.lookup("GCA_012295145.1")
        assert "2023_10" in record.providers["ensembl"]
        assert "2024_02" in record.providers["ensembl"]

    def test_no_gtf_assembly_is_indexed(self):
        """Assembly with no GTF/GFF3 is still indexed; absence is determined at render time."""
        record = self.manifest.lookup("GCA_888888888.1")
        assert record is not None


class TestEnsemblFtpManifestFromUrl:
    """Tests for EnsemblFtpManifest.from_url() — all network calls mocked."""

    def test_success(self):
        """A valid HTTP response is parsed and an instance is returned."""
        mock_response = MagicMock()
        mock_response.json.return_value = _FIXTURE_DATA
        with patch("requests.get", return_value=mock_response):
            manifest = EnsemblFtpManifest.from_url()
        assert len(manifest) > 0

    def test_timeout_raises_manifest_error(self):
        """A request timeout raises ManifestError (not a raw requests exception)."""
        import requests as req_mod

        with patch("requests.get", side_effect=req_mod.exceptions.Timeout()):
            with pytest.raises(ManifestError, match="Timed out"):
                EnsemblFtpManifest.from_url()

    def test_http_error_raises_manifest_error(self):
        """A non-200 HTTP response raises ManifestError."""
        import requests as req_mod

        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = req_mod.exceptions.HTTPError("503")
        with patch("requests.get", return_value=mock_response):
            with pytest.raises(ManifestError, match="Could not download"):
                EnsemblFtpManifest.from_url()

    def test_malformed_json_raises_manifest_error(self):
        """A response that is not valid JSON raises ManifestError."""
        mock_response = MagicMock()
        mock_response.json.side_effect = ValueError("No JSON object could be decoded")
        with patch("requests.get", return_value=mock_response):
            with pytest.raises(ManifestError, match="not valid JSON"):
                EnsemblFtpManifest.from_url()

    def test_invalid_structure_raises_manifest_error(self):
        """A valid JSON response with wrong schema raises ManifestError."""
        mock_response = MagicMock()
        mock_response.json.return_value = {"not_species": {}}
        with patch("requests.get", return_value=mock_response):
            with pytest.raises(ManifestError, match="'species'"):
                EnsemblFtpManifest.from_url()

    def test_exception_is_chained(self):
        """ManifestError wraps the original exception via exception chaining."""
        import requests as req_mod

        with patch("requests.get", side_effect=req_mod.exceptions.Timeout()):
            with pytest.raises(ManifestError) as exc_info:
                EnsemblFtpManifest.from_url()
        # The original cause must be preserved
        assert exc_info.value.__cause__ is not None

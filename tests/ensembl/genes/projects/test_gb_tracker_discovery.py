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
"""Tests for GbTrackerClient query logic — all database calls mocked.

Verifies that the SQL queries produced by ``fetch_by_identifier`` and
``fetch_project_pre_releases`` contain the correct status, last_attempt,
and release_type filters for discovering active unreleased genebuilds
including ``handed_over`` records.
"""

from unittest.mock import MagicMock, patch

import pytest

from ensembl.genes.projects.config import ProjectConfig
from ensembl.genes.projects.registry.gb_tracker import GbTrackerClient


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_client() -> GbTrackerClient:
    return GbTrackerClient(host="test-host", port=3306, user="test-user")


def _make_config(
    bioproject_scoping: list[str] | None = None,
    custom_group_scoping: list[str] | None = None,
) -> ProjectConfig:
    return ProjectConfig(
        project_name="test",
        bioproject_scoping=bioproject_scoping or ["CBP"],
        custom_group_scoping=custom_group_scoping,
    )


def _make_row(
    accession: str = "GCA_051175935.2",
    species_name: str = "Anthias nicholsi",
    asm_name: str = "fAntiNich1.1",
    annotation_method: str = "full_genebuild",
    busco_score: str | None = None,
    busco_lineage: str | None = None,
) -> dict:
    """Build a mock database row matching the SELECT column aliases."""
    return {
        "genome_uuid": None,
        "core_dbname": None,
        "assembly_accession": accession,
        "species_name": species_name,
        "taxon_id": 481325,
        "assembly_name": asm_name,
        "annotation_method": annotation_method,
        "busco_score": busco_score,
        "busco_lineage": busco_lineage,
    }


# ---------------------------------------------------------------------------
# Automatic project discovery: fetch_project_pre_releases
# ---------------------------------------------------------------------------


class TestProjectDiscovery:
    """Tests for ``fetch_project_pre_releases`` SQL query logic."""

    def test_discovery_pre_released_still_discoverable(self):
        """Active pre_released record is returned by fetch_project_pre_releases."""
        client = _make_client()
        config = _make_config()
        row = _make_row()

        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchall.return_value = [row]
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cursor)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)
        mock_conn.open = True

        with patch("pymysql.connect", return_value=mock_conn):
            results = client.fetch_project_pre_releases(config)

        assert len(results) == 1
        assert results[0].accession == "GCA_051175935.2"
        assert results[0].species_name == "Anthias nicholsi"
        assert results[0].is_released is False

    def test_discovery_handed_over_unreleased_discoverable(self):
        """Active handed_over + release_type='not_available' + last_attempt=1
        is returned by fetch_project_pre_releases (the query includes it)."""
        client = _make_client()
        config = _make_config()
        row = _make_row()

        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchall.return_value = [row]
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cursor)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)
        mock_conn.open = True

        with patch("pymysql.connect", return_value=mock_conn):
            results = client.fetch_project_pre_releases(config)

        assert len(results) == 1
        meta = results[0]
        assert meta.accession == "GCA_051175935.2"
        assert meta.annotation_method == "Ensembl Genebuild"
        assert meta.is_released is False

    def test_discovery_query_contains_status_filters(self):
        """The executed SQL contains the correct status, last_attempt, and
        release_type filters."""
        client = _make_client()
        config = _make_config()

        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchall.return_value = []
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cursor)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)
        mock_conn.open = True

        with patch("pymysql.connect", return_value=mock_conn):
            client.fetch_project_pre_releases(config)

        executed_sql = mock_cursor.execute.call_args[0][0]
        assert "'pre_released'" in executed_sql
        assert "'handed_over'" in executed_sql
        assert "last_attempt = 1" in executed_sql
        assert "release_type = 'not_available'" in executed_sql

    def test_discovery_no_scoping_returns_empty(self):
        """When no bioproject or custom_group scoping is configured, return []
        without executing a query."""
        client = _make_client()
        config = ProjectConfig(
            project_name="test",
            bioproject_scoping=None,
            custom_group_scoping=None,
        )

        with patch("pymysql.connect") as mock_connect:
            results = client.fetch_project_pre_releases(config)

        assert results == []
        mock_connect.assert_not_called()


# ---------------------------------------------------------------------------
# Explicit identifier lookup: fetch_by_identifier
# ---------------------------------------------------------------------------


class TestExplicitLookup:
    """Tests for ``fetch_by_identifier`` SQL query logic."""

    def test_explicit_lookup_handed_over_unreleased_succeeds(self):
        """Unreleased handed_over record resolves via fetch_by_identifier."""
        client = _make_client()
        row = _make_row()

        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchone.return_value = row
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cursor)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)
        mock_conn.open = True

        with patch("pymysql.connect", return_value=mock_conn):
            meta = client.fetch_by_identifier("GCA_051175935.2")

        assert meta is not None
        assert meta.accession == "GCA_051175935.2"
        assert meta.species_name == "Anthias nicholsi"
        assert meta.annotation_method == "Ensembl Genebuild"
        assert meta.is_released is False

    def test_explicit_lookup_in_progress_still_works(self):
        """in_progress record resolves via fetch_by_identifier (regression guard)."""
        client = _make_client()
        row = _make_row(
            accession="GCA_000000001.1",
            species_name="Test species",
            annotation_method="braker",
        )

        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchone.return_value = row
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cursor)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)
        mock_conn.open = True

        with patch("pymysql.connect", return_value=mock_conn):
            meta = client.fetch_by_identifier("GCA_000000001.1")

        assert meta is not None
        assert meta.annotation_method == "BRAKER2"

    def test_explicit_query_contains_status_filters(self):
        """The executed SQL contains the correct status, last_attempt, and
        release_type filters."""
        client = _make_client()

        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchone.return_value = None
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cursor)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)
        mock_conn.open = True

        with patch("pymysql.connect", return_value=mock_conn):
            client.fetch_by_identifier("GCA_051175935.2")

        executed_sql = mock_cursor.execute.call_args[0][0]
        assert "'in_progress'" in executed_sql
        assert "'pre_released'" in executed_sql
        assert "'handed_over'" in executed_sql
        assert "last_attempt = 1" in executed_sql
        assert "release_type = 'not_available'" in executed_sql

    def test_explicit_lookup_not_found_returns_none(self):
        """When no matching row is found, fetch_by_identifier returns None."""
        client = _make_client()

        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchone.return_value = None
        mock_conn.cursor.return_value.__enter__ = MagicMock(return_value=mock_cursor)
        mock_conn.cursor.return_value.__exit__ = MagicMock(return_value=False)
        mock_conn.open = True

        with patch("pymysql.connect", return_value=mock_conn):
            meta = client.fetch_by_identifier("GCA_999999999.1")

        assert meta is None

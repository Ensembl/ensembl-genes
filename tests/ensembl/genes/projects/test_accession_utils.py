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
"""Tests for accession_utils.accession_to_ftp_path."""

import pytest

from ensembl.genes.projects.accession_utils import accession_to_ftp_path


class TestAccessionToFtpPath:
    """Unit tests for accession_to_ftp_path."""

    # ------------------------------------------------------------------
    # Happy-path conversions
    # ------------------------------------------------------------------

    def test_gca_basic(self):
        """GCA accession produces the correct triplet path."""
        assert accession_to_ftp_path("GCA_922984935.2") == "GCA/922/984/935/2"

    def test_gcf_basic(self):
        """GCF accession produces the correct triplet path."""
        assert accession_to_ftp_path("GCF_000001405.40") == "GCF/000/001/405/40"

    def test_leading_zeroes_preserved(self):
        """Leading zeros in the nine-digit body are preserved as directory names."""
        assert accession_to_ftp_path("GCA_000001405.29") == "GCA/000/001/405/29"

    def test_multi_digit_version(self):
        """Multi-digit version numbers are placed as a single path component."""
        assert accession_to_ftp_path("GCF_021464435.1") == "GCF/021/464/435/1"

    def test_version_1(self):
        """Version 1 (single digit) is handled correctly."""
        assert accession_to_ftp_path("GCA_012295145.1") == "GCA/012/295/145/1"

    # ------------------------------------------------------------------
    # Error cases
    # ------------------------------------------------------------------

    def test_wrong_prefix_raises(self):
        """Unrecognised prefix raises ValueError."""
        with pytest.raises(ValueError, match="Malformed"):
            accession_to_ftp_path("GCX_922984935.2")

    def test_missing_prefix_raises(self):
        """Accession without GCA/GCF prefix raises ValueError."""
        with pytest.raises(ValueError):
            accession_to_ftp_path("922984935.2")

    def test_wrong_digit_count_raises(self):
        """Body with fewer than nine digits raises ValueError."""
        with pytest.raises(ValueError):
            accession_to_ftp_path("GCA_92298493.2")

    def test_too_many_digits_raises(self):
        """Body with more than nine digits raises ValueError."""
        with pytest.raises(ValueError):
            accession_to_ftp_path("GCA_9229849350.2")

    def test_missing_version_raises(self):
        """Accession with no version (no dot) raises ValueError."""
        with pytest.raises(ValueError):
            accession_to_ftp_path("GCA_922984935")

    def test_non_integer_version_raises(self):
        """Accession with a non-integer version raises ValueError."""
        with pytest.raises(ValueError):
            accession_to_ftp_path("GCA_922984935.abc")

    def test_empty_string_raises(self):
        """Empty string raises ValueError."""
        with pytest.raises(ValueError):
            accession_to_ftp_path("")

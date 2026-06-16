"""Tests for the haplotype resolver.

All NCBI API calls are mocked to avoid network dependencies.
Run with:
    pytest src/python/ensembl/genes/projects/test_haplotype_resolver.py -v
"""

from dataclasses import dataclass
from typing import List, Optional
from unittest.mock import MagicMock, patch

from ensembl.genes.projects.haplotype_resolver import (
    HaplotypeResolver,
    _extract_hap_group_key,
)


# ---------------------------------------------------------------------------
# Minimal GenomeMetadata stand-in
# ---------------------------------------------------------------------------


@dataclass
class _FakeMeta:
    accession: str
    species_name: str = "Test species"
    assembly_name: str = ""
    taxon_id: Optional[int] = 9606


# ---------------------------------------------------------------------------
# Mock NCBI Datasets API helper
# ---------------------------------------------------------------------------


def _mock_ncbi_response(accession_to_info):
    """Build a mock requests.get that returns NCBI Datasets-style JSON.

    accession_to_info: dict mapping accession → {biosample, sample_name, assembly_name, taxon_id}
    """
    reports = []
    for acc, info in accession_to_info.items():
        biosample_attrs = []
        if info.get("sample_name"):
            biosample_attrs.append(
                {"name": "sample_name", "value": info["sample_name"]}
            )

        report = {
            "accession": acc,
            "organism": {"tax_id": info.get("taxon_id", 9606)},
            "assembly_info": {
                "assembly_name": info.get("assembly_name", ""),
                "biosample": {
                    "accession": info.get("biosample", ""),
                    "attributes": biosample_attrs,
                },
            },
        }
        reports.append(report)

    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {"reports": reports}
    return mock_resp


def _no_ncbi():
    """Patch NCBI to return empty results."""
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {"reports": []}
    return patch(
        "ensembl.genes.projects.haplotype_resolver.requests.get",
        return_value=mock_resp,
    )


# ---------------------------------------------------------------------------
# Source 1: BioSample pairing
# ---------------------------------------------------------------------------


class TestBioSamplePairing:

    def test_two_assemblies_same_biosample(self):
        """Two assemblies sharing a BioSample should be paired."""
        genomes = [
            _FakeMeta(accession="GCA_001", assembly_name="asm_hap1"),
            _FakeMeta(accession="GCA_002", assembly_name="asm_hap2"),
            _FakeMeta(accession="GCA_003", assembly_name="other_species"),
        ]
        ncbi_info = {
            "GCA_001": {"biosample": "SAMN111", "assembly_name": "asm_hap1"},
            "GCA_002": {"biosample": "SAMN111", "assembly_name": "asm_hap2"},
            "GCA_003": {"biosample": "SAMN222", "assembly_name": "other"},
        }
        with patch(
            "ensembl.genes.projects.haplotype_resolver.requests.get",
            return_value=_mock_ncbi_response(ncbi_info),
        ):
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert pairs.get("GCA_001") == "GCA_002"
        assert pairs.get("GCA_002") == "GCA_001"
        assert "GCA_003" not in pairs

    def test_three_assemblies_same_biosample_not_paired(self):
        """Three assemblies with same BioSample — ambiguous, skip."""
        genomes = [
            _FakeMeta(accession="GCA_001"),
            _FakeMeta(accession="GCA_002"),
            _FakeMeta(accession="GCA_003"),
        ]
        ncbi_info = {
            "GCA_001": {"biosample": "SAMN111"},
            "GCA_002": {"biosample": "SAMN111"},
            "GCA_003": {"biosample": "SAMN111"},
        }
        with patch(
            "ensembl.genes.projects.haplotype_resolver.requests.get",
            return_value=_mock_ncbi_response(ncbi_info),
        ):
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        # Ambiguous — should not auto-pair
        assert len(pairs) == 0

    def test_no_biosample_no_crash(self):
        """Assemblies without BioSample data don't cause errors."""
        genomes = [
            _FakeMeta(accession="GCA_001"),
            _FakeMeta(accession="GCA_002"),
        ]
        ncbi_info = {
            "GCA_001": {"biosample": ""},
            "GCA_002": {"biosample": ""},
        }
        with patch(
            "ensembl.genes.projects.haplotype_resolver.requests.get",
            return_value=_mock_ncbi_response(ncbi_info),
        ):
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert len(pairs) == 0


# ---------------------------------------------------------------------------
# Source 2: Sample name pairing
# ---------------------------------------------------------------------------


class TestSampleNamePairing:

    def test_same_sample_name_pairs(self):
        """Assemblies with the same sample_name should be paired."""
        genomes = [
            _FakeMeta(accession="GCA_001", taxon_id=9606),
            _FakeMeta(accession="GCA_002", taxon_id=9606),
        ]
        ncbi_info = {
            "GCA_001": {
                "biosample": "",  # no biosample
                "sample_name": "individual_42",
                "taxon_id": 9606,
            },
            "GCA_002": {
                "biosample": "",
                "sample_name": "individual_42",
                "taxon_id": 9606,
            },
        }
        with patch(
            "ensembl.genes.projects.haplotype_resolver.requests.get",
            return_value=_mock_ncbi_response(ncbi_info),
        ):
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert pairs.get("GCA_001") == "GCA_002"
        assert pairs.get("GCA_002") == "GCA_001"


# ---------------------------------------------------------------------------
# Source 3: Assembly name heuristics
# ---------------------------------------------------------------------------


class TestAssemblyNameHeuristics:

    def test_hap1_hap2_naming(self):
        """Assembly names with hap1/hap2 should be paired."""
        genomes = [
            _FakeMeta(
                accession="GCA_001",
                species_name="Homo sapiens",
                assembly_name="HG002.hap1",
                taxon_id=9606,
            ),
            _FakeMeta(
                accession="GCA_002",
                species_name="Homo sapiens",
                assembly_name="HG002.hap2",
                taxon_id=9606,
            ),
        ]
        with _no_ncbi():
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert pairs.get("GCA_001") == "GCA_002"
        assert pairs.get("GCA_002") == "GCA_001"

    def test_maternal_paternal_naming(self):
        """Assembly names with mat/pat should be paired."""
        genomes = [
            _FakeMeta(
                accession="GCA_A",
                species_name="Species X",
                assembly_name="sX_maternal",
                taxon_id=1234,
            ),
            _FakeMeta(
                accession="GCA_B",
                species_name="Species X",
                assembly_name="sX_paternal",
                taxon_id=1234,
            ),
        ]
        with _no_ncbi():
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert pairs.get("GCA_A") == "GCA_B"

    def test_primary_alternate_naming(self):
        genomes = [
            _FakeMeta(
                accession="GCA_P",
                species_name="Species Y",
                assembly_name="spY_primary",
                taxon_id=5678,
            ),
            _FakeMeta(
                accession="GCA_Q",
                species_name="Species Y",
                assembly_name="spY_alternate",
                taxon_id=5678,
            ),
        ]
        with _no_ncbi():
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert pairs.get("GCA_P") == "GCA_Q"

    def test_different_species_not_paired(self):
        """Even with matching naming, different species must not pair."""
        genomes = [
            _FakeMeta(
                accession="GCA_X",
                species_name="Species A",
                assembly_name="common_hap1",
                taxon_id=111,
            ),
            _FakeMeta(
                accession="GCA_Y",
                species_name="Species B",
                assembly_name="common_hap2",
                taxon_id=222,
            ),
        ]
        with _no_ncbi():
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert len(pairs) == 0

    def test_no_haplotype_pattern_no_pair(self):
        """Assemblies without haplotype naming patterns are not paired."""
        genomes = [
            _FakeMeta(
                accession="GCA_M",
                species_name="Normal species",
                assembly_name="normal_assembly_v1",
                taxon_id=999,
            ),
            _FakeMeta(
                accession="GCA_N",
                species_name="Normal species",
                assembly_name="another_assembly_v2",
                taxon_id=999,
            ),
        ]
        with _no_ncbi():
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        assert len(pairs) == 0


# ---------------------------------------------------------------------------
# Helper function tests
# ---------------------------------------------------------------------------


class TestExtractHapGroupKey:

    def test_hap1(self):
        result = _extract_hap_group_key("HG002.hap1")
        assert result is not None
        key, label = result
        assert label == "hap1"

    def test_hap2(self):
        result = _extract_hap_group_key("HG002.hap2")
        assert result is not None
        key, label = result
        assert label == "hap2"

    def test_maternal(self):
        result = _extract_hap_group_key("sample_maternal")
        assert result is not None
        _, label = result
        assert label == "maternal"

    def test_paternal(self):
        result = _extract_hap_group_key("sample_paternal")
        assert result is not None
        _, label = result
        assert label == "paternal"

    def test_primary(self):
        result = _extract_hap_group_key("asm_primary")
        assert result is not None
        _, label = result
        assert label == "primary"

    def test_alternate(self):
        result = _extract_hap_group_key("asm_alternate")
        assert result is not None
        _, label = result
        assert label == "alternate"

    def test_same_group_key_for_pair(self):
        r1 = _extract_hap_group_key("HG002.hap1")
        r2 = _extract_hap_group_key("HG002.hap2")
        assert r1 is not None and r2 is not None
        assert r1[0] == r2[0]  # same group key
        assert r1[1] != r2[1]  # different labels

    def test_no_match(self):
        assert _extract_hap_group_key("normal_assembly_v3") is None
        assert _extract_hap_group_key("") is None
        assert _extract_hap_group_key(None) is None


# ---------------------------------------------------------------------------
# Priority: BioSample wins over naming
# ---------------------------------------------------------------------------


class TestPriority:

    def test_biosample_takes_priority(self):
        """BioSample pairing should take priority over naming patterns."""
        genomes = [
            _FakeMeta(
                accession="GCA_001",
                assembly_name="sample_hap1",
                taxon_id=9606,
            ),
            _FakeMeta(
                accession="GCA_002",
                assembly_name="sample_hap2",
                taxon_id=9606,
            ),
            _FakeMeta(
                accession="GCA_003",
                assembly_name="other_hap1",
                taxon_id=9606,
            ),
        ]
        # BioSample says 001 and 003 are from the same individual,
        # even though naming suggests 001+002
        ncbi_info = {
            "GCA_001": {"biosample": "SAMN_A", "assembly_name": "sample_hap1"},
            "GCA_002": {"biosample": "SAMN_B", "assembly_name": "sample_hap2"},
            "GCA_003": {"biosample": "SAMN_A", "assembly_name": "other_hap1"},
        }
        with patch(
            "ensembl.genes.projects.haplotype_resolver.requests.get",
            return_value=_mock_ncbi_response(ncbi_info),
        ):
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        # BioSample pairing wins: 001 <-> 003
        assert pairs.get("GCA_001") == "GCA_003"
        assert pairs.get("GCA_003") == "GCA_001"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:

    def test_empty_genome_list(self):
        with _no_ncbi():
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes([])
        assert pairs == {}

    def test_single_genome(self):
        with _no_ncbi():
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(
                [_FakeMeta(accession="GCA_ONLY")]
            )
        assert pairs == {}

    def test_ncbi_api_failure(self):
        """If NCBI API fails, naming fallback should still work."""
        genomes = [
            _FakeMeta(
                accession="GCA_X",
                species_name="Sp",
                assembly_name="x_hap1",
                taxon_id=1,
            ),
            _FakeMeta(
                accession="GCA_Y",
                species_name="Sp",
                assembly_name="x_hap2",
                taxon_id=1,
            ),
        ]
        mock_resp = MagicMock()
        mock_resp.status_code = 500
        with patch(
            "ensembl.genes.projects.haplotype_resolver.requests.get",
            return_value=mock_resp,
        ):
            resolver = HaplotypeResolver()
            pairs = resolver.find_alternate_haplotypes(genomes)

        # Should fall back to naming pattern
        assert pairs.get("GCA_X") == "GCA_Y"

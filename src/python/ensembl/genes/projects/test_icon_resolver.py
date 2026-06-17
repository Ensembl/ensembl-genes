"""Tests for the taxonomy-lineage-based icon resolver.

These tests mock NCBI responses to avoid network calls during CI.
Run with:
    pytest src/python/ensembl/genes/projects/test_icon_resolver.py -v
"""

import os
import tempfile
from dataclasses import dataclass, field
from typing import List, Optional
from unittest.mock import MagicMock, patch

import pytest

from ensembl.genes.projects.icon_resolver import IconResolver, _FALLBACK_ICON


# ---------------------------------------------------------------------------
# Minimal stand-in for GenomeMetadata (avoids import cycles in tests)
# ---------------------------------------------------------------------------


@dataclass
class _FakeMeta:
    accession: str = "GCA_TEST"
    taxon_id: Optional[int] = 1
    taxonomy_lineage: Optional[List[str]] = None
    busco_lineage: Optional[str] = None


# ---------------------------------------------------------------------------
# Helpers: synthetic NCBI Entrez XML responses
# ---------------------------------------------------------------------------


def _make_entrez_xml(lineage_names):
    """Build a minimal NCBI efetch XML with the given lineage (root → leaf)."""
    taxon_items = "\n".join(
        f"<Taxon><ScientificName>{name}</ScientificName></Taxon>"
        for name in lineage_names
    )
    return f"""<?xml version="1.0" ?>
<TaxaSet>
  <Taxon>
    <TaxId>0</TaxId>
    <ScientificName>Test species</ScientificName>
    <LineageEx>{taxon_items}</LineageEx>
  </Taxon>
</TaxaSet>"""


def _mock_response(xml_text, status_code=200):
    resp = MagicMock()
    resp.status_code = status_code
    resp.text = xml_text
    return resp


def _no_ncbi():
    """Patch NCBI to always return empty (simulating network failure)."""
    return patch(
        "ensembl.genes.projects.icon_resolver.requests.get",
        return_value=_mock_response("", status_code=500),
    )


def _with_ncbi(lineage_root_to_leaf):
    """Patch NCBI to return a specific lineage."""
    xml = _make_entrez_xml(lineage_root_to_leaf)
    return patch(
        "ensembl.genes.projects.icon_resolver.requests.get",
        return_value=_mock_response(xml),
    )


# ---------------------------------------------------------------------------
# Source 1: metadata DB lineage (highest priority)
# ---------------------------------------------------------------------------


class TestMetadataLineage:
    """When taxonomy_lineage is pre-populated, it should be used first."""

    def test_metadata_lineage_used_first(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=[
                    "Lepidoptera",
                    "Insecta",
                    "Arthropoda",
                    "Metazoa",
                ],
            )
            icon, matched, source = resolver.resolve_icon(meta)
            assert icon == "Lepidoptera.png"
            assert matched == "Lepidoptera"
            assert source == "metadata_lineage"

    def test_metadata_lineage_arthropod(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=["Arachnida", "Arthropoda", "Metazoa"],
            )
            icon, matched, source = resolver.resolve_icon(meta)
            assert icon == "Arthropods.png"
            assert source == "metadata_lineage"

    def test_metadata_lineage_bird(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=["Aves", "Vertebrata", "Chordata", "Metazoa"],
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Birds.png"
            assert source == "metadata_lineage"

    def test_metadata_lineage_plant(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=[
                    "Magnoliopsida",
                    "Streptophyta",
                    "Viridiplantae",
                ],
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Plants.png"
            assert source == "metadata_lineage"


# ---------------------------------------------------------------------------
# Source 2: NCBI Entrez lineage
# ---------------------------------------------------------------------------


class TestNCBILineage:
    """When metadata lineage is absent, NCBI Entrez should be used."""

    def test_ncbi_lepidoptera(self):
        lineage = ["Eukaryota", "Metazoa", "Arthropoda", "Insecta", "Lepidoptera"]
        with _with_ncbi(lineage):
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=42)
            icon, matched, source = resolver.resolve_icon(meta)
            assert icon == "Lepidoptera.png"
            assert source == "ncbi_lineage"

    def test_ncbi_mammal(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Mammalia"]
        with _with_ncbi(lineage):
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=42)
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"
            assert source == "ncbi_lineage"

    def test_ncbi_fish(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Chondrichthyes",
        ]
        with _with_ncbi(lineage):
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=42)
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Fish.png"
            assert source == "ncbi_lineage"

    def test_ncbi_caching(self):
        """NCBI should only be called once per taxon_id."""
        lineage = ["Eukaryota", "Metazoa", "Arthropoda"]
        with _with_ncbi(lineage) as mock_get:
            resolver = IconResolver(icons_file="/nonexistent")
            resolver.resolve_icon(_FakeMeta(taxon_id=42))
            resolver.resolve_icon(_FakeMeta(taxon_id=42))
            assert mock_get.call_count == 1


# ---------------------------------------------------------------------------
# Source 3: BUSCO lineage fallback
# ---------------------------------------------------------------------------


class TestBUSCOFallback:
    """When both metadata and NCBI fail, BUSCO lineage should be used."""

    def test_busco_insecta(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="insecta_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Arthropods.png"
            assert source == "busco_lineage"

    def test_busco_lepidoptera(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="lepidoptera_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Lepidoptera.png"
            assert source == "busco_lineage"

    def test_busco_mammalia(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="mammalia_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"
            assert source == "busco_lineage"

    def test_busco_aves(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="aves_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Birds.png"
            assert source == "busco_lineage"

    def test_busco_actinopterygii(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="actinopterygii_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Fish.png"
            assert source == "busco_lineage"

    def test_busco_fungi(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="ascomycota_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Fungi.png"
            assert source == "busco_lineage"

    def test_busco_plant(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="eudicots_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Plants.png"
            assert source == "busco_lineage"

    def test_busco_vertebrata(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="vertebrata_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Chordates.png"
            assert source == "busco_lineage"

    def test_busco_laurasiatheria_maps_to_mammals(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="laurasiatheria_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"
            assert source == "busco_lineage"

    def test_busco_rodentia_maps_to_mammals(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="rodentia_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"
            assert source == "busco_lineage"

    def test_busco_galloanserae_maps_to_birds(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="galloanserae_odb10")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Birds.png"
            assert source == "busco_lineage"


# ---------------------------------------------------------------------------
# Priority ordering between sources
# ---------------------------------------------------------------------------


class TestSourcePriority:
    """Metadata lineage should win over NCBI, which should win over BUSCO."""

    def test_metadata_lineage_overrides_ncbi(self):
        """Even if NCBI would give a different answer, metadata wins."""
        ncbi_lineage = ["Eukaryota", "Fungi", "Ascomycota"]
        with _with_ncbi(ncbi_lineage):
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxon_id=42,
                taxonomy_lineage=["Mammalia", "Vertebrata", "Chordata"],
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"
            assert source == "metadata_lineage"

    def test_ncbi_overrides_busco(self):
        """NCBI should be preferred over BUSCO when both are available."""
        lineage = ["Eukaryota", "Metazoa", "Arthropoda", "Insecta", "Lepidoptera"]
        with _with_ncbi(lineage):
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxon_id=42,
                busco_lineage="mammalia_odb10",  # contradicts NCBI
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Lepidoptera.png"
            assert source == "ncbi_lineage"

    def test_busco_used_when_ncbi_fails(self):
        """When NCBI fails, BUSCO should be the fallback."""
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxon_id=42,
                busco_lineage="insecta_odb10",
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Arthropods.png"
            assert source == "busco_lineage"


# ---------------------------------------------------------------------------
# Specificity: most-specific match wins
# ---------------------------------------------------------------------------


class TestSpecificity:

    def test_lepidoptera_beats_arthropoda(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=[
                    "Lepidoptera",
                    "Insecta",
                    "Arthropoda",
                    "Metazoa",
                ],
            )
            icon, matched, _ = resolver.resolve_icon(meta)
            assert icon == "Lepidoptera.png"
            assert matched == "Lepidoptera"

    def test_mammalia_beats_chordata(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=["Mammalia", "Vertebrata", "Chordata"],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"

    def test_aves_beats_chordata(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=["Aves", "Vertebrata", "Chordata"],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Birds.png"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:

    def test_no_taxon_id_no_busco(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage=None)
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == _FALLBACK_ICON
            assert source == "missing_taxon_id"

    def test_all_sources_empty(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=999)
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == _FALLBACK_ICON
            assert source == "fallback"


# ---------------------------------------------------------------------------
# icons.txt overrides
# ---------------------------------------------------------------------------


class TestIconsFileOverride:

    def test_override_from_icons_file(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("Mammalia\tCustomMammals.png\n")
            f.flush()
            icons_path = f.name

        try:
            with _no_ncbi():
                resolver = IconResolver(icons_file=icons_path)
                meta = _FakeMeta(
                    taxonomy_lineage=["Mammalia", "Vertebrata", "Chordata"],
                )
                icon, _, _ = resolver.resolve_icon(meta)
                assert icon == "CustomMammals.png"
        finally:
            os.unlink(icons_path)


# ---------------------------------------------------------------------------
# Realistic acceptance criteria species
# ---------------------------------------------------------------------------


class TestAcceptanceCriteria:
    """Validate the specific species from the acceptance criteria,
    using metadata lineage (primary) and BUSCO fallback."""

    def test_adelges_tsugae_arthropod(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_ADELGES",
                taxonomy_lineage=["Hemiptera", "Insecta", "Arthropoda", "Metazoa"],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Arthropods.png"

    def test_aix_sponsa_bird(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_AIX",
                taxonomy_lineage=["Aves", "Vertebrata", "Chordata", "Metazoa"],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Birds.png"

    def test_antilocapra_americana_mammal(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_ANTIC",
                taxonomy_lineage=["Mammalia", "Vertebrata", "Chordata", "Metazoa"],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"

    def test_carcharodon_carcharias_fish(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_CARCH",
                taxonomy_lineage=[
                    "Chondrichthyes",
                    "Vertebrata",
                    "Chordata",
                    "Metazoa",
                ],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Fish.png"

    def test_euthyatira_pudens_lepidoptera(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_EUTHY",
                taxonomy_lineage=[
                    "Lepidoptera",
                    "Insecta",
                    "Arthropoda",
                    "Metazoa",
                ],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Lepidoptera.png"

    def test_coprinus_comatus_fungi(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_COPRI",
                taxonomy_lineage=["Basidiomycota", "Fungi"],
            )
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Fungi.png"

    def test_plant_via_busco_when_no_lineage(self):
        """Plant should resolve via BUSCO when lineage is unavailable."""
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxon_id=None,
                busco_lineage="eudicots_odb10",
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Plants.png"
            assert source == "busco_lineage"

    def test_arthropod_via_busco_when_ncbi_fails(self):
        """Arthropod should resolve via BUSCO when NCBI is down."""
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxon_id=42,
                busco_lineage="insecta_odb10",
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Arthropods.png"
            assert source == "busco_lineage"


# ---------------------------------------------------------------------------
# BUSCO plant token coverage (eudicotyledons regression)
# ---------------------------------------------------------------------------


class TestBUSCOPlantTokens:
    """Ensure plant BUSCO dataset names resolve to Plants.png."""

    def test_eudicotyledons_odb12(self):
        """The exact token that caused the Arnica montana regression."""
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxon_id=None,
                busco_lineage="eudicotyledons_odb12",
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Plants.png"
            assert source == "busco_lineage"

    def test_eudicots_odb10(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="eudicots_odb10")
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Plants.png"

    def test_poales_odb12(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="poales_odb12")
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Plants.png"

    def test_brassicales_odb12(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="brassicales_odb12")
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Plants.png"

    def test_embryophyta_odb12(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="embryophyta_odb12")
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Plants.png"

    def test_asterales_odb12(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage="asterales_odb12")
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Plants.png"

    def test_unknown_metazoan_still_fallback(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, busco_lineage=None)
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Metazoa.png"
            assert source == "missing_taxon_id"


# ---------------------------------------------------------------------------
# Lineage normalisation (string input, whitespace)
# ---------------------------------------------------------------------------


class TestLineageNormalisation:
    """Ensure lineage data in various formats is handled correctly."""

    def test_semicolon_separated_string(self):
        """Lineage as semicolon-delimited string (some DB exports)."""
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage="Lepidoptera; Insecta; Arthropoda; Metazoa",
            )
            icon, matched, source = resolver.resolve_icon(meta)
            assert icon == "Lepidoptera.png"
            assert source == "metadata_lineage"

    def test_comma_separated_string(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage="Aves,Vertebrata,Chordata",
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Birds.png"
            assert source == "metadata_lineage"

    def test_whitespace_in_list(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                taxonomy_lineage=["  Mammalia  ", "Vertebrata", "Chordata"],
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Mammals.png"
            assert source == "metadata_lineage"

    def test_empty_string_lineage(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, taxonomy_lineage="")
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Metazoa.png"

    def test_none_lineage(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(taxon_id=None, taxonomy_lineage=None)
            icon, _, _ = resolver.resolve_icon(meta)
            assert icon == "Metazoa.png"


# ---------------------------------------------------------------------------
# Arnica montana specific regression test
# ---------------------------------------------------------------------------


class TestArnicaMontanaRegression:
    """The specific regression: Arnica montana with busco_lineage=eudicotyledons_odb12."""

    def test_arnica_montana_via_busco(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_965212395.1",
                taxon_id=None,
                busco_lineage="eudicotyledons_odb12",
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Plants.png"
            assert source == "busco_lineage"

    def test_arnica_montana_via_metadata_lineage(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_965212395.1",
                taxonomy_lineage=[
                    "Asteraceae",
                    "Asterales",
                    "Magnoliopsida",
                    "Streptophyta",
                    "Viridiplantae",
                ],
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == "Plants.png"
            assert source == "metadata_lineage"


# ---------------------------------------------------------------------------
# Regression: pre-release plant must not collapse to Metazoa when NCBI is down
# ---------------------------------------------------------------------------


class TestPreReleasePlantBuscoFallback:
    """Geum rivale regression: a pre-release genome with only taxon_id and a
    busco_lineage hint (no metadata lineage) must still resolve to Plants.png
    via BUSCO when the live NCBI lookup is unavailable."""

    def test_geum_rivale_busco_fallback_when_ncbi_down(self):
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_964205265.1",
                taxon_id=334465,
                taxonomy_lineage=None,
                busco_lineage="eudicotyledons_odb12",
            )
            icon, matched, source = resolver.resolve_icon(meta)
            assert icon == "Plants.png"
            assert source == "busco_lineage"

    def test_no_busco_no_lineage_ncbi_down_is_metazoa(self):
        """Without any usable source the fallback is still Metazoa (documents
        exactly the pre-fix failure mode)."""
        with _no_ncbi():
            resolver = IconResolver(icons_file="/nonexistent")
            meta = _FakeMeta(
                accession="GCA_964205265.1",
                taxon_id=334465,
                taxonomy_lineage=None,
                busco_lineage=None,
            )
            icon, _, source = resolver.resolve_icon(meta)
            assert icon == _FALLBACK_ICON
            assert source == "fallback"

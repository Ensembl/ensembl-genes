"""Tests for the taxonomy-lineage-based icon resolver.

These tests mock NCBI responses to avoid network calls during CI.
Run with:
    pytest src/python/ensembl/genes/projects/test_icon_resolver.py -v
"""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from ensembl.genes.projects.icon_resolver import IconResolver, _FALLBACK_ICON


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


# ---------------------------------------------------------------------------
# Test the built-in default rules
# ---------------------------------------------------------------------------


class TestDefaultRules:
    """Verify that the built-in rules select the correct icon for each
    major taxonomic group, using mocked NCBI lineage data."""

    def _resolve(self, lineage_root_to_leaf, taxon_id=1):
        """Helper: create a resolver with no icons.txt and mock NCBI."""
        xml = _make_entrez_xml(lineage_root_to_leaf)
        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.return_value = _mock_response(xml)
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            return resolver.resolve_icon(taxon_id, accession="GCA_TEST")

    def test_lepidoptera(self):
        lineage = ["Eukaryota", "Metazoa", "Arthropoda", "Insecta", "Lepidoptera"]
        assert self._resolve(lineage) == "Lepidoptera.png"

    def test_non_lepidoptera_insect(self):
        """Diptera should match Arthropods.png (via Diptera rule), not Metazoa."""
        lineage = ["Eukaryota", "Metazoa", "Arthropoda", "Insecta", "Diptera"]
        assert self._resolve(lineage) == "Arthropods.png"

    def test_arthropod_without_specific_order(self):
        """A generic arthropod (e.g. a spider) should match Arthropoda rule."""
        lineage = ["Eukaryota", "Metazoa", "Arthropoda", "Arachnida"]
        assert self._resolve(lineage) == "Arthropods.png"

    def test_mammalia(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Mammalia"]
        assert self._resolve(lineage) == "Mammals.png"

    def test_aves(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Aves"]
        assert self._resolve(lineage) == "Birds.png"

    def test_teleostei(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Actinopterygii",
            "Teleostei",
        ]
        assert self._resolve(lineage) == "Fish.png"

    def test_actinopterygii_without_teleostei(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Actinopterygii"]
        assert self._resolve(lineage) == "Fish.png"

    def test_chondrichthyes(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Chondrichthyes",
        ]
        assert self._resolve(lineage) == "Fish.png"

    def test_amphibia(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Amphibia"]
        assert self._resolve(lineage) == "Amphibians.png"

    def test_testudines(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Testudines",
        ]
        assert self._resolve(lineage) == "Reptiles.png"

    def test_squamata(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Squamata",
        ]
        assert self._resolve(lineage) == "Reptiles.png"

    def test_serpentes(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Squamata",
            "Serpentes",
        ]
        assert self._resolve(lineage) == "Reptiles.png"

    def test_crocodylia(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Crocodylia",
        ]
        assert self._resolve(lineage) == "Reptiles.png"

    def test_chordate_without_specific_class(self):
        """A chordate without a more specific mapping -> Chordates.png."""
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Tunicata"]
        assert self._resolve(lineage) == "Chordates.png"

    def test_viridiplantae(self):
        lineage = ["Eukaryota", "Viridiplantae", "Streptophyta", "Magnoliopsida"]
        assert self._resolve(lineage) == "Plants.png"

    def test_fungi(self):
        lineage = ["Eukaryota", "Fungi", "Ascomycota"]
        assert self._resolve(lineage) == "Fungi.png"

    def test_generic_metazoan(self):
        """A metazoan that is not an arthropod, chordate, plant, or fungus."""
        lineage = ["Eukaryota", "Metazoa", "Cnidaria"]
        assert self._resolve(lineage) == "Metazoa.png"

    def test_completely_unknown_lineage(self):
        """Lineage with no matching entries at all -> fallback."""
        lineage = ["Eukaryota", "SomethingUnknown"]
        assert self._resolve(lineage) == _FALLBACK_ICON


# ---------------------------------------------------------------------------
# Specificity: most-specific match wins
# ---------------------------------------------------------------------------


class TestSpecificity:
    """Ensure that the most-specific taxonomy match always wins."""

    def _resolve(self, lineage_root_to_leaf, taxon_id=1):
        xml = _make_entrez_xml(lineage_root_to_leaf)
        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.return_value = _mock_response(xml)
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            return resolver.resolve_icon(taxon_id, accession="GCA_TEST")

    def test_lepidoptera_beats_arthropoda(self):
        lineage = ["Eukaryota", "Metazoa", "Arthropoda", "Insecta", "Lepidoptera"]
        assert self._resolve(lineage) == "Lepidoptera.png"

    def test_mammalia_beats_chordata(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Mammalia"]
        assert self._resolve(lineage) == "Mammals.png"

    def test_aves_beats_chordata(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Aves"]
        assert self._resolve(lineage) == "Birds.png"

    def test_teleostei_beats_chordata(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Actinopterygii",
            "Teleostei",
        ]
        assert self._resolve(lineage) == "Fish.png"

    def test_amphibia_beats_chordata(self):
        lineage = ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Amphibia"]
        assert self._resolve(lineage) == "Amphibians.png"

    def test_testudines_beats_chordata(self):
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Testudines",
        ]
        assert self._resolve(lineage) == "Reptiles.png"


# ---------------------------------------------------------------------------
# Edge cases and fallbacks
# ---------------------------------------------------------------------------


class TestEdgeCases:

    def test_missing_taxon_id(self):
        resolver = IconResolver(icons_file="/nonexistent/icons.txt")
        assert resolver.resolve_icon(None, accession="GCA_NONE") == _FALLBACK_ICON

    def test_zero_taxon_id(self):
        resolver = IconResolver(icons_file="/nonexistent/icons.txt")
        assert resolver.resolve_icon(0, accession="GCA_ZERO") == _FALLBACK_ICON

    def test_ncbi_failure(self):
        """If NCBI returns an error, fall back gracefully."""
        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.return_value = _mock_response("", status_code=500)
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            assert resolver.resolve_icon(9999, accession="GCA_ERR") == _FALLBACK_ICON

    def test_ncbi_timeout(self):
        """Network timeout -> fall back gracefully."""
        import requests as real_requests

        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.side_effect = real_requests.Timeout("timeout")
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            assert resolver.resolve_icon(9999, accession="GCA_TMO") == _FALLBACK_ICON


# ---------------------------------------------------------------------------
# Lineage caching
# ---------------------------------------------------------------------------


class TestCaching:

    def test_lineage_is_cached(self):
        """The resolver should only call NCBI once per taxon_id."""
        xml = _make_entrez_xml(["Eukaryota", "Metazoa", "Arthropoda"])
        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.return_value = _mock_response(xml)
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            icon1 = resolver.resolve_icon(42, accession="GCA_1")
            icon2 = resolver.resolve_icon(42, accession="GCA_2")
            assert icon1 == icon2 == "Arthropods.png"
            # NCBI should only have been called once
            assert mock_get.call_count == 1


# ---------------------------------------------------------------------------
# icons.txt overrides
# ---------------------------------------------------------------------------


class TestIconsFileOverride:

    def test_override_from_icons_file(self):
        """Entries in icons.txt should override built-in defaults."""
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            # Override Mammalia to use a custom icon
            f.write("Mammalia\tCustomMammals.png\n")
            f.flush()
            icons_path = f.name

        xml = _make_entrez_xml(
            ["Eukaryota", "Metazoa", "Chordata", "Vertebrata", "Mammalia"]
        )
        try:
            with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
                mock_get.return_value = _mock_response(xml)
                resolver = IconResolver(icons_file=icons_path)
                assert resolver.resolve_icon(1) == "CustomMammals.png"
        finally:
            os.unlink(icons_path)

    def test_icons_file_adds_new_group(self):
        """icons.txt can map taxonomy names that are not in the defaults."""
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("Cnidaria\tJellyfish.png\n")
            f.flush()
            icons_path = f.name

        xml = _make_entrez_xml(["Eukaryota", "Metazoa", "Cnidaria"])
        try:
            with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
                mock_get.return_value = _mock_response(xml)
                resolver = IconResolver(icons_file=icons_path)
                assert resolver.resolve_icon(1) == "Jellyfish.png"
        finally:
            os.unlink(icons_path)


# ---------------------------------------------------------------------------
# Audit output
# ---------------------------------------------------------------------------


class TestAudit:

    def test_resolve_icon_with_audit_returns_matched_taxon(self):
        xml = _make_entrez_xml(
            ["Eukaryota", "Metazoa", "Arthropoda", "Insecta", "Lepidoptera"]
        )
        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.return_value = _mock_response(xml)
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            icon, matched = resolver.resolve_icon_with_audit(1, accession="GCA_AUD")
            assert icon == "Lepidoptera.png"
            assert matched == "Lepidoptera"

    def test_audit_missing_taxon_id(self):
        resolver = IconResolver(icons_file="/nonexistent/icons.txt")
        icon, matched = resolver.resolve_icon_with_audit(None, accession="GCA_X")
        assert icon == _FALLBACK_ICON
        assert matched == "missing_taxon_id"

    def test_audit_no_match(self):
        xml = _make_entrez_xml(["Eukaryota", "SomethingUnknown"])
        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.return_value = _mock_response(xml)
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            icon, matched = resolver.resolve_icon_with_audit(1, accession="GCA_NM")
            assert icon == _FALLBACK_ICON
            assert matched == "fallback"


# ---------------------------------------------------------------------------
# Realistic species examples
# ---------------------------------------------------------------------------


class TestRealisticSpecies:
    """Validate against the examples given in the acceptance criteria."""

    def _resolve(self, lineage_root_to_leaf):
        xml = _make_entrez_xml(lineage_root_to_leaf)
        with patch("ensembl.genes.projects.icon_resolver.requests.get") as mock_get:
            mock_get.return_value = _mock_response(xml)
            resolver = IconResolver(icons_file="/nonexistent/icons.txt")
            return resolver.resolve_icon(1)

    def test_euthyatira_pudens_lepidoptera(self):
        """Euthyatira pudens is a moth (Lepidoptera)."""
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Arthropoda",
            "Insecta",
            "Pterygota",
            "Neoptera",
            "Lepidoptera",
            "Drepanidae",
        ]
        assert self._resolve(lineage) == "Lepidoptera.png"

    def test_squalus_suckleyi_fish(self):
        """Squalus suckleyi is a shark (Chondrichthyes)."""
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Chordata",
            "Vertebrata",
            "Chondrichthyes",
            "Squaliformes",
        ]
        # Sharks map to Fish.png via the Chondrichthyes rule
        assert self._resolve(lineage) == "Fish.png"

    def test_generic_insect_not_metazoa(self):
        """An insect order not explicitly in the rules should still get
        Arthropods.png via the Arthropoda catch-all, not Metazoa.png."""
        lineage = [
            "Eukaryota",
            "Metazoa",
            "Arthropoda",
            "Insecta",
            "Odonata",  # dragonflies -- not in default rules
        ]
        assert self._resolve(lineage) == "Arthropods.png"

    def test_plant_example(self):
        lineage = [
            "Eukaryota",
            "Viridiplantae",
            "Streptophyta",
            "Magnoliopsida",
            "Asteraceae",
        ]
        assert self._resolve(lineage) == "Plants.png"

    def test_unmapped_metazoan(self):
        """A nematode should fall to Metazoa.png."""
        lineage = ["Eukaryota", "Metazoa", "Nematoda"]
        assert self._resolve(lineage) == "Metazoa.png"

"""
Taxonomy-lineage-based icon resolver for Ensembl project pages.

Resolves a species icon filename using multiple taxonomy data sources,
in priority order:

1. **Metadata DB lineage** (``GenomeMetadata.taxonomy_lineage``) — pre-fetched
   from the Ensembl metadata database during data gathering.
2. **NCBI Entrez lineage** — live lookup by ``taxon_id`` (cached per run).
3. **BUSCO lineage hint** — weak fallback parsing the ``busco_lineage``
   string (e.g. ``"insecta_odb10"`` → Arthropods.png).

The first source that produces a usable lineage is walked leaf → root;
the first matching taxonomy rule determines the icon.

Icon resolution must **never** infer icons from species names, sample
descriptions, or any free-text metadata field — taxonomy only.
"""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import requests
import xmltodict

if TYPE_CHECKING:
    from ensembl.genes.projects.models import GenomeMetadata

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default taxonomy → icon rules
# ---------------------------------------------------------------------------
# The lookup dict is checked against the lineage which is already ordered
# leaf → root, so the first hit is always the most-specific classification.

_DEFAULT_RULES: Dict[str, str] = {
    # --- Arthropoda: specific orders first ---
    "Lepidoptera": "Lepidoptera.png",
    "Hymenoptera": "Arthropods.png",
    "Diptera": "Arthropods.png",
    "Coleoptera": "Arthropods.png",
    "Hemiptera": "Arthropods.png",
    "Trichoptera": "Arthropods.png",
    "Neoptera": "Arthropods.png",
    "Endopterygota": "Arthropods.png",
    "Insecta": "Arthropods.png",
    "Arachnida": "Arthropods.png",
    "Myriapoda": "Arthropods.png",
    "Crustacea": "Arthropods.png",
    "Arthropoda": "Arthropods.png",
    # --- Vertebrata: specific classes first ---
    "Mammalia": "Mammals.png",
    "Aves": "Birds.png",
    "Galloanserae": "Birds.png",
    "Neoaves": "Birds.png",
    "Amphibia": "Amphibians.png",
    # Fish groups
    "Teleostei": "Fish.png",
    "Actinopterygii": "Fish.png",
    "Chondrichthyes": "Fish.png",
    "Actinistia": "Fish.png",
    "Dipnoi": "Fish.png",
    # Reptile groups
    "Testudines": "Reptiles.png",
    "Squamata": "Reptiles.png",
    "Serpentes": "Reptiles.png",
    "Crocodylia": "Reptiles.png",
    "Lepidosauria": "Reptiles.png",
    # Broad chordate catch-all
    "Vertebrata": "Chordates.png",
    "Chordata": "Chordates.png",
    # --- Kingdoms ---
    "Viridiplantae": "Plants.png",
    "Streptophyta": "Plants.png",
    "Embryophyta": "Plants.png",
    "Fungi": "Fungi.png",
    "Ascomycota": "Fungi.png",
    "Basidiomycota": "Fungi.png",
    # --- Broadest animal catch-all ---
    "Metazoa": "Metazoa.png",
}

# ---------------------------------------------------------------------------
# BUSCO lineage → taxonomy name mapping
# ---------------------------------------------------------------------------
# busco_lineage values look like "insecta_odb10", "mammalia_odb10", etc.
# We strip the "_odb*" suffix and match against known taxonomy terms.
# The keys are lowercase BUSCO tokens; values are the canonical taxonomy
# name to look up in the rules dict.

_BUSCO_TO_TAXON: Dict[str, str] = {
    # Arthropods
    "lepidoptera": "Lepidoptera",
    "hymenoptera": "Hymenoptera",
    "diptera": "Diptera",
    "coleoptera": "Coleoptera",
    "hemiptera": "Hemiptera",
    "trichoptera": "Trichoptera",
    "endopterygota": "Endopterygota",
    "insecta": "Insecta",
    "arachnida": "Arachnida",
    "arthropoda": "Arthropoda",
    # Vertebrates
    "mammalia": "Mammalia",
    "laurasiatheria": "Mammalia",
    "euarchontoglires": "Mammalia",
    "rodentia": "Mammalia",
    "primates": "Mammalia",
    "carnivora": "Mammalia",
    "glires": "Mammalia",
    "aves": "Aves",
    "galloanserae": "Aves",
    "neoaves": "Aves",
    "passeriformes": "Aves",
    "amphibia": "Amphibia",
    "actinopterygii": "Actinopterygii",
    "teleostei": "Teleostei",
    "chondrichthyes": "Chondrichthyes",
    "cyprinodontiformes": "Actinopterygii",
    "perciformes": "Actinopterygii",
    "tetrapoda": "Vertebrata",
    "vertebrata": "Vertebrata",
    "chordata": "Chordata",
    # Reptiles
    "testudines": "Testudines",
    "squamata": "Squamata",
    "crocodylia": "Crocodylia",
    "sauropsida": "Reptiles.png",  # direct icon, rare
    # Plants
    "viridiplantae": "Viridiplantae",
    "embryophyta": "Embryophyta",
    "eudicots": "Viridiplantae",
    "poales": "Viridiplantae",
    "brassicales": "Viridiplantae",
    "fabales": "Viridiplantae",
    "solanales": "Viridiplantae",
    "liliopsida": "Viridiplantae",
    "magnoliopsida": "Viridiplantae",
    # Fungi
    "fungi": "Fungi",
    "ascomycota": "Ascomycota",
    "basidiomycota": "Basidiomycota",
    "sordariomycetes": "Fungi",
    "eurotiomycetes": "Fungi",
    "saccharomycetes": "Fungi",
    "agaricomycetes": "Fungi",
    # Metazoa
    "metazoa": "Metazoa",
    "nematoda": "Metazoa",
    "mollusca": "Metazoa",
    "eukaryota": "Metazoa",
}

# Fallback when nothing matches at all
_FALLBACK_ICON = "Metazoa.png"


class IconResolver:
    """Resolves an icon filename using taxonomy lineage from multiple sources.

    Parameters
    ----------
    icons_file:
        Path to a project-level ``icons.txt`` override file.  Each line is
        ``<TaxonomyName>\\t<icon_filename>``.  Entries here take priority
        over the built-in defaults *for the same taxonomy name*.
    """

    def __init__(
        self,
        icons_file: Optional[str] = None,
    ) -> None:
        # Build the combined lookup: defaults first, then overrides
        self._lookup: Dict[str, str] = dict(_DEFAULT_RULES)

        # Load project-specific overrides from icons.txt
        if icons_file is None:
            icons_file = os.path.join(os.path.dirname(__file__), "icons.txt")
        if os.path.exists(icons_file):
            with open(icons_file, encoding="utf-8") as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) >= 2:
                        self._lookup[parts[0]] = parts[1]

        # Per-run lineage cache: taxon_id → list of ScientificName strings
        self._lineage_cache: Dict[int, List[str]] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def resolve_icon(self, meta: "GenomeMetadata") -> Tuple[str, str, str]:
        """Return the best icon filename for a genome.

        Tries lineage sources in order:
        1. ``meta.taxonomy_lineage`` (pre-fetched from metadata DB)
        2. NCBI Entrez lookup by ``meta.taxon_id`` (cached per run)
        3. BUSCO lineage hint from ``meta.busco_lineage``

        Returns
        -------
        (icon, matched_term, lineage_source)
            - *icon*: filename e.g. ``"Lepidoptera.png"``
            - *matched_term*: taxonomy name that triggered the match
            - *lineage_source*: one of ``"metadata_lineage"``,
              ``"ncbi_lineage"``, ``"busco_lineage"``,
              ``"missing_taxon_id"``, ``"fallback"``
        """
        accession = getattr(meta, "accession", "unknown")
        taxon_id = getattr(meta, "taxon_id", None)

        # Source 1: metadata DB lineage (already on the object)
        db_lineage = getattr(meta, "taxonomy_lineage", None) or []
        if db_lineage:
            icon = self._match_lineage(db_lineage)
            if icon != _FALLBACK_ICON:
                matched = self._find_matched_term(db_lineage)
                logger.debug(
                    "Icon resolved via metadata_lineage: taxon_id=%s "
                    "accession=%s matched=%s icon=%s",
                    taxon_id,
                    accession,
                    matched,
                    icon,
                )
                return icon, matched, "metadata_lineage"

        # Source 2: NCBI Entrez live lookup (cached)
        if taxon_id:
            ncbi_lineage = self._get_ncbi_lineage(taxon_id)
            if ncbi_lineage:
                icon = self._match_lineage(ncbi_lineage)
                if icon != _FALLBACK_ICON:
                    matched = self._find_matched_term(ncbi_lineage)
                    logger.debug(
                        "Icon resolved via ncbi_lineage: taxon_id=%d "
                        "accession=%s matched=%s icon=%s",
                        taxon_id,
                        accession,
                        matched,
                        icon,
                    )
                    return icon, matched, "ncbi_lineage"
                # NCBI returned lineage but nothing matched — still try
                # busco as last resort before accepting Metazoa
        elif not taxon_id:
            logger.warning(
                "taxon_id is missing for %s. Trying busco_lineage fallback.",
                accession,
            )

        # Source 3: BUSCO lineage hint (weak but reliable)
        busco_lineage = getattr(meta, "busco_lineage", None)
        icon_from_busco = self._resolve_from_busco(busco_lineage)
        if icon_from_busco:
            matched = self._busco_matched_term(busco_lineage)
            logger.debug(
                "Icon resolved via busco_lineage: taxon_id=%s "
                "accession=%s busco=%s matched=%s icon=%s",
                taxon_id,
                accession,
                busco_lineage,
                matched,
                icon_from_busco,
            )
            return icon_from_busco, matched, "busco_lineage"

        # All sources exhausted
        if not taxon_id:
            source = "missing_taxon_id"
        else:
            source = "fallback"
        logger.warning(
            "All lineage sources exhausted for accession=%s taxon_id=%s. "
            "Using fallback %s.",
            accession,
            taxon_id,
            _FALLBACK_ICON,
        )
        return _FALLBACK_ICON, source, source

    # ------------------------------------------------------------------
    # Lineage matching
    # ------------------------------------------------------------------

    def _match_lineage(self, lineage: List[str]) -> str:
        """Walk lineage leaf → root; return first matching icon."""
        for name in lineage:
            icon = self._lookup.get(name)
            if icon is not None:
                return icon
        return _FALLBACK_ICON

    def _find_matched_term(self, lineage: List[str]) -> str:
        """Return the first lineage name that matches a rule."""
        for name in lineage:
            if name in self._lookup:
                return name
        return "fallback"

    # ------------------------------------------------------------------
    # BUSCO lineage parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _normalise_busco_token(busco_lineage: Optional[str]) -> str:
        """Extract the taxonomy token from a BUSCO lineage string.

        E.g. ``"insecta_odb10"`` → ``"insecta"``,
             ``"mammalia_odb10"`` → ``"mammalia"``.
        """
        if not busco_lineage:
            return ""
        token = busco_lineage.strip().lower()
        # Strip _odb* suffix
        if "_odb" in token:
            token = token[: token.index("_odb")]
        # Strip any remaining underscores (some are e.g. "cyprinodontiformes_odb10")
        return token.strip("_")

    def _resolve_from_busco(self, busco_lineage: Optional[str]) -> Optional[str]:
        """Try to resolve an icon from the busco_lineage string."""
        token = self._normalise_busco_token(busco_lineage)
        if not token:
            return None
        taxon_name = _BUSCO_TO_TAXON.get(token)
        if taxon_name is None:
            return None
        # taxon_name might be a direct icon filename (e.g. "Reptiles.png")
        if taxon_name.endswith(".png"):
            return taxon_name
        return self._lookup.get(taxon_name)

    def _busco_matched_term(self, busco_lineage: Optional[str]) -> str:
        """Return the taxon name derived from busco_lineage."""
        token = self._normalise_busco_token(busco_lineage)
        if not token:
            return "busco_empty"
        return _BUSCO_TO_TAXON.get(token, token)

    # ------------------------------------------------------------------
    # NCBI Entrez lineage retrieval (cached)
    # ------------------------------------------------------------------

    def _get_ncbi_lineage(self, taxon_id: int) -> List[str]:
        """Return cached taxonomy lineage (leaf → root) for *taxon_id*."""
        if taxon_id in self._lineage_cache:
            return self._lineage_cache[taxon_id]

        lineage = self._fetch_lineage_from_ncbi(taxon_id)
        self._lineage_cache[taxon_id] = lineage
        return lineage

    @staticmethod
    def _fetch_lineage_from_ncbi(taxon_id: int) -> List[str]:
        """Fetch the full taxonomy lineage from NCBI Entrez.

        Returns a list of scientific names ordered leaf → root.
        """
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=taxonomy&id={taxon_id}&retmode=xml"
        )
        try:
            res = requests.get(url, timeout=15)
            if res.status_code != 200:
                logger.debug(
                    "NCBI Entrez returned HTTP %d for taxon_id=%d",
                    res.status_code,
                    taxon_id,
                )
                return []

            data = xmltodict.parse(res.text)
            taxa_set = data.get("TaxaSet", {}).get("Taxon", {})
            if isinstance(taxa_set, list):
                taxa_set = taxa_set[0]

            lineage = taxa_set.get("LineageEx", {}).get("Taxon", [])
            if isinstance(lineage, dict):
                lineage = [lineage]

            # NCBI LineageEx is root → leaf; reverse so most-specific first
            lineage = list(reversed(lineage))
            return [t.get("ScientificName") for t in lineage if "ScientificName" in t]
        except Exception as exc:  # pylint: disable=broad-exception-caught
            logger.debug(
                "Failed to fetch taxonomy lineage for taxon_id=%d: %s",
                taxon_id,
                exc,
            )
            return []

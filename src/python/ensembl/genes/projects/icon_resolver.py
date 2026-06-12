"""
Taxonomy-lineage-based icon resolver for Ensembl project pages.

Resolves a species icon filename from the NCBI taxonomy lineage,
walking from the most specific classification to the broadest and
selecting the first matching rule.  Rules are loaded from a built-in
default table (covering all available icon files) and optionally
overridden/extended by a project-level ``icons.txt`` file.

Icon resolution must **never** infer icons from species names, sample
descriptions, or any free-text metadata field -- taxonomy only.
"""

import logging
import os
from typing import Dict, List, Optional, Tuple

import requests
import xmltodict

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default taxonomy-to-icon rules
# ---------------------------------------------------------------------------
# Ordered most-specific first.  When the lineage is walked leaf-to-root the
# *first* matching entry wins, so Lepidoptera (→ Lepidoptera.png) beats
# Arthropoda (→ Arthropods.png) for butterflies.
#
# The ordering here is a tie-breaker for when *icons.txt* contains entries
# at the same specificity level -- the built-in priority list ensures that
# e.g. "Mammalia" always resolves before "Chordata" even if both appear in
# the lineage at the same time.

_DEFAULT_RULES: List[Tuple[str, str]] = [
    # --- Arthropoda: specific orders first ---
    ("Lepidoptera", "Lepidoptera.png"),
    ("Hymenoptera", "Arthropods.png"),
    ("Diptera", "Arthropods.png"),
    ("Coleoptera", "Arthropods.png"),
    ("Hemiptera", "Arthropods.png"),
    ("Trichoptera", "Arthropods.png"),
    ("Neoptera", "Arthropods.png"),
    ("Arthropoda", "Arthropods.png"),
    # --- Vertebrata: specific classes first ---
    ("Mammalia", "Mammals.png"),
    ("Aves", "Birds.png"),
    ("Amphibia", "Amphibians.png"),
    # Fish groups
    ("Teleostei", "Fish.png"),
    ("Actinopterygii", "Fish.png"),
    ("Chondrichthyes", "Fish.png"),
    ("Actinistia", "Fish.png"),
    ("Dipnoi", "Fish.png"),
    # Reptile groups
    ("Testudines", "Reptiles.png"),
    ("Squamata", "Reptiles.png"),
    ("Serpentes", "Reptiles.png"),
    ("Crocodylia", "Reptiles.png"),
    ("Lepidosauria", "Reptiles.png"),
    # Broad chordate catch-all
    ("Vertebrata", "Chordates.png"),
    ("Chordata", "Chordates.png"),
    # --- Kingdoms ---
    ("Viridiplantae", "Plants.png"),
    ("Streptophyta", "Plants.png"),
    ("Fungi", "Fungi.png"),
    # --- Broadest animal catch-all ---
    ("Metazoa", "Metazoa.png"),
]

# Lookup dict built from the ordered list above.  Because we walk the
# lineage leaf→root and check membership in this dict, the *lineage
# ordering* already ensures most-specific-wins.  The dict itself just
# maps name → icon filename.
_DEFAULT_LOOKUP: Dict[str, str] = {name: icon for name, icon in _DEFAULT_RULES}

# Fallback when nothing matches at all
_FALLBACK_ICON = "Metazoa.png"


class IconResolver:
    """Resolves an icon filename from a taxon ID using NCBI taxonomy lineage.

    Parameters
    ----------
    icons_file:
        Path to a project-level ``icons.txt`` override file.  Each line is
        ``<TaxonomyName>\\t<icon_filename>``.  Entries here take priority
        over the built-in defaults *for the same taxonomy name*.
    taxonomy_client:
        Reserved for future injection of a non-NCBI taxonomy source.
        Currently unused; lineage is always fetched via NCBI Entrez.
    """

    def __init__(
        self,
        icons_file: Optional[str] = None,
        taxonomy_client: object = None,
    ) -> None:
        # Build the combined lookup: defaults first, then overrides
        self._lookup: Dict[str, str] = dict(_DEFAULT_LOOKUP)

        # Load project-specific overrides from icons.txt
        self._overrides: Dict[str, str] = {}
        if icons_file is None:
            icons_file = os.path.join(os.path.dirname(__file__), "icons.txt")
        if os.path.exists(icons_file):
            with open(icons_file, encoding="utf-8") as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) >= 2:
                        self._overrides[parts[0]] = parts[1]
            # Overrides win over defaults for the same taxonomy name
            self._lookup.update(self._overrides)

        # Per-run lineage cache:  taxon_id → list of ScientificName strings
        # (leaf→root order, i.e. most specific first)
        self._lineage_cache: Dict[int, List[str]] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def resolve_icon(self, taxon_id: Optional[int], accession: str = "") -> str:
        """Return the best icon filename for *taxon_id*.

        Parameters
        ----------
        taxon_id:
            NCBI taxonomy ID.  If ``None`` or ``0``, the fallback icon
            is returned immediately.
        accession:
            Assembly accession, used only for logging context.

        Returns
        -------
        str
            Icon filename (e.g. ``"Lepidoptera.png"``).
        """
        if not taxon_id:
            logger.warning(
                "taxon_id is missing for %s. Falling back to %s.",
                accession or "unknown",
                _FALLBACK_ICON,
            )
            return _FALLBACK_ICON

        lineage = self._get_lineage(taxon_id)
        if not lineage:
            logger.warning(
                "Empty lineage for taxon_id=%d (%s). Falling back to %s.",
                taxon_id,
                accession,
                _FALLBACK_ICON,
            )
            return _FALLBACK_ICON

        # Walk leaf → root; first hit in the lookup wins
        for name in lineage:
            icon = self._lookup.get(name)
            if icon is not None:
                logger.debug(
                    "Icon resolved: taxon_id=%d accession=%s " "matched=%s icon=%s",
                    taxon_id,
                    accession,
                    name,
                    icon,
                )
                return icon

        # Nothing matched at all
        logger.debug(
            "No icon rule matched for taxon_id=%d (%s). " "Using fallback %s.",
            taxon_id,
            accession,
            _FALLBACK_ICON,
        )
        return _FALLBACK_ICON

    def resolve_icon_with_audit(
        self, taxon_id: Optional[int], accession: str = ""
    ) -> Tuple[str, str]:
        """Like :meth:`resolve_icon` but also returns the matched taxon name.

        Returns
        -------
        (icon, matched_taxon)
            ``matched_taxon`` is the taxonomy name that triggered the
            match, or ``"fallback"`` / ``"missing_taxon_id"`` if no
            rule matched.
        """
        if not taxon_id:
            logger.warning(
                "taxon_id is missing for %s. Falling back to %s.",
                accession or "unknown",
                _FALLBACK_ICON,
            )
            return _FALLBACK_ICON, "missing_taxon_id"

        lineage = self._get_lineage(taxon_id)
        if not lineage:
            logger.warning(
                "Empty lineage for taxon_id=%d (%s). Falling back to %s.",
                taxon_id,
                accession,
                _FALLBACK_ICON,
            )
            return _FALLBACK_ICON, "lineage_unavailable"

        for name in lineage:
            icon = self._lookup.get(name)
            if icon is not None:
                logger.debug(
                    "Icon resolved: taxon_id=%d accession=%s " "matched=%s icon=%s",
                    taxon_id,
                    accession,
                    name,
                    icon,
                )
                return icon, name

        logger.debug(
            "No icon rule matched for taxon_id=%d (%s). " "Using fallback %s.",
            taxon_id,
            accession,
            _FALLBACK_ICON,
        )
        return _FALLBACK_ICON, "fallback"

    # ------------------------------------------------------------------
    # Lineage retrieval (cached)
    # ------------------------------------------------------------------

    def _get_lineage(self, taxon_id: int) -> List[str]:
        """Return cached taxonomy lineage (leaf → root) for *taxon_id*."""
        if taxon_id in self._lineage_cache:
            return self._lineage_cache[taxon_id]

        lineage = self._fetch_lineage_from_ncbi(taxon_id)
        self._lineage_cache[taxon_id] = lineage
        return lineage

    @staticmethod
    def _fetch_lineage_from_ncbi(taxon_id: int) -> List[str]:
        """Fetch the full taxonomy lineage from NCBI Entrez.

        Returns a list of scientific names ordered leaf → root (most
        specific classification first).
        """
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=taxonomy&id={taxon_id}&retmode=xml"
        )
        try:
            res = requests.get(url, timeout=10)
            if res.status_code != 200:
                logger.warning(
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

            # NCBI LineageEx is root → leaf; reverse so most-specific is first
            lineage = list(reversed(lineage))
            return [t.get("ScientificName") for t in lineage if "ScientificName" in t]
        except Exception as exc:  # pylint: disable=broad-exception-caught
            logger.warning(
                "Failed to fetch taxonomy lineage for taxon_id=%d: %s",
                taxon_id,
                exc,
            )
            return []

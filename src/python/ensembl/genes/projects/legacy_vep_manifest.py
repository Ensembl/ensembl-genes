"""
Legacy VEP manifest reader for the old Ensembl organisms FTP structure.

The legacy manifest at:

    https://ftp.ebi.ac.uk/pub/ensemblorganisms/species.json

uses species-name-based directory paths and is no longer the primary source of
data.  This module reads it **only** to resolve VEP annotation file paths for
projects (currently only HPRC) where the new accession-based manifest
(``species.new_ftp_structure.json``) does not yet include VEP entries.

All other FTP data — genesets, genome sequences, homology, variation — continues
to be resolved exclusively from ``EnsemblFtpManifest``.

The ``species.json`` manifest stores paths to individual VEP files
(e.g. ``genes.gff3.bgz``).  The project page links to the *directory* that
contains those files, not to the file itself.  The directory is derived from
the primary VEP file path via :attr:`LegacyVepRecord.directory_path`.

Usage::

    # Production (once per run, injected into YamlRenderer):
    try:
        legacy = LegacyVepManifest.from_url()
    except LegacyVepManifestError as exc:
        logger.warning("Legacy VEP manifest unavailable: %s", exc)
        legacy = None

    # Tests (inject fixture data directly):
    legacy = LegacyVepManifest(fixture_data)
    record = legacy.lookup_vep("GCA_018852605.1", provider="ensembl")
    # Use record.directory_path (relative) or the full URL:
    # _manifest_url(EBI_FTP_BASE, record.directory_path)
"""

import logging
import posixpath
from dataclasses import dataclass, field
from typing import Any

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LEGACY_MANIFEST_URL = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/species.json"
EBI_FTP_BASE = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"

# The primary VEP file in every record observed in the live manifest.
# The companion index file (.csi) is not emitted in YAML.
VEP_PRIMARY_FILE = "genes.gff3.bgz"

# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class LegacyVepManifestError(RuntimeError):
    """Raised for any legacy manifest download, parse or validation failure."""


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------


@dataclass
class LegacyVepRecord:
    """One VEP annotation record extracted from the legacy manifest.

    ``accession``
        The assembly accession (e.g. ``"GCA_018852605.1"``).

    ``provider``
        The genebuild provider key (e.g. ``"ensembl"``).

    ``date_key``
        The annotation date key (e.g. ``"2022_07"``).

    ``vep_relative_path``
        The relative FTP path for the primary VEP file
        (``genes.gff3.bgz``), exactly as stored in the manifest.
        Prepend :data:`EBI_FTP_BASE` to form the full file URL, or use
        :attr:`directory_path` to obtain the containing release directory.

    ``all_vep_files``
        All VEP file entries for this provider/date, keyed by filename.
        Relative paths, as supplied by the manifest.
    """

    accession: str
    provider: str
    date_key: str
    vep_relative_path: str  # relative path for genes.gff3.bgz
    all_vep_files: dict[str, str] = field(default_factory=dict)

    @property
    def directory_path(self) -> str:
        """Relative FTP path for the dated VEP release directory.

        Derived from :attr:`vep_relative_path` by stripping the filename and
        appending a trailing slash.  For example::

            # vep_relative_path:
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/genes.gff3.bgz"
            # directory_path:
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/"

        Raises:
            LegacyVepManifestError: If the path has no parent directory
                (i.e. it is a bare filename with no directory component).
        """
        parent = posixpath.dirname(self.vep_relative_path.rstrip("/"))
        if not parent:
            raise LegacyVepManifestError(
                f"Cannot derive VEP directory from manifest path "
                f"{self.vep_relative_path!r}: no parent directory."
            )
        return f"{parent}/"

# ---------------------------------------------------------------------------
# URL joining helper
# ---------------------------------------------------------------------------


def _manifest_url(base_url: str, path: str) -> str:
    """Join a manifest relative path with a base URL.

    Handles any combination of leading/trailing slashes so that the result
    always has exactly one ``/`` between base and path.

    Examples::

        _manifest_url(
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/",
            "Homo_sapiens/GCA_018852605.1/vep/ensembl/geneset/2022_07/genes.gff3.bgz"
        )
        # -> "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/GCA_018852605.1/..."

    Args:
        base_url: Base URL, with or without a trailing slash.
        path: Relative path, with or without a leading slash.

    Returns:
        The joined URL as a string.
    """
    return f"{base_url.rstrip('/')}/{path.lstrip('/')}"


# ---------------------------------------------------------------------------
# Manifest reader
# ---------------------------------------------------------------------------


class LegacyVepManifest:
    """In-memory index of the legacy Ensembl organisms VEP manifest.

    The full ``species.json`` manifest is large (~16 MB) so this class indexes
    only the VEP-relevant fields and discards the rest.

    Use :meth:`from_url` in production.  Pass a pre-loaded dict in tests::

        manifest = LegacyVepManifest(fixture_data)
        record = manifest.lookup_vep("GCA_018852605.1", provider="ensembl")
    """

    MANIFEST_URL: str = LEGACY_MANIFEST_URL
    EBI_FTP_BASE: str = EBI_FTP_BASE

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, data: dict[str, Any]) -> None:
        """Build the in-memory VEP index from a parsed legacy manifest dict.

        The index maps ``accession -> list[LegacyVepRecord]``.  Multiple records
        per accession are allowed (different providers/dates).  Ambiguity is
        resolved at lookup time, not at build time.

        Args:
            data: The top-level manifest JSON object; must contain a
                ``"species"`` key.

        Raises:
            LegacyVepManifestError: If *data* is missing required keys or has
                an incompatible structure.
        """
        if not isinstance(data, dict) or "species" not in data:
            raise LegacyVepManifestError(
                "Legacy manifest JSON is missing the required top-level 'species' key."
            )

        # accession -> list of VEP records (one per provider/date)
        self._index: dict[str, list[LegacyVepRecord]] = {}

        try:
            self._build_index(data["species"])
        except LegacyVepManifestError:
            raise
        except Exception as exc:  # pylint: disable=broad-exception-caught
            raise LegacyVepManifestError(
                f"Failed to build legacy VEP manifest index: {exc}"
            ) from exc

        assemblies_with_vep = sum(1 for recs in self._index.values() if recs)
        logger.info(
            "Legacy VEP manifest loaded: %d assemblies with VEP records.",
            assemblies_with_vep,
        )

    @classmethod
    def from_url(cls, timeout: int = 30) -> "LegacyVepManifest":
        """Download the legacy manifest and return an instance.

        All network and JSON failures are wrapped in :exc:`LegacyVepManifestError`.

        Args:
            timeout: HTTP request timeout in seconds.

        Raises:
            LegacyVepManifestError: On download, parse or validation failure.
        """
        try:
            response = requests.get(cls.MANIFEST_URL, timeout=timeout)
            response.raise_for_status()
        except requests.exceptions.Timeout as exc:
            raise LegacyVepManifestError(
                f"Timed out downloading legacy VEP manifest from {cls.MANIFEST_URL}"
            ) from exc
        except requests.exceptions.RequestException as exc:
            raise LegacyVepManifestError(
                f"Could not download legacy VEP manifest from {cls.MANIFEST_URL}: {exc}"
            ) from exc

        try:
            data = response.json()
        except ValueError as exc:
            raise LegacyVepManifestError(
                f"Legacy VEP manifest response is not valid JSON: {exc}"
            ) from exc

        return cls(data)

    # ------------------------------------------------------------------
    # Index building
    # ------------------------------------------------------------------

    def _build_index(self, species_dict: dict[str, Any]) -> None:
        """Populate ``self._index`` from the manifest species dict.

        Only VEP-bearing records are stored; assemblies without any VEP data
        are skipped silently because the caller will treat ``None`` from
        :meth:`lookup_vep` as "not found".
        """
        for _, species_data in species_dict.items():
            assemblies = species_data.get("assemblies", {})
            if not isinstance(assemblies, dict):
                continue

            for accession, asm_data in assemblies.items():
                providers_data = asm_data.get("genebuild_providers", {})
                if not isinstance(providers_data, dict):
                    continue

                for provider, dates_data in providers_data.items():
                    if not isinstance(dates_data, dict):
                        continue
                    for date_key, date_data in dates_data.items():
                        vep_files = self._extract_vep_files(date_data)
                        if not vep_files:
                            continue  # no VEP data for this provider/date

                        primary = vep_files.get(VEP_PRIMARY_FILE)
                        if not primary:
                            logger.debug(
                                "VEP record for %s/%s/%s has no %r file; skipping.",
                                accession,
                                provider,
                                date_key,
                                VEP_PRIMARY_FILE,
                            )
                            continue

                        record = LegacyVepRecord(
                            accession=accession,
                            provider=provider,
                            date_key=date_key,
                            vep_relative_path=primary,
                            all_vep_files=dict(vep_files),
                        )
                        self._index.setdefault(accession, []).append(record)

    @staticmethod
    def _extract_vep_files(date_data: dict[str, Any]) -> dict[str, str]:
        """Extract the ``vep`` file dict from one provider-date entry.

        The VEP dict lives at::

            paths.genebuild.files.vep

        Returns an empty dict if any intermediate key is missing.
        """
        vep = (
            date_data.get("paths", {})
            .get("genebuild", {})
            .get("files", {})
            .get("vep", {})
        )
        return vep if isinstance(vep, dict) else {}

    # ------------------------------------------------------------------
    # Public lookup
    # ------------------------------------------------------------------

    def lookup_vep(
        self,
        accession: str,
        *,
        provider: str | None = None,
        annotation_date: str | None = None,
    ) -> LegacyVepRecord | None:
        """Look up a VEP record for *accession* from the legacy manifest.

        Resolution order when multiple records exist for an accession:

        1. Exact provider match (normalised lowercase).
        2. Exact or prefix-normalised annotation-date match
           (``"2022-07"`` matches ``"2022_07"``; only the first two components
           are compared so ``"2022-07-15"`` also matches ``"2022_07"``).
        3. Only one record remains after filtering → return it.
        4. Multiple records remain → return ``None`` with an ``ambiguous``
           indication (caller should log the ambiguity).

        Args:
            accession: Assembly accession, e.g. ``"GCA_018852605.1"``.
            provider: Provider key from the new manifest (used as a hint; not
                mandatory).
            annotation_date: Annotation date string from metadata, in any of
                ``YYYY_MM``, ``YYYY-MM``, or ``YYYY-MM-DD`` formats.

        Returns:
            A :class:`LegacyVepRecord` if exactly one match is found, or
            ``None`` if the accession is absent, has no VEP data, or is
            ambiguous.
        """
        candidates = list(self._index.get(accession, []))
        if not candidates:
            return None

        if len(candidates) == 1:
            return candidates[0]

        # Step 1: filter by provider
        if provider:
            prov_lower = provider.lower()
            by_provider = [c for c in candidates if c.provider.lower() == prov_lower]
            if by_provider:
                candidates = by_provider
            if len(candidates) == 1:
                return candidates[0]

        # Step 2: filter by date prefix (YYYY_MM)
        if annotation_date:
            date_prefix = annotation_date.replace("-", "_")[:7]  # "YYYY_MM"
            by_date = [
                c for c in candidates if c.date_key[:7] == date_prefix
            ]
            if by_date:
                candidates = by_date
            if len(candidates) == 1:
                return candidates[0]

        # Multiple records remain — ambiguous
        logger.info(
            "Ambiguous legacy VEP records for %s: %s candidates. "
            "Providers/dates: %s. No VEP URL will be emitted.",
            accession,
            len(candidates),
            [(c.provider, c.date_key) for c in candidates],
        )
        return None

    def is_ambiguous(
        self,
        accession: str,
        *,
        provider: str | None = None,
        annotation_date: str | None = None,
    ) -> bool:
        """Return ``True`` if *accession* has multiple unresolvable VEP records.

        This is a convenience method for building audit messages.  It mirrors
        the resolution logic in :meth:`lookup_vep` but returns a boolean rather
        than the resolved record.
        """
        candidates = list(self._index.get(accession, []))
        if len(candidates) <= 1:
            return False

        if provider:
            by_provider = [
                c for c in candidates if c.provider.lower() == provider.lower()
            ]
            if by_provider:
                candidates = by_provider
        if len(candidates) <= 1:
            return False

        if annotation_date:
            date_prefix = annotation_date.replace("-", "_")[:7]
            by_date = [c for c in candidates if c.date_key[:7] == date_prefix]
            if by_date:
                candidates = by_date
        return len(candidates) > 1

    def __len__(self) -> int:
        return sum(len(recs) for recs in self._index.values())

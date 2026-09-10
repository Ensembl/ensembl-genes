"""
Manifest reader for the new Ensembl accession-based FTP structure.

The species manifest at:

    https://ftp.ebi.ac.uk/pub/ensemblorganisms/species.json

is the authoritative index of what data is available on the new FTP site.
This module downloads it once, parses it, and exposes an efficient lookup API.

Usage::

    # Production — download once per run:
    try:
        manifest = EnsemblFtpManifest.from_url()
    except ManifestError as exc:
        logger.error("Manifest unavailable: %s", exc)
        manifest = None

    # Tests — inject fixture data directly:
    manifest = EnsemblFtpManifest(data)
    record = manifest.lookup("GCA_922984935.2")
"""

import logging
import re
from dataclasses import dataclass, field
from typing import Any

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MANIFEST_URL = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/species.json"
EBI_FTP_BASE = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"

# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class ManifestError(RuntimeError):
    """Raised for any manifest download, parse or validation failure."""


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class ProviderDateRecord:
    """One annotation release for a given assembly, provider, and date.

    ``date_key`` is the opaque manifest directory name — normally ``YYYY_MM``
    (e.g. ``"2022_11"``), but the format is not guaranteed.

    All file dictionaries map a filename (e.g. ``"genes.gtf.gz"``) to a
    *relative* FTP path such as
    ``"GCA/922/984/935/2/ensembl/2022_11/geneset/genes.gtf.gz"``.

    Genome/assembly files are **not** stored here; they belong to
    :class:`AssemblyRecord`.
    """

    date_key: str
    annotation_files: dict[str, str] = field(default_factory=dict)
    homology_files: dict[str, str] = field(default_factory=dict)
    variation_files: dict[str, str] = field(default_factory=dict)


@dataclass
class VariationRecord:
    """One candidate variation VCF release for an assembly.

    ``rel_path`` is the exact relative FTP path from the manifest
    (e.g. ``"GCA/018/469/665/2/ensembl/2025_08/variation/2025_12_15/variation.vcf.gz"``).

    ``variation_date`` is the variation release date string extracted from
    the path (e.g. ``"2025_12_15"``).

    ``parsed_date`` is the comparable integer tuple for sorting (e.g. ``(2025, 12, 15)``).

    ``genebuild_date`` is the containing genebuild/annotation date key (e.g. ``"2025_08"``).
    """

    rel_path: str
    variation_date: str
    parsed_date: tuple[int, ...]
    genebuild_date: str


@dataclass
class AssemblyRecord:
    """One assembly entry from the species manifest.

    ``assembly_genome_files`` holds genome-sequence files from the manifest's
    ``assembly.files.genome_sequences`` section (e.g. ``"softmasked.fa.bgz"``).
    These are assembly-level data, independent of any particular annotation
    provider or date.

    ``providers`` maps provider name → date_key → :class:`ProviderDateRecord`.
    """

    accession: str
    species_key: str  # top-level species name, used in error messages
    assembly_genome_files: dict[str, str] = field(default_factory=dict)
    providers: dict[str, dict[str, "ProviderDateRecord"]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Date normalisation helper
# ---------------------------------------------------------------------------

_DATE_RE = re.compile(r"^(\d{4})[_\-](\d{2})(?:[_\-](\d{2}))?$")


def _parse_manifest_date(date_str: str) -> tuple[int, ...] | None:
    """Normalise a date string to a comparable integer tuple.

    Supported input formats:
        ``YYYY_MM``, ``YYYY-MM``, ``YYYY_MM_DD``, ``YYYY-MM-DD``

    Returns a tuple such as ``(2022, 11)`` or ``(2023, 10, 18)``, or
    ``None`` if the string cannot be parsed or contains an out-of-range
    month/day value.

    Examples::

        _parse_manifest_date("2022_11")      # -> (2022, 11)
        _parse_manifest_date("2023-10-18")   # -> (2023, 10, 18)
        _parse_manifest_date("2024_99")      # -> None  (invalid month)
        _parse_manifest_date("bad")          # -> None
    """
    m = _DATE_RE.fullmatch(date_str)
    if not m:
        return None
    year, month = int(m.group(1)), int(m.group(2))
    if not 1 <= month <= 12:
        return None
    if m.group(3):
        day = int(m.group(3))
        if not 1 <= day <= 31:
            return None
        return (year, month, day)
    return (year, month)


_VARIATION_DATE_RE = re.compile(r"/variation/(\d{4}[_\-]\d{2}(?:[_\-]\d{2})?)")


def extract_variation_date_from_path(rel_path: str) -> str:
    """Extract the variation release date string from a variation file path.

    The variation release date lives inside the path after ``/variation/``, e.g.::

        extract_variation_date_from_path(
            "GCA/018/469/665/2/ensembl/2025_08/variation/2025_12_15/variation.vcf.gz"
        )
        # -> "2025_12_15"

        extract_variation_date_from_path(
            "GCA/018/506/975/2/ensembl/2025_08/variation/2025_12_15/2025_12_15/variation.vcf.gz"
        )
        # -> "2025_12_15"

    Returns an empty string if no standard date segment follows ``/variation/``.
    """
    m = _VARIATION_DATE_RE.search(rel_path)
    if m:
        return m.group(1)
    return ""


# ---------------------------------------------------------------------------
# Manifest reader
# ---------------------------------------------------------------------------


class EnsemblFtpManifest:
    """In-memory index of the new Ensembl FTP species manifest.

    One instance is created per generator run and shared across all genome
    lookups.  There is no on-disk cache; this object is held in memory for the
    duration of the process only.

    Use :meth:`from_url` in production and pass a pre-loaded dict in tests::

        manifest = EnsemblFtpManifest(fixture_data)
    """

    MANIFEST_URL: str = MANIFEST_URL
    EBI_FTP_BASE: str = EBI_FTP_BASE

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, data: dict[str, Any]) -> None:
        """Build the in-memory index from a parsed manifest dict.

        Args:
            data: The top-level manifest JSON object, which must contain a
                ``"species"`` key.

        Raises:
            ManifestError: If *data* is missing required keys, has an
                incompatible structure, or contains a duplicate accession.
        """
        if not isinstance(data, dict) or "species" not in data:
            raise ManifestError(
                "Manifest JSON is missing the required top-level 'species' key."
            )

        self._index: dict[str, AssemblyRecord] = {}
        try:
            self._build_index(data["species"])
        except ManifestError:
            raise
        except Exception as exc:  # pylint: disable=broad-exception-caught
            raise ManifestError(f"Failed to build manifest index: {exc}") from exc

        logger.info(
            "Manifest loaded: %d assemblies across %d species entries.",
            len(self._index),
            len(data["species"]),
        )

    @classmethod
    def from_url(cls, timeout: int = 30) -> "EnsemblFtpManifest":
        """Download the manifest from the canonical URL and return an instance.

        All network and JSON failures are wrapped in :exc:`ManifestError` so
        that callers need only catch one exception type.

        Args:
            timeout: HTTP request timeout in seconds.

        Raises:
            ManifestError: On any download, parse or validation failure.
        """
        try:
            response = requests.get(cls.MANIFEST_URL, timeout=timeout)
            response.raise_for_status()
        except requests.exceptions.Timeout as exc:
            raise ManifestError(
                f"Timed out downloading FTP manifest from {cls.MANIFEST_URL}"
            ) from exc
        except requests.exceptions.RequestException as exc:
            raise ManifestError(
                f"Could not download FTP manifest from {cls.MANIFEST_URL}: {exc}"
            ) from exc

        try:
            data = response.json()
        except ValueError as exc:
            raise ManifestError(
                f"FTP manifest response is not valid JSON: {exc}"
            ) from exc

        return cls(data)

    # ------------------------------------------------------------------
    # Index building
    # ------------------------------------------------------------------

    def _build_index(self, species_dict: dict[str, Any]) -> None:
        """Populate ``self._index`` from the manifest species dict."""
        for species_key, species_data in species_dict.items():
            assemblies = species_data.get("assemblies", {})
            if not isinstance(assemblies, dict):
                logger.warning(
                    "Species %r has unexpected 'assemblies' value; skipping.",
                    species_key,
                )
                continue

            for accession, asm_data in assemblies.items():
                if accession in self._index:
                    raise ManifestError(
                        f"Duplicate accession {accession!r} in manifest "
                        f"(already seen under species {self._index[accession].species_key!r}, "
                        f"now also under {species_key!r})."
                    )

                record = AssemblyRecord(
                    accession=accession,
                    species_key=species_key,
                )

                # Assembly-level genome files (softmasked, hardmasked, etc.)
                genome_seqs = (
                    asm_data.get("assembly", {})
                    .get("files", {})
                    .get("genome_sequences", {})
                )
                if isinstance(genome_seqs, dict):
                    record.assembly_genome_files = dict(genome_seqs)

                # Annotation providers → dates → files
                providers_data = asm_data.get("genebuild_providers", {})
                if isinstance(providers_data, dict):
                    for provider, dates_data in providers_data.items():
                        if not isinstance(dates_data, dict):
                            continue
                        record.providers[provider] = {}
                        for date_key, date_data in dates_data.items():
                            pdr = self._parse_provider_date_record(date_key, date_data)
                            record.providers[provider][date_key] = pdr

                self._index[accession] = record

    @staticmethod
    def _parse_provider_date_record(
        date_key: str, date_data: dict[str, Any]
    ) -> ProviderDateRecord:
        """Extract annotation, homology and variation file dicts from one date entry."""
        paths = date_data.get("paths", {})

        # Annotation files (geneset/)
        annotation_files: dict[str, str] = {}
        genebuild = paths.get("genebuild", {})
        ann = genebuild.get("files", {}).get("annotations", {})
        if isinstance(ann, dict):
            annotation_files = dict(ann)

        # Homology files
        homology_files: dict[str, str] = {}
        hom = paths.get("homologies", {}).get("files", {}).get("homology_data", {})
        if isinstance(hom, dict):
            homology_files = dict(hom)

        # Variation files
        variation_files: dict[str, str] = {}
        var = paths.get("short_variants", {}).get("files", {}).get("variation_data", {})
        if isinstance(var, dict):
            variation_files = dict(var)

        return ProviderDateRecord(
            date_key=date_key,
            annotation_files=annotation_files,
            homology_files=homology_files,
            variation_files=variation_files,
        )

    # ------------------------------------------------------------------
    # Public lookup
    # ------------------------------------------------------------------

    def lookup(self, accession: str) -> AssemblyRecord | None:
        """Look up an assembly by accession.

        Args:
            accession: An INSDC accession string, e.g. ``"GCA_922984935.2"``.

        Returns:
            The :class:`AssemblyRecord` for *accession*, or ``None`` if the
            accession is not present in the manifest.
        """
        return self._index.get(accession)

    def candidate_variation_vcfs(  # pylint: disable=too-many-locals
        self,
        accession: str,
        *,
        provider: str | None = None,
    ) -> list[VariationRecord]:
        """Return all usable variation VCF candidate records for *accession*,
        sorted newest variation release date first.

        Resolution rules:
        1. Look up accession. Return empty list if absent.
        2. Select provider using standard priority (exact match -> single provider -> 'ensembl').
        3. Enumerate all date records for that provider.
        4. For records with 'variation.vcf.gz', extract the variation release date from the path.
        5. Parse the date into a comparable tuple for sorting.
        6. Sort candidates by parsed variation date (newest first), with deterministic tie-breaking.
        """
        record = self._index.get(accession)
        if record is None:
            return []

        providers = record.providers
        if not providers:
            return []

        selected_provider: str | None = None
        if provider:
            prov_lower = provider.lower()
            if prov_lower in providers:
                selected_provider = prov_lower
        if selected_provider is None:
            if len(providers) == 1:
                selected_provider = next(iter(providers))
            elif "ensembl" in providers:
                selected_provider = "ensembl"
            else:
                return []

        provider_dates = providers[selected_provider]
        candidates: list[VariationRecord] = []
        for gb_date_key, pdr in provider_dates.items():
            vcf_path = pdr.variation_files.get("variation.vcf.gz")
            if not vcf_path:
                continue
            var_date = extract_variation_date_from_path(vcf_path)
            parsed = _parse_manifest_date(var_date) if var_date else None
            if parsed is None:
                # Fall back to genebuild date if variation date could not be parsed
                parsed = _parse_manifest_date(gb_date_key)
            sort_key = parsed if parsed is not None else (0,)
            candidates.append(
                VariationRecord(
                    rel_path=vcf_path,
                    variation_date=var_date or gb_date_key,
                    parsed_date=sort_key,
                    genebuild_date=gb_date_key,
                )
            )

        # Sort newest variation release first; tie-break deterministically
        candidates.sort(
            key=lambda c: (
                c.parsed_date,
                c.variation_date,
                c.genebuild_date,
                c.rel_path,
            ),
            reverse=True,
        )
        return candidates

    def lookup_variation_vcf(
        self,
        accession: str,
        *,
        provider: str | None = None,
    ) -> tuple[str, str] | None:
        """Return the manifest-provided relative path and variation release date
        for the newest usable variation VCF for *accession*, or None.

        The returned path is the *exact* value from the manifest (relative,
        no leading slash).  Callers must prepend :data:`EBI_FTP_BASE` and
        derive the containing directory via :func:`posixpath.dirname`.

        Args:
            accession: INSDC accession string.
            provider:  Optional provider hint (e.g. ``"ensembl"``).

        Returns:
            ``(relative_vcf_path, variation_date)`` or ``None``.
        """
        candidates = self.candidate_variation_vcfs(accession, provider=provider)
        if not candidates:
            return None
        best = candidates[0]
        return best.rel_path, best.variation_date

    def __len__(self) -> int:
        return len(self._index)

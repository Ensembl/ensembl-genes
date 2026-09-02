"""
Converts internal GenomeMetadata objects into specific project YAML schemas.
"""

import logging
import re
from typing import Any, Dict, List, Optional


from ensembl.genes.projects.accession_utils import accession_to_ftp_path
from ensembl.genes.projects.config import ProjectConfig
from ensembl.genes.projects.ftp_client import (
    check_beta_species_status,
    check_url_status,
)
from ensembl.genes.projects.ftp_manifest import (
    EBI_FTP_BASE,
    AssemblyRecord,
    EnsemblFtpManifest,
    ProviderDateRecord,
    _parse_manifest_date,
)
from ensembl.genes.projects.icon_resolver import IconResolver
from ensembl.genes.projects.legacy_vep_manifest import (
    LegacyVepManifest,
    _manifest_url,
)
from ensembl.genes.projects.models import GenomeMetadata

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Provider alias map: lower-cased metadata annotation_source -> manifest key.
# Only add entries that are verified against real metadata + manifest values.
# ---------------------------------------------------------------------------
_PROVIDER_ALIASES: dict[str, str] = {
    "braker2": "braker",
    "ensembl_braker": "braker",
}


class YamlRenderer:  # pylint: disable=too-few-public-methods
    """Renders GenomeMetadata into validated dictionary structures for YAML output."""

    def __init__(
        self,
        config: ProjectConfig,
        ftp_client=None,
        manifest: Optional[EnsemblFtpManifest] = None,
        legacy_vep_manifest: Optional[LegacyVepManifest] = None,
    ):
        self.config = config
        self.ftp_client = ftp_client
        self.manifest = manifest
        self.legacy_vep_manifest = legacy_vep_manifest
        self.icon_resolver = IconResolver()
        # Per-run cache of beta species availability, keyed by genome UUID.
        self._beta_status_cache: Dict[str, str] = {}
        # Per-run cache for FTP species name resolution.
        self._ftp_species_cache: Dict[str, str] = {}

    def _resolve_vep_url(
        self,
        meta: GenomeMetadata,
        ftp_assets: Dict[str, Any],
    ) -> tuple[str | None, str]:
        """Resolve a VEP annotation file URL for an HPRC genome.

        Resolution order:

        A. New manifest — reserved for a future ``vep`` section in
           ``species.new_ftp_structure.json``.  Currently that manifest does
           not expose VEP paths so this step is a no-op placeholder.

        B. Legacy manifest (``species.json``) — used when
           ``config.use_legacy_vep_fallback`` is True and the legacy manifest
           was successfully loaded.

           The lookup uses the exact assembly accession.  Provider and
           annotation date are used as tie-breakers when multiple records
           exist for the same accession.

           If the resolved path is a relative path it is joined with
           :data:`~legacy_vep_manifest.EBI_FTP_BASE`.  The resulting URL is
           validated via :func:`~ftp_client.check_url_status` before it is
           emitted.

        C. No result — ``variants_vep`` is omitted.

        Args:
            meta: Genome metadata for the current record.
            ftp_assets: The dict returned by ``_resolve_ftp_assets()``; used
                to extract the selected provider and resolved date.

        Returns:
            ``(url_or_none, audit_status)`` where *audit_status* is one of:

            * ``available_new_manifest`` — not currently reachable
            * ``available_legacy_manifest`` — directory URL from ``species.json``, validated
            * ``not_found`` — neither manifest contained a VEP path
            * ``legacy_manifest_unavailable`` — fallback enabled but manifest absent
            * ``legacy_url_unavailable`` — legacy directory found but URL check failed
            * ``ambiguous_legacy_record`` — multiple unresolvable VEP records
            * ``not_applicable`` — VEP fallback disabled for this project

        The emitted URL points to the dated release *directory* (e.g.
        ``.../vep/ensembl/geneset/2022_07/``), not to the
        ``genes.gff3.bgz`` file directly.
        """
        # Step A: new accession-based manifest.
        # The current new manifest schema has no VEP section; reserved for
        # when it is added in the future.  Do not speculate.
        # (If a "vep" key appears in ftp_assets in the future, handle it here.)

        # VEP fallback disabled for this project.
        if not self.config.use_legacy_vep_fallback:
            return None, "not_applicable"

        # Step B: legacy manifest.
        if self.legacy_vep_manifest is None:
            return None, "legacy_manifest_unavailable"

        resolved_provider = ftp_assets.get("resolved_provider")
        resolved_date = ftp_assets.get("resolved_date")

        # Check for ambiguity before attempting resolution.
        if self.legacy_vep_manifest.is_ambiguous(
            meta.accession,
            provider=resolved_provider,
            annotation_date=resolved_date,
        ):
            providers_dates = [
                f"{p}/{d}"
                for p, d in self.legacy_vep_manifest.candidate_provider_dates(
                    meta.accession
                )
            ]
            logger.info(
                "Ambiguous legacy VEP records for %s: %s. No VEP URL emitted.",
                meta.accession,
                providers_dates,
            )
            return None, "ambiguous_legacy_record"

        legacy_record = self.legacy_vep_manifest.lookup_vep(
            meta.accession,
            provider=resolved_provider,
            annotation_date=resolved_date,
        )

        if legacy_record is None:
            return None, "not_found"

        # Derive the dated release directory from the primary VEP file path.
        # The manifest stores a direct file path (genes.gff3.bgz); the project
        # page links to the containing directory so users can browse or download
        # any of the VEP files (the .csi companion index is also present there).
        vep_dir_path = legacy_record.directory_path
        if vep_dir_path.startswith(("http://", "https://")):
            # Guard for future schema changes where the manifest might supply
            # a full URL rather than a relative path.
            vep_directory_url = vep_dir_path
        else:
            vep_directory_url = _manifest_url(
                LegacyVepManifest.EBI_FTP_BASE, vep_dir_path
            )

        # Validate the directory URL before emitting.
        if not check_url_status(vep_directory_url):
            logger.info(
                "Legacy VEP directory for %s is not reachable: %s",
                meta.accession,
                vep_directory_url,
            )
            return None, "legacy_url_unavailable"

        return vep_directory_url, "available_legacy_manifest"

    def _check_beta_status(self, genome_uuid: str) -> str:
        """Cached wrapper around ``check_beta_species_status``."""
        if genome_uuid in self._beta_status_cache:
            return self._beta_status_cache[genome_uuid]
        status = check_beta_species_status(genome_uuid)
        self._beta_status_cache[genome_uuid] = status
        return status

    def _resolve_beta_link(self, meta: GenomeMetadata, target_released: bool):
        """Decide the ``beta_link`` value and its audit status.

        Only released genomes with a confirmed-usable beta species page get a
        real beta URL; everything else gets "Coming soon!". This avoids
        emitting links that resolve to the beta "species not recognised" page.

        Returns:
            Tuple[str, str]: (beta_link_value, beta_link_status) where status
            is one of: available, unavailable, error, skipped_prerelease,
            skipped_no_uuid.
        """
        if not target_released:
            return "Coming soon!", "skipped_prerelease"
        if not meta.genome_uuid or meta.genome_uuid == "unknown":
            return "Coming soon!", "skipped_no_uuid"
        status = self._check_beta_status(meta.genome_uuid)
        if status == "available":
            return (
                f"https://ensembl.org/species/{meta.genome_uuid}",
                "available",
            )
        # "unavailable" or "error" — do not emit a possibly-broken link.
        return "Coming soon!", status

    def render(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Dispatches to the correct schema renderer based on project config."""
        if self.config.schema_type == "hprc":
            return self._render_hprc(meta)
        if self.config.schema_type == "mouse":
            return self._render_mouse(meta)
        return self._render_standard(meta)

    # ------------------------------------------------------------------
    # Provider / date selection helpers (manifest-based)
    # ------------------------------------------------------------------

    def _resolve_provider(
        self,
        record: AssemblyRecord,
        meta: GenomeMetadata,
    ) -> tuple[str | None, str]:
        """Select the FTP provider for a manifest record.

        Selection order (see implementation plan §Decision 2):

        1. Exact normalised match of ``meta.annotation_source`` to a provider key.
        2. Alias map match.
        3. Single-provider shortcut.
        4. Explicit Ensembl preference (if ``config.prefer_ensembl_provider``).
        5. Ambiguity — no URL generated.

        Returns:
            (provider_key, provider_status) where provider_key is None on
            ambiguity.
        """
        providers = record.providers
        if not providers:
            return None, "ambiguous"

        # Step 1 & 2: normalised / alias match
        raw_source = (meta.annotation_source or "").lower().strip()
        if raw_source:
            # Resolve alias
            resolved_source = _PROVIDER_ALIASES.get(raw_source, raw_source)
            if resolved_source in providers:
                return resolved_source, (
                    "exact_match" if resolved_source == raw_source else "alias_match"
                )

        # Step 3: single provider
        if len(providers) == 1:
            return next(iter(providers)), "single_provider"

        # Step 4: explicit Ensembl preference
        if self.config.prefer_ensembl_provider and "ensembl" in providers:
            return "ensembl", "ensembl_preference"

        # Step 5: ambiguity
        logger.info(
            "Provider ambiguity for %s: annotation_source=%r, providers=%s. "
            "No released URL will be generated.",
            meta.accession,
            meta.annotation_source,
            sorted(providers.keys()),
        )
        return None, "ambiguous"

    @staticmethod
    def _resolve_date(
        provider_dates: dict[str, ProviderDateRecord],
        meta: GenomeMetadata,
    ) -> tuple[str | None, str]:
        """Select the annotation date from a provider's date records.

        Selection order:

        1. Normalised prefix match against ``meta.annotation_date``.
        2. Latest parseable date (largest tuple).
        3. All dates unparseable → return (None, "unparseable").

        Returns:
            (date_key, date_status) where date_key is None on failure.
        """
        if not provider_dates:
            return None, "unparseable"

        # Build a list of (tuple, date_key) for parseable dates
        parsed: list[tuple[tuple[int, ...], str]] = []
        for date_key in provider_dates:
            t = _parse_manifest_date(date_key)
            if t is None:
                logger.warning(
                    "Manifest date key %r could not be parsed; skipping.", date_key
                )
                continue
            parsed.append((t, date_key))

        if not parsed:
            return None, "unparseable"

        # Step 1: try to match meta.annotation_date
        if meta.annotation_date:
            meta_tuple = _parse_manifest_date(meta.annotation_date)
            if meta_tuple is not None:
                # Prefix match: manifest (YYYY, MM) matches metadata (YYYY, MM, *)
                meta_prefix = meta_tuple[:2]
                matches = [(t, dk) for (t, dk) in parsed if t[:2] == meta_prefix]
                if len(matches) == 1:
                    return matches[0][1], "exact_match"
                if len(matches) > 1:
                    # Multiple manifest dates match — pick latest
                    best = max(matches, key=lambda x: x[0])
                    return best[1], "latest_selected"

        # Step 2: latest date
        best = max(parsed, key=lambda x: x[0])
        return best[1], "latest_selected"

    # ------------------------------------------------------------------
    # Pre-release species name normalisation (unchanged)
    # ------------------------------------------------------------------

    def _normalise_species_for_ftp(self, species_name: str) -> List[str]:
        variants: List[str] = []

        # Ensure first letter is capitalized without lowercasing the rest
        if not species_name:
            return variants
        capitalized_species = species_name[0].upper() + species_name[1:]

        # A. Canonical Ensembl-style: replace spaces, hyphens, dots with underscores
        canonical = re.sub(r"[ \-\.]", "_", capitalized_species)
        # C. Collapse multiple underscores
        canonical = re.sub(r"_+", "_", canonical)
        variants.append(canonical)

        # B. Minimal normalisation (current behavior fallback): replace spaces with underscores only
        minimal = capitalized_species.replace(" ", "_")
        if minimal not in variants:
            variants.append(minimal)

        # C. Collapse multiple underscores for minimal
        minimal_collapsed = re.sub(r"_+", "_", minimal)
        if minimal_collapsed not in variants:
            variants.append(minimal_collapsed)

        # D. Lowercase variants (for pre-release/repeat paths if needed)
        lowercase_variants = [v.lower() for v in variants]
        for lv in lowercase_variants:
            if lv not in variants:
                variants.append(lv)

        # Preserve legacy fallbacks just in case
        legacy_v2 = minimal.replace(".", "")
        if legacy_v2 not in variants:
            variants.append(legacy_v2)

        legacy_v3 = minimal.replace(".", "_")
        if legacy_v3 not in variants:
            variants.append(legacy_v3)

        legacy_v4 = re.sub(r"_+", "_", legacy_v3)
        if legacy_v4 not in variants:
            variants.append(legacy_v4)

        # Reorder if we have a cached success for this species
        if (
            hasattr(self, "_ftp_species_cache")
            and species_name in self._ftp_species_cache
        ):
            known_good = self._ftp_species_cache[species_name]
            if known_good in variants:
                variants.remove(known_good)
            variants.insert(0, known_good)

        return variants

    # ------------------------------------------------------------------
    # FTP asset resolution
    # ------------------------------------------------------------------

    def _resolve_ftp_assets(  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
        self, meta: GenomeMetadata
    ) -> Dict[str, Any]:
        """Resolve FTP file URLs for one genome.

        Resolution order:
        1. Manifest-based released path (new accession-based structure).
        2. Pre-release fallback (unchanged species-name-based structure).

        Returns a dict that the schema renderers consume.  Internal
        ``__audit_*`` keys are present; they are stripped by the orchestration
        layer before YAML serialisation.
        """
        # ftp_species_name is still needed for pre-release and repeat-library paths.
        ftp_species_name_base = meta.species_name.capitalize().replace(" ", "_")
        variants = self._normalise_species_for_ftp(meta.species_name)

        metadata_date = (
            meta.annotation_date.replace("-", "_")
            if meta.annotation_date
            else "unknown_date"
        )

        audit_decision = "excluded"
        audit_reason = "No released or pre-release FTP assets found."
        manifest_status = "n/a"
        provider_status = "n/a"
        date_status = "n/a"

        # ------------------------------------------------------------------
        # 1. Try Released logic (manifest-based)
        # ------------------------------------------------------------------
        if meta.is_released and self.manifest is not None:
            manifest_record = self.manifest.lookup(meta.accession)

            if manifest_record is None:
                manifest_status = "not_found"
                logger.debug("Accession %s not found in manifest.", meta.accession)
            else:
                manifest_status = "found"

                # Select provider
                provider, provider_status = self._resolve_provider(
                    manifest_record, meta
                )

                if provider is None:
                    # Ambiguous — fall through to pre-release
                    audit_reason = (
                        f"Provider ambiguity for {meta.accession}: "
                        f"providers={sorted(manifest_record.providers.keys())}. "
                        "Falling back to pre-release."
                    )
                else:
                    # Select date
                    provider_dates = manifest_record.providers[provider]
                    date_key, date_status = self._resolve_date(provider_dates, meta)

                    if date_key is None:
                        audit_reason = (
                            f"No usable annotation date found for "
                            f"{meta.accession}/{provider} in manifest."
                        )
                    else:
                        pdr = provider_dates[date_key]

                        # Check for at least one required annotation file
                        ann = pdr.annotation_files
                        has_gtf = "genes.gtf.gz" in ann
                        has_gff3 = "genes.gff3.gz" in ann
                        if not has_gtf and not has_gff3:
                            audit_reason = (
                                f"Manifest record {meta.accession}/{provider}/{date_key} "
                                "contains neither genes.gtf.gz nor genes.gff3.gz."
                            )
                        else:
                            # Build released asset dict
                            try:
                                acc_path = accession_to_ftp_path(meta.accession)
                            except ValueError as exc:
                                audit_reason = f"Could not compute FTP path for {meta.accession}: {exc}"
                            else:
                                resolved_annotation_files = {
                                    fname: EBI_FTP_BASE + rel_path
                                    for fname, rel_path in ann.items()
                                }
                                resolved_homology_files = {
                                    fname: EBI_FTP_BASE + rel_path
                                    for fname, rel_path in pdr.homology_files.items()
                                }
                                resolved_variation_files = {
                                    fname: EBI_FTP_BASE + rel_path
                                    for fname, rel_path in pdr.variation_files.items()
                                }
                                resolved_genome_files = {
                                    fname: EBI_FTP_BASE + rel_path
                                    for fname, rel_path in manifest_record.assembly_genome_files.items()
                                }

                                audit_decision = "included_released"
                                audit_reason = "Found released FTP assets in manifest."

                                if provider_status == "exact_match":
                                    pass  # expected
                                elif provider_status in (
                                    "single_provider",
                                    "alias_match",
                                ):
                                    logger.info(
                                        "Provider %r selected for %s via %s.",
                                        provider,
                                        meta.accession,
                                        provider_status,
                                    )
                                elif provider_status == "ensembl_preference":
                                    logger.info(
                                        "Provider 'ensembl' selected for %s by "
                                        "explicit preference policy.",
                                        meta.accession,
                                    )

                                if date_status == "latest_selected":
                                    logger.info(
                                        "No metadata date match for %s/%s; "
                                        "using latest manifest date %r.",
                                        meta.accession,
                                        provider,
                                        date_key,
                                    )

                                return {
                                    "is_released": True,
                                    "ftp_species_name": ftp_species_name_base,
                                    "acc_ftp_path": acc_path,
                                    "resolved_provider": provider,
                                    "resolved_date": date_key,
                                    "annotation_files": resolved_annotation_files,
                                    "genome_files": resolved_genome_files,
                                    "homology_files": resolved_homology_files,
                                    "variation_files": resolved_variation_files,
                                    "audit_decision": audit_decision,
                                    "audit_reason": audit_reason,
                                    "__audit_manifest_status__": manifest_status,
                                    "__audit_provider_status__": provider_status,
                                    "__audit_date_status__": date_status,
                                }
        elif meta.is_released and self.manifest is None:
            manifest_status = "manifest_unavailable"
            audit_reason = (
                "Manifest unavailable; cannot resolve released FTP path. "
                "Attempting pre-release fallback."
            )

        # ------------------------------------------------------------------
        # 2. Try Pre-release Fallback (unchanged species-name-based logic)
        # ------------------------------------------------------------------
        resolved_pre_variant = None
        pre_release_urls: dict[str, str] = {}
        if self.ftp_client:
            for variant in variants:
                fb_gtf = self.ftp_client.check_pre_release_file(
                    variant, meta.accession, ".gtf.gz"
                )
                if fb_gtf:
                    resolved_pre_variant = variant
                    pre_release_urls["annotation_gtf"] = fb_gtf

                    fb_gff = self.ftp_client.check_pre_release_file(
                        variant, meta.accession, ".gff3.gz"
                    )
                    if not fb_gff:
                        fb_gff = self.ftp_client.check_pre_release_file(
                            variant, meta.accession, ".gff3"
                        )
                    if fb_gff:
                        pre_release_urls["annotation_gff3"] = fb_gff

                    fb_pep = self.ftp_client.check_pre_release_file(
                        variant, meta.accession, ".pep.fa.gz"
                    )
                    if fb_pep:
                        pre_release_urls["proteins"] = fb_pep

                    fb_cdna = self.ftp_client.check_pre_release_file(
                        variant, meta.accession, ".cdna.fa.gz"
                    )
                    if fb_cdna:
                        pre_release_urls["transcripts"] = fb_cdna

                    fb_soft = self.ftp_client.check_pre_release_file(
                        variant, meta.accession, ".dna.softmasked.fa.gz"
                    )
                    if fb_soft:
                        pre_release_urls["softmasked_genome"] = fb_soft
                    break

        if resolved_pre_variant:
            self._ftp_species_cache[meta.species_name] = resolved_pre_variant
            if resolved_pre_variant != ftp_species_name_base:
                logger.info(
                    "Resolved FTP species name: input=%r used=%r",
                    meta.species_name,
                    resolved_pre_variant,
                )

            return {
                "is_released": False,
                "ftp_species_name": resolved_pre_variant,
                "acc_ftp_path": None,
                "resolved_provider": None,
                "resolved_date": metadata_date,
                "annotation_files": {},
                "genome_files": {},
                "homology_files": {},
                "variation_files": {},
                "pre_release_urls": pre_release_urls,
                "audit_decision": "included_prerelease",
                "audit_reason": "Found pre-release FTP assets.",
                "__audit_manifest_status__": manifest_status,
                "__audit_provider_status__": provider_status,
                "__audit_date_status__": date_status,
            }

        return {
            "is_released": False,
            "ftp_species_name": ftp_species_name_base,
            "acc_ftp_path": None,
            "resolved_provider": None,
            "resolved_date": metadata_date,
            "annotation_files": {},
            "genome_files": {},
            "homology_files": {},
            "variation_files": {},
            "audit_decision": "excluded",
            "audit_reason": audit_reason,
            "__audit_manifest_status__": manifest_status,
            "__audit_provider_status__": provider_status,
            "__audit_date_status__": date_status,
        }

    # ------------------------------------------------------------------
    # Schema renderers
    # ------------------------------------------------------------------

    def _render_standard(  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
        self, meta: GenomeMetadata
    ) -> Dict[str, Any]:
        """Renders Schema A: Standard Projects (VGP, DToL, ERGA)"""
        doc: Dict[str, Any] = {}

        # species display name must only come from the scientific name -- never from strain,
        # sample description, habitat text, or any other free-text metadata field.
        # (strain is displayed separately in _render_mouse; it must never appear here.)
        doc["species"] = meta.species_name

        # Icon resolution via taxonomy lineage (see icon_resolver.py).
        # The resolver tries: metadata DB lineage → NCBI Entrez → BUSCO hint.
        # Results are cached per taxon_id for the run.
        icon, icon_matched_term, icon_source = self.icon_resolver.resolve_icon(meta)
        doc["image"] = icon
        doc["__audit_image_source__"] = f"{icon_source}:{icon_matched_term}"

        if self.config.scrape_ncbi_submitter and meta.assembly_submitter:
            doc["submitted_by"] = meta.assembly_submitter

        doc["accession"] = meta.accession
        doc["annotation_method"] = meta.annotation_method or "BRAKER2"

        ftp_resolution = self._resolve_ftp_assets(meta)
        target_released = ftp_resolution["is_released"]
        ftp_species_name = ftp_resolution["ftp_species_name"]

        doc["__audit_decision__"] = ftp_resolution["audit_decision"]
        doc["__audit_reason__"] = ftp_resolution["audit_reason"]
        doc["__audit_resolved_date__"] = ftp_resolution["resolved_date"]
        doc["__audit_manifest_status__"] = ftp_resolution.get(
            "__audit_manifest_status__", ""
        )
        doc["__audit_provider_status__"] = ftp_resolution.get(
            "__audit_provider_status__", ""
        )
        doc["__audit_date_status__"] = ftp_resolution.get("__audit_date_status__", "")

        if ftp_resolution["audit_decision"] == "excluded":
            return doc  # Returns only audit keys

        if target_released:
            meta.annotation_date = ftp_resolution["resolved_date"].replace("_", "-")
            ann = ftp_resolution.get("annotation_files", {})
            genome = ftp_resolution.get("genome_files", {})

            if "genes.gtf.gz" in ann:
                doc["annotation_gtf"] = ann["genes.gtf.gz"]
            if "genes.gff3.gz" in ann:
                doc["annotation_gff3"] = ann["genes.gff3.gz"]
            if "pep.fa.bgz" in ann:
                doc["proteins"] = ann["pep.fa.bgz"]
            if "cdna.fa.bgz" in ann:
                doc["transcripts"] = ann["cdna.fa.bgz"]
            if "softmasked.fa.bgz" in genome:
                doc["softmasked_genome"] = genome["softmasked.fa.bgz"]

            # repeat_library — uses different FTP hierarchy; ftp_species_name still correct
            repeat_species = ftp_species_name.lower()
            repeat_url = (
                "https://ftp.ebi.ac.uk/pub/databases/ensembl/"
                f"repeats/unfiltered_repeatmodeler/species/"
                f"{repeat_species}/{meta.accession}.repeatmodeler.fa"
            )
            if check_url_status(repeat_url):
                doc["repeat_library"] = repeat_url
        else:
            pre_urls = ftp_resolution.get("pre_release_urls", {})
            for k, v in pre_urls.items():
                doc[k] = v

        acc_path = ftp_resolution.get("acc_ftp_path")
        if target_released and acc_path:
            doc["ftp_dumps"] = EBI_FTP_BASE + acc_path + "/"
        else:
            doc["ftp_dumps"] = (
                "https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/"
                f"{ftp_species_name}/{meta.accession}/"
            )

        # Linkages
        if meta.is_on_rapid:
            doc["ensembl_link"] = f"https://rapid.ensembl.org/{meta.species_name}"
            doc["__audit_beta_status__"] = "skipped_rapid"
        elif self.config.allow_beta_urls:
            beta_link, beta_status = self._resolve_beta_link(meta, target_released)
            doc["beta_link"] = beta_link
            doc["__audit_beta_status__"] = beta_status

        # Quality
        if meta.busco_score:
            doc["busco_score"] = meta.busco_score
        if meta.busco_lineage:
            doc["busco_lineage"] = meta.busco_lineage

        # Relationships
        if meta.alternate_of:
            doc["alternate"] = meta.alternate_of

        return {k: v for k, v in doc.items() if v is not None}

    def _render_hprc(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema B: HPRC"""
        doc: Dict[str, Any] = {}

        doc["assembly"] = meta.assembly_name
        if meta.parent_of_origin:
            doc["parent_of_origin"] = meta.parent_of_origin

        doc["assembly_accession"] = meta.accession
        doc["assembly_link"] = (
            f"https://www.ebi.ac.uk/ena/browser/view/{meta.accession}"
        )

        if meta.assembly_submitter:
            doc["assembly_submitter"] = meta.assembly_submitter

        ftp_resolution = self._resolve_ftp_assets(meta)
        target_released = ftp_resolution["is_released"]
        ftp_species_name = ftp_resolution["ftp_species_name"]

        doc["__audit_decision__"] = ftp_resolution["audit_decision"]
        doc["__audit_reason__"] = ftp_resolution["audit_reason"]
        doc["__audit_resolved_date__"] = ftp_resolution["resolved_date"]
        doc["__audit_manifest_status__"] = ftp_resolution.get(
            "__audit_manifest_status__", ""
        )
        doc["__audit_provider_status__"] = ftp_resolution.get(
            "__audit_provider_status__", ""
        )
        doc["__audit_date_status__"] = ftp_resolution.get("__audit_date_status__", "")

        if ftp_resolution["audit_decision"] == "excluded":
            return doc  # Returns only audit keys

        if target_released:
            meta.annotation_date = ftp_resolution["resolved_date"].replace("_", "-")
            ann = ftp_resolution.get("annotation_files", {})

            if "genes.gtf.gz" in ann:
                doc["annotation_gtf"] = ann["genes.gtf.gz"]
            if "genes.gff3.gz" in ann:
                doc["annotation_gff3"] = ann["genes.gff3.gz"]
            if "pep.fa.bgz" in ann:
                doc["proteins"] = ann["pep.fa.bgz"]
            if "cdna.fa.bgz" in ann:
                doc["transcripts"] = ann["cdna.fa.bgz"]
        else:
            pre_urls = ftp_resolution.get("pre_release_urls", {})
            for k, v in pre_urls.items():
                doc[k] = v

        # VEP resolution: try the new manifest first (no VEP section yet),
        # then fall back to the legacy species.json manifest for HPRC.
        vep_url, vep_status = self._resolve_vep_url(meta, ftp_resolution)
        if vep_url:
            doc["variants_vep"] = vep_url
        doc["__audit_vep_status__"] = vep_status

        acc_path = ftp_resolution.get("acc_ftp_path")
        if target_released and acc_path:
            doc["ftp_dumps"] = EBI_FTP_BASE + acc_path + "/"
        else:
            doc["ftp_dumps"] = (
                "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
                f"{ftp_species_name}/{meta.accession}/"
            )

        beta_link, beta_status = self._resolve_beta_link(meta, target_released)
        doc["beta_link"] = beta_link
        doc["__audit_beta_status__"] = beta_status

        return {k: v for k, v in doc.items() if v is not None}

    def _render_mouse(self, meta: GenomeMetadata) -> Dict[str, Any]:
        """Renders Schema C: Mouse Genomes"""
        doc: Dict[str, Any] = {}

        doc["species"] = meta.species_name
        if meta.strain:
            doc["strain"] = meta.strain

        doc["accession"] = meta.accession

        ftp_resolution = self._resolve_ftp_assets(meta)
        target_released = ftp_resolution["is_released"]
        ftp_species_name = ftp_resolution["ftp_species_name"]

        doc["__audit_decision__"] = ftp_resolution["audit_decision"]
        doc["__audit_reason__"] = ftp_resolution["audit_reason"]
        doc["__audit_resolved_date__"] = ftp_resolution["resolved_date"]
        doc["__audit_manifest_status__"] = ftp_resolution.get(
            "__audit_manifest_status__", ""
        )
        doc["__audit_provider_status__"] = ftp_resolution.get(
            "__audit_provider_status__", ""
        )
        doc["__audit_date_status__"] = ftp_resolution.get("__audit_date_status__", "")

        if ftp_resolution["audit_decision"] == "excluded":
            return doc  # Returns only audit keys

        if target_released:
            meta.annotation_date = ftp_resolution["resolved_date"].replace("_", "-")
            ann = ftp_resolution.get("annotation_files", {})
            genome = ftp_resolution.get("genome_files", {})

            if "genes.gtf.gz" in ann:
                doc["annotation_gtf"] = ann["genes.gtf.gz"]
            if "genes.gff3.gz" in ann:
                doc["annotation_gff3"] = ann["genes.gff3.gz"]
            if "pep.fa.bgz" in ann:
                doc["proteins"] = ann["pep.fa.bgz"]
            if "cdna.fa.bgz" in ann:
                doc["transcripts"] = ann["cdna.fa.bgz"]
            if "softmasked.fa.bgz" in genome:
                doc["softmasked_genome"] = genome["softmasked.fa.bgz"]
        else:
            pre_urls = ftp_resolution.get("pre_release_urls", {})
            for k, v in pre_urls.items():
                doc[k] = v

        acc_path = ftp_resolution.get("acc_ftp_path")
        if target_released and acc_path:
            doc["ftp_dumps"] = EBI_FTP_BASE + acc_path + "/"
        else:
            doc["ftp_dumps"] = (
                "https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/"
                f"{ftp_species_name}/{meta.accession}/"
            )

        if self.config.allow_beta_urls:
            beta_link, beta_status = self._resolve_beta_link(meta, target_released)
            doc["beta_link"] = beta_link
            doc["__audit_beta_status__"] = beta_status

        if meta.alternate_of:
            doc["alternate"] = meta.alternate_of

        return {k: v for k, v in doc.items() if v is not None}

    def _build_rapid_ftp_url(self, meta: GenomeMetadata, resource_type: str) -> str:
        """Helper to build Rapid Release FTP URLs"""
        base = f"https://ftp.ensembl.org/pub/rapid-release/species/{meta.species_name}"
        return f"{base}/{resource_type}"

"""
Haplotype resolver for Ensembl project YAML generation.

Identifies alternate haplotype pairs among assemblies in the same project
dataset using NCBI Datasets metadata (BioSample accession, sample name,
assembly naming conventions).

Matching strategy (in priority order):
1. **BioSample accession** — assemblies sharing the same BioSample are from
   the same individual and are haplotype pairs.
2. **Sample/isolate name** — if BioSample is unavailable, assemblies with
   the same sample_name or isolate in the NCBI metadata are paired.
3. **Assembly name heuristics** — as a weak fallback, assembly names
   containing hap1/hap2, pat/mat, primary/alternate patterns are matched
   within the same species+taxon_id.

Only assemblies that are *both present* in the generated project dataset
are linked.  No external URLs are invented.
"""

from __future__ import annotations

import json
import logging
import re
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import requests

if TYPE_CHECKING:
    from ensembl.genes.projects.models import GenomeMetadata

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# NCBI Datasets batch lookup
# ---------------------------------------------------------------------------

_DATASETS_API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"


def _fetch_assembly_metadata_batch(
    accessions: List[str],
) -> Dict[str, Dict[str, str]]:
    """Fetch BioSample and sample metadata for a batch of accessions.

    Uses the NCBI Datasets REST API (POST endpoint for bulk lookup).

    Returns a dict keyed by accession with values::

        {
            "biosample": "SAMN12345678",
            "sample_name": "isolate_name",
            "assembly_name": "asm_name",
            "taxon_id": 12345,
        }

    Entries that could not be fetched are silently omitted.
    """
    if not accessions:
        return {}

    result: Dict[str, Dict] = {}

    # NCBI Datasets v2 API: POST with {"accessions": [...]}
    url = f"{_DATASETS_API_BASE}"
    # The v2 API accepts GET with comma-separated accessions for small batches.
    # For robustness we do chunked GET requests (100 at a time).
    chunk_size = 100
    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i : i + chunk_size]
        accs_param = ",".join(chunk)
        api_url = f"{url}/{accs_param}/dataset_report" "?filters.assembly_version=all"
        try:
            resp = requests.get(api_url, timeout=30)
            if resp.status_code != 200:
                logger.debug(
                    "NCBI Datasets API returned HTTP %d for batch starting %s",
                    resp.status_code,
                    chunk[0],
                )
                continue

            data = resp.json()
            for report in data.get("reports", []):
                acc = report.get("accession", "")
                if not acc:
                    continue

                assembly_info = report.get("assembly_info", {}) or {}
                biosample_info = assembly_info.get("biosample", {}) or {}

                # BioSample accession (e.g. "SAMN12345678")
                biosample_acc = (
                    biosample_info.get("accession", "")
                    if isinstance(biosample_info, dict)
                    else ""
                )

                # Sample name / isolate — NCBI stores these in biosample
                # attributes or at the top level
                sample_name = ""
                if isinstance(biosample_info, dict):
                    for attr in biosample_info.get("attributes", []):
                        attr_name = (attr.get("name") or "").lower()
                        if attr_name in ("sample_name", "isolate", "sample name"):
                            sample_name = attr.get("value", "")
                            break

                asm_name = assembly_info.get("assembly_name", "")
                taxon_id = (report.get("organism", {}) or {}).get("tax_id")

                result[acc] = {
                    "biosample": biosample_acc,
                    "sample_name": sample_name,
                    "assembly_name": asm_name,
                    "taxon_id": taxon_id,
                }

        except requests.RequestException as exc:
            logger.warning("NCBI Datasets API error for batch: %s", exc)
        except (json.JSONDecodeError, KeyError) as exc:
            logger.warning("NCBI Datasets API parse error: %s", exc)

    logger.info(
        "Fetched NCBI assembly metadata for %d / %d accessions.",
        len(result),
        len(accessions),
    )
    return result


# ---------------------------------------------------------------------------
# Assembly name haplotype patterns
# ---------------------------------------------------------------------------

# Pairs of (regex, group_key_extractor) for matching haplotype naming.
# group_key is a normalised identifier that should be identical for both
# haplotypes of the same individual.

_HAP_PATTERNS = [
    # hap1 / hap2
    (re.compile(r"^(.+?)[\._\-]?(hap[12]|haplotype[12])(.*)$", re.I), "hap"),
    # pat / mat / paternal / maternal
    (re.compile(r"^(.+?)[\._\-]?(pat(?:ernal)?|mat(?:ernal)?)(.*)$", re.I), "parent"),
    # primary / alternate
    (re.compile(r"^(.+?)[\._\-]?(primary|alternate|alt)(.*)$", re.I), "alt"),
]


def _extract_hap_group_key(assembly_name: str) -> Optional[Tuple[str, str]]:
    """Extract a grouping key and haplotype label from an assembly name.

    Returns (group_key, hap_label) or None if no pattern matches.
    The group_key should be identical for both haplotypes from the same
    individual; hap_label distinguishes them.
    """
    if not assembly_name:
        return None
    for pattern, _ in _HAP_PATTERNS:
        m = pattern.match(assembly_name)
        if m:
            prefix = m.group(1).rstrip("_.- ")
            hap_label = m.group(2).lower()
            suffix = m.group(3).lstrip("_.- ")
            group_key = f"{prefix}||{suffix}".lower()
            return group_key, hap_label
    return None


# ---------------------------------------------------------------------------
# Main resolver
# ---------------------------------------------------------------------------


class HaplotypeResolver:
    """Identifies alternate haplotype pairs within a project dataset.

    Only pairs assemblies that are both present in the provided genome
    list.  Never invents external links.
    """

    def find_alternate_haplotypes(
        self,
        genomes: List["GenomeMetadata"],
    ) -> Dict[str, str]:
        """Find alternate haplotype pairs among the given genomes.

        Parameters
        ----------
        genomes:
            All GenomeMetadata records that will appear in the final YAML
            (i.e. after deduplication/exclusion filtering).

        Returns
        -------
        Dict mapping accession → alternate accession.
        The mapping is bidirectional: if A→B then B→A.
        Only accessions present in *genomes* appear as values.
        """
        if len(genomes) < 2:
            return {}

        accessions = [g.accession for g in genomes]
        acc_set = set(accessions)

        # Build indexes
        by_accession: Dict[str, "GenomeMetadata"] = {g.accession: g for g in genomes}

        pairs: Dict[str, str] = {}

        # ------ Source 1: BioSample from NCBI Datasets API ------
        ncbi_meta = _fetch_assembly_metadata_batch(accessions)

        # Group by BioSample accession
        biosample_groups: Dict[str, List[str]] = {}
        for acc, info in ncbi_meta.items():
            bs = info.get("biosample", "")
            if bs and acc in acc_set:
                biosample_groups.setdefault(bs, []).append(acc)

        for bs, group_accs in biosample_groups.items():
            if len(group_accs) == 2:
                a, b = group_accs
                if a not in pairs and b not in pairs:
                    pairs[a] = b
                    pairs[b] = a
                    logger.info("Haplotype pair (biosample=%s): %s <-> %s", bs, a, b)
            elif len(group_accs) > 2:
                logger.debug(
                    "BioSample %s has %d assemblies; skipping auto-pair "
                    "(ambiguous).",
                    bs,
                    len(group_accs),
                )

        # ------ Source 2: Sample name / isolate ------
        sample_groups: Dict[str, List[str]] = {}
        for acc, info in ncbi_meta.items():
            sn = info.get("sample_name", "").strip()
            if sn and acc in acc_set and acc not in pairs:
                # Also require same taxon_id to avoid cross-species matches
                tid = info.get("taxon_id")
                key = f"{tid}||{sn}" if tid else sn
                sample_groups.setdefault(key, []).append(acc)

        for key, group_accs in sample_groups.items():
            unpaired = [a for a in group_accs if a not in pairs]
            if len(unpaired) == 2:
                a, b = unpaired
                pairs[a] = b
                pairs[b] = a
                logger.info("Haplotype pair (sample_name=%s): %s <-> %s", key, a, b)

        # ------ Source 3: Assembly name heuristics (weak fallback) ------
        # Group by species+taxon_id for safety, then by naming pattern
        species_groups: Dict[Tuple, List["GenomeMetadata"]] = {}
        for g in genomes:
            if g.accession in pairs:
                continue
            sp_key = (g.species_name, g.taxon_id)
            species_groups.setdefault(sp_key, []).append(g)

        for sp_key, sp_genomes in species_groups.items():
            if len(sp_genomes) < 2:
                continue
            hap_groups: Dict[str, List["GenomeMetadata"]] = {}
            for g in sp_genomes:
                # Use NCBI assembly name if available, fall back to meta
                asm_name = ""
                if g.accession in ncbi_meta:
                    asm_name = ncbi_meta[g.accession].get("assembly_name", "")
                if not asm_name:
                    asm_name = g.assembly_name or ""

                result = _extract_hap_group_key(asm_name)
                if result:
                    group_key, _label = result
                    hap_groups.setdefault(group_key, []).append(g)

            for gk, hap_genomes in hap_groups.items():
                unpaired = [g for g in hap_genomes if g.accession not in pairs]
                if len(unpaired) == 2:
                    a, b = unpaired[0].accession, unpaired[1].accession
                    pairs[a] = b
                    pairs[b] = a
                    logger.info(
                        "Haplotype pair (assembly_name pattern=%s): %s <-> %s",
                        gk,
                        a,
                        b,
                    )

        logger.info(
            "Haplotype resolution complete: %d pairs found for %d genomes.",
            len(pairs) // 2,
            len(genomes),
        )
        return pairs

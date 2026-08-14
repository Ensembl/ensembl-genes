"""Resolve NCBI RefSeq annotation-report URLs by assembly accession.

This module is intentionally isolated from the core loader. The NCBI lookup
can be removed later without changing GFF parsing or database loading logic.
"""

# The resolver supports multiple NCBI report schema variants and intentionally
# uses a small recursive lookup rather than coupling the loader to one schema.
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-boolean-expressions

from __future__ import annotations

import argparse
import logging
import re
from collections.abc import Mapping
from typing import Any, Protocol
from urllib.parse import quote

LOGGER = logging.getLogger(__name__)
NCBI_DATASETS_REPORT_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
    "{accession}/dataset_report"
)


class HttpResponse(Protocol):
    """Minimal response interface required by the resolver."""

    def raise_for_status(self) -> None:
        """Raise when the NCBI request failed."""

    def json(self) -> Mapping[str, Any]:
        """Return the decoded response payload."""


class HttpSession(Protocol):
    """Minimal session interface used for test injection."""

    def get(self, url: str, timeout: int) -> HttpResponse:
        """Perform an HTTP GET request."""


def _normalized_key(key: str) -> str:
    """Normalize API field names for matching across report versions."""

    return re.sub(r"[^a-z0-9]", "", key.lower())


def _find_annotation_report_url(
    value: Any,
    in_annotation_context: bool = False,
) -> str | None:
    """Find an annotation-report URL in a nested NCBI report payload."""

    if isinstance(value, Mapping):
        for key, child in value.items():
            normalized = _normalized_key(str(key))
            if (
                (
                    ("annotation" in normalized and "report" in normalized)
                    or (in_annotation_context and "report" in normalized)
                )
                and "url" in normalized
                and isinstance(child, str)
                and child
            ):
                return child
            found = _find_annotation_report_url(
                child,
                in_annotation_context or "annotation" in normalized,
            )
            if found:
                return found
    elif isinstance(value, list):
        for child in value:
            found = _find_annotation_report_url(child, in_annotation_context)
            if found:
                return found
    return None


def get_annotation_report_url(
    assembly_accession: str,
    *,
    session: HttpSession | None = None,
    timeout: int = 20,
) -> str:
    """Return NCBI's official annotation-report URL for a GCF accession."""

    accession = assembly_accession.strip()
    if not accession.upper().startswith("GCF_"):
        raise ValueError(
            "NCBI eukaryotic annotation reports require a GCF accession: "
            f"{assembly_accession}"
        )

    if session is None:
        try:
            import requests  # pylint: disable=import-outside-toplevel
        except ImportError as error:  # pragma: no cover - dependency issue
            raise ImportError(
                "The NCBI annotation resolver requires the 'requests' package."
            ) from error
        session = requests.Session()

    url = NCBI_DATASETS_REPORT_URL.format(accession=quote(accession, safe=""))
    response = session.get(url, timeout=timeout)
    response.raise_for_status()
    report_url = _find_annotation_report_url(response.json())
    if not report_url:
        raise ValueError(
            "NCBI did not return an annotation-report URL for " f"{assembly_accession}"
        )
    return report_url


def main() -> None:
    """Print the annotation-report URL for one GCF accession."""

    parser = argparse.ArgumentParser()
    parser.add_argument("assembly_accession")
    args = parser.parse_args()
    print(get_annotation_report_url(args.assembly_accession))


if __name__ == "__main__":  # pragma: no cover - CLI wrapper
    main()

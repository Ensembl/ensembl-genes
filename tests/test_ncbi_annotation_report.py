"""Tests for NCBI annotation provenance URL resolution."""

# pylint: disable=missing-class-docstring,missing-function-docstring
# pylint: disable=too-few-public-methods,wrong-import-position

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

SRC_PATH = Path(__file__).resolve().parents[1] / "src" / "python"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from ensembl.genes.ensembl_loading.ncbi_annotation_report import (
    NCBI_DATASETS_GENOME_URL,
    get_annotation_report_url,
)


class FakeResponse:
    def __init__(self, payload: dict[str, Any]) -> None:
        self.payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict[str, Any]:
        return self.payload


class FakeSession:
    def __init__(self, payload: dict[str, Any]) -> None:
        self.payload = payload
        self.url = ""

    def get(self, url: str, timeout: int) -> FakeResponse:
        del timeout
        self.url = url
        return FakeResponse(self.payload)


def test_annotation_report_url_uses_report_url_when_present() -> None:
    session = FakeSession(
        {"reports": [{"annotation_info": {"report_url": "https://example/report"}}]}
    )

    assert (
        get_annotation_report_url("GCF_000001.1", session=session)
        == "https://example/report"
    )


def test_annotation_report_url_falls_back_to_datasets_genome_page() -> None:
    session = FakeSession(
        {
            "reports": [
                {
                    "assembly_info": {"submitter": "WELLCOME SANGER INSTITUTE"},
                    "annotation_info": {
                        "name": "GCF_964188405.1-RS_2025_09",
                        "provider": "Vertebrate Genome Laboratory",
                    },
                }
            ]
        }
    )

    assert get_annotation_report_url("GCF_964188405.1", session=session) == (
        NCBI_DATASETS_GENOME_URL.format(accession="GCF_964188405.1")
    )

"""
Output layout of an experiment workspace (``output_dir`` in the configuration).

    <output_dir>/
      cache/checksums.json                 sha256 cache keyed by path, size, mtime
      prepared/<reference_id>/             derived reference inputs + preparation.json
      prepared/runs/<run_id>/              derived query inputs (only when preparation is configured)
      validation/validation.json           latest validation report
      comparisons/<run_id>/<policy>/       comparator outputs + run_record.json (status complete)
      comparisons/<run_id>/<policy>.attempts.jsonl   every attempt, including failures
      comparisons/_partial/                in-progress / failed attempt directories
      imported/import_report.json          what `benchmark import` matched, pilots, rejects
      independent/<run_id>/<policy>.*      independent analysis tables (computed or imported)
      dashboard/dashboard.sqlite           data read by the dashboard
      logs/

Source inputs are never written to.
"""

from __future__ import annotations

import json
from pathlib import Path

from ensembl.genes.annotation_qc.benchmark.config import Experiment
from ensembl.genes.annotation_qc.benchmark.provenance import ChecksumCache


class Workspace:
    def __init__(self, experiment: Experiment):
        self.experiment = experiment
        self.root = experiment.output_dir
        self.checksums = ChecksumCache(self.root / "cache" / "checksums.json")

    # directories -------------------------------------------------------------
    def prepared_dir(self, reference_id: str) -> Path:
        return self.root / "prepared" / reference_id

    def prepared_run_dir(self, run_id: str) -> Path:
        return self.root / "prepared" / "runs" / run_id

    def comparison_dir(self, run_id: str, policy_id: str) -> Path:
        return self.root / "comparisons" / run_id / policy_id

    def attempts_file(self, run_id: str, policy_id: str) -> Path:
        return self.root / "comparisons" / run_id / f"{policy_id}.attempts.jsonl"

    def partial_root(self) -> Path:
        return self.root / "comparisons" / "_partial"

    def independent_dir(self, run_id: str) -> Path:
        return self.root / "independent" / run_id

    def dashboard_db(self) -> Path:
        return self.root / "dashboard" / "dashboard.sqlite"

    def validation_file(self) -> Path:
        return self.root / "validation" / "validation.json"

    def import_report(self) -> Path:
        return self.root / "imported" / "import_report.json"

    # records -----------------------------------------------------------------
    def run_record(self, run_id: str, policy_id: str) -> dict | None:
        path = self.comparison_dir(run_id, policy_id) / "run_record.json"
        return read_json(path)

    def write_run_record(self, run_id: str, policy_id: str, record: dict) -> None:
        write_json(self.comparison_dir(run_id, policy_id) / "run_record.json", record)

    def append_attempt(self, run_id: str, policy_id: str, record: dict) -> None:
        path = self.attempts_file(run_id, policy_id)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a") as handle:
            handle.write(json.dumps(record, sort_keys=True, default=str) + "\n")

    def attempts(self, run_id: str, policy_id: str) -> list[dict]:
        path = self.attempts_file(run_id, policy_id)
        if not path.exists():
            return []
        return [
            json.loads(line) for line in path.read_text().splitlines() if line.strip()
        ]

    def preparation_record(self, reference_id: str) -> dict | None:
        return read_json(self.prepared_dir(reference_id) / "preparation.json")

    def run_preparation_record(self, run_id: str) -> dict | None:
        return read_json(self.prepared_run_dir(run_id) / "preparation.json")


def read_json(path: Path) -> dict | None:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None


def write_json(path: Path, data) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=1, sort_keys=True, default=str))
    tmp.replace(path)

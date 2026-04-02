#!/usr/bin/env python3
"""Benchmark harness for mapped-annotation truth sets."""

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class BenchmarkResult:
    case_name: str
    passed: bool
    metrics: Dict[str, float]
    failures: Dict[str, Tuple[float, float]]


def _parse_rate(value) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str) and value.endswith("%"):
        return float(value.rstrip("%")) / 100.0
    if isinstance(value, str):
        return float(value)
    return 0.0


def extract_core_metrics(report: Dict) -> Dict[str, float]:
    """Extract normalized metrics from pipeline JSON report."""
    summary = report.get("summary", {})
    validation = report.get("validation", {})
    metrics = {
        "gene_mapping_rate": _parse_rate(summary.get("genes", {}).get("mapping_rate", 0.0)),
        "transcript_mapping_rate": _parse_rate(summary.get("transcripts", {}).get("mapping_rate", 0.0)),
        "exon_mapping_rate": _parse_rate(summary.get("exons", {}).get("mapping_rate", 0.0)),
        "splice_site_validation_rate": _parse_rate(validation.get("splice_sites", {}).get("rate", 0.0)),
        "start_codon_validation_rate": _parse_rate(validation.get("start_codons", {}).get("rate", 0.0)),
        "stop_codon_validation_rate": _parse_rate(validation.get("stop_codons", {}).get("rate", 0.0)),
    }
    refinement = report.get("details", {}) or report.get("refinement", {})
    protein_qc = refinement.get("protein_qc", {})
    genes_checked = max(1, int(protein_qc.get("genes_checked", 0)))
    metrics["protein_qc_ok_rate"] = float(protein_qc.get("genes_ok", 0)) / genes_checked
    return metrics


def evaluate_case(case_name: str, report: Dict, thresholds: Dict[str, float]) -> BenchmarkResult:
    """Evaluate one benchmark case against threshold dictionary."""
    metrics = extract_core_metrics(report)
    failures: Dict[str, Tuple[float, float]] = {}
    for metric_name, min_expected in thresholds.items():
        observed = metrics.get(metric_name, 0.0)
        if observed < float(min_expected):
            failures[metric_name] = (observed, float(min_expected))
    return BenchmarkResult(
        case_name=case_name,
        passed=len(failures) == 0,
        metrics=metrics,
        failures=failures,
    )


def run_benchmarks(manifest: Dict, base_dir: Path) -> List[BenchmarkResult]:
    """Run benchmark evaluations from manifest."""
    results: List[BenchmarkResult] = []
    for case in manifest.get("cases", []):
        case_name = case["name"]
        report_path = Path(case["report_json"])
        if not report_path.is_absolute():
            report_path = base_dir / report_path
        with open(report_path, "r", encoding="utf-8") as handle:
            report = json.load(handle)
        thresholds = case.get("thresholds", {})
        results.append(evaluate_case(case_name, report, thresholds))
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate mapping benchmark reports against acceptance thresholds."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="JSON manifest with benchmark cases and thresholds.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional path to write benchmark summary JSON.",
    )
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    with open(manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    results = run_benchmarks(manifest, manifest_path.parent)
    overall_pass = all(r.passed for r in results)

    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(f"[{status}] {result.case_name}")
        if result.failures:
            for metric, (observed, expected) in sorted(result.failures.items()):
                print(f"  - {metric}: observed={observed:.4f}, expected>={expected:.4f}")

    summary = {
        "overall_pass": overall_pass,
        "results": [asdict(r) for r in results],
    }
    if args.output:
        out_path = Path(args.output)
        with open(out_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)

    raise SystemExit(0 if overall_pass else 1)


if __name__ == "__main__":
    main()

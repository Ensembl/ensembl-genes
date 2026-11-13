# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Orchestrate rescoring of stored confidence evidence without rerunning DIAMOND."""
from __future__ import annotations
import argparse
from pathlib import Path
from ensembl.genes.annotation_qc.config.high_confidence import CONFIDENCE_POLICY_VERSION
from ensembl.genes.annotation_qc.parsers.high_confidence import (
    parse_confidence_evidence,
    parse_json_manifest,
    parse_generated_gff,
    file_sha256,
)
from ensembl.genes.annotation_qc.metrics.high_confidence import rescore_evidence
from ensembl.genes.annotation_qc.reports.high_confidence import write_rescored_outputs
from ensembl.genes.annotation_qc.runners.high_confidence import implementation_hashes

COMPANIONS = [
    ".introns.tsv",
    ".cds.fa",
    ".prots.fa",
    ".prots.res",
    ".diamond.tsv",
    ".normalization.tsv",
    ".input_validation.tsv",
    ".translation_exceptions.tsv",
]


def add_arguments(parser):
    """Declare rescoring inputs and preserve originals by requiring a new destination."""
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--evidence", type=Path)
    source.add_argument(
        "--input-root", type=Path, help="Contains SPECIES/METHOD/*.evidence.tsv"
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--no-utr", "--no_utr", dest="no_utr", action="store_true")


def register(subparsers):
    """Register stored-evidence rescoring on annotation-qc."""
    parser = subparsers.add_parser(
        "high-confidence-rescore",
        help="Rescore stored evidence with confidence policy v2.",
    )
    add_arguments(parser)
    parser.set_defaults(func=run_cli)


def rescore(source, outdir, no_utr=False):
    """Load parsed tables, call the metric and delegate all output formatting."""
    source = source.resolve()
    if not source.name.endswith(".evidence.tsv"):
        raise ValueError("Input must end in .evidence.tsv")
    prefix = source.name[: -len(".evidence.tsv")]
    frame = parse_confidence_evidence(source)
    manifest_path = source.with_name(prefix + ".manifest.json")
    original = parse_json_manifest(manifest_path) if manifest_path.exists() else {}
    if original.get("processed_transcript_count", len(frame)) != len(frame):
        raise ValueError("Manifest transcript count does not match evidence")
    gff_path = source.with_name(prefix + ".annotations.gff3")
    gff = parse_generated_gff(gff_path) if gff_path.exists() else None
    result = rescore_evidence(frame, no_utr, gff)
    outdir = outdir.resolve()
    paths = {
        key: outdir / (prefix + suffix)
        for key, suffix in [
            ("evidence", ".evidence.tsv"),
            ("manifest", ".manifest.json"),
            ("changes", ".confidence_changes.tsv"),
        ]
    }
    if gff is not None:
        paths["gff"] = outdir / (prefix + ".annotations.gff3")
    companions = {
        outdir / (prefix + suffix): source.with_name(prefix + suffix)
        for suffix in COMPANIONS
        if source.with_name(prefix + suffix).is_file()
    }
    if any(p.exists() or p.is_symlink() for p in [*paths.values(), *companions]):
        raise ValueError("Output files already exist; choose a fresh directory")
    code_hash, hashes = implementation_hashes()
    manifest = dict(
        original,
        schema_version=1,
        script_sha256=code_hash,
        implementation_hashes=hashes,
        confidence_policy_version=CONFIDENCE_POLICY_VERSION,
        high_confidence_requires_intron_support=no_utr,
        options=dict(
            original.get("options", {}),
            no_utr=no_utr,
            output_prefix=str(outdir / prefix),
        ),
        processed_transcript_count=len(frame),
        selected_transcript_count=original.get("selected_transcript_count", len(frame)),
        confidence_counts=result["confidence_counts"],
        outputs={key: str(value) for key, value in paths.items()},
        rescoring={
            "source_evidence": str(source),
            "source_evidence_sha256": file_sha256(source),
            "source_manifest": str(manifest_path) if manifest_path.exists() else None,
            "extraction_script_sha256": original.get("script_sha256"),
            "source_confidence_policy_version": original.get(
                "confidence_policy_version", 1
            ),
            "old_confidence_counts": result["old_confidence_counts"],
            "class_changed_count": result["class_changed_count"],
            "unchanged_companions": "Absolute symlinks to the original run; retain original files.",
        },
    )
    write_rescored_outputs(paths, result, manifest, companions, gff)
    return manifest


def run_cli(args):
    """Handle single runs or preserve species/method layout in batch rescoring."""
    if args.evidence:
        manifest = rescore(args.evidence, args.output_dir, args.no_utr)
        print(manifest["confidence_counts"])
        return manifest
    root = args.input_root.resolve()
    if args.output_dir.resolve() == root or root in args.output_dir.resolve().parents:
        raise ValueError("Choose an output root outside the input root")
    paths = sorted(root.glob("*/*/*.evidence.tsv"))
    if not paths:
        raise ValueError("No SPECIES/METHOD/*.evidence.tsv inputs found")
    manifests = []
    for path in paths:
        manifest = rescore(
            path, args.output_dir / path.relative_to(root).parent, args.no_utr
        )
        print(f"{path}: {manifest['confidence_counts']}")
        manifests.append(manifest)
    return manifests


def main():
    """Support execution as an installed module."""
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    run_cli(parser.parse_args())


if __name__ == "__main__":
    main()

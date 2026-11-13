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
"""Run high-confidence annotation QC through parser, metric and report layers."""
from __future__ import annotations
import argparse
import hashlib
from collections import Counter
import json
import subprocess
import sys
from pathlib import Path

from ensembl.genes.annotation_qc.config.high_confidence import (
    DEFAULT_MITOCHONDRIAL_CONTIGS,
    CONFIDENCE_POLICY_VERSION,
    OUTPUT_SUFFIXES,
)
from ensembl.genes.annotation_qc.parsers.high_confidence import (
    Genome,
    parse_high_confidence_annotation,
    parse_canonical_assignments,
    parse_legacy_diamond_alignments,
    read_junctions,
    read_homology,
    read_transcript_region,
    file_sha256,
)
from ensembl.genes.annotation_qc.metrics.high_confidence import (
    TranslationExceptionError,
    biopython_api,
    select_representative_transcripts,
    validate_model_bounds,
    assess_transcript,
    score_diamond_alignments,
    build_confidence_tables,
    compute_run_diagnostics,
)
from ensembl.genes.annotation_qc.reports.high_confidence import (
    HighConfidenceReports,
    write_homology,
)


def add_arguments(parser):
    """Register scoring inputs, selection policies and output options."""
    parser.add_argument(
        "--annotation",
        "--gtf",
        dest="gtf",
        required=True,
        type=Path,
        help="GTF or GFF3 annotation (option name retained for compatibility)",
    )
    parser.add_argument(
        "--genome",
        required=True,
        type=Path,
        help="Plain regularly wrapped genome FASTA",
    )
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--canonical-map", type=Path)
    parser.add_argument(
        "--input-mode",
        choices=["auto", "strict"],
        default="auto",
        help="Auto normalizes missing biotypes/exons/phases; strict requires supplied data",
    )
    parser.add_argument(
        "--canonical-policy",
        choices=["auto", "strict", "longest-cds"],
        default="auto",
        help=(
            "Auto uses canonical tags or falls back to longest CDS; strict "
            "refuses ambiguous choices"
        ),
    )
    parser.add_argument(
        "--junctions", type=Path, help="Stranded intron read-count TSV; see README"
    )
    parser.add_argument("--min-junction-reads", type=int, default=1)
    parser.add_argument(
        "--no-utr",
        "--no_utr",
        dest="no_utr",
        action="store_true",
        help="Disable the UTR-based high route and require intron support for every high call",
    )
    parser.add_argument(
        "--mitochondrial-contigs",
        nargs="+",
        default=list(DEFAULT_MITOCHONDRIAL_CONTIGS),
        help="Exact contig IDs using vertebrate mitochondrial table 2 "
        "(default: MT chrM chrMT M); others use table 1 unless GTF declares a table",
    )
    parser.add_argument(
        "--translation-exception-policy",
        choices=["error", "skip"],
        default="error",
        help="Unresolved annotated translations: stop (default) or skip and report",
    )
    parser.add_argument(
        "--invalid-model-policy",
        choices=["error", "skip"],
        default="error",
        help="GTF/FASTA coordinate mismatches: stop with report (default), or skip affected models",
    )
    homology = parser.add_mutually_exclusive_group()
    homology.add_argument("--diamond-results", type=Path)
    homology.add_argument(
        "--diamond-script", type=Path, help="Existing get_diamond_coverage.py helper"
    )
    homology.add_argument(
        "--diamond-alignments",
        type=Path,
        help="Nine-column blastp output; reproduce helper metrics",
    )
    homology.add_argument(
        "--run-diamond",
        action="store_true",
        help="Run installed diamond and reproduce helper metrics",
    )
    parser.add_argument("--protein-db", type=Path)


def register(subparsers):
    """Register high-confidence on the shared annotation-qc command."""
    parser = subparsers.add_parser(
        "high-confidence",
        help="Assess gene confidence from annotation, genome and external evidence.",
    )
    add_arguments(parser)
    parser.set_defaults(func=run_cli)


def implementation_hashes():
    """Fingerprint every implementation layer, not just the runner file."""

    root = Path(__file__).parent.parent
    relative = [
        "config/high_confidence.py",
        "parsers/high_confidence.py",
        "parsers/_high_confidence_records.py",
        "metrics/high_confidence.py",
        "reports/high_confidence.py",
        "runners/high_confidence.py",
    ]
    hashes = {name: file_sha256(root / name) for name in relative}
    combined = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    return combined, hashes


# The orchestration/layout keeps related report operations together.
# pylint: disable-next=too-many-branches,too-many-locals,too-many-statements
def run(args):
    """Assemble the QC workflow; calculations and formatting are delegated."""
    if args.min_junction_reads < 1:
        raise ValueError("--min-junction-reads must be at least 1")
    if bool(args.diamond_script or args.run_diamond) != bool(args.protein_db):
        raise ValueError(
            "--protein-db is required with --diamond-script or --run-diamond only"
        )
    inputs = {
        key: value.resolve()
        for key, value in vars(args).items()
        if isinstance(value, Path) and key != "output_prefix"
    }
    for path in inputs.values():
        if not path.is_file():
            raise ValueError(f"Input file does not exist: {path}")
    prefix = str(args.output_prefix.resolve())
    paths = {key: Path(prefix + suffix) for key, suffix in OUTPUT_SUFFIXES.items()}
    if args.diamond_script or args.diamond_alignments or args.run_diamond:
        paths["homology"] = Path(prefix + ".prots.res")
    if args.run_diamond:
        paths["alignments"] = Path(prefix + ".diamond.tsv")
    if {path.resolve() for path in paths.values()} & set(inputs.values()):
        raise ValueError(
            "An output would overwrite an input file; choose another prefix"
        )
    _, _, biopython_version = biopython_api()
    annotation = parse_high_confidence_annotation(args.gtf, args.input_mode)
    policy = (
        "strict"
        if args.input_mode == "strict" and args.canonical_policy == "auto"
        else args.canonical_policy
    )
    genes, selected = select_representative_transcripts(
        annotation,
        parse_canonical_assignments(args.canonical_map),
        policy,
        args.input_mode,
    )
    junctions = read_junctions(args.junctions)
    reports = HighConfidenceReports(paths)
    results, exceptions, normalization = [], [], []
    genome = Genome(args.genome)
    try:
        issues, invalid_models = validate_model_bounds(selected, genome)
        reports.write_validation(issues, args.invalid_model_policy)
        if invalid_models and args.invalid_model_policy == "error":
            raise ValueError(
                f"GTF/FASTA coordinate mismatch in {len(invalid_models)} selected models; "
                f'check the assembly and contig bounds. See {paths["input_validation"]}. '
                "Use --invalid-model-policy skip only to explicitly exclude these models."
            )
        processable = [
            entry
            for entry in selected
            if (entry[0].id, entry[1].id) not in invalid_models
        ]
        if not processable:
            raise ValueError(
                f'No models remain after coordinate validation; see {paths["input_validation"]}'
            )
        for gene in genes.values():
            normalization.extend(
                # Keyword construction rejects duplicate export-schema keys.
                # pylint: disable-next=use-dict-literal
                dict(gene_id=gene.id, transcript_id=".", **note)
                for note in gene.normalization
            )
        for gene, tx, selection_source, biotype_source in processable:
            normalization.extend(
                # Keyword construction rejects duplicate export-schema keys.
                # pylint: disable-next=use-dict-literal
                dict(gene_id=gene.id, transcript_id=tx.id, **note)
                for note in tx.normalization
            )
            if selection_source == "longest_cds_fallback":
                normalization.append(
                    {
                        "gene_id": gene.id,
                        "transcript_id": tx.id,
                        "action": "representative_transcript_inferred",
                        "details": (
                            "No canonical assignment; selected longest CDS, then "
                            "exon length, then transcript ID"
                        ),
                    }
                )
            try:
                region = read_transcript_region(genome, tx)
                result = assess_transcript(
                    tx,
                    region,
                    junctions,
                    args.min_junction_reads,
                    mitochondrial_contigs=args.mitochondrial_contigs,
                    input_mode=args.input_mode,
                )
            except TranslationExceptionError as error:
                exceptions.append(
                    {
                        "gene_id": gene.id,
                        "transcript_id": tx.id,
                        "policy": args.translation_exception_policy,
                        "reason": str(error),
                    }
                )
                reports.write_extraction_diagnostics(normalization, exceptions)
                if args.translation_exception_policy == "error":
                    raise
                print(f"Skipped {tx.id}: {error}", file=sys.stderr)
                continue
            except ValueError as error:
                raise ValueError(
                    f"Gene {gene.id}, transcript {tx.id}: {error}"
                ) from error
            normalization.extend(
                # Keyword construction rejects duplicate export-schema keys.
                # pylint: disable-next=use-dict-literal
                dict(gene_id=gene.id, transcript_id=tx.id, **note)
                for note in result["normalization"]
            )
            result.update(
                gene=gene,
                tx=tx,
                selection_source=selection_source,
                biotype_source=biotype_source,
            )
            results.append(result)
    finally:
        genome.close()
    reports.write_extraction_diagnostics(normalization, exceptions)
    if not results:
        raise ValueError(
            f'No transcripts could be processed; see {paths["translation_exceptions"]}'
        )
    reports.write_sequences(results)
    result_path, alignment_path = args.diamond_results, args.diamond_alignments
    diamond_version = None
    if args.run_diamond:
        diamond_version = subprocess.run(
            ["diamond", "version"],
            text=True,
            capture_output=True,
            check=True,
        ).stdout.strip()
        paths["alignments"].unlink(missing_ok=True)
        subprocess.run(
            [
                "diamond",
                "blastp",
                "-d",
                str(args.protein_db.resolve()),
                "-q",
                str(paths["proteins"]),
                "--masking",
                "0",
                "-f",
                "6",
                "qseqid",
                "sseqid",
                "qstart",
                "qend",
                "slen",
                "pident",
                "nident",
                "length",
                "gaps",
                "-o",
                str(paths["alignments"]),
            ],
            check=True,
        )
        alignment_path = paths["alignments"]
    if args.diamond_script:
        paths["homology"].unlink(missing_ok=True)
        subprocess.run(
            [
                sys.executable,
                str(args.diamond_script.resolve()),
                str(paths["proteins"]),
                str(args.protein_db.resolve()),
                str(paths["homology"]),
            ],
            check=True,
        )
        result_path = paths["homology"]
    if alignment_path:
        alignments = parse_legacy_diamond_alignments(alignment_path)
        homology = score_diamond_alignments(
            alignments, {r["tx"].id: len(r["protein"]) for r in results}
        )
        result_path = paths["homology"]
        write_homology(result_path, homology)
    else:
        homology = read_homology(result_path, {tx.id for _, tx, _, _ in selected})
    tables = build_confidence_tables(
        results, homology, args.no_utr, result_path is not None
    )
    reports.write_evidence(tables)
    code_hash, hashes = implementation_hashes()
    manifest = {
        "schema_version": 1,
        "script_sha256": code_hash,
        "implementation_hashes": hashes,
        "confidence_policy_version": CONFIDENCE_POLICY_VERSION,
        "high_confidence_requires_intron_support": args.no_utr,
        "inputs": {
            key: {"path": str(path), "sha256": file_sha256(path)}
            for key, path in inputs.items()
        },
        "options": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
            if not callable(value)
        },
        "gene_count": len(genes),
        "selected_transcript_count": len(selected),
        "processed_transcript_count": len(results),
        "skipped_translation_exception_count": len(exceptions),
        "skipped_invalid_model_count": len(invalid_models),
        "input_validation_issue_count": len(issues),
        "biopython_version": biopython_version,
        "python_version": sys.version.split()[0],
        "effective_canonical_policy": policy,
        "confidence_counts": tables["confidence_counts"],
        "normalization_counts": dict(Counter((row["action"] for row in normalization))),
        "diamond_version": diamond_version,
        "homology_result": (
            None
            if result_path is None
            else {
                "path": str(result_path.resolve()),
                "sha256": file_sha256(result_path),
            }
        ),
        "outputs": {key: str(value) for key, value in paths.items()},
        "limitations": [
            (
                "Default UTR-aware scoring preserves the Perl branches; no_utr "
                "policy v2 requires intron support for every high call"
            ),
            (
                "Intronless models retain vacuous intron-support checks; absent "
                "RNA data cannot satisfy support-dependent rules for spliced "
                "models"
            ),
            "DIAMOND target coverage uses query span for legacy parity and may exceed 100 percent",
            (
                "The original ATG plus TAA/TAG/TGA completeness heuristic is "
                "independent of translation table"
            ),
            (
                "Auto-mode normalizations and exact CDS-junction support differ "
                "from some original database evidence exports"
            ),
        ],
    }
    manifest.update(compute_run_diagnostics(results))
    reports.write_manifest(manifest)
    return manifest


def run_cli(args):
    """Print workflow status and use a nonzero exit status on input/tool failures."""
    try:
        manifest = run(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
    print(
        f"Processed {manifest['processed_transcript_count']} genes: {manifest['confidence_counts']}"
    )
    return manifest


def main():
    """Support python -m ensembl.genes.annotation_qc.runners.high_confidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    run_cli(parser.parse_args())


if __name__ == "__main__":
    main()

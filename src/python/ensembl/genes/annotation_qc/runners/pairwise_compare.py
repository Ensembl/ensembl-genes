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
"""
Runner to compare a completed annotation (query) with a trusted reference.

Parses both annotations, harmonises sequence names, filters the reference,
classifies every reference and query gene and writes the summary, detail,
label, audit and manifest files. The reference is used for evaluation only.

Usage:
        annotation-qc pairwise-compare \
            --query finalise/consensus.gff3 --reference reference.gff3.gz \
            --genome genome.fa --evaluation-mode protein_coding --outdir comparison
"""

import argparse
import hashlib
import importlib.metadata
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    build_gene_models,
    classify_loci,
    gene_span_excess,
)
from ensembl.genes.annotation_qc.metrics.pairwise.matching import (
    one_to_one_exact_cds,
)
from ensembl.genes.annotation_qc.metrics.pairwise.selection import (
    EVALUATION_MODES,
    TRANSCRIPT_SELECTIONS,
    apply_evaluation_mode,
    drop_seqnames,
    filter_by_biotype,
    select_transcripts,
    subset_to_regions,
    summarise_annotation,
)
from ensembl.genes.annotation_qc.metrics.pairwise.summary import (
    query_transcript_labels,
    summarise_comparison,
    summarise_intron_support,
    summarise_span_excess,
)
from ensembl.genes.annotation_qc.parsers.annotation import (
    ANNOTATION_FORMATS,
    parse_annotation_for_comparison,
)
from ensembl.genes.annotation_qc.parsers.regions import load_regions_file, parse_region
from ensembl.genes.annotation_qc.parsers.seqnames import (
    apply_seqname_mapping,
    build_seqname_mapping,
    check_seqnames_against_genome,
)
from ensembl.genes.annotation_qc.parsers.sequence import parse_sequence_lengths
from ensembl.genes.annotation_qc.parsers.tool_output import (
    find_gmb_handover_manifest,
    read_evidence_attribution,
)
from ensembl.genes.annotation_qc.reports import pairwise_compare as report

GENOME_MISMATCH_POLICIES = ("error", "warn", "exclude")


class ComparisonInputError(ValueError):
    """Inputs are inconsistent; nothing has been compared."""


def _log(message: str) -> None:
    print(message, flush=True)


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: str | None) -> dict | None:
    if not path:
        return None
    return {
        "path": str(Path(path).resolve()),
        "size": os.path.getsize(path),
        "sha256": _sha256(path),
    }


def _code_version() -> dict:
    version = {"ensembl_genes": importlib.metadata.version("ensembl-genes")}
    source = Path(__file__).resolve().parent
    try:
        commit = subprocess.run(
            ["git", "-C", str(source), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(source), "status", "--porcelain", "--", "."],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        version.update(git_commit=commit, git_uncommitted_changes=bool(dirty))
    except (OSError, subprocess.CalledProcessError):
        version["git_commit"] = None
    for package in ("pandas", "pyranges1"):
        version[package] = importlib.metadata.version(package)
    version["python"] = platform.python_version()
    return version


def _gmb_provenance(query_path: str, query_sha256: str) -> dict | None:
    found = find_gmb_handover_manifest(query_path)
    if not found:
        return None
    manifest_path, manifest = found
    listed = {
        name: entry.get("sha256") for name, entry in manifest.get("outputs", {}).items()
    }
    matches = [name for name, digest in listed.items() if digest == query_sha256]
    return {
        "handover_manifest": manifest_path,
        "build_dir": manifest.get("build_dir"),
        "finalise_dir": manifest.get("output_dir"),
        "genome_path": manifest.get("genome_path"),
        "preset": manifest.get("preset"),
        "timestamp": manifest.get("timestamp"),
        "query_listed_as": matches[0] if matches else None,
        "query_matches_manifest": bool(matches),
    }


def _parse(path: str, format_hint: str, label: str):
    try:
        annotation = parse_annotation_for_comparison(path, format_hint)
    except ValueError as error:
        raise ComparisonInputError(f"{label} annotation {path}: {error}") from error
    diagnostics = annotation.attrs["parse_diagnostics"]
    _log(
        f"  {label}: {diagnostics['genes']} genes, {diagnostics['transcripts']} transcripts, "
        f"{diagnostics['exon_rows']} exons, {diagnostics['cds_rows']} CDS "
        f"(format: {diagnostics['format']})"
    )
    for key in (
        "namespaced_gene_ids",
        "synthesised_genes",
        "multi_parent_children",
        "cds_only_transcripts",
    ):
        if diagnostics.get(key):
            _log(f"  {label}: {key} = {diagnostics[key]}")
    return annotation, diagnostics


def _harmonise_seqnames(ref, query, args, mapping) -> tuple:
    """Apply the seqname map and check names/coordinates. Returns (ref, query, report)."""
    ref, ref_map = apply_seqname_mapping(ref, mapping)
    query, query_map = apply_seqname_mapping(query, mapping)
    checks = {"reference_mapping": ref_map, "query_mapping": query_map, "excluded": {}}

    shared = set(ref["Chromosome"]) & set(query["Chromosome"])
    if not shared:
        raise ComparisonInputError(
            "Reference and query share no sequence names after mapping "
            f"(reference e.g. {sorted(set(ref['Chromosome']))[:3]}, "
            f"query e.g. {sorted(set(query['Chromosome']))[:3]}). "
            "Supply --seqname-map (annotation name -> genome name) or --assembly-report."
        )

    if not args.genome:
        only_ref = sorted(set(ref["Chromosome"]) - shared)
        if only_ref:
            _log(
                f"  NOTE: {len(only_ref)} reference sequence(s) have no query features: {only_ref[:10]}"
            )
        return ref, query, checks

    lengths = parse_sequence_lengths(args.genome)
    for label, frame in (("reference", ref), ("query", query)):
        result = check_seqnames_against_genome(frame, lengths)
        checks[f"{label}_genome_check"] = result
        missing, beyond = result["missing_seqnames"], result["out_of_bounds_seqnames"]
        problems = []
        if missing:
            problems.append(
                f"{len(missing)} {label} sequence(s) are not in the genome FASTA "
                f"({sum(missing.values())} genes): {dict(list(missing.items())[:10])}"
            )
        if beyond:
            problems.append(
                f"{result['out_of_bounds_features']} {label} feature(s) on {len(beyond)} "
                f"sequence(s) end beyond the FASTA sequence ({sum(beyond.values())} genes on "
                f"them), e.g. {result['out_of_bounds_examples']}; the annotation may be on a "
                "different assembly version"
            )
        if not problems:
            continue
        message = ". ".join(problems)
        if args.genome_mismatch == "error":
            raise ComparisonInputError(
                message
                + ". Fix the genome or seqname map, or pass --genome-mismatch exclude "
                "to drop these sequences from the comparison (warn keeps them, as gmb-compare did)."
            )
        _log(f"  WARNING: {message}")
        if args.genome_mismatch == "exclude":
            checks["excluded"][label] = {**missing, **beyond}

    # Excluded sequences are removed from both annotations so that neither side
    # gains Missed or Novel genes from them.
    excluded = {name for names in checks["excluded"].values() for name in names}
    if excluded:
        ref, query = drop_seqnames(ref, excluded), drop_seqnames(query, excluded)
        _log(f"  Excluded sequences from both annotations: {sorted(excluded)}")
    return ref, query, checks


def _filter_reference(ref, args) -> tuple:
    before = summarise_annotation(ref)
    ref = apply_evaluation_mode(ref, args.evaluation_mode)
    gene_biotypes = _split(args.reference_gene_biotypes)
    transcript_biotypes = _split(args.reference_transcript_biotypes)
    if gene_biotypes or transcript_biotypes:
        ref = filter_by_biotype(ref, gene_biotypes, transcript_biotypes)
    ref = select_transcripts(ref, args.reference_transcript_selection)
    after = summarise_annotation(ref)
    audit = {
        "evaluation_mode": args.evaluation_mode,
        "transcript_selection": args.reference_transcript_selection,
        "gene_biotypes_filter": args.reference_gene_biotypes,
        "transcript_biotypes_filter": args.reference_transcript_biotypes,
        "pre_filter_genes": before["genes"],
        "pre_filter_transcripts": before["transcripts"],
        "post_filter_genes": after["genes"],
        "post_filter_transcripts": after["transcripts"],
        "post_filter_exons": after["exons"],
        "post_filter_cds": after["cds"],
        "cds_containing_transcripts": after["cds_containing_transcripts"],
        "canonical_transcripts": after["canonical_transcripts"],
        "gene_biotype_counts": after["gene_biotype_counts"],
        "transcript_biotype_counts": after["transcript_biotype_counts"],
    }
    _log(
        f"  Reference after filtering: {after['genes']} of {before['genes']} genes, "
        f"{after['transcripts']} of {before['transcripts']} transcripts"
    )
    return ref, audit


def _split(value: str | None) -> list[str] | None:
    return (
        [item.strip() for item in value.split(",") if item.strip()] if value else None
    )


def _regions(args) -> list:
    regions = [parse_region(text) for text in args.region or []]
    if args.regions_file:
        regions.extend(load_regions_file(args.regions_file))
    return regions


def run_pairwise_compare(args) -> dict:
    """
    Run the comparison described by parsed CLI arguments.
    Args:
            args: argparse.Namespace from this runner's parser
    Returns:
            The summary dict written to comparison_summary.json
    """
    started = time.time()
    for option in (
        "query",
        "reference",
        "genome",
        "seqname_map",
        "assembly_report",
        "evidence_attribution",
        "regions_file",
    ):
        path = getattr(args, option)
        if path and not os.path.isfile(path):
            raise ComparisonInputError(
                f"--{option.replace('_', '-')}: file not found: {path}"
            )
    regions = _regions(args)
    if args.plots_per_category:
        from ensembl.genes.annotation_qc.reports import (
            pairwise_plots,
        )  # optional matplotlib

        try:
            pairwise_plots.import_pyplot()
        except ImportError as error:
            raise ComparisonInputError(str(error)) from error

    _log("Loading annotations...")
    mapping = build_seqname_mapping(args.seqname_map, args.assembly_report)
    ref, ref_diag = _parse(args.reference, args.reference_format, "Reference")
    query, query_diag = _parse(args.query, args.query_format, "Query")
    ref, query, seqname_checks = _harmonise_seqnames(ref, query, args, mapping)

    _log("Filtering reference...")
    # The unfiltered reference only labels Novel query genes (e.g. pseudogene loci).
    ref_context = ref
    ref, audit = _filter_reference(ref, args)
    if seqname_checks["excluded"].get("reference"):
        audit["excluded_seqnames"] = seqname_checks["excluded"]["reference"]
    if args.query_transcript_selection != "all":
        query = select_transcripts(query, args.query_transcript_selection)
        _log(
            f"  Query transcripts kept ({args.query_transcript_selection}): "
            f"{int((query['Feature'] == 'transcript').sum())}"
        )
    if regions:
        ref = subset_to_regions(ref, regions)
        query = subset_to_regions(query, regions)
        ref_context = subset_to_regions(ref_context, regions)
        _log(f"  Restricted to {len(regions)} region(s)")

    _log("Classifying loci...")
    ref_models, query_models = build_gene_models(ref), build_gene_models(query)
    ref_results, query_results = classify_loci(
        ref_models, query_models, reference_context=build_gene_models(ref_context)
    )
    summary = summarise_comparison(ref_results, query_results)
    summary["intron_support"] = summarise_intron_support(ref_models, query_models)
    summary["cds_exact_one_to_one"], exact_pairs = one_to_one_exact_cds(
        ref_models, query_models
    )
    labels = query_transcript_labels(query_results)
    summary["filter_info"] = {
        key: audit[key]
        for key in (
            "evaluation_mode",
            "transcript_selection",
            "pre_filter_genes",
            "pre_filter_transcripts",
            "post_filter_genes",
            "post_filter_transcripts",
        )
    }
    for key, name in (
        ("gene_biotypes_filter", "gene_biotypes"),
        ("transcript_biotypes_filter", "transcript_biotypes"),
    ):
        if audit[key]:
            summary["filter_info"][name] = audit[key]
    if regions:
        summary["subset_regions"] = [str(region) for region in regions]
    summary["gene_span_checks"] = {
        "reference": summarise_span_excess(gene_span_excess(ref_models)),
        "query": summarise_span_excess(gene_span_excess(query_models)),
    }

    _log("Writing reports...")
    os.makedirs(args.outdir, exist_ok=True)
    report.write_filter_audit(audit, args.outdir)
    report.write_summary(summary, args.outdir)
    report.write_details(ref_results, query_results, args.outdir)
    report.write_transcript_labels(labels, args.outdir)
    report.write_split_merge(ref_results, query_results, args.outdir)
    report.write_exact_cds_pairs(exact_pairs, args.outdir)
    unlabelled = None
    if args.evidence_attribution:
        attribution = read_evidence_attribution(args.evidence_attribution)
        unlabelled = report.write_evidence_attribution_labels(
            attribution, labels, args.outdir
        )
        if unlabelled:
            _log(
                f"  WARNING: {unlabelled} evidence attribution rows have no comparison label"
            )
    if args.plots_per_category:
        pairwise_plots.write_locus_plots(
            ref_results, query_results, ref, query, args.outdir, args.plots_per_category
        )

    query_record = _file_record(args.query)
    manifest = {
        "command": "annotation-qc pairwise-compare",
        "argv": sys.argv,
        "created": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": round(time.time() - started, 1),
        "inputs": {
            "query": query_record,
            "reference": _file_record(args.reference),
            "genome": _file_record(args.genome),
            "seqname_map": _file_record(args.seqname_map),
            "assembly_report": _file_record(args.assembly_report),
            "evidence_attribution": _file_record(args.evidence_attribution),
            "regions_file": _file_record(args.regions_file),
        },
        "gmb_finalise": _gmb_provenance(args.query, query_record["sha256"]),
        "options": {key: value for key, value in vars(args).items() if key != "func"},
        "coordinates": "internally zero-based half-open; reported one-based inclusive",
        "parse_diagnostics": {"reference": ref_diag, "query": query_diag},
        "seqname_checks": seqname_checks,
        "evidence_attribution_unlabelled_rows": unlabelled,
        "code_version": _code_version(),
    }
    report.write_manifest(manifest, args.outdir)

    _log("\n" + report.format_console_summary(summary))
    _log(f"\nResults written to {args.outdir}")
    return summary


def _add_args(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--query",
        required=True,
        help="Completed annotation to evaluate (GFF3/GTF/Tiberius, may be .gz).",
    )
    required.add_argument(
        "--reference",
        required=True,
        help="Trusted reference annotation (GFF3/GTF, may be .gz).",
    )
    required.add_argument("--outdir", required=True, help="Output directory.")

    names = parser.add_argument_group("sequence names and coordinates")
    names.add_argument(
        "--genome",
        default=None,
        help="Genome FASTA (plain or gzip). Checks that every sequence name "
        "exists and every feature lies within its sequence.",
    )
    names.add_argument(
        "--seqname-map",
        default=None,
        help="Two-column TSV/CSV: annotation seqname -> genome seqname. Applied "
        "to reference and query; overrides --assembly-report.",
    )
    names.add_argument(
        "--assembly-report",
        default=None,
        help="NCBI assembly report mapping GenBank accessions to chromosome names.",
    )
    names.add_argument(
        "--genome-mismatch",
        choices=GENOME_MISMATCH_POLICIES,
        default="error",
        help="With --genome: what to do with annotation sequences that are absent "
        "from the FASTA or have features beyond its end. error (default) "
        "stops; exclude drops those sequences from both annotations; warn "
        "keeps them, as gmb-compare did.",
    )

    reference = parser.add_argument_group("reference selection")
    reference.add_argument(
        "--evaluation-mode",
        choices=EVALUATION_MODES,
        default="all",
        help="Reference preset (default: all).",
    )
    reference.add_argument(
        "--reference-transcript-selection",
        choices=TRANSCRIPT_SELECTIONS,
        default="all",
        help="Reference transcripts per gene (default: all).",
    )
    reference.add_argument(
        "--reference-gene-biotypes",
        default=None,
        help="Comma-separated reference gene biotypes to keep.",
    )
    reference.add_argument(
        "--reference-transcript-biotypes",
        default=None,
        help="Comma-separated reference transcript biotypes to keep.",
    )

    reference.add_argument(
        "--query-transcript-selection",
        choices=TRANSCRIPT_SELECTIONS,
        default="all",
        help="Query transcripts per gene (default: all). Use canonical with "
        "GMB consensus.canonical_annotated.gff3 to score only the "
        "Ensembl_canonical-tagged transcript.",
    )

    other = parser.add_argument_group("other options")
    other.add_argument(
        "--query-format",
        choices=("auto",) + ANNOTATION_FORMATS,
        default="auto",
        help="Query format (default: detect from content).",
    )
    other.add_argument(
        "--reference-format",
        choices=("auto",) + ANNOTATION_FORMATS,
        default="auto",
        help="Reference format (default: detect from content).",
    )
    other.add_argument(
        "--region",
        action="append",
        default=None,
        help="Restrict to genes overlapping seqname[:start-end] (one-based, "
        "comparison seqnames). Repeatable.",
    )
    other.add_argument(
        "--regions-file", default=None, help="File with one region per line."
    )
    other.add_argument(
        "--evidence-attribution",
        default=None,
        help="GMB build/evidence_attribution.tsv to label with query results.",
    )
    other.add_argument(
        "--plots-per-category",
        type=int,
        default=0,
        help="Locus plots per classification in qc/ (default: 0; needs matplotlib).",
    )


def register(subparsers):
    """Register this runner as a subcommand on the shared CLI subparsers object."""
    parser = subparsers.add_parser(
        "pairwise-compare",
        help="Compare a completed annotation with a trusted reference annotation.",
        description=(
            "Compare a completed annotation (query) with a trusted reference. Reports "
            "reference-side sensitivity (Exact/Partial/Structural/Missed, coordinate-exact "
            "CDS/exons, intron chains, splits) and query-side precision (matched, novel by "
            "overlap context, merges). The reference is used for evaluation only."
        ),
    )
    _add_args(parser)
    parser.set_defaults(func=_run)


def _run(args):
    try:
        return run_pairwise_compare(args)
    except ComparisonInputError as error:
        sys.exit(f"ERROR: {error}")


def main():
    """Entry point for standalone execution of the pairwise comparison runner."""
    parser = argparse.ArgumentParser(
        description="Compare an annotation with a reference."
    )
    _add_args(parser)
    _run(parser.parse_args())


if __name__ == "__main__":
    main()

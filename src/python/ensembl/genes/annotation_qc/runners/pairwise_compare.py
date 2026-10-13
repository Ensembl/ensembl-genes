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
import json
import os
import platform
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from ensembl.genes.annotation_qc.metrics.pairwise import contract
from ensembl.genes.annotation_qc.metrics.pairwise.classify import (
    RECIPROCAL_THRESHOLD,
    SPLIT_MERGE_MIN_SHARED_FRACTION,
    build_gene_models,
    classify_loci,
    gene_span_excess,
    transcripts_outside_gene,
)
from ensembl.genes.annotation_qc.metrics.pairwise.evaluation_audit import (
    ExclusionTracker,
    canonical_tag_count,
    selection_labels,
    transcript_ids,
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
from ensembl.genes.annotation_qc.metrics.pairwise.stop_codons import (
    GENETIC_CODE_STOPS,
    harmonise_stop_codons,
    parse_genetic_code_policy,
)
from ensembl.genes.annotation_qc.metrics.pairwise.summary import (
    query_transcript_labels,
    summarise_comparison,
    summarise_intron_support,
    summarise_span_excess,
    summarise_transcripts_outside,
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
from ensembl.genes.annotation_qc.parsers.sequence import (
    iter_sequences,
    parse_sequence_lengths,
)
from ensembl.genes.annotation_qc.parsers.tool_output import (
    find_gmb_handover_manifest,
    read_evidence_attribution,
)
from ensembl.genes.annotation_qc.reports import pairwise_compare as report
from ensembl.genes.annotation_qc.reports import pairwise_contract as contract_report

GENOME_MISMATCH_POLICIES = ("error", "warn", "exclude")
EXPORT_MODES = ("ids", "full")
_SCOPE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")


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


def _compression(path: str) -> str:
    """none, gzip or bgzf (gzip with the BGZF extra field), from the content."""
    with open(path, "rb") as handle:
        head = handle.read(18)
    if head[:2] != b"\x1f\x8b":
        return "none"
    return (
        "bgzf" if len(head) >= 14 and head[3] & 4 and head[12:14] == b"BC" else "gzip"
    )


def _file_record(path: str | None) -> dict | None:
    if not path:
        return None
    return {
        "path": str(Path(path).resolve()),
        "size": os.path.getsize(path),
        "sha256": _sha256(path),
        "compression": _compression(path),
    }


# Scheme 1 (contract 1.0 before this change) selected modules by import name
# (ensembl.genes.annotation_qc.*) and so missed modules that the annotation-qc
# entry point imports under top-level names (cli.py puts its directory on
# sys.path: runners.pairwise_compare, parsers.interpro, ...), including this
# runner. Scheme 2 selects every loaded module file inside the package.
FINGERPRINT_SCHEME = "annotation_qc_loaded_files_by_path/2"


def _implementation() -> dict:
    """
    Content fingerprint of the annotation_qc module files this run loaded.

    Modules are selected by file location (inside the annotation_qc package),
    whatever name they were imported under. Identifies the code even when it is
    not a clean Git checkout (uncommitted edits) or not a checkout at all (an
    installed wheel).
    """
    package = Path(__file__).resolve().parent.parent
    files = {}
    for module in list(sys.modules.values()):
        path = getattr(module, "__file__", None)
        if not path:
            continue
        resolved = Path(path).resolve()
        if resolved.is_relative_to(package):
            files[str(resolved.relative_to(package))] = _sha256(str(resolved))
    files = dict(sorted(files.items()))
    versions = {
        name: importlib.metadata.version(name) for name in ("pandas", "pyranges1")
    }
    versions["python"] = platform.python_version()
    text = json.dumps({"files": files, "versions": versions}, sort_keys=True)
    return {
        "scheme": FINGERPRINT_SCHEME,
        "fingerprint": hashlib.sha256(text.encode()).hexdigest(),
        "source_files": files,
        "versions": versions,
    }


def _code_version() -> dict:
    version = {"ensembl_genes": importlib.metadata.version("ensembl-genes")}
    # The whole annotation_qc package (this runner imports from all of it).
    source = Path(__file__).resolve().parent.parent
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


def _harmonise_seqnames(ref, query, args, mapping, lengths) -> tuple:
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


def _select(df, mode: str, label: str, tracker: ExclusionTracker) -> tuple:
    """Transcript selection with an explicit canonical check and labels."""
    if mode == "canonical" and canonical_tag_count(df) == 0:
        raise ComparisonInputError(
            f"{label}: canonical transcript selection was requested but no "
            "transcript carries the Ensembl_canonical tag, so every gene would "
            "silently fall back to the longest CDS. Choose all or longest_cds."
        )
    selected = select_transcripts(df, mode)
    labels = selection_labels(df, transcript_ids(selected), mode)
    tracker.step(selected, f"not_selected:{mode}")
    return selected, labels


def _filter_reference(ref, args, tracker: ExclusionTracker) -> tuple:
    """Evaluation mode, biotype filters and selection, each step tracked."""
    before = summarise_annotation(ref)
    mode_labels = {}
    if args.evaluation_mode == "canonical":
        # The canonical preset is a transcript selection.
        ref, mode_labels = _select(ref, "canonical", "reference", tracker)
    else:
        ref = apply_evaluation_mode(ref, args.evaluation_mode)
        tracker.step(ref, f"evaluation_mode:{args.evaluation_mode}")
    gene_biotypes = _split(args.reference_gene_biotypes)
    transcript_biotypes = _split(args.reference_transcript_biotypes)
    # filter_by_biotype applies the gene filter, then the transcript filter; the
    # two calls below give the same table and record which filter dropped what.
    if gene_biotypes:
        ref = filter_by_biotype(ref, gene_biotypes, None)
        tracker.step(ref, "gene_biotype_filter")
    if transcript_biotypes:
        ref = filter_by_biotype(ref, None, transcript_biotypes)
        tracker.step(ref, "transcript_biotype_filter")
    selection = args.reference_transcript_selection
    ref, labels = _select(ref, selection, "reference", tracker)
    if selection == "all" and mode_labels:
        labels = {tid: mode_labels.get(tid, "all") for tid in labels}
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
    return ref, audit, labels


def _stop_codon_policy(args):
    """Validate the stop-codon options; return the genetic-code policy or None."""
    if not args.add_stop_codon:
        if args.stop_codon_genetic_code != 1 or args.stop_codon_sequence_code:
            raise ComparisonInputError(
                "--stop-codon-genetic-code and --stop-codon-sequence-code need "
                "--add-stop-codon"
            )
        return None
    if not args.genome:
        raise ComparisonInputError("--add-stop-codon needs --genome")
    try:
        return parse_genetic_code_policy(
            args.stop_codon_genetic_code, args.stop_codon_sequence_code
        )
    except ValueError as error:
        raise ComparisonInputError(str(error)) from error


def _add_stop_codons(query, genome: str, policy) -> tuple:
    """Harmonise the query's stop-codon convention; returns (query, audit, summary)."""
    _log("Adding stop codons to the query CDS (--add-stop-codon)...")
    names = set(query.loc[query["Feature"] == "CDS", "Chromosome"])
    query, audit, summary = harmonise_stop_codons(
        query, iter_sequences(genome, names), policy
    )
    _log(
        f"  {summary['transcripts_with_cds']} query transcripts with CDS: "
        f"{summary['eligible']} with an established frame; {summary['extended']} "
        f"extended, {summary['unchanged']} unchanged, {summary['skipped']} skipped"
    )
    return query, audit, summary


def _split(value: str | None) -> list[str] | None:
    return (
        [item.strip() for item in value.split(",") if item.strip()] if value else None
    )


def _regions(args) -> list:
    regions = [parse_region(text) for text in args.region or []]
    if args.regions_file:
        regions.extend(load_regions_file(args.regions_file))
    return regions


def _analysis_scopes(args, ref, query, lengths, regions) -> dict[str, set]:
    """
    Named sequence subsets to report beside the whole comparison (scope "all").

    Every sequence must be verifiable: present in the genome FASTA when --genome
    is given (otherwise in one of the annotations), and inside the regions when
    the comparison is restricted. Scopes are applied to reference and query
    genes alike.
    """
    scopes: dict[str, list[str]] = {}
    for text in args.analysis_scope or []:
        name, _, sequences = text.partition("=")
        scopes.setdefault(name.strip(), []).extend(
            s.strip() for s in sequences.split(",") if s.strip()
        )
    if args.analysis_scopes_file:
        table = pd.read_csv(
            args.analysis_scopes_file, sep="\t", dtype=str, comment="#"
        ).fillna("")
        if not {"scope_id", "sequence"} <= set(table.columns):
            raise ComparisonInputError(
                "--analysis-scopes-file needs the columns scope_id and sequence"
            )
        for name, sequence in zip(table["scope_id"], table["sequence"]):
            scopes.setdefault(name.strip(), []).append(sequence.strip())
    known = (
        set(lengths) if lengths else set(ref["Chromosome"]) | set(query["Chromosome"])
    )
    region_names = {region.seqname for region in regions}
    out = {}
    for name, sequences in scopes.items():
        if not _SCOPE_ID.match(name) or name == contract.ALL_SCOPE:
            raise ComparisonInputError(
                f"analysis scope id {name!r} must be an identifier other than "
                f"{contract.ALL_SCOPE!r}"
            )
        if not sequences or len(set(sequences)) != len(sequences):
            raise ComparisonInputError(
                f"analysis scope {name}: list each sequence once, at least one"
            )
        unknown = sorted(set(sequences) - known)
        if unknown:
            where = "the genome FASTA" if lengths else "either annotation"
            raise ComparisonInputError(
                f"analysis scope {name}: sequences not in {where}: {unknown[:10]}"
            )
        outside = sorted(set(sequences) - region_names) if regions else []
        if outside:
            raise ComparisonInputError(
                f"analysis scope {name}: sequences outside --region/--regions-file: "
                f"{outside[:10]}"
            )
        out[name] = set(sequences)
    return out


def _provenance(args) -> dict | None:
    if not args.provenance_json:
        return None
    with open(args.provenance_json, encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ComparisonInputError("--provenance-json must hold a JSON object")
    return data


def _splice_available(args) -> bool:
    """Binned splice analyses need one reference transcript per gene."""
    return (
        args.reference_transcript_selection in ("longest_cds", "canonical")
        or args.evaluation_mode == "canonical"
    )


def _load_inputs(args, lengths):
    """Parse, map sequence names and (optionally) harmonise query stop codons."""
    _log("Loading annotations...")
    mapping = build_seqname_mapping(args.seqname_map, args.assembly_report)
    ref, ref_diag = _parse(args.reference, args.reference_format, "Reference")
    query, query_diag = _parse(args.query, args.query_format, "Query")
    parsed = {"reference": ref, "query": query}
    trackers = {role: ExclusionTracker(frame) for role, frame in parsed.items()}
    ref, query, seqname_checks = _harmonise_seqnames(ref, query, args, mapping, lengths)
    trackers["reference"].step(ref, "excluded_seqname")
    trackers["query"].step(query, "excluded_seqname")
    return {
        "mapping": mapping,
        "parsed": parsed,
        "trackers": trackers,
        "ref": ref,
        "query": query,
        "diagnostics": {"reference": ref_diag, "query": query_diag},
        "seqname_checks": seqname_checks,
    }


def run_pairwise_compare(  # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    args,
) -> dict:
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
        "analysis_scopes_file",
        "provenance_json",
    ):
        path = getattr(args, option)
        if path and not os.path.isfile(path):
            raise ComparisonInputError(
                f"--{option.replace('_', '-')}: file not found: {path}"
            )
    regions = _regions(args)
    stop_policy = _stop_codon_policy(args)
    provenance = _provenance(args)
    if args.plots_per_category:
        from ensembl.genes.annotation_qc.reports import (
            pairwise_plots,
        )  # optional matplotlib

        try:
            pairwise_plots.import_pyplot()
        except ImportError as error:
            raise ComparisonInputError(str(error)) from error

    lengths = parse_sequence_lengths(args.genome) if args.genome else None
    data = _load_inputs(args, lengths)
    ref, query = data["ref"], data["query"]
    trackers = data["trackers"]
    scopes = _analysis_scopes(args, ref, query, lengths, regions)
    stop: dict = {"audit": None, "transformed": {}}
    stop_summary = {"applied": False, "annotation": "query"}
    if stop_policy:
        query, audit_table, stop_summary = _add_stop_codons(
            query, args.genome, stop_policy
        )
        stop["audit"] = audit_table
        extended = audit_table[audit_table["status"] == "extended"]
        stop["transformed"] = dict.fromkeys(
            extended["transcript_id"], "stop_codon_added"
        )
    data["compared"] = {"reference": ref, "query": query}

    _log("Filtering reference...")
    # The unfiltered reference only labels Novel query genes (e.g. pseudogene loci).
    ref_context = ref
    ref, audit, ref_labels = _filter_reference(ref, args, trackers["reference"])
    if data["seqname_checks"]["excluded"].get("reference"):
        audit["excluded_seqnames"] = data["seqname_checks"]["excluded"]["reference"]
    query, query_labels = _select(
        query, args.query_transcript_selection, "query", trackers["query"]
    )
    if args.query_transcript_selection != "all":
        _log(
            f"  Query transcripts kept ({args.query_transcript_selection}): "
            f"{int((query['Feature'] == 'transcript').sum())}"
        )
    if regions:
        ref = subset_to_regions(ref, regions)
        query = subset_to_regions(query, regions)
        ref_context = subset_to_regions(ref_context, regions)
        trackers["reference"].step(ref, "outside_regions")
        trackers["query"].step(query, "outside_regions")
        _log(f"  Restricted to {len(regions)} region(s)")

    _log("Classifying loci...")
    ref_models, query_models = build_gene_models(ref), build_gene_models(query)
    trackers["reference"].finish(ref_models)
    trackers["query"].finish(query_models)
    ref_results, query_results = classify_loci(
        ref_models, query_models, reference_context=build_gene_models(ref_context)
    )
    summary = summarise_comparison(ref_results, query_results)
    summary["contract"] = {
        "id": contract.CONTRACT_ID,
        "version": contract.CONTRACT_VERSION,
        "definition_version": contract.DEFINITION_VERSION,
        "metrics_file": contract_report.METRICS_TSV,
    }
    summary["intron_support"] = summarise_intron_support(ref_models, query_models)
    summary["query_stop_codon_harmonisation"] = stop_summary
    summary["cds_exact_one_to_one"], exact_pairs = one_to_one_exact_cds(
        ref_models, query_models
    )
    contract.add_contract_columns(ref_results, query_results, exact_pairs, scopes)
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
    if scopes:
        summary["analysis_scopes"] = {sid: sorted(seqs) for sid, seqs in scopes.items()}
    summary["gene_span_checks"] = {
        role: {
            **summarise_span_excess(gene_span_excess(models)),
            "transcripts_outside": summarise_transcripts_outside(
                transcripts_outside_gene(models)
            ),
        }
        for role, models in (("reference", ref_models), ("query", query_models))
    }

    _log("Writing reports...")
    os.makedirs(args.outdir, exist_ok=True)
    report.write_filter_audit(audit, args.outdir)
    report.write_summary(summary, args.outdir)
    report.write_details(ref_results, query_results, args.outdir)
    report.write_transcript_labels(labels, args.outdir)
    report.write_split_merge(ref_results, query_results, args.outdir)
    report.write_exact_cds_pairs(exact_pairs, args.outdir)
    stop_outputs = None
    if stop["audit"] is not None:
        stop_outputs = {
            "audit": _file_record(
                report.write_stop_codon_audit(stop["audit"], args.outdir)
            ),
            "evaluated_query": _file_record(
                report.write_annotation_gtf(
                    query, args.outdir, report.HARMONISED_QUERY_GTF
                )
            ),
        }
    input_records = {
        "query": _file_record(args.query),
        "reference": _file_record(args.reference),
    }
    analyses = contract_report.write_contract_outputs(
        {
            "outdir": args.outdir,
            "export_evaluated_models": args.export_evaluated_models,
            "splice_available": _splice_available(args),
            "canonical_reference": "canonical"
            in (args.reference_transcript_selection, args.evaluation_mode),
        },
        data,
        {"reference": ref_models, "query": query_models},
        {"reference": ref_results, "query": query_results},
        exact_pairs,
        scopes,
        {
            "reference": {
                "mode": args.reference_transcript_selection,
                "labels": ref_labels,
            },
            "query": {"mode": args.query_transcript_selection, "labels": query_labels},
        },
        stop,
        input_records,
    )
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

    manifest = {
        "contract_id": contract.CONTRACT_ID,
        "contract_version": contract.CONTRACT_VERSION,
        "definition_version": contract.DEFINITION_VERSION,
        "command": "annotation-qc pairwise-compare",
        "argv": sys.argv,
        "created": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": round(time.time() - started, 1),
        "inputs": {
            **input_records,
            "genome": _file_record(args.genome),
            "seqname_map": _file_record(args.seqname_map),
            "assembly_report": _file_record(args.assembly_report),
            "evidence_attribution": _file_record(args.evidence_attribution),
            "regions_file": _file_record(args.regions_file),
            "analysis_scopes_file": _file_record(args.analysis_scopes_file),
            "provenance_json": _file_record(args.provenance_json),
        },
        "provenance": provenance,
        "gmb_finalise": _gmb_provenance(args.query, input_records["query"]["sha256"]),
        "options": {key: value for key, value in vars(args).items() if key != "func"},
        "policies": {
            "evaluation_mode": args.evaluation_mode,
            "reference_gene_biotypes": _split(args.reference_gene_biotypes),
            "reference_transcript_biotypes": _split(args.reference_transcript_biotypes),
            "reference_transcript_selection": args.reference_transcript_selection,
            "query_transcript_selection": args.query_transcript_selection,
            "genome_mismatch": args.genome_mismatch,
            "excluded_seqnames": data["seqname_checks"]["excluded"],
            "seqname_mapping_applied": {
                "reference": data["seqname_checks"]["reference_mapping"]["mapped"],
                "query": data["seqname_checks"]["query_mapping"]["mapped"],
            },
        },
        "transformations": (
            [
                {
                    "id": "stop_codon_added",
                    "side": "query",
                    "audit": report.STOP_AUDIT_TSV,
                    "summary": {k: v for k, v in stop_summary.items() if k != "policy"},
                    "policy": stop_summary.get("policy"),
                }
            ]
            if stop["audit"] is not None
            else []
        ),
        "scope": {
            "comparison_regions": [str(region) for region in regions] or None,
            "analysis_scopes": {sid: sorted(seqs) for sid, seqs in scopes.items()},
        },
        "thresholds": {
            "reciprocal_overlap": RECIPROCAL_THRESHOLD,
            "split_merge_min_shared_fraction": SPLIT_MERGE_MIN_SHARED_FRACTION,
        },
        "analyses": analyses,
        "metric_definitions": contract.metric_definitions(),
        "coordinates": "internally zero-based half-open; reported one-based inclusive",
        "parse_diagnostics": data["diagnostics"],
        "stop_codon_harmonisation": {**stop_summary, "outputs": stop_outputs},
        "seqname_checks": data["seqname_checks"],
        "evidence_attribution_unlabelled_rows": unlabelled,
        "code_version": _code_version(),
        "implementation": _implementation(),
    }
    manifest["outputs"] = contract_report.outputs_inventory(
        args.outdir, exclude=("comparison_manifest.json",)
    )
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

    stop = parser.add_argument_group("stop-codon harmonisation (optional; query only)")
    stop.add_argument(
        "--add-stop-codon",
        action="store_true",
        help="For query annotations whose CDS excludes the stop codon: extend each "
        "terminal CDS segment by 3 bp when the next three genomic bases are an "
        "in-frame stop codon and the frame is established from the sequence (ATG "
        "start, length a multiple of 3, no internal stop). Needs --genome. The "
        "reference is never changed. Writes stop_codon_harmonisation.tsv and "
        "query_evaluated.stop_harmonised.gtf. Default: compare original "
        "coordinates.",
    )
    stop.add_argument(
        "--stop-codon-genetic-code",
        type=int,
        default=1,
        metavar="TABLE",
        help="NCBI translation table for --add-stop-codon (default: 1, standard). "
        f"Supported: {', '.join(map(str, sorted(GENETIC_CODE_STOPS)))}. "
        "Organelle sequences (MT, chrM, Pt, ...) are skipped unless given a code "
        "with --stop-codon-sequence-code.",
    )
    stop.add_argument(
        "--stop-codon-sequence-code",
        action="append",
        default=None,
        metavar="SEQNAME=TABLE",
        help="Translation table for one sequence (comparison seqname), e.g. MT=2; "
        "TABLE none skips the sequence. Repeatable.",
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
        "--analysis-scope",
        action="append",
        default=None,
        metavar="ID=SEQ[,SEQ...]",
        help="Named sequence subset reported beside the whole comparison (scope "
        "'all') in metrics.tsv, the splice tables and scope_ids columns; applied "
        "to reference and query genes alike. Sequences must be in the genome "
        "FASTA (or an annotation without --genome). Repeatable.",
    )
    other.add_argument(
        "--analysis-scopes-file",
        default=None,
        help="TSV with columns scope_id and sequence (one row per sequence), "
        "same meaning as --analysis-scope.",
    )
    other.add_argument(
        "--provenance-json",
        default=None,
        help="JSON object copied verbatim into comparison_manifest.json "
        "('provenance'), e.g. the original source path of a renamed input, a "
        "predictor label or a run id. Never interpreted.",
    )
    other.add_argument(
        "--export-evaluated-models",
        choices=EXPORT_MODES,
        default="ids",
        help="evaluated_transcripts.tsv.gz content: ids (default) lists every "
        "transcript with its status and structure counts; full also writes "
        "one-based exon and CDS segments.",
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

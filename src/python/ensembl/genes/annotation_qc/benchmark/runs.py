"""
Run (or skip) pairwise comparisons with content-based cache keys.

A result is reused only when its cache key matches. The key covers:

    inputs          sha256 of the prepared query, prepared reference annotation,
                    prepared genome, scope regions content and seqname map
    preparation     the recipes that produced those prepared files (via their
                    output checksums and recorded source checksums)
    parameters      effective pairwise-compare options (config.comparator_parameters)
    implementation  content fingerprint of the comparator source files and the
                    python / pandas / pyranges1 versions

Comparisons run one at a time by default (a human whole-genome run needs up to
~7 GB RAM). Each attempt writes into ``comparisons/_partial/`` and is promoted to
``comparisons/<run_id>/<policy>/`` only after exit status 0 and all expected output
files exist; failures are appended to ``<policy>.attempts.jsonl`` and never
treated as results.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from ensembl.genes.annotation_qc.benchmark.config import (
    Experiment,
    Policy,
    Reference,
    Run,
    comparator_parameters,
)
from ensembl.genes.annotation_qc.benchmark.preparation import (
    prepare_reference,
    prepare_run,
)
from ensembl.genes.annotation_qc.benchmark.provenance import (
    current_implementation,
    environment,
    now,
    stable_hash,
)
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace

EXPECTED_OUTPUTS = (
    "comparison_summary.json",
    "comparison_details.tsv",
    "consensus_transcript_labels.tsv",
    "gene_splits.tsv",
    "gene_merges.tsv",
    "reference_filter_audit.json",
    "comparison_manifest.json",
)

STATUS_COMPLETE = "complete"
STATUS_FAILED = "failed"
STATUS_NOT_RUN = "not_run"
STATUS_STALE = "stale"
STATUS_UNAVAILABLE = "unavailable"


@dataclass
class Task:
    run: Run
    reference: Reference
    policy: Policy
    key: str
    key_components: dict
    argv_template: list[str]  # outdir appended at execution
    status: str
    reason: str


def key_components(inputs: dict, parameters: dict, implementation: dict) -> dict:
    return {
        "inputs": {
            k: inputs.get(k)
            for k in (
                "query_sha256",
                "reference_sha256",
                "genome_sha256",
                "scope_sha256",
                "seqname_map_sha256",
            )
        },
        "parameters": parameters,
        "implementation": {"fingerprint": implementation.get("fingerprint")},
    }


def cache_key(components: dict) -> str:
    return stable_hash(components)


def policy_available(policy: Policy, canonical_available: bool) -> tuple[bool, str]:
    if policy.reference_transcript_selection == "canonical" and not canonical_available:
        return False, (
            "reference has no Ensembl_canonical-tagged transcripts; canonical selection would silently "
            "fall back to longest CDS for every gene, so this policy is not run"
        )
    return True, ""


def canonical_available(ws: Workspace, prepared: dict) -> bool:
    from ensembl.genes.annotation_qc.benchmark.validation import cached_scan

    scan = cached_scan(ws, Path(prepared["annotation"]["output"]["path"]))
    return scan["canonical_tagged_records"] > 0


def plan(
    experiment: Experiment, ws: Workspace, run_ids=None, policy_ids=None, prepare=True
) -> list[Task]:
    """Work out, for every (run, policy), whether a valid result exists."""
    implementation = current_implementation()
    tasks = []
    prepared_refs: dict[str, dict] = {}
    for run in experiment.runs.values():
        if run_ids and run.run_id not in run_ids:
            continue
        reference = experiment.references[run.reference]
        if reference.id not in prepared_refs:
            prepared_refs[reference.id] = (
                prepare_reference(ws, reference)
                if prepare
                else ws.preparation_record(reference.id)
            )
        prep = prepared_refs[reference.id]
        if prep is None:
            raise RuntimeError(
                f"reference {reference.id} is not prepared; run `annotation-qc benchmark prepare`"
            )
        query = prepare_run(ws, run)
        inputs = {
            "query_sha256": query["output"]["sha256"],
            "query_path": query["output"]["path"],
            "reference_sha256": prep["annotation"]["output"]["sha256"],
            "reference_path": prep["annotation"]["output"]["path"],
            "genome_sha256": (
                prep["genome"]["output"]["sha256"] if prep.get("genome") else None
            ),
            "genome_path": (
                prep["genome"]["output"]["path"] if prep.get("genome") else None
            ),
            "scope_sha256": prep["scope"]["sha256"],
            "scope_path": prep["scope"]["regions_file"],
            "seqname_map_sha256": (prep.get("seqname_map") or {}).get("sha256"),
            "seqname_map_path": (prep.get("seqname_map") or {}).get("path"),
        }
        canonical = canonical_available(ws, prep)
        for policy in experiment.policies.values():
            if policy_ids and policy.id not in policy_ids:
                continue
            params = comparator_parameters(experiment, reference, run, policy)
            components = key_components(inputs, params, implementation)
            key = cache_key(components)
            argv = build_argv(inputs, params)
            ok, why = policy_available(policy, canonical)
            record = ws.run_record(run.run_id, policy.id)
            if not ok:
                status, reason = STATUS_UNAVAILABLE, why
            elif (
                record
                and record.get("status") == STATUS_COMPLETE
                and record.get("cache_key") == key
                and _outputs_present(record)
            ):
                status, reason = (
                    STATUS_COMPLETE,
                    "cached result matches inputs, parameters and implementation",
                )
            elif record and record.get("status") == STATUS_COMPLETE:
                status, reason = STATUS_STALE, describe_difference(
                    record.get("key_components") or {}, components
                )
            else:
                last = (ws.attempts(run.run_id, policy.id) or [None])[-1]
                if (
                    last
                    and last.get("status") == STATUS_FAILED
                    and last.get("cache_key") == key
                ):
                    status, reason = (
                        STATUS_FAILED,
                        f"last attempt failed (exit {last.get('exit_code')}); will retry",
                    )
                else:
                    status, reason = STATUS_NOT_RUN, "no result yet"
            tasks.append(
                Task(run, reference, policy, key, components, argv, status, reason)
            )
    return tasks


def describe_difference(old: dict, new: dict) -> str:
    changed = []
    for section in ("inputs", "parameters", "implementation"):
        a, b = old.get(section) or {}, new.get(section) or {}
        for k in sorted(set(a) | set(b)):
            if a.get(k) != b.get(k):
                changed.append(f"{section}.{k}")
    return "changed: " + (", ".join(changed) if changed else "cache key")


def _outputs_present(record: dict) -> bool:
    out = Path(record.get("output_dir", ""))
    return all((out / name).exists() for name in EXPECTED_OUTPUTS)


def build_argv(inputs: dict, params: dict) -> list[str]:
    argv = [
        sys.executable,
        "-m",
        "ensembl.genes.annotation_qc.runners.pairwise_compare",
        "--query",
        inputs["query_path"],
        "--reference",
        inputs["reference_path"],
    ]
    if inputs.get("genome_path"):
        argv += [
            "--genome",
            inputs["genome_path"],
            "--genome-mismatch",
            params["genome_mismatch"],
        ]
    argv += ["--evaluation-mode", params["evaluation_mode"]]
    if params["reference_gene_biotypes"]:
        argv += [
            "--reference-gene-biotypes",
            ",".join(params["reference_gene_biotypes"]),
        ]
    if params["reference_transcript_biotypes"]:
        argv += [
            "--reference-transcript-biotypes",
            ",".join(params["reference_transcript_biotypes"]),
        ]
    argv += [
        "--reference-transcript-selection",
        params["reference_transcript_selection"],
        "--query-transcript-selection",
        params["query_transcript_selection"],
        "--query-format",
        params["query_format"],
        "--reference-format",
        params["reference_format"],
    ]
    if inputs.get("seqname_map_path"):
        argv += ["--seqname-map", inputs["seqname_map_path"]]
    if inputs.get("scope_path"):
        argv += ["--regions-file", inputs["scope_path"]]
    if params["plots_per_category"]:
        argv += ["--plots-per-category", str(params["plots_per_category"])]
    return argv


def execute(ws: Workspace, task: Task, log=print) -> dict:
    """Run one comparison; promote output only on success."""
    final = ws.comparison_dir(task.run.run_id, task.policy.id)
    partial = (
        ws.partial_root() / f"{task.run.run_id}__{task.policy.id}__{int(time.time())}"
    )
    partial.mkdir(parents=True, exist_ok=True)
    argv = task.argv_template + ["--outdir", str(partial / "out")]
    started = now()
    t0 = time.time()
    with (
        open(partial / "stdout.log", "w") as out,
        open(partial / "stderr.log", "w") as err,
    ):
        proc = subprocess.Popen(argv, stdout=out, stderr=err)
        _, wait_status, usage = os.wait4(proc.pid, 0)  # per-child resource usage
        proc.returncode = os.waitstatus_to_exitcode(wait_status)
    peak = (
        usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024
    )  # bytes on macOS, KiB on Linux
    record = {
        "run_id": task.run.run_id,
        "policy_id": task.policy.id,
        "reference_id": task.reference.id,
        "cache_key": task.key,
        "key_components": task.key_components,
        "argv": argv,
        "started": started,
        "finished": now(),
        "wall_seconds": round(time.time() - t0, 1),
        "peak_rss_bytes": peak,
        "exit_code": proc.returncode,
        "environment": environment(),
        "implementation": current_implementation(),
        "source": "computed",
    }
    missing = [n for n in EXPECTED_OUTPUTS if not (partial / "out" / n).exists()]
    if proc.returncode != 0 or missing:
        record.update(
            status=STATUS_FAILED,
            partial_dir=str(partial),
            missing_outputs=missing,
            stderr_tail=(partial / "stderr.log").read_text()[-4000:],
        )
        ws.append_attempt(task.run.run_id, task.policy.id, record)
        log(
            f"  FAILED {task.run.run_id}/{task.policy.id} (exit {proc.returncode}); partial output kept in {partial}"
        )
        return record
    if final.exists():
        shutil.rmtree(final)
    final.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(partial / "out"), str(final))
    for name in ("stdout.log", "stderr.log"):
        shutil.move(str(partial / name), str(final / name))
    shutil.rmtree(partial, ignore_errors=True)
    record.update(status=STATUS_COMPLETE, output_dir=str(final))
    ws.write_run_record(task.run.run_id, task.policy.id, record)
    ws.append_attempt(
        task.run.run_id,
        task.policy.id,
        {k: v for k, v in record.items() if k not in ("environment", "implementation")},
    )
    log(f"  done {task.run.run_id}/{task.policy.id} in {record['wall_seconds']} s")
    return record


def run_missing(
    experiment: Experiment,
    ws: Workspace,
    run_ids=None,
    policy_ids=None,
    dry_run=False,
    force=False,
    log=print,
) -> list[Task]:
    tasks = plan(experiment, ws, run_ids, policy_ids)
    todo = [
        t
        for t in tasks
        if t.status in (STATUS_NOT_RUN, STATUS_STALE, STATUS_FAILED)
        or (force and t.status == STATUS_COMPLETE)
    ]
    for t in tasks:
        marker = "→ run" if t in todo else "  skip"
        log(f"{marker} {t.run.run_id:<28} {t.policy.id:<3} {t.status:<11} {t.reason}")
    if dry_run:
        return tasks
    for i, task in enumerate(todo, 1):
        log(f"[{i}/{len(todo)}] {task.run.run_id} policy {task.policy.id}")
        execute(ws, task, log)
    return tasks

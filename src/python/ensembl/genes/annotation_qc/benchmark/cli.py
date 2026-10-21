"""
``annotation-qc benchmark <command> --config EXPERIMENT.json``

    validate         check configuration and inputs (writes validation/validation.json)
    prepare          run recorded input preparation (decompress, retype) and derive scope
    import           register completed comparator results (results and pilot directories)
    run              run missing, stale, failed or unverified comparisons, sequentially
    independent      run the independent CDS-signature analysis (optional cross-check)
    build            build or update the dashboard dataset (SQLite)
    status           show every (run, policy) result status without running anything
    provenance       show the implementation identity of every result (registry check)
    register-build   write a registry entry for the build behind some manifests (git-checked)
    verify-headline  compare legacy dataset metrics with a headline_metrics.tsv (wide) table
"""

from __future__ import annotations

import json
import sqlite3
import sys

import pandas as pd

from ensembl.genes.annotation_qc.benchmark.config import ConfigError, load_experiment
from ensembl.genes.annotation_qc.benchmark.workspace import Workspace


def _experiment(args):
    try:
        experiment = load_experiment(args.config)
    except ConfigError as error:
        sys.exit(str(error))
    return experiment, Workspace(experiment)


def _read_only_plan(experiment, ws):
    """The plan without preparing anything: status and provenance never write records."""
    from ensembl.genes.annotation_qc.benchmark.preparation import PreparationError
    from ensembl.genes.annotation_qc.benchmark.runs import plan

    try:
        return plan(experiment, ws, prepare=False)
    except (PreparationError, RuntimeError) as error:
        sys.exit(str(error))


def cmd_validate(args):
    from ensembl.genes.annotation_qc.benchmark.preparation import (
        PreparationError,
        prepare_reference,
    )
    from ensembl.genes.annotation_qc.benchmark.validation import validate_experiment

    experiment, ws = _experiment(args)
    prepared = {}
    for ref_id in experiment.references:
        rec = ws.preparation_record(ref_id)
        if args.prepare:
            try:
                rec = prepare_reference(ws, experiment.references[ref_id])
            except PreparationError as error:
                print(f"ERROR preparing {ref_id}: {error}")
                rec = None
        if rec:
            prepared[ref_id] = rec
    report = validate_experiment(experiment, ws, prepared)
    for issue in report["issues"]:
        print(f"{issue['level'].upper():8} {issue['where']}: {issue['message']}")
    print("\nKnown comparator parser limitations (probed against this checkout):")
    for lim in report["parser_limitations"]:
        print(f"  {lim['title']}: {lim['status']} — {lim['handling']}")
    print(
        f"\nValidation {'passed' if report['ok'] else 'FAILED'}; report: {ws.validation_file()}"
    )
    sys.exit(0 if report["ok"] else 1)


def cmd_prepare(args):
    from ensembl.genes.annotation_qc.benchmark.preparation import (
        prepare_reference,
        prepare_run,
    )

    experiment, ws = _experiment(args)
    for ref in experiment.references.values():
        rec = prepare_reference(ws, ref, force=args.force)
        print(
            f"{ref.id}: annotation {rec['annotation']['method']} -> {rec['annotation']['output']['path']}"
        )
        if rec.get("genome"):
            print(
                f"{ref.id}: genome {rec['genome']['method']} -> {rec['genome']['output']['path']}"
            )
        print(
            f"{ref.id}: scope {rec['scope']['type']} ({len(rec['scope']['sequences'] or []) or 'all'} sequences)"
        )
    for run in experiment.runs.values():
        rec = prepare_run(ws, run, force=args.force)
        print(f"{run.run_id}: {rec['method']}")


def cmd_import(args):
    from ensembl.genes.annotation_qc.benchmark.importer import import_benchmark

    experiment, ws = _experiment(args)
    import_benchmark(experiment, ws, args.benchmark_dir)


def cmd_run(args):
    from ensembl.genes.annotation_qc.benchmark.runs import run_missing

    experiment, ws = _experiment(args)
    run_missing(
        experiment,
        ws,
        run_ids=args.runs,
        policy_ids=args.policies,
        dry_run=args.dry_run,
        force=args.force,
    )


def cmd_status(args):
    experiment, ws = _experiment(args)
    for task in _read_only_plan(experiment, ws):
        print(
            f"{task.reference.id:<24} {task.run.run_id:<28} {task.policy.id:<3} {task.status:<11} {task.reason}"
        )


def cmd_provenance(args):
    """Implementation identity of each result: (scheme, fingerprint) against the registries."""
    experiment, ws = _experiment(args)
    counts: dict = {}
    for task in _read_only_plan(experiment, ws):
        identity = (task.assessment or {}).get("identity") or {}
        status = identity.get("status", "no result")
        counts[status] = counts.get(status, 0) + 1
        print(
            f"{task.run.run_id:<28} {task.policy.id:<3} {task.status:<11} {status:<17} "
            f"{identity.get('scheme') or '-':<40} {(identity.get('fingerprint') or '-')[:16]:<16} "
            f"{identity.get('build_id') or ''} {identity.get('reason') or ''}"
        )
    print("; ".join(f"{v} {k}" for k, v in sorted(counts.items())))


def cmd_register_build(args):
    from pathlib import Path

    from ensembl.genes.annotation_qc.benchmark.registry import (
        RegistryError,
        register_build,
    )

    evidence = json.loads(args.evidence) if args.evidence else []
    try:
        entry = register_build(
            [Path(m) for m in args.manifest],
            Path(args.repo),
            args.commit,
            args.build_id,
            evidence,
            args.notes,
        )
    except RegistryError as error:
        sys.exit(f"not registered: {error}")
    out = Path(args.output)
    data = (
        json.loads(out.read_text())
        if out.exists()
        else {"registry_version": 1, "builds": []}
    )
    data["builds"] = [
        b
        for b in data["builds"]
        if (b["scheme"], b["fingerprint"]) != (entry["scheme"], entry["fingerprint"])
    ] + [entry]
    out.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")
    print(
        f"registered {entry['build_id']} ({entry['scheme']}, {entry['fingerprint'][:16]}…, "
        f"{len(entry['source_files'])} files) in {out}"
    )


def cmd_independent(args):
    from ensembl.genes.annotation_qc.benchmark.independent.cds_signatures import (
        run_independent,
    )
    from ensembl.genes.annotation_qc.benchmark.preparation import (
        prepare_reference,
        prepare_run,
    )
    from ensembl.genes.annotation_qc.benchmark.runs import (
        canonical_available,
        policy_available,
    )

    experiment, ws = _experiment(args)
    for ref in experiment.references.values():
        runs = [
            r
            for r in experiment.runs_for(ref.id)
            if not args.runs or r.run_id in args.runs
        ]
        if not runs:
            continue
        prep = prepare_reference(ws, ref)
        seqs = prep["scope"]["sequences"]
        canonical = canonical_available(ws, prep)
        policies = {
            p.id: (p.reference_transcript_selection, p.query_transcript_selection)
            for p in experiment.policies.values()
            if (not args.policies or p.id in args.policies)
            and policy_available(p, canonical)[0]
        }
        print(
            f"Independent analysis for {ref.id} ({len(runs)} runs, policies {list(policies)})"
        )
        run_independent(
            prep["annotation"]["output"]["path"],
            set(seqs) if seqs else None,
            experiment.evaluation,
            {r.run_id: prepare_run(ws, r)["output"]["path"] for r in runs},
            policies,
            ws.root / "independent",
        )


def cmd_build(args):
    from ensembl.genes.annotation_qc.benchmark.preparation import PreparationError
    from ensembl.genes.annotation_qc.benchmark.store import build_dashboard

    experiment, ws = _experiment(args)
    try:
        build_dashboard(experiment, ws)
    except PreparationError as error:
        sys.exit(str(error))


def cmd_verify_headline(args):
    """Every mapped column of every row of a headline table must equal the dashboard metric."""
    from ensembl.genes.annotation_qc.benchmark.legacy import (
        HEADLINE_COLUMNS,
        compare_headline_value,
    )

    experiment, ws = _experiment(args)
    mapping = json.loads(args.mapping) if args.mapping else {}
    headline = pd.read_csv(args.headline, sep="\t")
    con = sqlite3.connect(f"file:{ws.dashboard_db()}?mode=ro", uri=True)
    metrics = {
        (r, p, m): (n, d, v, s)
        for r, p, m, n, d, v, s in con.execute(
            "SELECT run_id, policy_id, metric_id, numerator, denominator, value, status FROM metrics "
            "WHERE scope_id='all'"
        )
    }
    runs = set(r for (r, _, _) in metrics)
    problems, checked, rows, unavailable = [], 0, 0, {}
    aliases = {
        alias: run.run_id for run in experiment.runs.values() for alias in run.aliases
    }
    for row in headline.to_dict("records"):
        label = f"{row['species']}_{row['tool']}"
        run_id = mapping.get("runs", {}).get(label) or aliases.get(label, label)
        policy_id = mapping.get("views", {}).get(row["view"], row["view"].split("_")[0])
        if run_id not in runs:
            problems.append(
                f"row {row['species']}/{row['view']}/{row['tool']}: run {run_id} has no metrics"
            )
            continue
        rows += 1
        for column, (metric_id, field) in HEADLINE_COLUMNS.items():
            if column not in row:
                continue
            got = metrics.get((run_id, policy_id, metric_id))
            if got is None:
                problems.append(f"{run_id}/{policy_id}: metric {metric_id} missing")
                continue
            if got[3] == "not_available":
                # not in the result's summary (e.g. an earlier dashboard derivation); never compared as 0
                unavailable[column] = unavailable.get(column, 0) + 1
                continue
            value = {"numerator": got[0], "denominator": got[1], "value": got[2]}[field]
            expected = row[column]
            checked += 1
            if not compare_headline_value(
                column, None if pd.isna(expected) else expected, value
            ):
                problems.append(
                    f"{run_id}/{policy_id} {column}: headline {expected} vs dashboard {value} ({got[3]})"
                )
    unmapped = [
        c
        for c in headline.columns
        if c not in HEADLINE_COLUMNS and c not in ("species", "view", "tool")
    ]
    print(
        f"Checked {checked} values in {rows} of {len(headline)} headline rows; {len(problems)} problems"
    )
    if unmapped:
        print(f"Columns without a dashboard metric: {unmapped}")
    if unavailable:
        print(
            "Mapped columns whose metric is not available for these results (not compared): "
            + ", ".join(f"{c} ({n} rows)" for c, n in sorted(unavailable.items()))
        )
    for problem in problems[:50]:
        print("  " + problem)
    sys.exit(0 if not problems and rows == len(headline) else 1)


def register(subparsers):
    parser = subparsers.add_parser(
        "benchmark",
        help="Reproducible annotation-comparison experiments (config, import, run, build).",
        description=__doc__,
        formatter_class=__import__("argparse").RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="benchmark_command", required=True)

    def add(name, func, help_text):
        p = sub.add_parser(name, help=help_text)
        p.add_argument(
            "--config",
            required=True,
            help="Experiment configuration (JSON, or YAML with PyYAML).",
        )
        p.set_defaults(func=func)
        return p

    p = add("validate", cmd_validate, "Validate configuration and inputs.")
    p.add_argument(
        "--prepare",
        action="store_true",
        help="Run preparation first and validate prepared files.",
    )
    p = add(
        "prepare",
        cmd_prepare,
        "Prepare inputs (recorded decompression/retyping, scope).",
    )
    p.add_argument("--force", action="store_true")
    p = add(
        "import", cmd_import, "Import completed comparisons from a benchmark directory."
    )
    p.add_argument(
        "--benchmark-dir",
        default=None,
        help="Defaults to import.benchmark_dir in the config.",
    )
    for name, func, text in (
        ("run", cmd_run, "Run missing/stale/failed comparisons sequentially."),
        ("independent", cmd_independent, "Run the independent cross-check analysis."),
    ):
        p = add(name, func, text)
        p.add_argument(
            "--runs", nargs="*", default=None, help="Restrict to these run_ids."
        )
        p.add_argument(
            "--policies", nargs="*", default=None, help="Restrict to these policy ids."
        )
        if name == "run":
            p.add_argument(
                "--dry-run", action="store_true", help="Show what would run."
            )
            p.add_argument(
                "--force", action="store_true", help="Rerun complete results too."
            )
    add("build", cmd_build, "Build or update the dashboard dataset.")
    add("status", cmd_status, "Show result status per run and policy.")
    add(
        "provenance",
        cmd_provenance,
        "Show the implementation identity (registry check) of every result.",
    )
    p = sub.add_parser(
        "register-build",
        help="Write a registry entry for the build behind manifests (file hashes checked against git).",
    )
    p.add_argument("--manifest", nargs="+", required=True)
    p.add_argument("--repo", required=True, help="Git checkout containing --commit.")
    p.add_argument("--commit", required=True)
    p.add_argument("--build-id", required=True)
    p.add_argument("--output", required=True, help="Registry JSON to create or update.")
    p.add_argument(
        "--evidence", default=None, help="JSON list of extra evidence records."
    )
    p.add_argument("--notes", default=None)
    p.set_defaults(func=cmd_register_build)
    p = add(
        "verify-headline",
        cmd_verify_headline,
        "Compare dashboard metrics with a headline_metrics.tsv.",
    )
    p.add_argument("--headline", required=True)
    p.add_argument(
        "--mapping",
        default=None,
        help='JSON {"runs": {"<species>_<tool>": run_id}, "views": {view: policy}}',
    )

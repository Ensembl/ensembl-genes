"""
Import completed pairwise-compare outputs without rerunning or copying them.

Every ``comparison_manifest.json`` below the configured directories is matched to
a configured (run, policy) by content, not by folder name:

    query      manifest query sha256 == prepared query sha256 of a configured run
    reference  manifest reference sha256 == prepared reference annotation sha256
    settings   manifest stop-codon options == the run's comparator_options, so an
               original-coordinate and a stop-harmonised result of the same file
               go to different runs (a manifest without the option predates it
               and compared original coordinates)
    policy     transcript selections in the manifest == the policy's

Directories:

    import.results_dirs   headline results (and import.benchmark_dir, the version-1
                          key, or --benchmark-dir)
    import.pilot_dirs     always pilots: listed in the import report, never shown
                          as headline results, whatever their scope

A manifest restricted with ``--region`` (or with a regions file that is not the
reference scope) is a pilot too. The status of each matched result comes from
benchmark/results.py: output checksums against the manifest, inputs, parameters
and definition version against the configuration, and the implementation
against the registries. Imported records never overwrite a computed result.

Independent-analysis tables are imported through ``import.independent_tables``, a
path pattern with ``{alias}`` and ``{policy}`` placeholders resolved with each run's
``aliases``; they are verification only.
"""

from __future__ import annotations

from pathlib import Path

from ensembl.genes.annotation_qc.benchmark.config import (
    Experiment,
    comparator_parameters,
)
from ensembl.genes.annotation_qc.benchmark.preparation import (
    prepare_reference,
    prepare_run,
)
from ensembl.genes.annotation_qc.benchmark.provenance import (
    implementation_from_snapshot,
    now,
)
from ensembl.genes.annotation_qc.benchmark import registry as reg
from ensembl.genes.annotation_qc.benchmark.results import (
    STATUS_COMPLETE,
    STATUS_LEGACY,
    STATUS_UNVERIFIED,
    assess,
    current_definition,
    key_components,
    manifest_inputs,
    manifest_parameters,
    manifest_stop_codon_options,
)
from ensembl.genes.annotation_qc.benchmark.runs import expected_inputs, load_registry
from ensembl.genes.annotation_qc.benchmark.workspace import (
    Workspace,
    read_json,
    write_json,
)

PREFERENCE = {STATUS_COMPLETE: 3, STATUS_UNVERIFIED: 2, STATUS_LEGACY: 1}


def find_manifests(root: Path, max_depth: int = 6) -> list[Path]:
    found = []
    root = Path(root)
    for path in sorted(root.rglob("comparison_manifest.json")):
        if len(path.relative_to(root).parts) <= max_depth + 1:
            found.append(path)
    return found


def _directories(experiment: Experiment, benchmark_dir: Path | None) -> list[tuple]:
    results = [Path(benchmark_dir)] if benchmark_dir else []
    if not results:
        if experiment.imports.get("benchmark_dir"):
            results.append(Path(experiment.imports["benchmark_dir"]))
        results += [Path(p) for p in experiment.imports.get("results_dirs") or []]
    pilots = [Path(p) for p in experiment.imports.get("pilot_dirs") or []]
    for path in results + pilots:
        if not path.is_dir():
            raise FileNotFoundError(f"result directory not found: {path}")
    if not results and not pilots:
        raise FileNotFoundError(
            "nothing to import: set import.results_dirs (or import.benchmark_dir) or pass --benchmark-dir"
        )
    return [(p, False) for p in results] + [(p, True) for p in pilots]


def _policy_for(experiment, ref, run, params):
    for pol in experiment.policies.values():
        expected = comparator_parameters(experiment, ref, run, pol)
        if all(
            expected[k] == params.get(k)
            for k in ("reference_transcript_selection", "query_transcript_selection")
        ):
            return pol
    return None


def _match_run(query_index, prepared, inputs, manifest):
    """The single configured run a manifest belongs to, or (None, reason)."""
    runs = [
        r
        for r in query_index.get(inputs["query_sha256"], [])
        if prepared[r.reference]["annotation"]["output"]["sha256"]
        == inputs["reference_sha256"]
    ]
    if not runs:
        return None, "query/reference checksums match no configured run"
    stop = manifest_stop_codon_options(manifest)
    same = [r for r in runs if r.comparator_options == stop]
    if len(same) == 1:
        return same[0], None
    if not same:
        return None, (
            f"stop-codon settings {stop or 'original coordinates'} match no configured run "
            f"with this query (runs: {', '.join(r.run_id + ' ' + str(r.comparator_options or 'original coordinates') for r in runs)})"
        )
    return None, (
        "several configured runs share this query, reference and settings: "
        f"{', '.join(r.run_id for r in same)}"
    )


def import_benchmark(  # pylint: disable=too-many-locals,too-many-statements
    experiment: Experiment, ws: Workspace, benchmark_dir: Path | None = None, log=print
) -> dict:
    directories = _directories(experiment, benchmark_dir)
    snapshot = experiment.imports.get("code_snapshot") or []
    registry = load_registry(experiment)
    definition = current_definition()

    prepared = {
        ref.id: prepare_reference(ws, ref) for ref in experiment.references.values()
    }
    query_index: dict[str, list] = {}
    for run in experiment.runs.values():
        rec = prepare_run(ws, run)
        query_index.setdefault(rec["output"]["sha256"], []).append(run)

    report = {
        "directories": [
            {"path": str(p), "kind": "pilots" if pilot else "results"}
            for p, pilot in directories
        ],
        "benchmark_dir": str(directories[0][0]),
        "imported": now(),
        "matched": [],
        "pilots": [],
        "unmatched": [],
        "conflicts": [],
        "code_snapshot": [str(p) for p in snapshot],
        "registries": registry.sources,
    }
    chosen: dict[tuple, dict] = {}
    for root, pilot_dir in directories:
        for manifest_path in find_manifests(root):
            manifest = read_json(manifest_path)
            out_dir = manifest_path.parent
            rel = str(out_dir.relative_to(root))
            if manifest is None:
                report["unmatched"].append(
                    {"dir": rel, "reason": "unreadable manifest"}
                )
                continue
            inputs = manifest_inputs(manifest)
            params = manifest_parameters(manifest)
            run, why = _match_run(query_index, prepared, inputs, manifest)
            if run is None:
                report["unmatched"].append(
                    {"dir": rel, "output_dir": str(out_dir), "reason": why}
                )
                continue
            ref = experiment.references[run.reference]
            prep = prepared[ref.id]
            policy = _policy_for(experiment, ref, run, params)
            entry = {
                "dir": rel,
                "output_dir": str(out_dir),
                "run_id": run.run_id,
                "reference_id": ref.id,
                "policy_id": policy.id if policy else None,
                "created": manifest.get("created"),
                "runtime_seconds": manifest.get("runtime_seconds"),
                "options_region": params["region"],
                "contract_version": manifest.get("contract_version"),
                "implementation": (
                    reg.label(reg.identify(manifest, registry))
                    if manifest.get("contract_version")
                    else "pre-contract"
                ),
                "provenance": manifest.get("provenance"),
            }
            if (
                pilot_dir
                or params["region"]
                or inputs["scope_sha256"] != prep["scope"]["sha256"]
            ):
                entry["reason"] = (
                    "configured pilot directory (never a headline result)"
                    if pilot_dir
                    else "region-restricted or scope differs from the reference scope "
                    "(pilot, not whole-genome)"
                )
                report["pilots"].append(entry)
                continue
            if policy is None:
                entry["reason"] = (
                    f"transcript selections {params['reference_transcript_selection']}/"
                    f"{params['query_transcript_selection']} match no policy"
                )
                report["unmatched"].append(entry)
                continue
            expected = key_components(
                expected_inputs(ws, run, prep),
                comparator_parameters(experiment, ref, run, policy),
                definition,
            )
            snapshot_identity = None
            if snapshot and not manifest.get("contract_version"):
                versions = {
                    k: (manifest.get("code_version") or {}).get(k)
                    for k in ("python", "pandas", "pyranges1")
                }
                snapshot_identity = implementation_from_snapshot(snapshot, versions)
            assessment = assess(
                out_dir, expected, registry, ws.checksums, snapshot_identity
            )
            status = assessment["status"]
            record = {
                "status": status,
                "source": "imported",
                "run_id": run.run_id,
                "policy_id": policy.id,
                "reference_id": ref.id,
                "output_dir": str(out_dir),
                "key_components": assessment.get("components"),
                "difference": (
                    None if status == STATUS_COMPLETE else assessment["reason"]
                ),
                "argv": manifest.get("argv"),
                "started": None,
                "finished": manifest.get("created"),
                "wall_seconds": manifest.get("runtime_seconds"),
                "peak_rss_bytes": None,
                "exit_code": 0,
                "manifest_code_version": manifest.get("code_version"),
                "identity": assessment.get("identity"),
                "snapshot_identity": snapshot_identity,
                "imported": now(),
            }
            k = (run.run_id, policy.id)
            if k in chosen:
                prev = chosen[k]
                report["conflicts"].append(
                    {
                        "run_id": run.run_id,
                        "policy_id": policy.id,
                        "dirs": [prev["output_dir"], str(out_dir)],
                    }
                )
                if (PREFERENCE.get(status, 0), record["finished"] or "") <= (
                    PREFERENCE.get(prev["status"], 0),
                    prev["finished"] or "",
                ):
                    continue
                report["matched"] = [
                    m
                    for m in report["matched"]
                    if not (
                        m["run_id"] == run.run_id and m.get("policy_id") == policy.id
                    )
                ]
            chosen[k] = record
            entry.update(status=status, difference=record["difference"])
            report["matched"].append(entry)

    for (run_id, policy_id), record in chosen.items():
        existing = ws.run_record(run_id, policy_id)
        if existing and existing.get("source") == "computed":
            log(
                f"  keep computed result for {run_id}/{policy_id} (not overwritten by import)"
            )
            continue
        ws.write_run_record(run_id, policy_id, record)
        ws.append_attempt(
            run_id,
            policy_id,
            {
                k: v
                for k, v in record.items()
                if k not in ("identity", "key_components", "snapshot_identity")
            },
        )

    _import_independent(experiment, ws, directories, report)
    write_json(ws.import_report(), report)
    by_status: dict[str, int] = {}
    for record in chosen.values():
        by_status[record["status"]] = by_status.get(record["status"], 0) + 1
    log(
        f"Imported {len(chosen)} results "
        f"({', '.join(f'{v} {k}' for k, v in sorted(by_status.items())) or 'none'}), "
        f"{len(report['pilots'])} pilots kept separate, {len(report['unmatched'])} unmatched, "
        f"{len(report['conflicts'])} conflicts"
    )
    for record in chosen.values():
        if record["status"] != STATUS_COMPLETE:
            log(
                f"  {record['status'].upper()} {record['run_id']}/{record['policy_id']}: "
                f"{record['difference']}"
            )
    return report


def _import_independent(experiment, ws, directories, report) -> None:
    pattern = experiment.imports.get("independent_tables")
    report["independent"] = []
    if not pattern:
        return
    root = directories[0][0]
    for run in experiment.runs.values():
        for policy in experiment.policies.values():
            for alias in run.aliases:
                path = root / pattern.format(
                    alias=alias, policy=policy.id, run_id=run.run_id
                )
                if path.exists():
                    write_json(
                        ws.independent_dir(run.run_id) / f"{policy.id}.import.json",
                        {
                            "source_path": str(path),
                            "sha256": ws.checksums.sha256(path),
                            "alias": alias,
                            "imported": now(),
                            "verified": None,
                        },
                    )
                    report["independent"].append(
                        {
                            "run_id": run.run_id,
                            "policy_id": policy.id,
                            "path": str(path),
                        }
                    )
                    break

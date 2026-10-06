"""
Import completed pairwise-compare outputs from an existing benchmark directory.

Nothing is rerun or copied. Every ``comparison_manifest.json`` below the benchmark
directory is matched to a configured (run, policy) by content, not by folder name:

    query      manifest query sha256 == prepared query sha256 of a configured run
    reference  manifest reference sha256 == prepared reference annotation sha256
    genome / scope / seqname map   compared through the cache key
    policy     effective comparator options in the manifest == policy parameters

Manifests restricted with ``--region`` (or whose regions file is not the
reference scope) are recorded as *pilots* and kept apart from whole-genome results.

The imported cache key is rebuilt from the manifest (input checksums, options) and
the implementation fingerprint of a supplied code snapshot. If it equals the key
the current configuration and checkout would produce, the result is
``complete`` and reused; otherwise it is recorded as ``stale`` with the difference.

Independent-analysis tables are imported through ``import.independent_tables``, a
path pattern with ``{alias}`` and ``{policy}`` placeholders resolved with each run's
``aliases``; they are verified against the comparator per-gene output at build time.
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
    current_implementation,
    implementation_from_snapshot,
    now,
)
from ensembl.genes.annotation_qc.benchmark.runs import (
    EXPECTED_OUTPUTS,
    STATUS_COMPLETE,
    STATUS_STALE,
    cache_key,
    describe_difference,
    key_components,
)
from ensembl.genes.annotation_qc.benchmark.workspace import (
    Workspace,
    read_json,
    write_json,
)


def _biotypes(value) -> list | None:
    if not value:
        return None
    items = value if isinstance(value, list) else str(value).split(",")
    return sorted(item.strip() for item in items if item.strip()) or None


def manifest_parameters(manifest: dict) -> dict:
    opts = manifest.get("options", {})
    return {
        "evaluation_mode": opts.get("evaluation_mode"),
        "reference_gene_biotypes": _biotypes(opts.get("reference_gene_biotypes")),
        "reference_transcript_biotypes": _biotypes(
            opts.get("reference_transcript_biotypes")
        ),
        "reference_transcript_selection": opts.get("reference_transcript_selection"),
        "query_transcript_selection": opts.get("query_transcript_selection"),
        "genome_mismatch": opts.get("genome_mismatch"),
        "query_format": opts.get("query_format"),
        "reference_format": opts.get("reference_format"),
        "plots_per_category": opts.get("plots_per_category") or 0,
        "region": opts.get("region") or None,
        "assembly_report": opts.get("assembly_report"),
    }


def manifest_inputs(manifest: dict) -> dict:
    inputs = manifest.get("inputs", {})
    sha = lambda key: (inputs.get(key) or {}).get("sha256")
    return {
        "query_sha256": sha("query"),
        "reference_sha256": sha("reference"),
        "genome_sha256": sha("genome"),
        "scope_sha256": sha("regions_file")
        or ("all" if not manifest.get("options", {}).get("region") else None),
        "seqname_map_sha256": sha("seqname_map"),
    }


def find_manifests(root: Path, max_depth: int = 4) -> list[Path]:
    found = []
    root = Path(root)
    for path in sorted(root.rglob("comparison_manifest.json")):
        if len(path.relative_to(root).parts) <= max_depth + 1:
            found.append(path)
    return found


def import_benchmark(
    experiment: Experiment, ws: Workspace, benchmark_dir: Path | None = None, log=print
) -> dict:
    benchmark_dir = Path(benchmark_dir or experiment.imports.get("benchmark_dir") or "")
    if not benchmark_dir.is_dir():
        raise FileNotFoundError(f"benchmark directory not found: {benchmark_dir}")
    snapshot = experiment.imports.get("code_snapshot") or []
    current = current_implementation()

    prepared = {
        ref.id: prepare_reference(ws, ref) for ref in experiment.references.values()
    }
    query_index: dict[str, list] = {}
    for run in experiment.runs.values():
        rec = prepare_run(ws, run)
        query_index.setdefault(rec["output"]["sha256"], []).append(run)

    report = {
        "benchmark_dir": str(benchmark_dir),
        "imported": now(),
        "matched": [],
        "pilots": [],
        "unmatched": [],
        "conflicts": [],
        "code_snapshot": [str(p) for p in snapshot],
    }
    chosen: dict[tuple, dict] = {}
    for manifest_path in find_manifests(benchmark_dir):
        manifest = read_json(manifest_path)
        out_dir = manifest_path.parent
        rel = str(out_dir.relative_to(benchmark_dir))
        if manifest is None:
            report["unmatched"].append({"dir": rel, "reason": "unreadable manifest"})
            continue
        inputs = manifest_inputs(manifest)
        params = manifest_parameters(manifest)
        runs = query_index.get(inputs["query_sha256"], [])
        runs = [
            r
            for r in runs
            if prepared[r.reference]["annotation"]["output"]["sha256"]
            == inputs["reference_sha256"]
        ]
        if not runs:
            report["unmatched"].append(
                {
                    "dir": rel,
                    "reason": "query/reference checksums match no configured run",
                }
            )
            continue
        run = runs[0]
        ref = experiment.references[run.reference]
        prep = prepared[ref.id]
        is_pilot = bool(params["region"]) or (
            inputs["scope_sha256"] != prep["scope"]["sha256"]
        )
        entry = {
            "dir": rel,
            "output_dir": str(out_dir),
            "run_id": run.run_id,
            "reference_id": ref.id,
            "created": manifest.get("created"),
            "runtime_seconds": manifest.get("runtime_seconds"),
            "options_region": params["region"],
        }
        if is_pilot:
            entry["reason"] = (
                "region-restricted or scope differs from the reference scope (pilot, not whole-genome)"
            )
            report["pilots"].append(entry)
            continue
        policy = None
        for pol in experiment.policies.values():
            expected = comparator_parameters(experiment, ref, run, pol)
            if all(
                expected[k] == params.get(k)
                for k in (
                    "reference_transcript_selection",
                    "query_transcript_selection",
                )
            ):
                policy = pol
                break
        if policy is None:
            entry["reason"] = (
                f"transcript selections {params['reference_transcript_selection']}/{params['query_transcript_selection']} match no policy"
            )
            report["unmatched"].append(entry)
            continue
        versions = {
            k: (manifest.get("code_version") or {}).get(k)
            for k in ("python", "pandas", "pyranges1")
        }
        implementation = (
            implementation_from_snapshot(snapshot, versions)
            if snapshot
            else {"fingerprint": None}
        )
        components = key_components(inputs, params, implementation)
        key = cache_key(components)
        expected_inputs = {
            "query_sha256": inputs["query_sha256"],
            "reference_sha256": prep["annotation"]["output"]["sha256"],
            "genome_sha256": (
                prep["genome"]["output"]["sha256"] if prep.get("genome") else None
            ),
            "scope_sha256": prep["scope"]["sha256"],
            "seqname_map_sha256": (prep.get("seqname_map") or {}).get("sha256"),
        }
        expected = key_components(
            expected_inputs,
            comparator_parameters(experiment, ref, run, policy),
            current,
        )
        expected_key = cache_key(expected)
        missing = [n for n in EXPECTED_OUTPUTS if not (out_dir / n).exists()]
        status = (
            STATUS_COMPLETE if key == expected_key and not missing else STATUS_STALE
        )
        record = {
            "status": status,
            "source": "imported",
            "run_id": run.run_id,
            "policy_id": policy.id,
            "reference_id": ref.id,
            "output_dir": str(out_dir),
            "cache_key": key,
            "key_components": components,
            "expected_key": expected_key,
            "difference": (
                None
                if key == expected_key
                else describe_difference(components, expected)
            ),
            "missing_outputs": missing,
            "argv": manifest.get("argv"),
            "started": None,
            "finished": manifest.get("created"),
            "wall_seconds": manifest.get("runtime_seconds"),
            "peak_rss_bytes": None,
            "exit_code": 0,
            "manifest_code_version": manifest.get("code_version"),
            "implementation": implementation,
            "imported": now(),
        }
        k = (run.run_id, policy.id)
        if k in chosen:
            report["conflicts"].append(
                {
                    "run_id": run.run_id,
                    "policy_id": policy.id,
                    "dirs": [chosen[k]["output_dir"], str(out_dir)],
                }
            )
            # Prefer a result whose key verifies (complete) over a stale one; then the newest.
            prev = chosen[k]
            if (status == STATUS_COMPLETE, record["finished"] or "") <= (
                prev["status"] == STATUS_COMPLETE,
                prev["finished"] or "",
            ):
                continue
            report["matched"] = [
                m
                for m in report["matched"]
                if not (m["run_id"] == run.run_id and m.get("policy_id") == policy.id)
            ]
        chosen[k] = record
        entry.update(
            policy_id=policy.id, status=status, difference=record["difference"]
        )
        report["matched"].append(entry)

    for (run_id, policy_id), record in chosen.items():
        existing = ws.run_record(run_id, policy_id)
        if (
            existing
            and existing.get("source") == "computed"
            and existing.get("status") == STATUS_COMPLETE
        ):
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
                if k not in ("implementation", "key_components")
            },
        )

    # Independent analysis tables (optional)
    pattern = experiment.imports.get("independent_tables")
    report["independent"] = []
    if pattern:
        for run in experiment.runs.values():
            for policy in experiment.policies.values():
                for alias in run.aliases:
                    path = benchmark_dir / pattern.format(
                        alias=alias, policy=policy.id, run_id=run.run_id
                    )
                    if path.exists():
                        target = (
                            ws.independent_dir(run.run_id) / f"{policy.id}.import.json"
                        )
                        write_json(
                            target,
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

    write_json(ws.import_report(), report)
    log(
        f"Imported {len(chosen)} whole-genome results, {len(report['pilots'])} pilots kept separate, "
        f"{len(report['unmatched'])} unmatched, {len(report['conflicts'])} conflicts"
    )
    stale = [r for r in chosen.values() if r["status"] != STATUS_COMPLETE]
    for record in stale:
        log(f"  STALE {record['run_id']}/{record['policy_id']}: {record['difference']}")
    return report

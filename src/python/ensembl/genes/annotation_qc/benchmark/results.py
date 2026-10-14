"""
Assess one comparator output directory against what the configuration expects.

Used by ``import`` and by run planning, so a result has one status whatever way
it entered the workspace (computed here or imported from elsewhere):

    complete     contract result; inputs, parameters and definition version match
                 the configuration; every output file matches the SHA-256 in its
                 manifest; the implementation is verified (a registered build) or
                 identical to this checkout
    unverified   as complete, but the implementation cannot be identified
                 (benchmark/registry.py). Shown, labelled "unverified
                 implementation", and recomputed by ``benchmark run``
    legacy       written before output contract 1.0 (no contract_version), with
                 matching inputs and parameters. Shown only through the explicit
                 legacy importer (legacy.py): identical-definition metrics only
    stale        anything else that does not match (the reason names what), or an
                 output whose checksum differs from the manifest
    failed / not_run / unavailable   as before (runs.py)

The cache key covers inputs (checksums), comparator parameters and the
*definition* (contract id, contract major version, definition version), not the
code fingerprint: code identity is assessed separately, so a verified historical
build stays valid when only the dashboard or display code changes.
"""

from __future__ import annotations

from pathlib import Path

from ensembl.genes.annotation_qc.benchmark import registry as reg
from ensembl.genes.annotation_qc.benchmark.config import stop_codon_options
from ensembl.genes.annotation_qc.benchmark.provenance import stable_hash
from ensembl.genes.annotation_qc.benchmark.workspace import read_json
from ensembl.genes.annotation_qc.metrics.pairwise import contract

STATUS_COMPLETE = "complete"
STATUS_UNVERIFIED = "unverified"
STATUS_LEGACY = "legacy"
STATUS_FAILED = "failed"
STATUS_NOT_RUN = "not_run"
STATUS_STALE = "stale"
STATUS_UNAVAILABLE = "unavailable"
DISPLAYED = (STATUS_COMPLETE, STATUS_UNVERIFIED, STATUS_LEGACY)

MANIFEST = "comparison_manifest.json"
LEGACY_OUTPUTS = (
    "comparison_summary.json",
    "comparison_details.tsv",
    "consensus_transcript_labels.tsv",
    "gene_splits.tsv",
    "gene_merges.tsv",
    "reference_filter_audit.json",
    MANIFEST,
)
CONTRACT_OUTPUTS = LEGACY_OUTPUTS + (
    "metrics.tsv",
    "evaluated_transcripts.tsv.gz",
    "per_sequence.tsv",
    "cds_exact_one_to_one_pairs.tsv",
)
SPLICE_OUTPUTS = ("splice_sites_by_position.tsv", "cds_introns_by_length.tsv")
STOP_CODON_OUTPUTS = (
    "stop_codon_harmonisation.tsv",
    "query_evaluated.stop_harmonised.gtf",
)
INPUT_KEYS = (
    "query_sha256",
    "reference_sha256",
    "genome_sha256",
    "scope_sha256",
    "seqname_map_sha256",
)


def expected_outputs(parameters: dict) -> tuple[str, ...]:
    """Files a successful contract comparison with these parameters must have written."""
    out: tuple[str, ...] = CONTRACT_OUTPUTS
    if parameters.get("reference_transcript_selection") != "all":
        out += SPLICE_OUTPUTS
    if parameters.get("add_stop_codon"):
        out += STOP_CODON_OUTPUTS
    return out


def current_definition() -> dict:
    return {
        "contract_id": contract.CONTRACT_ID,
        "contract_major": contract.CONTRACT_VERSION.split(".")[0],
        "definition_version": contract.DEFINITION_VERSION,
    }


def manifest_definition(manifest: dict) -> dict:
    version = manifest.get("contract_version")
    return {
        "contract_id": manifest.get("contract_id"),
        "contract_major": version.split(".")[0] if version else None,
        "definition_version": manifest.get("definition_version"),
    }


def is_contract(manifest: dict | None) -> bool:
    return bool(manifest and manifest.get("contract_version"))


def _biotypes(value) -> list | None:
    if not value:
        return None
    items = value if isinstance(value, list) else str(value).split(",")
    return sorted(item.strip() for item in items if item.strip()) or None


def manifest_stop_codon_options(manifest: dict) -> dict:
    """Stop-codon parameters recorded in a manifest ({} = original coordinates)."""
    opts = manifest.get("options", {})
    try:
        return stop_codon_options(
            opts.get("add_stop_codon", False),
            opts.get("stop_codon_genetic_code"),
            opts.get("stop_codon_sequence_code"),
        )
    except ValueError as error:
        return {"add_stop_codon": True, "unreadable": str(error)}


def manifest_parameters(manifest: dict) -> dict:
    opts = manifest.get("options", {})
    scopes = (manifest.get("scope") or {}).get("analysis_scopes") or {}
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
        "analysis_scopes": {k: sorted(v) for k, v in scopes.items()},
        **manifest_stop_codon_options(manifest),
    }


def manifest_inputs(manifest: dict) -> dict:
    inputs = manifest.get("inputs", {})

    def sha(key):
        return (inputs.get(key) or {}).get("sha256")

    return {
        "query_sha256": sha("query"),
        "reference_sha256": sha("reference"),
        "genome_sha256": sha("genome"),
        "scope_sha256": sha("regions_file")
        or ("all" if not manifest.get("options", {}).get("region") else None),
        "seqname_map_sha256": sha("seqname_map"),
    }


def key_components(inputs: dict, parameters: dict, definition: dict) -> dict:
    return {
        "inputs": {k: inputs.get(k) for k in INPUT_KEYS},
        "parameters": parameters,
        "definition": definition,
    }


def cache_key(components: dict) -> str:
    return stable_hash(components)


def describe_difference(old: dict, new: dict, ignore=()) -> str:
    changed = []
    for section in ("inputs", "parameters", "definition"):
        a, b = old.get(section) or {}, new.get(section) or {}
        for k in sorted(set(a) | set(b)):
            if f"{section}.{k}" not in ignore and a.get(k) != b.get(k):
                changed.append(f"{section}.{k}")
    return "changed: " + (", ".join(changed) if changed else "cache key")


def verify_outputs(out_dir: Path, manifest: dict, checksums) -> dict:
    """Every file listed in the manifest's ``outputs`` against its SHA-256."""
    listed = manifest.get("outputs") or {}
    if not listed:
        return {"status": "not_recorded", "files": 0}
    missing, differing = [], []
    for name, entry in sorted(listed.items()):
        path = Path(out_dir) / name
        if not path.exists():
            missing.append(name)
        elif checksums.sha256(path) != entry.get("sha256"):
            differing.append(name)
    return {
        "status": "verified" if not (missing or differing) else "mismatch",
        "files": len(listed),
        "missing": missing,
        "differing": differing,
    }


def assess(
    out_dir: Path,
    expected: dict,
    registry: reg.Registry,
    checksums,
    snapshot_identity: dict | None = None,
) -> dict:
    """
    Status of the result in ``out_dir`` for the expected key components.

    Returns {status, reason, components, identity, outputs_check, contract}.
    """
    out_dir = Path(out_dir)
    manifest = read_json(out_dir / MANIFEST)
    if manifest is None:
        return {"status": STATUS_STALE, "reason": f"{MANIFEST} missing or unreadable"}
    components = key_components(
        manifest_inputs(manifest),
        manifest_parameters(manifest),
        manifest_definition(manifest),
    )
    base = {
        "components": components,
        "contract": is_contract(manifest),
        "contract_version": manifest.get("contract_version"),
        "definition_version": manifest.get("definition_version"),
        "provenance": manifest.get("provenance"),
        "created": manifest.get("created"),
        "runtime_seconds": manifest.get("runtime_seconds"),
    }
    if not is_contract(manifest):
        ignore = ("parameters.analysis_scopes",) + tuple(
            f"definition.{k}" for k in current_definition()
        )
        diff = describe_difference(components, expected, ignore)
        missing = [n for n in LEGACY_OUTPUTS if not (out_dir / n).exists()]
        identity = {
            "status": reg.UNVERIFIED,
            "reason": "pre-contract manifest: no implementation fingerprint",
            **({"code_snapshot": snapshot_identity} if snapshot_identity else {}),
        }
        if diff == "changed: cache key" and not missing:
            return {
                **base,
                "status": STATUS_LEGACY,
                "reason": "written before output contract 1.0; shown through the legacy importer "
                "(identical-definition metrics only)",
                "identity": identity,
                "outputs_check": {"status": "not_recorded"},
            }
        return {
            **base,
            "status": STATUS_STALE,
            "identity": identity,
            "reason": "; ".join(
                filter(
                    None,
                    [
                        diff if diff != "changed: cache key" else None,
                        f"outputs missing: {', '.join(missing)}" if missing else None,
                    ],
                )
            ),
        }
    identity = reg.identify(manifest, registry)
    if components != expected:
        return {
            **base,
            "status": STATUS_STALE,
            "identity": identity,
            "reason": describe_difference(components, expected),
        }
    check = verify_outputs(out_dir, manifest, checksums)
    missing = [
        n
        for n in expected_outputs(components["parameters"])
        if n not in (manifest.get("outputs") or {}) and n != MANIFEST
    ]
    if check["status"] != "verified" or missing:
        problems = []
        if check.get("missing"):
            problems.append(f"files missing: {', '.join(check['missing'])}")
        if check.get("differing"):
            problems.append(
                f"files differ from their manifest SHA-256: {', '.join(check['differing'])}"
            )
        if check["status"] == "not_recorded":
            problems.append("the manifest lists no output checksums")
        if missing:
            problems.append(f"expected outputs not written: {', '.join(missing)}")
        return {
            **base,
            "status": STATUS_STALE,
            "identity": identity,
            "outputs_check": check,
            "reason": "; ".join(problems),
        }
    if identity["status"] == reg.UNVERIFIED:
        return {
            **base,
            "status": STATUS_UNVERIFIED,
            "identity": identity,
            "outputs_check": check,
            "reason": f"unverified implementation: {identity['reason']}",
        }
    return {
        **base,
        "status": STATUS_COMPLETE,
        "identity": identity,
        "outputs_check": check,
        "reason": f"{reg.label(identity)}; {check['files']} outputs match their manifest checksums",
    }

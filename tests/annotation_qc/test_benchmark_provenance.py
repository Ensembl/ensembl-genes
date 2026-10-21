"""Implementation registry, legacy import mapping and exact model reconstruction (SYNTHETIC fixtures)."""

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from ensembl.genes.annotation_qc.benchmark import legacy
from ensembl.genes.annotation_qc.benchmark import registry as reg
from ensembl.genes.annotation_qc.benchmark.models import ModelsUnavailable, rebuild
from ensembl.genes.annotation_qc.metrics.pairwise import contract


def _impl(files, scheme=None, versions=None):
    versions = versions or {
        "pandas": "3.0.6",
        "pyranges1": "1.4.5",
        "python": "3.12.14",
    }
    out = {
        "fingerprint": reg.recompute_fingerprint(files, versions),
        "source_files": files,
        "versions": versions,
    }
    if scheme:
        out["scheme"] = scheme
    return out


FILES = {"metrics/pairwise/classify.py": "a" * 64, "parsers/annotation.py": "b" * 64}


def _registry(tmp_path, entries):
    path = tmp_path / "registry.json"
    path.write_text(json.dumps({"builds": entries}))
    registry = reg.Registry()
    registry.add_file(path)
    return registry


def _entry(build_id, impl, scheme):
    return {
        "build_id": build_id,
        "scheme": scheme,
        **{k: impl[k] for k in ("fingerprint", "source_files", "versions")},
    }


# --------------------------------------------------------------- registry
def test_packaged_registry_identifies_the_contract_build():
    registry = reg.Registry.load()
    builds = {(b["scheme"], b["build_id"]) for b in registry.public()}
    assert (reg.SCHEME_BY_NAME, "d445054") in builds and (
        reg.SCHEME_BY_PATH,
        "f0f838d",
    ) in builds
    d445054 = next(b for b in registry.builds.values() if b["build_id"] == "d445054")
    manifest = {
        "implementation": {
            k: d445054[k] for k in ("fingerprint", "source_files", "versions")
        }
    }
    identity = reg.identify(manifest, registry)  # no "scheme" key: scheme 1
    assert (identity["status"], identity["build_id"], identity["scheme_recorded"]) == (
        "verified",
        "d445054",
        False,
    )
    assert any(e["type"] == "wheel" for e in identity["evidence"])


def test_identity_is_keyed_by_scheme_and_fingerprint(tmp_path):
    scheme1 = _impl(FILES)
    registry = _registry(tmp_path, [_entry("old", scheme1, reg.SCHEME_BY_NAME)])
    assert reg.identify({"implementation": scheme1}, registry)["status"] == "verified"
    # the same fingerprint under another scheme is a different build: never compared across schemes
    other = {**scheme1, "scheme": reg.SCHEME_BY_PATH}
    identity = reg.identify({"implementation": other}, registry)
    assert (
        identity["status"] == "unverified"
        and "differ from this checkout" in identity["reason"]
    )
    # an unknown scheme is unsupported
    assert (
        "unsupported scheme"
        in reg.identify({"implementation": {**scheme1, "scheme": "x/9"}}, registry)[
            "reason"
        ]
    )


@pytest.mark.parametrize(
    "change, reason",
    [
        (lambda i: i.update(fingerprint="0" * 64), "does not recompute"),
        (
            lambda i: i["source_files"].update(
                {"metrics/pairwise/classify.py": "c" * 64}
            ),
            "does not recompute",
        ),
        (lambda i: i.pop("fingerprint"), "pre-contract"),
    ],
)
def test_inconsistent_or_missing_implementations_are_unverified(
    tmp_path, change, reason
):
    impl = _impl(dict(FILES))
    registry = _registry(
        tmp_path, [_entry("old", _impl(dict(FILES)), reg.SCHEME_BY_NAME)]
    )
    change(impl)
    identity = reg.identify({"implementation": impl}, registry)
    assert identity["status"] == "unverified" and reason in identity["reason"]
    assert reg.label(identity) == "unverified implementation"


def test_registered_pair_with_other_files_is_unverified(tmp_path):
    registered = _impl(dict(FILES))
    registry = _registry(tmp_path, [_entry("old", registered, reg.SCHEME_BY_NAME)])
    # forge a manifest that claims the registered fingerprint with other files: recompute fails first
    forged = {**registered, "source_files": {**FILES, "x.py": "d" * 64}}
    assert reg.identify({"implementation": forged}, registry)["status"] == "unverified"
    # a registry entry whose files do not reproduce its fingerprint is refused
    with pytest.raises(reg.RegistryError, match="does not match its file hashes"):
        registry.add(
            {
                **_entry("again", registered, reg.SCHEME_BY_NAME),
                "source_files": {"y.py": "e" * 64},
                "fingerprint": registered["fingerprint"],
            },
            "test",
        )


def test_scheme2_identical_to_this_checkout_is_current_code():
    rel = "metrics/pairwise/contract.py"
    digest = hashlib.sha256((reg.PACKAGE_DIR / rel).read_bytes()).hexdigest()
    impl = _impl({rel: digest}, scheme=reg.SCHEME_BY_PATH)
    identity = reg.identify({"implementation": impl}, reg.Registry())
    assert (
        identity["status"] == "current_checkout"
        and reg.label(identity) == "current code (this checkout)"
    )
    # scheme 1 manifests are never matched against the checkout
    assert (
        reg.identify({"implementation": _impl({rel: digest})}, reg.Registry())["status"]
        == "unverified"
    )


def test_register_build_checks_every_file_against_git(tmp_path):
    repo = tmp_path / "repo"
    pkg = repo / reg.PACKAGE_REPO_PATH / "metrics"
    pkg.mkdir(parents=True)
    (pkg / "m.py").write_text("x = 1\n")
    for args in (
        ["init", "-q"],
        ["add", "."],
        ["-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "c"],
    ):
        subprocess.run(["git", "-C", str(repo), *args], check=True)
    good = hashlib.sha256(b"x = 1\n").hexdigest()
    manifest = tmp_path / "m.json"
    manifest.write_text(
        json.dumps(
            {
                "contract_version": "1.0",
                "implementation": _impl({"metrics/m.py": good}, reg.SCHEME_BY_PATH),
            }
        )
    )
    entry = reg.register_build([manifest], repo, "HEAD", "build1")
    assert (
        entry["scheme"] == reg.SCHEME_BY_PATH
        and entry["evidence"][0]["type"] == "git_objects"
    )
    manifest.write_text(
        json.dumps(
            {"implementation": _impl({"metrics/m.py": "f" * 64}, reg.SCHEME_BY_PATH)}
        )
    )
    with pytest.raises(reg.RegistryError, match="files differ"):
        reg.register_build([manifest], repo, "HEAD", "build1")


# --------------------------------------------------------------- legacy import
def _summary(**extra):
    return {
        "total_reference_genes": 10,
        "total_consensus_genes": 0,
        "sensitivity": {
            "locus_detected_count": 4,
            "missed_count": 6,
            "strand_mismatch_count": 0,
            "exon_coordinate_exact_count": 1,
            "exact_match_count": 1,
        },
        "sensitivity_cds": {
            "cds_coordinate_exact_count": 2,
            "cds_exact_match_count": 2,
            "cds_intron_chain_recovered": 3,
            "multi_segment_cds_reference_genes": 8,
        },
        "specificity": {
            "consensus_cds_coordinate_exact_count": 0,
            "novel_consensus_count": 0,
            "strand_mismatch_consensus_count": 0,
        },
        "split_merge": {
            "gene_split_count": 0,
            "gene_merge_count": 0,
            "one_to_one_reference_genes": 3,
        },
        "novel_categories": {},
        "locus_detection_exonic": {"span_only_detected_count": 0},
        "intron_support": {
            "cds": {"shared_introns": 5, "reference_introns": 20, "query_introns": 0}
        },
        **extra,
    }


def test_legacy_maps_identical_definitions_only():
    rows = {r["metric_id"]: r for r in legacy.legacy_metric_rows(_summary())}
    assert set(rows) == set(contract.METRICS)  # every contract metric has a row
    exact = rows["cds_coordinate_exact_reference"]
    assert (
        exact["numerator"],
        exact["denominator"],
        exact["value"],
        exact["source"],
    ) == (2, 10, 0.2, "legacy_summary")
    # zero denominators are undefined, never 0
    q = rows["cds_coordinate_exact_query"]
    assert (q["status"], q["reason"], q["value"]) == (
        "undefined",
        "zero_denominator",
        None,
    )
    assert rows["query_unique_cds_introns_supported"]["status"] == "undefined"
    # contract-only analyses and the pre-matching summary's one-to-one fields are not available
    for mid in (
        "cds_outcome_exact",
        "cds_intron_chain_any_pair",
        "exact_cds_one_to_one_f1",
        "locus_recovery_cds_overlap",
    ):
        assert (
            rows[mid]["status"] == "not_available" and rows[mid]["value"] is None
        ), mid
    with_block = {
        r["metric_id"]: r
        for r in legacy.legacy_metric_rows(
            _summary(
                cds_exact_one_to_one={"true_positives": 2, "exact_gene_pairs": 3},
                total_consensus_genes=6,
            )
        )
    }
    f1 = with_block["exact_cds_one_to_one_f1"]
    assert (f1["numerator"], f1["denominator"], f1["value"]) == (4, 16, 0.25)


def test_historical_estimate_never_fills_one_to_one_ids():
    details = pd.DataFrame(
        {
            "source": ["reference"] * 3,
            "gene_id": ["r1", "r2", "r3"],
            "cds_coordinate_exact": ["True", "True", "False"],
            "cds_matched_id": ["q1", "q1", ""],
        }
    )
    rows = {r["metric_id"]: r for r in legacy.historical_rows(details, 3, 2)}
    assert set(rows) == set(legacy.HISTORICAL_IDS)
    assert (
        rows["historical_exact_tp_best_pair"]["numerator"] == 1
    )  # min(2 reference, 1 query)
    assert all(r["source"] == "historical" for r in rows.values())


# --------------------------------------------------------------- model reconstruction
def _parsed():
    rows = [
        ("1", 0, 100, "+", "gene", "g1", ""),
        ("1", 0, 100, "+", "transcript", "g1", "t1"),
        ("1", 0, 30, "+", "exon", "g1", "t1"),
        ("1", 60, 90, "+", "exon", "g1", "t1"),
        ("1", 0, 30, "+", "CDS", "g1", "t1"),
        ("1", 60, 90, "+", "CDS", "g1", "t1"),
        ("1", 0, 100, "+", "transcript", "g1", "t2"),
        ("1", 0, 100, "+", "exon", "g1", "t2"),
    ]
    df = pd.DataFrame(
        rows,
        columns=[
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "Feature",
            "gene_id",
            "transcript_id",
        ],
    )
    df["source_feature"] = df.Feature.map({"exon": "inferred_exon"}).fillna("x")
    df["tags"] = ""
    return df


def _export(transformed=""):
    return pd.DataFrame(
        {
            "transcript_id": ["t1", "t2"],
            "gene_id": ["g1", "g1"],
            "original_gene_id": ["g1", "g1"],
            "original_transcript_id": ["t1", "t2"],
            "chrom": ["1", "1"],
            "strand": ["+", "+"],
            "start": ["1", "1"],
            "end": ["93" if transformed else "90", "100"],
            "transcript_biotype": ["protein_coding", "lncRNA"],
            "status": ["evaluated", "excluded"],
            "exclusion_reason": ["", "transcript_biotype_filter"],
            "selected_by": ["all", ""],
            "exon_count": ["2", "1"],
            "cds_segment_count": ["2", "0"],
            "cds_length": ["63" if transformed else "60", "0"],
            "transformed": [transformed, ""],
        }
    )


def test_rebuild_keeps_exported_transcripts_without_reselection():
    models = rebuild(_parsed(), _export())
    assert list(models.transcript_id) == [
        "t1",
        "t2",
    ]  # excluded isoforms kept for display
    assert models.set_index("transcript_id").at["t1", "cds"] == "1-30,61-90"


def test_rebuild_applies_the_audit_and_refuses_without_it():
    audit = pd.DataFrame(
        {
            "transcript_id": ["t1"],
            "status": ["extended"],
            "terminal_cds_start": ["61"],
            "terminal_cds_end": ["90"],
            "new_terminal_cds_start": ["61"],
            "new_terminal_cds_end": ["93"],
        }
    )
    models = rebuild(_parsed(), _export("stop_codon_added"), audit)
    t1 = models.set_index("transcript_id").loc["t1"]
    assert (t1.cds, t1.exons, t1.transformed) == (
        "1-30,61-93",
        "1-30,61-93",
        "stop_codon_added",
    )
    with pytest.raises(ModelsUnavailable, match="audit file is absent"):
        rebuild(_parsed(), _export("stop_codon_added"), None)


def test_rebuild_refuses_counts_that_differ_from_the_export():
    bad = _export()
    bad.loc[0, "cds_length"] = "61"
    with pytest.raises(ModelsUnavailable, match="cds_length"):
        rebuild(_parsed(), bad)

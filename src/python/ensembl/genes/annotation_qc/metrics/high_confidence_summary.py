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
"""Pure confidence-evidence summaries from normalized DataFrames."""
from __future__ import annotations
from collections import Counter
import itertools
import math
from statistics import mean, median
import pandas as pd
from ensembl.genes.annotation_qc.config.high_confidence import (
    FLAGS,
    CONFIDENCE_POLICY_VERSION,
    STATES,
    CLASSES,
    BINARY,
    NUMERIC,
    CATEGORICAL,
    MISSING,
)
from ensembl.genes.annotation_qc.metrics.high_confidence import confidence_rating


def state(value):
    """Map a binary evidence value to pass, fail or unknown."""
    if pd.isna(value):
        return "unknown"
    if value in ("1", 1, True):
        return "pass"
    if value in ("0", 0, False):
        return "fail"
    if value in MISSING:
        return "unknown"
    raise ValueError(f"Expected 1/0/. binary evidence; received {value!r}")


def number(value):
    """Parse an optional finite numeric evidence value."""
    if pd.isna(value) or value in MISSING:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite numeric evidence: {value}")
    return result


def pct(n, total):
    """Return percentage with an explicit empty value for zero denominator."""
    return round(100 * n / total, 4) if total else ""


def quantile(values, q):
    """Interpolate a quantile from sorted observed values."""
    if not values:
        return ""
    index = (len(values) - 1) * q
    lo, hi = math.floor(index), math.ceil(index)
    return values[lo] + (values[hi] - values[lo]) * (index - lo)


def distribution(rows, model, field, values, group="all"):
    """Append descriptive statistics for an observed metric."""
    counts = Counter(values)
    return [
        {
            "model": model,
            "group": group,
            "evidence": field,
            "category": k,
            "count": v,
            "denominator": len(rows),
            "percent": pct(v, len(rows)),
        }
        for k, v in sorted(counts.items())
    ]


# Keep the established reconstruction/scoring flow together for parity.
# pylint: disable-next=too-many-branches,too-many-locals,too-many-statements
def summarize_confidence(
    frame, model, manifest=None, intron_frame=None, evidence_hash=None
):
    """Compute evidence distributions and mode comparisons from parsed tables only."""
    fields, rows = list(frame.columns), frame.to_dict("records")
    for row in rows:
        for key, value in row.items():
            row[key] = (
                "."
                if value is None or pd.isna(value)
                else ("1" if value else "0") if isinstance(value, bool) else str(value)
            )
    path = model
    required = [
        "gene_id",
        "transcript_id",
        "confidence",
        "confidence_rule",
        "homology_status",
    ] + FLAGS
    missing = set(required) - set(fields)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")
    ids = [r["transcript_id"] for r in rows]
    if any(not tx or tx in MISSING for tx in ids) or len(set(ids)) != len(ids):
        raise ValueError(f"{path}: missing or duplicate transcript IDs")
    if any(r["confidence"] not in CLASSES for r in rows):
        raise ValueError(f"{path}: unexpected confidence label")
    if any(r["homology_status"] not in {"hit", "no_hit", "not_measured"} for r in rows):
        raise ValueError(f"{path}: unexpected homology_status")
    binary = [f for f in BINARY if f in fields]
    numeric = [f for f in NUMERIC if f in fields]
    for r in rows:
        for f in binary:
            state(r[f])
        for f in numeric:
            number(r[f])
    warnings = []
    if not rows:
        warnings.append("No transcript rows; percentages and statistics are undefined.")
    if len(set(r["gene_id"] for r in rows)) != len(rows):
        warnings.append(
            "Multiple rows per gene; percentages use transcripts, not genes."
        )
    if manifest is None:
        warnings.append(
            "Run manifest unavailable; scoring policy and input comparability cannot be verified."
        )
    elif manifest.get("processed_transcript_count") != len(rows):
        warnings.append(
            "Manifest processed-transcript count differs from evidence table."
        )
    if (
        manifest
        and manifest.get("options", {}).get("no_utr") is True
        and manifest.get("confidence_policy_version") != CONFIDENCE_POLICY_VERSION
    ):
        warnings.append(
            (
                "This no_utr run predates the RNA-required confidence policy v2; "
                "rescore evidence before making current-policy figures."
            )
        )
    if not {"intron_count", "cds_intron_count"} <= set(fields):
        warnings.append(
            "Intron counts unavailable; intronless strata cannot be recovered."
        )

    annotated = []
    for r in rows:
        categories = {s: [f for f in FLAGS if state(r[f]) == s] for s in STATES}
        annotated.append(
            # Keyword construction rejects duplicate export-schema keys.
            # pylint: disable-next=use-dict-literal
            dict(
                model=model,
                **r,
                evidence_pattern="|".join(r[f] for f in FLAGS),
                **{f"{s}_flag_count": len(categories[s]) for s in STATES},
                **{f"{s}_flags": ";".join(categories[s]) for s in STATES},
                support_applicability=(
                    "no_cds_introns"
                    if r.get("cds_intron_count") == "0"
                    else (
                        "has_cds_introns"
                        if number(r.get("cds_intron_count", ".")) is not None
                        else "unknown"
                    )
                ),
            )
        )
    added = (
        ["evidence_pattern"]
        + [f"{s}_flag_count" for s in STATES]
        + [f"{s}_flags" for s in STATES]
        + ["support_applicability"]
    )

    flag_rows, category_rows, numeric_rows, patterns, joint, strata = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    groups = [("all", rows)] + [
        (c, [r for r in rows if r["confidence"] == c]) for c in CLASSES
    ]
    if "cds_intron_count" in fields:
        groups += [
            ("no_cds_introns", [r for r in rows if number(r["cds_intron_count"]) == 0]),
            (
                "has_cds_introns",
                [
                    r
                    for r in rows
                    if number(r["cds_intron_count"]) is not None
                    and number(r["cds_intron_count"]) > 0
                ],
            ),
        ]
    for group, subset in groups:
        for f in binary:
            counts = Counter(state(r[f]) for r in subset)
            measured = counts["pass"] + counts["fail"]
            for s in STATES:
                flag_rows.append(
                    {
                        "model": model,
                        "group": group,
                        "flag": f,
                        "state": s,
                        "count": counts[s],
                        "denominator": len(subset),
                        "percent": pct(counts[s], len(subset)),
                        "measured_count": measured,
                        "pass_percent_measured": pct(counts["pass"], measured),
                    }
                )
        for f in CATEGORICAL:
            if f in fields:
                category_rows.extend(
                    distribution(
                        subset,
                        model,
                        f,
                        [r[f] if r[f] not in MISSING else "unknown" for r in subset],
                        group,
                    )
                )
        for f in numeric:
            values = sorted(v for r in subset if (v := number(r[f])) is not None)
            numeric_rows.append(
                {
                    "model": model,
                    "group": group,
                    "metric": f,
                    "n": len(values),
                    "missing": len(subset) - len(values),
                    "mean": mean(values) if values else "",
                    "min": values[0] if values else "",
                    "q25": quantile(values, 0.25),
                    "median": median(values) if values else "",
                    "q75": quantile(values, 0.75),
                    "p95": quantile(values, 0.95),
                    "max": values[-1] if values else "",
                }
            )
    for (pattern, confidence), count in sorted(
        Counter((tuple(r[f] for f in FLAGS), r["confidence"]) for r in rows).items()
    ):
        patterns.append(
            # Keyword construction rejects duplicate export-schema keys.
            # pylint: disable-next=use-dict-literal
            dict(
                model=model,
                **dict(zip(FLAGS, pattern)),
                confidence=confidence,
                count=count,
                percent=pct(count, len(rows)),
            )
        )
    for a, b in itertools.combinations(FLAGS, 2):
        counts = Counter((state(r[a]), state(r[b])) for r in rows)
        for sa, sb in itertools.product(STATES, repeat=2):
            joint.append(
                {
                    "model": model,
                    "flag_a": a,
                    "flag_b": b,
                    "state_a": sa,
                    "state_b": sb,
                    "count": counts[sa, sb],
                    "percent": pct(counts[sa, sb], len(rows)),
                }
            )
    for field in [
        "seqname",
        "homology_status",
        "utr_annotation_source",
        "biotype_source",
        "cds_phase_source",
    ]:
        if field in fields:
            counts = Counter((r[field], r["confidence"]) for r in rows)
            totals = Counter(r[field] for r in rows)
            for category in sorted(totals):
                for c in CLASSES:
                    strata.append(
                        {
                            "model": model,
                            "stratum": field,
                            "category": category,
                            "confidence": c,
                            "count": counts[category, c],
                            "denominator": totals[category],
                            "percent": pct(counts[category, c], totals[category]),
                        }
                    )
    for field in ["intron_count", "cds_intron_count"]:
        if field in fields:

            def buckets(r, field=field):
                return (
                    "unknown"
                    if number(r[field]) is None
                    else (
                        "0"
                        if number(r[field]) == 0
                        else (
                            "1"
                            if number(r[field]) == 1
                            else "2-5" if number(r[field]) <= 5 else "6+"
                        )
                    )
                )

            totals = Counter(buckets(r) for r in rows)
            counts = Counter((buckets(r), r["confidence"]) for r in rows)
            for category, total in sorted(totals.items()):
                for c in CLASSES:
                    strata.append(
                        {
                            "model": model,
                            "stratum": field,
                            "category": category,
                            "confidence": c,
                            "count": counts[category, c],
                            "denominator": total,
                            "percent": pct(counts[category, c], total),
                        }
                    )

    sensitivity, transitions = [], Counter()
    can_score = all(f in fields for f in ["error", "score"])
    for r in rows:
        if not can_score:
            break
        evidence = {
            f: {"pass": True, "fail": False, "unknown": None}[state(r[f])]
            for f in FLAGS
        }
        hit = (
            None
            if r["homology_status"] != "hit"
            else {f: number(r[f]) for f in ["error", "score"]}
        )
        if hit is not None and any(v is None for v in hit.values()):
            raise ValueError(
                f"{path}: hit with missing error/score for {r['transcript_id']}"
            )
        standard, standard_rule = confidence_rating(evidence, hit)
        no_utr, no_utr_rule = confidence_rating(evidence, hit, no_utr=True)
        sensitivity.append(
            {
                "model": model,
                "gene_id": r["gene_id"],
                "transcript_id": r["transcript_id"],
                "reported_confidence": r["confidence"],
                "standard_confidence": standard,
                "no_utr_confidence": no_utr,
                "standard_rule": standard_rule,
                "no_utr_rule": no_utr_rule,
                "changed": int(standard != no_utr),
            }
        )
        transitions[standard, no_utr] += 1
    transition_rows = [
        {
            "model": model,
            "standard_confidence": a,
            "no_utr_confidence": b,
            "count": transitions[a, b],
            "denominator": len(sensitivity),
            "percent": pct(transitions[a, b], len(sensitivity)),
        }
        for a, b in itertools.product(CLASSES, repeat=2)
    ]
    if (
        can_score
        and manifest
        and isinstance(manifest.get("options", {}).get("no_utr"), bool)
    ):
        policy = (
            "no_utr_confidence"
            if manifest["options"]["no_utr"]
            else "standard_confidence"
        )
        mismatches = sum(r["reported_confidence"] != r[policy] for r in sensitivity)
        if mismatches:
            warnings.append(
                f"{mismatches} reported labels disagree with recomputed manifest scoring policy."
            )
    elif not can_score:
        warnings.append("Missing error/score columns; UTR sensitivity not computed.")

    tables = {
        "flag_distribution": (
            flag_rows,
            [
                "model",
                "group",
                "flag",
                "state",
                "count",
                "denominator",
                "percent",
                "measured_count",
                "pass_percent_measured",
            ],
        ),
        "category_distribution": (
            category_rows,
            [
                "model",
                "group",
                "evidence",
                "category",
                "count",
                "denominator",
                "percent",
            ],
        ),
        "numeric_summary": (
            numeric_rows,
            [
                "model",
                "group",
                "metric",
                "n",
                "missing",
                "mean",
                "min",
                "q25",
                "median",
                "q75",
                "p95",
                "max",
            ],
        ),
        "evidence_patterns": (
            patterns,
            ["model"] + FLAGS + ["confidence", "count", "percent"],
        ),
        "flag_pairs": (
            joint,
            ["model", "flag_a", "flag_b", "state_a", "state_b", "count", "percent"],
        ),
        "confidence_by_stratum": (
            strata,
            [
                "model",
                "stratum",
                "category",
                "confidence",
                "count",
                "denominator",
                "percent",
            ],
        ),
        "utr_sensitivity": (
            sensitivity,
            [
                "model",
                "gene_id",
                "transcript_id",
                "reported_confidence",
                "standard_confidence",
                "no_utr_confidence",
                "standard_rule",
                "no_utr_rule",
                "changed",
            ],
        ),
        "utr_transitions": (
            transition_rows,
            [
                "model",
                "standard_confidence",
                "no_utr_confidence",
                "count",
                "denominator",
                "percent",
            ],
        ),
    }

    intron_summary = []
    if intron_frame is not None:
        intron_fields, introns = list(intron_frame.columns), intron_frame.to_dict(
            "records"
        )
        introns_path = model
        required_intron = {
            "transcript_id",
            "is_cds_intron",
            "canonical",
            "regular_size",
            "supported",
            "splice_motif",
            "length",
            "read_count",
        }
        if not required_intron <= set(intron_fields):
            raise ValueError(
                f"{introns_path}: missing intron columns "
                f"{sorted(required_intron-set(intron_fields))}"
            )
        id_set = set(ids)
        if any(r["transcript_id"] not in id_set for r in introns):
            raise ValueError(f"{introns_path}: transcript absent from evidence table")
        for group, subset in [
            ("all_introns", introns),
            ("cds_introns", [r for r in introns if r["is_cds_intron"] == "1"]),
        ]:
            for f in ["canonical", "regular_size", "supported", "splice_motif"]:
                vals = [state(r[f]) if f != "splice_motif" else r[f] for r in subset]
                intron_summary.extend(distribution(subset, model, f, vals, group))
        tables["intron_distribution"] = (
            intron_summary,
            [
                "model",
                "group",
                "evidence",
                "category",
                "count",
                "denominator",
                "percent",
            ],
        )

    counts = Counter(r["confidence"] for r in rows)
    # Keyword construction rejects duplicate export-schema keys.
    # pylint: disable-next=use-dict-literal
    overview = dict(
        model=model,
        transcripts=len(rows),
        genes=len(set(r["gene_id"] for r in rows)),
        **{c: counts[c] for c in CLASSES},
        **{c + "_percent": pct(counts[c], len(rows)) for c in CLASSES},
        no_utr=(manifest or {}).get("options", {}).get("no_utr", "unknown"),
        utr_sensitive_transcripts=sum(r["changed"] for r in sensitivity),
        evidence_sha256=evidence_hash,
    )
    provenance = {
        "model": model,
        "run_manifest": manifest,
        "warnings": warnings,
        "interpretation": [
            "Unit is one selected transcript per row, not all predicted isoforms.",
            (
                "Unknown flags are separate from failures. Percentages use all "
                "rows unless marked measured."
            ),
            "Intronless transcripts pass intron checks vacuously; inspect intron-count strata.",
            "Pair and pattern counts describe co-occurrence, not independent evidence.",
            "DIAMOND error/score are helper metrics, not calibrated probabilities.",
            "Confidence is a heuristic; these summaries do not measure annotation accuracy.",
        ],
    }
    return {
        "model": model,
        "rows": rows,
        "overview": overview,
        "tables": tables,
        "provenance": provenance,
        "transcripts": annotated,
        "transcript_fields": ["model"] + fields + added,
    }


def routes(results):
    """Count matched confidence rules from previously scored evidence."""
    counts = Counter()
    for result in results:
        for row in result["rows"]:
            rule = row["confidence_rule"]
            if rule == "incomplete_cds_or_no_homology":
                incomplete = row["transcript_cds_complete"] != "1"
                missing = row["homology_status"] != "hit"
                rule = (
                    "Incomplete CDS + no homology hit"
                    if incomplete and missing
                    else "Incomplete CDS" if incomplete else "No homology hit"
                )
            elif rule == "no_high_or_medium_rule":
                if (
                    row["transcript_regular_intron_size"] != "1"
                    and row["transcript_canonical_splicing"] != "1"
                ):
                    rule = "Both intron checks fail"
                elif float(row["error"]) > 10:
                    rule = "Coverage error >10"
                elif row["transcript_intron_support"] != "1" and (
                    float(row["error"]) > 5 or float(row["score"]) < 80
                ):
                    rule = (
                        "RNA unmeasured"
                        if row["transcript_intron_support"] in {".", "", "NA"}
                        else "RNA unsupported"
                    ) + "; error >5 or score <80"
            counts[result["model"], rule, row["confidence"]] += 1
    return counts


def workflow_paths(results):
    """Count sequential decisions, validating them against the scoring implementation."""
    paths = Counter()
    for result in results:
        no_utr = result["overview"]["no_utr"]
        for row in result["rows"]:
            e = {
                f: True if row[f] == "1" else False if row[f] == "0" else None
                for f in FLAGS
            }
            hit = (
                {k: float(row[k]) for k in ["error", "score"]}
                if row["homology_status"] == "hit"
                else None
            )
            confidence = (
                confidence_rating(e, hit, no_utr=no_utr)[0]
                if isinstance(no_utr, bool)
                else row["confidence"]
            )
            if confidence != row["confidence"] or (
                confidence in {"high", "medium"}
                and (not e["transcript_cds_complete"] or hit is None)
            ):
                raise ValueError(
                    f"{result['model']}/{row['transcript_id']}: "
                    "reported confidence does not match scoring policy; "
                    "rescore evidence with annotation-qc high-confidence-rescore "
                    "before regenerating figures"
                )
            if not e["transcript_cds_complete"] or hit is None:
                path = ("method:" + result["model"], "required:no", "final:low")
            elif confidence == "high":
                path = (
                    "method:" + result["model"],
                    "required:yes",
                    "high:yes",
                    "final:high",
                )
            elif confidence == "medium":
                path = (
                    "method:" + result["model"],
                    "required:yes",
                    "high:no",
                    "medium:yes",
                    "final:medium",
                )
            else:
                path = (
                    "method:" + result["model"],
                    "required:yes",
                    "high:no",
                    "medium:no",
                    "final:low",
                )
            paths[result["model"], path, confidence] += 1
    return paths


def prepare_confidence_plots(results):
    """Validate routes and prepare diagram counts before the report layer runs."""
    prepared = []
    for result in results:
        result = dict(result)
        result["plot_flags"] = [
            f
            for f in FLAGS
            if not (
                result["overview"]["no_utr"] is True
                and f in {"transcript_five_prime_utr", "transcript_three_prime_utr"}
            )
        ]
        result["sankey_routes"] = routes([result])
        result["sankey_paths"] = workflow_paths([result])
        prepared.append(result)
    return prepared

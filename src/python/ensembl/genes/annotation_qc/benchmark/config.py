"""
Experiment configuration (format versions 1 and 2).

Version 2 adds optional keys to version 1; every version-1 file loads unchanged.
Version-2 additions (all optional):

    runs[].run_label        separates several runs of one tool ("GMB run 2")
    runs[].source           provenance of the compared file: original_path (as
                            recorded, never resolved or renamed), sha256,
                            recorded_in; validation checks the sha256
    runs[].legacy_ids       earlier run ids of the same prediction (linkage only)
    runs[].input_audit      JSON audit of the prediction file (e.g. which sequences
                            the predictor recorded as processed); shown, never used
                            to change a value
    import.results_dirs     directories of completed comparator results (headline)
    import.pilot_dirs       directories whose results are always pilots
    provenance.registries   extra implementation registries (benchmark/registry.py)
    display.preferred_tools tools listed first, in this order, before the ranked rest
    runs[].display_role     "focus" or "primary_comparator": the dataset of interest
                            and the main dataset it is compared with (presentation)

A configuration is a JSON file (or YAML, when PyYAML is installed through the
``dashboard`` extra). Relative paths are resolved against the directory that
contains the configuration file. See ``templates/experiment_template.json`` and
the README for every key.

Identity rules:
    reference ``id`` and run ``run_id`` are the unique identifiers. ``tool`` is a
    display/grouping label only, so several runs (models, versions, settings) of
    the same tool can coexist. Policy keys (e.g. "A", "B", "C") are identifiers too.

Comparator options per run (``runs[].comparator_options``):
    ``add_stop_codon`` (default false: original coordinates),
    ``stop_codon_genetic_code`` (NCBI table, default 1) and
    ``stop_codon_sequence_code`` ({sequence: table or "none"}) map to the
    pairwise-compare options of the same names. They are part of the run's
    parameters, command and cache key, so an original-coordinate run and a
    stop-harmonised run of the same file are distinct results.

Presentation only (never part of a cache key or comparison), except analysis scopes:
    ``runs[].companion_of``         a supplementary run (e.g. the same file with
                                    stop-codon harmonisation) is listed directly after
                                    this run and is excluded from the ordering
    ``runs[].display_role``         ``focus`` (the dataset of interest) or
                                    ``primary_comparator``; such runs are listed
                                    first (focus, then primary comparator), get
                                    reserved colours and are drawn on top, with
                                    heavier lines and their own marker shapes. A role
                                    whose run is absent is simply not shown
    ``display.pinned_runs``         run ids listed first, in this order
    ``display.preferred_tools``     tools listed next, in this order; the other runs
                                    follow in descending order of
                                    ``display.order_metric`` (default: one-to-one
                                    exact-CDS recall) for the selected policy and scope
    ``references[].analysis_scopes`` named sequence subsets ({id, label, sequences,
                                    description}); passed to the comparator as
                                    ``--analysis-scope`` (so part of the parameters and
                                    cache key) and read back from metrics.tsv.
                                    ``scope.label`` names the whole compared scope.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

CONFIG_VERSION = 2
CONFIG_VERSIONS = (1, 2)
SELECTIONS = ("all", "longest_cds", "canonical")
EVALUATION_MODES = ("all", "protein_coding", "cds_only", "canonical")
GENOME_MISMATCH = ("error", "warn", "exclude")
FORMATS = ("auto", "gff3", "gtf", "tiberius")
SCOPE_TYPES = ("all", "reference_sequence_regions", "sequences", "regions_file")
PREPARATION_STEPS = ("decompress", "retype_transcript_types")
_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")
ROLES = ("primary", "supplementary")
DISPLAY_ROLES = ("focus", "primary_comparator")
ORDER_METRICS = (
    "exact_cds_one_to_one_recall",
    "exact_cds_one_to_one_f1",
    "exact_cds_one_to_one_precision",
    "locus_recovery_cds_overlap",
)
# Display ids of earlier dashboard versions -> contract metric ids (same definitions).
ORDER_METRIC_ALIASES = {
    "exact_recall_one_to_one": "exact_cds_one_to_one_recall",
    "exact_f1_one_to_one": "exact_cds_one_to_one_f1",
    "exact_precision_one_to_one": "exact_cds_one_to_one_precision",
}
RESERVED_SCOPE_IDS = ("whole_assembly", "all")
COMPARATOR_OPTIONS = (
    "add_stop_codon",
    "stop_codon_genetic_code",
    "stop_codon_sequence_code",
)

DEFAULT_POLICIES = {
    "A": {
        "label": "Any eligible reference CDS isoform",
        "description": "A reference gene counts once when any eligible isoform is reproduced; all query transcripts are used.",
        "reference_transcript_selection": "all",
        "query_transcript_selection": "all",
    },
    "B": {
        "label": "Longest CDS per gene",
        "description": "Longest total CDS per gene on both sides (ties: first transcript in file order), chosen independently.",
        "reference_transcript_selection": "longest_cds",
        "query_transcript_selection": "longest_cds",
    },
    "C": {
        "label": "Canonical reference transcript",
        "description": "Reference transcript tagged Ensembl_canonical (falls back to longest CDS per gene without a tag); query longest CDS. Unavailable for references without canonical tags.",
        "reference_transcript_selection": "canonical",
        "query_transcript_selection": "longest_cds",
        "optional": True,
    },
}


class ConfigError(ValueError):
    """The configuration is invalid; ``problems`` lists every issue found."""

    def __init__(self, problems: list[str]):
        super().__init__(
            "Invalid experiment configuration:\n  - " + "\n  - ".join(problems)
        )
        self.problems = problems


@dataclass
class Preparation:
    steps: list[dict] = field(
        default_factory=list
    )  # [{"name": "decompress"}, {"name": "retype_transcript_types", "types": [...]}]
    reuse_existing: Path | None = None

    def recipe(self) -> list[dict]:
        return [dict(step) for step in self.steps]


@dataclass
class Reference:
    id: str
    species: str
    assembly: str
    annotation: Path
    genome: Path | None
    assembly_accession: str | None = None
    annotation_source: str | None = None
    annotation_release: str | None = None
    annotation_format: str = "auto"
    annotation_preparation: Preparation = field(default_factory=Preparation)
    genome_preparation: Preparation = field(default_factory=Preparation)
    scope: dict = field(default_factory=lambda: {"type": "all"})
    seqname_map: Path | None = None
    notes: str | None = None
    analysis_scopes: list[dict] = field(default_factory=list)
    raw: dict = field(default_factory=dict)


@dataclass
class Run:
    run_id: str
    reference: str
    tool: str
    annotation: Path
    tool_version: str | None = None
    model: str | None = None
    format: str = "auto"
    preparation: Preparation = field(default_factory=Preparation)
    provenance: str | None = None
    caveats: list[str] = field(default_factory=list)
    aliases: list[str] = field(default_factory=list)
    color: str | None = None
    synthetic: bool = False
    label: str | None = None  # display label; defaults to tool (+ run_label)
    role: str = "primary"  # or "supplementary"
    comparator_options: dict = field(default_factory=dict)  # normalised, see below
    companion_of: str | None = None  # presentation: listed after this run
    run_label: str | None = None  # e.g. "GMB run 2": several runs of one tool
    display_role: str | None = None  # presentation: focus / primary_comparator
    source: dict = field(default_factory=dict)  # original_path, sha256, recorded_in
    legacy_ids: dict = field(default_factory=dict)
    input_audit: Path | None = None
    raw: dict = field(default_factory=dict)

    @property
    def display_label(self) -> str:
        if self.label:
            return self.label
        if self.run_label:
            if self.tool.lower() in self.run_label.lower():
                return self.run_label
            return f"{self.tool} · {self.run_label}"
        return self.tool

    @property
    def variant(self) -> str:
        """Coordinate treatment of the query: original or stop_harmonised."""
        return "stop_harmonised" if self.comparator_options.get("add_stop_codon") else "original"


@dataclass
class Policy:
    id: str
    label: str
    description: str
    reference_transcript_selection: str
    query_transcript_selection: str
    optional: bool = False
    plots_per_category: int = 0
    role: str = "primary"  # or "supplementary"


@dataclass
class Experiment:
    path: Path
    base_dir: Path
    id: str
    title: str
    description: str
    synthetic: bool
    output_dir: Path
    evaluation: dict
    policies: dict[str, Policy]
    references: dict[str, Reference]
    runs: dict[str, Run]
    imports: dict
    raw: dict
    display: dict = field(default_factory=dict)
    registries: list[Path] = field(default_factory=list)
    config_version: int = CONFIG_VERSION

    def runs_for(self, reference_id: str) -> list[Run]:
        return [run for run in self.runs.values() if run.reference == reference_id]

    def to_public(self) -> dict:
        """Configuration as written (paths as strings) for provenance and the UI."""
        return self.raw


def _load_text(path: Path) -> dict:
    text = path.read_text()
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml  # optional: pip install "ensembl-genes[dashboard]"
        except ImportError as error:
            raise ConfigError(
                [
                    f"{path}: YAML configuration needs PyYAML (pip install 'ensembl-genes[dashboard]'), or use JSON"
                ]
            ) from error
        return yaml.safe_load(text)
    return json.loads(text)


def _path(
    base: Path, value, problems: list[str], what: str, required=True
) -> Path | None:
    if value in (None, ""):
        if required:
            problems.append(f"{what}: path is required")
        return None
    if not isinstance(value, str):
        problems.append(f"{what}: path must be a string")
        return None
    path = Path(value).expanduser()
    return (base / path).resolve() if not path.is_absolute() else path


def _preparation(base: Path, value, problems: list[str], what: str) -> Preparation:
    if not value:
        return Preparation()
    if not isinstance(value, dict):
        problems.append(f"{what}: must be an object with 'steps'")
        return Preparation()
    steps = []
    for raw in value.get("steps", []):
        step = {"name": raw} if isinstance(raw, str) else dict(raw)
        name = step.get("name")
        if name not in PREPARATION_STEPS:
            problems.append(
                f"{what}: unknown preparation step {name!r} (allowed: {', '.join(PREPARATION_STEPS)})"
            )
            continue
        if name == "retype_transcript_types":
            types = step.get("types")
            if not types or not all(isinstance(t, str) for t in types):
                problems.append(
                    f"{what}: retype_transcript_types needs a non-empty 'types' list"
                )
            step["to"] = step.get("to", "transcript")
        steps.append(step)
    reuse = _path(
        base,
        value.get("reuse_existing"),
        problems,
        f"{what}.reuse_existing",
        required=False,
    )
    return Preparation(steps=steps, reuse_existing=reuse)


def stop_codon_options(add_stop_codon, genetic_code=None, sequence_code=None) -> dict:
    """
    Canonical stop-codon parameters, shared by configurations and manifests.

    Returns {} when stop-codon harmonisation is off, so original-coordinate
    results keep the parameters (and cache keys) they had before the option
    existed. Otherwise the table and per-sequence overrides as a sorted
    ``SEQNAME=TABLE`` list ("none" = no table; that sequence is skipped).
    Raises:
            ValueError for unsupported tables or malformed overrides
    """
    from ensembl.genes.annotation_qc.metrics.pairwise.stop_codons import (
        parse_genetic_code_policy,
    )

    if not add_stop_codon:
        return {}
    if isinstance(sequence_code, dict):
        sequence_code = [f"{name}={table}" for name, table in sequence_code.items()]
    policy = parse_genetic_code_policy(
        1 if genetic_code is None else genetic_code, list(sequence_code or [])
    )
    return {
        "add_stop_codon": True,
        "stop_codon_genetic_code": policy.default,
        "stop_codon_sequence_code": sorted(
            f"{name}={'none' if table is None else table}"
            for name, table in policy.overrides.items()
        )
        or None,
    }


def _comparator_options(spec, problems: list[str], what: str) -> dict:
    spec = spec or {}
    if not isinstance(spec, dict):
        problems.append(f"{what}: must be an object")
        return {}
    unknown = sorted(set(spec) - set(COMPARATOR_OPTIONS))
    if unknown:
        problems.append(
            f"{what}: unknown option(s) {unknown} (allowed: {', '.join(COMPARATOR_OPTIONS)})"
        )
    add = spec.get("add_stop_codon", False)
    if not isinstance(add, bool):
        problems.append(f"{what}.add_stop_codon must be true or false")
        return {}
    if not add and (
        spec.get("stop_codon_genetic_code") is not None
        or spec.get("stop_codon_sequence_code")
    ):
        problems.append(
            f"{what}: stop_codon_genetic_code / stop_codon_sequence_code need add_stop_codon: true"
        )
    try:
        return stop_codon_options(
            add,
            spec.get("stop_codon_genetic_code"),
            spec.get("stop_codon_sequence_code"),
        )
    except ValueError as error:
        problems.append(f"{what}: {error}")
        return {}


def _analysis_scopes(spec, problems: list[str], what: str) -> list[dict]:
    scopes, seen = [], set()
    for i, raw in enumerate(spec or []):
        where = f"{what}[{i}]"
        if not isinstance(raw, dict):
            problems.append(f"{where}: must be an object")
            continue
        sid = raw.get("id", "")
        if not _ID.match(str(sid)) or sid in RESERVED_SCOPE_IDS:
            problems.append(
                f"{where}.id must be an identifier other than {RESERVED_SCOPE_IDS}"
            )
            continue
        if sid in seen:
            problems.append(f"{where}: duplicate scope id {sid!r}")
            continue
        seen.add(sid)
        seqs = raw.get("sequences")
        if (
            not isinstance(seqs, list)
            or not seqs
            or not all(isinstance(x, str) and x for x in seqs)
        ):
            problems.append(f"{where}.sequences must be a non-empty list of sequence names")
            continue
        if len(set(seqs)) != len(seqs):
            problems.append(f"{where}.sequences lists a sequence twice")
        scopes.append(
            {
                "id": sid,
                "label": raw.get("label", sid),
                "sequences": list(seqs),
                "description": raw.get("description", ""),
                "role": "supplementary",
            }
        )
    return scopes


def _display(spec, runs: dict, problems: list[str]) -> dict:
    spec = spec or {}
    if not isinstance(spec, dict):
        problems.append("display: must be an object")
        return {}
    pinned = spec.get("pinned_runs") or []
    for run_id in pinned:
        if run_id not in runs:
            problems.append(f"display.pinned_runs: {run_id!r} is not a configured run_id")
    tools = spec.get("preferred_tools") or []
    if not isinstance(tools, list) or not all(isinstance(t, str) for t in tools):
        problems.append("display.preferred_tools must be a list of tool names")
        tools = []
    known = {run.tool for run in runs.values()}
    for tool in tools:
        if tool not in known:
            problems.append(f"display.preferred_tools: {tool!r} is not the tool of any run")
    order = spec.get("order_metric", "exact_cds_one_to_one_recall")
    order = ORDER_METRIC_ALIASES.get(order, order)
    if order not in ORDER_METRICS:
        problems.append(f"display.order_metric must be one of {ORDER_METRICS}")
    return {"pinned_runs": list(pinned), "preferred_tools": list(tools), "order_metric": order}


def load_experiment(path: str | Path) -> Experiment:
    """
    Load and statically validate an experiment configuration.

    File existence and content checks belong to ``validation.validate_experiment``.
    Raises:
            ConfigError listing every structural problem
    """
    path = Path(path).resolve()
    try:
        raw = _load_text(path)
    except (OSError, json.JSONDecodeError) as error:
        raise ConfigError([f"{path}: cannot read configuration: {error}"]) from error
    problems: list[str] = []
    base = path.parent
    if raw.get("config_version") not in CONFIG_VERSIONS:
        problems.append(
            f"config_version must be one of {CONFIG_VERSIONS} (got {raw.get('config_version')!r})"
        )

    meta = raw.get("experiment") or {}
    exp_id = meta.get("id", "")
    if not _ID.match(str(exp_id)):
        problems.append("experiment.id is required (letters, digits, '_', '-', '.')")
    output_dir = _path(base, raw.get("output_dir"), problems, "output_dir")

    evaluation = {
        "evaluation_mode": "protein_coding",
        "reference_gene_biotypes": [],
        "reference_transcript_biotypes": [],
        "genome_mismatch": "error",
        **(raw.get("evaluation") or {}),
    }
    if evaluation["evaluation_mode"] not in EVALUATION_MODES:
        problems.append(f"evaluation.evaluation_mode must be one of {EVALUATION_MODES}")
    if evaluation["genome_mismatch"] not in GENOME_MISMATCH:
        problems.append(f"evaluation.genome_mismatch must be one of {GENOME_MISMATCH}")
    for key in ("reference_gene_biotypes", "reference_transcript_biotypes"):
        if not isinstance(evaluation[key], list):
            problems.append(f"evaluation.{key} must be a list")

    policies = {}
    for pid, spec in (raw.get("policies") or DEFAULT_POLICIES).items():
        if pid.startswith("_"):  # comments
            continue
        if not _ID.match(pid):
            problems.append(f"policy id {pid!r} is not a valid identifier")
            continue
        spec = {**DEFAULT_POLICIES.get(pid, {}), **(spec or {})}
        for key in ("reference_transcript_selection", "query_transcript_selection"):
            if spec.get(key) not in SELECTIONS:
                problems.append(f"policies.{pid}.{key} must be one of {SELECTIONS}")
        policies[pid] = Policy(
            id=pid,
            label=spec.get("label", pid),
            description=spec.get("description", ""),
            reference_transcript_selection=spec.get(
                "reference_transcript_selection", "all"
            ),
            query_transcript_selection=spec.get("query_transcript_selection", "all"),
            optional=bool(
                spec.get(
                    "optional",
                    spec.get("reference_transcript_selection") == "canonical",
                )
            ),
            plots_per_category=int(spec.get("plots_per_category", 0)),
            role=spec.get("role", "primary"),
        )
        if policies[pid].role not in ROLES:
            problems.append(f"policies.{pid}.role must be one of {ROLES}")

    references = {}
    for i, spec in enumerate(raw.get("references") or []):
        what = f"references[{i}]"
        rid = spec.get("id", "")
        if not _ID.match(str(rid)):
            problems.append(f"{what}.id is required and must be an identifier")
            continue
        if rid in references:
            problems.append(f"duplicate reference id {rid!r}")
            continue
        scope = spec.get("scope") or {"type": "all"}
        if scope.get("type") not in SCOPE_TYPES:
            problems.append(f"{what}.scope.type must be one of {SCOPE_TYPES}")
        if scope.get("type") == "sequences" and not scope.get("sequences"):
            problems.append(f"{what}.scope.sequences must list sequence names")
        if scope.get("type") == "regions_file":
            scope = {
                **scope,
                "path": str(
                    _path(base, scope.get("path"), problems, f"{what}.scope.path")
                ),
            }
        fmt = spec.get("annotation_format", "auto")
        if fmt not in FORMATS:
            problems.append(f"{what}.annotation_format must be one of {FORMATS}")
        references[rid] = Reference(
            id=rid,
            species=spec.get("species")
            or problems.append(f"{what}.species is required")
            or "",
            assembly=spec.get("assembly")
            or problems.append(f"{what}.assembly is required")
            or "",
            annotation=_path(
                base, spec.get("annotation"), problems, f"{what}.annotation"
            ),
            genome=_path(
                base, spec.get("genome"), problems, f"{what}.genome", required=False
            ),
            assembly_accession=spec.get("assembly_accession"),
            annotation_source=spec.get("annotation_source"),
            annotation_release=spec.get("annotation_release"),
            annotation_format=fmt,
            annotation_preparation=_preparation(
                base,
                spec.get("annotation_preparation"),
                problems,
                f"{what}.annotation_preparation",
            ),
            genome_preparation=_preparation(
                base,
                spec.get("genome_preparation"),
                problems,
                f"{what}.genome_preparation",
            ),
            scope=scope,
            seqname_map=_path(
                base,
                spec.get("seqname_map"),
                problems,
                f"{what}.seqname_map",
                required=False,
            ),
            notes=spec.get("notes"),
            analysis_scopes=_analysis_scopes(
                spec.get("analysis_scopes"), problems, f"{what}.analysis_scopes"
            ),
            raw=spec,
        )
    if not references:
        problems.append("at least one reference is required")

    runs = {}
    for i, spec in enumerate(raw.get("runs") or []):
        what = f"runs[{i}]"
        run_id = spec.get("run_id", "")
        if not _ID.match(str(run_id)):
            problems.append(f"{what}.run_id is required and must be an identifier")
            continue
        if run_id in runs:
            problems.append(
                f"duplicate run_id {run_id!r} (tool name is not an identifier; give each run its own run_id)"
            )
            continue
        if spec.get("reference") not in references:
            problems.append(
                f"{what}.reference {spec.get('reference')!r} is not a configured reference id"
            )
        if not spec.get("tool"):
            problems.append(f"{what}.tool is required")
        fmt = spec.get("format", "auto")
        if fmt not in FORMATS:
            problems.append(f"{what}.format must be one of {FORMATS}")
        options = _comparator_options(
            spec.get("comparator_options"), problems, f"{what}.comparator_options"
        )
        reference = references.get(spec.get("reference"))
        if options.get("add_stop_codon") and reference and reference.genome is None:
            problems.append(
                f"{what}: add_stop_codon needs a genome on reference {reference.id}"
            )
        role = spec.get("role", "primary")
        if role not in ROLES:
            problems.append(f"{what}.role must be one of {ROLES}")
        display_role = spec.get("display_role")
        if display_role is not None and display_role not in DISPLAY_ROLES:
            problems.append(f"{what}.display_role must be one of {DISPLAY_ROLES}")
        source = spec.get("source") or {}
        if not isinstance(source, dict) or (
            source.get("sha256") and not re.fullmatch(r"[0-9a-f]{64}", str(source["sha256"]))
        ):
            problems.append(f"{what}.source must be an object; sha256 must be 64 hex digits")
            source = {}
        legacy_ids = spec.get("legacy_ids") or {}
        if not isinstance(legacy_ids, dict):
            problems.append(f"{what}.legacy_ids must be an object")
            legacy_ids = {}
        runs[run_id] = Run(
            run_id=run_id,
            reference=spec.get("reference", ""),
            tool=spec.get("tool", ""),
            annotation=_path(
                base, spec.get("annotation"), problems, f"{what}.annotation"
            ),
            tool_version=spec.get("tool_version"),
            model=spec.get("model"),
            format=fmt,
            preparation=_preparation(
                base, spec.get("preparation"), problems, f"{what}.preparation"
            ),
            provenance=spec.get("provenance"),
            caveats=list(spec.get("caveats") or []),
            aliases=list(spec.get("aliases") or []),
            color=spec.get("color"),
            synthetic=bool(spec.get("synthetic", meta.get("synthetic", False))),
            label=spec.get("label"),
            role=role,
            comparator_options=options,
            companion_of=spec.get("companion_of"),
            run_label=spec.get("run_label"),
            display_role=display_role,
            source=dict(source),
            legacy_ids=dict(legacy_ids),
            input_audit=_path(
                base, spec.get("input_audit"), problems, f"{what}.input_audit", required=False
            ),
            raw=spec,
        )
    _check_run_consistency(runs, problems)
    for run in runs.values():
        if run.companion_of is None:
            continue
        other = runs.get(run.companion_of)
        if other is None or other.reference != run.reference or other is run:
            problems.append(
                f"run {run.run_id}: companion_of must name another run on the same reference"
            )
        elif run.role != "supplementary":
            problems.append(f"run {run.run_id}: companion_of needs role 'supplementary'")

    imports = dict(raw.get("import") or {})
    if imports.get("benchmark_dir"):
        imports["benchmark_dir"] = _path(
            base, imports["benchmark_dir"], problems, "import.benchmark_dir"
        )
    for key in ("results_dirs", "pilot_dirs"):
        values = imports.get(key) or []
        if not isinstance(values, list):
            problems.append(f"import.{key} must be a list of directories")
            values = []
        imports[key] = [_path(base, p, problems, f"import.{key}") for p in values]
    registries = [
        _path(base, p, problems, "provenance.registries")
        for p in ((raw.get("provenance") or {}).get("registries") or [])
    ]
    for key in ("code_snapshot", "context_documents", "paper_analogous"):
        imports[key] = [
            _path(base, p, problems, f"import.{key}") for p in imports.get(key) or []
        ]

    display = _display(raw.get("display"), runs, problems)

    if problems:
        raise ConfigError(problems)
    return Experiment(
        path=path,
        base_dir=base,
        id=exp_id,
        title=meta.get("title", exp_id),
        description=meta.get("description", ""),
        synthetic=bool(meta.get("synthetic", False)),
        output_dir=output_dir,
        evaluation=evaluation,
        policies=policies,
        references=references,
        runs=runs,
        imports=imports,
        raw=raw,
        display=display,
        registries=[p for p in registries if p],
        config_version=raw.get("config_version", CONFIG_VERSION),
    )


def _check_run_consistency(runs: dict, problems: list[str]) -> None:
    """Settings that would make two runs indistinguishable or a companion meaningless."""
    seen: dict[tuple, str] = {}
    for run in runs.values():
        key = (
            run.reference,
            str(run.annotation),
            json.dumps(run.comparator_options, sort_keys=True),
        )
        if key in seen:
            problems.append(
                f"runs {seen[key]} and {run.run_id} compare the same annotation with the same "
                "settings on the same reference; their results could not be told apart"
            )
        else:
            seen[key] = run.run_id
        other = runs.get(run.companion_of) if run.companion_of else None
        if other is not None and other.annotation != run.annotation:
            problems.append(
                f"run {run.run_id}: companion_of {other.run_id} must compare the same annotation "
                "file (a companion is another treatment of the same predictions)"
            )
        if other is not None and other.comparator_options == run.comparator_options:
            problems.append(
                f"run {run.run_id}: companion_of {other.run_id} has identical comparator options"
            )


def analysis_scope_parameters(reference: Reference) -> dict:
    """Configured analysis scopes as the comparator records them (sorted sequences)."""
    return {
        scope["id"]: sorted(scope["sequences"]) for scope in reference.analysis_scopes
    }


def comparator_parameters(
    experiment: Experiment, reference: Reference, run: Run, policy: Policy
) -> dict:
    """
    Effective pairwise-compare parameters for one (run, policy); used in cache keys.

    Stop-codon keys appear only when harmonisation is on (see stop_codon_options),
    so the parameters of original-coordinate runs are unchanged.
    """
    ev = experiment.evaluation
    return {
        "evaluation_mode": ev["evaluation_mode"],
        "reference_gene_biotypes": sorted(ev["reference_gene_biotypes"]) or None,
        "reference_transcript_biotypes": sorted(ev["reference_transcript_biotypes"])
        or None,
        "reference_transcript_selection": policy.reference_transcript_selection,
        "query_transcript_selection": policy.query_transcript_selection,
        "genome_mismatch": ev["genome_mismatch"],
        "query_format": run.format,
        "reference_format": reference.annotation_format,
        "plots_per_category": policy.plots_per_category,
        "region": None,
        "assembly_report": None,
        "analysis_scopes": analysis_scope_parameters(reference),
        **run.comparator_options,
    }

"""
Experiment configuration (format version 1).

A configuration is a JSON file (or YAML, when PyYAML is installed through the
``dashboard`` extra). Relative paths are resolved against the directory that
contains the configuration file. See ``templates/experiment_template.json`` and
the README for every key.

Identity rules:
    reference ``id`` and run ``run_id`` are the unique identifiers. ``tool`` is a
    display/grouping label only, so several runs (models, versions, settings) of
    the same tool can coexist. Policy keys (e.g. "A", "B", "C") are identifiers too.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

CONFIG_VERSION = 1
SELECTIONS = ("all", "longest_cds", "canonical")
EVALUATION_MODES = ("all", "protein_coding", "cds_only", "canonical")
GENOME_MISMATCH = ("error", "warn", "exclude")
FORMATS = ("auto", "gff3", "gtf", "tiberius")
SCOPE_TYPES = ("all", "reference_sequence_regions", "sequences", "regions_file")
PREPARATION_STEPS = ("decompress", "retype_transcript_types")
_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")

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
    raw: dict = field(default_factory=dict)


@dataclass
class Policy:
    id: str
    label: str
    description: str
    reference_transcript_selection: str
    query_transcript_selection: str
    optional: bool = False
    plots_per_category: int = 0


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
    if raw.get("config_version") != CONFIG_VERSION:
        problems.append(
            f"config_version must be {CONFIG_VERSION} (got {raw.get('config_version')!r})"
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
        )

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
            raw=spec,
        )

    imports = dict(raw.get("import") or {})
    if imports.get("benchmark_dir"):
        imports["benchmark_dir"] = _path(
            base, imports["benchmark_dir"], problems, "import.benchmark_dir"
        )
    for key in ("code_snapshot", "context_documents", "paper_analogous"):
        imports[key] = [
            _path(base, p, problems, f"import.{key}") for p in imports.get(key) or []
        ]

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
    )


def comparator_parameters(
    experiment: Experiment, reference: Reference, run: Run, policy: Policy
) -> dict:
    """Effective pairwise-compare parameters for one (run, policy); used in cache keys."""
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
    }

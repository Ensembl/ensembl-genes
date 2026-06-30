#!/usr/bin/env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
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
"""Compare old metadata output with registry-backed metadata output.

This script is intended for eHive checking analyses. It exits with:
  0 when the generated metadata JSON values match
  1 when both scripts ran but generated different metadata values
  2 when either script could not be run or its JSON could not be parsed
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

SQL_STRING = r"'(?P<{name}>(?:''|[^'])*)'"
INSERT_RE = re.compile(
    r"^INSERT\s+IGNORE\s+INTO\s+meta\s*"
    r"\(\s*species_id\s*,\s*meta_key\s*,\s*meta_value\s*\)\s*"
    r"VALUES\s*\(\s*(?P<species_id>\d+)\s*,\s*"
    + SQL_STRING.format(name="meta_key")
    + r"\s*,\s*"
    + SQL_STRING.format(name="meta_value")
    + r"\s*\)\s*;\s*$",
    re.IGNORECASE,
)
UPDATE_RE = re.compile(
    r"^UPDATE\s+meta\s+SET\s+meta_value\s*=\s*"
    + SQL_STRING.format(name="meta_value")
    + r"\s+WHERE\s+meta_key\s*=\s*"
    + SQL_STRING.format(name="meta_key")
    + r"\s*;\s*$",
    re.IGNORECASE,
)
DELETE_RE = re.compile(
    r"^DELETE\s+from\s+meta\s+WHERE\s+meta_key\s*=\s*"
    + SQL_STRING.format(name="meta_key")
    + r"\s*;\s*$",
    re.IGNORECASE,
)
USE_RE = re.compile(r"^USE\s+\S+\s*;\s*$", re.IGNORECASE)
# method and start_date are set by the pipeline and so it is not expected that the old script see these values
# assembly.date used to come from ENA_LAST_UPDATE but the registry gets assembly date from NCBI directly which we
# consider as authoritative. it is therefore ignored if there is ambiguity.
# These keys are still logged for record keeping but they do not cause the script to fail.
IGNORED_METADATA_KEYS = {"genebuild.method", "genebuild.start_date", "assembly.date"}


class TruthCheckError(Exception):
    """Base error for expected checker failures."""

    exit_code = 2

    def __init__(self, message: str, report: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.report = report or {}


class TruthMismatchError(TruthCheckError):
    """Raised when old and registry metadata outputs differ."""

    exit_code = 1


@dataclass(frozen=True)
class ScriptRun:
    """Result paths and process output for one metadata script run."""

    name: str
    command: List[str]
    output_dir: Path
    sql_path: Path
    json_path: Optional[Path]
    stdout: str
    stderr: str


def sql_unescape(value: str) -> str:
    """Unescape SQL single-quoted string content."""
    return value.replace("''", "'")


def output_name(db_name: str, production_name: Optional[str]) -> str:
    """Return the metadata script output basename."""
    return production_name or db_name


def tail_text(text: str, max_lines: int = 40) -> str:
    """Return a bounded process-output tail for readable error reports."""
    lines = text.splitlines()
    return "\n".join(lines[-max_lines:])


def parse_sql_patch(sql_path: Path) -> List[Dict[str, str]]:
    """Parse generated SQL into normalized metadata actions."""
    if not sql_path.exists():
        raise TruthCheckError(
            f"Expected SQL output was not created: {sql_path}",
            {"sql_path": str(sql_path)},
        )

    actions: List[Dict[str, str]] = []
    with open(sql_path, "r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or USE_RE.match(line):
                continue

            insert_match = INSERT_RE.match(line)
            if insert_match:
                actions.append(
                    {
                        "action": "insert",
                        "meta_key": sql_unescape(insert_match.group("meta_key")),
                        "meta_value": sql_unescape(insert_match.group("meta_value")),
                        "species_id": insert_match.group("species_id"),
                    }
                )
                continue

            update_match = UPDATE_RE.match(line)
            if update_match:
                actions.append(
                    {
                        "action": "update",
                        "meta_key": sql_unescape(update_match.group("meta_key")),
                        "meta_value": sql_unescape(update_match.group("meta_value")),
                    }
                )
                continue

            delete_match = DELETE_RE.match(line)
            if delete_match:
                actions.append(
                    {
                        "action": "delete",
                        "meta_key": sql_unescape(delete_match.group("meta_key")),
                        "meta_value": "",
                    }
                )
                continue

            raise TruthCheckError(
                f"Could not parse SQL generated by {sql_path} at line {line_number}",
                {
                    "sql_path": str(sql_path),
                    "line_number": line_number,
                    "line": line,
                },
            )

    return sorted(
        actions,
        key=lambda item: (
            item["meta_key"],
            item["action"],
            item.get("meta_value", ""),
            item.get("species_id", ""),
        ),
    )


def actions_by_key(
    actions: Iterable[Dict[str, str]]
) -> Dict[str, List[Dict[str, str]]]:
    """Group SQL actions by meta_key for clearer diffs."""
    grouped: Dict[str, List[Dict[str, str]]] = {}
    for action in actions:
        grouped.setdefault(action["meta_key"], []).append(action)
    return grouped


def compare_actions(
    old_actions: List[Dict[str, str]], registry_actions: List[Dict[str, str]]
) -> Dict[str, Any]:
    """Return a structured difference report."""
    old_by_key = actions_by_key(old_actions)
    registry_by_key = actions_by_key(registry_actions)
    all_keys = sorted(set(old_by_key) | set(registry_by_key))

    missing_from_registry = []
    extra_in_registry = []
    mismatched = []

    for meta_key in all_keys:
        old_entries = old_by_key.get(meta_key, [])
        registry_entries = registry_by_key.get(meta_key, [])
        if old_entries == registry_entries:
            continue
        if old_entries and not registry_entries:
            missing_from_registry.append(
                {"meta_key": meta_key, "old": old_entries, "registry": []}
            )
        elif registry_entries and not old_entries:
            extra_in_registry.append(
                {"meta_key": meta_key, "old": [], "registry": registry_entries}
            )
        else:
            mismatched.append(
                {
                    "meta_key": meta_key,
                    "old": old_entries,
                    "registry": registry_entries,
                }
            )

    return {
        "missing_from_registry": missing_from_registry,
        "extra_in_registry": extra_in_registry,
        "mismatched": mismatched,
        "different_key_count": (
            len(missing_from_registry) + len(extra_in_registry) + len(mismatched)
        ),
    }


def load_metadata_json(json_path: Optional[Path]) -> Dict[str, Any]:
    """Load generated metadata JSON and return its metadata dictionary."""
    if json_path is None:
        raise TruthCheckError("Expected JSON output path was not configured")
    if not json_path.exists():
        raise TruthCheckError(
            f"Expected JSON output was not created: {json_path}",
            {"json_path": str(json_path)},
        )

    try:
        with open(json_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        raise TruthCheckError(
            f"Could not parse JSON generated by {json_path}",
            {
                "json_path": str(json_path),
                "line_number": exc.lineno,
                "column": exc.colno,
                "message": exc.msg,
            },
        ) from exc

    metadata = payload.get("metadata")
    if not isinstance(metadata, dict):
        raise TruthCheckError(
            f"JSON generated by {json_path} does not contain a metadata object",
            {"json_path": str(json_path)},
        )

    return {key: "" if value is None else str(value) for key, value in metadata.items()}


def compare_metadata(
    old_metadata: Dict[str, str], registry_metadata: Dict[str, str]
) -> Dict[str, Any]:
    """Return a structured difference report for generated metadata values."""
    all_keys = sorted(set(old_metadata) | set(registry_metadata))

    missing_from_registry = []
    extra_in_registry = []
    mismatched = []
    ignored = []

    for meta_key in all_keys:
        old_has_key = meta_key in old_metadata
        registry_has_key = meta_key in registry_metadata
        old_value = old_metadata.get(meta_key)
        registry_value = registry_metadata.get(meta_key)

        if meta_key in IGNORED_METADATA_KEYS:
            if old_has_key != registry_has_key or old_value != registry_value:
                ignored.append(
                    {
                        "meta_key": meta_key,
                        "old": old_value if old_has_key else None,
                        "registry": registry_value if registry_has_key else None,
                    }
                )
            continue

        if old_has_key and registry_has_key and old_value == registry_value:
            continue
        if old_has_key and not registry_has_key:
            missing_from_registry.append(
                {"meta_key": meta_key, "old": old_value, "registry": None}
            )
        elif registry_has_key and not old_has_key:
            extra_in_registry.append(
                {"meta_key": meta_key, "old": None, "registry": registry_value}
            )
        else:
            mismatched.append(
                {
                    "meta_key": meta_key,
                    "old": old_value,
                    "registry": registry_value,
                }
            )

    return {
        "missing_from_registry": missing_from_registry,
        "extra_in_registry": extra_in_registry,
        "mismatched": mismatched,
        "ignored": ignored,
        "ignored_keys": sorted(IGNORED_METADATA_KEYS),
        "different_key_count": (
            len(missing_from_registry) + len(extra_in_registry) + len(mismatched)
        ),
    }


def script_env(repo_root: Path) -> Dict[str, str]:
    """Build an environment that can import local ensembl modules."""
    env = os.environ.copy()
    src_python = str(repo_root / "src" / "python")
    existing_pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        src_python
        if not existing_pythonpath
        else src_python + os.pathsep + existing_pythonpath
    )
    return env


def run_command(
    name: str,
    command: List[str],
    output_dir: Path,
    sql_path: Path,
    json_path: Optional[Path],
    repo_root: Path,
) -> ScriptRun:
    """Run a metadata script and convert failures to eHive-friendly errors."""
    process = subprocess.run(
        command,
        cwd=repo_root,
        env=script_env(repo_root),
        check=False,
        capture_output=True,
        text=True,
    )
    run = ScriptRun(
        name=name,
        command=command,
        output_dir=output_dir,
        sql_path=sql_path,
        json_path=json_path,
        stdout=process.stdout,
        stderr=process.stderr,
    )
    if process.returncode != 0:
        raise TruthCheckError(
            f"{name} metadata script failed with exit code {process.returncode}",
            {
                "script": name,
                "exit_code": process.returncode,
                "command": command,
                "stdout_tail": tail_text(process.stdout),
                "stderr_tail": tail_text(process.stderr),
            },
        )
    return run


def build_old_command(
    args: argparse.Namespace, output_dir: Path, json_path: Path
) -> List[str]:
    """Build the old core_meta_data.py command."""
    command = [
        sys.executable,
        str(args.old_script),
        "-d",
        args.db_name,
        "-s",
        args.host,
        "-p",
        str(args.port),
        "-t",
        args.team,
        "-o",
        str(output_dir),
        "--json_output",
        str(json_path),
    ]
    if args.production_name:
        command.extend(["-n", args.production_name])
    if args.verbose_scripts:
        command.append("-v")
    return command


def build_registry_command(
    args: argparse.Namespace, output_dir: Path, json_path: Path
) -> List[str]:
    """Build the registry-backed core_meta_reg.py command."""
    command = [
        sys.executable,
        str(args.registry_script),
        "-d",
        args.db_name,
        "-s",
        args.host,
        "-p",
        str(args.port),
        "-t",
        args.team,
        "-o",
        str(output_dir),
        "--json_output",
        str(json_path),
    ]
    if args.production_name:
        command.extend(["-n", args.production_name])
    if args.assembly_accession:
        command.extend(["--assembly_accession", args.assembly_accession])
    if args.registry_host:
        command.extend(["--registry_host", args.registry_host])
    if args.registry_port:
        command.extend(["--registry_port", str(args.registry_port)])
    if args.registry_user:
        command.extend(["--registry_user", args.registry_user])
    if args.registry_password:
        command.extend(["--registry_password", args.registry_password])
    if args.registry_db:
        command.extend(["--registry_db", args.registry_db])
    if args.genebuilder:
        command.extend(["--genebuilder", args.genebuilder])
    if args.db_user:
        command.extend(["--db_user", args.db_user])
    if args.db_password:
        command.extend(["--db_password", args.db_password])
    if args.verbose_scripts:
        command.append("-v")
    return command


def write_report(path: Optional[Path], report: Dict[str, Any]) -> None:
    """Write a JSON report if requested."""
    if not path:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


def resolve_work_dir(
    args: argparse.Namespace,
) -> tempfile.TemporaryDirectory[str] | None:
    """Create a tempdir unless the user provided persistent output."""
    if args.work_dir:
        Path(args.work_dir).mkdir(parents=True, exist_ok=True)
        return None
    if args.keep_outputs:
        raise TruthCheckError("--keep_outputs requires --work_dir")
    try:
        return tempfile.TemporaryDirectory(prefix="metadata_truth_check_")
    except OSError as exc:
        raise TruthCheckError(
            "Could not create a temporary directory. Provide --work_dir pointing "
            "to a writable location for eHive job output.",
            {"error": str(exc)},
        ) from exc


def run_truth_check(args: argparse.Namespace) -> Dict[str, Any]:
    """Run both scripts, parse JSON outputs, and compare generated metadata."""
    repo_root = args.repo_root
    temp_dir = resolve_work_dir(args)
    try:
        base_dir = Path(args.work_dir) if args.work_dir else Path(temp_dir.name)
        base_dir.mkdir(parents=True, exist_ok=True)

        basename = output_name(args.db_name, args.production_name)
        old_sql = base_dir / f"{basename}.sql"
        old_json = base_dir / f"{basename}.json"
        registry_sql = base_dir / f"{basename}_registry.sql"
        registry_json = base_dir / f"{basename}_registry.json"

        old_run = run_command(
            "old",
            build_old_command(args, base_dir, old_json),
            base_dir,
            old_sql,
            old_json,
            repo_root,
        )
        registry_run = run_command(
            "registry",
            build_registry_command(args, base_dir, registry_json),
            base_dir,
            registry_sql,
            registry_json,
            repo_root,
        )

        old_metadata = load_metadata_json(old_run.json_path)
        registry_metadata = load_metadata_json(registry_run.json_path)
        differences = compare_metadata(old_metadata, registry_metadata)

        report = {
            "status": "ok" if differences["different_key_count"] == 0 else "failed",
            "database": args.db_name,
            "production_name": args.production_name,
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "old_sql": str(old_sql),
            "old_json": str(old_json),
            "registry_sql": str(registry_sql),
            "registry_json": str(registry_json),
            "old_metadata_count": len(old_metadata),
            "registry_metadata_count": len(registry_metadata),
            "differences": differences,
        }
        write_report(args.report_json, report)

        if differences["different_key_count"]:
            raise TruthMismatchError(
                "Metadata truth check failed: old and registry JSON metadata differ",
                report,
            )

        return report
    finally:
        if temp_dir is not None:
            temp_dir.cleanup()


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the command line parser."""
    repo_root = Path(__file__).resolve().parents[5]
    metadata_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Compare old core metadata output with registry-backed output"
    )
    parser.add_argument("-d", "--db_name", required=True, help="Core database name")
    parser.add_argument("-s", "--host", required=True, help="Core database host")
    parser.add_argument("-p", "--port", required=True, type=int, help="Core port")
    parser.add_argument(
        "-n", "--production_name", help="species.production_name for collection DBs"
    )
    parser.add_argument(
        "-t",
        "--team",
        required=True,
        type=lambda value: value.capitalize(),
        help="Team responsible for the database",
    )
    parser.add_argument(
        "-a",
        "--assembly_accession",
        help="Assembly accession override for the registry-backed script",
    )
    parser.add_argument("--registry_host", default=os.environ.get("GBS1"))
    parser.add_argument("--registry_port", default=os.environ.get("GBP1"), type=int)
    parser.add_argument("--registry_user", default="ensro")
    parser.add_argument("--registry_password", default="")
    parser.add_argument("--registry_db", default="gb_assembly_metadata")
    parser.add_argument("--genebuilder")
    parser.add_argument(
        "--db_user",
        default="ensro",
        help="Core DB user for registry script. The old script always uses ensro.",
    )
    parser.add_argument("--db_password", default="")
    parser.add_argument(
        "--old_script",
        type=Path,
        default=metadata_dir / "core_meta_data.py",
        help="Path to old metadata script",
    )
    parser.add_argument(
        "--registry_script",
        type=Path,
        default=metadata_dir / "core_meta_reg.py",
        help="Path to registry-backed metadata script",
    )
    parser.add_argument(
        "--repo_root",
        type=Path,
        default=repo_root,
        help="Repository root used as subprocess working directory",
    )
    parser.add_argument(
        "--work_dir",
        type=Path,
        help="Directory for generated outputs. Defaults to a temporary directory.",
    )
    parser.add_argument(
        "--keep_outputs",
        action="store_true",
        help="Keep output files in --work_dir for debugging",
    )
    parser.add_argument(
        "--report_json",
        type=Path,
        help="Optional path to write the comparison report JSON",
    )
    parser.add_argument(
        "--verbose_scripts",
        action="store_true",
        help="Pass -v through to both metadata scripts",
    )
    return parser


def main() -> None:
    """Entry point with compact eHive-readable failures."""
    parser = build_arg_parser()
    args = parser.parse_args()

    try:
        report = run_truth_check(args)
    except TruthMismatchError as exc:
        print("TRUTH_CHECK_FAILED: " + str(exc), file=sys.stderr)
        print(json.dumps(exc.report, indent=2, sort_keys=True), file=sys.stderr)
        raise SystemExit(exc.exit_code) from None
    except TruthCheckError as exc:
        print("TRUTH_CHECK_SCRIPT_ERROR: " + str(exc), file=sys.stderr)
        if exc.report:
            print(json.dumps(exc.report, indent=2, sort_keys=True), file=sys.stderr)
        raise SystemExit(exc.exit_code) from None

    print(
        "TRUTH_CHECK_OK: old and registry metadata JSON outputs match "
        f"for {args.db_name}; compared {report['old_metadata_count']} metadata values"
    )


if __name__ == "__main__":
    main()

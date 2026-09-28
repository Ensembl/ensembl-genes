#!/usr/bin/env python3

"""
Collect compact runtime software provenance.

The configuration key is treated only as a logical pipeline name.
The value of the "software" field is used to identify the actual software.

For example:

    "dust": {
        "software": "dustmasker"
    }

produces:

    "dustmasker": {
        "version": "...",
        "package_manager": "spack"
    }

Output schema:

{
    "software": {
        "<software executable/file name>": {
            "version": "...",
            "package_manager": "spack|homebrew|git|container|unknown",
            "git": {
                "tag": "...",
                "commit": "...",
                "dirty": false
            },
            "container": {
                "path": "...",
                "last_modified": "..."
            }
        }
    },
    "python_packages": {
        "<package>": {
            "version": "...",
            "package_manager": "git",
            "git": {
                "tag": "...",
                "commit": "...",
                "dirty": null
            }
        }
    }
}

Package-manager detection precedence:

    container
    spack
    homebrew/linuxbrew
    git
    executable-reported version
    unknown
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any

VERSION_COMMANDS = (
    ["--version"],
    ["-version"],
    ["-v"],
)

VERSION_RE = re.compile(
    r"(?<![A-Za-z0-9])" r"v?\d+(?:\.\d+)+" r"(?:[-+_.][A-Za-z0-9.-]+)?"
)

CONTAINER_SUFFIXES = {
    ".sif",
    ".img",
    ".simg",
}


# ============================================================================
# Generic helpers
# ============================================================================


def run_command(
    command: list[str],
    timeout: int = 10,
) -> dict[str, Any]:
    """Run a command and return stdout/stderr/return code."""
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
            check=False,
        )

        return {
            "returncode": result.returncode,
            "stdout": result.stdout.strip(),
            "stderr": result.stderr.strip(),
        }

    except (OSError, subprocess.SubprocessError) as exc:
        return {
            "returncode": None,
            "stdout": "",
            "stderr": "",
            "error": str(exc),
        }


def resolve_path(value: str) -> Path | None:
    """Resolve an existing filesystem path and follow symlinks."""
    path = Path(value).expanduser()

    if not path.exists():
        return None

    try:
        return path.resolve()
    except OSError:
        return path.absolute()


def resolve_executable(value: str) -> Path | None:
    """
    Resolve either:
        - an existing filesystem path
        - a command available on PATH
    """
    path = resolve_path(value)

    if path is not None:
        return path

    resolved = shutil.which(value)

    if resolved:
        try:
            return Path(resolved).resolve()
        except OSError:
            return Path(resolved)

    return None


def extract_version(text: str) -> str | None:
    """Extract a likely version number from arbitrary command output."""
    if not text:
        return None

    match = VERSION_RE.search(text)

    return match.group(0) if match else None


def load_json(path: Path) -> Any | None:
    """Load JSON safely."""
    try:
        with path.open(encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None


# ============================================================================
# Naming
# ============================================================================


def software_name_from_config(
    software_value: str,
) -> str:
    """
    Get the actual software name from the "software" value.

    Examples:

        "dustmasker"                -> "dustmasker"
        "/bin/samtools"             -> "samtools"
        "/path/test_cpc2.sif"       -> "test_cpc2.sif"
        "/path/eponine-scan.jar"    -> "eponine-scan.jar"
    """
    return Path(software_value).name


# ============================================================================
# Container handling
# ============================================================================


def is_container_path(
    software_value: str,
) -> bool:
    """Determine whether the configured software is a container image."""
    return Path(software_value).suffix.lower() in CONTAINER_SUFFIXES


def get_container_metadata(
    software_value: str,
) -> dict[str, Any]:
    """
    Return the container path and filesystem modification time.

    The path is retained even when the image does not currently exist.
    """
    path = Path(software_value).expanduser()

    if not path.exists():
        return {
            "path": str(path),
            "last_modified": None,
        }

    try:
        resolved = path.resolve()
    except OSError:
        resolved = path

    try:
        modified = datetime.fromtimestamp(
            resolved.stat().st_mtime,
            tz=timezone.utc,
        ).isoformat()
    except OSError:
        modified = None

    return {
        "path": str(resolved),
        "last_modified": modified,
    }


# ============================================================================
# Spack
# ============================================================================


def find_spack_spec(
    executable: Path,
) -> Path | None:
    """
    Find the nearest .spack/spec.json above the executable.

    This works with Spack paths containing __spack_path_placeholder__.
    """
    current = executable.parent

    while True:
        spec_path = current / ".spack" / "spec.json"

        if spec_path.is_file():
            return spec_path

        parent = current.parent

        if parent == current:
            return None

        current = parent


def get_spack_version(
    executable: Path,
) -> str | None:
    """Read the package version from Spack's spec.json."""
    spec_path = find_spack_spec(executable)

    if spec_path is None:
        return None

    spec = load_json(spec_path)

    if not isinstance(spec, dict):
        return None

    version = spec.get("version")

    if isinstance(version, str):
        return version

    package = spec.get("package")

    if isinstance(package, dict):
        version = package.get("version")

        if isinstance(version, str):
            return version

    return None


# ============================================================================
# Homebrew / Linuxbrew
# ============================================================================


def get_homebrew_metadata(
    executable: Path,
) -> tuple[str | None, str | None]:
    """
    Detect Homebrew/Linuxbrew from the resolved path.

    Supports both:

        .../Cellar/samtools/1.16.1/bin/samtools

    and:

        .../linuxbrew/opt/eponine/...
    """
    parts = executable.parts

    # ------------------------------------------------------------------
    # Standard Homebrew Cellar layout
    # ------------------------------------------------------------------

    try:
        cellar_index = parts.index("Cellar")
    except ValueError:
        cellar_index = -1

    if cellar_index >= 0:
        if len(parts) > cellar_index + 2:
            package = parts[cellar_index + 1]
            version = parts[cellar_index + 2]

            return package, version

    # ------------------------------------------------------------------
    # Linuxbrew / Homebrew opt layout
    #
    # .../linuxbrew/opt/eponine/libexec/eponine-scan.jar
    #
    # We can identify the package from "opt/<package>".
    # The version is then obtained from brew when possible.
    # ------------------------------------------------------------------

    for index, part in enumerate(parts):
        if part == "opt" and index + 1 < len(parts):
            package = parts[index + 1]

            brew = shutil.which("brew")

            if brew:
                result = run_command(
                    [
                        brew,
                        "list",
                        "--versions",
                        package,
                    ]
                )

                if result["returncode"] == 0 and result["stdout"]:
                    fields = result["stdout"].split()

                    if len(fields) >= 2:
                        return (
                            package,
                            fields[1],
                        )

            # We know the package manager even if the version
            # cannot be retrieved.
            return package, None

    return None, None


# ============================================================================
# Git
# ============================================================================


def find_git_root(
    path: Path,
) -> Path | None:
    """Return the containing Git repository, if there is one."""
    result = run_command(
        [
            "git",
            "-C",
            str(path),
            "rev-parse",
            "--show-toplevel",
        ]
    )

    if result["returncode"] != 0:
        return None

    if not result["stdout"]:
        return None

    try:
        return Path(result["stdout"]).resolve()
    except OSError:
        return Path(result["stdout"])


def get_git_metadata(
    repo: Path,
) -> dict[str, Any]:
    """Return compact Git provenance."""
    git_data: dict[str, Any] = {
        "tag": None,
        "commit": None,
        "dirty": None,
    }

    # ------------------------------------------------------------------
    # Exact commit
    # ------------------------------------------------------------------

    result = run_command(
        [
            "git",
            "-C",
            str(repo),
            "rev-parse",
            "HEAD",
        ]
    )

    if result["returncode"] == 0:
        git_data["commit"] = result["stdout"]

    # ------------------------------------------------------------------
    # Exact tags pointing to HEAD
    # ------------------------------------------------------------------

    result = run_command(
        [
            "git",
            "-C",
            str(repo),
            "tag",
            "--points-at",
            "HEAD",
        ]
    )

    if result["returncode"] == 0 and result["stdout"]:
        tags = result["stdout"].splitlines()

        if len(tags) == 1:
            git_data["tag"] = tags[0]
        elif tags:
            git_data["tag"] = tags

    # ------------------------------------------------------------------
    # Dirty state
    # ------------------------------------------------------------------

    result = run_command(
        [
            "git",
            "-C",
            str(repo),
            "status",
            "--porcelain=v1",
        ]
    )

    if result["returncode"] == 0:
        git_data["dirty"] = bool(result["stdout"])

    return git_data


# ============================================================================
# Executable version
# ============================================================================


def get_executable_version(
    executable: Path,
    config: dict[str, Any],
) -> str | None:
    """
    Try to get a version from the executable itself.

    This is only used if no package-manager-derived version was found.
    """
    custom = config.get("version_command")

    if custom is not None:
        if not isinstance(custom, list):
            return None

        commands = [
            [
                str(executable),
                *[str(item) for item in custom],
            ]
        ]
    else:
        commands = [
            [
                str(executable),
                *args,
            ]
            for args in VERSION_COMMANDS
        ]

    for command in commands:
        result = run_command(command)

        if result["returncode"] != 0:
            continue

        output = "\n".join(
            value
            for value in (
                result["stdout"],
                result["stderr"],
            )
            if value
        )

        version = extract_version(output)

        if version:
            return version

    return None


# ============================================================================
# Python Git packages
# ============================================================================


def read_direct_url(
    distribution: metadata.Distribution,
) -> dict[str, Any] | None:
    """Read PEP 610 direct_url.json."""
    try:
        text = distribution.read_text("direct_url.json")
    except (OSError, FileNotFoundError):
        return None

    if not text:
        return None

    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return None

    return data if isinstance(data, dict) else None


def collect_python_git_packages() -> dict[str, Any]:
    """
    Find installed Python distributions with Git provenance.
    """
    packages: dict[str, Any] = {}

    for distribution in metadata.distributions():

        name = distribution.metadata.get("Name") or distribution.name

        if not name:
            continue

        direct_url = read_direct_url(distribution)

        if direct_url is None:
            continue

        vcs_info = direct_url.get("vcs_info")

        # --------------------------------------------------------------
        # Standard pip Git/VCS installation
        # --------------------------------------------------------------

        if isinstance(vcs_info, dict) and vcs_info.get("vcs") == "git":
            packages[name] = {
                "version": distribution.version,
                "package_manager": "git",
                "git": {
                    "tag": vcs_info.get("requested_revision"),
                    "commit": vcs_info.get("commit_id"),
                    "dirty": None,
                },
            }

            continue

        # --------------------------------------------------------------
        # Editable local Git installation
        # --------------------------------------------------------------

        url = direct_url.get("url")

        if isinstance(url, str) and url.startswith("file://"):
            repo_path = Path(url[7:])

            if repo_path.exists():
                repo = find_git_root(repo_path)

                if repo is not None:
                    packages[name] = {
                        "version": distribution.version,
                        "package_manager": "git",
                        "git": get_git_metadata(repo),
                    }

    return dict(
        sorted(
            packages.items(),
            key=lambda item: item[0].lower(),
        )
    )


# ============================================================================
# One software entry
# ============================================================================


def collect_software(
    config: dict[str, Any],
) -> dict[str, Any]:
    """
    Collect compact provenance for one software config entry.

    The logical configuration key is deliberately not passed here.
    The actual "software" field determines the software being measured.
    """
    software_value = config.get("software")

    result: dict[str, Any] = {
        "version": None,
        "package_manager": "unknown",
    }

    if not software_value:
        return result

    software_value = str(software_value)

    # ------------------------------------------------------------------
    # Container
    #
    # Check BEFORE resolve_executable().
    # ------------------------------------------------------------------

    if is_container_path(software_value):
        result["package_manager"] = "container"
        result["container"] = get_container_metadata(software_value)

        return result

    # ------------------------------------------------------------------
    # Resolve executable
    # ------------------------------------------------------------------

    executable = resolve_executable(software_value)

    if executable is None:
        return result

    # ------------------------------------------------------------------
    # Spack
    # ------------------------------------------------------------------

    spack_version = get_spack_version(executable)

    if spack_version is not None:
        result["version"] = spack_version
        result["package_manager"] = "spack"

    # ------------------------------------------------------------------
    # Homebrew / Linuxbrew
    # ------------------------------------------------------------------

    if result["package_manager"] == "unknown":
        _, brew_version = get_homebrew_metadata(executable)

        # A Homebrew installation is identified even if brew cannot
        # provide a version.
        brew_package, brew_version = get_homebrew_metadata(executable)

        if brew_package is not None:
            result["package_manager"] = "homebrew"

            if brew_version is not None:
                result["version"] = brew_version

    # ------------------------------------------------------------------
    # Git
    #
    # Git is ONLY a fallback after Spack/Homebrew/container.
    # ------------------------------------------------------------------

    if result["package_manager"] == "unknown":
        git_root = find_git_root(executable.parent)

        if git_root is not None:
            git_data = get_git_metadata(git_root)

            result["package_manager"] = "git"
            result["git"] = git_data

            # An exact tag is useful as a version.
            tag = git_data.get("tag")

            if isinstance(tag, str):
                result["version"] = tag

    # ------------------------------------------------------------------
    # Executable-reported version
    # ------------------------------------------------------------------

    if result["version"] is None:
        result["version"] = get_executable_version(
            executable,
            config,
        )

    return result


# ============================================================================
# Complete collection
# ============================================================================


def collect_all(
    config: dict[str, Any],
) -> dict[str, Any]:
    """
    Collect all software.

    The output key is taken from the actual "software" value, NOT the
    logical configuration key.

    Example:

        "dust": {
            "software": "dustmasker"
        }

    becomes:

        "dustmasker": {...}
    """
    software: dict[str, Any] = {}

    for _, entry in config.items():

        if not isinstance(entry, dict):
            continue

        software_value = entry.get("software")

        if not software_value:
            continue

        software_name = software_name_from_config(str(software_value))

        software[software_name] = collect_software(entry)

    return {
        "software": software,
        "python_packages": (collect_python_git_packages()),
    }


# ============================================================================
# Command line
# ============================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description=("Collect compact software provenance")
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Software configuration JSON",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output JSON file",
    )

    args = parser.parse_args()

    with open(
        args.config,
        encoding="utf-8",
    ) as handle:
        config = json.load(handle)

    if not isinstance(config, dict):
        raise ValueError("Software config must be a JSON object")

    result = collect_all(config)

    with open(
        args.output,
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            result,
            handle,
            indent=2,
        )
        handle.write("\n")


if __name__ == "__main__":
    main()

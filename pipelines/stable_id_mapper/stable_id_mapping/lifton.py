"""LiftOn execution helpers."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

DEFAULT_LIFTON_FEATURE_TYPES = (
    "CDS",
    "exon",
    "five_prime_UTR",
    "gene",
    "mRNA",
    "three_prime_UTR",
)
GENERATED_FEATURE_TYPES_FILE_NAME = "lifton_feature_types.txt"


@dataclass(frozen=True)
class LiftonRunConfig:
    ref_gff: Path
    ref_fasta: Path
    target_fasta: Path
    output_gff: Path
    threads: int = 8
    executable: str = "lifton"
    feature_types: tuple[str, ...] = DEFAULT_LIFTON_FEATURE_TYPES
    feature_types_file: Optional[Path] = None
    extra_args: tuple[str, ...] = field(default_factory=tuple)
    work_dir: Optional[Path] = None

    def validate_inputs(self) -> None:
        if self.threads < 1:
            raise ValueError("threads must be >= 1")
        if self.feature_types_file is None and not self.feature_types:
            raise ValueError("At least one feature type is required")
        for path in (self.ref_gff, self.ref_fasta, self.target_fasta):
            if not path.exists():
                raise FileNotFoundError(path)
        if self.feature_types_file is not None and not self.feature_types_file.exists():
            raise FileNotFoundError(self.feature_types_file)
        if self.work_dir is not None and not self.work_dir.exists():
            raise FileNotFoundError(self.work_dir)

    def feature_types_arg(self) -> str:
        if self.feature_types_file is not None:
            return str(self.feature_types_file)
        return str(self.generated_feature_types_file)

    @property
    def generated_feature_types_file(self) -> Path:
        return self.output_gff.parent / GENERATED_FEATURE_TYPES_FILE_NAME


def prepare_lifton_feature_types_file(config: LiftonRunConfig) -> Optional[Path]:
    """Create the LiftOn feature-type file when inline feature types are configured."""
    config.validate_inputs()
    if config.feature_types_file is not None:
        return None
    feature_types_file = config.generated_feature_types_file
    feature_types_file.parent.mkdir(parents=True, exist_ok=True)
    feature_types_file.write_text(
        "".join(f"{feature_type}\n" for feature_type in config.feature_types),
        encoding="utf-8",
    )
    return feature_types_file


def build_lifton_command(config: LiftonRunConfig) -> list[str]:
    config.validate_inputs()
    return [
        config.executable,
        "-f",
        config.feature_types_arg(),
        "-t",
        str(config.threads),
        "-g",
        str(config.ref_gff),
        "-o",
        str(config.output_gff),
        *config.extra_args,
        str(config.ref_fasta),
        str(config.target_fasta),
    ]


def find_lifton_executable(executable: str) -> str:
    if Path(executable).parent != Path("."):
        if Path(executable).exists():
            return executable
        raise FileNotFoundError(f"LiftOn executable not found: {executable}")
    resolved = shutil.which(executable)
    if resolved is None:
        raise FileNotFoundError(
            f"LiftOn executable {executable!r} was not found on PATH"
        )
    return resolved


def run_lifton(
    config: LiftonRunConfig,
    stdout=None,
    stderr=None,
) -> subprocess.CompletedProcess:
    config.validate_inputs()
    find_lifton_executable(config.executable)
    config.output_gff.parent.mkdir(parents=True, exist_ok=True)
    prepare_lifton_feature_types_file(config)
    command = build_lifton_command(config)
    completed = subprocess.run(
        command,
        cwd=config.work_dir,
        check=True,
        stdout=stdout,
        stderr=stderr,
    )
    if not config.output_gff.exists():
        raise RuntimeError(
            f"LiftOn completed but expected output was not created: {config.output_gff}"
        )
    return completed


def format_command(command: Sequence[str]) -> str:
    return " ".join(_shell_quote(part) for part in command)


def _shell_quote(value: str) -> str:
    if value and all(char.isalnum() or char in "-_./:=,+" for char in value):
        return value
    return "'" + value.replace("'", "'\"'\"'") + "'"

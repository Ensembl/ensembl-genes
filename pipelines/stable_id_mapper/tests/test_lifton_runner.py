from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stable_id_mapping.lifton import (
    LiftonRunConfig,
    build_lifton_command,
    find_lifton_executable,
    format_command,
)


def touch(path: Path) -> None:
    path.write_text("", encoding="utf-8")


def make_config(tmp_path: Path, **kwargs) -> LiftonRunConfig:
    ref_gff = tmp_path / "ref.gff3"
    ref_fasta = tmp_path / "ref.fa"
    target_fasta = tmp_path / "target.fa"
    for path in (ref_gff, ref_fasta, target_fasta):
        touch(path)
    defaults = {
        "ref_gff": ref_gff,
        "ref_fasta": ref_fasta,
        "target_fasta": target_fasta,
        "output_gff": tmp_path / "lifton.gff3",
        "threads": 4,
        "executable": "lifton",
        "feature_types": ("gene", "mRNA", "exon", "CDS"),
        "extra_args": ("--debug",),
    }
    defaults.update(kwargs)
    return LiftonRunConfig(**defaults)


def test_build_lifton_command(tmp_path: Path) -> None:
    config = make_config(tmp_path)

    command = build_lifton_command(config)

    assert command == [
        "lifton",
        "-f",
        "gene,mRNA,exon,CDS",
        "-t",
        "4",
        "-g",
        str(config.ref_gff),
        "-o",
        str(config.output_gff),
        "--debug",
        str(config.ref_fasta),
        str(config.target_fasta),
    ]


def test_feature_types_file_overrides_inline_feature_types(tmp_path: Path) -> None:
    feature_types_file = tmp_path / "feature_types.txt"
    feature_types_file.write_text("gene\nmRNA\n", encoding="utf-8")
    config = make_config(
        tmp_path,
        feature_types_file=feature_types_file,
        feature_types=("ignored",),
    )

    command = build_lifton_command(config)

    assert command[2] == str(feature_types_file)


def test_lifton_config_requires_inputs(tmp_path: Path) -> None:
    config = LiftonRunConfig(
        ref_gff=tmp_path / "missing.gff3",
        ref_fasta=tmp_path / "missing.fa",
        target_fasta=tmp_path / "missing_target.fa",
        output_gff=tmp_path / "out.gff3",
    )

    with pytest.raises(FileNotFoundError):
        build_lifton_command(config)


def test_missing_lifton_executable_error_is_clear() -> None:
    with pytest.raises(FileNotFoundError, match="was not found on PATH"):
        find_lifton_executable("definitely-not-lifton")


def test_format_command_quotes_spaces() -> None:
    command = ["lifton", "-g", "/tmp/ref genes.gff3", "/tmp/ref.fa", "/tmp/tar.fa"]

    assert format_command(command) == "lifton -g '/tmp/ref genes.gff3' /tmp/ref.fa /tmp/tar.fa"


from __future__ import annotations

import json
import os
from typing import Any, Dict

import pandas as pd

from .utils import ensure_dir


def _can_parquet() -> bool:
    try:
        import pyarrow  # noqa: F401

        return True
    except Exception:
        return False


def write_df(df: pd.DataFrame, path: str, fmt: str = "tsv") -> None:
    ensure_dir(os.path.dirname(path))
    fmt = fmt.lower()
    if fmt == "parquet":
        if not _can_parquet():
            fmt = "tsv"
        else:
            df.to_parquet(path, index=False)
            return
    if fmt == "tsv":
        df.to_csv(path, sep="\t", index=False)
    elif fmt == "csv":
        df.to_csv(path, index=False)
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def append_tsv(df: pd.DataFrame, path: str, header_if_new: bool = True) -> None:
    ensure_dir(os.path.dirname(path))
    write_header = header_if_new and not os.path.exists(path)
    df.to_csv(path, sep="\t", index=False, mode="a", header=write_header)


def write_manifest(manifest: Dict[str, Any], path: str) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)

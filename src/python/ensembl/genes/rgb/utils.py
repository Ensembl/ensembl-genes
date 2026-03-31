from __future__ import annotations

import hashlib
import json
import os
import re
from datetime import datetime
from typing import Any, Dict


def now_utc_stamp() -> str:
    return datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")


def sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def make_run_id(core_db: str, layer_db: str) -> str:
    ts = now_utc_stamp()
    tag = sanitize(f"{core_db}__{layer_db}")[:80]
    return f"{ts}__{tag}"


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def sha256_of_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_json(obj: Dict[str, Any], path: str) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)

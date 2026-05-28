# id_manager.py
from __future__ import annotations
import re
from typing import Dict, Tuple

# Match the *core* part only (no namespace), e.g. ENSG00000123456 or ENSCGRG00015017959.3
CORE_RE = re.compile(r"^(?P<prefix>[A-Za-z]+)(?P<num>\d+)(?:\.(?P<ver>\d+))?$")

def split_ns(stable_id: str) -> tuple[str, str]:
    """
    Split optional namespace. Returns (namespace, core).
    Examples:
      'gene:ENSMUSG00000000001.3' -> ('gene', 'ENSMUSG00000000001.3')
      'ENSG00000141510.16'        -> ('',     'ENSG00000141510.16')
    """
    if ":" in stable_id:
        ns, core = stable_id.split(":", 1)
        return ns, core
    return "", stable_id

def join_ns(ns: str, core: str) -> str:
    """Reattach namespace if present."""
    return f"{ns}:{core}" if ns else core

class IDManager:
    """
    Tracks maximum numeric part per *core prefix* (e.g. ENSG, ENSMUSG, ENSBBBG)
    across ref+target. Provides:
      - increment_version(stable_id)  # preserves namespace/prefix and width
      - new_id(prefix, ns='')         # generates ns:prefix#########.1
    Also records a mapping log for reporting.
    """
    def __init__(self, ref: Dict, tgt: Dict) -> None:
        self.maxnum: Dict[str, int] = {}
        self.width: Dict[str, int] = {}
        # Scan both reference and target to learn max numbers and widths
        def _touch_any(sid: str) -> None:
            ns, core = split_ns(sid)
            m = CORE_RE.match(core)
            if not m:
                return
            pre, num_s = m["prefix"], m["num"]
            num = int(num_s)
            self.maxnum[pre] = max(self.maxnum.get(pre, 0), num)
            self.width[pre] = max(self.width.get(pre, 0), len(num_s))

        for gene_map in (ref, tgt):
            for gid, g in gene_map.items():
                _touch_any(g.attrs.get("ID", gid))
                for tid, t in g.transcripts.items():
                    _touch_any(t.attrs.get("ID", t.tid))

        # old_id -> (new_id, status: 'mapped'|'missing'|'new')
        self.mapping: Dict[str, Tuple[str, str]] = {}

    @staticmethod
    def parse(stable_id: str) -> Tuple[str, int, int, str]:
        """
        Parse possibly-namespaced ID -> (prefix, number, version, namespace).
        Version defaults to 1 if absent.
        """
        ns, core = split_ns(stable_id)
        m = CORE_RE.match(core)
        if not m:
            raise ValueError(f"Bad stable ID: {stable_id}")
        pre = m["prefix"]
        num = int(m["num"])
        ver = int(m["ver"]) if m["ver"] else 1
        return pre, num, ver, ns

    def increment_version(self, stable_id: str) -> str:
        """
        Bump version while preserving namespace, prefix, and numeric width.
        If no version present, assumes .1 -> returns .2.
        """
        pre, num, ver, ns = self.parse(stable_id)
        w = self.width.get(pre, len(str(num)))
        core = f"{pre}{num:0{w}d}.{ver+1}"
        return join_ns(ns, core)

    def new_id(self, prefix: str, ns: str = "") -> str:
        """
        Create a new ID with the learned width for this prefix and version .1.
        Namespace is applied if provided (e.g., ns='gene' or 'transcript').
        """
        nxt = self.maxnum.get(prefix, 0) + 1
        self.maxnum[prefix] = nxt
        w = max(self.width.get(prefix, len(str(nxt))), len(str(nxt)))
        core = f"{prefix}{nxt:0{w}d}.1"
        return join_ns(ns, core)

    def set(self, old_id: str, new_id: str, status: str) -> None:
        self.mapping[old_id] = (new_id, status)

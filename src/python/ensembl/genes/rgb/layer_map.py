from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import requests


RAW_URL_DEFAULT = (
    "https://raw.githubusercontent.com/Ensembl/ensembl-analysis/refs/heads/main/"
    "modules/Bio/EnsEMBL/Analysis/Hive/Config/LayerAnnotationStatic.pm"
)


@dataclass
class LayerEntry:
    layer_id: str  # e.g., LAYER1
    biotypes: List[str]


@dataclass
class LayerMap:
    config_key: str  # e.g., mammals_basic
    layers: List[LayerEntry]

    def biotype_to_tier(self) -> Dict[str, Tuple[str, int]]:
        m: Dict[str, Tuple[str, int]] = {}
        for idx, le in enumerate(self.layers, start=1):
            for bt in le.biotypes:
                m[bt] = (le.layer_id, idx)
        return m


def fetch_layer_static(url: str = RAW_URL_DEFAULT) -> str:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text


def parse_config_block(text: str, config_key: str) -> LayerMap:
    # Find the array assigned to config_key => [ ... ]
    # Use a naive bracket counter for square brackets starting at the '[' after the key
    key_pattern = re.compile(r"\b" + re.escape(config_key) + r"\s*=>\s*\[")
    m = key_pattern.search(text)
    if not m:
        raise ValueError(f"config_key '{config_key}' not found in LayerAnnotationStatic.pm")
    start = m.end()  # position after '['
    depth = 1
    i = start
    n = len(text)
    while i < n and depth > 0:
        ch = text[i]
        if ch == '[':
            depth += 1
        elif ch == ']':
            depth -= 1
        i += 1
    block = text[start : i - 1]

    # Extract layer entries: { ID => 'LAYERX', ... BIOTYPES => [ 'a', 'b', ... ] }
    # We will iterate over each hash entry between { and }
    layers: List[LayerEntry] = []
    for mobj in re.finditer(r"\{([^{}]+)\}", block, flags=re.S):
        chunk = mobj.group(1)
        idm = re.search(r"ID\s*=>\s*'([^']+)'", chunk)
        bm = re.search(r"BIOTYPES\s*=>\s*\[(.*?)\]", chunk, flags=re.S)
        if not idm or not bm:
            continue
        layer_id = idm.group(1).strip()
        # Split quoted tokens inside BIOTYPES list
        biotypes = [t.strip() for t in re.findall(r"'([^']+)'", bm.group(1))]
        layers.append(LayerEntry(layer_id=layer_id, biotypes=biotypes))

    if not layers:
        raise ValueError(f"No layers parsed for '{config_key}'")
    return LayerMap(config_key=config_key, layers=layers)


def to_yaml(mapping: Dict[str, Tuple[str, int]]) -> str:
    # Produce a simple YAML: biotype: {tier: LAYERX, tier_index: N}
    lines = ["# Auto-generated from LayerAnnotationStatic.pm", "# biotype -> tier mapping", ""]
    for bt, (tier, idx) in sorted(mapping.items()):
        lines.append(f"{bt}: {{tier: {tier}, tier_index: {idx}}}")
    return "\n".join(lines) + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description="Parse LayerAnnotationStatic.pm to a YAML biotype->tier map")
    p.add_argument("--config_key", default="mammals_basic")
    p.add_argument("--layer_static_url", default=RAW_URL_DEFAULT)
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    text = fetch_layer_static(args.layer_static_url)
    lm = parse_config_block(text, args.config_key)
    mapping = lm.biotype_to_tier()
    yaml_text = to_yaml(mapping)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(yaml_text)
    print(f"[layer-map] wrote {len(mapping)} biotypes for '{args.config_key}' → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

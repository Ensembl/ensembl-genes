"""Stable-ID parsing and allocation helpers."""

from __future__ import annotations

import argparse
from dataclasses import dataclass

from .models import Feature


@dataclass(frozen=True)
class StableIdRange:
    prefix: str
    start: int
    end: int
    width: int


@dataclass
class IdAllocator:
    prefix: str
    next_value: int
    end: int
    width: int
    reserved: set[str]

    def allocate(self) -> str:
        while self.next_value <= self.end:
            candidate = f"{self.prefix}{self.next_value:0{self.width}d}"
            self.next_value += 1
            if candidate in self.reserved:
                continue
            self.reserved.add(candidate)
            return candidate
        raise RuntimeError(
            f"Ran out of IDs for {self.prefix}:{self.next_value}-{self.end}"
        )


def parse_id_range(value: str) -> StableIdRange:
    try:
        prefix, interval = value.split(":", 1)
        start_text, end_text = interval.split("-", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Bad range {value!r}; expected PREFIX:START-END"
        ) from exc
    if not prefix:
        raise argparse.ArgumentTypeError("Range prefix cannot be empty")
    if not start_text.isdigit() or not end_text.isdigit():
        raise argparse.ArgumentTypeError("Range START and END must be integers")
    start = int(start_text)
    end = int(end_text)
    if end < start:
        raise argparse.ArgumentTypeError("Range END must be >= START")
    return StableIdRange(prefix, start, end, max(len(start_text), len(end_text)))


def make_allocator(
    id_range: StableIdRange,
    reserved_ids: set[str],
) -> IdAllocator:
    return IdAllocator(
        prefix=id_range.prefix,
        next_value=id_range.start,
        end=id_range.end,
        width=id_range.width,
        reserved=set(reserved_ids),
    )


def collect_reserved_ids(*feature_sets: dict[str, dict[str, Feature]]) -> set[str]:
    reserved: set[str] = set()
    for features_by_type in feature_sets:
        for features in features_by_type.values():
            reserved.update(features)
    return reserved


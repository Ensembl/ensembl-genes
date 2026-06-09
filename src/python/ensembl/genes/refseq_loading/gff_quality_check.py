"""Post-load quality checks for GFF imports into Ensembl core databases."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Mapping, Protocol, Sequence

try:  # Support both package imports and direct same-directory imports.
    from .gff_models import ParsedAnnotation
    from .gff_source_config import GffSourceConfig
except ImportError:  # pragma: no cover - used when run beside this file.
    from gff_models import ParsedAnnotation  # type: ignore
    from gff_source_config import GffSourceConfig  # type: ignore


LOGGER = logging.getLogger(__name__)


class QualityCursor(Protocol):
    """Cursor protocol required by the post-load quality checker."""

    def execute(self, operation: str, params: Any | None = None) -> Any:
        """Execute one SQL statement."""
        ...

    def fetchall(self) -> Sequence[tuple[Any, ...]]:
        """Fetch all result rows."""
        ...


@dataclass(frozen=True)
class QualityMetric:
    """One expected-vs-observed quality metric."""

    name: str
    expected: int
    observed: int
    missing: int = 0
    unexpected: int = 0
    mismatched: int = 0

    @property
    def passed(self) -> bool:
        """Return True when this metric has no discrepancy."""

        return self.missing == 0 and self.unexpected == 0 and self.mismatched == 0


@dataclass
class CoreLoadQualityReport:
    """Summary of a post-load comparison between parsed GFF and core rows."""

    source_name: str
    metrics: list[QualityMetric] = field(default_factory=list)
    missing_examples: dict[str, list[str]] = field(default_factory=dict)
    unexpected_examples: dict[str, list[str]] = field(default_factory=dict)
    mismatch_examples: dict[str, list[str]] = field(default_factory=dict)
    gene_biotypes: dict[str, int] = field(default_factory=dict)
    transcript_biotypes: dict[str, int] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        """Return True when every quality metric passed."""

        return all(metric.passed for metric in self.metrics)

    def failure_summary(self) -> str:
        """Return a compact summary of all quality check discrepancies."""

        missing = sum(metric.missing for metric in self.metrics)
        unexpected = sum(metric.unexpected for metric in self.metrics)
        mismatched = sum(metric.mismatched for metric in self.metrics)
        return (
            f"missing={missing} unexpected={unexpected} mismatched={mismatched}"
        )

    def render_lines(self) -> list[str]:
        """Return printable/loggable report lines."""

        status = "PASS" if self.passed else "FAIL"
        lines = [f"GFF core quality check: {status}", f"source: {self.source_name}"]
        lines.append("expected vs core rows:")
        for metric in self.metrics:
            lines.append(
                "  "
                f"{metric.name}: expected={metric.expected} observed={metric.observed} "
                f"missing={metric.missing} unexpected={metric.unexpected} "
                f"mismatched={metric.mismatched}"
            )

        if self.gene_biotypes:
            lines.append(
                "gene biotypes: "
                + ", ".join(
                    f"{biotype}={count}"
                    for biotype, count in sorted(self.gene_biotypes.items())
                )
            )
        if self.transcript_biotypes:
            lines.append(
                "transcript biotypes: "
                + ", ".join(
                    f"{biotype}={count}"
                    for biotype, count in sorted(self.transcript_biotypes.items())
                )
            )

        for label, examples in sorted(self.missing_examples.items()):
            if examples:
                lines.append(f"missing {label} examples: {', '.join(examples[:10])}")
        for label, examples in sorted(self.unexpected_examples.items()):
            if examples:
                lines.append(
                    f"unexpected {label} examples: {', '.join(examples[:10])}"
                )
        for label, examples in sorted(self.mismatch_examples.items()):
            if examples:
                lines.append(f"mismatched {label} examples: {', '.join(examples[:10])}")

        return lines


def _range_from_ids(ids: Sequence[int]) -> tuple[int, int] | None:
    """Return the inclusive range for a non-empty list of allocated IDs."""

    if not ids:
        return None
    return min(ids), max(ids)


def _fetch_id_to_values(
    cursor: QualityCursor,
    table_name: str,
    id_column: str,
    value_columns: Sequence[str],
    id_values: Sequence[int],
) -> dict[int, tuple[Any, ...]]:
    """Fetch rows from a table for a contiguous set of allocated IDs."""

    id_range = _range_from_ids(id_values)
    if id_range is None:
        return {}

    selected_columns = ", ".join((id_column, *value_columns))
    cursor.execute(
        f"SELECT {selected_columns} FROM {table_name} "
        f"WHERE {id_column} BETWEEN %s AND %s",
        id_range,
    )
    return {int(row[0]): tuple(row[1:]) for row in cursor.fetchall()}


def _fetch_exon_transcript_pairs(
    cursor: QualityCursor,
    exon_ids: Sequence[int],
    transcript_ids: Sequence[int],
) -> set[tuple[int, int]]:
    """Fetch exon_transcript links for allocated exon and transcript IDs."""

    exon_range = _range_from_ids(exon_ids)
    transcript_range = _range_from_ids(transcript_ids)
    if exon_range is None or transcript_range is None:
        return set()

    cursor.execute(
        "SELECT exon_id, transcript_id FROM exon_transcript "
        "WHERE exon_id BETWEEN %s AND %s "
        "AND transcript_id BETWEEN %s AND %s",
        (*exon_range, *transcript_range),
    )
    return {(int(row[0]), int(row[1])) for row in cursor.fetchall()}


def _fetch_translations(
    cursor: QualityCursor,
    transcript_ids: Sequence[int],
) -> dict[int, str | None]:
    """Fetch translations attached to allocated transcript IDs."""

    transcript_range = _range_from_ids(transcript_ids)
    if transcript_range is None:
        return {}

    cursor.execute(
        "SELECT transcript_id, stable_id FROM translation "
        "WHERE transcript_id BETWEEN %s AND %s",
        transcript_range,
    )
    return {int(row[0]): row[1] for row in cursor.fetchall()}


def expected_translation_stable_ids(annotation: ParsedAnnotation) -> dict[str, str]:
    """Return transcript ID to expected translation stable ID for loadable CDS."""

    expected: dict[str, str] = {}
    for transcript_id, cds_list in annotation.cds_segments.items():
        transcript = annotation.transcripts.get(transcript_id)
        if transcript is None or not cds_list or not transcript.exons:
            continue

        if transcript.strand == 1:
            translation_start_pos = min(cds.start for cds in cds_list)
            translation_end_pos = max(cds.end for cds in cds_list)
        else:
            translation_start_pos = max(cds.end for cds in cds_list)
            translation_end_pos = min(cds.start for cds in cds_list)

        start_found = any(
            exon.start <= translation_start_pos <= exon.end for exon in transcript.exons
        )
        end_found = any(
            exon.start <= translation_end_pos <= exon.end for exon in transcript.exons
        )
        if start_found and end_found:
            expected[transcript_id] = transcript.protein_id or f"{transcript_id}_prot"

    return expected


def _count_biotypes(records: Sequence[str]) -> dict[str, int]:
    """Count biotype strings in stable sorted form."""

    counts: dict[str, int] = {}
    for biotype in records:
        counts[biotype] = counts.get(biotype, 0) + 1
    return dict(sorted(counts.items()))


def run_core_load_quality_check(
    cursor: QualityCursor,
    annotation: ParsedAnnotation,
    gene_id_map: Mapping[str, int],
    transcript_id_map: Mapping[str, int],
    exon_id_map: Mapping[tuple[str, int], int],
    source_config: GffSourceConfig,
) -> CoreLoadQualityReport:
    """Compare rows inserted into core against the parsed GFF annotation."""

    report = CoreLoadQualityReport(source_name=source_config.name)
    report.gene_biotypes = _count_biotypes(
        [gene.biotype for gene in annotation.genes.values()]
    )
    report.transcript_biotypes = _count_biotypes(
        [transcript.biotype for transcript in annotation.transcripts.values()]
    )

    expected_gene_by_db_id = {
        db_id: (gene.stable_id, gene.biotype)
        for gene_id, db_id in gene_id_map.items()
        for gene in (annotation.genes[gene_id],)
    }
    actual_gene_by_db_id = _fetch_id_to_values(
        cursor,
        "gene",
        "gene_id",
        ("stable_id", "biotype"),
        list(expected_gene_by_db_id),
    )
    _add_id_value_metric(
        report,
        "genes",
        expected_gene_by_db_id,
        actual_gene_by_db_id,
    )

    expected_transcript_by_db_id = {
        db_id: (transcript.stable_id, transcript.biotype)
        for transcript_id, db_id in transcript_id_map.items()
        for transcript in (annotation.transcripts[transcript_id],)
    }
    actual_transcript_by_db_id = _fetch_id_to_values(
        cursor,
        "transcript",
        "transcript_id",
        ("stable_id", "biotype"),
        list(expected_transcript_by_db_id),
    )
    _add_id_value_metric(
        report,
        "transcripts",
        expected_transcript_by_db_id,
        actual_transcript_by_db_id,
    )

    expected_exon_by_db_id: dict[int, tuple[str | None]] = {}
    for transcript_id, transcript in annotation.transcripts.items():
        sorted_exons = sorted(
            transcript.exons,
            key=lambda exon: exon.start,
            reverse=(transcript.strand == -1),
        )
        for rank, exon in enumerate(sorted_exons, start=1):
            db_id = exon_id_map[(transcript_id, rank)]
            expected_exon_by_db_id[db_id] = (exon.stable_id,)
    actual_exon_by_db_id = _fetch_id_to_values(
        cursor,
        "exon",
        "exon_id",
        ("stable_id",),
        list(expected_exon_by_db_id),
    )
    _add_id_value_metric(
        report,
        "exons",
        expected_exon_by_db_id,
        actual_exon_by_db_id,
        ignore_expected_none=True,
    )

    expected_exon_transcript_pairs = {
        (exon_id, transcript_id_map[transcript_id])
        for (transcript_id, _rank), exon_id in exon_id_map.items()
    }
    actual_exon_transcript_pairs = _fetch_exon_transcript_pairs(
        cursor,
        list(expected_exon_by_db_id),
        list(expected_transcript_by_db_id),
    )
    missing_pairs = expected_exon_transcript_pairs - actual_exon_transcript_pairs
    unexpected_pairs = actual_exon_transcript_pairs - expected_exon_transcript_pairs
    report.metrics.append(
        QualityMetric(
            name="exon_transcript_links",
            expected=len(expected_exon_transcript_pairs),
            observed=len(actual_exon_transcript_pairs),
            missing=len(missing_pairs),
            unexpected=len(unexpected_pairs),
        )
    )
    if missing_pairs:
        report.missing_examples["exon_transcript_links"] = [
            f"exon_id={exon_id}:transcript_id={transcript_id}"
            for exon_id, transcript_id in sorted(missing_pairs)[:10]
        ]
    if unexpected_pairs:
        report.unexpected_examples["exon_transcript_links"] = [
            f"exon_id={exon_id}:transcript_id={transcript_id}"
            for exon_id, transcript_id in sorted(unexpected_pairs)[:10]
        ]

    expected_translation_by_transcript_id = expected_translation_stable_ids(annotation)
    expected_translation_by_db_transcript_id = {
        transcript_id_map[transcript_id]: stable_id
        for transcript_id, stable_id in expected_translation_by_transcript_id.items()
        if transcript_id in transcript_id_map
    }
    actual_translation_by_db_transcript_id = _fetch_translations(
        cursor,
        list(expected_transcript_by_db_id),
    )
    _add_id_value_metric(
        report,
        "translations",
        {
            db_id: (stable_id,)
            for db_id, stable_id in expected_translation_by_db_transcript_id.items()
        },
        {
            db_id: (stable_id,)
            for db_id, stable_id in actual_translation_by_db_transcript_id.items()
        },
    )

    return report


def _add_id_value_metric(
    report: CoreLoadQualityReport,
    name: str,
    expected: Mapping[int, tuple[Any, ...]],
    actual: Mapping[int, tuple[Any, ...]],
    ignore_expected_none: bool = False,
) -> None:
    """Add a metric comparing expected DB IDs and selected column values."""

    missing_ids = sorted(set(expected) - set(actual))
    unexpected_ids = sorted(set(actual) - set(expected))
    mismatched_ids: list[int] = []
    for db_id in sorted(set(expected) & set(actual)):
        expected_values = expected[db_id]
        actual_values = actual[db_id]
        if ignore_expected_none and all(value is None for value in expected_values):
            continue
        if expected_values != actual_values:
            mismatched_ids.append(db_id)

    report.metrics.append(
        QualityMetric(
            name=name,
            expected=len(expected),
            observed=len(actual),
            missing=len(missing_ids),
            unexpected=len(unexpected_ids),
            mismatched=len(mismatched_ids),
        )
    )
    if missing_ids:
        report.missing_examples[name] = [str(db_id) for db_id in missing_ids[:10]]
    if unexpected_ids:
        report.unexpected_examples[name] = [
            str(db_id) for db_id in unexpected_ids[:10]
        ]
    if mismatched_ids:
        report.mismatch_examples[name] = [
            f"{db_id}: expected={expected[db_id]} observed={actual[db_id]}"
            for db_id in mismatched_ids[:10]
        ]


def emit_quality_report(
    report: CoreLoadQualityReport,
    logger: logging.Logger | None = None,
) -> None:
    """Print and log a post-load quality report."""

    log = logger or LOGGER
    for line in report.render_lines():
        print(line)
        log.info(line)

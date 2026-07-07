#!/usr/bin/env python3
"""
Generate SQL stable-ID updates from projected/mapped GFF3 output.

This remains as a compatibility entrypoint. The implementation lives in the
stable_id_mapping package so the same decision/SQL logic can be reused by the
future Lifton-first pipeline.
"""

from __future__ import annotations

import argparse
import tempfile
import textwrap
from pathlib import Path
from typing import Optional

from stable_id_mapping.config import StableIdEventConfig, default_backup_prefix
from stable_id_mapping.ids import parse_id_range
from stable_id_mapping.models import Decision
from stable_id_mapping.pipeline import run_pipeline


def write_test_file(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip(), encoding="utf-8")


def assert_has_decision(
    decisions: list[Decision],
    feature_type: str,
    action: str,
    current_stable_id: Optional[str] = None,
    old_stable_id: Optional[str] = None,
    new_stable_id: Optional[str] = None,
    new_version: Optional[int] = None,
) -> None:
    for decision in decisions:
        if decision.feature_type != feature_type or decision.action != action:
            continue
        if current_stable_id is not None and decision.current_stable_id != current_stable_id:
            continue
        if old_stable_id is not None and decision.old_stable_id != old_stable_id:
            continue
        if new_stable_id is not None and decision.new_stable_id != new_stable_id:
            continue
        if new_version is not None and decision.new_version != new_version:
            continue
        return
    raise AssertionError(
        "Missing decision "
        f"type={feature_type} action={action} current={current_stable_id} "
        f"old={old_stable_id} new={new_stable_id} new_version={new_version}"
    )


def run_tests() -> None:
    old_g1 = "ENSNWIG00000000001"
    old_g2 = "ENSNWIG00000000002"
    old_g3 = "ENSNWIG00000000003"
    old_t1 = "ENSNWIT00000000001"
    old_t2 = "ENSNWIT00000000002"
    old_t3 = "ENSNWIT00000000003"
    target_g1 = "ENSNWIG00000090001"
    target_g2 = "ENSNWIG00000090002"
    target_g3 = "ENSNWIG00000090003"
    target_t1 = "ENSNWIT00000090001"
    target_t2 = "ENSNWIT00000090002"
    target_t3 = "ENSNWIT00000090003"

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        ref_gff = tmp / "ref.gff3"
        target_gff = tmp / "target.gff3"
        mapped_gff = tmp / "mapped.gff3"
        report = tmp / "report.txt"
        output_sql = tmp / "out.sql"
        output_tsv = tmp / "out.tsv"

        write_test_file(
            ref_gff,
            f"""
            ##gff-version 3
            chr1\ttest\tregion\t1\t9999\t.\t.\t.\tID=region:chr1
            chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:{old_g1}.1
            chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:{old_t1}.1;Parent=gene:{old_g1}.1
            chr1\ttest\tCDS\t150\t450\t.\t+\t0\tParent=transcript:{old_t1}.1
            chr1\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene:{old_g2}.1
            chr1\ttest\tmRNA\t1000\t1500\t.\t+\t.\tID=transcript:{old_t2}.1;Parent=gene:{old_g2}.1
            chr1\ttest\tgene\t2000\t2500\t.\t+\t.\tID=gene:{old_g3}.1
            chr1\ttest\tmRNA\t2000\t2500\t.\t+\t.\tID=transcript:{old_t3}.1;Parent=gene:{old_g3}.1
            chr1\ttest\tCDS\t2050\t2450\t.\t+\t0\tParent=transcript:{old_t3}.1
            """,
        )
        write_test_file(
            target_gff,
            f"""
            ##gff-version 3
            chrT\ttest\tregion\t1\t9999\t.\t.\t.\tID=region:chrT
            chrT\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:{target_g1}.1
            chrT\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:{target_t1}.1;Parent=gene:{target_g1}.1
            chrT\ttest\tCDS\t150\t450\t.\t+\t0\tParent=transcript:{target_t1}.1
            chrT\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene:{target_g2}.1
            chrT\ttest\tmRNA\t1000\t1500\t.\t+\t.\tID=transcript:{target_t2}.1;Parent=gene:{target_g2}.1
            chrT\ttest\tgene\t3000\t3500\t.\t+\t.\tID=gene:{target_g3}.1
            chrT\ttest\tmRNA\t3000\t3500\t.\t+\t.\tID=transcript:{target_t3}.1;Parent=gene:{target_g3}.1
            chrT\ttest\tCDS\t3050\t3450\t.\t+\t0\tParent=transcript:{target_t3}.1
            """,
        )
        write_test_file(
            mapped_gff,
            f"""
            ##gff-version 3
            chrT\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:{old_g1}.1
            chrT\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:{old_t1}.1;Parent=gene:{old_g1}.1
            chrT\ttest\texon\t100\t500\t.\t+\t.\tID=exon:{old_t1}.exon1;Parent=transcript:{old_t1}.1
            chrT\ttest\tCDS\t150\t450\t.\t+\t0\tParent=transcript:{old_t1}.1
            chrT\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene:{old_g3}.1
            chrT\ttest\tmRNA\t1000\t1500\t.\t+\t.\tID=transcript:{old_t3}.1;Parent=gene:{old_g3}.1
            chrT\ttest\tCDS\t1050\t1450\t.\t+\t0\tParent=transcript:{old_t3}.1
            """,
        )
        write_test_file(
            report,
            f"""
            Stable ID mapper report

            Missing gene IDs:
              gene:{old_g2}
            """,
        )

        config = StableIdEventConfig(
            ref_gff=ref_gff,
            target_gff=target_gff,
            mapped_gff=mapped_gff,
            report=report,
            db_name="synthetic_core",
            mapping_session_id=42,
            gene_range=parse_id_range("ENSNWIG:5000001-5000010"),
            transcript_range=parse_id_range("ENSNWIT:5000001-5000010"),
            translation_range=parse_id_range("ENSNWIP:5000001-5000010"),
            output_sql=output_sql,
            output_tsv=output_tsv,
            include_translations=True,
            dry_run=True,
            backup_prefix="stable_id_mapper_backup_test",
            batch_size=2,
            replace_events_for_session=False,
            min_overlap=0.10,
        )
        decisions = run_pipeline(config)

        expected_counts = {
            ("gene", "mapped"): 2,
            ("gene", "missing"): 1,
            ("gene", "new"): 1,
            ("transcript", "mapped"): 2,
            ("transcript", "missing"): 1,
            ("transcript", "new"): 1,
            ("translation", "mapped"): 1,
            ("translation", "missing"): 1,
            ("translation", "new"): 1,
        }
        observed_counts: dict[tuple[str, str], int] = {}
        for decision in decisions:
            key = (decision.feature_type, decision.action)
            observed_counts[key] = observed_counts.get(key, 0) + 1
        if observed_counts != expected_counts:
            raise AssertionError(
                f"Unexpected decision counts: {observed_counts} != {expected_counts}"
            )

        assert_has_decision(
            decisions,
            "gene",
            "mapped",
            current_stable_id=target_g1,
            old_stable_id=old_g1,
            new_stable_id=old_g1,
            new_version=2,
        )
        assert_has_decision(decisions, "gene", "missing", old_stable_id=old_g2)
        assert_has_decision(decisions, "gene", "new", current_stable_id=target_g3)
        assert_has_decision(
            decisions,
            "transcript",
            "mapped",
            current_stable_id=target_t1,
            old_stable_id=old_t1,
            new_stable_id=old_t1,
            new_version=2,
        )
        assert_has_decision(
            decisions,
            "translation",
            "mapped",
            current_stable_id=target_t1,
            old_stable_id=old_t1,
            new_stable_id=old_t1,
            new_version=2,
        )

        sql_text = output_sql.read_text(encoding="utf-8")
        if "USE `synthetic_core`;" not in sql_text:
            raise AssertionError("SQL should include db_name when provided")
        if "CREATE TEMPORARY TABLE tmp_stable_id_mapper_decisions" not in sql_text:
            raise AssertionError("Dry-run decision table was not written")
        if "ROLLBACK;" not in sql_text:
            raise AssertionError("Dry-run SQL should end with ROLLBACK")
        tsv_header = output_tsv.read_text(encoding="utf-8").splitlines()[0]
        if not tsv_header.startswith("type\taction\tcurrent_stable_id"):
            raise AssertionError("TSV header is incorrect")

    print("All tests passed.")


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Ensembl core SQL from stable-ID mapped GFF3 output."
    )
    parser.add_argument("--test", action="store_true", help="run built-in synthetic tests")
    parser.add_argument("--ref-gff", type=Path)
    parser.add_argument("--target-gff", type=Path)
    parser.add_argument("--mapped-gff", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--db-name")
    parser.add_argument("--mapping-session-id", type=int)
    parser.add_argument("--gene-range", type=parse_id_range)
    parser.add_argument("--transcript-range", type=parse_id_range)
    parser.add_argument("--translation-range", type=parse_id_range)
    parser.add_argument("--output-sql", type=Path)
    parser.add_argument("--output-tsv", type=Path)
    parser.add_argument("--include-translations", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--backup-prefix", default=default_backup_prefix())
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--replace-events-for-session", action="store_true")
    parser.add_argument("--min-overlap", type=float, default=0.10)

    args = parser.parse_args(argv)
    if args.test:
        return args

    required = [
        ("ref_gff", "--ref-gff"),
        ("target_gff", "--target-gff"),
        ("mapped_gff", "--mapped-gff"),
        ("report", "--report"),
        ("mapping_session_id", "--mapping-session-id"),
        ("gene_range", "--gene-range"),
        ("transcript_range", "--transcript-range"),
        ("translation_range", "--translation-range"),
        ("output_sql", "--output-sql"),
    ]
    missing = [flag for attr, flag in required if getattr(args, attr) is None]
    if missing:
        parser.error("the following arguments are required: " + ", ".join(missing))
    if args.batch_size < 1:
        parser.error("--batch-size must be >= 1")
    if args.min_overlap < 0:
        parser.error("--min-overlap must be >= 0")
    return args


def config_from_args(args: argparse.Namespace) -> StableIdEventConfig:
    return StableIdEventConfig(
        ref_gff=args.ref_gff,
        target_gff=args.target_gff,
        mapped_gff=args.mapped_gff,
        report=args.report,
        db_name=args.db_name,
        mapping_session_id=args.mapping_session_id,
        gene_range=args.gene_range,
        transcript_range=args.transcript_range,
        translation_range=args.translation_range,
        output_sql=args.output_sql,
        output_tsv=args.output_tsv,
        include_translations=args.include_translations,
        dry_run=args.dry_run,
        backup_prefix=args.backup_prefix,
        batch_size=args.batch_size,
        replace_events_for_session=args.replace_events_for_session,
        min_overlap=args.min_overlap,
    )


def main() -> None:
    args = parse_args()
    if args.test:
        run_tests()
        return
    run_pipeline(config_from_args(args))


if __name__ == "__main__":
    main()


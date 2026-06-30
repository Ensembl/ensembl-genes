"""Tests for the Ensembl GFF/GTF loading package."""

# pylint: disable=missing-function-docstring,missing-class-docstring
# pylint: disable=too-many-arguments,too-many-instance-attributes
# pylint: disable=too-many-branches,too-many-statements
# pylint: disable=wrong-import-position

from __future__ import annotations

import importlib
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pytest

SRC_PATH = Path(__file__).resolve().parents[1] / "src" / "python"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from ensembl.genes.ensembl_loading import (
    gff_cli,
    gff_core_loader,
    gff_repeat_loader,
)
from ensembl.genes.ensembl_loading.gff_quality_check import (
    expected_translation_stable_ids,
)
from ensembl.genes.ensembl_loading.gff_source_config import (
    available_source_configs,
    get_source_config,
)
from ensembl.genes.ensembl_loading.refseq_conversion import (
    convert_gff_to_ensembl,
    default_gff_output_path,
    load_refseq_name_map,
)
from ensembl.genes.ensembl_loading.refseq_ncbi import (
    accession_subdir,
    parse_assembly_summary,
)


def feature_row(
    seqid: str,
    feature_type: str,
    start: int,
    end: int,
    strand: str,
    attributes: str,
    source: str = "ensembl",
    phase: str = ".",
) -> str:
    return "\t".join(
        (
            seqid,
            source,
            feature_type,
            str(start),
            str(end),
            ".",
            strand,
            phase,
            attributes,
        )
    )


def write_lines(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n")
    return path


def test_public_modules_import_and_cli_exposes_new_package() -> None:
    modules = (
        "gff_cli",
        "gff_core_loader",
        "gff_models",
        "gff_quality_check",
        "gff_repeat_loader",
        "gff_source_config",
        "refseq_constants",
        "refseq_conversion",
        "refseq_models",
        "refseq_ncbi",
    )

    for module in modules:
        imported = importlib.import_module(f"ensembl.genes.ensembl_loading.{module}")
        assert imported.__name__.endswith(module)

    args = gff_cli.build_parser().parse_args(
        [
            "load-features",
            "annotation.gtf",
            "--db-name",
            "example_core",
            "--db-host",
            "mysql-host",
            "--db-user",
            "ensro",
            "--db-password",
            "secret",
            "--db-port",
            "3306",
            "--source",
            "anno_gtf",
        ]
    )

    assert args.func is gff_cli.run_load_features
    assert args.source == "anno_gtf"
    assert "anno_gtf" in available_source_configs()
    assert "ncrna_gtf" in available_source_configs()

    ncrna_args = gff_cli.build_parser().parse_args(
        [
            "load-features",
            "rfam_output/annotation.gtf",
            "--db-name",
            "example_core",
            "--db-host",
            "mysql-host",
            "--db-user",
            "ensro",
            "--db-password",
            "secret",
            "--db-port",
            "3306",
            "--source",
            "ncrna_gtf",
        ]
    )

    assert ncrna_args.func is gff_cli.run_load_features
    assert ncrna_args.source == "ncrna_gtf"

    repeat_args = gff_cli.build_parser().parse_args(
        [
            "load-single-line-features",
            "dust_output/annotation.gtf",
            "--analysis-name",
            "dust",
            "--db-name",
            "example_core",
            "--db-host",
            "mysql-host",
            "--db-user",
            "ensro",
            "--db-password",
            "secret",
            "--db-port",
            "3306",
        ]
    )

    assert repeat_args.func is gff_cli.run_load_single_line_features
    assert repeat_args.analysis_name == "dust"

    anno_output_args = gff_cli.build_parser().parse_args(
        [
            "load-anno-output",
            "GCA_000000000.1",
            "--db-name",
            "example_core",
            "--db-host",
            "mysql-host",
            "--db-user",
            "ensro",
            "--db-password",
            "secret",
            "--db-port",
            "3306",
        ]
    )

    assert anno_output_args.func is gff_cli.run_load_anno_output
    assert anno_output_args.output_dir == Path("GCA_000000000.1")


def test_source_specific_annotation_parsers(tmp_path: Path) -> None:
    anno_gtf = write_lines(
        tmp_path / "annotation.gtf",
        [
            feature_row(
                "17",
                "transcript",
                50,
                80,
                "+",
                'gene_id "gene1"; transcript_id "tx2";',
            ),
            feature_row(
                "17",
                "exon",
                50,
                80,
                "+",
                'gene_id "gene1"; transcript_id "tx2"; exon_number "1";',
            ),
            feature_row(
                "17",
                "transcript",
                100,
                500,
                "+",
                'gene_id "gene1"; transcript_id "tx1"; '
                'translation_coords "100:150:11:400:500:51";',
            ),
            feature_row(
                "17",
                "exon",
                100,
                150,
                "+",
                'gene_id "gene1"; transcript_id "tx1"; exon_number "1";',
            ),
            feature_row(
                "17",
                "exon",
                200,
                300,
                "+",
                'gene_id "gene1"; transcript_id "tx1"; exon_number "2";',
            ),
            feature_row(
                "17",
                "exon",
                400,
                500,
                "+",
                'gene_id "gene1"; transcript_id "tx1"; exon_number "3";',
            ),
        ],
    )
    ncrna_gtf = write_lines(
        tmp_path / "rfam_annotation.gtf",
        [
            feature_row(
                "17",
                "transcript",
                900,
                950,
                "+",
                'gene_id "rfam_gene1"; transcript_id "rfam_gene1.t1"; '
                'biotype "snoRNA";',
                source="Rfam",
            ),
            feature_row(
                "17",
                "exon",
                900,
                950,
                "+",
                'gene_id "rfam_gene1"; transcript_id "rfam_gene1.t1";',
                source="Rfam",
            ),
            feature_row(
                "17",
                "transcript",
                1000,
                1100,
                "-",
                'gene_id "rfam_gene2"; transcript_id "rfam_gene2.t1"; biotype "rRNA";',
                source="Rfam",
            ),
            feature_row(
                "17",
                "exon",
                1000,
                1100,
                "-",
                'gene_id "rfam_gene2"; transcript_id "rfam_gene2.t1";',
                source="Rfam",
            ),
            feature_row(
                "17",
                "transcript",
                1200,
                1250,
                "+",
                'gene_id "rfam_gene3"; transcript_id "rfam_gene3.t1";',
                source="Rfam",
            ),
            feature_row(
                "17",
                "exon",
                1200,
                1250,
                "+",
                'gene_id "rfam_gene3"; transcript_id "rfam_gene3.t1";',
                source="Rfam",
            ),
        ],
    )
    ensembl_gff = write_lines(
        tmp_path / "ensembl.gff3",
        [
            "##gff-version 3",
            feature_row(
                "20",
                "gene",
                100,
                900,
                "+",
                "ID=gene:ENSG000001;Name=ABC1;biotype=protein_coding;"
                "gene_id=ENSG000001",
            ),
            feature_row(
                "20",
                "mRNA",
                100,
                900,
                "+",
                "ID=transcript:ENST000001;Parent=gene:ENSG000001;"
                "Name=ABC1-201;biotype=protein_coding;"
                "transcript_id=ENST000001",
            ),
            feature_row(
                "20",
                "exon",
                100,
                200,
                "+",
                "Parent=transcript:ENST000001;Name=ENSE000001;"
                "exon_id=ENSE000001;rank=1",
            ),
            feature_row(
                "20",
                "CDS",
                120,
                200,
                "+",
                "Parent=transcript:ENST000001;protein_id=ENSP000001",
                phase="0",
            ),
        ],
    )
    refseq_gff = write_lines(
        tmp_path / "refseq.gff3",
        [
            "##gff-version 3",
            feature_row(
                "chr1",
                "gene",
                100,
                500,
                "+",
                "ID=gene-GeneA;Name=GeneA;gbkey=Gene",
                source="RefSeq",
            ),
            feature_row(
                "chr1",
                "mRNA",
                100,
                500,
                "+",
                "ID=rna-TxA;Parent=gene-GeneA;Name=TxA;gbkey=mRNA",
                source="RefSeq",
            ),
            feature_row(
                "chr1",
                "exon",
                100,
                200,
                "+",
                "ID=exon-TxA-1;Parent=rna-TxA;gbkey=mRNA",
                source="RefSeq",
            ),
            feature_row(
                "chr1",
                "CDS",
                150,
                200,
                "+",
                "ID=cds-ProtA;Parent=rna-TxA;protein_id=ProtA",
                source="RefSeq",
                phase="0",
            ),
        ],
    )

    anno = gff_core_loader.prepare_annotation_for_load(
        anno_gtf,
        source_config=get_source_config("anno_gtf"),
    )
    ncrna_config = get_source_config("ncrna_gtf")
    ncrna = gff_core_loader.prepare_annotation_for_load(
        ncrna_gtf,
        source_config=ncrna_config,
    )
    ensembl = gff_core_loader.prepare_annotation_for_load(
        ensembl_gff,
        source_config=get_source_config("ensembl"),
    )
    refseq = gff_core_loader.prepare_annotation_for_load(
        refseq_gff,
        source_config=get_source_config("refseq"),
    )

    assert anno.genes["gene1"].start == 50
    assert anno.genes["gene1"].end == 500
    assert anno.transcripts["tx1"].biotype == "protein_coding"
    assert anno.transcripts["tx2"].biotype == "not_set"
    assert [(cds.start, cds.end, cds.phase) for cds in anno.cds_segments["tx1"]] == [
        (110, 150, "0"),
        (200, 300, "."),
        (400, 450, "."),
    ]

    assert ncrna_config.source_label == "ensembl"
    assert ncrna_config.analysis_logic_name == "ncrna"
    assert ncrna.genes["rfam_gene1"].biotype == "snoRNA"
    assert ncrna.transcripts["rfam_gene1.t1"].biotype == "snoRNA"
    assert ncrna.genes["rfam_gene2"].biotype == "rRNA"
    assert ncrna.transcripts["rfam_gene2.t1"].biotype == "rRNA"
    assert ncrna.genes["rfam_gene3"].biotype == "misc_RNA"
    assert ncrna.transcripts["rfam_gene3.t1"].biotype == "misc_RNA"
    assert {
        transcript_id: [(exon.start, exon.end) for exon in transcript.exons]
        for transcript_id, transcript in ncrna.transcripts.items()
    } == {
        "rfam_gene1.t1": [(900, 950)],
        "rfam_gene2.t1": [(1000, 1100)],
        "rfam_gene3.t1": [(1200, 1250)],
    }

    assert set(ensembl.genes) == {"ENSG000001"}
    assert set(ensembl.transcripts) == {"ENST000001"}
    assert ensembl.transcripts["ENST000001"].protein_id == "ENSP000001"
    assert expected_translation_stable_ids(ensembl) == {"ENST000001": "ENSP000001"}

    assert set(refseq.genes) == {"GeneA"}
    assert set(refseq.transcripts) == {"TxA"}
    assert refseq.transcripts["TxA"].stable_id == "TxA"
    assert expected_translation_stable_ids(refseq) == {"TxA": "TxA_prot"}


def test_refseq_discovery_and_conversion_helpers(tmp_path: Path) -> None:
    assert accession_subdir("GCF_000001635.27") == "000/001/635"
    assert default_gff_output_path("GCF_000001635.27_genomic.gff.gz").name == (
        "GCF_000001635.27_genomic_ensembl.gff3"
    )
    assert (
        gff_core_loader.derive_core_db_name(
            "Mus musculus",
            "GCF_000001635.27",
        )
        == "mus_musculus_core_000001635_27"
    )
    assert (
        gff_core_loader.derive_core_db_name(
            "Scientific name",
            "GCF_037462849.1",
            source_config=get_source_config("refseq"),
        )
        == "scientific_name_gcf037462849v1_rs_core_114_1"
    )

    report = write_lines(
        tmp_path / "assembly_report.txt",
        [
            "# header",
            "chrOne\tassembled-molecule\t1\tChromosome\tna\tna\tNC_000001.11",
            "scaffoldA\tunlocalized-scaffold\tna\tna\tna\tna\tNW_000001.1",
        ],
    )
    assert load_refseq_name_map(report) == {
        "NC_000001.11": "1",
        "NW_000001.1": "scaffoldA",
    }

    gff = write_lines(
        tmp_path / "input.gff3",
        [
            "##gff-version 3",
            "##sequence-region NC_000001.11 1 1000",
            feature_row(
                "NC_000001.11",
                "gene",
                10,
                50,
                "+",
                "ID=gene-GeneA",
                source="RefSeq",
            ),
        ],
    )
    converted = convert_gff_to_ensembl(gff, report, tmp_path / "converted.gff3")
    converted_text = converted.read_text(encoding="utf-8")
    assert "##sequence-region 1 1 1000" in converted_text
    assert "\n1\tRefSeq\tgene\t10\t50" in converted_text

    columns = ["x"] * 20
    columns[0] = "GCF_000001635.27"
    columns[5] = "10090"
    columns[7] = "Mus musculus strain example"
    columns[10] = "latest"
    columns[15] = "GRCm39"
    columns[19] = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/"
        "GCF_000001635.27_GRCm39"
    )
    records = parse_assembly_summary("\t".join(columns), "vertebrate_mammalian")
    assert len(records) == 1
    assert records[0].species_name == "Mus musculus"
    assert records[0].paths.gff_url.endswith("_genomic.gff.gz")


class FakeCoreCursor:
    def __init__(self) -> None:
        self.lastrowid = 0
        self.fetchone_result: tuple[Any, ...] | None = None
        self.fetchall_result: list[tuple[Any, ...]] = []
        self.seq_regions = {"17": 101}
        self.genes: dict[int, tuple[str, str]] = {}
        self.transcripts: dict[int, tuple[str, str]] = {}
        self.exons: dict[int, tuple[str | None]] = {}
        self.exon_transcripts: set[tuple[int, int]] = set()
        self.translations: dict[int, str | None] = {}
        self.analysis_inserted: tuple[Any, ...] | None = None

    def execute(self, operation: str, params: Sequence[Any] | None = None) -> None:
        sql = " ".join(operation.lower().split())
        if sql.startswith("select coord_system_id from coord_system"):
            self.fetchone_result = (1,)
        elif sql.startswith("select analysis_id from analysis"):
            self.fetchone_result = None
        elif sql.startswith("insert into analysis"):
            assert params is not None
            self.analysis_inserted = tuple(params)
            self.lastrowid = 7
        elif sql.startswith("select seq_region_id from seq_region"):
            assert params is not None
            seq_name = params[0]
            seq_region_id = self.seq_regions.get(seq_name)
            self.fetchone_result = (seq_region_id,) if seq_region_id else None
        elif sql.startswith("select coalesce(max(gene_id)"):
            self.fetchone_result = (1,)
        elif sql.startswith("select coalesce(max(transcript_id)"):
            self.fetchone_result = (1,)
        elif sql.startswith("select coalesce(max(exon_id)"):
            self.fetchone_result = (1,)
        elif sql.startswith("insert into gene"):
            assert params is not None
            self.genes[int(params[0])] = (params[7], params[5])
        elif sql.startswith("insert into transcript"):
            assert params is not None
            self.transcripts[int(params[0])] = (
                params[8],
                params[6],
            )
        elif sql.startswith("insert into exon "):
            assert params is not None
            self.exons[int(params[0])] = (params[7],)
        elif sql.startswith("insert into exon_transcript"):
            assert params is not None
            self.exon_transcripts.add((int(params[0]), int(params[1])))
        elif sql.startswith("insert into translation"):
            assert params is not None
            self.lastrowid = len(self.translations) + 1
            self.translations[int(params[0])] = params[5]
        elif sql.startswith("update transcript set canonical_translation_id"):
            pass
        elif sql.startswith("select gene_id, stable_id, biotype from gene"):
            self.fetchall_result = [
                (db_id, *values) for db_id, values in sorted(self.genes.items())
            ]
        elif sql.startswith("select transcript_id, stable_id, biotype from transcript"):
            self.fetchall_result = [
                (db_id, *values) for db_id, values in sorted(self.transcripts.items())
            ]
        elif sql.startswith("select exon_id, stable_id from exon"):
            self.fetchall_result = [
                (db_id, *values) for db_id, values in sorted(self.exons.items())
            ]
        elif sql.startswith("select exon_id, transcript_id from exon_transcript"):
            self.fetchall_result = sorted(self.exon_transcripts)
        elif sql.startswith("select transcript_id, stable_id from translation"):
            self.fetchall_result = list(sorted(self.translations.items()))
        else:  # pragma: no cover - makes unsupported SQL obvious in failures.
            raise AssertionError(f"Unexpected SQL: {operation}")

    def executemany(self, operation: str, seq_params: list[Any]) -> None:
        for params in seq_params:
            self.execute(operation, params)

    def fetchone(self) -> tuple[Any, ...] | None:
        return self.fetchone_result

    def fetchall(self) -> list[tuple[Any, ...]]:
        return self.fetchall_result


class FakeCoreConnection:
    def __init__(self) -> None:
        self.cursor_instance = FakeCoreCursor()
        self.committed = False
        self.rolled_back = False
        self.closed = False

    def cursor(self) -> FakeCoreCursor:
        return self.cursor_instance

    def commit(self) -> None:
        self.committed = True

    def rollback(self) -> None:
        self.rolled_back = True

    def close(self) -> None:
        self.closed = True


class FakeSingleLineFeatureCursor:
    def __init__(self) -> None:
        self.lastrowid = 0
        self.fetchone_result: tuple[Any, ...] | None = None
        self.seq_regions = {"17": 101}
        self.analysis_inserted: tuple[Any, ...] | None = None
        self.repeat_consensus: dict[int, tuple[Any, ...]] = {}
        self.repeat_features: list[tuple[Any, ...]] = []
        self.simple_features: list[tuple[Any, ...]] = []

    def execute(self, operation: str, params: Sequence[Any] | None = None) -> None:
        sql = " ".join(operation.lower().split())
        if sql.startswith("select coord_system_id from coord_system"):
            self.fetchone_result = (1,)
        elif sql.startswith("select analysis_id from analysis"):
            self.fetchone_result = None
        elif sql.startswith("insert into analysis"):
            assert params is not None
            self.analysis_inserted = tuple(params)
            self.lastrowid = 7
        elif sql.startswith("select seq_region_id from seq_region"):
            assert params is not None
            seq_region_id = self.seq_regions.get(params[0])
            self.fetchone_result = (seq_region_id,) if seq_region_id else None
        elif sql.startswith("select repeat_consensus_id from repeat_consensus"):
            self.fetchone_result = None
        elif sql.startswith("insert into repeat_consensus"):
            assert params is not None
            self.lastrowid = len(self.repeat_consensus) + 1
            self.repeat_consensus[self.lastrowid] = tuple(params)
        elif sql.startswith("insert into repeat_feature"):
            assert params is not None
            self.repeat_features.append(tuple(params))
        elif sql.startswith("insert into simple_feature"):
            assert params is not None
            self.simple_features.append(tuple(params))
        else:  # pragma: no cover - makes unsupported SQL obvious in failures.
            raise AssertionError(f"Unexpected SQL: {operation}")

    def fetchone(self) -> tuple[Any, ...] | None:
        return self.fetchone_result


class FakeSingleLineFeatureConnection:
    def __init__(self) -> None:
        self.cursor_instance = FakeSingleLineFeatureCursor()
        self.committed = False
        self.rolled_back = False
        self.closed = False

    def cursor(self) -> FakeSingleLineFeatureCursor:
        return self.cursor_instance

    def commit(self) -> None:
        self.committed = True

    def rollback(self) -> None:
        self.rolled_back = True

    def close(self) -> None:
        self.closed = True


def test_single_line_feature_loader_loads_repeats_and_simple_features(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repeat_gtf = write_lines(
        tmp_path / "dust.gtf",
        [
            feature_row(
                "17",
                "low_complexity",
                10,
                20,
                "+",
                'repeat_name "Low_complexity"; repeat_class "dust"; '
                'repeat_type "Dust"; repeat_consensus "NNN"; score "42";',
                source="dust",
            ),
        ],
    )
    repeat_connection = FakeSingleLineFeatureConnection()
    monkeypatch.setattr(
        gff_repeat_loader,
        "connect_mysql",
        lambda **_: repeat_connection,
    )

    repeat_summary = gff_repeat_loader.load_single_line_features_to_core(
        gtf_path=repeat_gtf,
        analysis_name="dust",
        db_name="existing_core",
        db_host="mysql-host",
        db_user="ensro",
        db_password="secret",
        db_port=3306,
    )

    assert repeat_summary == {
        "repeat_features": 1,
        "repeat_consensus": 1,
        "simple_features": 0,
    }
    assert repeat_connection.committed is True
    assert repeat_connection.rolled_back is False
    assert repeat_connection.closed is True
    assert repeat_connection.cursor_instance.analysis_inserted == ("dust", "Anno")
    assert repeat_connection.cursor_instance.repeat_consensus == {
        1: ("Low_complexity", "dust", "Dust", "NNN")
    }
    assert repeat_connection.cursor_instance.repeat_features == [
        (101, 10, 20, 1, 1, 11, 1, 7, 42.0)
    ]

    simple_gtf = write_lines(
        tmp_path / "cpg.gtf",
        [
            feature_row(
                "17",
                "CpG",
                30,
                40,
                "-",
                'feature_id "cpg1";',
                source="cpg",
            ),
        ],
    )
    simple_connection = FakeSingleLineFeatureConnection()
    monkeypatch.setattr(
        gff_repeat_loader,
        "connect_mysql",
        lambda **_: simple_connection,
    )

    simple_summary = gff_repeat_loader.load_single_line_features_to_core(
        gtf_path=simple_gtf,
        analysis_name="cpg",
        db_name="existing_core",
        db_host="mysql-host",
        db_user="ensro",
        db_password="secret",
        db_port=3306,
    )

    assert simple_summary == {
        "repeat_features": 0,
        "repeat_consensus": 0,
        "simple_features": 1,
    }
    assert simple_connection.committed is True
    assert simple_connection.cursor_instance.analysis_inserted == ("cpg", "Anno")
    assert simple_connection.cursor_instance.simple_features == [
        (101, 30, 40, -1, "", 7, 0.0)
    ]


def test_load_anno_output_runs_available_loaders(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_dir = tmp_path / "GCA_000000000.1"
    for relative_dir in (
        "annotation_output",
        "rfam_output",
        "trnascan_output",
        "dust_output",
        "cpg_output",
    ):
        (output_dir / relative_dir).mkdir(parents=True)

    write_lines(
        output_dir / "annotation_output" / "annotation.gtf",
        [
            feature_row(
                "17",
                "transcript",
                10,
                20,
                "+",
                'gene_id "gene1"; transcript_id "tx1";',
            )
        ],
    )
    write_lines(
        output_dir / "rfam_output" / "annotation.gtf",
        [
            feature_row(
                "17",
                "transcript",
                30,
                40,
                "+",
                'gene_id "rfam1"; transcript_id "rfam1.t1"; biotype "snRNA";',
                source="Rfam",
            )
        ],
    )
    (output_dir / "trnascan_output" / "annotation.gtf").write_text(
        "",
        encoding="utf-8",
    )
    write_lines(
        output_dir / "dust_output" / "annotation.gtf",
        [
            feature_row(
                "17",
                "low_complexity",
                50,
                60,
                "+",
                'repeat_name "dust";',
                source="dust",
            )
        ],
    )
    write_lines(
        output_dir / "cpg_output" / "annotation.gtf",
        [
            feature_row(
                "17",
                "CpG",
                70,
                80,
                "+",
                'feature_id "cpg1";',
                source="cpg",
            )
        ],
    )

    calls: list[tuple[str, str, str]] = []

    def fake_load_gff_features_to_core(**kwargs: Any) -> dict[str, int]:
        gff_path = Path(kwargs["gff_path"])
        calls.append(
            (
                "feature",
                str(gff_path.relative_to(output_dir)),
                kwargs["source_config"].name,
            )
        )
        return {"genes": 2, "transcripts": 3, "cds_transcript_groups": 1}

    def fake_load_single_line_features_to_core(**kwargs: Any) -> dict[str, int]:
        gtf_path = Path(kwargs["gtf_path"])
        calls.append(
            (
                "single_line",
                str(gtf_path.relative_to(output_dir)),
                kwargs["analysis_name"],
            )
        )
        return {
            "repeat_features": 1 if kwargs["analysis_name"] == "dust" else 0,
            "repeat_consensus": 1 if kwargs["analysis_name"] == "dust" else 0,
            "simple_features": 1 if kwargs["analysis_name"] == "cpg" else 0,
        }

    monkeypatch.setattr(
        gff_cli,
        "load_gff_features_to_core",
        fake_load_gff_features_to_core,
    )
    monkeypatch.setattr(
        gff_cli,
        "load_single_line_features_to_core",
        fake_load_single_line_features_to_core,
    )

    args = gff_cli.build_parser().parse_args(
        [
            "load-anno-output",
            str(output_dir),
            "--db-name",
            "existing_core",
            "--db-host",
            "mysql-host",
            "--db-user",
            "ensro",
            "--db-password",
            "secret",
            "--db-port",
            "3306",
        ]
    )

    assert args.func(args) == 0
    assert calls == [
        ("feature", "annotation_output/annotation.gtf", "anno_gtf"),
        ("feature", "rfam_output/annotation.gtf", "ncrna_gtf"),
        ("single_line", "dust_output/annotation.gtf", "dust"),
        ("single_line", "cpg_output/annotation.gtf", "cpg"),
    ]


def test_anno_gtf_loads_into_existing_core_with_quality_check(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    gtf = write_lines(
        tmp_path / "annotation.gtf",
        [
            feature_row(
                "17",
                "transcript",
                100,
                300,
                "+",
                'gene_id "gene1"; transcript_id "tx1"; '
                'translation_coords "100:150:1:200:300:51";',
            ),
            feature_row(
                "17",
                "exon",
                100,
                150,
                "+",
                'gene_id "gene1"; transcript_id "tx1"; exon_number "1";',
            ),
            feature_row(
                "17",
                "exon",
                200,
                300,
                "+",
                'gene_id "gene1"; transcript_id "tx1"; exon_number "2";',
            ),
        ],
    )
    connection = FakeCoreConnection()
    monkeypatch.setattr(gff_core_loader, "connect_mysql", lambda **_: connection)

    summary = gff_core_loader.load_gff_features_to_core(
        gff_path=gtf,
        db_name="existing_core",
        db_host="mysql-host",
        db_user="ensro",
        db_password="secret",
        db_port=3306,
        source_config=get_source_config("anno_gtf"),
    )

    assert summary == {
        "genes": 1,
        "transcripts": 1,
        "cds_transcript_groups": 1,
    }
    assert connection.committed is True
    assert connection.rolled_back is False
    assert connection.closed is True
    assert connection.cursor_instance.analysis_inserted == ("ensembl", "Anno_GTF")
    assert connection.cursor_instance.genes == {1: ("gene1", "protein_coding")}
    assert connection.cursor_instance.transcripts == {1: ("tx1", "protein_coding")}
    assert connection.cursor_instance.exon_transcripts == {(1, 1), (2, 1)}
    assert connection.cursor_instance.translations == {1: "tx1_prot"}
    assert "GFF core quality check: PASS" in capsys.readouterr().out

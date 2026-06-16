from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import pytest

from ensembl.genes.ensembl_loading import gff_cli, gff_core_loader
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
        "gff_source_config",
        "refseq_constants",
        "refseq_conversion",
        "refseq_models",
        "refseq_ncbi",
    )

    for module in modules:
        imported = importlib.import_module(
            f"ensembl.genes.ensembl_loading.{module}"
        )
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
                'gene_id "rfam_gene2"; transcript_id "rfam_gene2.t1"; '
                'biotype "rRNA";',
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
    assert [
        (cds.start, cds.end, cds.phase) for cds in anno.cds_segments["tx1"]
    ] == [
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
    assert expected_translation_stable_ids(ensembl) == {
        "ENST000001": "ENSP000001"
    }

    assert set(refseq.genes) == {"GeneA"}
    assert set(refseq.transcripts) == {"TxA"}
    assert refseq.transcripts["TxA"].stable_id == "TxA"
    assert expected_translation_stable_ids(refseq) == {"TxA": "TxA_prot"}


def test_refseq_discovery_and_conversion_helpers(tmp_path: Path) -> None:
    assert accession_subdir("GCF_000001635.27") == "000/001/635"
    assert default_gff_output_path("GCF_000001635.27_genomic.gff.gz").name == (
        "GCF_000001635.27_genomic_ensembl.gff3"
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
    assert "##sequence-region 1 1 1000" in converted.read_text()
    assert "\n1\tRefSeq\tgene\t10\t50" in converted.read_text()

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

    def execute(self, operation: str, params: Any | None = None) -> None:
        sql = " ".join(operation.lower().split())
        if sql.startswith("select coord_system_id from coord_system"):
            self.fetchone_result = (1,)
        elif sql.startswith("select analysis_id from analysis"):
            self.fetchone_result = None
        elif sql.startswith("insert into analysis"):
            self.analysis_inserted = tuple(params)
            self.lastrowid = 7
        elif sql.startswith("select seq_region_id from seq_region"):
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
            self.genes[int(params[0])] = (params[7], params[5])
        elif sql.startswith("insert into transcript"):
            self.transcripts[int(params[0])] = (params[8], params[6])
        elif sql.startswith("insert into exon "):
            self.exons[int(params[0])] = (params[7],)
        elif sql.startswith("insert into exon_transcript"):
            self.exon_transcripts.add((int(params[0]), int(params[1])))
        elif sql.startswith("insert into translation"):
            self.lastrowid = len(self.translations) + 1
            self.translations[int(params[0])] = params[5]
        elif sql.startswith("update transcript set canonical_translation_id"):
            pass
        elif sql.startswith("select gene_id, stable_id, biotype from gene"):
            self.fetchall_result = [
                (db_id, *values) for db_id, values in sorted(self.genes.items())
            ]
        elif sql.startswith(
            "select transcript_id, stable_id, biotype from transcript"
        ):
            self.fetchall_result = [
                (db_id, *values)
                for db_id, values in sorted(self.transcripts.items())
            ]
        elif sql.startswith("select exon_id, stable_id from exon"):
            self.fetchall_result = [
                (db_id, *values) for db_id, values in sorted(self.exons.items())
            ]
        elif sql.startswith("select exon_id, transcript_id from exon_transcript"):
            self.fetchall_result = sorted(self.exon_transcripts)
        elif sql.startswith("select transcript_id, stable_id from translation"):
            self.fetchall_result = [
                (db_id, stable_id)
                for db_id, stable_id in sorted(self.translations.items())
            ]
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

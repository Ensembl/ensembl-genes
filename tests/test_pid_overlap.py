from pathlib import Path
from ensembl.genes.stable_id.pid_overlap import extract_proteins_from_gff3
import pandas as pd


def test_extract_proteins_from_gff3(tmp_path):
    fake_file: Path = tmp_path / "dummy.gff3"

    fake_file.write_text(
        "chr1\tDUMMY\tCDS\t100\t200\t.\t+\t.\tID=CDS:ENSP00000354687;\n"
    )
    result_df = extract_proteins_from_gff3(str(fake_file))
    assert isinstance(result_df, pd.DataFrame)
    assert len(result_df) == 1


def test_extract_proteins_from_gff3_merges_records(tmp_path):
    fake_file: Path = tmp_path / "dummy.gff3"

    fake_file.write_text(
        "chr1\tDUMMY\tCDS\t100\t200\t.\t+\t.\tID=CDS:ENSP00000354687;\n"
        "chr1\tDUMMY\tCDS\t50\t250\t.\t+\t.\tID=CDS:ENSP00000354687;\n"
    )
    result_df = extract_proteins_from_gff3(str(fake_file))
    assert result_df.iloc[0]["Start"] == 50
    assert result_df.iloc[0]["End"] == 250

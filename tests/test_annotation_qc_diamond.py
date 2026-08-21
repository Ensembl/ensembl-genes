from ensembl.genes.annotation_qc.parsers.diamond import parse_diamond_hits


def test_parse_diamond_hits_normalizes_gffread_transcript_prefix(tmp_path):
    hits_file = tmp_path / "diamond.tsv"
    hits_file.write_text(
        "transcript:ENSHGLT00000000139\tQ2KHR2.1\t91.9\t1364\t89\t1449\t1449\t93.9\t1\t1363\t1363\t0.0\t2320\n"
        "ENSHGLT00000000143\tQ62940.1\t87.4\t174\t1\t174\t201\t86.6\t108\t280\t887\t3.89e-102\t315\n"
    )

    hits = parse_diamond_hits(str(hits_file))

    assert hits["qseqid"].tolist() == [
        "ENSHGLT00000000139",
        "ENSHGLT00000000143",
    ]

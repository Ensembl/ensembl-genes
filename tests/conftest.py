import sys
from pathlib import Path

# Make metrics/parsers/runners importable as top-level packages.
sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "src/python/ensembl/genes/annotation-qc"),
)

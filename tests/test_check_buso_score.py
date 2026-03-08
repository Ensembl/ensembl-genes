import pytest
from ensembl.genes.metrics.check_busco_score import extract_completeness


@pytest.mark.parametrize("busco_str, expected", [
    ("C:99.2%[S:98.4%,D:0.8%],F:0.6%,M:0.2%,n:3640", 99.2),
    ("C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:100", 100.0),
    ("C:0.0%[S:0.0%,D:0.0%],F:50.0%,M:50.0%,n:10", 0.0),
    ("C:100%[S:100%,D:0%],F:0%,M:0%,n:10", 100)
])
def test_extract_completeness_valid(busco_str, expected):
    assert extract_completeness(busco_str) == expected


def test_extract_completeness_invalid():
    """Test that invalid BUSCO strings raise ValueError."""
    busco_str = "invalid_busco_string"
    with pytest.raises(ValueError, match=f"Invalid BUSCO format: {busco_str}"):
        extract_completeness(busco_str)


def test_extract_completeness_empty():
    """Test that empty BUSCO strings raise ValueError."""
    busco_str = ""
    with pytest.raises(ValueError, match=f"Invalid BUSCO format: {busco_str}"):
        extract_completeness(busco_str)

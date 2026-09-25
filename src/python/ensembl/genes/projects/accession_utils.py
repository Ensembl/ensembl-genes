"""Utilities for converting INSDC assembly accessions to FTP directory paths."""

import re

# Matches GCA_NNNNNNNNN.V or GCF_NNNNNNNNN.V (exactly 9 digits, integer version).
_ACCESSION_RE = re.compile(r"^(GCA|GCF)_(\d{9})\.(\d+)$")


def accession_to_ftp_path(accession: str) -> str:
    """Convert an INSDC assembly accession to its new EBI FTP triplet path.

    The nine-digit accession body is split into three groups of three digits to
    keep directory sizes manageable.  Both GCA and GCF prefixes are supported,
    and leading zeroes in the digit body are preserved.

    Examples::

        accession_to_ftp_path("GCA_922984935.2")   # -> "GCA/922/984/935/2"
        accession_to_ftp_path("GCF_000001405.40")  # -> "GCF/000/001/405/40"
        accession_to_ftp_path("GCA_000001405.29")  # -> "GCA/000/001/405/29"

    Args:
        accession: An INSDC assembly accession string of the form
            ``GCA_NNNNNNNNN.V`` or ``GCF_NNNNNNNNN.V``.

    Returns:
        The FTP path component, e.g. ``GCA/922/984/935/2``.

    Raises:
        ValueError: If *accession* does not match the expected format (wrong
            prefix, wrong number of digits, or missing version).
    """
    m = _ACCESSION_RE.match(accession)
    if not m:
        raise ValueError(
            f"Malformed assembly accession {accession!r}. "
            "Expected format: GCA_NNNNNNNNN.V or GCF_NNNNNNNNN.V "
            "(nine-digit body, integer version after the dot)."
        )
    prefix, digits, version = m.group(1), m.group(2), m.group(3)
    return f"{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{version}"

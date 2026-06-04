"""NCBI assembly metadata helpers for annotation QC."""

import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile
from typing import Dict, Optional


def require_exe(name: str) -> str:
    """
    Ensure that a required executable is available on the system PATH.

    This function checks for the presence of an external command-line tool
    using ``shutil.which`` and raises a RuntimeError if the executable cannot
    be found. It is intended to fail fast before invoking subprocesses that
    depend on external tools.

    Parameters
    ----------
    name : str
        Name of the executable to locate (e.g. ``datasets``).

    Returns
    -------
    str
        Absolute path to the located executable.

    Raises
    ------
    RuntimeError
        If the executable is not found on the system PATH.
    """

    path = shutil.which(name)
    if not path:
        raise RuntimeError(
            f"Required executable '{name}' not found on PATH. "
            f"Install NCBI Datasets CLI and try again."
        )
    return path


def ncbi_assembly_stats(
    assembly_accession: str, outdir: Optional[str] = None
) -> Dict[str, object]:
    """
    Retrieve assembly-level statistics from NCBI using the Datasets CLI.

    This function downloads an assembly package for the given accession using
    the NCBI Datasets command-line tool, extracts the embedded
    ``assembly_data_report.jsonl`` file, and parses key assembly and organism
    statistics into a flat dictionary suitable for reporting.

    The function requires the ``datasets`` executable to be available on the
    system PATH and will raise an error if the command fails or if expected
    files are missing from the downloaded package.

    Parameters
    ----------
    assembly_accession : str
        NCBI assembly accession (e.g. ``GCF_000001405.40``).
    outdir : Optional[str], default=None
        Optional output directory for debugging artifacts. If provided, the
        downloaded ZIP file, extracted contents, command invocation, and
        parsed results are preserved for inspection.

    Returns
    -------
    Dict[str, object]
        Dictionary of assembly statistics, including:
        - Assembly accession, name, and release date
        - Scientific name and taxon ID
        - Contig N50
        - Total assembly length and total gap length
        - Chromosome and component sequence counts
        - GC percentage

        Numeric fields may be integers or empty strings depending on data
        availability in the NCBI report.

    Raises
    ------
    RuntimeError
        If the ``datasets`` executable is not available, the download command
        fails, the expected JSONL report cannot be found, or the report is
        empty or malformed.

    Notes
    -----
    - The NCBI Datasets JSON report may encode numeric fields as strings; this
      function attempts to normalize those values where possible.
    - Only the first record in ``assembly_data_report.jsonl`` is read, which
      is the expected behavior for single-accession downloads.
    - Some assembly metrics (e.g. spanned gaps) may not have a
      direct equivalent in the NCBI report and are left blank if unavailable.
    """
    require_exe("datasets")

    debug_dir = None
    if outdir:
        debug_dir = os.path.join(outdir, "ncbi_debug")
        os.makedirs(debug_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        zip_path = os.path.join(td, f"{assembly_accession}.ncbi_genome.zip")

        cmd = [
            "datasets",
            "download",
            "genome",
            "accession",
            assembly_accession,
            "--filename",
            zip_path,
        ]
        p = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )

        if debug_dir:
            with open(
                os.path.join(debug_dir, "datasets_cmd.txt"), "w", encoding="utf-8"
            ) as fh:
                fh.write(" ".join(cmd) + "\n")
            with open(
                os.path.join(debug_dir, "datasets_stdout.txt"), "w", encoding="utf-8"
            ) as fh:
                fh.write(p.stdout or "")
            with open(
                os.path.join(debug_dir, "datasets_stderr.txt"), "w", encoding="utf-8"
            ) as fh:
                fh.write(p.stderr or "")

        if p.returncode != 0:
            raise RuntimeError(f"datasets failed:\n{p.stderr}\nCMD: {' '.join(cmd)}")

        # Persist zip + unzip for inspection
        extracted_root = None
        if debug_dir:
            shutil.copy2(
                zip_path,
                os.path.join(debug_dir, f"{assembly_accession}.ncbi_genome.zip"),
            )
            extracted_root = os.path.join(debug_dir, "unzipped")
            if os.path.exists(extracted_root):
                shutil.rmtree(extracted_root)
            os.makedirs(extracted_root, exist_ok=True)
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(extracted_root)
            with zipfile.ZipFile(zip_path, "r") as zf:
                with open(
                    os.path.join(debug_dir, "zip_contents.txt"), "w", encoding="utf-8"
                ) as fh:
                    for n in zf.namelist():
                        fh.write(n + "\n")
        else:
            extracted_root = os.path.join(td, "unzipped")
            os.makedirs(extracted_root, exist_ok=True)
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(extracted_root)

        # Find assembly_data_report.jsonl inside extracted package
        jsonl_path = None
        for root, _, files in os.walk(extracted_root):
            for fn in files:
                if fn == "assembly_data_report.jsonl":
                    jsonl_path = os.path.join(root, fn)
                    break
            if jsonl_path:
                break

        if not jsonl_path or not os.path.exists(jsonl_path):
            raise RuntimeError(
                "Could not find assembly_data_report.jsonl after unzipping datasets package"
            )

        # Read first JSON object (JSONL => one object per line; usually just one line for accession download)
        with open(jsonl_path, "rt", encoding="utf-8", errors="replace") as fh:
            first = fh.readline().strip()
        if not first:
            raise RuntimeError("assembly_data_report.jsonl was empty")

        obj = json.loads(first)

        # Helper to navigate nested dicts safely
        def dig(d, path, default=""):
            cur = d
            for k in path:
                if not isinstance(cur, dict):
                    return default
                cur = cur.get(k)
                if cur is None:
                    return default
            return cur

        paired_acc = dig(obj, ["pairedAccession"], default="")  # e.g. GCF_907164915.1

        asm_name = dig(obj, ["assemblyInfo", "assemblyName"])
        asm_date = dig(obj, ["assemblyInfo", "releaseDate"])
        sci_name = dig(obj, ["organism", "organismName"])
        tax_id = dig(obj, ["organism", "taxId"])

        contig_n50 = dig(obj, ["assemblyStats", "contigN50"])
        total_len = dig(obj, ["assemblyStats", "totalSequenceLength"])
        ungapped_len = dig(obj, ["assemblyStats", "totalUngappedLength"])

        # total_len / ungapped_len may be strings in JSON
        def to_int(x):
            if x is None:
                return None
            if isinstance(x, int):
                return x
            s = str(x).strip()
            if not s:
                return None
            try:
                return int(s)
            except ValueError:
                return None

        total_len_i = to_int(total_len)
        ungapped_len_i = to_int(ungapped_len)

        total_gap_length = ""
        if (
            total_len_i is not None
            and ungapped_len_i is not None
            and total_len_i >= ungapped_len_i
        ):
            total_gap_length = str(total_len_i - ungapped_len_i)

        chromosome_count = dig(obj, ["assemblyStats", "totalNumberOfChromosomes"])
        component_count = dig(obj, ["assemblyStats", "numberOfComponentSequences"])
        gc_percent = dig(obj, ["assemblyStats", "gcPercent"])

        # NCBI report doesn’t always include a field that exactly equals “spanned gaps”
        # we can at least try the explicit key if present
        spanned_gaps = dig(
            obj, ["assemblyStats", "gapsBetweenScaffoldsCount"], default=""
        )

        result = {
            "assembly_accession": dig(obj, ["accession"], default=assembly_accession)
            or assembly_accession,
            "assembly_name": asm_name,
            "assembly_date": asm_date,
            "scientific_name": sci_name,
            "taxon_id": tax_id,
            "contig_n50": contig_n50,
            "total_length": total_len_i if total_len_i is not None else total_len,
            "total_gap_length": total_gap_length,
            "spanned_gaps": spanned_gaps,
            "molecule_count": chromosome_count,
            "component_count": component_count,
            "gc_percent": gc_percent,
        }

        if debug_dir:
            with open(
                os.path.join(debug_dir, "parsed_assembly_stats.json"),
                "w",
                encoding="utf-8",
            ) as fh:
                json.dump(result, fh, indent=2, sort_keys=True)

        # If "spanned gaps" is missing from datasets JSONL, fall back to legacy stats report (might want to fallback on this for other values, too... already added top_level_count)
        legacy = {}
        if not result.get("spanned_gaps"):
            report_acc = paired_acc or result["assembly_accession"]
            legacy = fetch_and_parse_ncbi_assembly_stats_report(
                accession=report_acc,
                assembly_name=result.get("assembly_name") or asm_name,
                outdir=outdir,
            )
            if legacy.get("spanned_gaps"):
                result["spanned_gaps"] = legacy["spanned_gaps"]
        if not result.get("toplevel_sequences"):
            if legacy.get("top_level_count"):
                result["toplevel_sequences"] = legacy["top_level_count"]

        return result


def ncbi_stats_report_url(accession: str, assembly_name: str) -> str:
    """
    Build the NCBI genomes/all URL for the legacy *_assembly_stats.txt report.

    Example:
      accession = GCF_907164915.1
      assembly_name = dImpGla2.1

    URL:
      https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/907/164/915/GCF_907164915.1_dImpGla2.1/GCF_907164915.1_dImpGla2.1_assembly_stats.txt
    """
    m = re.match(r"^(GCF|GCA)_(\d+)\.(\d+)$", accession)
    if not m:
        raise ValueError(f"Unexpected accession format: {accession}")

    prefix, digits, version = m.group(1), m.group(2), m.group(3)

    # NCBI groups digits into 3s: 907/164/915
    chunks = [digits[i : i + 3] for i in range(0, len(digits), 3)]
    if len(chunks) < 3:
        # pad defensively, though NCBI accessions are normally 9 digits here
        chunks = (chunks + ["000", "000", "000"])[:3]

    base = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
    dirname = f"{accession}_{assembly_name}"
    filename = f"{accession}_{assembly_name}_assembly_stats.txt"
    return f"{base}/{prefix}/{chunks[0]}/{chunks[1]}/{chunks[2]}/{dirname}/{filename}"


def fetch_and_parse_ncbi_assembly_stats_report(
    accession: str,
    assembly_name: str,
    outdir: Optional[str] = None,
) -> Dict[str, object]:
    """
    Download and parse NCBI legacy *_assembly_stats.txt.

    Robust parsing:
      - does not care about tabs vs spaces
      - matches 'all ... all ... all ... all <statistic> <value>'
    """
    url = ncbi_stats_report_url(accession, assembly_name)

    debug_dir = None
    if outdir:
        debug_dir = os.path.join(outdir, "ncbi_debug")
        os.makedirs(debug_dir, exist_ok=True)

    if debug_dir:
        with open(
            os.path.join(debug_dir, "ncbi_stats_report_url.txt"), "w", encoding="utf-8"
        ) as fh:
            fh.write(url + "\n")

    local_path = (
        os.path.join(debug_dir, f"{accession}_{assembly_name}_assembly_stats.txt")
        if debug_dir
        else os.path.join(
            tempfile.gettempdir(), f"{accession}_{assembly_name}_assembly_stats.txt"
        )
    )

    # Download (raise on failure)
    urllib.request.urlretrieve(url, local_path)

    out: Dict[str, object] = {
        "spanned_gaps": "",
        "top_level_count": "",
        "sex": "",
        "date": "",
    }

    # Regex for the “all ... statistic value” lines
    # Example: all  all  all  all  spanned-gaps  856
    stat_re = re.compile(r"^all\s+all\s+all\s+all\s+(?P<stat>\S+)\s+(?P<val>\S+)\s*$")

    with open(local_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")

            m = re.match(r"^\#\s*Sex\:\s+(.+)", line)
            if m:
                out["sex"] = m.group(1).strip()
                continue

            m = re.match(r"^\#\s*Date\:\s+([\d\-]+)", line)
            if m:
                out["date"] = m.group(1).strip()
                continue

            m = stat_re.match(line)
            if not m:
                continue

            stat = m.group("stat")
            val = m.group("val")

            # Take leading integer if present
            m2 = re.search(r"\d+", val)
            if not m2:
                continue
            ival = m2.group(0)

            if stat == "spanned-gaps":
                out["spanned_gaps"] = ival
            elif stat == "top-level-count":
                out["top_level_count"] = ival

    # Always write parsed output for inspection
    if debug_dir:
        with open(
            os.path.join(debug_dir, "parsed_ncbi_assembly_stats_report.json"),
            "w",
            encoding="utf-8",
        ) as fh:
            json.dump(
                {"url": url, "local_path": local_path, "parsed": out}, fh, indent=2
            )

    return out

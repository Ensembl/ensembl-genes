# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

#!/usr/bin/env python3
from __future__ import annotations

import urllib.request
import argparse
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any

# ----------------------------
# GFF3 parsing helpers
# ----------------------------


def open_maybe_gzip(path: str):
    """
    Open a text file that may be gzip-compressed.

    If the file path ends with '.gz', the file is opened using gzip in text
    mode. Otherwise, it is opened as a regular UTF-8 text file. In both cases,
    decoding errors are replaced to avoid parse failures on malformed input.

    Parameters
    ----------
    path : str
        Path to the input file, optionally ending with '.gz'.

    Returns
    -------
    TextIO
        A file-like object opened in text read mode.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_gff3_attrs(s: str) -> Dict[str, str]:
    """
    Parse the GFF3 attributes column into a dictionary.

    The GFF3 attributes field consists of semicolon-separated key-value pairs,
    where keys and values are separated by '='. Attributes without an '=' are
    treated as keys with empty string values.

    Parameters
    ----------
    s : str
        Raw attributes string from the 9th column of a GFF3 record.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping attribute keys to their corresponding values.
    """
    d: Dict[str, str] = {}
    for part in s.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = ""
    return d


def interval_len(start: int, end: int) -> int:
    """
    Compute the length of a genomic interval using GFF3 coordinates.

    GFF3 uses 1-based, inclusive coordinates, so the length is calculated as
    (end - start + 1).

    Parameters
    ----------
    start : int
        Start coordinate (1-based, inclusive).
    end : int
        End coordinate (1-based, inclusive).

    Returns
    -------
    int
        Length of the interval.
    """
    return end - start + 1


def compute_introns(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Compute intron intervals from a list of exon coordinates.

    Exons are assumed to be provided as (start, end) tuples using 1-based,
    inclusive coordinates. The function sorts exons by start position and
    returns the gaps between consecutive exons as introns. Overlapping or
    adjacent exons do not produce introns.

    Parameters
    ----------
    exons : List[Tuple[int, int]]
        List of exon intervals as (start, end) tuples.

    Returns
    -------
    List[Tuple[int, int]]
        List of intron intervals as (start, end) tuples. Returns an empty list
        if fewer than two exons are provided.
    """
    if len(exons) < 2:
        return []
    ex = sorted(exons, key=lambda x: x[0])
    introns: List[Tuple[int, int]] = []
    prev_s, prev_e = ex[0]
    for s, e in ex[1:]:
        a, b = prev_e + 1, s - 1
        if a <= b:
            introns.append((a, b))
        if e > prev_e:
            prev_e = e
    return introns


# ----------------------------
# Data model
# ----------------------------


@dataclass
class Transcript:
    """
    Representation of a transcript and its genomic structure.

    Coordinates follow GFF3 conventions (1-based, inclusive). Exons and CDS
    segments are stored as genomic intervals and are assumed to be on the same
    strand.

    Attributes
    ----------
    tid : str
        Transcript identifier.
    gid : str
        Parent gene identifier.
    strand : str
        Strand of the transcript ('+' or '-').
    exons : List[Tuple[int, int]]
        List of exon intervals as (start, end) tuples.
    cds : List[Tuple[int, int]]
        List of coding sequence (CDS) intervals as (start, end) tuples.
    is_canonical : bool
        Whether this transcript is designated as the canonical isoform.
    """

    tid: str
    gid: str
    strand: str
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds: List[Tuple[int, int]] = field(default_factory=list)
    is_canonical: bool = False

    def exon_spliced_len(self) -> int:
        """
        Return the total spliced length of the transcript exons.

        This is the sum of exon lengths using GFF3 (1-based, inclusive)
        coordinates.

        Returns
        -------
        int
            Total spliced exon length.
        """
        return sum(interval_len(s, e) for s, e in self.exons)

    def cds_spliced_len(self) -> int:
        """
        Return the total spliced length of the coding sequence (CDS).

        The length is computed as the sum of all CDS segment lengths, matching
        the common convention used for translated sequence length calculations.

        Returns
        -------
        int
            Total spliced CDS length.
        """
        # Perl API uses translateable_seq length, i.e. summed CDS segment lengths
        return sum(interval_len(s, e) for s, e in self.cds)

    def introns(self) -> List[Tuple[int, int]]:
        """
        Compute intron intervals from the transcript's exon structure.

        Introns are defined as the genomic gaps between consecutive exons after
        sorting exons by start coordinate.

        Returns
        -------
        List[Tuple[int, int]]
            List of intron intervals as (start, end) tuples.
        """
        return compute_introns(self.exons)

    def is_coding(self) -> bool:
        """
        Indicate whether the transcript contains a coding sequence.

        Returns
        -------
        bool
            True if one or more CDS intervals are present, False otherwise.
        """
        return len(self.cds) > 0


@dataclass
class Gene:
    """
    Representation of a gene and its associated transcripts.

    Gene coordinates follow GFF3 conventions (1-based, inclusive) and typically
    span the full genomic extent of all associated transcripts.

    Attributes
    ----------
    gid : str
        Gene identifier.
    start : int
        Genomic start coordinate of the gene.
    end : int
        Genomic end coordinate of the gene.
    strand : str
        Strand of the gene ('+' or '-').
    biotype : str, optional
        Gene biotype (e.g. 'protein_coding', 'lncRNA').
    transcripts : Dict[str, Transcript]
        Mapping of transcript IDs to Transcript objects.
    """

    gid: str
    start: int
    end: int
    strand: str
    biotype: str = ""
    transcripts: Dict[str, Transcript] = field(default_factory=dict)

    def length(self) -> int:
        """
        Return the genomic length of the gene.

        The length is computed using GFF3 (1-based, inclusive) coordinates.

        Returns
        -------
        int
            Length of the gene in base pairs.
        """
        return interval_len(self.start, self.end)


# ----------------------------
# Canonical transcript handling
# ----------------------------


def pick_canonical_transcript(g: Gene) -> Optional[Transcript]:
    """
    Select a canonical transcript for a gene.

    If one or more transcripts are explicitly flagged as canonical, the
    canonical transcript is chosen from that subset. If multiple flagged
    transcripts exist, the one with the greatest spliced CDS length is
    preferred, with spliced exon length used as a secondary tie-breaker.

    If no transcript is explicitly marked as canonical, a fallback selection
    is performed over all transcripts using the same criteria. This behavior
    mirrors typical downstream conventions when no canonical annotation
    (e.g. Ensembl_canonical) is present.

    Parameters
    ----------
    g : Gene
        Gene object containing one or more transcripts.

    Returns
    -------
    Optional[Transcript]
        The selected canonical Transcript, or None if the gene has no
        transcripts.
    """
    if not g.transcripts:
        return None
    flagged = [t for t in g.transcripts.values() if t.is_canonical]
    if flagged:
        # if multiple, choose most CDS then most exon length
        return sorted(
            flagged,
            key=lambda t: (t.cds_spliced_len(), t.exon_spliced_len()),
            reverse=True,
        )[0]
    txs = list(g.transcripts.values())
    return sorted(
        txs, key=lambda t: (t.cds_spliced_len(), t.exon_spliced_len()), reverse=True
    )[0]


# ----------------------------
# Biotype grouping (i.e. get_Biotype->biotype_group from Perl API)
# For Ensembl GFF3, biotype=protein_coding etc. is present.
# Noncoding subclasses are not always explicit at gene level.
# ----------------------------

_SMALL_NC = {"miRNA", "rRNA", "snRNA", "snoRNA", "tRNA", "scaRNA", "sRNA", "ncRNA"}
_LONG_NC = {"lncRNA", "lincRNA"}


def biotype_group(g: Gene) -> str:
    """
    Match Ensembl GFF3 conventions for this dataset:
    - protein_coding genes are coding (gene biotype=protein_coding)
    - ncRNA_gene genes are non-coding; subtype from biotype:
        * lncRNA -> lnoncoding
        * small structured RNAs -> snoncoding
        * everything else -> mnoncoding
    - pseudogene if biotype contains 'pseudogene' or gene feature type says so (handled upstream if you store it)
    """
    bt = (g.biotype or "").strip()

    # Pseudogene
    if "pseudogene" in bt.lower():
        return "pseudogene"

    # Coding
    if bt == "protein_coding":
        return "coding"

    # If it's not protein_coding but has CDS, treat as coding (rare but safe)
    if any(t.cds for t in g.transcripts.values()):
        return "coding"

    # Non-coding subtype from biotype
    if bt == "lncRNA":
        return "lnoncoding"

    # Treat these as "small ncRNAs"
    small = {"snoRNA", "snRNA", "miRNA", "rRNA", "tRNA", "scaRNA", "sRNA"}
    if bt in small:
        return "snoncoding"

    # Everything else noncoding goes to misc
    return "mnoncoding"


# ----------------------------
# Parse Ensembl-style GFF3
# ----------------------------


def parse_ensembl_gff3(gff_path: str) -> Dict[str, Gene]:
    """
    Parse an Ensembl-style GFF3 file into Gene and Transcript objects.

    This parser is tailored to Ensembl GFF3 conventions and extracts genes,
    transcripts, exons, and CDS features into an in-memory data model. All
    coordinates are interpreted using GFF3 semantics (1-based, inclusive).

    Supported gene feature types include protein-coding and non-coding genes
    (e.g. gene, ncRNA_gene, pseudogene). Supported transcript feature types
    include mRNA, transcript, and common Ensembl RNA biotypes (lncRNA, miRNA,
    snoRNA, etc.).

    Canonical transcripts are detected using the Ensembl-specific
    ``tag=Ensembl_canonical`` attribute when present.

    Notes
    -----
    - Transcript-to-gene relationships are resolved via the ``Parent``
      attribute. Ensembl-style prefixes (e.g. ``gene:`` or ``transcript:``)
      are preserved as-is.
    - If a transcript feature appears before its corresponding gene feature,
      a temporary gene record is created and later reused.
    - Exons and CDS features are attached to transcripts using the same
      ``Parent`` resolution logic.
    - Features with missing or malformed identifiers are skipped silently.

    Parameters
    ----------
    gff_path : str
        Path to the Ensembl GFF3 file. May be plain text or gzip-compressed
        (``.gz``).

    Returns
    -------
    Dict[str, Gene]
        Dictionary mapping gene IDs to populated Gene objects, each containing
        its associated Transcript objects with exon and CDS coordinates.
    """
    genes: Dict[str, Gene] = {}
    transcript_to_gene: Dict[str, str] = {}

    with open_maybe_gzip(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start_s, end_s, score, strand, phase, attrs_s = parts
            start, end = int(start_s), int(end_s)
            attrs = parse_gff3_attrs(attrs_s)

            if ftype in {
                "gene",
                "ncRNA_gene",
                "pseudogene",
                "transposable_element_gene",
            }:
                gid = attrs.get("ID", "")
                if not gid:
                    continue
                bt = (
                    attrs.get("biotype", "")
                    or attrs.get("gene_biotype", "")
                    or attrs.get("gene_type", "")
                )
                genes[gid] = Gene(
                    gid=gid, start=start, end=end, strand=strand, biotype=bt
                )

            elif ftype in {
                "mRNA",
                "transcript",
                "snoRNA",
                "snRNA",
                "miRNA",
                "rRNA",
                "tRNA",
                "lnc_RNA",
                "lncRNA",
                "ncRNA",
                "antisense_RNA",
                "sRNA",
                "scaRNA",
            }:
                tid = attrs.get("ID", "")
                parent = attrs.get("Parent", "")
                if not tid or not parent:
                    continue
                # Ensembl uses Parent=gene:XXXX, IDs include prefixes like gene:/transcript:
                gid = parent.split(",")[0]
                transcript_to_gene[tid] = gid

                g = genes.get(gid)
                if g is None:
                    # gene may appear later, create stub with transcript span for now
                    g = Gene(gid=gid, start=start, end=end, strand=strand, biotype="")
                    genes[gid] = g
                else:
                    # keep real gene coords from gene feature
                    pass

                t = g.transcripts.get(tid)
                if t is None:
                    t = Transcript(tid=tid, gid=gid, strand=strand)
                    g.transcripts[tid] = t

                # canonical from tag=Ensembl_canonical
                tag = attrs.get("tag", "")
                if isinstance(tag, str) and "Ensembl_canonical" in tag:
                    t.is_canonical = True

                # transcript biotype can exist, but Perl API groups at gene level
                # t_biotype = attrs.get("biotype","")

            elif ftype == "exon":
                parent = attrs.get("Parent", "")
                if not parent:
                    continue
                tid = parent.split(",")[0]
                gid_opt = transcript_to_gene.get(tid)
                if gid_opt is None:
                    continue
                gid = gid_opt

                g = genes.get(gid)
                if not g:
                    continue

                t = g.transcripts.get(tid)
                if not t:
                    t = Transcript(tid=tid, gid=gid, strand=strand)
                    g.transcripts[tid] = t
                t.exons.append((start, end))

            elif ftype == "CDS":
                parent = attrs.get("Parent", "")
                if not parent:
                    continue
                tid = parent.split(",")[0]
                gid_opt = transcript_to_gene.get(tid)
                if gid_opt is None:
                    continue
                gid = gid_opt

                g = genes.get(gid)
                if not g:
                    continue

                t = g.transcripts.get(tid)
                if not t:
                    t = Transcript(tid=tid, gid=gid, strand=strand)
                    g.transcripts[tid] = t
                t.cds.append((start, end))

    return genes


# ----------------------------
# Metric computation
# ----------------------------


def process_coding_genes(
    genes: List[Gene], scientific_name: str
) -> Tuple[Dict[str, object], int]:
    """
    Compute summary statistics for a set of coding genes.

    This function aggregates gene-, transcript-, exon-, CDS-, and intron-level
    metrics across the provided genes, following conventions commonly used in
    Ensembl/Perl-based genome statistics pipelines. All coordinate calculations
    use GFF3 semantics (1-based, inclusive).

    Metrics are computed across all transcripts belonging to the input genes,
    while certain values (e.g. sequence length, coding sequence length) are
    derived specifically from each gene's canonical transcript when available.

    Parameters
    ----------
    genes : List[Gene]
        List of Gene objects to process. These are assumed to represent
        protein-coding genes and to contain populated Transcript, exon, and CDS
        data.
    scientific_name : str
        Scientific name of the organism, included verbatim in the output
        summary.

    Returns
    -------
    Tuple[Dict[str, object], int]
        A tuple containing:

        1. A dictionary of computed summary metrics, including:
           - Gene counts and genomic spans
           - Transcript and coding transcript counts
           - Exon, CDS (coding exon), and intron counts and average lengths
           - Per-gene and per-transcript averages
           - Canonical transcript–based sequence statistics

           Missing or undefined values are represented as empty strings to
           match downstream reporting expectations.

        2. An integer representing the total spliced CDS length summed across
           canonical coding transcripts only. This mirrors the Perl
           ``translateable_seq`` length aggregation used in some pipelines.

    Notes
    -----
    - Canonical transcripts are selected using ``pick_canonical_transcript``.
    - All exons and introns are counted across *all* transcripts of coding
      genes, not just canonical transcripts.
    - "Coding exons" correspond to individual CDS segments, not exon features.
    - Genes with no transcripts or transcripts with no CDS contribute only to
      applicable metrics.
    - If the input gene list is empty, all numeric metrics are returned as empty
      strings or zero as appropriate.
    """
    if not genes:
        return (
            {
                "Scientific name": scientific_name,
                "Coding genes": 0,
                "Average genomic span": "",
                "Average sequence length": "",
                "Average CDS length": "",
                "Shortest gene": "",
                "Longest gene": "",
                "Total transcripts": "",
                "Coding transcripts": "",
                "Transcripts per gene": "",
                "Coding transcripts per gene": "",
                "Total exons": "",
                "Total coding exons": "",
                "Average exon length": "",
                "Average coding exon length": "",
                "Average exons per transcript": "",
                "Average coding exons per coding transcript": "",
                "Total introns": "",
                "Average intron length": "",
            },
            0,
        )

    gene_count = len(genes)
    shortest_gene = 10**18
    longest_gene = 0

    total_gene_len = 0
    total_canonical_tx_len = 0

    total_transcripts = 0
    coding_transcripts = 0

    total_exons = 0
    total_exon_len = 0

    total_cds_segments = 0
    total_cds_len = 0

    total_introns = 0
    total_intron_len = 0

    # From Perl API: sum canonical coding transcript translateable_seq length
    total_coding_sequence_length = 0

    for g in genes:
        gl = g.length()
        shortest_gene = min(shortest_gene, gl)
        longest_gene = max(longest_gene, gl)
        total_gene_len += gl

        canonical = pick_canonical_transcript(g)
        if canonical:
            total_canonical_tx_len += canonical.exon_spliced_len()

        for t in g.transcripts.values():
            total_transcripts += 1

            # From Perl API: count exons and introns across ALL transcripts for coding genes
            total_exons += len(t.exons)
            total_exon_len += sum(interval_len(s, e) for s, e in t.exons)

            intrs = t.introns()
            total_introns += len(intrs)
            total_intron_len += sum(interval_len(s, e) for s, e in intrs)

            # coding transcript defined by having CDS (translation)
            if t.is_coding():
                coding_transcripts += 1
                cds_len = t.cds_spliced_len()
                total_cds_len += cds_len
                total_cds_segments += len(t.cds)

                if canonical and t.tid == canonical.tid:
                    total_coding_sequence_length += cds_len

    out = {
        "Scientific name": scientific_name,
        "Coding genes": gene_count,
        "Average genomic span": round(total_gene_len / gene_count, 2)
        if gene_count
        else "",
        "Average sequence length": round(total_canonical_tx_len / gene_count, 2)
        if gene_count
        else "",
        "Average CDS length": round(total_cds_len / coding_transcripts, 2)
        if coding_transcripts
        else "",
        "Shortest gene": shortest_gene if shortest_gene < 10**18 else "",
        "Longest gene": longest_gene,
        "Total transcripts": total_transcripts,
        "Coding transcripts": coding_transcripts,
        "Transcripts per gene": round(total_transcripts / gene_count, 2)
        if gene_count
        else "",
        "Coding transcripts per gene": round(coding_transcripts / gene_count, 2)
        if gene_count
        else "",
        "Total exons": total_exons,
        # From Perl API: "coding exons" counts CDS segments (get_all_CDS)
        "Total coding exons": total_cds_segments,
        "Average exon length": round(total_exon_len / total_exons, 2)
        if total_exons
        else "",
        "Average coding exon length": round(total_cds_len / total_cds_segments, 2)
        if total_cds_segments
        else "",
        "Average exons per transcript": round(total_exons / total_transcripts, 2)
        if total_transcripts
        else "",
        "Average coding exons per coding transcript": round(
            total_cds_segments / coding_transcripts, 2
        )
        if coding_transcripts
        else "",
        "Total introns": total_introns,
        "Average intron length": (
            round(total_intron_len / total_introns, 2) if total_introns else ""
        ),
    }
    if total_introns == 0:
        out["Average intron length"] = ""
    return out, total_coding_sequence_length


def process_non_coding_genes(
    genes: List[Gene], scientific_name: str
) -> Dict[str, object]:
    """
    Compute summary statistics for a set of non-coding genes.

    This function aggregates gene-, transcript-, exon-, and intron-level metrics
    for non-coding genes, following conventions commonly used in genome
    annotation summary pipelines. All coordinate calculations use GFF3
    semantics (1-based, inclusive).

    Non-coding genes are additionally categorized into small non-coding,
    long non-coding, or miscellaneous non-coding groups based on their biotype,
    as determined by ``biotype_group``.

    Parameters
    ----------
    genes : List[Gene]
        List of non-coding Gene objects to process. These are assumed to contain
        populated Transcript and exon data.
    scientific_name : str
        Scientific name of the organism, included verbatim in the output
        summary.

    Returns
    -------
    Dict[str, object]
        Dictionary of computed summary metrics, including:
        - Counts of non-coding genes by category (small, long, miscellaneous)
        - Genomic span statistics
        - Transcript counts and per-gene averages
        - Exon and intron counts and average lengths
        - Canonical transcript–based sequence length statistics

        Missing or undefined values are represented as empty strings to match
        downstream reporting expectations.

    Notes
    -----
    - Canonical transcripts are selected using ``pick_canonical_transcript``.
    - Exons and introns are counted across *all* transcripts for each gene,
      not just canonical transcripts.
    - Coding sequence (CDS) features are ignored by design.
    - If the input gene list is empty, all numeric metrics are returned as empty
      strings or zero as appropriate.
    """
    if not genes:
        return {
            "Scientific name": scientific_name,
            "Non-coding genes": 0,
            "Small non-coding genes": "",
            "Long non-coding genes": "",
            "Misc non-coding genes": "",
            "Average genomic span": "",
            "Average sequence length": "",
            "Shortest gene": "",
            "Longest gene": "",
            "Total transcripts": "",
            "Transcripts per gene": "",
            "Total exons": "",
            "Average exon length": "",
            "Average exons per transcript": "",
            "Total introns": "",
            "Average intron length": "",
        }

    gene_count = len(genes)
    shortest_gene = 10**18
    longest_gene = 0

    total_gene_len = 0
    total_canonical_tx_len = 0

    sn = ln = mn = 0
    total_transcripts = 0
    total_exons = 0
    total_exon_len = 0
    total_introns = 0
    total_intron_len = 0

    for g in genes:
        gl = g.length()
        shortest_gene = min(shortest_gene, gl)
        longest_gene = max(longest_gene, gl)
        total_gene_len += gl

        grp = biotype_group(g)
        if grp == "snoncoding":
            sn += 1
        elif grp == "lnoncoding":
            ln += 1
        else:
            mn += 1

        canonical = pick_canonical_transcript(g)
        if canonical:
            total_canonical_tx_len += canonical.exon_spliced_len()

        for t in g.transcripts.values():
            total_transcripts += 1
            total_exons += len(t.exons)
            total_exon_len += sum(interval_len(s, e) for s, e in t.exons)

            intrs = t.introns()
            total_introns += len(intrs)
            total_intron_len += sum(interval_len(s, e) for s, e in intrs)

    out = {
        "Scientific name": scientific_name,
        "Non-coding genes": gene_count,
        "Small non-coding genes": sn,
        "Long non-coding genes": ln,
        "Misc non-coding genes": mn,
        "Average genomic span": round(total_gene_len / gene_count, 2)
        if gene_count
        else "",
        "Average sequence length": round(total_canonical_tx_len / gene_count, 2)
        if gene_count
        else "",
        "Shortest gene": shortest_gene if shortest_gene < 10**18 else "",
        "Longest gene": longest_gene,
        "Total transcripts": total_transcripts,
        "Transcripts per gene": round(total_transcripts / gene_count, 2)
        if gene_count
        else "",
        "Total exons": total_exons,
        "Average exon length": round(total_exon_len / total_exons, 2)
        if total_exons
        else "",
        "Average exons per transcript": round(total_exons / total_transcripts, 2)
        if total_transcripts
        else "",
        "Total introns": total_introns,
        "Average intron length": (
            round(total_intron_len / total_introns, 2) if total_introns else ""
        ),
    }
    if total_introns == 0:
        out["Average intron length"] = ""
    return out


def process_pseudogenes(genes: List[Gene], scientific_name: str) -> Dict[str, object]:
    """
    Compute summary statistics for a set of pseudogenes.

    This function aggregates gene-, transcript-, exon-, and intron-level metrics
    for pseudogenes using conventions commonly applied in genome annotation
    summary pipelines. All coordinate calculations use GFF3 semantics
    (1-based, inclusive).

    Pseudogenes are treated as non-coding for the purposes of these statistics:
    coding sequence (CDS) features are ignored, and no biotype sub-classification
    is performed.

    Parameters
    ----------
    genes : List[Gene]
        List of pseudogene Gene objects to process. These are assumed to contain
        populated Transcript and exon data.
    scientific_name : str
        Scientific name of the organism, included verbatim in the output
        summary.

    Returns
    -------
    Dict[str, object]
        Dictionary of computed summary metrics, including:
        - Pseudogene counts
        - Genomic span statistics
        - Transcript counts and per-gene averages
        - Exon and intron counts and average lengths
        - Canonical transcript–based sequence length statistics

        Missing or undefined values are represented as empty strings to match
        downstream reporting expectations.

    Notes
    -----
    - Canonical transcripts are selected using ``pick_canonical_transcript``.
    - Exons and introns are counted across *all* transcripts, not just canonical
      transcripts.
    - If the input gene list is empty, all numeric metrics are returned as empty
      strings or zero as appropriate.
    """
    if not genes:
        return {
            "Scientific name": scientific_name,
            "Pseudogenes": 0,
            "Average genomic span": "",
            "Average sequence length": "",
            "Shortest gene": "",
            "Longest gene": "",
            "Total transcripts": "",
            "Transcripts per gene": "",
            "Total exons": "",
            "Average exon length": "",
            "Average exons per transcript": "",
            "Total introns": "",
            "Average intron length": "",
        }

    gene_count = len(genes)
    shortest_gene = 10**18
    longest_gene = 0

    total_gene_len = 0
    total_canonical_tx_len = 0
    total_transcripts = 0
    total_exons = 0
    total_exon_len = 0
    total_introns = 0
    total_intron_len = 0

    for g in genes:
        gl = g.length()
        shortest_gene = min(shortest_gene, gl)
        longest_gene = max(longest_gene, gl)
        total_gene_len += gl

        canonical = pick_canonical_transcript(g)
        if canonical:
            total_canonical_tx_len += canonical.exon_spliced_len()

        for t in g.transcripts.values():
            total_transcripts += 1
            total_exons += len(t.exons)
            total_exon_len += sum(interval_len(s, e) for s, e in t.exons)

            intrs = t.introns()
            total_introns += len(intrs)
            total_intron_len += sum(interval_len(s, e) for s, e in intrs)

    out = {
        "Scientific name": scientific_name,
        "Pseudogenes": gene_count,
        "Average genomic span": round(total_gene_len / gene_count, 2)
        if gene_count
        else "",
        "Average sequence length": round(total_canonical_tx_len / gene_count, 2)
        if gene_count
        else "",
        "Shortest gene": shortest_gene if shortest_gene < 10**18 else "",
        "Longest gene": longest_gene,
        "Total transcripts": total_transcripts,
        "Transcripts per gene": round(total_transcripts / gene_count, 2)
        if gene_count
        else "",
        "Total exons": total_exons,
        "Average exon length": round(total_exon_len / total_exons, 2)
        if total_exons
        else "",
        "Average exons per transcript": round(total_exons / total_transcripts, 2)
        if total_transcripts
        else "",
        "Total introns": total_introns,
        "Average intron length": (
            round(total_intron_len / total_introns, 2) if total_introns else ""
        ),
    }
    if total_introns == 0:
        out["Average intron length"] = ""
    return out


# ----------------------------
# Assembly stats (using NCBI Datasets)
# ----------------------------


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


# ----------------------------
# Fallback to legacy stats file
# if missing stats form Datasets
# ----------------------------


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


# ----------------------------
# Output
# ----------------------------


def write_tsv(path: str, headers: List[str], row: Dict[str, object]) -> None:
    """
    Write a single-row TSV file with a specified header order.

    This function writes a tab-separated values (TSV) file containing a header
    row followed by exactly one data row. Output values are written in the
    order specified by ``headers``. Missing values (None) are written as empty
    fields.

    Parameters
    ----------
    path : str
        Path to the output TSV file.
    headers : List[str]
        Ordered list of column names to write as the header row.
    row : Dict[str, object]
        Mapping from column names to values. Values are converted to strings;
        keys missing from the mapping or explicitly set to None result in empty
        fields in the output.

    Returns
    -------
    None
    """
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        fh.write(
            "\t".join("" if row.get(h) is None else str(row.get(h)) for h in headers)
            + "\n"
        )


def write_json_from_tsv(tsv_path: str, json_path: str) -> None:
    """
    Convert a TSV file into a JSON array of row objects.

    The input TSV file is expected to contain a single header row followed by
    one or more data rows. Each data row is converted into a dictionary keyed
    by the header fields, and all rows are written as a JSON array.

    Empty fields in the TSV are preserved as empty strings in the resulting
    JSON, matching TSV semantics and avoiding implicit type coercion.

    Parameters
    ----------
    tsv_path : str
        Path to the input TSV file.
    json_path : str
        Path to the output JSON file.

    Returns
    -------
    None
    """

    rows: List[Dict[str, Any]] = []

    with open(tsv_path, "rt", encoding="utf-8") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if header is None:
                header = fields
                continue

            row = {}
            for k, v in zip(header, fields):
                # Keep empty strings as empty strings (matches TSV semantics)
                row[k] = v
            rows.append(row)

    with open(json_path, "wt", encoding="utf-8") as out:
        json.dump(rows, out, indent=2, sort_keys=False)


def write_meta_sql(
    out_path: str,
    dbname: str,
    species_id: str,
    coding_stats: Dict[str, object],
    noncoding_stats: Dict[str, object],
    pseudogene_stats: Dict[str, object],
    assembly_row: Dict[str, object],
) -> None:
    """
    Write SQL statements to populate the Ensembl-style `meta` table with
    gene build and assembly statistics.

    This function generates a SQL script containing `INSERT IGNORE` statements
    for numeric metadata values derived from coding, non-coding, pseudogene,
    and assembly statistics. The output is intended to be executed against a
    database that follows the Ensembl schema conventions.

    Only values that appear to be numeric (integers or floats, including
    negatives) are written. Empty values, missing keys, and non-numeric values
    are silently skipped.

    The resulting file:
      * Optionally begins with a `USE <dbname>;` statement
      * Contains one `INSERT IGNORE INTO meta (...)` statement per valid value
      * Never overwrites existing rows in the `meta` table

    Parameters
    ----------
    out_path : str
        Path to the output `.sql` file to be written. Any existing file at this
        path will be overwritten.

    dbname : str
        Name of the target database. If provided (non-empty), a `USE <dbname>;`
        statement is written at the top of the file. If empty, no `USE`
        statement is emitted.

    species_id : str
        Numeric species identifier used in the `meta.species_id` column.
        This value is written verbatim into the SQL statements and must be
        non-empty. A `ValueError` is raised if it is missing.

    coding_stats : Dict[str, object]
        Mapping of human-readable column names to values for coding gene
        statistics (e.g. "Coding genes", "Average exon length"). Only the keys
        referenced internally by this function are consulted.

    noncoding_stats : Dict[str, object]
        Mapping of human-readable column names to values for non-coding gene
        statistics.

    pseudogene_stats : Dict[str, object]
        Mapping of human-readable column names to values for pseudogene
        statistics.

    assembly_row : Dict[str, object]
        Mapping of human-readable column names to values for genome assembly
        statistics (e.g. N50, total genome length, GC percentage).

    Returns
    -------
    None
        This function does not return a value. Its sole side effect is writing
        a SQL script to `out_path`.

    Raises
    ------
    ValueError
        If `species_id` is empty or evaluates to False.

    Notes
    -----
    * Values are normalized by converting them to strings and stripping
      surrounding whitespace.
    * A value is considered valid only if it matches the regular expression
      for a signed integer or floating-point number.
    * SQL values are emitted without quotes, assuming numeric semantics.
    * The function intentionally uses `INSERT IGNORE` to avoid conflicts with
      existing meta keys already present in the database.

    Example
    -------
    >>> write_meta_sql(
    ...     out_path="meta.sql",
    ...     dbname="homo_sapiens_core",
    ...     species_id="1",
    ...     coding_stats={"Coding genes": 20000},
    ...     noncoding_stats={},
    ...     pseudogene_stats={},
    ...     assembly_row={"Contig N50": 50000000},
    ... )
    """

    if not species_id:
        raise ValueError("species_id is required to write meta SQL")

    coding_map = {
        "genebuild.stats.coding_genes": "Coding genes",
        "genebuild.stats.average_genomic_span": "Average genomic span",
        "genebuild.stats.average_sequence_length": "Average sequence length",
        "genebuild.stats.average_cds_length": "Average CDS length",
        "genebuild.stats.shortest_gene_length": "Shortest gene",
        "genebuild.stats.longest_gene_length": "Longest gene",
        "genebuild.stats.total_transcripts": "Total transcripts",
        "genebuild.stats.coding_transcripts": "Coding transcripts",
        "genebuild.stats.transcripts_per_gene": "Transcripts per gene",
        "genebuild.stats.coding_transcripts_per_gene": "Coding transcripts per gene",
        "genebuild.stats.total_exons": "Total exons",
        "genebuild.stats.total_coding_exons": "Total coding exons",
        "genebuild.stats.average_exon_length": "Average exon length",
        "genebuild.stats.average_coding_exon_length": "Average coding exon length",
        "genebuild.stats.average_coding_exons_per_transcript": "Average exons per transcript",
        "genebuild.stats.average_coding_exons_per_coding_transcript": "Average coding exons per coding transcript",
        "genebuild.stats.total_introns": "Total introns",
        "genebuild.stats.average_intron_length": "Average intron length",
    }

    noncoding_map = {
        "genebuild.stats.nc_non_coding_genes": "Non-coding genes",
        "genebuild.stats.nc_small_non_coding_genes": "Small non-coding genes",
        "genebuild.stats.nc_long_non_coding_genes": "Long non-coding genes",
        "genebuild.stats.nc_misc_non_coding_genes": "Misc non-coding genes",
        "genebuild.stats.nc_average_genomic_span": "Average genomic span",
        "genebuild.stats.nc_average_sequence_length": "Average sequence length",
        "genebuild.stats.nc_shortest_gene_length": "Shortest gene",
        "genebuild.stats.nc_longest_gene_length": "Longest gene",
        "genebuild.stats.nc_total_transcripts": "Total transcripts",
        "genebuild.stats.nc_transcripts_per_gene": "Transcripts per gene",
        "genebuild.stats.nc_total_exons": "Total exons",
        "genebuild.stats.nc_average_exon_length": "Average exon length",
        "genebuild.stats.nc_average_exons_per_transcript": "Average exons per transcript",
        "genebuild.stats.nc_total_introns": "Total introns",
        "genebuild.stats.nc_average_intron_length": "Average intron length",
    }

    pseudogene_map = {
        "genebuild.stats.ps_pseudogenes": "Pseudogenes",
        "genebuild.stats.ps_average_genomic_span": "Average genomic span",
        "genebuild.stats.ps_average_sequence_length": "Average sequence length",
        "genebuild.stats.ps_shortest_gene_length": "Shortest gene",
        "genebuild.stats.ps_longest_gene_length": "Longest gene",
        "genebuild.stats.ps_total_transcripts": "Total transcripts",
        "genebuild.stats.ps_transcripts_per_gene": "Transcripts per gene",
        "genebuild.stats.ps_total_exons": "Total exons",
        "genebuild.stats.ps_average_exon_length": "Average exon length",
        "genebuild.stats.ps_average_exons_per_transcript": "Average exons per transcript",
        "genebuild.stats.ps_total_introns": "Total introns",
        "genebuild.stats.ps_average_intron_length": "Average intron length",
    }

    assembly_map = {
        "assembly.stats.contig_n50": "Contig N50",
        "assembly.stats.total_genome_length": "Total genome length",
        "assembly.stats.total_coding_sequence_length": "Total coding sequence length",
        "assembly.stats.total_gap_length": "Total gap length",
        "assembly.stats.spanned_gaps": "Spanned gaps",
        "assembly.stats.chromosomes": "Chromosomes",
        "assembly.stats.toplevel_sequences": "Toplevel sequences",
        "assembly.stats.component_sequences": "Component sequences",
        "assembly.stats.gc_percentage": "% GC",
    }

    def normalize_value(v) -> str:
        if v is None:
            return ""
        s = str(v).strip()
        return s

    def is_empty(v: str) -> bool:
        return v == ""

    # Keep values if they look numeric, else we skip
    num_re = re.compile(r"^-?\d+(\.\d+)?$")

    def emit(fh, meta_key: str, raw_val: str) -> None:
        if is_empty(raw_val):
            return
        if not num_re.match(raw_val):
            return
        fh.write(
            "INSERT IGNORE INTO meta (species_id, meta_key, meta_value) "
            f"VALUES({species_id}, '{meta_key}', {raw_val});\n"
        )

    with open(out_path, "w", encoding="utf-8") as fh:
        if dbname:
            fh.write(f"USE {dbname};\n")

        for mk, col in coding_map.items():
            emit(fh, mk, normalize_value(coding_stats.get(col, "")))

        for mk, col in noncoding_map.items():
            emit(fh, mk, normalize_value(noncoding_stats.get(col, "")))

        for mk, col in pseudogene_map.items():
            emit(fh, mk, normalize_value(pseudogene_stats.get(col, "")))

        for mk, col in assembly_map.items():
            emit(fh, mk, normalize_value(assembly_row.get(col, "")))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Compute Ensembl-like metrics from Ensembl GFF3 + NCBI datasets assembly stats."
    )
    ap.add_argument(
        "--gff3", required=True, help="Path to Ensembl GFF3 (.gff3 or .gff3.gz)"
    )
    ap.add_argument("--outdir", required=True, help="Output directory for TSVs")
    ap.add_argument(
        "--assembly_accession",
        required=True,
        help="GCA_/GCF_ accession for NCBI Datasets stats",
    )
    ap.add_argument(
        "--scientific_name",
        default="",
        help="Override scientific name (else from NCBI report if available)",
    )
    ap.add_argument(
        "--taxon_id",
        default="",
        help="Override taxonomy id (else from NCBI report if available)",
    )
    ap.add_argument("--strain", default="", help="Strain/Breed/Cultivar (optional)")
    ap.add_argument(
        "--sex",
        default="",
        help="Sex (optional override; otherwise not always present in NCBI report)",
    )
    ap.add_argument(
        "--keep_ncbi",
        action="store_true",
        help="Keep NCBI datasets zip and extracted reports in outdir/ncbi_debug",
    )
    ap.add_argument(
        "--ncbi_debug_dir",
        default="ncbi_debug",
        help="Subdirectory under outdir to store NCBI datasets artifacts",
    )
    ap.add_argument(
        "--json",
        action="store_true",
        help="Also write JSON versions of all TSV output files",
    )
    ap.add_argument(
        "--sql",
        action="store_true",
        help="Also write a SQL file with INSERT IGNORE statements into meta table",
    )
    ap.add_argument("--dbname", default="", help="Database name, required for SQL file")
    ap.add_argument(
        "--species_id",
        type=int,
        default=1,
        help="Species ID for meta inserts (default: 1; required if --sql is set)",
    )
    ap.add_argument(
        "--sql_filename",
        default="",
        help="Optional explicit SQL output path (default: outdir/stats_<assembly_accession>.sql)",
    )

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    genes = parse_ensembl_gff3(args.gff3)

    coding: List[Gene] = []
    noncoding: List[Gene] = []
    pseudo: List[Gene] = []

    for g in genes.values():
        grp = biotype_group(g)
        if grp == "coding":
            coding.append(g)
        elif grp == "pseudogene":
            pseudo.append(g)
        else:
            noncoding.append(g)

    asm = ncbi_assembly_stats(args.assembly_accession, outdir=args.outdir)

    scientific_name: str = str(args.scientific_name or asm.get("scientific_name") or "")
    taxon_id: str = str(args.taxon_id or asm.get("taxon_id") or "")

    coding_stats, total_coding_seq_len = process_coding_genes(coding, scientific_name)
    noncoding_stats = process_non_coding_genes(noncoding, scientific_name)
    pseudogene_stats = process_pseudogenes(pseudo, scientific_name)

    # Assembly headers
    assembly_headers = [
        "Scientific name",
        "Sex",
        "Breed/Cultivar/Strain",
        "Taxonomy id",
        "Assembly name",
        "Assembly accession",
        "Assembly date",
        "Contig N50",
        "Total genome length",
        "Total coding sequence length",
        "Total gap length",
        "Spanned gaps",
        "Chromosomes",
        "Toplevel sequences",
        "Component sequences",
        "% GC",
    ]

    # Note: some metrics are not reliably present in NCBI assembly report, keep these blank (should add a compute from FASTA option for GC%)
    assembly_row = {
        "Scientific name": scientific_name,
        "Sex": args.sex or "",
        "Breed/Cultivar/Strain": args.strain or "",
        "Taxonomy id": taxon_id,
        "Assembly name": asm.get("assembly_name") or "",
        "Assembly accession": asm.get("assembly_accession") or args.assembly_accession,
        "Assembly date": asm.get("assembly_date") or "",
        "Contig N50": asm.get("contig_n50") or "",
        "Total genome length": asm.get("total_length") or "",
        "Total coding sequence length": total_coding_seq_len,
        "Total gap length": asm.get("total_gap_length") or "",
        "Spanned gaps": asm.get("spanned_gaps") or "",
        "Chromosomes": asm.get("molecule_count") or "",
        "Toplevel sequences": asm.get("toplevel_sequences") or "",
        "Component sequences": asm.get("component_count") or "",
        "% GC": asm.get("gc_percent") or "",
    }

    coding_headers = [
        "Scientific name",
        "Coding genes",
        "Average genomic span",
        "Average sequence length",
        "Average CDS length",
        "Shortest gene",
        "Longest gene",
        "Total transcripts",
        "Coding transcripts",
        "Transcripts per gene",
        "Coding transcripts per gene",
        "Total exons",
        "Total coding exons",
        "Average exon length",
        "Average coding exon length",
        "Average exons per transcript",
        "Average coding exons per coding transcript",
        "Total introns",
        "Average intron length",
    ]
    noncoding_headers = [
        "Scientific name",
        "Non-coding genes",
        "Small non-coding genes",
        "Long non-coding genes",
        "Misc non-coding genes",
        "Average genomic span",
        "Average sequence length",
        "Shortest gene",
        "Longest gene",
        "Total transcripts",
        "Transcripts per gene",
        "Total exons",
        "Average exon length",
        "Average exons per transcript",
        "Total introns",
        "Average intron length",
    ]
    pseudogene_headers = [
        "Scientific name",
        "Pseudogenes",
        "Average genomic span",
        "Average sequence length",
        "Shortest gene",
        "Longest gene",
        "Total transcripts",
        "Transcripts per gene",
        "Total exons",
        "Average exon length",
        "Average exons per transcript",
        "Total introns",
        "Average intron length",
    ]

    write_tsv(
        os.path.join(args.outdir, "coding_stats.tsv"), coding_headers, coding_stats
    )
    write_tsv(
        os.path.join(args.outdir, "noncoding_stats.tsv"),
        noncoding_headers,
        noncoding_stats,
    )
    write_tsv(
        os.path.join(args.outdir, "pseudogene_stats.tsv"),
        pseudogene_headers,
        pseudogene_stats,
    )
    write_tsv(
        os.path.join(args.outdir, "assembly_stats.tsv"), assembly_headers, assembly_row
    )

    if args.json:
        write_json_from_tsv(
            os.path.join(args.outdir, "coding_stats.tsv"),
            os.path.join(args.outdir, "coding_stats.json"),
        )

        write_json_from_tsv(
            os.path.join(args.outdir, "noncoding_stats.tsv"),
            os.path.join(args.outdir, "noncoding_stats.json"),
        )

        write_json_from_tsv(
            os.path.join(args.outdir, "pseudogene_stats.tsv"),
            os.path.join(args.outdir, "pseudogene_stats.json"),
        )

        write_json_from_tsv(
            os.path.join(args.outdir, "assembly_stats.tsv"),
            os.path.join(args.outdir, "assembly_stats.json"),
        )

    if args.sql:
        if not args.species_id:
            raise SystemExit("--species_id is required when using --sql")

        sql_path = args.sql_filename or os.path.join(
            args.outdir, f"stats_{args.assembly_accession}.sql"
        )

        write_meta_sql(
            out_path=sql_path,
            dbname=args.dbname,
            species_id=str(args.species_id),
            coding_stats=coding_stats,
            noncoding_stats=noncoding_stats,
            pseudogene_stats=pseudogene_stats,
            assembly_row=assembly_row,
        )

    print("[done]", file=sys.stderr)
    print(f"coding genes: {len(coding)}", file=sys.stderr)
    print(f"noncoding genes: {len(noncoding)}", file=sys.stderr)
    print(f"pseudogenes: {len(pseudo)}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

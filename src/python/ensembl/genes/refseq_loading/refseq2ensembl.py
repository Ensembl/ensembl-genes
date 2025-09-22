import re
import os
import csv
import requests
import threading
import gzip
import mysql.connector
from typing import Optional, Union, Dict, List

# Define the major RefSeq clade groups (taxonomic divisions)
NCBI_GROUPS = [
    "vertebrate_mammalian", "vertebrate_other", "plant", "invertebrate", 
    "fungi", "protozoa"#,"archaea", "bacteria", "viral"
]

def list_available_annotations(
    base_dir: str = "refseq_data",
    group_filter: Optional[str] = None,
    max_print: Optional[int] = 50,
    return_dict: bool = False,
    output_tsv: Optional[str] = "refseq_annotation_metadata.tsv"
) -> Union[None, Dict[str, List[str]]]:
    """
    List RefSeq genome annotations grouped by clade and export metadata to TSV.

    Parameters
    ----------
    base_dir : str
        Local base directory for RefSeq downloads.
    group_filter : str or None
        Restrict to a single clade (e.g. 'vertebrate_mammalian') if set.
    max_print : int or None
        Max species to print per clade. Use None to print all.
    return_dict : bool
        Return dict of clade -> species summary lines.
    output_tsv : str or None
        If set, write TSV metadata file with all assemblies.

    Returns
    -------
    Dict[str, List[str]] or None
        Species summary per group if return_dict=True.
    """
    listings: Dict[str, List[str]] = {}
    metadata_rows = []

    groups_to_check = [group_filter] if group_filter else NCBI_GROUPS

    for group in groups_to_check:
        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{group}/assembly_summary.txt"
        try:
            res = requests.get(url, timeout=15)
            res.raise_for_status()
        except Exception as e:
            print(f"[ERROR] Failed to fetch summary for '{group}': {e}")
            continue

        lines = res.text.splitlines()
        species_seen = []
        clade_list = []

        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 20:
                continue

            asm_acc = cols[0]                        # e.g. GCF_000001635.27
            asm_name = cols[15]                      # e.g. GRCm39
            ftp_path = cols[19]                      # ends in /GCF_xxx_ASMname/
            species_name = cols[7]                   # e.g. Mus musculus
            species_binomial = " ".join(species_name.split()[:2])
            taxon_id = cols[5]

            acc_main = asm_acc.split("_")[1].split(".")[0].zfill(9)
            subdir = "/".join([acc_main[i:i+3] for i in range(0, 9, 3)])
            asm_dir = os.path.join(base_dir, "GCF", subdir, asm_acc)
            downloaded = os.path.exists(asm_dir)

            # Filenames from FTP path
            ftp_base = ftp_path.rstrip("/").split("/")[-1]
            gff_url = f"{ftp_path}/{ftp_base}_genomic.gff.gz"
            fa_url = f"{ftp_path}/{ftp_base}_genomic.fna.gz"
            rpt_url = f"{ftp_path}/{ftp_base}_assembly_report.txt"

            # Local paths (only approximate — might not exist)
            gff_local = os.path.join(asm_dir, f"{ftp_base}_genomic.gff.gz")
            fa_local = os.path.join(asm_dir, f"{ftp_base}_genomic.fna.gz")
            rpt_local = os.path.join(asm_dir, f"{ftp_base}_assembly_report.txt")

            metadata_rows.append({
                "group": group,
                "species_name": species_binomial,
                "taxon_id": taxon_id,
                "assembly_accession": asm_acc,
                "assembly_name": asm_name,
                "ftp_path": ftp_path,
                "gff3_ftp": gff_url,
                "fasta_ftp": fa_url,
                "assembly_report_ftp": rpt_url,
                "gff3_local": gff_local,
                "fasta_local": fa_local,
                "assembly_report_local": rpt_local,
                "downloaded": "Downloaded" if downloaded else "Not downloaded"
            })

            clade_list.append(f"{species_binomial} [{'Downloaded' if downloaded else 'Not downloaded'}]")

        # Sort species list alphabetically
        listings[group] = sorted(clade_list)

    # Sort metadata: by group, then species name
    metadata_rows.sort(key=lambda r: (r["group"], r["species_name"]))

    # Write to TSV
    if output_tsv:
        fieldnames = [
            "group", "species_name", "taxon_id", "assembly_accession", "assembly_name",
            "ftp_path", "gff3_ftp", "fasta_ftp", "assembly_report_ftp",
            "gff3_local", "fasta_local", "assembly_report_local",
            "downloaded"
        ]
        with open(output_tsv, "w", newline="") as tsvfile:
            writer = csv.DictWriter(tsvfile, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            for row in metadata_rows:
                writer.writerow(row)
        print(f"\n[✓] Metadata written to {output_tsv}\n")

    # Print summary
    for group, species_list in listings.items():
        print(f"**{group}**: {len(species_list)} species")
        if max_print is None or len(species_list) <= max_print:
            for entry in species_list:
                print(" - " + entry)
        else:
            for entry in species_list[:max_print]:
                print(" - " + entry)
            print(f" ... and {len(species_list) - max_print} more species ...")

    return listings if return_dict else None



def download_annotations(
    base_dir: str = "refseq_data",
    species_name: Optional[str] = None,
    assembly_acc: Optional[str] = None,
    group: Optional[str] = None
):
    """
    Download GFF3, FASTA, and assembly report for a given RefSeq assembly or species.

    Parameters
    ----------
    base_dir : str
        Local root for storing downloaded files.
    species_name : str or None
        Scientific name (e.g. "Mus musculus"). Optional.
    assembly_acc : str or None
        RefSeq accession (e.g. "GCF_000001635.27"). Recommended.
    group : str or None
        RefSeq clade (e.g. "vertebrate_mammalian"). Optional for batch mode.
    """
    def download_assembly(asm_acc: str, ftp_path: str):
        # Parse numeric portion of accession for folder path
        acc_main = asm_acc.split("_")[1].split(".")[0].zfill(9)
        subdir = "/".join([acc_main[i:i+3] for i in range(0, 9, 3)])
        asm_dir = os.path.join(base_dir, "GCF", subdir, asm_acc)
        os.makedirs(asm_dir, exist_ok=True)

        # Assembly name prefix used in file names
        ftp_base = ftp_path.rstrip("/").split("/")[-1]
        gff_url = f"{ftp_path}/{ftp_base}_genomic.gff.gz"
        fa_url = f"{ftp_path}/{ftp_base}_genomic.fna.gz"
        rpt_url = f"{ftp_path}/{ftp_base}_assembly_report.txt"

        # Local paths
        gff_local = os.path.join(asm_dir, os.path.basename(gff_url))
        fa_local = os.path.join(asm_dir, os.path.basename(fa_url))
        rpt_local = os.path.join(asm_dir, os.path.basename(rpt_url))

        for url, dest in [(gff_url, gff_local), (fa_url, fa_local), (rpt_url, rpt_local)]:
            if os.path.exists(dest):
                print(f"  [✓] Already exists: {dest}")
                continue
            print(f"  [↓] Downloading {url} → {dest}")
            try:
                resp = requests.get(url, stream=True, timeout=30)
                resp.raise_for_status()
                with open(dest, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                if os.path.getsize(dest) == 0:
                    raise IOError(f"Downloaded file is empty: {dest}")
            except Exception as e:
                print(f"  [✗] Failed to download {url}: {e}")

    targets = []

    # If assembly accession is given, locate FTP path in all groups
    if assembly_acc:
        found = False
        for grp in (NCBI_GROUPS if not group else [group]):
            url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{grp}/assembly_summary.txt"
            try:
                res = requests.get(url, timeout=20)
                res.raise_for_status()
            except Exception:
                continue
            for line in res.text.splitlines():
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 20:
                    continue
                if cols[0] == assembly_acc:
                    ftp_path = cols[19]  # correct column for FTP path
                    targets.append((assembly_acc, ftp_path))
                    found = True
                    break
            if found:
                break
        if not found:
            raise ValueError(f"Assembly {assembly_acc} not found in RefSeq summaries.")

    # Fallback: locate by species name
    elif species_name:
        found = False
        for grp in (NCBI_GROUPS if not group else [group]):
            url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{grp}/assembly_summary.txt"
            try:
                res = requests.get(url, timeout=20)
                res.raise_for_status()
            except Exception:
                continue
            for line in res.text.splitlines():
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 20:
                    continue
                organism = cols[7]
                if organism.lower().startswith(species_name.lower()):
                    ftp_path = cols[19]
                    asm_acc = cols[0]
                    targets.append((asm_acc, ftp_path))
                    found = True
                    break
            if found:
                break
        if not found:
            raise ValueError(f"Species {species_name} not found in RefSeq summaries.")

    # Batch mode: all assemblies in a group
    elif group:
        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{group}/assembly_summary.txt"
        try:
            res = requests.get(url, timeout=30)
            res.raise_for_status()
        except Exception as e:
            raise RuntimeError(f"Failed to fetch group summary for {group}: {e}")
        for line in res.text.splitlines():
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 20:
                continue
            if cols[10] != "latest":
                continue
            asm_acc = cols[0]
            ftp_path = cols[19]
            targets.append((asm_acc, ftp_path))
    else:
        raise ValueError("Must provide at least one of: assembly_acc, species_name, or group")

    # Download in threads (one per assembly, safe for 1–2 at a time)
    threads = []
    for asm_acc, ftp_path in targets:
        print(f"\n→ Preparing to download {asm_acc}")
        t = threading.Thread(target=download_assembly, args=(asm_acc, ftp_path))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()



def convert_fna_headers(fna_in_path, assembly_report_path, fna_out_path):
    """
    Convert FASTA headers in a genome .fna file to Ensembl-style names
    based on the NCBI assembly report mapping.

    Parameters:
    - fna_in_path: path to input RefSeq-style FASTA
    - assembly_report_path: path to the NCBI assembly report
    - fna_out_path: path to write the converted FASTA
    """

    # Parse the assembly report to build accession → assigned-name map
    accession_to_name = {}
    with open(assembly_report_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 3:
                refseq_acc = cols[6]  # GenBank-Accn or RefSeq-Accn
                assigned = cols[2]    # Assigned-Molecule
                if refseq_acc != "na" and assigned != "na":
                    accession_to_name[refseq_acc] = assigned

    # Process the FASTA
    with open(fna_in_path, "r") as fin, open(fna_out_path, "w") as fout:
        current_name = None
        for line in fin:
            if line.startswith(">"):
                acc = line[1:].split()[0]  # e.g., "NC_000067.7"
                name = accession_to_name.get(acc, acc)
                fout.write(f">{name}\n")
            else:
                fout.write(line)

    print(f"Converted headers written to: {fna_out_path}")



def convert_gff_to_ensembl(gff_path, assembly_report_path, output_path=None, chrom_filter=None):
    """
    Convert a RefSeq-style GFF3 file to Ensembl-style:
     - Replace sequence IDs (RefSeq accessions) with Ensembl region names (chromosome numbers, etc.).
     - Optionally filter to only include certain sequences (e.g., one chromosome for testing).
    If chrom_filter is provided (as a set of desired sequence names or refseq accessions), only those will be included.
    """
    # Load assembly report to build mapping from RefSeq accession -> Ensembl name
    refseq_to_name = {}
    with open(assembly_report_path, "r") as rpt:
        for line in rpt:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            # The assembly report might be tab-delimited with columns:
            # [0] Sequence-Name, [1] Sequence-Role, [2] Assigned-Molecule, [3] Assigned-Molecule-Location/Type,
            # [4] GenBank-Accn, [5] Relationship, [6] RefSeq-Accn, [7] Assembly-Unit, [8] Sequence-Length, [9] UCSC-style-name (if present).
            # Identify the indices by searching header if needed.
            # Here we assume standard columns per NCBI format:contentReference[oaicite:11]{index=11}:contentReference[oaicite:12]{index=12}.
            seq_name = parts[0]
            assigned = parts[2]
            refseq_acc = parts[6]
            # Choose name: use assigned molecule if available, else sequence name
            if assigned != "na":
                new_name = assigned
            else:
                new_name = seq_name
            refseq_to_name[refseq_acc] = new_name
    # Prepare output file path
    if output_path is None:
        # Default output path: replace .gff or .gff.gz with _ensembl.gff3
        base = os.path.splitext(gff_path)[0]
        output_path = base + "_ensembl.gff3"
    # Open GFF (might be gzip) and output
    gff_open = gzip.open if gff_path.endswith(".gz") else open
    with gff_open(gff_path, "rt") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                # Replace sequence IDs in sequence-region directives if present
                if line.lower().startswith("##sequence-region"):
                    # Format: ##sequence-region <seqID> <start> <end>
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        seq_id = parts[1]
                        if seq_id in refseq_to_name:
                            new_id = refseq_to_name[seq_id]
                            parts[1] = new_id
                            line = " ".join(parts) + "\n"
                # Other meta lines we leave unchanged
                outfile.write(line)
            else:
                # GFF feature line
                cols = line.rstrip().split("\t")
                if len(cols) < 9:
                    continue  # skip malformed line
                seq_id = cols[0]
                # Filter by chromosome if needed (for testing scenarios)
                if chrom_filter:
                    # chrom_filter could be a set of allowed seq names (after mapping) or refseq IDs
                    # We accept either in set by checking against both original and new name.
                    allow = False
                    if seq_id in chrom_filter:
                        allow = True
                    elif seq_id in refseq_to_name and refseq_to_name[seq_id] in chrom_filter:
                        allow = True
                    if not allow:
                        continue
                # Replace sequence ID if mapping exists
                if seq_id in refseq_to_name:
                    cols[0] = refseq_to_name[seq_id]
                else:
                    # If no mapping (should not happen for RefSeq sequences in annotation), skip or leave as is.
                    # We'll leave it as is but warn.
                    # print(f"Warning: No mapping for sequence {seq_id}, leaving as is.")
                    pass
                # Write the modified line
                outfile.write("\t".join(cols) + "\n")
    print(f"Converted GFF saved to {output_path}")
    return output_path


def load_seq_regions_from_fna(fna_path, cur, cs_id):
    """
    Load seq_region and dna entries from a converted .fna FASTA file.

    Parameters:
    - fna_path: path to the converted genome FASTA (chrom names already standardized)
    - cur: active MySQL cursor
    - cs_id: coord_system_id for 'primary_assembly'

    Returns:
    - dict: seq_region_name -> seq_region_id
    """
    seq_region_ids = {}
    name = None
    seq_lines = []

    def store_sequence(name, sequence):
        length = len(sequence)
        cur.execute(
            "SELECT seq_region_id FROM seq_region WHERE name = %s AND coord_system_id = %s",
            (name, cs_id)
        )
        row = cur.fetchone()
        if row:
            sr_id = row[0]
            already_exists = True
        else:
            cur.execute(
                "INSERT INTO seq_region (name, coord_system_id, length) VALUES (%s, %s, %s)",
                (name, cs_id, length)
            )
            sr_id = cur.lastrowid
            already_exists = False

        if not already_exists:
            cur.execute(
                "INSERT INTO dna (seq_region_id, sequence) VALUES (%s, %s)",
                (sr_id, sequence)
            )
            # Add toplevel attribute
            cur.execute(
                "INSERT INTO seq_region_attrib (seq_region_id, attrib_type_id, value) VALUES (%s, 6, '')",
                (sr_id,)
            )
        seq_region_ids[name] = sr_id

    with open(fna_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    store_sequence(name, ''.join(seq_lines))
                name = line[1:].split()[0]  # use first token after ">"
                seq_lines = []
            else:
                seq_lines.append(line.upper())
        if name is not None:
            store_sequence(name, ''.join(seq_lines))

    return seq_region_ids


def normalize_id(raw_id: str) -> str:
    """Strip RefSeq prefixes to create stable IDs."""
    return (
        raw_id.replace("gene-", "")
        .replace("rna-", "")
        .replace("cds-", "")
        .replace("exon-", "")
    )


def load_to_ensembl_core(
    converted_gff_path: str,
    converted_fna_path: str,
    assembly_report_path: str,
    species_name: str,
    assembly_accession: str,
    db_host: str,
    db_user: str,
    db_password: str,
    db_port: int = 3306,
    schema_sql_path: Optional[str] = None,
):
    """Load a RefSeq-converted GFF3 into an empty Ensembl core DB.

    Edge-case logic ensures **every** gene ↔ transcript ↔ exon chain exists:
    ─ any transcript w/o exons gets a dummy exon
    ─ any transcript/gene gaps are patched with dummy features
    ─ genes with no transcripts obtain a dummy T/E pair
    Numeric IDs start at 1 for gene / transcript / exon so canonical links
    are valid at insertion time.
    """

    # 0 ─ Connect & prepare schema
    genus, species = species_name.split()[:2]
    db_name = f"{genus.lower()}_{species.lower()}_core_{assembly_accession.split('_')[1].replace('.', '_')}"
    conn = mysql.connector.connect(
        host=db_host, user=db_user, password=db_password, port=db_port, autocommit=False
    )
    cur = conn.cursor()
    cur.execute(f"CREATE DATABASE IF NOT EXISTS {db_name}")
    cur.execute(f"USE {db_name}")

    # Load schema if provided
    if schema_sql_path:
        with open(schema_sql_path) as fh:
            raw = fh.read()
        clean = re.sub(r"/\*\*.*?\*/", "", raw, flags=re.DOTALL)
        for stmt in filter(None, map(str.strip, clean.split(";"))):
            cur.execute(stmt)

    # 1 ─ coord_system, meta, analysis
    cur.execute("INSERT INTO coord_system (name,version,rank,attrib) VALUES ('primary_assembly','',1,'default_version,sequence_level')")
    cs_id = cur.lastrowid
    cur.executemany(
        "INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,%s,%s)",
        [("species.scientific_name", species_name), ("assembly.accession", assembly_accession), ("assembly.name", assembly_accession),
         ("genebuild.level", "toplevel"), ("transcriptbuild.level", "toplevel"), ("exonbuild.level", "toplevel")],
    )
    cur.execute("INSERT INTO analysis (logic_name,created,program) VALUES ('refseq_import',NOW(),'NCBI_RefSeq')")
    analysis_id = cur.lastrowid

    # 2 ─ Load seq_regions
    seq_region_ids = {}
    seq_region_ids = load_seq_regions_from_fna(converted_fna_path, cur, cs_id)

    # 3 ─ Parse GFF3
    genes = {}
    transcripts = {}
    cds_segments = {}

    with open(converted_gff_path) as gff:
        for ln in gff:
            if ln.startswith('#'):
                continue
            seq, src, ftype, s, e, sc, strand, phase, attrs = ln.rstrip().split("\t")
            start, end = int(s), int(e)
            strand_val = 1 if strand == '+' else -1
            ad = dict(kv.split('=', 1) for kv in attrs.split(';') if '=' in kv)

            def resolve_biotype(ftype, ad, fid):
                gbkey = ad.get('gbkey', '')
                pseudo_flag = ad.get('pseudo', 'false') == 'true'
                transcript_biotype = ad.get('transcript_biotype', None)
                gene_biotype = ad.get('gene_biotype', None)

                # Immunoglobulin/T-cell segments remap based on gbkey
                if gbkey.endswith('_segment'):
                    segment_type = gbkey[0]  # D, J, V, C
                    return f"IG_{segment_type}_gene"

                # Transcribed pseudogene heuristic
                if 'Transcribed_Pseudogene' in gbkey:
                    return 'transcribed_pseudogene'

                # Explicit pseudo flag
                if pseudo_flag:
                    return 'pseudogene'

                # If transcript_biotype is explicitly provided, use it
                if transcript_biotype:
                    return transcript_biotype

                # For transcripts, use gbkey as biotype if relevant
                if ftype in ('mRNA', 'transcript', 'lnc_RNA', 'snRNA', 'rRNA', 'snoRNA', 'ncRNA',
                             'antisense_RNA', 'scRNA', 'piRNA', 'siRNA', 'tRNA', 'pseudogenic_tRNA', 'C_region', 'precursor_RNA'):
                    if ftype == 'mRNA':
                        return 'protein_coding'
                    elif ftype == 'C_region':
                        return 'IG_C_gene'
                    else:
                        return ftype

                if ftype == 'exon':
                    if gbkey in ('ncRNA','precursor_RNA'):
                        return gbkey
                    elif gbkey == 'C_region':
                        return 'IG_C_gene'

                # For genes, prefer gene_biotype if available
                if gene_biotype:
                    return gene_biotype

                # Default fallback
                print("DEFAULTING: " + fid + " " + ftype + " " + gbkey)
                return 'protein_coding'

            if ftype == 'gene':
                gid = normalize_id(ad.get('ID', ''))
                dbxref_geneid = next((x.split(':')[1] for x in ad.get('Dbxref', '').split(',') if x.startswith('GeneID:')), None)
                biotype = resolve_biotype(ftype, ad, gid)

                genes[gid] = {
                    'seq_name': seq,
                    'start': start,
                    'end': end,
                    'strand': strand_val,
                    'biotype': biotype,
                    'stable_id': gid,  # Always use ID= for stable_id
                    'xref_geneid': dbxref_geneid,  # Optional reference
                    'name': ad.get('Name') or ad.get('gene') or gid
                }
            elif ftype in ("mRNA", "transcript", "lnc_RNA", "snRNA", "rRNA", "snoRNA", "ncRNA", "antisense_RNA", "scRNA",
                           "piRNA", "siRNA", "tRNA", "pseudogenic_tRNA", "D_gene_segment", "V_gene_segment", "J_gene_segment", "C_gene_segment", "C_region"):
                tid = normalize_id(ad.get('ID', ''))
                pg = normalize_id(ad.get('Parent', ''))
                biotype = resolve_biotype(ftype, ad, tid)

                transcripts[tid] = {
                    'gene_id': pg,
                    'seq_name': seq,
                    'start': start,
                    'end': end,
                    'strand': strand_val,
                    'biotype': biotype,
                    'stable_id': ad.get('Name') or tid,
                    'exons': []
                }

            elif ftype == 'exon':
                pt = normalize_id(ad.get('Parent', ''))
                if pt not in transcripts:
                    # Create dummy gene/transcript if missing
                    if pt not in genes:
                        biotype = resolve_biotype(ftype, ad, pt)
                        genes[pt] = {
                            'seq_name': seq,
                            'start': start,
                            'end': end,
                            'strand': strand_val,
                            'biotype': biotype,
                            'stable_id': pt,
                            'name': pt
                        }
                    dt = f"{pt}_dTx"
                    transcripts.setdefault(dt, {
                        'gene_id': pt,
                        'seq_name': seq,
                        'start': start,
                        'end': end,
                        'strand': strand_val,
                        'biotype': biotype,
                        'stable_id': dt,
                        'exons': []
                    })
                    pt = dt
                transcripts[pt]['exons'].append({'start': start, 'end': end, 'strand': strand_val, 'phase': None, 'end_phase': None})

            elif ftype == 'CDS':
                pt = normalize_id(ad.get('Parent', ''))
                cds_segments.setdefault(pt, []).append((start, end, strand_val, phase))

    # 4 ─ Reconciliation
    for t in transcripts.values():
        if not t['exons']:
            t['exons']=[{'start':t['start'],'end':t['end'],'strand':t['strand'],'phase':None,'end_phase':None}]
    for tid,t in list(transcripts.items()):
        gid=t['gene_id']
        if gid not in genes:
            genes[gid]={'seq_name':t['seq_name'],'start':t['start'],'end':t['end'],'strand':t['strand'],'biotype':t['biotype'],'stable_id':gid,'name':gid}
    for gid,g in genes.items():
        if not any(t['gene_id']==gid for t in transcripts.values()):
            dt=f"{gid}_dTx"
            transcripts[dt]={'gene_id':gid,'seq_name':g['seq_name'],'start':g['start'],'end':g['end'],'strand':g['strand'],'biotype':g['biotype'],'stable_id':dt,'exons':[{'start':g['start'],'end':g['end'],'strand':g['strand'],'phase':None,'end_phase':None}]}


    # ───────────────────────────────────────────────────────────────
    # 4.5 ─ Compute exon.phase & exon.end_phase using CDS (strand‑aware)
    # ───────────────────────────────────────────────────────────────
    for tid, cds_list in cds_segments.items():
        tx = transcripts.get(tid)
        if not tx:
            continue

        strand = tx['strand']

        # Exons in *transcript* order (5' → 3' of transcript)
        exons_tx_order = sorted(tx['exons'], key=lambda x: x['start'], reverse=(strand == -1))

        # CDS segments in transcript order as well
        # cds_list items are (start, end, strand_val, phase_str)
        if strand == 1:
            ordered_cds = sorted(cds_list, key=lambda c: c[0])              # by genomic start
        else:
            ordered_cds = sorted(cds_list, key=lambda c: c[1], reverse=True)  # by genomic end desc

        # First usable GFF3 phase (0/1/2) from the first CDS segment in transcript order
        first_phase = None
        for cs, ce, _s, ph in ordered_cds:
            if ph is not None and ph != '.':
                try:
                    first_phase = int(ph)
                except ValueError:
                    pass
                break

        # Initialise all exons as non-coding
        for ex in tx['exons']:
            ex['phase'] = -1
            ex['end_phase'] = -1

        coding_bases = 0  # total coding bases traversed so far (mod 3 drives phase)
        first_coding_exon_seen = False

        # For each exon in transcript order, work out how many CDS bases overlap it
        for ex in exons_tx_order:
            exon_start, exon_end = ex['start'], ex['end']
            exon_coding_len = 0

            for cs, ce, _s, _ph in ordered_cds:
                ov_start = max(exon_start, cs)
                ov_end   = min(exon_end, ce)
                if ov_start <= ov_end:
                    exon_coding_len += (ov_end - ov_start + 1)

            if exon_coding_len == 0:
                # Non-coding exon, already set to -1/-1
                continue

            # Coding exon: assign phase
            if not first_coding_exon_seen and first_phase is not None:
                ex['phase'] = first_phase
            else:
                ex['phase'] = coding_bases % 3

            coding_bases += exon_coding_len
            ex['end_phase'] = coding_bases % 3
            first_coding_exon_seen = True


    # 4.75 Biotype Overrides: Fix transcript biotypes based on parent gene
    gene_biotype_override_map = {
        'V_segment': 'IG_V_gene',
        'D_segment': 'IG_D_gene',
        'J_segment': 'IG_J_gene',
        'C_segment': 'IG_C_gene',
        'C_region': 'IG_C_gene',
        # Add any other overrides you need here
    }

    for gid, gene in genes.items():
        original_biotype = gene['biotype']
        if original_biotype in gene_biotype_override_map:
            gene['biotype'] = gene_biotype_override_map[original_biotype]


    transcript_biotype_override_map = {
        'lncRNA': 'lncRNA',
        'antisense_RNA': 'antisense_RNA',
        'pseudogene': 'pseudogene',
        'transcribed_pseudogene': 'transcribed_pseudogene',
        'rRNA': 'rRNA',
        'snRNA': 'snRNA',
        'snoRNA': 'snoRNA',
        'tRNA': 'tRNA',
        'miRNA': 'miRNA',
        'ncRNA':'ncRNA',
        'misc_RNA': 'misc_RNA',
        'IG_V_gene': 'IG_V_gene',
        'IG_D_gene': 'IG_D_gene',
        'IG_J_gene': 'IG_J_gene',
        'IG_C_gene': 'IG_C_gene',
        # Add more overrides as needed
    }

    for tid, tx in transcripts.items():
        gene_biotype = genes.get(tx['gene_id'], {}).get('biotype')
        if gene_biotype in transcript_biotype_override_map:
            tx['biotype'] = transcript_biotype_override_map[gene_biotype]

    # ───────────────────────────────────────────────────────────────
    # 5 ─ Allocate numeric IDs
    # ───────────────────────────────────────────────────────────────
    gene_id_map = {}
    tx_id_map   = {}
    exon_id_map = {}

    ng, nt, ne = 1, 1, 1
    for gid in sorted(genes):
        gene_id_map[gid] = ng
        ng += 1
    for tid in sorted(transcripts):
        tx_id_map[tid] = nt
        nt += 1

    # ───────────────────────────────────────────────────────────────
    # 6 ─ Insert genes
    # ───────────────────────────────────────────────────────────────
    for gid, g in genes.items():
        sr_id = seq_region_ids[g['seq_name']]
        # first transcript for this gene (there will always be one after reconciliation)
        first_tx = next((t for t in transcripts if transcripts[t]['gene_id'] == gid), None)
        canonical_tx = tx_id_map.get(first_tx)
        cur.execute(
            """INSERT INTO gene
               (gene_id, seq_region_id, seq_region_start, seq_region_end,
                seq_region_strand, biotype, analysis_id, stable_id,
                display_xref_id, source, canonical_transcript_id)
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,NULL,%s,%s)""",
            (gene_id_map[gid], sr_id, g['start'], g['end'], g['strand'],
             g['biotype'], analysis_id, g['stable_id'], 'refseq', canonical_tx)
        )

    # ───────────────────────────────────────────────────────────────
    # 7 ─ Insert transcripts & exons (rank reversed on − strand)
    # ───────────────────────────────────────────────────────────────
    # We'll keep a per-transcript map from (start,end) -> exon_id to reuse later
    per_tx_coord_to_eid = {}

    for tid, tx in transcripts.items():
        sr_id = seq_region_ids[tx['seq_name']]
        cur.execute(
            """INSERT INTO transcript
               (transcript_id, gene_id, seq_region_id, seq_region_start, seq_region_end,
                seq_region_strand, biotype, analysis_id, stable_id, display_xref_id, source)
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,NULL,%s)""",
            (tx_id_map[tid], gene_id_map[tx['gene_id']], sr_id,
             tx['start'], tx['end'], tx['strand'], tx['biotype'], analysis_id,
             tx['stable_id'], 'refseq')
        )

        # rank in transcript order (reverse for negative strand)
        sorted_exons = sorted(tx['exons'], key=lambda x: x['start'], reverse=(tx['strand'] == -1))
        coord_to_eid = {}
        for rank, ex in enumerate(sorted_exons, start=1):
            eid = ne
            ne += 1
            phase_val     = 0 if ex['phase']     in (None, -1) else ex['phase']
            end_phase_val = 0 if ex['end_phase'] in (None, -1) else ex['end_phase']
            cur.execute(
                """INSERT INTO exon
                   (exon_id, seq_region_id, seq_region_start, seq_region_end,
                    seq_region_strand, phase, end_phase, stable_id)
                   VALUES (%s,%s,%s,%s,%s,%s,%s,NULL)""",
                (eid, sr_id, ex['start'], ex['end'], ex['strand'], phase_val, end_phase_val)
            )
            exon_id_map[(tid, rank)] = eid
            coord_to_eid[(ex['start'], ex['end'])] = eid

            cur.execute(
                "INSERT INTO exon_transcript (exon_id, transcript_id, rank) VALUES (%s,%s,%s)",
                (eid, tx_id_map[tid], rank)
            )

        per_tx_coord_to_eid[tid] = coord_to_eid

    # ───────────────────────────────────────────────────────────────
    # 8 ─ Insert translations (strand-aware start/end & correct ranks)
    # ───────────────────────────────────────────────────────────────
    for tid, cds_list in cds_segments.items():
        db_tx_id = tx_id_map.get(tid)
        if not db_tx_id:
            continue

        strand = transcripts[tid]['strand']
        exons  = transcripts[tid]['exons']

        # CDS extremes in genomic space
        if strand == 1:
            trans_start_pos = min(cs for cs, ce, _s, _p in cds_list)
            trans_end_pos   = max(ce for cs, ce, _s, _p in cds_list)
        else:
            # translation starts at the most 5' genomic coord (which is the larger number on - strand)
            trans_start_pos = max(ce for cs, ce, _s, _p in cds_list)
            trans_end_pos   = min(cs for cs, ce, _s, _p in cds_list)

        # Exons in genomic order (ascending) to find the coords, but
        # ranks were stored in transcript order, so we'll translate later.
        genomic_sorted = sorted(exons, key=lambda e: e['start'])

        start_rank = end_rank = None
        start_off  = end_off  = None

        # helper to convert genomic exon index to stored rank
        # ranks array we used during insert:
        inserted_order = sorted(exons, key=lambda x: x['start'], reverse=(strand == -1))
        rank_map = {(ex['start'], ex['end']): r for r, ex in enumerate(inserted_order, start=1)}

        for ex in genomic_sorted:
            if ex['start'] <= trans_start_pos <= ex['end']:
                # offset is within-exon, 1-based
                start_off = (trans_start_pos - ex['start'] + 1) if strand == 1 else (ex['end'] - trans_start_pos + 1)
                start_rank = rank_map[(ex['start'], ex['end'])]
            if ex['start'] <= trans_end_pos <= ex['end']:
                end_off = (trans_end_pos - ex['start'] + 1) if strand == 1 else (ex['end'] - trans_end_pos + 1)
                end_rank = rank_map[(ex['start'], ex['end'])]

        if start_rank and end_rank and start_off and end_off:
            pid = transcripts[tid].get("protein_id") or f"{tid}_prot"
            cur.execute(
                """INSERT INTO translation
                     (transcript_id, start_exon_id, end_exon_id,
                      seq_start, seq_end, stable_id)
                   VALUES (%s,%s,%s,%s,%s,%s)""",
                (
                    db_tx_id,
                    exon_id_map[(tid, start_rank)],
                    exon_id_map[(tid, end_rank)],
                    start_off,
                    end_off,
                    pid,
                ),
            )
            tr_id = cur.lastrowid
            cur.execute(
                "UPDATE transcript SET canonical_translation_id = %s WHERE transcript_id = %s",
                (tr_id, db_tx_id),
            )

    conn.commit()
    conn.close()

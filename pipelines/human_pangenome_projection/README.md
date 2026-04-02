# Pangenome Annotation Mapping System

A three-stage system for accurately transferring genome annotations from a reference assembly to a target assembly in the context of a species pangenome.

## Overview

This tool maps gene annotations (GFF3 format) from a reference genome to a closely related target assembly. Designed for human pangenome use cases (e.g., GRCh38 → CHM13), it:

1. **Stage 1 - Synteny Detection**: Uses minimap2 whole-genome alignment to identify syntenic blocks
2. **Stage 2 - Feature Mapping**: Projects gene coordinates using exact cs tag parsing with structural validation
3. **Stage 3 - Refinement**: Resolves conflicts, analyzes copy number changes, identifies gaps

Recent redevelopment additions:
- Synteny-aware paralog reassignment for overlapping clusters
- Boundary/codon/splice rescue for method-induced edge errors
- Protein-level QC (identity/coverage/internal-stop checks)
- Second-pass missed-locus recovery from syntenic gap neighborhoods
- Deterministic CNV labeling with ambiguity tags
- Per-gene audit traces in JSON outputs

## Installation

```bash
# Clone repository
cd antigravity_pangenome_mapping_opus

# Install dependencies
pip install -r requirements.txt

# Ensure minimap2 and samtools are available
which minimap2 samtools
```

### Dependencies

- Python 3.9+
- minimap2 (for genome alignment)
- samtools (for FASTA indexing)
- pysam, biopython (Python packages)

## Quick Start

```bash
# Download test data (chr20, ~200MB total)
./scripts/download_test_data.sh data/

# Run mapping
python cli.py map \
    --ref-fasta data/GRCh38.chr20.fa \
    --ref-gff data/gencode.v44.chr20.gff3 \
    --target-fasta data/CHM13v2.chr20.fa \
    --output-gff data/output/chm13.chr20.mapped.gff3 \
    --output-stats data/output/chr20.stats.json \
    --chromosomes chr20 \
    --threads 4
```

## CLI Commands

### `map` - Full Pipeline

```bash
python cli.py map \
    --ref-fasta REFERENCE.fa \
    --ref-gff REFERENCE.gff3 \
    --target-fasta TARGET.fa \
    --output-gff OUTPUT.gff3 \
    --output-stats stats.json \
    [--chromosomes chr1 chr2 ...] \
    [--threads 4] \
    [--min-identity 0.95]
```

Production presets:
```bash
# Very close assemblies (e.g. human pangenome haplotypes)
python cli.py map ... --preset-profile ultra-close

# Same-species but less-close assemblies
python cli.py map ... --preset-profile less-close
```

### `synteny` - Synteny Detection Only

```bash
python cli.py synteny \
    --ref-fasta REFERENCE.fa \
    --target-fasta TARGET.fa \
    --output alignment.paf
```

### `validate` - Validate Mapped Annotations

```bash
python cli.py validate \
    --gff mapped.gff3 \
    --target-fasta TARGET.fa \
    --output validation.json
```

### `summary` - Display Statistics

```bash
python cli.py summary --stats stats.json
```

## Output

### GFF3 Output

Output GFF3 follows Ensembl/GENCODE format with additional provenance attributes:

```
mapped_from=<original_feature_id>
mapping_identity=<0.0-1.0>
mapping_status=<mapped|partial|unmapped>
```

### Statistics (JSON)

Includes:
- Gene/transcript/exon mapping rates
- Splice site validation rates
- Start/stop codon validation rates
- Mapping identity distribution
- Conflict resolution summary
- Copy number changes
- Performance metrics
- Protein QC summary and per-gene audit traces (`*.audit.json`)

## Algorithm

### Stage 1: Synteny Detection

```
Reference FASTA ─┬─► minimap2 -cx asm5 --cs=long ─► PAF file
Target FASTA ────┘                                     │
                                                       ▼
                                        Parse & merge syntenic blocks
                                                       │
                                                       ▼
                                              SyntenicMap (indexed)
```

### Stage 2: Feature Mapping

For each gene:
1. Find overlapping syntenic blocks
2. Project coordinates using cs tag (exact, not linear interpolation)
3. Map exons, CDS, UTRs maintaining hierarchy
4. Validate exon order preservation

### Stage 3: Refinement

1. **Conflict Resolution**: When multiple genes map to overlapping regions (that didn't overlap in reference), select best based on identity and structure
2. **Copy Number**: Track genes with 0 (contraction) or >1 (expansion) target copies
3. **Gap Detection**: Identify large unmapped regions in target
4. **Missed-Locus Recovery**: Re-attempt unmapped/partial genes using local anchor search around syntenic gaps

## Validation

The system validates:
- **Splice sites**: GT-AG (or GC-AG) canonical dinucleotides
- **Start codons**: ATG at CDS start
- **Stop codons**: TAA/TAG/TGA at CDS end  
- **CDS frame**: Length divisible by 3
- **Exon order**: Preserved after mapping

## Testing

```bash
# Run unit tests
pytest tests/ -v

# Run benchmark threshold checks (after producing per-case stats.json files)
python scripts/benchmark_harness.py --manifest benchmarks/manifest.example.json

# Run on test data (chr20)
./scripts/download_test_data.sh data/
python cli.py map --ref-fasta data/GRCh38.chr20.fa ...
```

## File Structure

```
├── cli.py                 # Command-line interface
├── pipeline.py            # Main orchestration
├── src/
│   ├── config.py          # Configuration and constants
│   ├── models.py          # Data models (Gene, Transcript, etc.)
│   ├── gff3_parser.py     # GFF3 input parsing
│   ├── gff3_writer.py     # GFF3 output generation
│   └── fasta_handler.py   # Indexed FASTA access
├── stages/
│   ├── synteny.py         # Stage 1: Synteny detection
│   ├── coordinate_projection.py  # cs tag parsing
│   ├── mapping.py         # Stage 2: Feature mapping
│   └── refinement.py      # Stage 3: Conflict resolution
├── validation/
│   ├── structural.py      # Splice/codon validation
│   ├── protein.py         # Protein-level QC
│   └── statistics.py      # Metrics and reporting
├── scripts/
│   ├── download_test_data.sh
│   └── benchmark_harness.py
├── benchmarks/
│   └── manifest.example.json
├── tests/
│   ├── test_gff3_parser.py
│   └── test_coordinate_projection.py
└── requirements.txt
```

## License

MIT License

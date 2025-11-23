# gene-id-resolver

[![PyPI version](https://badge.fury.io/py/gene-id-resolver.svg)](https://badge.fury.io/py/gene-id-resolver)
[![Python versions](https://img.shields.io/pypi/pyversions/gene-id-resolver.svg)](https://pypi.org/project/gene-id-resolver/)
[![Downloads](https://static.pepy.tech/badge/gene-id-resolver)](https://pepy.tech/project/gene-id-resolver)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17547077.svg)](https://doi.org/10.5281/zenodo.17547077)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Convert gene IDs without the headaches.**

Handle deprecated symbols, ambiguous mappings, and silent failures transparently. Know exactly what happened to every gene in your list.

> Born from the frustration of inconsistent gene mappings

**Key Features:**
- ✅ **Multi-species support:** Human, mouse, and rat genomes
- ✅ **Comprehensive deprecated gene handling:** 145,000+ mappings (HGNC, MGI, RGD)
- ✅ Smart ambiguity resolution (don't silently drop genes)
- ✅ **Database updates:** Stay current with latest Ensembl releases
- ✅ Offline operation (no API limits, cached downloads)
- ✅ Full audit trail (know why conversions fail)
- ✅ Works with TSV, CSV, and plain text files

## Installation

```bash
pip install gene-id-resolver
```

## Quick Start

**Command line:**
```bash
# One-time setup
gene-resolver init

# Convert genes (with professional colored output)
gene-resolver convert TP53 BRCA1 EGFR --to-type ensembl

# Output example:
# ┌─ Conversion Results ──────────────────────────────┐
# │ Genome: GRCh38     │ Strategy: primary  │
# ├───────────────────────────────────────────────────┤
# │ ✓ Successful: 3 │ ⚠ Ambiguous: 0 │ ✗ Failed: 0 │
# └───────────────────────────────────────────────────┘
#
# [SUCCESSFUL CONVERSIONS]
#   TP53 → ENSG00000141510
#   BRCA1 → ENSG00000012048
#   EGFR → ENSG00000146648
```

**Python API:**
```python
from pathlib import Path
from gene_id_resolver.core.resolver import GeneResolver

resolver = GeneResolver(Path("./data"))
result = resolver.convert(
    gene_ids=["TP53", "BRCA1", "SEPT4"],
    from_type="symbol",
    to_type="ensembl",
    genome_build="GRCh38"  # Updated default
)

print(f"Converted: {len(result.successful)} genes")
for gene, mapping in result.successful.items():
    print(f"{gene} → {mapping.identifiers.ensembl_id}")
```

## ✨ New in v0.2.0: Performance & Developer Features

### 🚀 Parallel Processing
Process large gene lists faster with automatic parallelization:
```bash
# Automatic parallel processing for >10 genes
gene-resolver convert gene_list.txt --parallel --max-workers 8

# Force sequential processing
gene-resolver convert genes.txt --sequential
```

### 🔍 Enhanced Debugging
Get detailed conversion information:
```bash
gene-resolver convert TP53 BRCA1 --debug
# Shows: conversion paths, fallback attempts, resolver config
```

### 📊 Structured Output
Export results as JSON for programmatic use:
```bash
gene-resolver convert genes.txt --format json > results.json
```

### ⚙️ Configuration Files
Set defaults for your team or project:
```yaml
# gene_resolver.yaml or ~/.gene_resolver.yaml
species: mus_musculus
genome-build: GRCm39
parallel: true
auto-correct: true
```

### 🛡️ Data Integrity
Automatic validation of downloaded files with corruption detection and recovery.

If you've ever spent hours tracking down why your gene list shrank after ID conversion, or discovered that "MEF2C" maps to multiple Ensembl IDs, this tool is for you. Gene ID conversion seems simple until you encounter:

- **Deprecated symbols** (SEPT4 → SEPTIN4, ENO1L1 → ENO1)
- **Ambiguous genes** (genes with multiple genomic locations)
- **Silent failures** (tools that drop genes without telling you why)
- **Reproducibility issues** (different tools, different results)

This tool handles all these cases transparently, so you know exactly what happened to every gene in your list.

## Detailed Usage

### Command Line Interface

Initialize the database (one-time setup):

```bash
gene-resolver init
```

Convert gene IDs:

```bash
# Convert gene symbols to Ensembl IDs
gene-resolver convert TP53 BRCA1 EGFR --to-type ensembl

# Convert from a file
gene-resolver convert-file genes.txt --output results.csv

# Check database status
gene-resolver status
```

### Python API

```python
from pathlib import Path
from gene_id_resolver.core.resolver import GeneResolver

# Initialize resolver
resolver = GeneResolver(Path("./data"))

# Convert gene symbols to Ensembl IDs
result = resolver.convert(
    gene_ids=["TP53", "BRCA1", "EGFR"],
    from_type="symbol",
    to_type="ensembl",
    genome_build="GRCh38",
    ambiguity_strategy="primary"
)

# Check results
print(f"Successful: {len(result.successful)}")
print(f"Failed: {len(result.failed)}")
print(f"Ambiguous: {len(result.ambiguous)}")

# Get converted IDs
for gene, mapping in result.successful.items():
    print(f"{gene} → {mapping.identifiers.ensembl_id}")

resolver.close()
```

## Features

### Deprecated Gene Handling

Automatically corrects outdated gene symbols:

```bash
gene-resolver convert SEPT4 ENO1L1 CDKN2 --to-type ensembl
# Auto-corrected: SEPT4 → SEPTIN4
# Auto-corrected: ENO1L1 → ENO1
# Auto-corrected: CDKN2 → CDKN2A
```

The tool uses HGNC's complete gene history dataset (58,000+ mappings) to handle symbol changes.

### Multi-Species Support

Work with human, mouse, and rat genomes with species-specific deprecated gene support:

```bash
# Initialize mouse database
gene-resolver init --species mus_musculus

# Initialize rat database
gene-resolver init --species rattus_norvegicus --release 109

# Convert mouse genes
gene-resolver convert Actb Gapdh --to-type ensembl --from-type symbol
```

**Deprecated gene coverage:**
- **Human:** 58,083 mappings from HGNC
- **Mouse:** 77,100 mappings from MGI (Mouse Genome Informatics)
- **Rat:** 10,751 mappings from RGD (Rat Genome Database)

Each species automatically downloads and uses its specific gene history database.

### Database Updates

Keep your gene mappings current with the latest Ensembl releases:

```bash
# Check for updates
gene-resolver update --check

# Update to latest release
gene-resolver update

# Update to specific release
gene-resolver update --release 115

# Update mouse database
gene-resolver update --species mus_musculus
```

Updates are safe - your current database is backed up automatically before updating.

### Ambiguity Resolution

When genes map to multiple locations, you decide what to do:

```bash
# Strategy 1: Prefer protein-coding genes (recommended)
gene-resolver convert AMBIGUOUS_GENE --ambiguity-strategy primary

# Strategy 2: Return all matches for manual review
gene-resolver convert AMBIGUOUS_GENE --ambiguity-strategy all

# Strategy 3: Use first match found
gene-resolver convert AMBIGUOUS_GENE --ambiguity-strategy first

# Strategy 4: Treat as failure
gene-resolver convert AMBIGUOUS_GENE --ambiguity-strategy fail
```

### File Processing

Process entire gene lists from CSV, TSV, or plain text files:

```bash
# From a text file (one gene per line)
gene-resolver convert-file genes.txt --output results.csv

# From a CSV (specify column)
gene-resolver convert-file data.csv --column 0 --output results.csv

# With custom settings
gene-resolver convert-file genes.txt \
    --from-type symbol \
    --to-type ensembl \
    --genome-build GRCh38 \
    --ambiguity-strategy primary \
    --output results.csv
```

### Genome Build Support

- **Human**: hg38/GRCh38, hg19/GRCh37
- **Mouse**: GRCm39, GRCm38
- **Rat**: Rnor_6.0

```bash
gene-resolver convert Tp53 --genome-build GRCm39 --to-type ensembl
```

### Multi-Species Support

Work with human, mouse, or rat genes - each species has its own optimized database:

```bash
# Initialize databases for different species
gene-resolver init --species homo_sapiens --release 109  # Human (default)
gene-resolver init --species mus_musculus --release 109  # Mouse
gene-resolver init --species rattus_norvegicus --release 109  # Rat

# Convert mouse genes
gene-resolver init --species mus_musculus --data-dir ./mouse_data
gene-resolver convert Actb Gapdh --genome-build GRCm39 \
    --data-dir ./mouse_data --to-type ensembl
```

**Python API:**
```python
# Work with mouse genes
mouse_resolver = GeneResolver(Path("./mouse_data"))
mouse_resolver.initialize_database(species="mus_musculus", release="109")

result = mouse_resolver.convert(
    gene_ids=["Actb", "Gapdh"],
    from_type="symbol",
    to_type="ensembl",
    genome_build="GRCm39"
)
```

### Database Updates

Keep your gene annotations up-to-date with Ensembl's latest releases:

```bash
# Check for available updates
gene-resolver update --check

# Update to latest Ensembl release
gene-resolver update

# Update to specific release
gene-resolver update --release 115

# Update mouse database
gene-resolver update --species mus_musculus --data-dir ./mouse_data
```

**Python API:**
```python
# Check for updates
update_info = resolver.check_for_updates()
print(f"Current: Ensembl {update_info['current']['ensembl_release']}")

# Perform update
success = resolver.update_database(target_release="115")
```

> **Note:** Updates download new Ensembl data (~100MB) and rebuild the database. Your current database is automatically backed up with a `.backup` extension.

## Database Information

The tool uses Ensembl annotations (release 109+, auto-updatable) with comprehensive deprecated gene support:

**Gene Coverage:**
- **Human:** 62,710 genes (GRCh38) + 58,083 deprecated mappings (HGNC)
- **Mouse:** 57,010 genes (GRCm39) + 77,100 deprecated mappings (MGI)
- **Rat:** 30,560 genes (mRatBN7.2) + 10,751 deprecated mappings (RGD)

**Total: 145,934 deprecated gene mappings across all species**

**Storage & Performance:**
- Database: `./data/genes.db` (~18 MB per species)
- Gene history: `./data/downloads/gene_history_*.txt.gz` (1-14 MB per species)
- Downloads cached locally - no repeated downloads
- Offline operation after initialization
- Update to latest Ensembl releases anytime

## Advanced Usage

### Python: Batch Processing

```python
from pathlib import Path
from gene_id_resolver.core.resolver_enhanced import EnhancedGeneResolver

# Enhanced resolver with auto-correction
resolver = EnhancedGeneResolver(Path("./data"))

# Convert with automatic deprecated gene correction
result = resolver.convert_with_correction(
    gene_ids=["TP53", "SEPT4", "BRCA1"],
    from_type="symbol",
    to_type="ensembl",
    genome_build="GRCh38",
    ambiguity_strategy="primary"
)

# Check which genes were auto-corrected
if hasattr(result, 'corrections_applied'):
    for old, new in result.corrections_applied.items():
        print(f"Corrected: {old} → {new}")

resolver.close()
```

### Python: Detailed Results

```python
# Examine successful conversions
for gene_id, mapping in result.successful.items():
    print(f"Gene: {gene_id}")
    print(f"  Ensembl ID: {mapping.identifiers.ensembl_id}")
    print(f"  Symbol: {mapping.identifiers.gene_symbol}")
    print(f"  Entrez ID: {mapping.identifiers.entrez_id}")
    print(f"  Biotype: {mapping.biotype}")
    print(f"  Location: {mapping.chromosome}:{mapping.start}-{mapping.end}")

# Check ambiguous genes
for gene_id, mappings in result.ambiguous.items():
    print(f"\n{gene_id} has {len(mappings)} possible matches:")
    for i, mapping in enumerate(mappings):
        print(f"  {i+1}. {mapping.identifiers.ensembl_id} ({mapping.biotype})")

# Examine failures
for failed_gene in result.failed:
    resolution = result.ambiguity_resolutions.get(failed_gene, "not found")
    print(f"Failed: {failed_gene} - {resolution}")
```

## Use Cases

### Cross-Species Comparative Studies
Work seamlessly across human, mouse, and rat genomes with species-specific deprecated gene handling.

### Cancer Research
Convert gene lists from publications to your preferred ID system while tracking deprecated symbols.

### Multi-Omics Integration
Ensure consistent gene identifiers across RNA-seq, proteomics, and methylation datasets.

### Database Maintenance
Keep your gene annotations current with automatic Ensembl updates - no manual downloads needed.

### Reproducible Pipelines
Version-controlled annotations mean your conversions are reproducible across time and platforms.

### Quality Control
Audit trails show exactly which genes failed conversion and why.

## Technical Details

**Built with:**
- Ensembl REST API and FTP (gene annotations)
- SQLite (local database)
- HGNC, MGI, RGD (species-specific gene history)
- Python 3.8+ (pandas, click, tqdm)

**Data sources:**
- **Human:** Ensembl + HGNC (Human Genome Organisation)
- **Mouse:** Ensembl + MGI (Mouse Genome Informatics)
- **Rat:** Ensembl + RGD (Rat Genome Database)

**Design principles:**
- Offline-first (no API rate limits)
- Transparent failures (no silent drops)
- Comprehensive testing (96.9% coverage)
- Clean CLI and Python API

## Contributing

Found a bug or have a feature request? Please open an issue on [GitHub](https://github.com/tahagill/gene-id-resolver/issues).

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```
Ahmad, T. (2025). gene-id-resolver: Production-grade gene ID conversion with deprecated gene support.
https://doi.org/10.5281/zenodo.17547077
```

**BibTeX:**
```bibtex
@software{ahmad2025geneidresolver,
  author       = {Ahmad, Taha},
  title        = {gene-id-resolver: Production-grade gene ID conversion with deprecated gene support},
  year         = 2025,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17547077},
  url          = {https://doi.org/10.5281/zenodo.17547077}
}
```

## Acknowledgments

- Gene annotations from [Ensembl](https://www.ensembl.org/)
- Deprecated gene mappings from [HGNC](https://www.genenames.org/)
- Inspired by the countless hours spent debugging gene ID mismatches

---

**Developed by [Taha Ahmad](https://github.com/tahagill)**
- Bioinformatics-optimized .gitignore
- Modular package structure
- Core dependencies for data engineering phase"


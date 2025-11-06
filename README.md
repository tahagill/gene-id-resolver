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
- ✅ Auto-corrects 58,000+ deprecated gene symbols (SEPT4 → SEPTIN4)
- ✅ Smart ambiguity resolution (don't silently drop genes)
- ✅ Offline database (no API limits, 62,710 genes)
- ✅ Full audit trail (know why conversions fail)
- ✅ Works with human, mouse, and rat genomes

## Installation

```bash
pip install gene-id-resolver
```

## Quick Start

**Command line:**
```bash
# One-time setup
gene-resolver init

# Convert genes
gene-resolver convert TP53 BRCA1 EGFR --to-type ensembl
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
    genome_build="hg38"
)

print(f"Converted: {len(result.successful)} genes")
for gene, mapping in result.successful.items():
    print(f"{gene} → {mapping.identifiers.ensembl_id}")
```

## Why This Tool?

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
    gene_ids=["TP53", "BRCA1", "SEPT4"],  # SEPT4 is deprecated
    from_type="symbol",
    to_type="ensembl",
    genome_build="hg38",
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
    --genome-build hg38 \
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

## Database Information

The tool uses Ensembl annotations (release 109) with:
- **62,710 genes** from human, mouse, and rat genomes
- **58,083 deprecated gene mappings** from HGNC
- **Offline operation** - no API calls needed after initialization

Database is stored locally in `./data/genes.db` (~18 MB) and `./data/gene_history.csv` (~4 MB).

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
    genome_build="hg38",
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

### Cancer Research
Convert gene lists from publications to your preferred ID system while tracking deprecated symbols.

### Multi-Omics Integration
Ensure consistent gene identifiers across RNA-seq, proteomics, and methylation datasets.

### Reproducible Pipelines
Version-controlled annotations mean your conversions are reproducible across time and platforms.

### Quality Control
Audit trails show exactly which genes failed conversion and why.

## Technical Details

**Built with:**
- Ensembl REST API and FTP (data source)
- SQLite (local database)
- HGNC gene history (deprecated symbol tracking)
- Python 3.8+ (pandas, click, tqdm)

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


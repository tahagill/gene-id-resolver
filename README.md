# gene-id-resolver

[![PyPI version](https://badge.fury.io/py/gene-id-resolver.svg)](https://badge.fury.io/py/gene-id-resolver)
[![Python versions](https://img.shields.io/pypi/pyversions/gene-id-resolver.svg)](https://pypi.org/project/gene-id-resolver/)
[![Downloads](https://static.pepy.tech/badge/gene-id-resolver)](https://pepy.tech/project/gene-id-resolver)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**The definitive gene ID conversion tool for reproducible bioinformatics**

> Born from the frustration of inconsistent gene mappings in multi-omics research

## Features

- 🔄 **Reproducible conversions** - version-controlled gene annotations
- 🕵️ **Ambiguity intelligence** - don't silently drop problematic genes  
- 📊 **Complete audit trails** - know exactly why each conversion succeeded/failed
- 💾 **Offline-first** - no API limits, no network dependencies

## Quick Start

```bash
pip install gene-id-resolver

from gene_id_resolver import Resolver

resolver = Resolver(genome_build="hg38", annotation_version="109")
result = resolver.convert(["MEF2C", "BRCA1"], from_type="symbol", to_type="ensembl")
print(result.audit_report)  


For Researchers
This tool solves the critical problem of gene identifier inconsistencies that plague reproducible bioinformatics analysis.

---

## 🔧 **GIT INITIALIZATION - PROPER VERSION CONTROL**

```bash
# Initialize Git
git init

# Add everything except ignored files
git add .

# Professional first commit
git commit -m "feat: initial project structure

- Modern pyproject.toml configuration
- Bioinformatics-optimized .gitignore
- Modular package structure
- Core dependencies for data engineering phase"


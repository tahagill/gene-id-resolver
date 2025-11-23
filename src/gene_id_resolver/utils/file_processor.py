"""
Utilities for processing gene files.
"""

import csv
from pathlib import Path
from typing import List, Tuple
import logging

logger = logging.getLogger(__name__)


def read_genes_from_file(file_path: Path, column: int = 0) -> List[str]:
    """
    Read gene IDs from various file formats.
    
    Supports:
    - CSV files (specify column)
    - TSV files (specify column) 
    - Plain text (one gene per line)
    - BED files (extract gene names from 4th column)
    
    Args:
        file_path: Path to input file
        column: Column index containing gene IDs (0-based)
    
    Returns:
        List of gene IDs
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    # Detect file type by extension
    if file_path.suffix.lower() in ['.csv']:
        return _read_csv(file_path, column)
    elif file_path.suffix.lower() in ['.tsv', '.tab']:
        return _read_tsv(file_path, column)
    elif file_path.suffix.lower() in ['.bed']:
        return _read_bed(file_path)
    else:
        # Assume plain text (one gene per line)
        return _read_plain_text(file_path)


def _read_csv(file_path: Path, column: int) -> List[str]:
    """Read genes from CSV file."""
    genes = []
    with open(file_path, 'r', newline='') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if row and len(row) > column:
                gene = row[column].strip()
                
                # Remove inline comments
                if '#' in gene:
                    gene = gene.split('#')[0].strip()
                
                # Skip empty genes AND header row (first line)
                if gene and not gene.startswith('#') and i > 0:  # ← ADDED i > 0
                    genes.append(gene)
    logger.info(f"Read {len(genes)} genes from CSV file: {file_path}")
    return genes


def _read_tsv(file_path: Path, column: int) -> List[str]:
    """Read genes from TSV file."""
    genes = []
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if row and len(row) > column:
                gene = row[column].strip()
                
                # Remove inline comments
                if '#' in gene:
                    gene = gene.split('#')[0].strip()
                
                # Skip empty genes AND header row (first line)
                if gene and not gene.startswith('#') and i > 0:  # ← ADDED i > 0
                    genes.append(gene)
    logger.info(f"Read {len(genes)} genes from TSV file: {file_path}")
    return genes


def _read_bed(file_path: Path) -> List[str]:
    """Read gene names from BED file (4th column)."""
    genes = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):  # Skip comments and empty lines
                parts = line.split('\t')
                if len(parts) >= 4:  # BED format: chrom, start, end, name, ...
                    gene = parts[3].strip()
                    if gene:
                        genes.append(gene)
    logger.info(f"Read {len(genes)} genes from BED file: {file_path}")
    return genes


def _read_plain_text(file_path: Path) -> List[str]:
    """Read genes from plain text file (one per line)."""
    genes = []
    with open(file_path, 'r') as f:
        for line in f:
            gene = line.strip()
            
            # Skip empty lines and full-line comments
            if not gene or gene.startswith('#'):
                continue
            
            # Remove inline comments (e.g., "SEPT4  # comment")
            if '#' in gene:
                gene = gene.split('#')[0].strip()
            
            if gene:  # Add only if not empty after comment removal
                genes.append(gene)
    logger.info(f"Read {len(genes)} genes from text file: {file_path}")
    return genes


def save_results_to_file(result, output_path: Path, input_genes: List[str] = None):
    """
    Save conversion results to a CSV file.
    
    Args:
        result: ConversionResult object
        output_path: Path to output file
        input_genes: Original input genes (for preserving order)
    """
    output_path = Path(output_path)
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['input_gene', 'ensembl_id', 'gene_symbol', 'entrez_id', 'status'])
        
        # Process genes in input order if provided
        genes_to_process = input_genes if input_genes else list(result.successful.keys()) + result.failed
        
        for gene in genes_to_process:
            if gene in result.successful:
                mapping = result.successful[gene]
                writer.writerow([
                    gene,
                    mapping.identifiers.ensembl_id or '',
                    mapping.identifiers.gene_symbol or '', 
                    mapping.identifiers.entrez_id or '',
                    'success'
                ])
            elif gene in result.ambiguous:
                # For ambiguous genes, write the first match
                mapping = result.ambiguous[gene][0]
                writer.writerow([
                    gene,
                    mapping.identifiers.ensembl_id or '',
                    mapping.identifiers.gene_symbol or '',
                    mapping.identifiers.entrez_id or '',
                    'ambiguous'
                ])
            else:
                writer.writerow([gene, '', '', '', 'failed'])
    
    logger.info(f"Saved conversion results to: {output_path}")
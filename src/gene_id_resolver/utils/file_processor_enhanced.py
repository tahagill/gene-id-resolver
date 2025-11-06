"""
Enhanced file processor with better error handling and encoding detection.
"""

import csv
import chardet  # You'll need to pip install chardet
from pathlib import Path
from typing import List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)

def detect_encoding(file_path: Path) -> str:
    """Detect file encoding automatically."""
    with open(file_path, 'rb') as f:
        raw_data = f.read(10000)  # Sample first 10KB
        result = chardet.detect(raw_data)
        encoding = result['encoding'] or 'utf-8'
        logger.info(f"Detected encoding: {encoding} (confidence: {result['confidence']:.2f})")
        return encoding

def read_genes_from_file_enhanced(file_path: Path, column: int = 0, 
                                has_header: bool = True) -> Tuple[List[str], List[str]]:
    """
    Enhanced file reading with better error handling.
    
    Returns:
        Tuple of (genes, warnings)
    """
    file_path = Path(file_path)
    genes = []
    warnings = []
    
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    try:
        encoding = detect_encoding(file_path)
        
        with open(file_path, 'r', encoding=encoding, newline='') as f:
            # Peek at first line to detect format
            first_line = f.readline()
            f.seek(0)  # Reset to beginning
            
            is_tsv = '\t' in first_line
            delimiter = '\t' if is_tsv else ','
            
            reader = csv.reader(f, delimiter=delimiter)
            
            for i, row in enumerate(reader):
                if not row:  # Skip empty rows
                    continue
                    
                if has_header and i == 0:
                    logger.info(f"Header detected: {row}")
                    continue
                
                if len(row) <= column:
                    warnings.append(f"Row {i+1}: No column {column}, available columns: {len(row)}")
                    continue
                
                gene = row[column].strip()
                
                # Skip comments and empty genes
                if not gene or gene.startswith('#'):
                    continue
                
                # Clean common issues
                gene = _clean_gene_name(gene)
                
                if gene:
                    genes.append(gene)
                else:
                    warnings.append(f"Row {i+1}: Empty gene after cleaning")
        
    except UnicodeDecodeError as e:
        logger.error(f"Encoding error in {file_path}: {e}")
        # Fall back to simple text reading
        genes = _read_with_fallback(file_path)
        warnings.append("Used fallback reading due to encoding issues")
    
    logger.info(f"Read {len(genes)} genes from {file_path} with {len(warnings)} warnings")
    return genes, warnings

def _clean_gene_name(gene: str) -> str:
    """Clean common gene name issues."""
    # Remove quotes, extra spaces, etc.
    gene = gene.strip('"\'').strip()
    
    # Remove version numbers (e.g., "Gene.1" -> "Gene")
    if '.' in gene and gene.split('.')[-1].isdigit():
        gene = gene.rsplit('.', 1)[0]
    
    return gene

def _read_with_fallback(file_path: Path) -> List[str]:
    """Fallback reading for problematic files."""
    genes = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.append(gene)
    return genes
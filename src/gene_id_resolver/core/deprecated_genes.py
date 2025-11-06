"""
Smart deprecated gene detection using Ensembl's official gene history.
"""

import sqlite3
import gzip
from pathlib import Path
from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)

class SmartDeprecatedGeneHandler:
    """Uses Ensembl's gene history for comprehensive deprecated gene detection."""
    
    def __init__(self, db_path: Path, data_dir: Path):
        self.db_path = db_path
        self.data_dir = data_dir
        self._deprecated_map: Dict[str, str] = {}
        self._load_gene_history()
    
    def _load_gene_history(self):
        """Load deprecated genes from HGNC gene symbol history file."""
        gene_history_file = self.data_dir / "downloads" / "gene_history_homo_sapiens_109.txt.gz"
        
        if not gene_history_file.exists():
            logger.warning("Gene history file not found. Using fallback mappings.")
            self._load_fallback_mappings()
            return
        
        try:
            import csv
            
            with gzip.open(gene_history_file, 'rt') as f:
                reader = csv.DictReader(f, delimiter='\t')
                
                for row in reader:
                    # HGNC format: current approved symbol is in 'symbol' column
                    # Previous symbols are in 'prev_symbol' column (pipe-separated)
                    # Alias symbols are in 'alias_symbol' column (pipe-separated)
                    
                    current_symbol = row.get('symbol', '').strip()
                    
                    if not current_symbol:
                        continue
                    
                    # Parse previous symbols
                    prev_symbols = row.get('prev_symbol', '')
                    if prev_symbols:
                        for old_symbol in prev_symbols.split('|'):
                            old_symbol = old_symbol.strip()
                            if old_symbol and old_symbol != current_symbol:
                                self._deprecated_map[old_symbol.upper()] = current_symbol
                    
                    # Also map common aliases
                    alias_symbols = row.get('alias_symbol', '')
                    if alias_symbols:
                        for alias in alias_symbols.split('|'):
                            alias = alias.strip()
                            if alias and alias != current_symbol:
                                # Only map if not already mapped (previous symbols take priority)
                                if alias.upper() not in self._deprecated_map:
                                    self._deprecated_map[alias.upper()] = current_symbol
            
            logger.info(f"Loaded {len(self._deprecated_map)} deprecated/alias gene mappings from HGNC")
            
        except Exception as e:
            logger.error(f"Failed to load gene history: {e}")
            self._load_fallback_mappings()
    
    def _load_fallback_mappings(self):
        """Fallback to common known deprecated mappings."""
        common_deprecated = {
            "SEPT4": "SEPTIN4",
            "ENO1L1": "ENO1",
            "CDKN2": "CDKN2A",  # Common ambiguous symbol
            "CDKN2B": "CDKN2B",  # Keep for reference
        }
        self._deprecated_map = common_deprecated
        logger.info(f"Using fallback mappings: {len(self._deprecated_map)} genes")
    
    def get_current_symbol(self, old_symbol: str) -> Optional[str]:
        """Get current symbol for deprecated gene."""
        return self._deprecated_map.get(old_symbol.upper())
    
    def is_deprecated(self, symbol: str) -> bool:
        """Check if a gene symbol is deprecated."""
        return symbol.upper() in self._deprecated_map
    
    def suggest_corrections(self, failed_genes: List[str]) -> Dict[str, str]:
        """Suggest corrections for failed gene conversions."""
        suggestions = {}
        for gene in failed_genes:
            current = self.get_current_symbol(gene)
            if current:
                suggestions[gene] = current
        return suggestions
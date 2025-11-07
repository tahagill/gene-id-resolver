"""
Smart deprecated gene detection using official gene history databases.
Supports human (HGNC), mouse (MGI), and rat (RGD).
"""

import sqlite3
import gzip
from pathlib import Path
from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)

class SmartDeprecatedGeneHandler:
    """Uses species-specific gene history databases for comprehensive deprecated gene detection."""
    
    def __init__(self, db_path: Path, data_dir: Path, species: str = "homo_sapiens"):
        self.db_path = db_path
        self.data_dir = data_dir
        self.species = species
        self._deprecated_map: Dict[str, str] = {}
        self._load_gene_history()
    
    def _load_gene_history(self):
        """Load deprecated genes from species-specific gene history file."""
        gene_history_file = self.data_dir / "downloads" / f"gene_history_{self.species}_109.txt.gz"
        
        if not gene_history_file.exists():
            logger.warning(f"Gene history file not found for {self.species}. Using fallback mappings.")
            self._load_fallback_mappings()
            return
        
        try:
            if self.species == "homo_sapiens":
                self._load_hgnc_history(gene_history_file)
            elif self.species == "mus_musculus":
                self._load_mgi_history(gene_history_file)
            elif self.species == "rattus_norvegicus":
                self._load_rgd_history(gene_history_file)
            else:
                logger.warning(f"Unknown species: {self.species}")
                self._load_fallback_mappings()
                
        except Exception as e:
            logger.error(f"Failed to load gene history: {e}")
            self._load_fallback_mappings()
    
    def _load_hgnc_history(self, file_path: Path):
        """Load human gene history from HGNC."""
        import csv
        
        with gzip.open(file_path, 'rt') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
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
                            if alias.upper() not in self._deprecated_map:
                                self._deprecated_map[alias.upper()] = current_symbol
        
        logger.info(f"Loaded {len(self._deprecated_map)} deprecated/alias gene mappings from HGNC")
    
    def _load_mgi_history(self, file_path: Path):
        """Load mouse gene history from MGI."""
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if line.startswith('#') or line.startswith('MGI Accession'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 14:
                    continue
                
                # Column 7 (index 6): Marker Symbol
                # Column 8 (index 7): Status (W = withdrawn, O = official)
                # Column 14 (index 13): Current Marker Symbol (if withdrawn)
                
                status = fields[7].strip()
                if status == 'W':  # Withdrawn
                    old_symbol = fields[6].strip()
                    new_symbol = fields[13].strip() if len(fields) > 13 else None
                    
                    if old_symbol and new_symbol:
                        self._deprecated_map[old_symbol.upper()] = new_symbol
        
        logger.info(f"Loaded {len(self._deprecated_map)} withdrawn gene mappings from MGI")
    
    def _load_rgd_history(self, file_path: Path):
        """Load rat gene history from RGD."""
        import csv
        
        with gzip.open(file_path, 'rt') as f:
            # Skip comment lines starting with #
            for line in f:
                if line.startswith('SPECIES'):
                    # This is the header line
                    break
            
            # Now read the rest as tab-delimited
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 7:
                    continue
                
                # Column 0: SPECIES
                # Column 2: OLD_GENE_SYMBOL
                # Column 6: NEW_GENE_SYMBOL
                
                species = fields[0].strip()
                if species != 'rat':
                    continue
                
                old_symbol = fields[2].strip()
                new_symbol = fields[6].strip() if len(fields) > 6 else ''
                
                if old_symbol and new_symbol:
                    self._deprecated_map[old_symbol.upper()] = new_symbol
        
        logger.info(f"Loaded {len(self._deprecated_map)} obsolete gene mappings from RGD")
    
    def _load_fallback_mappings(self):
        """Fallback to common known deprecated mappings for each species."""
        if self.species == "homo_sapiens":
            common_deprecated = {
                "SEPT4": "SEPTIN4",
                "ENO1L1": "ENO1",
                "CDKN2": "CDKN2A",
            }
        elif self.species == "mus_musculus":
            common_deprecated = {
                "GSTM7": "0610005A07RIK",
                "LYPD2": "0610005K03RIK",
            }
        elif self.species == "rattus_norvegicus":
            common_deprecated = {
                "ABL1": "ABL1_MAPPED",
                "ACAA": "ACAA1A",
            }
        else:
            common_deprecated = {}
        
        self._deprecated_map = {k.upper(): v for k, v in common_deprecated.items()}
        logger.info(f"Using fallback mappings for {self.species}: {len(self._deprecated_map)} genes")
    
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
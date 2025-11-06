"""
Main gene ID resolution engine.

This ties together the database and data sources to provide
a clean API for converting between gene identifier systems.
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional
from .ambiguity import AmbiguityStrategy, resolve_ambiguity
from .models import GeneIdentifier, GeneMapping, ConversionResult
from .database import GeneDatabase
from ..datasources.ensembl import EnsemblDownloader

logger = logging.getLogger(__name__)


class GeneResolver:
    """
    Main class for gene ID conversion operations.
    
    This is the primary interface that users will interact with.
    It manages the database and coordinates between different
    data sources.
    """
    
    def __init__(self, data_dir: Path = Path("data")):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(exist_ok=True)
        
        # Initialize components
        self.db = GeneDatabase(self.data_dir / "genes.db")
        self.downloader = EnsemblDownloader(self.data_dir / "downloads")
        
        logger.info(f"GeneResolver initialized with data directory: {self.data_dir}")
    
    def initialize_database(self, species: str = "homo_sapiens", 
                            release: str = "109",
                            genome_build: str = None) -> int:
        """
        Download Ensembl data and populate the database.
        """
        logger.info(f"Initializing database with Ensembl {release} {species}")
        
        # Download the GTF file
        gtf_path = self.downloader.download_gtf(species, release, genome_build)
        
        # NEW: Download gene history file for deprecated genes
        try:
            history_path = self.downloader.download_gene_history(species, release)
            logger.info(f"Gene history downloaded: {history_path}")
        except Exception as e:
            logger.warning(f"Could not download gene history: {e}")
        
        # Parse and load into database
        gene_count = 0
        for gene_mapping in self.downloader.parse_gtf(gtf_path):
            success = self.db.insert_mapping(gene_mapping)
            if success:
                gene_count += 1
        
        logger.info(f"Database initialized with {gene_count} genes")
        return gene_count
        
    def convert(self, gene_ids: List[str], from_type: str, to_type: str, 
            genome_build: str = "hg38",
            ambiguity_strategy: str = "primary") -> ConversionResult:  # ← CHANGE PARAMETER NAME
        """
        Convert gene IDs from one system to another.
        
        Args:
            gene_ids: List of gene identifiers to convert
            from_type: Input ID type (symbol, ensembl, entrez)
            to_type: Output ID type (symbol, ensembl, entrez) 
            genome_build: Genome build (GRCh38, GRCh37, GRCm39, etc.)
            ambiguity_strategy: How to handle ambiguous mappings  # ← UPDATED
                            (primary, first, all, fail)
        """
        # Convert string to enum
        try:
            strategy_enum = AmbiguityStrategy(ambiguity_strategy)  # ← DIFFERENT VARIABLE NAME
        except ValueError:
            raise ValueError(f"Invalid ambiguity strategy: {ambiguity_strategy}. "
                        f"Valid options: {[s.value for s in AmbiguityStrategy]}")
        
        result = ConversionResult(
            input_type=from_type,
            output_type=to_type,
            resolver_config={
                "genome_build": genome_build,
                "ambiguity_strategy": ambiguity_strategy  # ← UPDATED
            }
        )
        
        logger.info(f"Converting {len(gene_ids)} genes from {from_type} to {to_type} "
                f"(build: {genome_build}, ambiguity: {ambiguity_strategy})")  # ← UPDATED
        
        for gene_id in gene_ids:
            mappings = self._find_gene_mappings(gene_id, from_type, genome_build)
            
            if not mappings:
                result.add_failed(gene_id)
                result.ambiguity_resolutions[gene_id] = "failed"
                
            elif len(mappings) == 1:
                result.add_success(gene_id, mappings[0])
                result.ambiguity_resolutions[gene_id] = "unique"
                
            else:
                # AMBIGUOUS CASE - Apply resolution strategy
                resolved_mappings = resolve_ambiguity(mappings, strategy_enum)  # ← USE strategy_enum
                
                if strategy_enum == AmbiguityStrategy.ALL:
                    # Return all matches for manual resolution
                    result.add_ambiguous(gene_id, mappings)
                    result.ambiguity_resolutions[gene_id] = "unresolved"
                    
                elif resolved_mappings:
                    # Strategy resolved the ambiguity
                    result.add_success(gene_id, resolved_mappings[0])
                    if strategy_enum == AmbiguityStrategy.PRIMARY:
                        result.ambiguity_resolutions[gene_id] = "resolved_primary"
                    else:
                        result.ambiguity_resolutions[gene_id] = "resolved_first"
                else:
                    # Strategy says to fail (AmbiguityStrategy.FAIL)
                    result.add_failed(gene_id)
                    result.ambiguity_resolutions[gene_id] = "failed_ambiguous"
        
        logger.info(f"Conversion complete: {len(result.successful)} successful, "
                f"{len(result.ambiguous)} ambiguous, {len(result.failed)} failed")
        
        return result

    def _find_gene_mappings(self, gene_id: str, id_type: str, genome_build: str) -> List[GeneMapping]:
        """
        Find gene mappings for a given identifier and genome build.
        """
        # Normalize the ID type
        id_type = id_type.lower().strip()
        
        try:
            if id_type in ['ensembl', 'ensembl_id', 'ensemblid']:
                return self.db.find_by_ensembl(gene_id, genome_build)
            
            elif id_type in ['symbol', 'gene_symbol', 'genesymbol', 'name']:
                return self.db.find_by_symbol(gene_id, genome_build)
            
            elif id_type in ['entrez', 'entrez_id', 'entrezid', 'ncbi']:
                return self.db.find_by_entrez(gene_id, genome_build)
            
            else:
                raise ValueError(f"Unsupported ID type: {id_type}. "
                               f"Valid types: symbol, ensembl, entrez")
                
        except ValueError:
            # Re-raise ValueError (unsupported ID type)
            raise
        except Exception as e:
            logger.error(f"Error searching for {id_type}: {gene_id} (build: {genome_build}): {e}")
            return []
    def close(self):
        """Clean up resources."""
        self.db.close()
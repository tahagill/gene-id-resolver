"""
Main gene ID resolution engine.

This ties together the database and data sources to provide
a clean API for converting between gene identifier systems.
"""

import logging
import concurrent.futures
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
        
        # Set default genome builds
        if not genome_build:
            defaults = {
                "homo_sapiens": "GRCh38",
                "mus_musculus": "GRCm39",
                "rattus_norvegicus": "mRatBN7.2"
            }
            genome_build = defaults.get(species, "GRCh38")
        
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
        for gene_mapping in self.downloader.parse_gtf(gtf_path, genome_build, release):
            success = self.db.insert_mapping(gene_mapping)
            if success:
                gene_count += 1
        
        # Update metadata
        self.db.update_metadata(
            ensembl_release=release,
            species=species,
            genome_build=genome_build
        )
        
        logger.info(f"Database initialized with {gene_count} genes")
        return gene_count
        
    def convert(self, gene_ids: List[str], from_type: str, to_type: str, 
            genome_build: str = "hg38",
            ambiguity_strategy: str = "primary",
            parallel: bool = True,
            max_workers: int = 4) -> ConversionResult:  # ← CHANGE PARAMETER NAME
        """
        Convert gene IDs from one system to another.
        
        Args:
            gene_ids: List of gene identifiers to convert
            from_type: Input ID type (symbol, ensembl, entrez)
            to_type: Output ID type (symbol, ensembl, entrez) 
            genome_build: Genome build (GRCh38, GRCh37, GRCm39, etc.)
            ambiguity_strategy: How to handle ambiguous mappings  # ← UPDATED
                            (primary, first, all, fail)
            parallel: Whether to use parallel processing for large lists
            max_workers: Maximum number of parallel workers
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
        
        # Use parallel processing for large gene lists
        if parallel and len(gene_ids) > 10:
            return self._convert_parallel(gene_ids, from_type, to_type, genome_build, 
                                        strategy_enum, max_workers)
        else:
            return self._convert_sequential(gene_ids, from_type, to_type, genome_build, 
                                          strategy_enum)

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
    
    def check_for_updates(self) -> Dict:
        """Check if database updates are available."""
        from .updater import DatabaseUpdater
        updater = DatabaseUpdater(self.data_dir / "genes.db")
        return updater.check_for_updates()
    
    def update_database(self, target_release: str, species: str = "homo_sapiens",
                       genome_build: str = None) -> bool:
        """Update database to a new Ensembl release."""
        from .updater import DatabaseUpdater
        updater = DatabaseUpdater(self.data_dir / "genes.db")
        return updater.incremental_update(target_release, species, genome_build)
    
    def _convert_sequential(self, gene_ids: List[str], from_type: str, to_type: str,
                           genome_build: str, strategy_enum: AmbiguityStrategy) -> ConversionResult:
        """Convert genes sequentially (original implementation)."""
        result = ConversionResult(
            input_type=from_type,
            output_type=to_type,
            resolver_config={
                "genome_build": genome_build,
                "ambiguity_strategy": strategy_enum.value
            }
        )
        
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
                resolved_mappings = resolve_ambiguity(mappings, strategy_enum)
                
                if strategy_enum == AmbiguityStrategy.ALL:
                    result.add_ambiguous(gene_id, mappings)
                    result.ambiguity_resolutions[gene_id] = "unresolved"
                    
                elif resolved_mappings:
                    result.add_success(gene_id, resolved_mappings[0])
                    if strategy_enum == AmbiguityStrategy.PRIMARY:
                        result.ambiguity_resolutions[gene_id] = "resolved_primary"
                    else:
                        result.ambiguity_resolutions[gene_id] = "resolved_first"
                else:
                    result.add_failed(gene_id)
                    result.ambiguity_resolutions[gene_id] = "failed_ambiguous"
        
        return result
    
    def _convert_parallel(self, gene_ids: List[str], from_type: str, to_type: str,
                         genome_build: str, strategy_enum: AmbiguityStrategy, 
                         max_workers: int) -> ConversionResult:
        """Convert genes in parallel using thread pools."""
        result = ConversionResult(
            input_type=from_type,
            output_type=to_type,
            resolver_config={
                "genome_build": genome_build,
                "ambiguity_strategy": strategy_enum.value,
                "parallel": True,
                "max_workers": max_workers
            }
        )

        # STEP 1: Fetch all gene mappings in main thread (SQLite-safe)
        logger.info(f"Fetching mappings for {len(gene_ids)} genes...")
        gene_mappings = {}
        for gene_id in gene_ids:
            mappings = self._find_gene_mappings(gene_id, from_type, genome_build)
            gene_mappings[gene_id] = mappings

        # STEP 2: Process mappings in parallel (no database calls)
        logger.info(f"Processing mappings with {max_workers} workers...")

        # Split genes into batches for parallel processing
        batch_size = max(1, len(gene_ids) // max_workers)
        gene_batches = [gene_ids[i:i + batch_size]
                       for i in range(0, len(gene_ids), batch_size)]

        # Process batches in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all batches for processing
            future_to_batch = {
                executor.submit(self._process_mappings_batch, batch, gene_mappings, strategy_enum): batch
                for batch in gene_batches
            }

            # Collect results as they complete
            for future in concurrent.futures.as_completed(future_to_batch):
                batch_result = future.result()

                # Merge batch results into main result
                result.successful.update(batch_result.successful)
                result.ambiguous.update(batch_result.ambiguous)
                result.failed.extend(batch_result.failed)
                result.ambiguity_resolutions.update(batch_result.ambiguity_resolutions)

        logger.info(f"Parallel conversion complete: {len(result.successful)} successful, "
                   f"{len(result.ambiguous)} ambiguous, {len(result.failed)} failed")

        return result
    
    def _process_gene_batch(self, gene_batch: List[str], from_type: str, 
                           genome_build: str, strategy_enum: AmbiguityStrategy) -> ConversionResult:
        """Process a batch of genes (called by parallel executor)."""
        # Create a temporary result for this batch
        batch_result = ConversionResult(
            input_type=from_type,
            output_type="temp",  # Not used in batch processing
            resolver_config={}
        )
        
        for gene_id in gene_batch:
            mappings = self._find_gene_mappings(gene_id, from_type, genome_build)
            
            if not mappings:
                batch_result.add_failed(gene_id)
                batch_result.ambiguity_resolutions[gene_id] = "failed"
                
            elif len(mappings) == 1:
                batch_result.add_success(gene_id, mappings[0])
                batch_result.ambiguity_resolutions[gene_id] = "unique"
                
            else:
                # AMBIGUOUS CASE - Apply resolution strategy
                resolved_mappings = resolve_ambiguity(mappings, strategy_enum)
                
                if strategy_enum == AmbiguityStrategy.ALL:
                    batch_result.add_ambiguous(gene_id, mappings)
                    batch_result.ambiguity_resolutions[gene_id] = "unresolved"
                    
                elif resolved_mappings:
                    batch_result.add_success(gene_id, resolved_mappings[0])
                    if strategy_enum == AmbiguityStrategy.PRIMARY:
                        batch_result.ambiguity_resolutions[gene_id] = "resolved_primary"
                    else:
                        batch_result.ambiguity_resolutions[gene_id] = "resolved_first"
                else:
                    batch_result.add_failed(gene_id)
                    batch_result.ambiguity_resolutions[gene_id] = "failed_ambiguous"
        
        return batch_result

    def _process_mappings_batch(self, gene_batch: List[str], gene_mappings: Dict[str, List[GeneMapping]],
                               strategy_enum: AmbiguityStrategy) -> ConversionResult:
        """Process a batch of gene mappings in parallel (no database calls)."""
        batch_result = ConversionResult()

        for gene_id in gene_batch:
            mappings = gene_mappings.get(gene_id, [])

            if not mappings:
                batch_result.failed.append(gene_id)
                continue

            # Handle ambiguity resolution
            if len(mappings) > 1:
                resolved_mapping = self._resolve_ambiguity(mappings, strategy_enum)
                if resolved_mapping:
                    batch_result.successful[gene_id] = resolved_mapping
                    batch_result.ambiguity_resolutions[gene_id] = resolved_mapping
                else:
                    batch_result.ambiguous[gene_id] = mappings
            else:
                # Single mapping - use it directly
                mapping = mappings[0]
                batch_result.successful[gene_id] = mapping

        return batch_result

    def close(self):
        """Clean up resources."""
        self.db.close()
"""
SQLite database layer for storing gene mappings.

Using SQLite because it's fast, reliable, and doesn't require a server.
We optimize for read performance since most operations are queries.
"""

import sqlite3
import json
import logging
from pathlib import Path
from typing import List, Optional, Dict, Any

from .models import GeneMapping, GeneIdentifier, GenomicCoordinates

logger = logging.getLogger(__name__)


class GeneDatabase:
    """
    Handles all database operations for gene mappings.
    
    We use SQLite's WAL mode for better concurrent read performance.
    The indexes are tuned for the types of queries we'll run most often.
    """
    
    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.connection = None
        self._setup_database()
    
    def _setup_database(self):
        """Initialize the database with proper settings and tables."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Connect and configure for performance
        self.connection = sqlite3.connect(self.db_path)
        self.connection.row_factory = sqlite3.Row  # Get dict-like rows
        
        # Performance tweaks - makes a big difference with large datasets
        self.connection.execute("PRAGMA journal_mode=WAL")  # Better for concurrent reads
        self.connection.execute("PRAGMA synchronous=NORMAL")  # Good balance of speed/safety
        self.connection.execute("PRAGMA cache_size=-64000")  # 64MB cache
        self.connection.execute("PRAGMA foreign_keys=ON")
        
        self._create_tables()
        logger.info(f"Database ready at {self.db_path}")
    
    def _create_tables(self):
        """Create the main table and indexes."""
        
        # Main table with all gene mappings
        # Using separate columns for each ID type makes querying faster
        self.connection.execute("""
            CREATE TABLE IF NOT EXISTS gene_mappings (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                
                -- Different gene identifier systems
                ensembl_id TEXT,
                gene_symbol TEXT, 
                entrez_id TEXT,
                uniprot_id TEXT,
                
                -- Genomic location
                chromosome TEXT,
                start INTEGER,
                end INTEGER,
                strand TEXT,
                
                -- Gene metadata
                biotype TEXT,
                description TEXT,
                aliases_json TEXT,  -- Store list as JSON
                
                -- Version info for reproducibility
                genome_build TEXT NOT NULL,
                annotation_version TEXT NOT NULL,
                data_source TEXT NOT NULL,
                import_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                
                -- Status flags
                is_current BOOLEAN DEFAULT TRUE,
                is_primary BOOLEAN DEFAULT TRUE,
                
                -- Unique hash to prevent duplicates
                content_hash TEXT UNIQUE,
                
                -- Must have at least one identifier
                CHECK (ensembl_id IS NOT NULL OR gene_symbol IS NOT NULL OR 
                       entrez_id IS NOT NULL OR uniprot_id IS NOT NULL)
            )
        """)
        
        # Create indexes on columns we'll query frequently
        # These make a huge difference with 50K+ gene records
        index_queries = [
            "CREATE INDEX IF NOT EXISTS idx_ensembl ON gene_mappings(ensembl_id)",
            "CREATE INDEX IF NOT EXISTS idx_symbol ON gene_mappings(gene_symbol)",
            "CREATE INDEX IF NOT EXISTS idx_entrez ON gene_mappings(entrez_id)", 
            "CREATE INDEX IF NOT EXISTS idx_content_hash ON gene_mappings(content_hash)",
            "CREATE INDEX IF NOT EXISTS idx_build_version ON gene_mappings(genome_build, annotation_version)",
        ]
        
        for query in index_queries:
            self.connection.execute(query)
        
        self.connection.commit()
    
    def insert_mapping(self, mapping: GeneMapping) -> bool:
        """
        Insert a gene mapping, skipping duplicates.
        
        We use INSERT OR REPLACE so if the same gene appears in multiple
        data sources, the last one inserted wins.
        """
        try:
            self.connection.execute("""
                INSERT OR REPLACE INTO gene_mappings (
                    ensembl_id, gene_symbol, entrez_id, uniprot_id,
                    chromosome, start, end, strand,
                    biotype, description, aliases_json,
                    genome_build, annotation_version, data_source,
                    is_current, is_primary, content_hash
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                mapping.identifiers.ensembl_id,
                mapping.identifiers.gene_symbol, 
                mapping.identifiers.entrez_id,
                mapping.identifiers.uniprot_id,
                mapping.coordinates.chromosome if mapping.coordinates else None,
                mapping.coordinates.start if mapping.coordinates else None,
                mapping.coordinates.end if mapping.coordinates else None, 
                mapping.coordinates.strand if mapping.coordinates else None,
                mapping.biotype,
                mapping.description,
                json.dumps(mapping.aliases),  # Store list as JSON string
                mapping.genome_build,
                mapping.annotation_version,
                mapping.data_source,
                mapping.is_current,
                mapping.is_primary,
                mapping.content_hash
            ))
            
            self.connection.commit()
            return True
            
        except sqlite3.IntegrityError:
            # This usually means a duplicate content_hash - we can skip it
            logger.debug(f"Duplicate gene mapping: {mapping.identifiers.primary_id}")
            return False
        except Exception as e:
            logger.error(f"Failed to insert gene mapping: {e}")
            return False
    
    def close(self):
        """Clean up database connection."""
        if self.connection:
            self.connection.close()
    
    # Support 'with' statement for automatic cleanup
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def find_by_ensembl(self, ensembl_id: str, genome_build: str) -> List[GeneMapping]:
        """Find gene mappings by Ensembl ID and genome build."""
        if not genome_build:
            raise ValueError("genome_build is required for queries")
        
        cursor = self.connection.execute("""
            SELECT * FROM gene_mappings 
            WHERE ensembl_id = ? AND genome_build = ? AND is_current = TRUE
        """, (ensembl_id, genome_build))
        
        return [self._row_to_mapping(row) for row in cursor.fetchall()]

    def find_by_symbol(self, symbol: str, genome_build: str) -> List[GeneMapping]:
        """Find gene mappings by gene symbol and genome build."""
        if not genome_build:
            raise ValueError("genome_build is required for queries")
        
        cursor = self.connection.execute("""
            SELECT * FROM gene_mappings 
            WHERE gene_symbol = ? AND genome_build = ? AND is_current = TRUE
        """, (symbol, genome_build))
        
        return [self._row_to_mapping(row) for row in cursor.fetchall()]

    def find_by_entrez(self, entrez_id: str, genome_build: str) -> List[GeneMapping]:
        """Find gene mappings by Entrez ID and genome build."""
        if not genome_build:
            raise ValueError("genome_build is required for queries")
        
        cursor = self.connection.execute("""
            SELECT * FROM gene_mappings 
            WHERE entrez_id = ? AND genome_build = ? AND is_current = TRUE
        """, (entrez_id, genome_build))
        
        return [self._row_to_mapping(row) for row in cursor.fetchall()]

    def _row_to_mapping(self, row) -> GeneMapping:
        """Convert a database row to a GeneMapping object."""
        import json
        
        # Reconstruct identifiers
        identifiers = GeneIdentifier(
            ensembl_id=row['ensembl_id'],
            gene_symbol=row['gene_symbol'],
            entrez_id=row['entrez_id'],
            uniprot_id=row['uniprot_id']
        )
        
        # Reconstruct coordinates if available
        coordinates = None
        if row['chromosome'] and row['start'] and row['end']:
            coordinates = GenomicCoordinates(
                chromosome=row['chromosome'],
                start=row['start'],
                end=row['end'],
                strand=row['strand'],
                genome_build=row['genome_build']
            )
        
        # Reconstruct the full mapping
        mapping = GeneMapping(
            identifiers=identifiers,
            coordinates=coordinates,
            biotype=row['biotype'],
            description=row['description'],
            aliases=json.loads(row['aliases_json']) if row['aliases_json'] else [],
            genome_build=row['genome_build'],
            annotation_version=row['annotation_version'],
            data_source=row['data_source'],
            is_current=bool(row['is_current']),
            is_primary=bool(row['is_primary'])
        )
        
        return mapping
"""
Database versioning and update management.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


class DatabaseVersion:
    """Manages database version information and updates."""
    
    def __init__(self, db_path: Path):
        self.db_path = db_path
    
    def get_current_version(self) -> Optional[Dict[str, Any]]:
        """Get current database version info."""
        import sqlite3
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Check if metadata table exists
            cursor.execute("""
                SELECT name FROM sqlite_master 
                WHERE type='table' AND name='database_metadata'
            """)
            if not cursor.fetchone():
                return None
            
            # Get version info
            cursor.execute("""
                SELECT ensembl_release, species, genome_build, last_updated 
                FROM database_metadata 
                WHERE id = 1
            """)
            result = cursor.fetchone()
            conn.close()
            
            if result:
                return {
                    'ensembl_release': result[0],
                    'species': result[1],
                    'genome_build': result[2],
                    'last_updated': result[3]
                }
            return None
            
        except Exception as e:
            logger.error(f"Error getting database version: {e}")
            return None
    
    def set_version(self, ensembl_release: str, species: str, genome_build: str):
        """Set database version information."""
        import sqlite3
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Create metadata table if it doesn't exist
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS database_metadata (
                    id INTEGER PRIMARY KEY CHECK (id = 1),
                    ensembl_release TEXT NOT NULL,
                    species TEXT NOT NULL,
                    genome_build TEXT NOT NULL,
                    last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Insert or update version info
            cursor.execute("""
                INSERT OR REPLACE INTO database_metadata 
                (id, ensembl_release, species, genome_build, last_updated)
                VALUES (1, ?, ?, ?, ?)
            """, (ensembl_release, species, genome_build, datetime.now()))
            
            conn.commit()
            conn.close()
            logger.info(f"Database version set: {species} release {ensembl_release} ({genome_build})")
            
        except Exception as e:
            logger.error(f"Error setting database version: {e}")
            raise
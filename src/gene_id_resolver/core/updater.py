"""
Database update management for new Ensembl releases.
"""

import sqlite3
from pathlib import Path
from typing import Optional, Dict, Any
import logging
from datetime import datetime
import urllib.request
import re

logger = logging.getLogger(__name__)

class DatabaseUpdater:
    """Manages database updates and version migration."""
    
    ENSEMBL_FTP_BASE = "http://ftp.ensembl.org/pub"
    
    def __init__(self, db_path: Path):
        self.db_path = db_path
    
    def get_latest_ensembl_release(self) -> Optional[str]:
        """Fetch the latest Ensembl release number from FTP server."""
        try:
            # The current_README file contains the latest release info
            url = f"{self.ENSEMBL_FTP_BASE}/current_README"
            with urllib.request.urlopen(url, timeout=10) as response:
                content = response.read().decode('utf-8')
                # Look for "Ensembl Release XXX"
                match = re.search(r'Ensembl Release (\d+)', content)
                if match:
                    return match.group(1)
            
            # Fallback: check directory listing
            url = f"{self.ENSEMBL_FTP_BASE}/"
            with urllib.request.urlopen(url, timeout=10) as response:
                content = response.read().decode('utf-8')
                # Find all release-XXX directories
                releases = re.findall(r'release-(\d+)', content)
                if releases:
                    return str(max(int(r) for r in releases))
            
            return None
            
        except Exception as e:
            logger.warning(f"Could not fetch latest Ensembl release: {e}")
            return None
    
    def check_for_updates(self) -> Dict[str, Any]:
        """Check if database needs updates."""
        current_version = self._get_current_version()
        
        if not current_version:
            return {'status': 'no_version_info', 'message': 'No version information found'}
        
        latest_release = self.get_latest_ensembl_release()
        
        if not latest_release:
            return {
                'status': 'check_failed',
                'current': current_version,
                'message': 'Could not check for updates (network issue)'
            }
        
        if current_version['ensembl_release'] == latest_release:
            return {
                'status': 'current', 
                'current': current_version,
                'message': 'Database is up to date'
            }
        else:
            return {
                'status': 'update_available',
                'current': current_version,
                'latest_release': latest_release,
                'message': f"Update available: {current_version['ensembl_release']} → {latest_release}"
            }
    
    def _get_current_version(self) -> Optional[Dict[str, Any]]:
        """Get current database version."""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                SELECT ensembl_release, species, genome_build, last_updated 
                FROM database_metadata WHERE id = 1
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
            logger.error(f"Error getting current version: {e}")
            return None
    
    def incremental_update(self, target_release: str, species: str = "homo_sapiens", 
                          genome_build: str = None) -> bool:
        """Perform incremental update to target release."""
        from ..datasources.ensembl import EnsemblDownloader
        from .database import GeneDatabase
        
        logger.info(f"Updating database to Ensembl release {target_release}")
        
        try:
            # Get current version info
            current_version = self._get_current_version()
            if not current_version:
                logger.error("Cannot update: no version information found")
                return False
            
            # Check if we're already at target version
            if current_version['ensembl_release'] == target_release:
                logger.info(f"Already at release {target_release}")
                return True
            
            # For now, do a full rebuild rather than incremental
            # (incremental would be complex - need to diff genes between releases)
            logger.info("Performing full database rebuild with new release...")
            
            # Download new data
            downloads_dir = self.db_path.parent / "downloads"
            downloads_dir.mkdir(exist_ok=True)
            
            downloader = EnsemblDownloader(downloads_dir)
            gtf_path = downloader.download_gtf(species, target_release, genome_build)
            
            # Backup current database
            backup_path = self.db_path.with_suffix('.db.backup')
            import shutil
            shutil.copy2(self.db_path, backup_path)
            logger.info(f"Created backup: {backup_path}")
            
            # Rebuild database with new data
            # First, clear existing data
            conn = sqlite3.connect(self.db_path)
            conn.execute("DELETE FROM gene_mappings")
            conn.commit()
            conn.close()
            
            # Now load new data
            db = GeneDatabase(self.db_path)
            gene_count = 0
            genome_build_final = genome_build or current_version.get('genome_build', 'GRCh38')
            for gene_mapping in downloader.parse_gtf(gtf_path, genome_build_final, target_release):
                if db.insert_mapping(gene_mapping):
                    gene_count += 1
            
            # Update metadata
            db.update_metadata(
                ensembl_release=target_release,
                species=species,
                genome_build=genome_build_final
            )
            db.close()
            
            logger.info(f"Database updated successfully! {gene_count} genes loaded")
            return True
            
        except Exception as e:
            logger.error(f"Update failed: {e}")
            # Try to restore backup if it exists
            backup_path = self.db_path.with_suffix('.db.backup')
            if backup_path.exists():
                import shutil
                shutil.copy2(backup_path, self.db_path)
                logger.info("Restored database from backup")
            return False
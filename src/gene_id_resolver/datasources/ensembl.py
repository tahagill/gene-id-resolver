"""
Download and parse gene data from Ensembl.

Ensembl is our primary data source because it has the most consistent
gene annotations and versioning. We'll download their GTF files which
contain both gene coordinates and identifiers.
"""

import gzip
import logging
from pathlib import Path
from typing import List, Iterator
import urllib.request
from urllib.error import URLError

import pandas as pd
from tqdm import tqdm

from ..core.models import GeneIdentifier, GenomicCoordinates, GeneMapping

logger = logging.getLogger(__name__)


class EnsemblDownloader:
    """
    Handles downloading and parsing Ensembl gene data.
    
    We use their GTF files because they contain both gene annotations
    and cross-references to other databases like Entrez and UniProt.
    """
    
    ENSEMBL_FTP_BASE = "http://ftp.ensembl.org/pub"
    
    def __init__(self, download_dir: Path):
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)
    
    def get_download_url(self, species: str = "homo_sapiens", 
                        release: str = "109", 
                        genome_build: str = None) -> str:
        """
        Construct the download URL for Ensembl GTF file.
        
        Supports multiple species and genome builds.
        """
        # Default genome builds for common species
        species_config = {
            "homo_sapiens": {"name": "Homo_sapiens", "build": "GRCh38"},
            "mus_musculus": {"name": "Mus_musculus", "build": "GRCm39"},  # mm10 equivalent
            "rattus_norvegicus": {"name": "Rattus_norvegicus", "build": "mRatBN7.2"},
        }
        
        if species not in species_config:
            raise ValueError(f"Unsupported species: {species}. Supported: {list(species_config.keys())}")
        
        config = species_config[species]
        if genome_build:
            config["build"] = genome_build
            
        filename = f"{config['name']}.{config['build']}.{release}.gtf.gz"
        url = (f"{self.ENSEMBL_FTP_BASE}/release-{release}/gtf/"
               f"{species}/{filename}")
        return url
    
    def download_gtf(self, species: str = "homo_sapiens", 
                    release: str = "109",
                    genome_build: str = None) -> Path:
        """
        Download Ensembl GTF file with progress tracking.
        
        Returns the path to the downloaded file. Uses tqdm to show
        download progress since these files can be 50-100MB.
        """
        url = self.get_download_url(species, release, genome_build)
        local_path = self.download_dir / f"ensembl_{species}_{release}.gtf.gz"
        
        if local_path.exists():
            logger.info(f"Using cached file: {local_path}")
            return local_path
        
        logger.info(f"Downloading Ensembl release {release} for {species}...")
        logger.info(f"URL: {url}")
        
        try:
            # Use urllib with progress tracking
            def update_progress(block_num, block_size, total_size):
                if total_size > 0:
                    percent = min(100, (block_num * block_size * 100) // total_size)
                    print(f"\rDownload progress: {percent}%", end='', flush=True)
            
            urllib.request.urlretrieve(url, local_path, update_progress)
            print()  # New line after progress bar
            
            logger.info(f"Download completed: {local_path}")
            return local_path
            
        except URLError as e:
            logger.error(f"Failed to download from {url}: {e}")
            raise
    
    def parse_gtf(self, gtf_path: Path, genome_build: str = None, 
                  release: str = None) -> Iterator[GeneMapping]:
        """
        Parse Ensembl GTF file and yield GeneMapping objects.
        Handles both gzipped and plain text files.
        """
        logger.info(f"Parsing GTF file: {gtf_path}")
        
        # Extract metadata from filename if not provided
        # Filename format: Homo_sapiens.GRCh38.109.gtf.gz or ensembl_homo_sapiens_109.gtf.gz
        if not genome_build or not release:
            filename = gtf_path.stem.replace('.gtf', '')
            parts = filename.split('.')
            if len(parts) >= 3 and not genome_build:
                genome_build = parts[1]  # e.g., GRCh38
            if len(parts) >= 3 and not release:
                release = parts[2]  # e.g., 109
            # Fallback to defaults if parsing failed
            if not genome_build:
                genome_build = "GRCh38"
            if not release:
                release = "109"
        
        # Detect if file is gzipped or plain text
        is_gzipped = False
        with open(gtf_path, 'rb') as f:
            magic_number = f.read(2)
            is_gzipped = (magic_number == b'\x1f\x8b')  # Gzip magic number
        
        # Use pandas for efficient parsing
        column_names = ['seqname', 'source', 'feature', 'start', 'end', 
                    'score', 'strand', 'frame', 'attributes']
        
        try:
            if is_gzipped:
                # Read compressed GTF
                df = pd.read_csv(gtf_path, compression='gzip', 
                            sep='\t', comment='#', header=None,
                            names=column_names, low_memory=False)
            else:
                # Read plain text GTF
                df = pd.read_csv(gtf_path, sep='\t', comment='#', header=None,
                            names=column_names, low_memory=False)
        
        except Exception as e:
            logger.error(f"Failed to read GTF file {gtf_path}: {e}")
            raise
        
        # Filter only gene entries - other features like transcripts 
        # and exons would create duplicates
        gene_df = df[df['feature'] == 'gene']
        total_genes = len(gene_df)
        
        logger.info(f"Found {total_genes} gene entries to process")
        
        # Parse each gene entry with progress tracking
        for _, row in tqdm(gene_df.iterrows(), total=total_genes, 
                        desc="Parsing genes"):
            try:
                gene_mapping = self._parse_gene_row(row, genome_build, release)
                if gene_mapping:
                    yield gene_mapping
            except Exception as e:
                # Log but continue parsing - don't fail on single bad entry
                logger.debug(f"Failed to parse gene row: {e}")
                continue
    
    def _parse_gene_row(self, row, genome_build: str = "GRCh38", 
                       release: str = "109") -> GeneMapping:
        """
        Parse a single GTF row into a GeneMapping object.
        """
        # Parse the attribute column
        attributes = self._parse_attributes(row['attributes'])
        
        # Extract core identifiers
        gene_id = attributes.get('gene_id', '').split('.')[0]
        gene_name = attributes.get('gene_name')
        entrez_id = attributes.get('entrezgene_id')
        
        # Skip if no basic identifiers
        if not gene_id and not gene_name:
            return None
        
        # Create genomic coordinates
        coordinates = GenomicCoordinates(
            chromosome=str(row['seqname']),
            start=int(row['start']),
            end=int(row['end']),
            strand=row['strand'],
            genome_build=genome_build
        )
        
        # Create gene identifiers
        identifiers = GeneIdentifier(
            ensembl_id=gene_id,
            gene_symbol=gene_name,
            entrez_id=entrez_id
        )
        
        # Create the complete gene mapping
        mapping = GeneMapping(
            identifiers=identifiers,
            coordinates=coordinates,
            biotype=attributes.get('gene_biotype'),
            description=attributes.get('description', ''),
            genome_build=genome_build,
            annotation_version=release,
            data_source="ensembl"
        )
        
        return mapping
    
    def _parse_attributes(self, attribute_string: str) -> dict:
        """
        Parse GTF attribute column into a dictionary.
        
        The attribute format is: key1 "value1"; key2 "value2";
        We need to handle the quotes and trailing semicolons.
        """
        attributes = {}
        
        # Split by semicolon and clean up each pair
        pairs = attribute_string.strip().split(';')
        for pair in pairs:
            pair = pair.strip()
            if not pair:
                continue
            
            # Split into key and value
            if ' ' in pair:
                key, value = pair.split(' ', 1)
                # Remove quotes from value
                value = value.strip('"')
                attributes[key] = value
        
        return attributes

    def download_gene_history(self, species: str = "homo_sapiens", release: str = "109") -> Path:
        """
        Download gene symbol history file for deprecated gene tracking.
        
        For human genes, we use HGNC's comprehensive gene symbol history.
        For mouse genes, we use MGI's marker list with withdrawn symbols.
        For rat genes, we use RGD's obsolete gene IDs.
        """
        local_path = self.download_dir / f"gene_history_{species}_{release}.txt.gz"
        
        if local_path.exists():
            logger.info(f"Using cached gene history: {local_path}")
            return local_path
        
        if species == "homo_sapiens":
            # HGNC provides comprehensive gene symbol history via their REST API
            url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
            logger.info(f"Downloading HGNC gene symbol history...")
            
        elif species == "mus_musculus":
            # MGI provides marker list including withdrawn symbols
            url = "http://www.informatics.jax.org/downloads/reports/MRK_List1.rpt"
            logger.info(f"Downloading MGI marker list with withdrawn symbols...")
            
        elif species == "rattus_norvegicus":
            # RGD provides obsolete gene IDs with mappings
            url = "https://download.rgd.mcw.edu/data_release/GENES_OBSOLETE_IDS.txt"
            logger.info(f"Downloading RGD obsolete gene IDs...")
            
        else:
            logger.info(f"Gene history not available for {species}, using fallback mappings")
            return None
        
        try:
            import gzip
            import urllib.request
            
            logger.info(f"URL: {url}")
            
            # Download uncompressed, then compress locally
            temp_path = self.download_dir / f"gene_history_{species}_{release}.txt"
            
            def update_progress(block_num, block_size, total_size):
                if total_size > 0:
                    percent = min(100, (block_num * block_size * 100) // total_size)
                    print(f"\rDownload progress: {percent}%", end='', flush=True)
            
            urllib.request.urlretrieve(url, temp_path, update_progress)
            print()  # New line after progress bar
            
            # Compress it
            with open(temp_path, 'rb') as f_in:
                with gzip.open(local_path, 'wb') as f_out:
                    f_out.writelines(f_in)
            
            # Remove temp file
            temp_path.unlink()
            
            logger.info(f"Gene history download completed: {local_path}")
            return local_path
            
        except Exception as e:
            logger.warning(f"Failed to download gene history: {e}")
            logger.warning("Deprecated gene detection will use fallback mappings only")
            return None
"""
Core data structures for gene ID resolution.

This is where we define how gene information is stored and handled.
Using dataclasses here makes the code cleaner and less error-prone.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any
from datetime import datetime
import hashlib


@dataclass(frozen=True)
class GeneIdentifier:
    """
    Represents a gene across different naming systems.
    
    Frozen so we can use these as dictionary keys safely.
    Having this as a separate class prevents ID mixing bugs.
    """
    ensembl_id: Optional[str] = None
    gene_symbol: Optional[str] = None
    entrez_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    
    def __post_init__(self):
        # At least one ID must be provided - no empty identifiers
        if not any([self.ensembl_id, self.gene_symbol, self.entrez_id, self.uniprot_id]):
            raise ValueError("Gene must have at least one identifier")
    
    @property
    def primary_id(self) -> str:
        """Get the most reliable ID, preferring Ensembl as it's most stable."""
        return self.ensembl_id or self.gene_symbol or self.entrez_id or self.uniprot_id
    
    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'ensembl_id': self.ensembl_id,
            'gene_symbol': self.gene_symbol,
            'entrez_id': self.entrez_id,
            'uniprot_id': self.uniprot_id
        }

@dataclass
class GenomicCoordinates:
    """
    Standardized genomic location.
    
    Using 1-based coordinates because that's what Ensembl and most 
    bioinformatics tools use. This avoids off-by-one errors.
    """
    chromosome: str
    start: int  # 1-based position
    end: int    
    strand: str  # '+' or '-'
    genome_build: str  # hg38, hg19, mm10, etc.
    
    def __post_init__(self):
        # Validate we have sane genomic coordinates
        if self.start < 1:
            raise ValueError(f"Start position can't be less than 1: {self.start}")
        if self.end < self.start:
            raise ValueError(f"End position {self.end} can't be before start {self.start}")
        
        # Convert chromosome to string and standardize naming
        # Some files use numbers (1), others use strings (chr1)
        chromosome_str = str(self.chromosome)
        if not chromosome_str.startswith('chr'):
            object.__setattr__(self, 'chromosome', f'chr{chromosome_str}')


@dataclass
class GeneMapping:
    """
    Complete gene information from a specific data source.
    
    This is our main data structure that gets stored in the database.
    The content_hash lets us detect duplicates across data sources.
    """
    identifiers: GeneIdentifier
    coordinates: Optional[GenomicCoordinates] = None
    
    # Gene metadata
    biotype: Optional[str] = None  # protein_coding, lncRNA, etc.
    description: Optional[str] = None
    aliases: List[str] = field(default_factory=list)  # Other symbols this gene is known by
    
    # Version info - critical for reproducibility
    genome_build: str = "hg38"
    annotation_version: str = "109"  # Ensembl release version
    data_source: str = "ensembl"  # Where this mapping came from
    import_date: datetime = field(default_factory=datetime.now)
    
    # For handling genes that map to multiple locations
    is_current: bool = True  # Not retired/deprecated
    is_primary: bool = True  # Main mapping when multiple exist
    
    def __post_init__(self):
        # Create a unique fingerprint of this gene version
        # This helps us avoid duplicates when loading from multiple sources
        content_str = f"{self.identifiers.primary_id}|{self.genome_build}|{self.annotation_version}"
        content_hash = hashlib.md5(content_str.encode()).hexdigest()
        object.__setattr__(self, 'content_hash', content_hash)


@dataclass
class ConversionResult:
    """
    Tracks what happened during ID conversion.
    
    Most tools just return a dict and you don't know why some IDs failed.
    This gives complete transparency about the conversion process.
    """
    successful: Dict[str, GeneMapping] = field(default_factory=dict)
    ambiguous: Dict[str, List[GeneMapping]] = field(default_factory=dict)  # 1-to-many mappings
    failed: List[str] = field(default_factory=list)  # IDs we couldn't convert
    
    # Audit information
    input_type: Optional[str] = None
    output_type: Optional[str] = None
    resolver_config: Dict[str, Any] = field(default_factory=dict)
    conversion_timestamp: datetime = field(default_factory=datetime.now)
    
    # Track ambiguity resolution
    ambiguity_resolutions: Dict[str, str] = field(default_factory=dict)

    def add_success(self, input_id: str, mapping: GeneMapping):
        """Record a successful conversion."""
        self.successful[input_id] = mapping
    
    def add_ambiguous(self, input_id: str, mappings: List[GeneMapping]):
        """Record when one ID maps to multiple genes."""
        self.ambiguous[input_id] = mappings
    
    def add_failed(self, input_id: str):
        """Record IDs that couldn't be converted."""
        self.failed.append(input_id)
    
    @property
    def success_rate(self) -> float:
        """Quick way to see how well the conversion went."""
        total = len(self.successful) + len(self.ambiguous) + len(self.failed)
        return len(self.successful) / total if total > 0 else 0.0
    
    def generate_audit_report(self) -> Dict[str, Any]:
        """Generate a comprehensive audit report."""
        return {
            'summary': {
                'total_inputs': (len(self.successful) + 
                               len(self.ambiguous) + 
                               len(self.failed)),
                'successful': len(self.successful),
                'ambiguous': len(self.ambiguous),
                'failed': len(self.failed),
                'success_rate': self.success_rate,
                'timestamp': self.conversion_timestamp.isoformat(),
                'ambiguity_resolutions': self.ambiguity_resolutions
            },
            'resolver_config': self.resolver_config,
            'ambiguous_details': {
                input_id: [m.identifiers.to_dict() for m in mappings]
                for input_id, mappings in self.ambiguous.items()
            },
            'failed_details': self.failed
        }
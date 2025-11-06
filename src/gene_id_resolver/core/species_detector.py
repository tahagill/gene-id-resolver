"""
Detect and handle species-specific gene symbols.
"""

from typing import Dict, List, Optional
from .models import GeneMapping
import logging

logger = logging.getLogger(__name__)

class SpeciesDetector:
    """Detect and handle species-specific gene conflicts."""
    
    # Common genes that exist in multiple species
    MULTI_SPECIES_GENES = {
        "MEF2C", "TP53", "BRCA1", "BRCA2", "MYC", "ACTB", "GAPDH",
        "TUBB", "HSP90", "EGFR", "KRAS", "PTEN", "AKT1", "MAPK1"
    }
    
    def __init__(self):
        self.species_indicators = {
            "human": {"GRCh38", "GRCh37", "hg38", "hg19"},
            "mouse": {"GRCm39", "GRCm38", "mm10", "mm9"},
            "rat": {"mRatBN7.2", "rn6", "rn5"}
        }
    
    def detect_species_conflict(self, mappings: List[GeneMapping]) -> Optional[str]:
        """Detect if mappings span multiple species."""
        if not mappings or len(mappings) < 2:
            return None
        
        genome_builds = set()
        for mapping in mappings:
            if mapping.coordinates:
                genome_builds.add(mapping.genome_build)
        
        if len(genome_builds) > 1:
            return f"Multiple genome builds detected: {', '.join(genome_builds)}"
        
        return None
    
    def is_multi_species_gene(self, symbol: str) -> bool:
        """Check if a gene symbol commonly exists in multiple species."""
        return symbol.upper() in self.MULTI_SPECIES_GENES
    
    def filter_by_species(self, mappings: List[GeneMapping], target_species: str) -> List[GeneMapping]:
        """Filter mappings to target species."""
        target_indicators = self.species_indicators.get(target_species.lower(), set())
        
        filtered = []
        for mapping in mappings:
            if any(indicator in mapping.genome_build for indicator in target_indicators):
                filtered.append(mapping)
        
        if filtered and len(filtered) < len(mappings):
            logger.info(f"Filtered {len(mappings) - len(filtered)} non-{target_species} mappings")
        
        return filtered
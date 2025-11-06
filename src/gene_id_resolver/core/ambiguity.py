"""
Enhanced ambiguity resolution with automatic conflict detection.
"""

from enum import Enum
from typing import List, Dict, Tuple
from .models import GeneMapping
import logging

logger = logging.getLogger(__name__)

class AmbiguityStrategy(Enum):
    FAIL = "fail"
    FIRST = "first"  
    PRIMARY = "primary"
    ALL = "all"
    AUTO = "auto"  # NEW: Smart automatic resolution

def resolve_ambiguity(mappings: List[GeneMapping], strategy: AmbiguityStrategy) -> List[GeneMapping]:
    """Enhanced ambiguity resolution with conflict detection."""
    if not mappings:
        return []
    
    # NEW: Detect and log ambiguity type
    ambiguity_type = _classify_ambiguity(mappings)
    logger.debug(f"Ambiguity type for {mappings[0].identifiers.gene_symbol}: {ambiguity_type}")
    
    if strategy == AmbiguityStrategy.FAIL:
        return []
    elif strategy == AmbiguityStrategy.FIRST:
        return [mappings[0]]
    elif strategy == AmbiguityStrategy.PRIMARY:
        return _resolve_primary(mappings)
    elif strategy == AmbiguityStrategy.ALL:
        return mappings
    elif strategy == AmbiguityStrategy.AUTO:
        return _resolve_auto(mappings, ambiguity_type)
    else:
        raise ValueError(f"Unknown ambiguity strategy: {strategy}")

def _classify_ambiguity(mappings: List[GeneMapping]) -> str:
    """Classify the type of ambiguity for better resolution."""
    if len(mappings) == 1:
        return "unique"
    
    # Check if it's biotype ambiguity (like SEPT4 vs SEPTIN4)
    biotypes = set(m.biotype for m in mappings if m.biotype)
    if len(biotypes) > 1:
        return "biotype_conflict"
    
    # Check if it's genomic location ambiguity
    locations = set((m.coordinates.chromosome, m.coordinates.start) for m in mappings if m.coordinates)
    if len(locations) > 1:
        return "location_conflict"
    
    # Check if it's species ambiguity (same symbol, different species)
    genome_builds = set(m.genome_build for m in mappings)
    if len(genome_builds) > 1:
        return "species_conflict"
    
    return "unknown_conflict"

def _resolve_auto(mappings: List[GeneMapping], ambiguity_type: str) -> List[GeneMapping]:
    """Smart automatic resolution based on ambiguity type."""
    if ambiguity_type == "biotype_conflict":
        # Prefer protein-coding genes
        protein_coding = [m for m in mappings if m.biotype == "protein_coding"]
        if protein_coding:
            logger.info(f"Auto-resolved biotype conflict: selected protein-coding gene")
            return [protein_coding[0]]
    
    elif ambiguity_type == "location_conflict":
        # For same gene in different locations, prefer canonical location
        # This is complex - for now, fall back to primary strategy
        return _resolve_primary(mappings)
    
    elif ambiguity_type == "species_conflict":
        # Prefer human genes if available
        human_genes = [m for m in mappings if "GRCh" in m.genome_build]
        if human_genes:
            logger.info(f"Auto-resolved species conflict: selected human gene")
            return [human_genes[0]]
    
    # Fall back to primary strategy
    return _resolve_primary(mappings)

def _resolve_primary(mappings: List[GeneMapping]) -> List[GeneMapping]:
    """Enhanced primary resolution with better logging."""
    protein_coding = [m for m in mappings if m.biotype == "protein_coding"]
    
    if protein_coding:
        if len(protein_coding) > 1:
            logger.warning(f"Multiple protein-coding genes found, selecting first")
        return [protein_coding[0]]
    else:
        logger.warning(f"No protein-coding genes found, falling back to first match")
        return [mappings[0]]
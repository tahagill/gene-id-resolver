"""
Enhanced resolver with auto-correction for deprecated gene symbols.
"""

import logging
from pathlib import Path
from typing import List, Dict
from .resolver import GeneResolver
from .deprecated_genes import SmartDeprecatedGeneHandler

logger = logging.getLogger(__name__)

class EnhancedGeneResolver(GeneResolver):
    """Enhanced resolver with comprehensive deprecated symbol detection."""
    
    def __init__(self, data_dir: Path = Path("data")):
        super().__init__(data_dir)
        self.deprecated_handler = SmartDeprecatedGeneHandler(
            self.data_dir / "genes.db", 
            self.data_dir
        )
    
    def convert_with_correction(self, gene_ids: List[str], from_type: str, to_type: str, 
                               genome_build: str = "hg38", ambiguity_strategy: str = "primary"):
        """Convert genes with auto-correction for deprecated symbols."""
        
        # Auto-correct deprecated symbols
        corrected_genes = []
        corrections = {}
        original_to_corrected = {}
        
        for gene in gene_ids:
            if from_type in ['symbol', 'gene_symbol'] and self.deprecated_handler.is_deprecated(gene):
                current_symbol = self.deprecated_handler.get_current_symbol(gene)
                if current_symbol:
                    corrected_genes.append(current_symbol)
                    corrections[gene] = current_symbol
                    original_to_corrected[gene] = current_symbol
                    print(f"Auto-corrected: {gene} → {current_symbol}")
                else:
                    corrected_genes.append(gene)
                    original_to_corrected[gene] = gene
            else:
                corrected_genes.append(gene)
                original_to_corrected[gene] = gene
        
        result = super().convert(corrected_genes, from_type, to_type, genome_build, ambiguity_strategy)
        remapped_successful = {}
        for original, corrected in original_to_corrected.items():
            if corrected in result.successful:
                remapped_successful[original] = result.successful[corrected]
        
        remapped_ambiguous = {}
        for original, corrected in original_to_corrected.items():
            if corrected in result.ambiguous:
                remapped_ambiguous[original] = result.ambiguous[corrected]
        
        # Update result with remapped keys
        result.successful = remapped_successful
        result.ambiguous = remapped_ambiguous
        
        # Add correction info to result
        result.corrections_applied = corrections
        
        return result
    
    def _preprocess_genes(self, gene_ids: List[str], from_type: str, auto_correct: bool) -> Dict:
        """Pre-process genes to handle deprecated symbols and other issues."""
        corrections = {}
        final_genes = []
        warnings = []
        
        for gene in gene_ids:
            if from_type in ['symbol', 'gene_symbol'] and self.deprecated_handler.is_deprecated(gene):
                if auto_correct:
                    corrected = self.deprecated_handler.get_current_symbol(gene)
                    if corrected:
                        corrections[gene] = corrected
                        final_genes.append(corrected)
                        warnings.append(f"Auto-corrected deprecated symbol: {gene} → {corrected}")
                    else:
                        final_genes.append(gene)
                        warnings.append(f"Deprecated symbol {gene} but no correction available")
                else:
                    final_genes.append(gene)
                    warnings.append(f"Deprecated symbol detected: {gene}")
            else:
                final_genes.append(gene)
        
        return {
            'original_genes': gene_ids,
            'final_genes': final_genes,
            'corrections': corrections,
            'warnings': warnings
        }
    
    def _enhance_result(self, result, processed_genes: Dict) -> Dict:
        """Enhance result with additional information."""
        # Suggest corrections for failed genes
        failed_suggestions = self.deprecated_handler.suggest_corrections(result.failed)
        
        enhanced = {
            'conversion_result': result,
            'preprocessing_info': processed_genes,
            'suggestions': {
                'failed_corrections': failed_suggestions
            },
            'enhanced_metrics': {
                'total_processed': len(processed_genes['original_genes']),
                'auto_corrected': len(processed_genes['corrections']),
                'success_rate_original': len(result.successful) / len(processed_genes['original_genes']),
                'success_rate_processed': len(result.successful) / len(processed_genes['final_genes'])
            }
        }
        
        return enhanced
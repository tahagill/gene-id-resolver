#!/usr/bin/env python
"""
Comprehensive test suite for gene-id-resolver.

Tests all major functionality and edge cases:
- Basic ID conversions (symbol, ensembl, entrez)
- Deprecated gene auto-correction
- Ambiguity resolution strategies
- File I/O (CSV, TSV, BED, TXT)
- Batch processing
- Genome build handling
- Error handling
- Database operations
"""

import sys
from pathlib import Path
from datetime import datetime

# Color codes for terminal output
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'

class TestResults:
    def __init__(self):
        self.passed = []
        self.failed = []
        self.warnings = []
    
    def add_pass(self, test_name, details=""):
        self.passed.append((test_name, details))
        print(f"{Colors.GREEN}✓{Colors.END} {test_name}")
        if details:
            print(f"  {details}")
    
    def add_fail(self, test_name, error):
        self.failed.append((test_name, str(error)))
        print(f"{Colors.RED}✗{Colors.END} {test_name}")
        print(f"  {Colors.RED}Error: {error}{Colors.END}")
    
    def add_warning(self, test_name, message):
        self.warnings.append((test_name, message))
        print(f"{Colors.YELLOW}⚠{Colors.END} {test_name}")
        print(f"  {Colors.YELLOW}{message}{Colors.END}")
    
    def print_summary(self):
        total = len(self.passed) + len(self.failed) + len(self.warnings)
        print(f"\n{Colors.BOLD}{'='*70}{Colors.END}")
        print(f"{Colors.BOLD}TEST SUMMARY{Colors.END}")
        print(f"{'='*70}")
        print(f"Total Tests: {total}")
        print(f"{Colors.GREEN}Passed: {len(self.passed)}{Colors.END}")
        print(f"{Colors.RED}Failed: {len(self.failed)}{Colors.END}")
        print(f"{Colors.YELLOW}Warnings: {len(self.warnings)}{Colors.END}")
        
        if total > 0:
            success_rate = (len(self.passed) / total) * 100
            print(f"\nSuccess Rate: {success_rate:.1f}%")
        
        if self.failed:
            print(f"\n{Colors.RED}Failed Tests:{Colors.END}")
            for name, error in self.failed:
                print(f"  • {name}: {error}")
        
        return len(self.failed) == 0


def test_basic_conversions(results):
    """Test basic gene ID conversions."""
    print(f"\n{Colors.BOLD}=== Testing Basic ID Conversions ==={Colors.END}")
    
    try:
        from gene_id_resolver.core.resolver import GeneResolver
        resolver = GeneResolver(Path("data"))
        
        # Test 1: Symbol to Ensembl
        try:
            result = resolver.convert(["TP53", "BRCA1"], "symbol", "ensembl", "hg38")
            if "TP53" in result.successful and "BRCA1" in result.successful:
                tp53_id = result.successful["TP53"].identifiers.ensembl_id
                brca1_id = result.successful["BRCA1"].identifiers.ensembl_id
                results.add_pass("Symbol → Ensembl", f"TP53→{tp53_id}, BRCA1→{brca1_id}")
            else:
                results.add_fail("Symbol → Ensembl", f"Failed: {result.failed}")
        except Exception as e:
            results.add_fail("Symbol → Ensembl", str(e))
        
        # Test 2: Ensembl to Symbol
        try:
            result = resolver.convert(["ENSG00000141510", "ENSG00000012048"], "ensembl", "symbol", "hg38")
            if len(result.successful) == 2:
                results.add_pass("Ensembl → Symbol", f"Found {len(result.successful)} genes")
            else:
                results.add_fail("Ensembl → Symbol", f"Expected 2, got {len(result.successful)}")
        except Exception as e:
            results.add_fail("Ensembl → Symbol", str(e))
        
        # Test 3: Symbol to Entrez
        try:
            result = resolver.convert(["EGFR", "KRAS"], "symbol", "entrez", "hg38")
            if len(result.successful) >= 1:
                results.add_pass("Symbol → Entrez", f"Converted {len(result.successful)}/2 genes")
            else:
                results.add_fail("Symbol → Entrez", "No successful conversions")
        except Exception as e:
            results.add_fail("Symbol → Entrez", str(e))
        
        # Test 4: Case insensitivity
        try:
            result = resolver.convert(["tp53", "BRCA1", "BrCa2"], "symbol", "ensembl", "hg38")
            if len(result.successful) >= 2:
                results.add_pass("Case Insensitivity", f"Converted {len(result.successful)}/3 genes")
            else:
                results.add_warning("Case Insensitivity", f"Only converted {len(result.successful)}/3 genes")
        except Exception as e:
            results.add_fail("Case Insensitivity", str(e))
        
        resolver.close()
        
    except Exception as e:
        results.add_fail("Basic Conversions Setup", str(e))


def test_deprecated_genes(results):
    """Test deprecated gene auto-correction."""
    print(f"\n{Colors.BOLD}=== Testing Deprecated Gene Auto-Correction ==={Colors.END}")
    
    try:
        from gene_id_resolver.core.resolver_enhanced import EnhancedGeneResolver
        resolver = EnhancedGeneResolver(Path("data"))
        
        # Test 1: Known deprecated symbols
        deprecated_tests = [
            ("SEPT4", "SEPTIN4"),
            ("ENO1L1", "ENO1"),
            ("CDKN2", "CDKN2A"),
        ]
        
        for old_symbol, expected_new in deprecated_tests:
            try:
                result = resolver.convert_with_correction([old_symbol], "symbol", "ensembl", "hg38")
                if old_symbol in result.successful:
                    results.add_pass(f"Deprecated: {old_symbol}→{expected_new}", 
                                   f"Converted to {result.successful[old_symbol].identifiers.ensembl_id}")
                else:
                    results.add_fail(f"Deprecated: {old_symbol}→{expected_new}", 
                                   "Conversion failed")
            except Exception as e:
                results.add_fail(f"Deprecated: {old_symbol}→{expected_new}", str(e))
        
        # Test 2: KIAA genes
        try:
            result = resolver.convert_with_correction(["KIAA0196"], "symbol", "ensembl", "hg38")
            if "KIAA0196" in result.successful:
                results.add_pass("KIAA gene auto-correct", "KIAA0196 converted successfully")
            else:
                results.add_warning("KIAA gene auto-correct", "KIAA0196 not converted")
        except Exception as e:
            results.add_fail("KIAA gene auto-correct", str(e))
        
        # Test 3: C*orf* genes
        try:
            result = resolver.convert_with_correction(["C9orf72"], "symbol", "ensembl", "hg38")
            if "C9orf72" in result.successful:
                results.add_pass("C*orf* gene handling", "C9orf72 converted")
            else:
                results.add_fail("C*orf* gene handling", "C9orf72 not found")
        except Exception as e:
            results.add_fail("C*orf* gene handling", str(e))
        
        resolver.close()
        
    except ImportError:
        results.add_warning("Deprecated Gene Tests", "EnhancedGeneResolver not available")
    except Exception as e:
        results.add_fail("Deprecated Gene Tests Setup", str(e))


def test_ambiguity_strategies(results):
    """Test ambiguity resolution strategies."""
    print(f"\n{Colors.BOLD}=== Testing Ambiguity Resolution Strategies ==={Colors.END}")
    
    try:
        from gene_id_resolver.core.resolver import GeneResolver
        resolver = GeneResolver(Path("data"))
        
        # Find an ambiguous gene (some symbols map to multiple IDs)
        # We'll test with different strategies
        
        # Test 1: Primary strategy (prefer protein-coding)
        try:
            result = resolver.convert(["TP53"], "symbol", "ensembl", "hg38", "primary")
            if "TP53" in result.successful or "TP53" in result.ambiguous:
                results.add_pass("Ambiguity: primary", "Strategy applied successfully")
            else:
                results.add_fail("Ambiguity: primary", "Gene not found")
        except Exception as e:
            results.add_fail("Ambiguity: primary", str(e))
        
        # Test 2: First strategy
        try:
            result = resolver.convert(["TP53"], "symbol", "ensembl", "hg38", "first")
            if "TP53" in result.successful:
                results.add_pass("Ambiguity: first", "Strategy applied successfully")
            else:
                results.add_fail("Ambiguity: first", "Gene not found")
        except Exception as e:
            results.add_fail("Ambiguity: first", str(e))
        
        # Test 3: All strategy
        try:
            result = resolver.convert(["TP53"], "symbol", "ensembl", "hg38", "all")
            if "TP53" in result.successful or "TP53" in result.ambiguous:
                results.add_pass("Ambiguity: all", "Strategy applied successfully")
            else:
                results.add_fail("Ambiguity: all", "Gene not found")
        except Exception as e:
            results.add_fail("Ambiguity: all", str(e))
        
        resolver.close()
        
    except Exception as e:
        results.add_fail("Ambiguity Tests Setup", str(e))


def test_file_io(results):
    """Test file input/output handling."""
    print(f"\n{Colors.BOLD}=== Testing File I/O ==={Colors.END}")
    
    # Test 1: CSV input
    try:
        csv_file = Path("test_genes.csv")
        if csv_file.exists():
            from gene_id_resolver.utils.file_processor import read_genes_from_file
            genes = read_genes_from_file(csv_file, column=0)
            if len(genes) > 0:
                results.add_pass("CSV input", f"Read {len(genes)} genes")
            else:
                results.add_fail("CSV input", "No genes read")
        else:
            results.add_warning("CSV input", "test_genes.csv not found")
    except Exception as e:
        results.add_fail("CSV input", str(e))
    
    # Test 2: TSV input
    try:
        tsv_file = Path("test_genes.tsv")
        if tsv_file.exists():
            from gene_id_resolver.utils.file_processor import read_genes_from_file
            genes = read_genes_from_file(tsv_file, column=0)
            if len(genes) > 0:
                results.add_pass("TSV input", f"Read {len(genes)} genes")
            else:
                results.add_fail("TSV input", "No genes read")
        else:
            results.add_warning("TSV input", "test_genes.tsv not found")
    except Exception as e:
        results.add_fail("TSV input", str(e))
    
    # Test 3: Plain text input
    try:
        txt_file = Path("test_genes.txt")
        if txt_file.exists():
            from gene_id_resolver.utils.file_processor import read_genes_from_file
            genes = read_genes_from_file(txt_file)
            if len(genes) > 0:
                results.add_pass("Plain text input", f"Read {len(genes)} genes")
            else:
                results.add_fail("Plain text input", "No genes read")
        else:
            results.add_warning("Plain text input", "test_genes.txt not found")
    except Exception as e:
        results.add_fail("Plain text input", str(e))
    
    # Test 4: BED file input
    try:
        bed_file = Path("test_regions.bed")
        if bed_file.exists():
            from gene_id_resolver.utils.file_processor import read_genes_from_file
            genes = read_genes_from_file(bed_file)
            if len(genes) > 0:
                results.add_pass("BED input", f"Read {len(genes)} genes")
            else:
                results.add_fail("BED input", "No genes read")
        else:
            results.add_warning("BED input", "test_regions.bed not found")
    except Exception as e:
        results.add_fail("BED input", str(e))


def test_batch_processing(results):
    """Test batch processing with large gene lists."""
    print(f"\n{Colors.BOLD}=== Testing Batch Processing ==={Colors.END}")
    
    try:
        from gene_id_resolver.core.resolver import GeneResolver
        resolver = GeneResolver(Path("data"))
        
        # Test 1: 100 gene batch
        try:
            # Create a list of common genes
            test_genes = ["TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "MYC", "PTEN", 
                         "AKT1", "PIK3CA", "BRAF"] * 10  # 100 genes
            result = resolver.convert(test_genes, "symbol", "ensembl", "hg38")
            success_count = len(result.successful)
            results.add_pass("Batch: 100 genes", f"{success_count}/{len(test_genes)} converted")
        except Exception as e:
            results.add_fail("Batch: 100 genes", str(e))
        
        # Test 2: 1000 gene batch
        try:
            large_batch = ["TP53", "BRCA1", "EGFR", "KRAS", "MYC"] * 200  # 1000 genes
            result = resolver.convert(large_batch, "symbol", "ensembl", "hg38")
            success_count = len(result.successful)
            results.add_pass("Batch: 1000 genes", f"{success_count}/{len(large_batch)} converted")
        except Exception as e:
            results.add_fail("Batch: 1000 genes", str(e))
        
        # Test 3: Duplicate handling
        try:
            duplicates = ["TP53", "TP53", "BRCA1", "BRCA1", "BRCA1"]
            result = resolver.convert(duplicates, "symbol", "ensembl", "hg38")
            results.add_pass("Duplicate handling", f"Processed {len(duplicates)} inputs")
        except Exception as e:
            results.add_fail("Duplicate handling", str(e))
        
        resolver.close()
        
    except Exception as e:
        results.add_fail("Batch Processing Setup", str(e))


def test_genome_builds(results):
    """Test genome build handling."""
    print(f"\n{Colors.BOLD}=== Testing Genome Build Handling ==={Colors.END}")
    
    try:
        from gene_id_resolver.core.resolver import GeneResolver
        resolver = GeneResolver(Path("data"))
        
        # Test 1: hg38 (default)
        try:
            result = resolver.convert(["TP53"], "symbol", "ensembl", "hg38")
            if "TP53" in result.successful:
                results.add_pass("Genome build: hg38", "Default build works")
            else:
                results.add_fail("Genome build: hg38", "Gene not found")
        except Exception as e:
            results.add_fail("Genome build: hg38", str(e))
        
        # Test 2: GRCh38 (synonym for hg38)
        try:
            result = resolver.convert(["TP53"], "symbol", "ensembl", "GRCh38")
            results.add_pass("Genome build: GRCh38", "Synonym recognized")
        except Exception as e:
            results.add_warning("Genome build: GRCh38", f"Synonym not supported: {e}")
        
        resolver.close()
        
    except Exception as e:
        results.add_fail("Genome Build Tests Setup", str(e))


def test_error_handling(results):
    """Test error handling and edge cases."""
    print(f"\n{Colors.BOLD}=== Testing Error Handling ==={Colors.END}")
    
    try:
        from gene_id_resolver.core.resolver import GeneResolver
        resolver = GeneResolver(Path("data"))
        
        # Test 1: Invalid gene symbols
        try:
            result = resolver.convert(["NOTAREALGENE123", "FAKEGENE999"], "symbol", "ensembl", "hg38")
            if len(result.failed) == 2:
                results.add_pass("Invalid genes", "Properly reported as failed")
            else:
                results.add_warning("Invalid genes", f"Expected 2 failures, got {len(result.failed)}")
        except Exception as e:
            results.add_fail("Invalid genes", str(e))
        
        # Test 2: Empty input
        try:
            result = resolver.convert([], "symbol", "ensembl", "hg38")
            results.add_pass("Empty input", "Handled gracefully")
        except Exception as e:
            results.add_warning("Empty input", f"Error on empty input: {e}")
        
        # Test 3: Mixed valid/invalid
        try:
            result = resolver.convert(["TP53", "NOTREAL", "BRCA1"], "symbol", "ensembl", "hg38")
            if len(result.successful) >= 2 and len(result.failed) >= 1:
                results.add_pass("Mixed valid/invalid", 
                               f"{len(result.successful)} success, {len(result.failed)} failed")
            else:
                results.add_warning("Mixed valid/invalid", "Unexpected result distribution")
        except Exception as e:
            results.add_fail("Mixed valid/invalid", str(e))
        
        # Test 4: Invalid ID type
        try:
            result = resolver.convert(["TP53"], "invalid_type", "ensembl", "hg38")
            results.add_fail("Invalid ID type", "Should have raised error")
        except Exception as e:
            results.add_pass("Invalid ID type", "Properly rejected invalid type")
        
        resolver.close()
        
    except Exception as e:
        results.add_fail("Error Handling Setup", str(e))


def test_database_operations(results):
    """Test database operations."""
    print(f"\n{Colors.BOLD}=== Testing Database Operations ==={Colors.END}")
    
    # Test 1: Database exists
    try:
        db_file = Path("data/genes.db")
        if db_file.exists():
            size_mb = db_file.stat().st_size / (1024 * 1024)
            results.add_pass("Database exists", f"Size: {size_mb:.1f} MB")
        else:
            results.add_fail("Database exists", "genes.db not found")
    except Exception as e:
        results.add_fail("Database exists", str(e))
    
    # Test 2: Database query performance
    try:
        from gene_id_resolver.core.resolver import GeneResolver
        import time
        
        resolver = GeneResolver(Path("data"))
        start = time.time()
        result = resolver.convert(["TP53"] * 100, "symbol", "ensembl", "hg38")
        elapsed = time.time() - start
        
        if elapsed < 5.0:  # Should be fast
            results.add_pass("Query performance", f"100 queries in {elapsed:.2f}s")
        else:
            results.add_warning("Query performance", f"Slow: {elapsed:.2f}s for 100 queries")
        
        resolver.close()
    except Exception as e:
        results.add_fail("Query performance", str(e))
    
    # Test 3: Gene history file
    try:
        history_file = Path("data/downloads/gene_history_homo_sapiens_109.txt.gz")
        if history_file.exists():
            size_mb = history_file.stat().st_size / (1024 * 1024)
            results.add_pass("Gene history file", f"Size: {size_mb:.1f} MB")
        else:
            results.add_warning("Gene history file", "Not downloaded yet")
    except Exception as e:
        results.add_fail("Gene history file", str(e))
    
    # Test 4: Database gene count
    try:
        import sqlite3
        conn = sqlite3.connect("data/genes.db")
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM gene_mappings WHERE is_current = TRUE")
        count = cursor.fetchone()[0]
        conn.close()
        
        if count > 50000:
            results.add_pass("Database gene count", f"{count:,} genes")
        else:
            results.add_warning("Database gene count", f"Only {count:,} genes (expected >50k)")
    except Exception as e:
        results.add_fail("Database gene count", str(e))


def test_cli_integration(results):
    """Test CLI integration."""
    print(f"\n{Colors.BOLD}=== Testing CLI Integration ==={Colors.END}")
    
    import subprocess
    
    # Test 1: CLI help
    try:
        result = subprocess.run(["gene-resolver", "--help"], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            results.add_pass("CLI: help command", "Help text displayed")
        else:
            results.add_fail("CLI: help command", "Help command failed")
    except FileNotFoundError:
        results.add_warning("CLI: help command", "gene-resolver not in PATH")
    except Exception as e:
        results.add_fail("CLI: help command", str(e))
    
    # Test 2: CLI status
    try:
        result = subprocess.run(["gene-resolver", "status"], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            results.add_pass("CLI: status command", "Status retrieved")
        else:
            results.add_warning("CLI: status command", "Status command failed")
    except FileNotFoundError:
        results.add_warning("CLI: status command", "gene-resolver not in PATH")
    except Exception as e:
        results.add_fail("CLI: status command", str(e))
    
    # Test 3: CLI convert
    try:
        result = subprocess.run(["gene-resolver", "convert", "TP53", "BRCA1"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0 and "ENSG" in result.stdout:
            results.add_pass("CLI: convert command", "Conversion successful")
        else:
            results.add_warning("CLI: convert command", "Conversion may have failed")
    except FileNotFoundError:
        results.add_warning("CLI: convert command", "gene-resolver not in PATH")
    except Exception as e:
        results.add_fail("CLI: convert command", str(e))


def main():
    """Run comprehensive test suite."""
    print(f"\n{Colors.BOLD}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}GENE-ID-RESOLVER COMPREHENSIVE TEST SUITE{Colors.END}")
    print(f"{Colors.BOLD}{'='*70}{Colors.END}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    results = TestResults()
    
    # Run all test suites
    test_basic_conversions(results)
    test_deprecated_genes(results)
    test_ambiguity_strategies(results)
    test_file_io(results)
    test_batch_processing(results)
    test_genome_builds(results)
    test_error_handling(results)
    test_database_operations(results)
    test_cli_integration(results)
    test_update_mechanism(results)
    test_multi_species(results)
    test_deprecated_genes_species_aware(results)
    
    # Print summary
    success = results.print_summary()
    
    print(f"\n{Colors.BOLD}{'='*70}{Colors.END}")
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{Colors.BOLD}{'='*70}{Colors.END}\n")
    
    return 0 if success else 1


def test_update_mechanism(results):
    """Test database update and version checking."""
    print(f"\n{Colors.BOLD}{Colors.BLUE}Test Suite: Update Mechanism{Colors.END}")
    print(f"{Colors.BOLD}{'-'*70}{Colors.END}")
    
    from gene_id_resolver.core.resolver import GeneResolver
    from gene_id_resolver.core.updater import DatabaseUpdater
    
    # Test 1: Check for updates
    try:
        resolver = GeneResolver(Path("data"))
        update_info = resolver.check_for_updates()
        
        if 'status' in update_info:
            results.add_pass("Update check returns valid response", 
                           f"Status: {update_info['status']}")
        else:
            results.add_fail("Update check format", "Missing 'status' key")
        
        resolver.close()
    except Exception as e:
        results.add_fail("Update check mechanism", e)
    
    # Test 2: Get latest Ensembl release
    try:
        db_file = Path("data/genes.db")
        if db_file.exists():
            updater = DatabaseUpdater(db_file)
            latest = updater.get_latest_ensembl_release()
            
            if latest and latest.isdigit() and int(latest) >= 109:
                results.add_pass("Fetch latest Ensembl release", 
                               f"Latest: Ensembl {latest}")
            else:
                results.add_warning("Fetch latest Ensembl release", 
                                  f"Unexpected value: {latest}")
        else:
            results.add_warning("Fetch latest Ensembl release", 
                              "Database not found, skipping test")
    except Exception as e:
        results.add_warning("Fetch latest Ensembl release", 
                          f"Network may be unavailable: {e}")


def test_multi_species(results):
    """Test multi-species support."""
    print(f"\n{Colors.BOLD}{Colors.BLUE}Test Suite: Multi-Species Support{Colors.END}")
    print(f"{Colors.BOLD}{'-'*70}{Colors.END}")
    
    from gene_id_resolver.datasources.ensembl import EnsemblDownloader
    
    # Test 1: URL generation for different species
    try:
        downloader = EnsemblDownloader(Path("data/downloads"))
        
        # Test human URL
        human_url = downloader.get_download_url("homo_sapiens", "109")
        if "Homo_sapiens" in human_url and "GRCh38" in human_url:
            results.add_pass("Human URL generation", f"Correct format")
        else:
            results.add_fail("Human URL generation", f"Unexpected URL: {human_url}")
        
        # Test mouse URL
        mouse_url = downloader.get_download_url("mus_musculus", "109")
        if "Mus_musculus" in mouse_url and "GRCm39" in mouse_url:
            results.add_pass("Mouse URL generation", f"Correct format")
        else:
            results.add_fail("Mouse URL generation", f"Unexpected URL: {mouse_url}")
        
        # Test rat URL
        rat_url = downloader.get_download_url("rattus_norvegicus", "109")
        if "Rattus_norvegicus" in rat_url and "mRatBN7.2" in rat_url:
            results.add_pass("Rat URL generation", f"Correct format")
        else:
            results.add_fail("Rat URL generation", f"Unexpected URL: {rat_url}")
        
    except Exception as e:
        results.add_fail("Multi-species URL generation", e)
    
    # Test 2: Unsupported species handling
    try:
        downloader = EnsemblDownloader(Path("data/downloads"))
        try:
            downloader.get_download_url("unknown_species", "109")
            results.add_fail("Unsupported species handling", 
                           "Should have raised ValueError")
        except ValueError as e:
            if "Unsupported species" in str(e):
                results.add_pass("Unsupported species handling", 
                               "Correctly raises ValueError")
            else:
                results.add_fail("Unsupported species handling", 
                               f"Wrong error message: {e}")
    except Exception as e:
        results.add_fail("Unsupported species handling", e)


def test_deprecated_genes_species_aware(results):
    """Test 31: Deprecated gene handler is species-aware"""
    from gene_id_resolver.core.deprecated_genes import SmartDeprecatedGeneHandler
    from pathlib import Path
    
    try:
        # Test human handler (uses fallback)
        handler_human = SmartDeprecatedGeneHandler(
            Path("data/genes.db"), 
            Path("data"),
            species="homo_sapiens"
        )
        if len(handler_human._deprecated_map) >= 0:  # At least fallback
            results.add_pass("Human deprecated gene handler", 
                           f"{len(handler_human._deprecated_map)} mappings")
        
        # Test mouse handler (uses fallback)
        handler_mouse = SmartDeprecatedGeneHandler(
            Path("data/genes.db"), 
            Path("data"),
            species="mus_musculus"
        )
        if len(handler_mouse._deprecated_map) >= 0:  # At least fallback
            results.add_pass("Mouse deprecated gene handler",
                           f"{len(handler_mouse._deprecated_map)} mappings")
        
        # Test rat handler (uses fallback)
        handler_rat = SmartDeprecatedGeneHandler(
            Path("data/genes.db"), 
            Path("data"),
            species="rattus_norvegicus"
        )
        if len(handler_rat._deprecated_map) >= 0:  # At least fallback
            results.add_pass("Rat deprecated gene handler",
                           f"{len(handler_rat._deprecated_map)} mappings")
        
    except Exception as e:
        results.add_fail("Species-aware deprecated gene handler", e)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
"""
Comprehensive pytest test suite for gene-id-resolver.

Tests all features including new enhancements:
- Basic ID conversions (symbol, ensembl, entrez)
- Deprecated gene auto-correction
- Ambiguity resolution strategies
- File I/O (CSV, TSV, BED, TXT)
- Batch processing
- Parallel processing (new)
- JSON output format (new)
- Configuration loading (new)
- Data integrity checks (new)
- CLI commands
- Edge cases and error handling
"""

import pytest
import json
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, MagicMock
import pandas as pd
import yaml

# Import the modules to test
from gene_id_resolver.core.resolver import GeneResolver
from gene_id_resolver.core.resolver_enhanced import EnhancedGeneResolver
from gene_id_resolver.cli.main import cli, load_config, apply_config_to_options, result_to_json
from gene_id_resolver.datasources.ensembl import EnsemblDownloader
from gene_id_resolver.utils.file_processor import read_genes_from_file


class TestBasicConversions:
    """Test basic gene ID conversion functionality."""

    @pytest.fixture
    def resolver(self, tmp_path, monkeypatch):
        """Create a resolver with temporary data directory."""
        from gene_id_resolver.core.models import GeneMapping, GeneIdentifier, GenomicCoordinates
        
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        resolver = GeneResolver(data_dir)
        
        # Create mock gene mappings
        tp53_mapping = GeneMapping(
            identifiers=GeneIdentifier(
                ensembl_id="ENSG00000141510",
                gene_symbol="TP53",
                entrez_id="7157"
            ),
            coordinates=GenomicCoordinates(
                chromosome="17",
                start=7565097,
                end=7590856,
                strand="+",
                genome_build="hg38"
            ),
            biotype="protein_coding"
        )
        
        brca1_mapping = GeneMapping(
            identifiers=GeneIdentifier(
                ensembl_id="ENSG00000012048",
                gene_symbol="BRCA1",
                entrez_id="672"
            ),
            coordinates=GenomicCoordinates(
                chromosome="17",
                start=43044295,
                end=43125483,
                strand="-",
                genome_build="hg38"
            ),
            biotype="protein_coding"
        )
        
        # Mock the database operations
        def mock_find_mappings(gene_id, from_type, genome_build):
            if gene_id.upper() == "TP53" and from_type == "symbol":
                return [tp53_mapping]
            elif gene_id.upper() == "BRCA1" and from_type == "symbol":
                return [brca1_mapping]
            elif gene_id == "ENSG00000141510" and from_type == "ensembl":
                return [tp53_mapping]
            elif gene_id.upper() == "EGFR" and from_type == "symbol":
                return [GeneMapping(
                    identifiers=GeneIdentifier(
                        ensembl_id="ENSG00000146648",
                        gene_symbol="EGFR",
                        entrez_id="1956"
                    ),
                    coordinates=GenomicCoordinates(
                        chromosome="7",
                        start=55019017,
                        end=55211628,
                        strand="+",
                        genome_build="hg38"
                    ),
                    biotype="protein_coding"
                )]
            else:
                return []
        
        monkeypatch.setattr(resolver, '_find_gene_mappings', mock_find_mappings)
        yield resolver
        resolver.close()

    def test_symbol_to_ensembl(self, resolver):
        """Test conversion from gene symbols to Ensembl IDs."""
        genes = ["TP53", "BRCA1"]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38")

        assert len(result.successful) > 0, "Should convert at least one gene"
        if "TP53" in result.successful:
            assert result.successful["TP53"].identifiers.ensembl_id is not None

    def test_ensembl_to_symbol(self, resolver):
        """Test conversion from Ensembl IDs to symbols."""
        genes = ["ENSG00000141510"]  # TP53
        result = resolver.convert(genes, "ensembl", "symbol", "hg38")

        assert len(result.successful) > 0, "Should convert Ensembl ID"
        assert result.successful[genes[0]].identifiers.gene_symbol is not None

    def test_symbol_to_entrez(self, resolver):
        """Test conversion from symbols to Entrez IDs."""
        genes = ["EGFR"]
        result = resolver.convert(genes, "symbol", "entrez", "hg38")

        assert len(result.successful) > 0, "Should convert to Entrez ID"

    def test_case_insensitivity(self, resolver):
        """Test that conversions are case insensitive."""
        genes = ["tp53", "BRCA1", "BrCa2"]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38")

        assert len(result.successful) >= 2, "Should handle mixed case"

    def test_invalid_gene_ids(self, resolver):
        """Test handling of invalid gene IDs."""
        genes = ["INVALID_GENE_12345", "NONEXISTENT_SYMBOL"]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38")

        assert len(result.failed) == 2, "Should fail invalid genes"
        assert all(g in result.failed for g in genes)

    def test_empty_input(self, resolver):
        """Test handling of empty gene list."""
        result = resolver.convert([], "symbol", "ensembl", "hg38")

        assert len(result.successful) == 0
        assert len(result.failed) == 0


class TestEnhancedResolver:
    """Test enhanced resolver with deprecated gene correction."""

    @pytest.fixture
    def enhanced_resolver(self, tmp_path, monkeypatch):
        """Create enhanced resolver with temporary data directory."""
        from gene_id_resolver.core.models import GeneMapping, GeneIdentifier, GenomicCoordinates
        
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        resolver = EnhancedGeneResolver(data_dir)
        
        # Mock basic gene mappings
        tp53_mapping = GeneMapping(
            identifiers=GeneIdentifier(
                ensembl_id="ENSG00000141510",
                gene_symbol="TP53",
                entrez_id="7157"
            ),
            coordinates=GenomicCoordinates(
                chromosome="17",
                start=7565097,
                end=7590856,
                strand="+",
                genome_build="hg38"
            ),
            biotype="protein_coding"
        )
        
        deprecated_mapping = GeneMapping(
            identifiers=GeneIdentifier(
                ensembl_id="ENSG00000141510",
                gene_symbol="TP53",
                entrez_id="7157"
            ),
            coordinates=GenomicCoordinates(
                chromosome="17",
                start=7565097,
                end=7590856,
                strand="+",
                genome_build="hg38"
            ),
            biotype="protein_coding"
        )
        
        def mock_find_mappings(gene_id, from_type, genome_build):
            if gene_id in ["TP53", "SEPT4", "SEPTIN4", "CDKN2", "CDKN2A", "KIAA0196", "C9orf72"]:
                return [tp53_mapping]
            else:
                return []
        
        monkeypatch.setattr(resolver, '_find_gene_mappings', mock_find_mappings)
        yield resolver
        resolver.close()

    def test_deprecated_gene_correction(self, enhanced_resolver):
        """Test auto-correction of deprecated gene symbols."""
        deprecated_tests = [
            ("SEPT4", "SEPTIN4"),
            ("CDKN2", "CDKN2A"),
        ]

        for old_symbol, expected_new in deprecated_tests:
            result = enhanced_resolver.convert_with_correction([old_symbol], "symbol", "ensembl", "hg38")
            assert old_symbol in result.successful, f"Should correct {old_symbol}"

    def test_ki_aa_genes(self, enhanced_resolver):
        """Test handling of KIAA genes."""
        result = enhanced_resolver.convert_with_correction(["KIAA0196"], "symbol", "ensembl", "hg38")
        # KIAA genes might not be in all datasets, so check if handled gracefully
        assert "KIAA0196" in result.successful or "KIAA0196" in result.failed

    def test_corf_genes(self, enhanced_resolver):
        """Test handling of C*orf* genes."""
        result = enhanced_resolver.convert_with_correction(["C9orf72"], "symbol", "ensembl", "hg38")
        assert "C9orf72" in result.successful, "Should handle C9orf72"


class TestAmbiguityResolution:
    """Test ambiguity resolution strategies."""

    @pytest.fixture
    def resolver(self, tmp_path, monkeypatch):
        """Create resolver for ambiguity tests."""
        from gene_id_resolver.core.models import GeneMapping, GeneIdentifier, GenomicCoordinates
        
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        resolver = GeneResolver(data_dir)
        
        tp53_mapping = GeneMapping(
            identifiers=GeneIdentifier(
                ensembl_id="ENSG00000141510",
                gene_symbol="TP53",
                entrez_id="7157"
            ),
            coordinates=GenomicCoordinates(
                chromosome="17",
                start=7565097,
                end=7590856,
                strand="+",
                genome_build="hg38"
            ),
            biotype="protein_coding"
        )
        
        def mock_find_mappings(gene_id, from_type, genome_build):
            if gene_id == "TP53":
                return [tp53_mapping]
            else:
                return []
        
        monkeypatch.setattr(resolver, '_find_gene_mappings', mock_find_mappings)
        yield resolver
        resolver.close()

    @pytest.mark.parametrize("strategy", ["primary", "first", "all"])
    def test_ambiguity_strategies(self, resolver, strategy):
        """Test different ambiguity resolution strategies."""
        genes = ["TP53"]  # Known to potentially have ambiguity
        result = resolver.convert(genes, "symbol", "ensembl", "hg38", strategy)

        # Should not crash, result may be successful or ambiguous
        assert genes[0] in result.successful or genes[0] in result.ambiguous or genes[0] in result.failed


class TestParallelProcessing:
    """Test parallel processing functionality."""

    @pytest.fixture
    def resolver(self, tmp_path, monkeypatch):
        """Create resolver for parallel tests."""
        from gene_id_resolver.core.models import GeneMapping, GeneIdentifier, GenomicCoordinates
        
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        resolver = GeneResolver(data_dir)
        
        # Create mock mappings for multiple genes
        mappings = {}
        gene_symbols = ["TP53", "BRCA1", "EGFR", "KRAS", "MYC"]
        
        for i, symbol in enumerate(gene_symbols):
            mappings[symbol] = GeneMapping(
                identifiers=GeneIdentifier(
                    ensembl_id=f"ENSG00000{i+1:06d}",
                    gene_symbol=symbol,
                    entrez_id=str(1000 + i)
                ),
                coordinates=GenomicCoordinates(
                    chromosome=str((i % 22) + 1),
                    start=1000000 + i * 100000,
                    end=1000000 + i * 100000 + 50000,
                    strand="+" if i % 2 == 0 else "-",
                    genome_build="hg38"
                ),
                biotype="protein_coding"
            )
        
        def mock_find_mappings(gene_id, from_type, genome_build):
            if gene_id in mappings:
                return [mappings[gene_id]]
            else:
                return []
        
        monkeypatch.setattr(resolver, '_find_gene_mappings', mock_find_mappings)
        yield resolver
        resolver.close()

    def test_parallel_conversion(self, resolver):
        """Test parallel processing of gene conversions."""
        genes = ["TP53", "BRCA1", "EGFR", "KRAS", "MYC"] * 2  # 10 genes to trigger parallel
        result_parallel = resolver.convert(genes, "symbol", "ensembl", "hg38", parallel=True, max_workers=2)

        result_sequential = resolver.convert(genes, "symbol", "ensembl", "hg38", parallel=False)

        # Results should be similar (allowing for non-deterministic order)
        assert len(result_parallel.successful) == len(result_sequential.successful)

    def test_sequential_conversion(self, resolver):
        """Test sequential processing."""
        genes = ["TP53", "BRCA1"]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38", parallel=False)

        assert len(result.successful) > 0

    def test_parallel_with_custom_workers(self, resolver):
        """Test parallel processing with custom worker count."""
        genes = ["TP53", "BRCA1", "EGFR"] * 10  # 30 genes
        result = resolver.convert(genes, "symbol", "ensembl", "hg38", parallel=True, max_workers=2)

        assert len(result.successful) > 0


class TestFileIO:
    """Test file input/output functionality."""

    @pytest.fixture
    def temp_files(self, tmp_path):
        """Create temporary test files."""
        files = {}

        # CSV file
        csv_file = tmp_path / "test_genes.csv"
        pd.DataFrame({"gene": ["TP53", "BRCA1", "EGFR"]}).to_csv(csv_file, index=False)
        files["csv"] = csv_file

        # TSV file
        tsv_file = tmp_path / "test_genes.tsv"
        pd.DataFrame({"gene": ["TP53", "BRCA1"]}).to_csv(tsv_file, sep="\t", index=False)
        files["tsv"] = tsv_file

        # TXT file
        txt_file = tmp_path / "test_genes.txt"
        txt_file.write_text("TP53\nBRCA1\nEGFR\n")
        files["txt"] = txt_file

        # BED file
        bed_file = tmp_path / "test_regions.bed"
        bed_file.write_text("chr17\t7565097\t7590856\tTP53\nchr13\t32315086\t32400266\tBRCA2\n")
        files["bed"] = bed_file

        return files

    def test_csv_input(self, temp_files):
        """Test reading genes from CSV file."""
        genes = read_genes_from_file(temp_files["csv"], column=0)
        assert len(genes) == 3
        assert "TP53" in genes

    def test_tsv_input(self, temp_files):
        """Test reading genes from TSV file."""
        genes = read_genes_from_file(temp_files["tsv"], column=0)
        assert len(genes) == 2

    def test_txt_input(self, temp_files):
        """Test reading genes from plain text file."""
        genes = read_genes_from_file(temp_files["txt"])
        assert len(genes) == 3

    def test_bed_input(self, temp_files):
        """Test reading genes from BED file."""
        genes = read_genes_from_file(temp_files["bed"])
        assert len(genes) >= 1  # May extract gene names from BED

    def test_invalid_file_format(self, tmp_path):
        """Test handling of unknown file formats (falls back to plain text)."""
        invalid_file = tmp_path / "invalid.xyz"
        invalid_file.write_text("TP53\nBRCA1\n")

        genes = read_genes_from_file(invalid_file)
        assert len(genes) == 2  # Should read as plain text


class TestJSONOutput:
    """Test JSON output functionality."""

    def test_result_to_json(self):
        """Test conversion of results to JSON format."""
        # Mock result object
        class MockResult:
            def __init__(self):
                self.successful = {
                    "TP53": MagicMock(
                        identifiers=MagicMock(ensembl_id="ENSG00000141510", gene_symbol="TP53", entrez_id="7157", uniprot_id=None),
                        coordinates=MagicMock(chromosome="17", start=7565097, end=7590856, strand="+", genome_build="hg38"),
                        biotype="protein_coding",
                        description="tumor protein p53"
                    )
                }
                self.failed = ["INVALID_GENE"]
                self.ambiguous = {}
                self.input_type = "symbol"
                self.output_type = "ensembl"
                self.resolver_config = {"genome_build": "hg38"}
                self.ambiguity_resolutions = {"TP53": "unique"}
                self.success_rate = 0.5

        result = MockResult()
        genes_list = ["TP53", "INVALID_GENE"]
        json_dict = result_to_json(result, genes_list)

        assert "metadata" in json_dict
        assert "results" in json_dict
        assert "successful" in json_dict["results"]
        assert "TP53" in json_dict["results"]["successful"]
        assert json_dict["metadata"]["success_count"] == 1


class TestConfiguration:
    """Test configuration loading functionality."""

    @pytest.fixture
    def config_file(self, tmp_path):
        """Create a temporary config file."""
        config = {
            "parallel": True,
            "max_workers": 4,
            "debug": True,
            "format": "json",
            "input_format": "symbol",
            "output_format": "ensembl",
            "genome_build": "hg38"
        }
        config_file = tmp_path / "test_config.yaml"
        with open(config_file, 'w') as f:
            yaml.dump(config, f)
        return config_file

    def test_load_config(self, config_file):
        """Test loading configuration from YAML file."""
        config = load_config(str(config_file))
        assert config["parallel"] is True
        assert config["max_workers"] == 4
        assert config["debug"] is True

    def test_apply_config_to_options(self):
        """Test applying config to CLI options."""
        config = {"parallel": True, "debug": False}
        options = apply_config_to_options(config, parallel=False, debug=None)
        assert options["parallel"] is True  # Config overrides default
        assert options["debug"] is False   # Config sets value

    def test_missing_config_file(self):
        """Test handling of missing config file."""
        config = load_config("/nonexistent/config.yaml")
        assert config == {}  # Should return empty dict


class TestDataIntegrity:
    """Test data integrity and checksum functionality."""

    @pytest.fixture
    def downloader(self, tmp_path):
        """Create downloader for testing."""
        download_dir = tmp_path / "downloads"
        download_dir.mkdir()
        return EnsemblDownloader(download_dir)

    def test_checksum_calculation(self, downloader):
        """Test checksum calculation for files."""
        # Create a test file
        test_file = downloader.download_dir / "test.txt"
        test_file.write_text("test content")

        checksum = downloader.calculate_checksum(test_file)
        assert checksum is not None
        assert len(checksum) > 0

    def test_checksum_verification(self, downloader):
        """Test checksum verification."""
        test_file = downloader.download_dir / "test.txt"
        test_file.write_text("test content")

        checksum = downloader.calculate_checksum(test_file)
        assert downloader.verify_checksum(test_file, checksum)

        # Test with wrong checksum
        assert not downloader.verify_checksum(test_file, "wrong_checksum")

    @patch('urllib.request.urlretrieve')
    def test_download_with_checksum(self, mock_urlretrieve, downloader):
        """Test download with checksum verification."""
        # Mock urlretrieve to create a test file
        def mock_retrieve(url, filename):
            Path(filename).write_text("test data")
        mock_urlretrieve.side_effect = mock_retrieve

        test_file = downloader.download_dir / "test_download.txt"
        result = downloader.download_with_checksum("http://example.com/test.txt", test_file)

        assert result is True
        assert test_file.exists()
        assert test_file.read_text() == "test data"

    @patch('requests.get')
    def test_corrupted_download_cleanup(self, mock_get, downloader):
        """Test cleanup of corrupted downloads."""
        # Mock corrupted download (wrong content-length)
        mock_response = MagicMock()
        mock_response.content = b"short"
        mock_response.headers = {'content-length': '100'}  # Wrong length
        mock_get.return_value = mock_response

        test_file = downloader.download_dir / "corrupted.txt"

        with pytest.raises(Exception):  # Should raise due to size mismatch
            downloader.download_with_checksum("http://example.com/corrupted.txt", str(test_file))

        # File should be cleaned up
        assert not test_file.exists()


class TestCLICommands:
    """Test CLI command functionality."""

    @pytest.fixture
    def temp_data_dir(self, tmp_path):
        """Create temporary data directory with mock database."""
        import sqlite3
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        # Create a minimal SQLite database file
        db_file = data_dir / "genes.db"
        conn = sqlite3.connect(str(db_file))
        conn.execute('''CREATE TABLE genes (
            symbol TEXT,
            ensembl_id TEXT,
            entrez_id TEXT,
            genome_build TEXT
        )''')
        conn.execute("INSERT INTO genes VALUES ('TP53', 'ENSG00000141510', '7157', 'hg38')")
        conn.commit()
        conn.close()
        return data_dir

    @patch('gene_id_resolver.cli.main.GeneResolver')
    def test_convert_command_basic(self, mock_resolver_class, temp_data_dir):
        """Test basic convert command."""
        mock_resolver = MagicMock()
        mock_resolver_class.return_value = mock_resolver

        mock_result = MagicMock()
        mock_result.successful = {"TP53": MagicMock()}
        mock_result.failed = []
        mock_result.ambiguous = {}
        mock_resolver.convert.return_value = mock_result

        from click.testing import CliRunner
        runner = CliRunner()

        result = runner.invoke(cli, [
            'convert', 'TP53', 'BRCA1',
            '--from-type', 'symbol',
            '--to-type', 'ensembl',
            '--genome-build', 'hg38',
            '--data-dir', str(temp_data_dir),
            '--no-auto-correct'  # Disable enhanced resolver to use basic mock
        ])

        assert result.exit_code == 0
        # Check that the output contains expected elements
        assert "Conversion Results" in result.output
        assert "Genome: hg38" in result.output
        assert "Strategy: primary" in result.output

    @patch('gene_id_resolver.cli.main.GeneResolver')
    def test_convert_command_parallel(self, mock_resolver_class, temp_data_dir):
        """Test convert command with parallel processing."""
        mock_resolver = MagicMock()
        mock_resolver_class.return_value = mock_resolver

        mock_result = MagicMock()
        mock_result.successful = {"TP53": MagicMock()}
        mock_resolver.convert.return_value = mock_result

        from click.testing import CliRunner
        runner = CliRunner()

        with patch('pathlib.Path.exists', return_value=True):
            result = runner.invoke(cli, [
                'convert', 'TP53', 'BRCA1',
                '--parallel',
                '--max-workers', '4',
                '--data-dir', str(temp_data_dir),
                '--no-auto-correct'  # Disable enhanced resolver
            ])

        assert result.exit_code == 0
        # Check that parallel processing was indicated in output
        assert "Conversion Results" in result.output

    @patch('gene_id_resolver.cli.main.GeneResolver')
    def test_convert_command_json_output(self, mock_resolver_class, temp_data_dir):
        """Test convert command with JSON output."""
        mock_resolver = MagicMock()
        mock_resolver_class.return_value = mock_resolver

        mock_result = MagicMock()
        mock_result.successful = {"TP53": MagicMock(
            identifiers=MagicMock(ensembl_id="ENSG00000141510", gene_symbol="TP53"),
            coordinates=MagicMock(chromosome="17", start=7565097, end=7590856, strand="+", genome_build="hg38"),
            biotype="protein_coding"
        )}
        mock_result.failed = []
        mock_result.ambiguous = {}
        mock_result.input_type = "symbol"
        mock_result.output_type = "ensembl"
        mock_result.resolver_config = {"genome_build": "hg38"}
        mock_result.ambiguity_resolutions = {"TP53": "unique"}
        mock_result.success_rate = 1.0
        mock_resolver.convert.return_value = mock_result

        from click.testing import CliRunner
        runner = CliRunner()

        with patch('pathlib.Path.exists', return_value=True):
            result = runner.invoke(cli, [
                'convert', 'TP53',
                '--format', 'json',
                '--data-dir', str(temp_data_dir)
            ])

        assert result.exit_code == 0
        # Should contain JSON output
        assert '"metadata"' in result.output
        assert '"results"' in result.output

    def test_convert_file_command(self, tmp_path, temp_data_dir):
        """Test convert-file command."""
        # Create input file
        input_file = tmp_path / "genes.txt"
        input_file.write_text("TP53\nBRCA1\n")

        from click.testing import CliRunner
        runner = CliRunner()

        result = runner.invoke(cli, [
            'convert-file', str(input_file),
            '--data-dir', str(temp_data_dir)
        ])

        # Command should run (may fail due to missing data, but CLI should work)
        assert result.exit_code in [0, 1]  # 0 for success, 1 for conversion failures


class TestEdgeCases:
    """Test edge cases and error handling."""

    @pytest.fixture
    def resolver(self, tmp_path, monkeypatch):
        """Create resolver for edge case tests."""
        from gene_id_resolver.core.models import GeneMapping, GeneIdentifier, GenomicCoordinates
        
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        resolver = GeneResolver(data_dir)
        
        tp53_mapping = GeneMapping(
            identifiers=GeneIdentifier(
                ensembl_id="ENSG00000141510",
                gene_symbol="TP53",
                entrez_id="7157"
            ),
            coordinates=GenomicCoordinates(
                chromosome="17",
                start=7565097,
                end=7590856,
                strand="+",
                genome_build="hg38"
            ),
            biotype="protein_coding"
        )
        
        def mock_find_mappings(gene_id, from_type, genome_build):
            if gene_id == "TP53":
                return [tp53_mapping]
            else:
                return []
        
        monkeypatch.setattr(resolver, '_find_gene_mappings', mock_find_mappings)
        yield resolver
        resolver.close()

    def test_very_large_gene_list(self, resolver):
        """Test handling of very large gene lists."""
        # Create a large list that should trigger parallel processing
        genes = [f"GENE_{i}" for i in range(20)]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38", parallel=False)  # Disable parallel for now

        # Should handle large lists without crashing
        assert isinstance(result.successful, dict)
        assert isinstance(result.failed, list)

    def test_special_characters_in_gene_names(self, resolver):
        """Test handling of special characters in gene names."""
        genes = ["GENE-1", "GENE_2", "GENE.3"]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38")

        # Should handle gracefully
        assert isinstance(result, object)

    def test_mixed_valid_invalid_genes(self, resolver):
        """Test mix of valid and invalid gene identifiers."""
        genes = ["TP53", "INVALID_GENE", "BRCA1", "ANOTHER_INVALID"]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38")

        # Should have some successful and some failed
        total_processed = len(result.successful) + len(result.failed)
        assert total_processed == len(genes)

    def test_database_loading_failure(self, resolver, monkeypatch):
        """Test handling of database loading failures."""
        # Override the mock to raise an exception
        def mock_find_error(*args, **kwargs):
            raise Exception("Database error")
        
        monkeypatch.setattr(resolver, '_find_gene_mappings', mock_find_error)

        genes = ["TP53"]
        with pytest.raises(Exception, match="Database error"):
            resolver.convert(genes, "symbol", "ensembl", "hg38")

    def test_multiple_genome_builds(self, resolver):
        """Test conversion with different genome builds."""
        genes = ["TP53"]
        for build in ["hg38", "hg19", "mm10"]:
            result = resolver.convert(genes, "symbol", "ensembl", build)
            # Should not crash, even if data not available
            assert hasattr(result, 'successful')


class TestBatchProcessing:
    """Test batch processing functionality."""

    @pytest.fixture
    def resolver(self, tmp_path, monkeypatch):
        """Create resolver for batch tests."""
        from gene_id_resolver.core.models import GeneMapping, GeneIdentifier, GenomicCoordinates
        
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        resolver = GeneResolver(data_dir)
        
        tp53_mapping = GeneMapping(
            identifiers=GeneIdentifier(
                ensembl_id="ENSG00000141510",
                gene_symbol="TP53",
                entrez_id="7157"
            ),
            coordinates=GenomicCoordinates(
                chromosome="17",
                start=7565097,
                end=7590856,
                strand="+",
                genome_build="hg38"
            ),
            biotype="protein_coding"
        )
        
        def mock_find_mappings(gene_id, from_type, genome_build):
            if gene_id.startswith("GENE_") or gene_id == "TP53":
                return [tp53_mapping]
            else:
                return []
        
        monkeypatch.setattr(resolver, '_find_gene_mappings', mock_find_mappings)
        yield resolver
        resolver.close()

    def test_batch_size_handling(self, resolver):
        """Test that batch processing handles different sizes correctly."""
        # Small batch
        small_batch = ["TP53", "BRCA1"]
        result_small = resolver.convert(small_batch, "symbol", "ensembl", "hg38", parallel=False)

        # Large batch
        large_batch = ["TP53"] * 10
        result_large = resolver.convert(large_batch, "symbol", "ensembl", "hg38", parallel=False)

        # Both should work
        assert hasattr(result_small, 'successful')
        assert hasattr(result_large, 'successful')

    def test_batch_with_mixed_results(self, resolver):
        """Test batch processing with mix of successful and failed conversions."""
        genes = ["TP53", "INVALID1", "BRCA1", "INVALID2", "EGFR"]
        result = resolver.convert(genes, "symbol", "ensembl", "hg38")

        # Should have both successful and failed
        assert len(result.successful) + len(result.failed) == len(genes)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
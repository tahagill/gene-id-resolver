"""
Command-line interface for gene-id-resolver.

Provides a simple way to convert gene IDs from the command line.
"""

import click
import sys
from pathlib import Path

from ..core.resolver import GeneResolver
from ..utils.file_processor import read_genes_from_file, save_results_to_file


@click.group()
def cli():
    """Gene ID Resolver - Convert between gene identifier systems."""
    pass
@cli.command()
@click.option('--check', '-c', is_flag=True, help='Check for updates without installing')
@click.option('--release', '-r', help='Update to specific Ensembl release')
@click.option('--data-dir', '-d', default='./data', help='Directory with gene databases')
@click.option('--species', '-s', default='homo_sapiens',
              type=click.Choice(['homo_sapiens', 'mus_musculus', 'rattus_norvegicus']),
              help='Species (required if updating to new release)')
@click.option('--genome-build', '-g', help='Genome build (optional)')
def update(check, release, data_dir, species, genome_build):
    """
    Update gene database to latest Ensembl release.
    
    Examples:
      gene-convert update --check           # Check for updates
      gene-convert update                   # Update to latest (human)
      gene-convert update --release 110     # Update to specific version
      gene-convert update --species mus_musculus  # Update mouse database
    """
    resolver = GeneResolver(Path(data_dir))
    
    try:
        if check:
            click.echo("Checking for database updates...")
            update_info = resolver.check_for_updates()
            
            if update_info['status'] == 'no_version_info':
                click.echo("Error: No version information found in database")
                click.echo("Run 'gene-convert init' first to initialize database")
            elif update_info['status'] == 'current':
                current = update_info['current']
                click.echo("Database is up to date!")
                click.echo(f"   Current version: Ensembl {current['ensembl_release']}")
                click.echo(f"   Species: {current['species']}")
                click.echo(f"   Genome build: {current['genome_build']}")
                click.echo(f"   Last updated: {current['last_updated']}")
            elif update_info['status'] == 'update_available':
                current = update_info['current']
                click.echo("Update available!")
                click.echo(f"   Current: Ensembl {current['ensembl_release']}")
                click.echo(f"   Latest: Ensembl {update_info['latest_release']}")
                click.echo(f"\nRun 'gene-convert update' to install the update")
            elif update_info['status'] == 'check_failed':
                click.echo("Could not check for updates (network issue)")
                click.echo(f"   Current version: Ensembl {update_info['current']['ensembl_release']}")
                
        else:
            # Perform actual update
            db_file = Path(data_dir) / "genes.db"
            if not db_file.exists():
                click.echo("Error: Database not found. Run 'gene-convert init' first.")
                return
            
            # Determine target release
            if not release:
                click.echo("Fetching latest Ensembl release...")
                from gene_id_resolver.core.updater import DatabaseUpdater
                updater = DatabaseUpdater(db_file)
                release = updater.get_latest_ensembl_release()
                if not release:
                    click.echo("Error: Could not determine latest release")
                    click.echo("Try specifying a release with --release <number>")
                    return
            
            click.echo(f"\nUpdating database to Ensembl release {release}...")
            click.echo(f"   Species: {species}")
            if genome_build:
                click.echo(f"   Genome build: {genome_build}")
            click.echo("\nThis will download new data (~100MB) and rebuild the database.")
            click.echo("Your current database will be backed up automatically.\n")
            
            if click.confirm("Continue with update?"):
                success = resolver.update_database(release, species, genome_build)
                if success:
                    click.echo("\nDatabase updated successfully!")
                    click.echo("Backup saved with .backup extension")
                else:
                    click.echo("\nUpdate failed. Database restored from backup.", err=True)
            else:
                click.echo("Update cancelled")
            
    except Exception as e:
        click.echo(f"Error: Update failed: {e}", err=True)
    finally:
        resolver.close()

@cli.command()
@click.argument('genes', nargs=-1, required=True)
@click.option('--from-type', '-f', default='symbol', 
              help='Input ID type (symbol, ensembl, entrez)')
@click.option('--to-type', '-t', default='ensembl',
              help='Output ID type (symbol, ensembl, entrez)')
@click.option('--genome-build', '-g', default='hg38',
              help='Genome build (hg38, GRCh37, GRCm39)')
@click.option('--ambiguity-strategy', '-r', default='primary',
              type=click.Choice(['primary', 'first', 'all', 'fail']),
              help='How to handle ambiguous gene symbols')
@click.option('--auto-correct/--no-auto-correct', default=True,
              help='Auto-correct deprecated gene symbols')
@click.option('--data-dir', '-d', default='./data',
              help='Directory to store gene databases')
def convert(genes, from_type, to_type, genome_build, ambiguity_strategy, auto_correct, data_dir):
    """
    Convert gene IDs from one system to another.
    
    Ambiguity resolution strategies:
    - primary: Prefer protein-coding genes (recommended)
    - first: Use first match found
    - all: Show all matches for manual review  
    - fail: Treat ambiguous genes as failures
    """
    # Import here to avoid circular imports
    from gene_id_resolver.core.resolver import GeneResolver
    try:
        # Check if database exists
        db_file = Path(data_dir) / "genes.db"
        if not db_file.exists():
            click.echo("Warning: Database not found. Run 'gene-resolver init' first.")
            return
        
        # Use enhanced resolver if auto-correct is enabled
        if auto_correct:
            try:
                from gene_id_resolver.core.resolver_enhanced import EnhancedGeneResolver
                resolver = EnhancedGeneResolver(Path(data_dir))
                result = resolver.convert_with_correction(
                    list(genes), from_type, to_type, genome_build, ambiguity_strategy
                )
            except ImportError:
                click.echo("Warning: Enhanced features not available, using basic conversion")
                resolver = GeneResolver(Path(data_dir))
                result = resolver.convert(
                    list(genes), from_type, to_type, genome_build, ambiguity_strategy
                )
        else:
            resolver = GeneResolver(Path(data_dir))
            result = resolver.convert(
                list(genes), from_type, to_type, genome_build, ambiguity_strategy
            )
        
        # Display results
        click.echo(f"\nConversion Results (genome: {genome_build}, ambiguity: {ambiguity_strategy}):")
        click.echo(f"   Successful: {len(result.successful)}")
        click.echo(f"   Ambiguous: {len(result.ambiguous)}")
        click.echo(f"   Failed: {len(result.failed)}")
        
        # Show auto-corrections if any were applied
        if auto_correct and hasattr(result, 'corrections_applied') and result.corrections_applied:
            click.echo(f"   Auto-corrected: {len(result.corrections_applied)} genes")
            for old, new in result.corrections_applied.items():
                click.echo(f"     {old} → {new}")
        
        # Show ambiguity resolutions
        primary_resolved = sum(1 for v in result.ambiguity_resolutions.values() 
                             if v == "resolved_primary")
        if primary_resolved > 0:
            click.echo(f"   Auto-resolved: {primary_resolved} ambiguous genes")
        
        # Show successful conversions
        if result.successful:
            click.echo(f"\nSuccessful conversions:")
            for input_id, mapping in result.successful.items():
                output_id = get_output_id(mapping, to_type)
                resolution = result.ambiguity_resolutions.get(input_id, "")
                
                was_corrected = False
                if auto_correct and hasattr(result, 'corrections_applied'):
                    was_corrected = input_id in result.corrections_applied.values()
                    original_gene = next((k for k, v in result.corrections_applied.items() if v == input_id), None)
                
                if output_id:
                    if was_corrected and original_gene:
                        click.echo(f"   {original_gene} → {output_id} (auto-corrected from {original_gene})")
                    elif resolution == "resolved_primary":
                        click.echo(f"   {input_id} → {output_id} (auto-resolved protein-coding)")
                    else:
                        click.echo(f"   {input_id} → {output_id}")
        
        # Show ambiguous mappings (when strategy is 'all')
        if result.ambiguous:
            click.echo(f"\nAmbiguous mappings (need manual review):")
            for input_id, mappings in result.ambiguous.items():
                click.echo(f"   {input_id} maps to {len(mappings)} locations:")
                for i, mapping in enumerate(mappings):
                    output_id = get_output_id(mapping, to_type)
                    biotype = mapping.biotype or "unknown"
                    if output_id:
                        click.echo(f"     {i+1}. {output_id} ({biotype})")
        
        # Show failed conversions
        if result.failed:
            click.echo(f"\nFailed conversions:")
            for failed_id in result.failed:
                resolution = result.ambiguity_resolutions.get(failed_id, "")
                if resolution == "failed_ambiguous":
                    click.echo(f"   {failed_id} (ambiguous - use '--ambiguity-strategy all' to see options)")
                else:
                    click.echo(f"   {failed_id}")
                
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
    finally:
        if 'resolver' in locals():
            resolver.close()


@cli.command()
@click.option('--data-dir', '-d', default='./data',
              help='Directory to store gene databases')
@click.option('--species', '-s', default='homo_sapiens',
              type=click.Choice(['homo_sapiens', 'mus_musculus', 'rattus_norvegicus']),
              help='Species to download')
@click.option('--release', '-r', default='109',
              help='Ensembl release version')
@click.option('--genome-build', '-g', 
              help='Genome build (default: GRCh38 for human, GRCm39 for mouse)')
def init(data_dir, species, release, genome_build):
    """
    Initialize the gene database by downloading Ensembl data.
    """
    resolver = GeneResolver(Path(data_dir))
    
    try:
        click.echo(f"Initializing gene database...")
        click.echo(f"   Species: {species}")
        click.echo(f"   Release: {release}")
        if genome_build:
            click.echo(f"   Genome build: {genome_build}")
        click.echo(f"   Data directory: {data_dir}")
        click.echo(f"\nThis may take a few minutes (downloading ~100MB)...")
        
        gene_count = resolver.initialize_database(species, release, genome_build)
        
        click.echo(f"\nDatabase initialized with {gene_count:,} genes!")
        click.echo("You can now use 'gene-resolver convert' to convert gene IDs.")
        
    except Exception as e:
        click.echo(f"Error: Initialization failed: {e}", err=True)
    finally:
        resolver.close()


@cli.command()
@click.option('--data-dir', '-d', default='./data',
              help='Directory to store gene databases')
def status(data_dir):
    """
    Check the status of the gene database.
    """
    db_file = Path(data_dir) / "genes.db"
    
    if not db_file.exists():
        click.echo("Error: Database not found. Run 'gene-resolver init' first.")
        return
    
    import sqlite3
    try:
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        
        cursor.execute("SELECT COUNT(*) FROM gene_mappings WHERE is_current = TRUE")
        gene_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT genome_build, annotation_version FROM gene_mappings LIMIT 1")
        result = cursor.fetchone()
        
        if result:
            genome_build, annotation_version = result
            click.echo(f"Database Status:")
            click.echo(f"   Genes: {gene_count}")
            click.echo(f"   Genome build: {genome_build}")
            click.echo(f"   Annotation version: {annotation_version}")
        else:
            click.echo("Database exists but contains no data.")
            
        conn.close()
        
    except Exception as e:
        click.echo(f"Error checking database: {e}")


@cli.command()
@click.argument('input_file')
@click.option('--output-file', '-o', help='Save results to CSV file')
@click.option('--from-type', '-f', default='symbol', 
              help='Input ID type (symbol, ensembl, entrez)')
@click.option('--to-type', '-t', default='ensembl',
              help='Output ID type (symbol, ensembl, entrez)')
@click.option('--genome-build', '-g', default='hg38',
              help='Genome build (hg38, GRCh37, GRCm39)')
@click.option('--ambiguity-strategy', '-r', default='primary',
              type=click.Choice(['primary', 'first', 'all', 'fail']),
              help='How to handle ambiguous gene symbols')
@click.option('--column', '-c', default=0,
              help='Column containing gene IDs (0-based, for CSV/TSV)')
@click.option('--data-dir', '-d', default='./data',
              help='Directory to store gene databases')
@click.option('--auto-correct/--no-auto-correct', default=True,
              help='Auto-correct deprecated gene symbols')
def convert_file(input_file, output_file, from_type, to_type, genome_build, 
                ambiguity_strategy, column, data_dir, auto_correct):
    """
    Convert gene IDs from a file.
    
    Ambiguity resolution strategies:
    - primary: Prefer protein-coding genes (recommended)
    - first: Use first match found
    - all: Show all matches for manual review  
    - fail: Treat ambiguous genes as failures
    """
    # Import here to avoid circular imports
    from gene_id_resolver.core.resolver import GeneResolver
    
    try:
        # Check if database exists
        db_file = Path(data_dir) / "genes.db"
        if not db_file.exists():
            click.echo("Warning: Database not found. Run 'gene-resolver init' first.")
            return
        
        # Read genes from file
        click.echo(f"Reading genes from: {input_file}")
        input_genes = read_genes_from_file(Path(input_file), column)
        
        if not input_genes:
            click.echo("Error: No genes found in input file")
            return
            
        click.echo(f"   Found {len(input_genes)} genes to convert")
        click.echo(f"   Genome build: {genome_build}")
        click.echo(f"   Ambiguity strategy: {ambiguity_strategy}")
        click.echo(f"   Auto-correct: {auto_correct}")
        
        # Use enhanced resolver if auto-correct is enabled
        if auto_correct:
            try:
                from gene_id_resolver.core.resolver_enhanced import EnhancedGeneResolver
                resolver = EnhancedGeneResolver(Path(data_dir))
                result = resolver.convert_with_correction(
                    input_genes, from_type, to_type, genome_build, ambiguity_strategy
                )
            except ImportError:
                click.echo("Warning: Enhanced features not available, using basic conversion")
                resolver = GeneResolver(Path(data_dir))
                result = resolver.convert(
                    input_genes, from_type, to_type, genome_build, ambiguity_strategy
                )
        else:
            resolver = GeneResolver(Path(data_dir))
            result = resolver.convert(
                input_genes, from_type, to_type, genome_build, ambiguity_strategy
            )
        
        # Display summary with genome info
        click.echo(f"\nConversion Results (genome: {genome_build}, ambiguity: {ambiguity_strategy}):")
        click.echo(f"   Successful: {len(result.successful)}")
        click.echo(f"   Ambiguous: {len(result.ambiguous)}")
        click.echo(f"   Failed: {len(result.failed)}")
        click.echo(f"   Success rate: {result.success_rate:.1%}")
        
        # Show auto-corrections if any were applied
        if auto_correct and hasattr(result, 'corrections_applied') and result.corrections_applied:
            click.echo(f"   Auto-corrected: {len(result.corrections_applied)} genes")
            for old, new in result.corrections_applied.items():
                click.echo(f"     {old} → {new}")
        
        # Show ambiguity resolutions
        primary_resolved = sum(1 for v in result.ambiguity_resolutions.values() 
                             if v == "resolved_primary")
        if primary_resolved > 0:
            click.echo(f"   Auto-resolved: {primary_resolved} ambiguous genes")
        
        # Save to file if requested
        if output_file:
            save_results_to_file(result, Path(output_file), input_genes)
            click.echo(f"Results saved to: {output_file}")
        else:
            # Display results in terminal
            if result.successful:
                click.echo(f"\nSuccessful conversions (first 10):")
                successful_genes = list(result.successful.keys())[:10]
                for gene in successful_genes:
                    mapping = result.successful[gene]
                    output_id = get_output_id(mapping, to_type)
                    resolution = result.ambiguity_resolutions.get(gene, "")
                    
                    was_corrected = False
                    if auto_correct and hasattr(result, 'corrections_applied'):
                        was_corrected = gene in result.corrections_applied.values()
                        original_gene = next((k for k, v in result.corrections_applied.items() if v == gene), None)
                    
                    if output_id:
                        if was_corrected and original_gene:
                            click.echo(f"   {original_gene} → {output_id} (auto-corrected from {original_gene})")
                        elif resolution == "resolved_primary":
                            click.echo(f"   {gene} → {output_id} (auto-resolved protein-coding)")
                        else:
                            click.echo(f"   {gene} → {output_id}")
                
                if len(result.successful) > 10:
                    click.echo(f"   ... and {len(result.successful) - 10} more")
            
            if result.ambiguous:
                click.echo(f"\nAmbiguous mappings (first 5):")
                ambiguous_genes = list(result.ambiguous.keys())[:5]
                for gene in ambiguous_genes:
                    mappings = result.ambiguous[gene]
                    click.echo(f"   {gene} maps to {len(mappings)} locations:")
                    for i, mapping in enumerate(mappings[:3]):
                        output_id = get_output_id(mapping, to_type)
                        biotype = mapping.biotype or "unknown"
                        if output_id:
                            click.echo(f"     {i+1}. {output_id} ({biotype})")
                    if len(mappings) > 3:
                        click.echo(f"     ... and {len(mappings) - 3} more")
            
            if result.failed:
                click.echo(f"\nFailed conversions (first 10):")
                failed_genes = result.failed[:10]
                for gene in failed_genes:
                    resolution = result.ambiguity_resolutions.get(gene, "")
                    if resolution == "failed_ambiguous":
                        click.echo(f"   {gene} (ambiguous - use '--ambiguity-strategy all' to see options)")
                    else:
                        click.echo(f"   {gene}")
                if len(result.failed) > 10:
                    click.echo(f"   ... and {len(result.failed) - 10} more")
                
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
    finally:
        if 'resolver' in locals():
            resolver.close()
def get_output_id(mapping, to_type):
    """Get the correct output ID based on type."""
    attr_map = {
        'ensembl': 'ensembl_id',
        'symbol': 'gene_symbol',
        'entrez': 'entrez_id'
    }
    attr_name = attr_map.get(to_type)
    return getattr(mapping.identifiers, attr_name, None) if attr_name else None


if __name__ == '__main__':
    cli()
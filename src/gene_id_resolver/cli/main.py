"""
Command-line interface for gene-id-resolver.

Provides a simple way to convert gene IDs from the command line.
"""

import click
from pathlib import Path
import json
import yaml

from ..core.resolver import GeneResolver
from ..utils.file_processor import read_genes_from_file, save_results_to_file


def load_config(config_file=None):
    """Load configuration from YAML file."""
    if not config_file:
        # Look for default config files
        default_configs = [
            Path.home() / '.gene_resolver.yaml',
            Path.home() / '.gene_resolver.yml',
            Path.cwd() / 'gene_resolver.yaml',
            Path.cwd() / 'gene_resolver.yml'
        ]
        for config_path in default_configs:
            if config_path.exists():
                config_file = config_path
                break

    if not config_file or not Path(config_file).exists():
        return {}

    try:
        with open(config_file, 'r') as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        click.echo(f"Warning: Could not load config file {config_file}: {e}", err=True)
        return {}


def apply_config_to_options(config, **kwargs):
    """Apply configuration defaults to CLI options."""
    result = {}
    for key, value in kwargs.items():
        if value is None or (isinstance(value, bool) and not value):
            # Use config default if CLI option not provided
            config_key = key.replace('_', '-')
            if config_key in config:
                result[key] = config[config_key]
            else:
                result[key] = value
        else:
            # CLI option takes precedence
            result[key] = value
    return result


def result_to_json(result, genes_list):
    """Convert conversion result to detailed JSON format."""
    def gene_mapping_to_dict(mapping):
        return {
            "identifiers": {
                "ensembl_id": mapping.identifiers.ensembl_id,
                "gene_symbol": mapping.identifiers.gene_symbol,
                "entrez_id": mapping.identifiers.entrez_id,
                "uniprot_id": mapping.identifiers.uniprot_id
            },
            "coordinates": {
                "chromosome": mapping.coordinates.chromosome,
                "start": mapping.coordinates.start,
                "end": mapping.coordinates.end,
                "strand": mapping.coordinates.strand,
                "genome_build": mapping.coordinates.genome_build
            } if mapping.coordinates else None,
            "gene_type": getattr(mapping, 'gene_type', None),
            "description": getattr(mapping, 'description', None)
        }

    successful = {}
    for gene, mapping in result.successful.items():
        successful[gene] = gene_mapping_to_dict(mapping)

    ambiguous = {}
    for gene, mappings in result.ambiguous.items():
        ambiguous[gene] = [gene_mapping_to_dict(m) for m in mappings]

    return {
        "metadata": {
            "input_type": result.input_type,
            "output_type": result.output_type,
            "resolver_config": result.resolver_config,
            "timestamp": result.timestamp.isoformat() if hasattr(
                result,
                'timestamp') else None,
            "input_gene_count": len(genes_list),
            "success_count": len(
                result.successful),
            "ambiguous_count": len(
                result.ambiguous),
            "failed_count": len(
                result.failed),
            "success_rate": result.success_rate},
        "results": {
            "successful": successful,
            "ambiguous": ambiguous,
            "failed": list(
                result.failed)},
        "ambiguity_resolutions": result.ambiguity_resolutions,
        "corrections_applied": getattr(
            result,
            'corrections_applied',
            {}) if hasattr(
            result,
            'corrections_applied') else {}}


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
                click.echo(
                    f"   Current version: Ensembl {update_info['current']['ensembl_release']}")

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
@click.option('--genome-build', '-g', default='GRCh38',
              help='Genome build (GRCh38, GRCh37, GRCm39)')
@click.option('--ambiguity-strategy', '-r', default='primary',
              type=click.Choice(['primary', 'first', 'all', 'fail']),
              help='How to handle ambiguous gene symbols')
@click.option('--parallel/--sequential', default=True,
              help='Use parallel processing for large gene lists')
@click.option('--max-workers', default=4, type=int,
              help='Maximum number of parallel workers')
@click.option('--debug', is_flag=True,
              help='Show detailed conversion information and debugging output')
@click.option('--format', type=click.Choice(['table', 'json']), default='table',
              help='Output format (table or json)')
@click.option('--config', type=click.Path(exists=True),
              help='Configuration file (YAML format)')
@click.option('--auto-correct/--no-auto-correct', default=True,
              help='Auto-correct deprecated gene symbols')
@click.option('--data-dir', '-d', default='./data',
              help='Directory to store gene databases')
def convert(
        genes,
        from_type,
        to_type,
        genome_build,
        ambiguity_strategy,
        parallel,
        max_workers,
        debug,
        format,
        config,
        auto_correct,
        data_dir):
    """
    Convert gene IDs from one system to another.

    Ambiguity resolution strategies:
    - primary: Prefer protein-coding genes (recommended)
    - first: Use first match found
    - all: Show all matches for manual review
    - fail: Treat ambiguous genes as failures
    """
    # Load configuration
    config_data = load_config(config)

    # Apply config defaults
    options = apply_config_to_options(config_data,
                                      genome_build=genome_build,
                                      ambiguity_strategy=ambiguity_strategy,
                                      parallel=parallel,
                                      max_workers=max_workers,
                                      auto_correct=auto_correct,
                                      data_dir=data_dir
                                      )

    # Update variables with config defaults
    genome_build = options['genome_build'] or 'GRCh38'
    ambiguity_strategy = options['ambiguity_strategy'] or 'primary'
    parallel = options['parallel'] if options['parallel'] is not None else True
    max_workers = options['max_workers'] or 4
    auto_correct = options['auto_correct'] if options['auto_correct'] is not None else True
    data_dir = options['data_dir'] or './data'

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
                    list(genes),
                    from_type,
                    to_type,
                    genome_build,
                    ambiguity_strategy,
                    parallel,
                    max_workers)
            except ImportError:
                click.echo("Warning: Enhanced features not available, using basic conversion")
                resolver = GeneResolver(Path(data_dir))
                result = resolver.convert(
                    list(genes),
                    from_type,
                    to_type,
                    genome_build,
                    ambiguity_strategy,
                    parallel,
                    max_workers)
        else:
            resolver = GeneResolver(Path(data_dir))
            result = resolver.convert(
                list(genes),
                from_type,
                to_type,
                genome_build,
                ambiguity_strategy,
                parallel,
                max_workers)

        # Display results
        if format == 'json':
            json_output = result_to_json(result, list(genes))
            click.echo(json.dumps(json_output, indent=2))
        else:
            # Professional table format with colors
            success_count = len(result.successful)
            ambiguous_count = len(result.ambiguous)
            failed_count = len(result.failed)

            # Format with ANSI colors (direct escape codes for guaranteed support)
            success_text = f"\033[32m✓ Successful: {success_count}\033[0m"  # Green
            ambiguous_text = f"\033[33m⚠ Ambiguous: {ambiguous_count}\033[0m"  # Yellow
            failed_text = f"\033[31m✗ Failed: {failed_count}\033[0m"  # Red

            # Calculate dynamic width based on content (excluding ANSI codes)
            plain_success = f"✓ Successful: {success_count}"
            plain_ambiguous = f"⚠ Ambiguous: {ambiguous_count}"
            plain_failed = f"✗ Failed: {failed_count}"
            content_width = len(f"{plain_success} │ {plain_ambiguous} │ {plain_failed}")
            box_width = max(59, content_width + 4)  # Minimum 59, or content + padding

            click.echo(f"\n┌─ Conversion Results {'─' * (box_width - 25)}─┐")
            click.echo(f"│ Genome: {genome_build:<10} │ Strategy: {ambiguity_strategy:<8} │")
            click.echo(f"├{'─' * box_width}┤")
            print(f"│ {success_text} │ {ambiguous_text} │ {failed_text} │")  # Use print for colors
            click.echo(f"└{'─' * box_width}┘")

            if debug:
                click.echo(f"\n\033[34m[DEBUG INFO]\033[0m")  # Blue
                click.echo(f"  Parallel: {parallel} | Workers: {max_workers} | Auto-correct: {auto_correct}")
                click.echo(f"  Config: {result.resolver_config}")

                if result.successful:
                    click.echo(f"\n\033[32m[SUCCESSFUL CONVERSIONS]\033[0m")  # Green
                    for gene, mapping in list(result.successful.items())[:5]:
                        target_id = mapping.identifiers.ensembl_id or mapping.identifiers.entrez_id
                        click.echo(f"  {gene} → \033[32m{target_id}\033[0m")  # Green
                    if len(result.successful) > 5:
                        click.echo(f"  ... and {len(result.successful) - 5} more")

                if result.ambiguous:
                    click.echo(f"\n\033[33m[AMBIGUOUS GENES]\033[0m")  # Yellow
                    for gene, mappings in list(result.ambiguous.items())[:3]:
                        ensembl_ids = [m.identifiers.ensembl_id for m in mappings if m.identifiers.ensembl_id]
                        click.echo(f"  {gene} → {len(mappings)} matches: {ensembl_ids[:3]}{'...' if len(ensembl_ids) > 3 else ''}")

                if result.failed:
                    click.echo(f"\n\033[31m[FAILED GENES]\033[0m")  # Red
                    for gene in list(result.failed)[:5]:
                        click.echo(f"  {gene} → \033[31mNo matches found\033[0m")  # Red
                    if len(result.failed) > 5:
                        click.echo(f"  ... and {len(result.failed) - 5} more")

                if hasattr(result, 'ambiguity_resolutions') and result.ambiguity_resolutions:
                    click.echo(f"\n\033[36m[AMBIGUITY RESOLUTIONS]\033[0m")  # Cyan
                    resolutions = list(result.ambiguity_resolutions.items())[:5]
                    for gene, resolution in resolutions:
                        click.echo(f"  {gene}: {resolution}")

            # Show auto-corrections if any were applied
        if auto_correct and hasattr(result, 'corrections_applied') and result.corrections_applied:
            click.echo(f"  Auto-corrected: {len(result.corrections_applied)} genes")
            for old, new in result.corrections_applied.items():
                click.echo(f"    {old} → {new}")

        # Show ambiguity resolutions
        primary_resolved = sum(1 for v in result.ambiguity_resolutions.values()
                               if v == "resolved_primary")
        if primary_resolved > 0:
            click.echo(f"  Auto-resolved: {primary_resolved} ambiguous genes")

        # Show successful conversions
        if result.successful:
            click.echo(f"\n[SUCCESSFUL CONVERSIONS]")
            for input_id, mapping in result.successful.items():
                output_id = get_output_id(mapping, to_type)
                resolution = result.ambiguity_resolutions.get(input_id, "")

                was_corrected = False
                if auto_correct and hasattr(result, 'corrections_applied'):
                    was_corrected = input_id in result.corrections_applied.values()
                    original_gene = next(
                        (k for k, v in result.corrections_applied.items() if v == input_id), None)

                if output_id:
                    if was_corrected and original_gene:
                        click.echo(
                            f"  {original_gene} → {output_id} (auto-corrected from {original_gene})")
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
                        click.echo(f"     {i + 1}. {output_id} ({biotype})")

        # Show failed conversions
        if result.failed:
            click.echo(f"\nFailed conversions:")
            for failed_id in result.failed:
                resolution = result.ambiguity_resolutions.get(failed_id, "")
                if resolution == "failed_ambiguous":
                    click.echo(
                        f"   {failed_id} (ambiguous - use '--ambiguity-strategy all' to see options)")
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
@click.option('--genome-build', '-g', default='GRCh38',
              help='Genome build (GRCh38, GRCh37, GRCm39)')
@click.option('--ambiguity-strategy', '-r', default='primary',
              type=click.Choice(['primary', 'first', 'all', 'fail']),
              help='How to handle ambiguous gene symbols')
@click.option('--parallel/--sequential', default=True,
              help='Use parallel processing for large gene lists')
@click.option('--max-workers', default=4, type=int,
              help='Maximum number of parallel workers')
@click.option('--debug', is_flag=True,
              help='Show detailed conversion information and debugging output')
@click.option('--format', type=click.Choice(['table', 'json']), default='table',
              help='Output format (table or json)')
@click.option('--config', type=click.Path(exists=True),
              help='Configuration file (YAML format)')
@click.option('--column', '-c', default=0,
              help='Column containing gene IDs (0-based, for CSV/TSV)')
@click.option('--data-dir', '-d', default='./data',
              help='Directory to store gene databases')
@click.option('--auto-correct/--no-auto-correct', default=True,
              help='Auto-correct deprecated gene symbols')
def convert_file(
        input_file,
        output_file,
        from_type,
        to_type,
        genome_build,
        ambiguity_strategy,
        parallel,
        max_workers,
        debug,
        format,
        config,
        column,
        data_dir,
        auto_correct):
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

    # Load configuration
    config_data = load_config(config)

    # Apply config defaults
    options = apply_config_to_options(config_data,
                                      genome_build=genome_build,
                                      ambiguity_strategy=ambiguity_strategy,
                                      parallel=parallel,
                                      max_workers=max_workers,
                                      auto_correct=auto_correct,
                                      data_dir=data_dir,
                                      column=column
                                      )

    # Update variables with config defaults
    genome_build = options['genome_build'] or 'GRCh38'
    ambiguity_strategy = options['ambiguity_strategy'] or 'primary'
    parallel = options['parallel'] if options['parallel'] is not None else True
    max_workers = options['max_workers'] or 4
    auto_correct = options['auto_correct'] if options['auto_correct'] is not None else True
    data_dir = options['data_dir'] or './data'
    column = options['column'] if options['column'] is not None else 0

    try:
        # Check if database exists
        db_file = Path(data_dir) / "genes.db"
        if not db_file.exists():
            click.echo("Warning: Database not found. Run 'gene-resolver init' first.")
            return

        # Read genes from file
        click.echo(f"[INPUT] Reading genes from: {input_file}")
        input_genes = read_genes_from_file(Path(input_file), column)

        if not input_genes:
            click.echo("[ERROR] No genes found in input file")
            return

        click.echo(f"[INFO] Found {len(input_genes)} genes to convert")
        click.echo(f"[INFO] Genome build: {genome_build}")
        click.echo(f"[INFO] Ambiguity strategy: {ambiguity_strategy}")
        click.echo(f"[INFO] Auto-correct: {auto_correct}")

        # Use enhanced resolver if auto-correct is enabled
        if auto_correct:
            try:
                from gene_id_resolver.core.resolver_enhanced import EnhancedGeneResolver
                resolver = EnhancedGeneResolver(Path(data_dir))
                result = resolver.convert_with_correction(
                    input_genes,
                    from_type,
                    to_type,
                    genome_build,
                    ambiguity_strategy,
                    parallel,
                    max_workers)
            except ImportError:
                click.echo("Warning: Enhanced features not available, using basic conversion")
                resolver = GeneResolver(Path(data_dir))
                result = resolver.convert(
                    input_genes,
                    from_type,
                    to_type,
                    genome_build,
                    ambiguity_strategy,
                    parallel,
                    max_workers)
        else:
            resolver = GeneResolver(Path(data_dir))
            result = resolver.convert(
                input_genes,
                from_type,
                to_type,
                genome_build,
                ambiguity_strategy,
                parallel,
                max_workers)

        # Display summary with genome info
        if format == 'json':
            json_output = result_to_json(result, input_genes)
            click.echo(json.dumps(json_output, indent=2))
        else:
            # Professional table format with colors
            success_count = len(result.successful)
            ambiguous_count = len(result.ambiguous)
            failed_count = len(result.failed)

            # Format with ANSI colors (direct escape codes for guaranteed support)
            success_text = f"\033[32m✓ Successful: {success_count}\033[0m"  # Green
            ambiguous_text = f"\033[33m⚠ Ambiguous: {ambiguous_count}\033[0m"  # Yellow
            failed_text = f"\033[31m✗ Failed: {failed_count}\033[0m"  # Red

            # Calculate dynamic width based on content (excluding ANSI codes)
            plain_success = f"✓ Successful: {success_count}"
            plain_ambiguous = f"⚠ Ambiguous: {ambiguous_count}"
            plain_failed = f"✗ Failed: {failed_count}"
            content_width = len(f"{plain_success} │ {plain_ambiguous} │ {plain_failed}")
            box_width = max(59, content_width + 4)  # Minimum 59, or content + padding

            click.echo(f"\n┌─ Conversion Results {'─' * (box_width - 25)}─┐")
            click.echo(f"│ Genome: {genome_build:<10} │ Strategy: {ambiguity_strategy:<8} │")
            click.echo(f"├{'─' * box_width}┤")
            print(f"│ {success_text} │ {ambiguous_text} │ {failed_text} │")  # Use print for colors
            click.echo(f"│ Success Rate: {result.success_rate:.1%}{' ' * (box_width - len(f'Success Rate: {result.success_rate:.1%}') - 1)} │")
            click.echo(f"└{'─' * box_width}┘")

            if debug:
                click.echo(f"\n\033[34m[DEBUG INFO]\033[0m")  # Blue
                click.echo(f"  Parallel: {parallel} | Workers: {max_workers} | Auto-correct: {auto_correct}")
                click.echo(f"  Input: {input_file} | Output: {output_file or 'None'} | Column: {column}")
                click.echo(f"  Config: {result.resolver_config}")

                if result.successful:
                    click.echo(f"\n\033[32m[SUCCESSFUL CONVERSIONS]\033[0m")  # Green
                    for gene, mapping in list(result.successful.items())[:5]:
                        target_id = mapping.identifiers.ensembl_id or mapping.identifiers.entrez_id
                        click.echo(f"  {gene} → \033[32m{target_id}\033[0m")  # Green
                    if len(result.successful) > 5:
                        click.echo(f"  ... and {len(result.successful) - 5} more")

                if result.ambiguous:
                    click.echo(f"\n\033[33m[AMBIGUOUS GENES]\033[0m")  # Yellow
                    for gene, mappings in list(result.ambiguous.items())[:3]:
                        ensembl_ids = [m.identifiers.ensembl_id for m in mappings if m.identifiers.ensembl_id]
                        click.echo(f"  {gene} → {len(mappings)} matches: {ensembl_ids[:3]}{'...' if len(ensembl_ids) > 3 else ''}")

                if result.failed:
                    click.echo(f"\n\033[31m[FAILED GENES]\033[0m")  # Red
                    for gene in list(result.failed)[:5]:
                        click.echo(f"  {gene} → \033[31mNo matches found\033[0m")  # Red
                    if len(result.failed) > 5:
                        click.echo(f"  ... and {len(result.failed) - 5} more")

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
            click.echo(f"[OUTPUT] Results saved to: {output_file}")
        else:
            # Display results in terminal
            if result.successful:
                click.echo(f"\n\033[32m[SUCCESSFUL CONVERSIONS]\033[0m")  # Green
                successful_genes = list(result.successful.keys())[:10]
                for gene in successful_genes:
                    mapping = result.successful[gene]
                    output_id = get_output_id(mapping, to_type)
                    resolution = result.ambiguity_resolutions.get(gene, "")

                    was_corrected = False
                    if auto_correct and hasattr(result, 'corrections_applied'):
                        was_corrected = gene in result.corrections_applied.values()
                        original_gene = next(
                            (k for k, v in result.corrections_applied.items() if v == gene), None)

                    if output_id:
                        if was_corrected and original_gene:
                            click.echo(
                                f"  {original_gene} → \033[32m{output_id}\033[0m (auto-corrected from {original_gene})")
                        elif resolution == "resolved_primary":
                            click.echo(f"   {gene} → \033[32m{output_id}\033[0m (auto-resolved protein-coding)")
                        else:
                            click.echo(f"   {gene} → \033[32m{output_id}\033[0m")

                if len(result.successful) > 10:
                    click.echo(f"   ... and {len(result.successful) - 10} more")

            if result.ambiguous:
                click.echo(f"\n\033[33mAmbiguous mappings (first 5):\033[0m")  # Yellow
                ambiguous_genes = list(result.ambiguous.keys())[:5]
                for gene in ambiguous_genes:
                    mappings = result.ambiguous[gene]
                    click.echo(f"   {gene} maps to {len(mappings)} locations:")
                    for i, mapping in enumerate(mappings[:3]):
                        output_id = get_output_id(mapping, to_type)
                        biotype = mapping.biotype or "unknown"
                        if output_id:
                            click.echo(f"     {i + 1}. \033[33m{output_id}\033[0m ({biotype})")  # Yellow
                    if len(mappings) > 3:
                        click.echo(f"     ... and {len(mappings) - 3} more")

            if result.failed:
                click.echo(f"\n\033[31mFailed conversions (first 10):\033[0m")  # Red
                failed_genes = result.failed[:10]
                for gene in failed_genes:
                    resolution = result.ambiguity_resolutions.get(gene, "")
                    if resolution == "failed_ambiguous":
                        click.echo(
                            f"   \033[31m{gene}\033[0m (ambiguous - use '--ambiguity-strategy all' to see options)")
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

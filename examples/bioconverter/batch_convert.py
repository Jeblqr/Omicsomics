#!/usr/bin/env python3
"""
Batch Bioconverter Example
===========================

Converts all files of a specific format in a directory.

Usage:
    python batch_convert.py <input_dir> <output_dir> [--from-format FORMAT] [--to-format FORMAT]

Example:
    python batch_convert.py vcf_files/ csv_output/
    python batch_convert.py data/ converted/ --from-format vcf --to-format csv
"""

import sys
import argparse
from pathlib import Path


def batch_convert(
    input_dir, output_dir, from_format="vcf", to_format="csv", verbose=True
):
    """
    Convert all files in input_dir to output_dir

    Args:
        input_dir: Input directory path
        output_dir: Output directory path
        from_format: Source file format
        to_format: Target file format
        verbose: Print progress information
    """

    try:
        from bioconverter import Converter
    except ImportError:
        print("Error: bioconverter is not installed", file=sys.stderr)
        print("\nInstall it with:", file=sys.stderr)
        print("  pip install bioconverter", file=sys.stderr)
        sys.exit(1)

    # Convert to Path objects
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    # Check if input directory exists
    if not input_path.exists():
        print(f"Error: Input directory '{input_dir}' not found", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Find all input files
    pattern = f"*.{from_format}"
    files = list(input_path.glob(pattern))

    if not files:
        print(f"No .{from_format} files found in {input_dir}")
        return

    if verbose:
        print(f"Found {len(files)} .{from_format} files")
        print(f"Converting to .{to_format} format...")

    # Try to use tqdm for progress bar
    try:
        from tqdm import tqdm

        files_iter = tqdm(files, desc="Converting")
    except ImportError:
        files_iter = files
        if verbose:
            print("(Install tqdm for progress bar: pip install tqdm)")

    # Convert files
    converter = Converter()
    success_count = 0
    error_count = 0

    for input_file in files_iter:
        output_file = output_path / f"{input_file.stem}.{to_format}"

        try:
            converter.convert(str(input_file), str(output_file))
            success_count += 1

            if verbose and "tqdm" not in sys.modules:
                print(f"  ✓ {input_file.name} -> {output_file.name}")

        except Exception as e:
            error_count += 1
            if verbose:
                print(f"\n  ✗ Error converting {input_file.name}: {e}")

    # Summary
    if verbose:
        print(f"\n{'='*60}")
        print(f"Conversion Summary:")
        print(f"  Success: {success_count}/{len(files)}")
        print(f"  Errors:  {error_count}/{len(files)}")
        print(f"  Output directory: {output_path.absolute()}")
        print(f"{'='*60}")


def main():
    """Main function with argument parsing"""

    parser = argparse.ArgumentParser(
        description="Batch convert bioinformatics files using bioconverter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python batch_convert.py vcf_files/ csv_output/
  python batch_convert.py data/ converted/ --from-format fastq --to-format fasta
  python batch_convert.py input/ output/ --from-format vcf --to-format json
        """,
    )

    parser.add_argument("input_dir", help="Input directory containing files to convert")

    parser.add_argument("output_dir", help="Output directory for converted files")

    parser.add_argument(
        "--from-format", default="vcf", help="Source file format (default: vcf)"
    )

    parser.add_argument(
        "--to-format", default="csv", help="Target file format (default: csv)"
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress output"
    )

    args = parser.parse_args()

    # Run batch conversion
    batch_convert(
        args.input_dir,
        args.output_dir,
        from_format=args.from_format,
        to_format=args.to_format,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()

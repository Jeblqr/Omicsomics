#!/usr/bin/env python3
"""
Validate and Convert Example
==============================

Demonstrates file validation before conversion using bioconverter.

Usage:
    python validate_and_convert.py <input_file> <output_file> [--skip-validation]

Example:
    python validate_and_convert.py input.vcf output.csv
"""

import sys
import argparse
from pathlib import Path


def validate_file(file_path, format_type=None):
    """
    Validate a file format

    Args:
        file_path: Path to file
        format_type: Expected format type (auto-detected if None)

    Returns:
        tuple: (is_valid, errors, detected_format)
    """

    try:
        from bioconverter import validate_file, detect_format
    except ImportError:
        print("Error: bioconverter is not installed", file=sys.stderr)
        print("\nInstall it with:", file=sys.stderr)
        print("  pip install bioconverter", file=sys.stderr)
        sys.exit(1)

    # Detect format if not provided
    if format_type is None:
        try:
            format_type = detect_format(file_path)
            print(f"Detected format: {format_type}")
        except Exception as e:
            print(f"Warning: Could not detect format: {e}")
            return False, [str(e)], None

    # Validate file
    try:
        is_valid, errors = validate_file(file_path, format_type=format_type)
        return is_valid, errors, format_type
    except Exception as e:
        return False, [str(e)], format_type


def convert_with_validation(input_file, output_file, skip_validation=False):
    """
    Convert file with validation

    Args:
        input_file: Input file path
        output_file: Output file path
        skip_validation: Skip validation step
    """

    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)

    try:
        from bioconverter import Converter
    except ImportError:
        print("Error: bioconverter is not installed", file=sys.stderr)
        print("\nInstall it with:", file=sys.stderr)
        print("  pip install bioconverter", file=sys.stderr)
        sys.exit(1)

    print(f"{'='*60}")
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")
    print(f"{'='*60}\n")

    # Validation step
    if not skip_validation:
        print("Step 1: Validating input file...")
        is_valid, errors, detected_format = validate_file(input_file)

        if is_valid:
            print("✓ File is valid\n")
        else:
            print("✗ Validation failed:")
            for error in errors:
                print(f"  - {error}")

            response = input("\nContinue anyway? (y/n): ")
            if response.lower() != "y":
                print("Conversion cancelled")
                sys.exit(1)
            print()
    else:
        print("Step 1: Skipping validation\n")

    # Conversion step
    print("Step 2: Converting file...")

    try:
        converter = Converter()
        converter.convert(input_file, output_file)

        print("✓ Conversion successful\n")

        # Show output file information
        output_path = Path(output_file)
        if output_path.exists():
            size_kb = output_path.stat().st_size / 1024
            print(f"Output file details:")
            print(f"  Path: {output_path.absolute()}")
            print(f"  Size: {size_kb:.2f} KB")

        print(f"\n{'='*60}")
        print("Conversion completed successfully!")
        print(f"{'='*60}")

    except Exception as e:
        print(f"✗ Conversion failed: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main function with argument parsing"""

    parser = argparse.ArgumentParser(
        description="Convert bioinformatics files with validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python validate_and_convert.py input.vcf output.csv
  python validate_and_convert.py reads.fastq reads.fasta
  python validate_and_convert.py data.gff data.bed --skip-validation
        """,
    )

    parser.add_argument("input_file", help="Input file to convert")

    parser.add_argument("output_file", help="Output file path")

    parser.add_argument(
        "--skip-validation", action="store_true", help="Skip validation step"
    )

    args = parser.parse_args()

    # Run conversion with validation
    convert_with_validation(
        args.input_file, args.output_file, skip_validation=args.skip_validation
    )


if __name__ == "__main__":
    main()

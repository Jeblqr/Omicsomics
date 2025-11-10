#!/usr/bin/env python3
"""
Simple Bioconverter Example
===========================

Demonstrates basic file format conversion using bioconverter.

Usage:
    python simple_convert.py <input_file> <output_file>

Example:
    python simple_convert.py input.vcf output.csv
"""

import sys
from pathlib import Path


def main():
    """Main conversion function"""

    # Check arguments
    if len(sys.argv) != 3:
        print("Usage: python simple_convert.py <input_file> <output_file>")
        print("\nExample:")
        print("  python simple_convert.py input.vcf output.csv")
        print("  python simple_convert.py reads.fastq reads.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)

    try:
        # Import bioconverter
        from bioconverter import Converter

        print(f"Converting {input_file} to {output_file}...")

        # Create converter instance
        converter = Converter()

        # Perform conversion
        converter.convert(input_file, output_file)

        print(f"✓ Successfully converted to {output_file}")

        # Show output file size
        output_path = Path(output_file)
        if output_path.exists():
            size_kb = output_path.stat().st_size / 1024
            print(f"  Output size: {size_kb:.2f} KB")

    except ImportError as e:
        print("Error: bioconverter is not installed", file=sys.stderr)
        print("\nInstall it with:", file=sys.stderr)
        print("  pip install bioconverter", file=sys.stderr)
        sys.exit(1)

    except Exception as e:
        print(f"Error during conversion: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Advanced Bioconverter Usage
============================

Demonstrates advanced features including:
- Format detection
- Listing available converters
- Multiple format conversions
- Error handling
- Statistics tracking

Usage:
    python advanced_usage.py
"""

import sys
from pathlib import Path
from typing import List, Dict, Any


def print_section(title: str):
    """Print a section header"""
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


def demo_format_detection():
    """Demonstrate format detection"""
    print_section("Format Detection")

    try:
        from bioconverter import detect_format
    except ImportError:
        print("Error: bioconverter is not installed")
        return

    # Example files to detect
    test_files = [
        "../../test_data/sample_genomics.vcf",
        "../../test_data/sample.fastq",
        "../../test_data/sample.fasta",
    ]

    for file_path in test_files:
        if Path(file_path).exists():
            try:
                format_type = detect_format(file_path)
                print(f"✓ {Path(file_path).name:30} -> {format_type}")
            except Exception as e:
                print(f"✗ {Path(file_path).name:30} -> Error: {e}")
        else:
            print(f"- {Path(file_path).name:30} -> File not found (example)")


def demo_list_converters():
    """Demonstrate listing available converters"""
    print_section("Available Converters")

    try:
        from bioconverter import list_converters
    except ImportError:
        print("Error: bioconverter is not installed")
        return

    # Get all converters
    try:
        converters = list_converters()

        print(f"Total converters available: {len(converters)}\n")

        # Group by source format
        from collections import defaultdict

        by_format = defaultdict(list)

        for conv in converters:
            by_format[conv["from"]].append(conv["to"])

        # Print grouped converters
        for from_format in sorted(by_format.keys()):
            to_formats = ", ".join(sorted(by_format[from_format]))
            print(f"{from_format:15} -> {to_formats}")

    except Exception as e:
        print(f"Error: {e}")


def demo_multiple_conversions():
    """Demonstrate converting one file to multiple formats"""
    print_section("Multiple Format Conversions")

    try:
        from bioconverter import Converter
    except ImportError:
        print("Error: bioconverter is not installed")
        return

    # Input file
    input_file = "../../test_data/sample_genomics.vcf"

    if not Path(input_file).exists():
        print(f"Example input file not found: {input_file}")
        print("(This would work with an actual VCF file)")
        return

    # Output formats
    output_formats = ["csv", "json", "tsv"]

    print(f"Converting {Path(input_file).name} to multiple formats:\n")

    converter = Converter()

    for fmt in output_formats:
        output_file = f"output_example.{fmt}"

        try:
            converter.convert(input_file, output_file)

            # Get file size
            size_kb = Path(output_file).stat().st_size / 1024
            print(f"  ✓ {fmt:10} -> {output_file:30} ({size_kb:.2f} KB)")

            # Clean up example file
            Path(output_file).unlink()

        except Exception as e:
            print(f"  ✗ {fmt:10} -> Error: {e}")


def demo_error_handling():
    """Demonstrate error handling"""
    print_section("Error Handling")

    try:
        from bioconverter import Converter
    except ImportError:
        print("Error: bioconverter is not installed")
        return

    converter = Converter()

    # Test cases
    test_cases = [
        {
            "input": "nonexistent_file.vcf",
            "output": "output.csv",
            "description": "Non-existent input file",
        },
        {
            "input": "/root/protected.vcf",
            "output": "output.csv",
            "description": "Permission denied (example)",
        },
    ]

    for test in test_cases:
        print(f"Test: {test['description']}")

        try:
            converter.convert(test["input"], test["output"])
            print(f"  ✓ Unexpected success\n")

        except FileNotFoundError as e:
            print(f"  ✓ Caught FileNotFoundError: {e}\n")

        except PermissionError as e:
            print(f"  ✓ Caught PermissionError: {e}\n")

        except Exception as e:
            print(f"  ✓ Caught {type(e).__name__}: {e}\n")


def demo_conversion_pipeline():
    """Demonstrate a conversion pipeline"""
    print_section("Conversion Pipeline")

    try:
        from bioconverter import Converter
    except ImportError:
        print("Error: bioconverter is not installed")
        return

    print("Pipeline: VCF -> CSV -> JSON\n")

    # Check if we have a test file
    input_file = "../../test_data/sample_genomics.vcf"

    if not Path(input_file).exists():
        print("Example input file not found")
        print("(This would work with an actual VCF file)\n")
        print("Pipeline steps:")
        print("  1. VCF -> CSV (tabular format)")
        print("  2. CSV -> JSON (structured format)")
        print("\nThis allows multi-step transformations!")
        return

    converter = Converter()

    try:
        # Step 1: VCF to CSV
        print("Step 1: Converting VCF to CSV...")
        csv_file = "temp_output.csv"
        converter.convert(input_file, csv_file)
        print(f"  ✓ Created {csv_file}\n")

        # Step 2: CSV to JSON
        print("Step 2: Converting CSV to JSON...")
        json_file = "temp_output.json"
        converter.convert(csv_file, json_file)
        print(f"  ✓ Created {json_file}\n")

        print("Pipeline completed successfully!")

        # Clean up
        for temp_file in [csv_file, json_file]:
            if Path(temp_file).exists():
                Path(temp_file).unlink()

    except Exception as e:
        print(f"  ✗ Pipeline failed: {e}")


def demo_statistics():
    """Demonstrate tracking conversion statistics"""
    print_section("Conversion Statistics")

    try:
        from bioconverter import Converter
    except ImportError:
        print("Error: bioconverter is not installed")
        return

    print("Tracking multiple conversions:\n")

    # Simulate multiple conversions
    conversions = [
        ("file1.vcf", "file1.csv"),
        ("file2.vcf", "file2.csv"),
        ("file3.vcf", "file3.csv"),
    ]

    converter = Converter()
    stats = {"total": len(conversions), "success": 0, "failed": 0, "errors": []}

    for input_file, output_file in conversions:
        try:
            # This would fail with non-existent files (for demo)
            converter.convert(input_file, output_file)
            stats["success"] += 1
        except Exception as e:
            stats["failed"] += 1
            stats["errors"].append((input_file, str(e)))

    # Print statistics
    print(f"Total conversions attempted: {stats['total']}")
    print(f"Successful: {stats['success']}")
    print(f"Failed: {stats['failed']}")

    if stats["errors"]:
        print(f"\nErrors:")
        for file_name, error in stats["errors"]:
            print(f"  - {file_name}: {error[:50]}...")


def main():
    """Run all demonstrations"""

    print(
        """
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║              Bioconverter - Advanced Usage Examples                 ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
    """
    )

    # Check if bioconverter is installed
    try:
        import bioconverter

        print(f"Bioconverter version: {bioconverter.__version__}\n")
    except ImportError:
        print("\n⚠ Bioconverter is not installed!\n")
        print("Install it with:")
        print("  pip install bioconverter\n")
        sys.exit(1)

    # Run demonstrations
    demos = [
        ("Format Detection", demo_format_detection),
        ("List Converters", demo_list_converters),
        ("Multiple Conversions", demo_multiple_conversions),
        ("Error Handling", demo_error_handling),
        ("Conversion Pipeline", demo_conversion_pipeline),
        ("Statistics", demo_statistics),
    ]

    for name, demo_func in demos:
        try:
            demo_func()
        except Exception as e:
            print(f"\nError in {name} demo: {e}")

    print_section("All Demonstrations Complete")
    print("See individual demo functions for detailed examples.\n")


if __name__ == "__main__":
    main()

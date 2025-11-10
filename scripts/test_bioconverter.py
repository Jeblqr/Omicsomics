#!/usr/bin/env python3
"""
Bioconverter Integration Test
==============================

Tests bioconverter functionality with the Omicsomics test data.

Usage:
    python test_bioconverter.py
"""

import sys
from pathlib import Path


# Colors for terminal output
class Colors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"


def print_header(text):
    """Print a section header"""
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*70}")
    print(f"  {text}")
    print(f"{'='*70}{Colors.ENDC}\n")


def print_success(text):
    """Print success message"""
    print(f"{Colors.OKGREEN}✓ {text}{Colors.ENDC}")


def print_error(text):
    """Print error message"""
    print(f"{Colors.FAIL}✗ {text}{Colors.ENDC}")


def print_info(text):
    """Print info message"""
    print(f"{Colors.OKCYAN}ℹ {text}{Colors.ENDC}")


def test_installation():
    """Test if bioconverter is installed"""
    print_header("Test 1: Check Bioconverter Installation")

    try:
        import bioconverter

        version = bioconverter.__version__
        print_success(f"Bioconverter is installed (version {version})")
        return True
    except ImportError:
        print_error("Bioconverter is not installed")
        print_info("Install with: pip install bioconverter")
        return False


def test_import_modules():
    """Test importing bioconverter modules"""
    print_header("Test 2: Import Bioconverter Modules")

    try:
        from bioconverter import (
            Converter,
            list_converters,
            detect_format,
            validate_file,
        )

        print_success("Successfully imported all modules")
        return True
    except ImportError as e:
        print_error(f"Failed to import modules: {e}")
        return False


def test_list_converters():
    """Test listing available converters"""
    print_header("Test 3: List Available Converters")

    try:
        from bioconverter import list_converters

        converters = list_converters()
        print_success(f"Found {len(converters)} converters")

        # Group by format
        from collections import defaultdict

        by_format = defaultdict(list)

        for conv in converters[:10]:  # Show first 10
            by_format[conv["from"]].append(conv["to"])

        print_info("Sample converters:")
        for from_format in sorted(by_format.keys()):
            to_formats = ", ".join(sorted(by_format[from_format]))
            print(f"  {from_format:15} -> {to_formats}")

        return True
    except Exception as e:
        print_error(f"Failed to list converters: {e}")
        return False


def test_format_detection():
    """Test format detection"""
    print_header("Test 4: Format Detection")

    try:
        from bioconverter import detect_format

        # Test files in test_data directory
        test_files = [
            "../test_data/sample_genomics.vcf",
            "../test_data/sample.fastq",
            "../test_data/sample.fasta",
        ]

        detected_count = 0
        for file_path in test_files:
            if Path(file_path).exists():
                try:
                    format_type = detect_format(file_path)
                    print_success(f"{Path(file_path).name:30} -> {format_type}")
                    detected_count += 1
                except Exception as e:
                    print_error(f"{Path(file_path).name:30} -> {e}")
            else:
                print_info(f"{Path(file_path).name:30} -> File not found (skipped)")

        return detected_count > 0

    except Exception as e:
        print_error(f"Format detection failed: {e}")
        return False


def test_vcf_conversion():
    """Test VCF to CSV conversion"""
    print_header("Test 5: VCF to CSV Conversion")

    try:
        from bioconverter import Converter

        input_file = "../test_data/sample_genomics.vcf"
        output_file = "test_output.csv"

        if not Path(input_file).exists():
            print_info(f"Test file not found: {input_file}")
            print_info("Creating a minimal test VCF...")

            # Create minimal VCF for testing
            minimal_vcf = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t50\tPASS\tDP=100
chr1\t200\t.\tT\tC\t60\tPASS\tDP=120
"""
            input_file = "minimal_test.vcf"
            with open(input_file, "w") as f:
                f.write(minimal_vcf)

            print_success("Created minimal test VCF")

        # Convert
        converter = Converter()
        converter.convert(input_file, output_file)

        # Check output
        if Path(output_file).exists():
            size_kb = Path(output_file).stat().st_size / 1024
            print_success(f"Conversion successful: {output_file} ({size_kb:.2f} KB)")

            # Show first few lines
            with open(output_file, "r") as f:
                lines = f.readlines()[:5]

            print_info(f"First {len(lines)} lines of output:")
            for line in lines:
                print(f"    {line.rstrip()}")

            # Clean up
            Path(output_file).unlink()
            if input_file == "minimal_test.vcf":
                Path(input_file).unlink()

            return True
        else:
            print_error("Output file not created")
            return False

    except Exception as e:
        print_error(f"Conversion failed: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_validation():
    """Test file validation"""
    print_header("Test 6: File Validation")

    try:
        from bioconverter import validate_file

        # Create a test VCF
        test_vcf = "test_validation.vcf"
        with open(test_vcf, "w") as f:
            f.write(
                """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t50\tPASS\tDP=100
"""
            )

        is_valid, errors = validate_file(test_vcf, format_type="vcf")

        if is_valid:
            print_success("File validation passed")
        else:
            print_error("File validation failed:")
            for error in errors:
                print(f"  - {error}")

        # Clean up
        Path(test_vcf).unlink()

        return True

    except Exception as e:
        print_error(f"Validation test failed: {e}")
        return False


def test_batch_conversion():
    """Test batch conversion"""
    print_header("Test 7: Batch Conversion")

    try:
        from bioconverter import Converter

        # Create test files
        test_dir = Path("test_batch_input")
        output_dir = Path("test_batch_output")

        test_dir.mkdir(exist_ok=True)
        output_dir.mkdir(exist_ok=True)

        # Create 3 test VCF files
        for i in range(1, 4):
            vcf_file = test_dir / f"test{i}.vcf"
            with open(vcf_file, "w") as f:
                f.write(
                    f"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t{i}00\t.\tA\tG\t50\tPASS\tDP={i}00
"""
                )

        print_info(f"Created 3 test VCF files in {test_dir}")

        # Batch convert
        converter = Converter()
        files = list(test_dir.glob("*.vcf"))

        success_count = 0
        for vcf_file in files:
            csv_file = output_dir / f"{vcf_file.stem}.csv"
            try:
                converter.convert(str(vcf_file), str(csv_file))
                success_count += 1
                print_success(f"Converted {vcf_file.name} -> {csv_file.name}")
            except Exception as e:
                print_error(f"Failed to convert {vcf_file.name}: {e}")

        print_success(f"Batch conversion: {success_count}/{len(files)} successful")

        # Clean up
        import shutil

        shutil.rmtree(test_dir)
        shutil.rmtree(output_dir)

        return success_count == len(files)

    except Exception as e:
        print_error(f"Batch conversion test failed: {e}")
        return False


def test_error_handling():
    """Test error handling"""
    print_header("Test 8: Error Handling")

    try:
        from bioconverter import Converter

        converter = Converter()

        # Test 1: Non-existent file
        try:
            converter.convert("nonexistent.vcf", "output.csv")
            print_error("Should have raised FileNotFoundError")
            return False
        except FileNotFoundError:
            print_success("Correctly handled FileNotFoundError")
        except Exception as e:
            print_info(f"Raised {type(e).__name__}: {e}")

        # Test 2: Invalid format
        invalid_file = "invalid.txt"
        with open(invalid_file, "w") as f:
            f.write("This is not a valid VCF file\n")

        try:
            converter.convert(invalid_file, "output.csv", from_format="vcf")
            print_info("Conversion attempted with invalid file")
        except Exception as e:
            print_success(f"Correctly handled invalid format: {type(e).__name__}")
        finally:
            Path(invalid_file).unlink()

        return True

    except Exception as e:
        print_error(f"Error handling test failed: {e}")
        return False


def main():
    """Run all tests"""

    print(
        f"""
{Colors.HEADER}{Colors.BOLD}╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║              Bioconverter Integration Test Suite                    ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝{Colors.ENDC}
"""
    )

    tests = [
        ("Installation Check", test_installation),
        ("Module Import", test_import_modules),
        ("List Converters", test_list_converters),
        ("Format Detection", test_format_detection),
        ("VCF Conversion", test_vcf_conversion),
        ("File Validation", test_validation),
        ("Batch Conversion", test_batch_conversion),
        ("Error Handling", test_error_handling),
    ]

    results = []

    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print_error(f"Test '{name}' crashed: {e}")
            results.append((name, False))

    # Summary
    print_header("Test Summary")

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for name, result in results:
        status = (
            f"{Colors.OKGREEN}PASS{Colors.ENDC}"
            if result
            else f"{Colors.FAIL}FAIL{Colors.ENDC}"
        )
        print(f"  {name:30} [{status}]")

    print(f"\n{Colors.BOLD}Total: {passed}/{total} tests passed{Colors.ENDC}")

    if passed == total:
        print(f"{Colors.OKGREEN}{Colors.BOLD}\n✓ All tests passed!{Colors.ENDC}\n")
        return 0
    else:
        print(f"{Colors.WARNING}{Colors.BOLD}\n⚠ Some tests failed{Colors.ENDC}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())

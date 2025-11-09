#!/usr/bin/env python3
"""
Test script for proteomics and metabolomics converters
"""
import asyncio
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from app.converters.base import ConverterFactory
from app.schemas.unified_format import OmicsType


async def test_proteomics_converter():
    """Test proteomics converter with sample data."""
    print("=" * 60)
    print("Testing Proteomics Converter")
    print("=" * 60)

    try:
        converter = ConverterFactory.create(OmicsType.PROTEOMICS)
        print(f"✓ Proteomics converter created: {converter.__class__.__name__}")
        print(f"  Omics type: {converter.omics_type}")
        print(f"  Supported formats: {converter.supported_inputs}")

        # Test with a sample MGF format (text-based, easy to create)
        test_mgf_content = """BEGIN IONS
TITLE=Test Spectrum 1
RTINSECONDS=100.5
PEPMASS=500.25
CHARGE=2+
147.113 1000.0
175.119 2000.0
END IONS

BEGIN IONS
TITLE=Test Spectrum 2
RTINSECONDS=150.8
PEPMASS=600.30
CHARGE=3+
200.100 1500.0
250.150 2500.0
300.200 3000.0
END IONS
"""

        # Create temporary test file
        test_file = Path("/tmp/test_proteomics.mgf")
        test_file.write_text(test_mgf_content)

        # Test conversion to unified
        unified_data = await converter.to_unified(
            file_path=test_file,
            sample_id="test_proteomics_sample",
            source_format="mgf",
            organism="Homo sapiens",
        )

        print(f"\n✓ Converted MGF to unified format:")
        print(f"  Sample ID: {unified_data.metadata.sample_id}")
        print(f"  Omics type: {unified_data.metadata.omics_type}")
        print(f"  Source format: {unified_data.metadata.source_format}")
        print(f"  Total spectra: {len(unified_data.records)}")
        print(f"  Headers: {unified_data.headers}")
        print(
            f"  First record: {unified_data.records[0].id if unified_data.records else 'N/A'}"
        )

        # Test conversion from unified to CSV
        output_file = Path("/tmp/test_proteomics_output.csv")
        result_path = await converter.from_unified(
            unified_data=unified_data, target_format="csv", output_path=output_file
        )

        print(f"\n✓ Converted unified to CSV:")
        print(f"  Output file: {result_path}")
        print(f"  File exists: {result_path.exists()}")
        if result_path.exists():
            content_preview = result_path.read_text()[:200]
            print(f"  Content preview: {content_preview}...")

        # Cleanup
        test_file.unlink()
        if result_path.exists():
            result_path.unlink()

        print("\n✓ Proteomics converter test PASSED\n")
        return True

    except Exception as e:
        print(f"\n✗ Proteomics converter test FAILED: {e}\n")
        import traceback

        traceback.print_exc()
        return False


async def test_metabolomics_converter():
    """Test metabolomics converter with sample data."""
    print("=" * 60)
    print("Testing Metabolomics Converter")
    print("=" * 60)

    try:
        converter = ConverterFactory.create(OmicsType.METABOLOMICS)
        print(f"✓ Metabolomics converter created: {converter.__class__.__name__}")
        print(f"  Omics type: {converter.omics_type}")
        print(f"  Supported formats: {converter.supported_inputs}")

        # Test with a sample CSV format
        test_csv_content = """metabolite_id,compound_name,mz,retention_time,intensity,annotation
M001,Glucose,180.063,5.2,1000000,Hexose
M002,Lactate,90.032,3.5,500000,Organic acid
M003,Alanine,89.047,2.8,750000,Amino acid
M004,Citrate,192.027,4.1,1200000,TCA cycle
"""

        # Create temporary test file
        test_file = Path("/tmp/test_metabolomics.csv")
        test_file.write_text(test_csv_content)

        # Test conversion to unified
        unified_data = await converter.to_unified(
            file_path=test_file,
            sample_id="test_metabolomics_sample",
            source_format="csv",
            organism="Homo sapiens",
        )

        print(f"\n✓ Converted CSV to unified format:")
        print(f"  Sample ID: {unified_data.metadata.sample_id}")
        print(f"  Omics type: {unified_data.metadata.omics_type}")
        print(f"  Source format: {unified_data.metadata.source_format}")
        print(f"  Total metabolites: {len(unified_data.records)}")
        print(f"  Headers: {unified_data.headers}")
        print(
            f"  First record: {unified_data.records[0].id if unified_data.records else 'N/A'}"
        )
        if unified_data.records:
            print(f"  First record data: {unified_data.records[0].values}")

        # Test conversion from unified to TSV
        output_file = Path("/tmp/test_metabolomics_output.tsv")
        result_path = await converter.from_unified(
            unified_data=unified_data, target_format="tsv", output_path=output_file
        )

        print(f"\n✓ Converted unified to TSV:")
        print(f"  Output file: {result_path}")
        print(f"  File exists: {result_path.exists()}")
        if result_path.exists():
            content_preview = result_path.read_text()[:200]
            print(f"  Content preview: {content_preview}...")

        # Cleanup
        test_file.unlink()
        if result_path.exists():
            result_path.unlink()

        print("\n✓ Metabolomics converter test PASSED\n")
        return True

    except Exception as e:
        print(f"\n✗ Metabolomics converter test FAILED: {e}\n")
        import traceback

        traceback.print_exc()
        return False


async def main():
    """Run all converter tests."""
    print("\n" + "=" * 60)
    print("PROTEOMICS & METABOLOMICS CONVERTER TESTS")
    print("=" * 60 + "\n")

    results = []

    # Test proteomics converter
    results.append(await test_proteomics_converter())

    # Test metabolomics converter
    results.append(await test_metabolomics_converter())

    # Summary
    print("=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    print(f"Total tests: {len(results)}")
    print(f"Passed: {sum(results)}")
    print(f"Failed: {len(results) - sum(results)}")

    if all(results):
        print("\n✓ All tests PASSED!")
        return 0
    else:
        print("\n✗ Some tests FAILED!")
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

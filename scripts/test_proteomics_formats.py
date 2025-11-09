#!/usr/bin/env python3
"""
Test script for advanced proteomics format parsing.

Tests RAW file metadata extraction and enhanced mzML/mzXML parsing with pyteomics.
"""

import asyncio
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from app.converters.proteomics import ProteomicsConverter


async def test_raw_metadata():
    """Test RAW file metadata extraction."""
    print("\n" + "=" * 80)
    print("TEST 1: RAW File Metadata Extraction")
    print("=" * 80)

    converter = ProteomicsConverter()

    # Create a dummy RAW file for testing
    test_raw = Path("/tmp/test_sample.raw")
    test_raw.write_bytes(b"MOCK_RAW_FILE_DATA" * 1000)  # 18KB mock file

    try:
        result = await converter.to_unified(
            file_path=test_raw,
            sample_id="test_001",
            source_format="raw",
            organism="Homo sapiens",
        )

        print(f"\n✓ Successfully extracted RAW file metadata")
        print(f"  - File: {result.metadata.custom_fields.get('file_name')}")
        print(f"  - Size: {result.metadata.custom_fields.get('file_size_mb')} MB")
        print(f"  - Records: {len(result.records)}")
        print(f"\nMetadata records:")
        for record in result.records:
            print(f"  • {record.values['property']}: {record.values['value']}")

        print(f"\nConversion command:")
        print(f"  {result.metadata.custom_fields.get('conversion_command')}")

        return True
    except Exception as e:
        print(f"\n✗ RAW metadata extraction failed: {e}")
        return False
    finally:
        test_raw.unlink()


async def test_mzml_parsing():
    """Test mzML parsing with pyteomics vs XML fallback."""
    print("\n" + "=" * 80)
    print("TEST 2: mzML Parsing (Pyteomics vs XML Fallback)")
    print("=" * 80)

    converter = ProteomicsConverter()

    # Create a minimal valid mzML file
    test_mzml = Path("/tmp/test_sample.mzML")
    mzml_content = """<?xml version="1.0" encoding="UTF-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1">
  <cvList count="2">
    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"/>
    <cv id="UO" fullName="Unit Ontology"/>
  </cvList>
  <run id="test_run">
    <spectrumList count="2">
      <spectrum id="scan=1" index="0" defaultArrayLength="5">
        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="120.5" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>
        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="1234.5"/>
        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="524.3"/>
        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="98765.4"/>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="0">
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
            <binary></binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
            <binary></binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
      <spectrum id="scan=2" index="1" defaultArrayLength="10">
        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>
        <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="121.2" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>
        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="2345.6"/>
        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="262.15"/>
        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="45678.9"/>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="0">
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>
            <binary></binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>
            <binary></binary>
          </binaryDataArray>
        </binaryDataArrayList>
        <precursorList count="1">
          <precursor>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="524.3"/>
                <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2"/>
              </selectedIon>
            </selectedIonList>
          </precursor>
        </precursorList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>"""

    test_mzml.write_text(mzml_content)

    try:
        result = await converter.to_unified(
            file_path=test_mzml,
            sample_id="test_002",
            source_format="mzml",
            organism="Homo sapiens",
        )

        parser_used = result.metadata.custom_fields.get("parser", "unknown")
        print(f"\n✓ Successfully parsed mzML file")
        print(f"  - Parser: {parser_used}")
        print(f"  - Total spectra: {result.statistics.get('total_spectra', 0)}")
        print(f"  - MS1 spectra: {result.statistics.get('ms1_spectra', 0)}")
        print(f"  - MS2 spectra: {result.statistics.get('ms2_spectra', 0)}")

        if result.records:
            print(f"\nFirst spectrum details:")
            first = result.records[0]
            for key, value in first.values.items():
                if value:  # Only show non-empty values
                    print(f"  • {key}: {value}")

        if parser_used == "pyteomics":
            print("\n🎉 Pyteomics is available - using enhanced parsing!")
        else:
            print("\n⚠️  Using XML fallback parser (pyteomics not available)")
            print("   Install with: pip install pyteomics>=4.6.0")

        return True
    except Exception as e:
        print(f"\n✗ mzML parsing failed: {e}")
        import traceback

        traceback.print_exc()
        return False
    finally:
        test_mzml.unlink()


async def test_mgf_parsing():
    """Test MGF format parsing."""
    print("\n" + "=" * 80)
    print("TEST 3: MGF (Mascot Generic Format) Parsing")
    print("=" * 80)

    converter = ProteomicsConverter()

    # Create a test MGF file
    test_mgf = Path("/tmp/test_sample.mgf")
    mgf_content = """BEGIN IONS
TITLE=Spectrum 1
RTINSECONDS=120.5
PEPMASS=524.3 1234.5
CHARGE=2+
100.0 500.0
200.0 1000.0
524.3 1234.5
END IONS

BEGIN IONS
TITLE=Spectrum 2
RTINSECONDS=121.2
PEPMASS=262.15 2345.6
CHARGE=2+
150.0 800.0
262.15 2345.6
400.0 600.0
END IONS
"""

    test_mgf.write_text(mgf_content)

    try:
        result = await converter.to_unified(
            file_path=test_mgf,
            sample_id="test_003",
            source_format="mgf",
            organism="Homo sapiens",
        )

        print(f"\n✓ Successfully parsed MGF file")
        print(f"  - Total spectra: {result.statistics.get('total_spectra', 0)}")
        print(f"  - Records: {len(result.records)}")

        if result.records:
            print(f"\nFirst spectrum:")
            first = result.records[0]
            for key, value in first.values.items():
                print(f"  • {key}: {value}")

        return True
    except Exception as e:
        print(f"\n✗ MGF parsing failed: {e}")
        return False
    finally:
        test_mgf.unlink()


async def test_format_comparison():
    """Compare different formats."""
    print("\n" + "=" * 80)
    print("TEST 4: Format Feature Comparison")
    print("=" * 80)

    formats = {
        "RAW": {
            "full_spectrum_data": False,
            "metadata_extraction": True,
            "requires_conversion": True,
            "recommended_converter": "ThermoRawFileParser",
        },
        "mzML": {
            "full_spectrum_data": True,
            "metadata_extraction": True,
            "requires_conversion": False,
            "parser": "pyteomics (preferred) / XML fallback",
        },
        "mzXML": {
            "full_spectrum_data": True,
            "metadata_extraction": True,
            "requires_conversion": False,
            "parser": "pyteomics (preferred) / XML fallback",
        },
        "MGF": {
            "full_spectrum_data": True,
            "metadata_extraction": True,
            "requires_conversion": False,
            "parser": "Text-based",
        },
        "CSV/TSV": {
            "full_spectrum_data": False,
            "metadata_extraction": True,
            "requires_conversion": False,
            "parser": "Pandas/CSV",
        },
    }

    print("\nFormat Support Matrix:\n")
    print(
        f"{'Format':<12} {'Spectra':<10} {'Metadata':<10} {'Direct Parse':<13} {'Parser/Tool'}"
    )
    print("-" * 80)

    for fmt, features in formats.items():
        spectra = "✓" if features.get("full_spectrum_data") else "Metadata"
        metadata = "✓" if features.get("metadata_extraction") else "✗"
        direct = "✗ Convert" if features.get("requires_conversion") else "✓ Yes"
        parser = features.get("parser") or features.get("recommended_converter", "")

        print(f"{fmt:<12} {spectra:<10} {metadata:<10} {direct:<13} {parser}")

    return True


async def main():
    """Run all proteomics format tests."""
    print("\n" + "=" * 80)
    print("ADVANCED PROTEOMICS FORMAT PARSING TESTS")
    print("=" * 80)
    print("\nTesting enhanced proteomics format support:")
    print("  - RAW file metadata extraction")
    print("  - mzML/mzXML parsing with pyteomics")
    print("  - MGF format support")
    print("  - Format comparison")

    results = []

    results.append(await test_raw_metadata())
    results.append(await test_mzml_parsing())
    results.append(await test_mgf_parsing())
    results.append(await test_format_comparison())

    # Summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)

    passed = sum(results)
    total = len(results)

    print(f"\nTests passed: {passed}/{total}")

    if passed == total:
        print("\n✅ All tests passed!")
    else:
        print(f"\n⚠️  {total - passed} test(s) failed")

    print("\n" + "=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print(
        """
1. Install pyteomics for enhanced parsing:
   pip install pyteomics>=4.6.0

2. Install ThermoRawFileParser for RAW conversion:
   wget https://github.com/compomics/ThermoRawFileParser/releases/latest/download/ThermoRawFileParser.zip
   unzip ThermoRawFileParser.zip

3. Convert RAW files to mzML:
   ThermoRawFileParser -i sample.raw -o output_dir -f 2

4. Upload mzML files via API:
   curl -X POST http://localhost:8000/api/v1/data/upload \\
     -H "Authorization: Bearer $TOKEN" \\
     -F "file=@sample.mzML" \\
     -F "project_id=1" \\
     -F "omics_type=proteomics" \\
     -F "process_file=true"

See docs/PROTEOMICS_FORMATS.md for detailed documentation.
"""
    )

    return 0 if passed == total else 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)

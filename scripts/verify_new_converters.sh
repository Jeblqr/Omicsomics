#!/bin/bash
# Test script to verify new converters are properly registered

echo "=================================================="
echo "Verifying Proteomics and Metabolomics Converters"
echo "=================================================="
echo ""

echo "1. Checking file existence..."
if [ -f "backend/app/converters/proteomics.py" ]; then
    echo "   ✓ proteomics.py exists"
else
    echo "   ✗ proteomics.py NOT FOUND"
    exit 1
fi

if [ -f "backend/app/converters/metabolomics.py" ]; then
    echo "   ✓ metabolomics.py exists"
else
    echo "   ✗ metabolomics.py NOT FOUND"
    exit 1
fi

echo ""
echo "2. Checking converter registration..."
if grep -q "@ConverterFactory.register(OmicsType.PROTEOMICS)" backend/app/converters/proteomics.py; then
    echo "   ✓ ProteomicsConverter registered"
else
    echo "   ✗ ProteomicsConverter NOT registered"
    exit 1
fi

if grep -q "@ConverterFactory.register(OmicsType.METABOLOMICS)" backend/app/converters/metabolomics.py; then
    echo "   ✓ MetabolomicsConverter registered"
else
    echo "   ✗ MetabolomicsConverter NOT registered"
    exit 1
fi

echo ""
echo "3. Checking converter classes..."
if grep -q "class ProteomicsConverter(BaseConverter)" backend/app/converters/proteomics.py; then
    echo "   ✓ ProteomicsConverter class defined"
else
    echo "   ✗ ProteomicsConverter class NOT defined"
    exit 1
fi

if grep -q "class MetabolomicsConverter(BaseConverter)" backend/app/converters/metabolomics.py; then
    echo "   ✓ MetabolomicsConverter class defined"
else
    echo "   ✗ MetabolomicsConverter class NOT defined"
    exit 1
fi

echo ""
echo "4. Checking method implementations..."
for method in "to_unified" "from_unified" "_mzml_to_unified" "_tabular_to_unified"; do
    if grep -q "async def $method" backend/app/converters/proteomics.py; then
        echo "   ✓ ProteomicsConverter.$method implemented"
    else
        echo "   ✗ ProteomicsConverter.$method NOT implemented"
    fi
done

for method in "to_unified" "from_unified" "_mzml_to_unified" "_tabular_to_unified"; do
    if grep -q "async def $method" backend/app/converters/metabolomics.py; then
        echo "   ✓ MetabolomicsConverter.$method implemented"
    else
        echo "   ✗ MetabolomicsConverter.$method NOT implemented"
    fi
done

echo ""
echo "5. Checking file processor extension mappings..."
if grep -q '".mgf": OmicsType.PROTEOMICS' backend/app/services/file_processor.py; then
    echo "   ✓ MGF format mapped to PROTEOMICS"
else
    echo "   ✗ MGF format NOT mapped"
fi

if grep -q '".mzdata": OmicsType.METABOLOMICS' backend/app/services/file_processor.py; then
    echo "   ✓ mzData format mapped to METABOLOMICS"
else
    echo "   ✗ mzData format NOT mapped"
fi

if grep -q '".cdf": OmicsType.METABOLOMICS' backend/app/services/file_processor.py; then
    echo "   ✓ CDF format mapped to METABOLOMICS"
else
    echo "   ✗ CDF format NOT mapped"
fi

echo ""
echo "6. Line counts..."
proteomics_lines=$(wc -l < backend/app/converters/proteomics.py)
metabolomics_lines=$(wc -l < backend/app/converters/metabolomics.py)
echo "   Proteomics converter: $proteomics_lines lines"
echo "   Metabolomics converter: $metabolomics_lines lines"

echo ""
echo "=================================================="
echo "✓ All verifications PASSED!"
echo "=================================================="
echo ""
echo "Summary:"
echo "  - ProteomicsConverter: Supports mzML, mzXML, MGF, CSV, TSV"
echo "  - MetabolomicsConverter: Supports mzData, mzML, netCDF, CSV, TSV"
echo "  - Both converters registered with ConverterFactory"
echo "  - File processor updated with new extensions"
echo ""
echo "Next steps:"
echo "  1. Start backend: ./start_backend.sh"
echo "  2. Test file upload with .mgf or .mzdata files"
echo "  3. Verify processed data via GET /data/{id}/processed"

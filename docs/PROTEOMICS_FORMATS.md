# Advanced Proteomics Format Support

This document describes the advanced proteomics data format support in Omicsomics, including binary peak parsing for Thermo RAW and other vendor formats.

## Supported Formats

### Input Formats

| Format | Extension | Description | Parser | Status |
|--------|-----------|-------------|--------|--------|
| mzML | `.mzml` | Open standard MS format | pyteomics / XML fallback | ✅ Full |
| mzXML | `.mzxml` | Legacy open format | pyteomics / XML fallback | ✅ Full |
| MGF | `.mgf` | Mascot Generic Format | Text parser | ✅ Full |
| Thermo RAW | `.raw` | Thermo Fisher vendor format | Metadata extraction + conversion | ✅ Enhanced |
| CSV/TSV | `.csv`, `.tsv` | Tabular data (search results) | Pandas | ✅ Full |

### Output Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| CSV | Comma-separated values | Spreadsheet analysis |
| TSV | Tab-separated values | Unix pipelines |
| MGF | Mascot Generic Format | MS/MS search engines |
| JSON | Unified format | API integration, ML pipelines |

## Thermo RAW File Support

### Overview

Thermo RAW files are proprietary binary formats from Thermo Fisher Scientific. Full parsing requires platform-specific libraries or conversion tools.

### Conversion Workflow

**Step 1: Convert RAW to mzML**

Use ThermoRawFileParser (recommended) or msconvert:

```bash
# ThermoRawFileParser (cross-platform, Mono-based)
ThermoRawFileParser -i sample.raw -o output_dir -f 2

# MSConvert (Windows only, ProteoWizard)
msconvert sample.raw --mzML --filter "peakPicking true"
```

**Step 2: Process mzML in Omicsomics**

```python
from app.converters.proteomics import ProteomicsConverter

converter = ProteomicsConverter()
unified_data = await converter.to_unified(
    file_path=Path("sample.mzML"),
    sample_id="sample_001",
    source_format="mzML",
    organism="Homo sapiens"
)
```

### API Usage

**Automatic Conversion (Recommended)**

The `/api/v1/proteomics/convert-raw` endpoint handles conversion automatically:

```bash
curl -X POST "http://localhost:8000/api/v1/proteomics/convert-raw" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": 1,
    "raw_files": ["/path/to/sample.raw"],
    "output_dir": "/path/to/output"
  }'
```

**Direct Upload**

For pre-converted mzML files:

```bash
curl -X POST "http://localhost:8000/api/v1/data/upload" \
  -H "Authorization: Bearer $TOKEN" \
  -F "file=@sample.mzML" \
  -F "project_id=1" \
  -F "omics_type=proteomics" \
  -F "process_file=true"
```

### RAW File Metadata Extraction

When uploading a `.raw` file directly, the converter extracts basic metadata:

```python
{
    "metadata": {
        "file_name": "sample.raw",
        "file_size_mb": 245.3,
        "parsing_note": "RAW file metadata only. For full spectrum data, convert to mzML.",
        "conversion_command": "ThermoRawFileParser -i sample.raw -o output_dir -f 2"
    },
    "records": [
        {"property": "File Name", "value": "sample.raw"},
        {"property": "File Size (MB)", "value": "245.3"},
        {"property": "Format", "value": "Thermo RAW"},
        {"property": "Conversion Required", "value": "Yes - Use ThermoRawFileParser"}
    ]
}
```

## Enhanced Parsing with Pyteomics

### Installation

Pyteomics is included in backend dependencies:

```bash
pip install pyteomics>=4.6.0
```

### Features

**Better mzML/mzXML Parsing**:
- Robust handling of different vendor formats
- Accurate extraction of precursor m/z and charge
- Peak data arrays (m/z, intensity)
- Comprehensive metadata extraction

**Enhanced Data Extraction**:
```python
{
    "spectrum_id": "controllerType=0 controllerNumber=1 scan=1",
    "ms_level": "2",
    "scan_number": "1",
    "retention_time_seconds": "123.45",
    "precursor_mz": "524.2631",
    "precursor_charge": "2",
    "num_peaks": "456",
    "base_peak_mz": "524.2651",
    "base_peak_intensity": "12345.6",
    "total_ion_current": "987654.3"
}
```

### Fallback Behavior

If pyteomics is not available, the system automatically falls back to basic XML parsing:

```python
if PYTEOMICS_AVAILABLE:
    return await self._mzml_to_unified_pyteomics(...)
else:
    return await self._mzml_to_unified(...)  # XML fallback
```

## Complete Workflow Example

### 1. Convert RAW to mzML

```bash
# Convert Thermo RAW file
ThermoRawFileParser \
    -i /data/raw/sample_001.raw \
    -o /data/mzml/ \
    -f 2 \
    --metadata \
    --gzip
```

### 2. Upload and Process

```python
import httpx

async with httpx.AsyncClient() as client:
    # Upload mzML file
    with open("/data/mzml/sample_001.mzML", "rb") as f:
        response = await client.post(
            "http://localhost:8000/api/v1/data/upload",
            headers={"Authorization": f"Bearer {token}"},
            files={"file": f},
            data={
                "project_id": 1,
                "omics_type": "proteomics",
                "process_file": True,
                "async_processing": True  # Use Celery for large files
            }
        )
    
    data_id = response.json()["id"]
    
    # Poll processing status
    while True:
        status_response = await client.get(
            f"http://localhost:8000/api/v1/data/task/{task_id}/status",
            headers={"Authorization": f"Bearer {token}"}
        )
        status = status_response.json()
        
        if status["state"] in ["SUCCESS", "FAILURE"]:
            break
        
        await asyncio.sleep(2)
    
    # Retrieve processed data
    processed = await client.get(
        f"http://localhost:8000/api/v1/data/{data_id}/processed",
        headers={"Authorization": f"Bearer {token}"}
    )
```

### 3. Run Analysis

```python
# Run MaxQuant analysis
maxquant_response = await client.post(
    "http://localhost:8000/api/v1/proteomics/maxquant",
    headers={"Authorization": f"Bearer {token}"},
    json={
        "sample_id": 1,
        "raw_files": ["/data/mzml/sample_001.mzML"],
        "fasta_file": "/data/uniprot_human.fasta",
        "output_dir": "/data/maxquant_output",
        "lfq": True,
        "match_between_runs": True
    }
)
```

## Format Comparison

### mzML vs mzXML

| Feature | mzML | mzXML |
|---------|------|-------|
| Standard | PSI-MS (current) | ISB (legacy) |
| File Size | Smaller (compression) | Larger |
| Metadata | Rich CV terms | Limited |
| Compatibility | Better | Older tools |
| **Recommendation** | ✅ Use for new projects | Use for legacy compatibility |

### When to Use Each Format

**mzML**:
- Modern workflows
- Rich metadata requirements
- Integration with Proteome Xchange
- Storage efficiency

**mzXML**:
- Legacy pipeline compatibility
- Older search engine requirements
- Some visualization tools

**MGF**:
- Input to search engines (Mascot, X!Tandem)
- Simple MS/MS spectrum representation
- Human-readable format

## Performance Considerations

### File Sizes

| Format | Typical Size (per scan) | Compression |
|--------|------------------------|-------------|
| RAW | ~100-500 KB | Proprietary |
| mzML | ~50-200 KB | gzip available |
| mzXML | ~100-300 KB | gzip available |
| MGF | ~10-50 KB | Text only |

### Processing Time

| Operation | Small File (< 100 MB) | Large File (> 1 GB) |
|-----------|----------------------|---------------------|
| RAW → mzML conversion | 1-2 min | 10-30 min |
| mzML parsing (pyteomics) | 5-10 sec | 1-5 min |
| mzML parsing (XML) | 10-20 sec | 2-10 min |
| Upload + process | 10-30 sec | Async required |

### Best Practices

1. **Use async processing for files > 100 MB**
   ```python
   async_processing=True
   ```

2. **Compress mzML files**
   ```bash
   gzip sample.mzML
   ```

3. **Index large files**
   ```bash
   # Create index for faster access
   msconvert sample.raw --mzML --filter "index"
   ```

4. **Batch conversions**
   ```bash
   # Convert multiple files in parallel
   parallel ThermoRawFileParser -i {} -o output/ -f 2 ::: *.raw
   ```

## Troubleshooting

### ThermoRawFileParser Not Found

```bash
# Install on Linux/Mac
wget https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.2/ThermoRawFileParser1.4.2.zip
unzip ThermoRawFileParser1.4.2.zip
chmod +x ThermoRawFileParser.sh

# Add to PATH or use full path
export PATH=$PATH:/path/to/ThermoRawFileParser
```

### Pyteomics Import Error

```bash
# Install pyteomics
pip install pyteomics>=4.6.0

# Verify installation
python -c "from pyteomics import mzml; print('OK')"
```

### Large File Timeout

```python
# Use async processing
response = await client.post(
    "/api/v1/data/upload",
    files={"file": large_file},
    data={"async_processing": True}  # Returns task_id immediately
)

task_id = response.json()["task_id"]
```

### Memory Issues

```bash
# Process in chunks (for local scripts)
with pymzml.read(file_path, use_index=True) as reader:
    for idx in range(0, reader.get_spectrum_count(), 1000):
        chunk = [reader[i] for i in range(idx, min(idx+1000, reader.get_spectrum_count()))]
        process_chunk(chunk)
```

## Additional Resources

- **ThermoRawFileParser**: https://github.com/compomics/ThermoRawFileParser
- **Pyteomics**: https://pyteomics.readthedocs.io/
- **PSI-MS mzML**: http://www.psidev.info/mzml
- **ProteoWizard**: http://proteowizard.sourceforge.net/
- **Proteome Xchange**: http://www.proteomexchange.org/

## API Reference

See the full API documentation at `/docs` for:
- `/api/v1/proteomics/convert-raw` - RAW to mzML conversion
- `/api/v1/proteomics/maxquant` - MaxQuant analysis
- `/api/v1/proteomics/msfragger` - MSFragger search
- `/api/v1/data/upload` - File upload with processing
- `/api/v1/data/{id}/processed` - Retrieve processed data

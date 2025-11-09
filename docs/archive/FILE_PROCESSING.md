# File Processing Module

## Overview

The file processing module automatically processes uploaded data files and converts them to a unified format. This enables:

- **Consistent Data Format**: All omics data stored in a standardized JSON structure
- **Format Detection**: Automatic detection of file type and omics category
- **Multi-Format Support**: Support for various genomics, transcriptomics, proteomics formats
- **Dual Storage**: Both raw and processed versions are stored
- **Easy Analysis**: Processed data can be directly used in pipelines

---

## Architecture

### Components

1. **FileProcessor** (`backend/app/services/file_processor.py`)

   - Main service for file processing
   - Auto-detects omics type and format
   - Coordinates conversion to unified format

2. **Converters** (`backend/app/converters/`)

   - **base.py**: Base converter class and factory
   - **genomics.py**: VCF, BED, GTF/GFF converters
   - **transcriptomics.py**: CSV, TSV, Excel converters
   - Additional converters can be easily added

3. **Unified Format** (`backend/app/schemas/unified_format.py`)
   - Standard JSON schema for all omics data
   - Includes metadata, records, and statistics

---

## Supported Formats

### Genomics

- **Input**: VCF, VCF.GZ, BED, GTF, GFF, GFF3, FASTA, FASTQ
- **Output**: VCF, BED, PLINK

### Transcriptomics

- **Input**: CSV, TSV, TXT, XLSX, XLS (expression matrices, count tables)
- **Output**: CSV, TSV, DESeq2 format

### Proteomics

- **Input**: mzML, mzXML, MGF (Mascot Generic Format), CSV, TSV
- **Output**: CSV, TSV, MGF
- **Features**:
  - Extracts spectrum information (m/z, intensity, retention time)
  - Supports tabular protein identification results
  - MS/MS spectra parsing

### Metabolomics

- **Input**: mzData, mzML, netCDF/CDF, CSV, TSV
- **Output**: CSV, TSV
- **Features**:
  - Mass spectrometry data extraction
  - Metabolite identification and quantification tables
  - Note: Full netCDF support requires netCDF4 library

---

## API Usage

### Upload with Processing

**Endpoint**: `POST /data/upload`

**Parameters**:

- `project_id` (required): Project ID
- `file` (required): File to upload
- `process_file` (optional, default=true): Enable/disable processing
- `sample_id` (optional): Sample identifier
- `organism` (optional): Organism name (e.g., "Homo sapiens")
- `reference_genome` (optional): Reference genome (e.g., "hg38")

**Example (curl)**:

```bash
curl -X POST "http://localhost:8000/data/upload" \
  -H "Authorization: Bearer $TOKEN" \
  -F "project_id=1" \
  -F "file=@data.vcf" \
  -F "process_file=true" \
  -F "organism=Homo sapiens" \
  -F "reference_genome=hg38"
```

**Response**:

```json
{
  "id": 123,
  "filename": "data.vcf",
  "size": 1234567,
  "project_id": 1,
  "metadata_": {
    "processed_file_id": 124,
    "processing_info": {
      "omics_type": "genomics",
      "source_format": "vcf",
      "sample_id": "sample_1_data",
      "converted": true,
      "record_count": 5000
    },
    "omics_type": "genomics"
  },
  "processing": {
    "processed": true,
    "processed_file_id": 124,
    "processed_filename": "data.vcf.processed.json",
    "omics_type": "genomics",
    "source_format": "vcf",
    "converted": true,
    "record_count": 5000
  }
}
```

---

### Get Processed Data

**Endpoint**: `GET /data/{datafile_id}/processed`

**Example**:

```bash
curl -H "Authorization: Bearer $TOKEN" \
  "http://localhost:8000/data/123/processed"
```

**Response**:

```json
{
  "datafile_id": 123,
  "processed_file_id": 124,
  "unified_data": {
    "metadata": {
      "omics_type": "genomics",
      "source_format": "vcf",
      "sample_id": "sample_1_data",
      "organism": "Homo sapiens",
      "reference_genome": "hg38",
      "processing_steps": [],
      "created_at": "2025-01-09T10:00:00Z"
    },
    "headers": ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
    "records": [
      {
        "id": "chr1:12345:A:G",
        "values": {
          "CHROM": "chr1",
          "POS": 12345,
          "REF": "A",
          "ALT": "G",
          "QUAL": 99.9
        }
      }
    ],
    "statistics": {
      "total_variants": 5000,
      "file_size_bytes": 1234567
    }
  },
  "processing_info": {
    "converted": true,
    "record_count": 5000
  }
}
```

---

## Unified Data Format

### Structure

```json
{
  "metadata": {
    "omics_type": "genomics|transcriptomics|proteomics|metabolomics",
    "source_format": "vcf|bed|csv|tsv|...",
    "sample_id": "sample_identifier",
    "organism": "Homo sapiens",
    "reference_genome": "hg38",
    "processing_steps": [
      {
        "step_name": "quality_control",
        "tool": "fastqc",
        "version": "0.11.9",
        "parameters": {}
      }
    ],
    "created_at": "2025-01-09T10:00:00Z",
    "custom_fields": {}
  },
  "headers": ["column1", "column2", ...],
  "records": [
    {
      "id": "unique_id",
      "values": {
        "column1": "value1",
        "column2": 123
      }
    }
  ],
  "statistics": {
    "total_records": 5000,
    "file_size_bytes": 1234567
  }
}
```

---

## How It Works

### Upload Flow

```
1. User uploads file
   ↓
2. Raw file stored (encrypted in S3)
   ↓
3. FileProcessor detects format and omics type
   ↓
4. Appropriate converter is selected
   ↓
5. File converted to unified format
   ↓
6. Processed version stored as JSON
   ↓
7. Metadata links raw ↔ processed files
   ↓
8. Both files returned to user
```

### Processing Logic

```python
# Detect omics type
omics_type = FileProcessor.detect_omics_type(filename, content)

# Detect format
format = FileProcessor.detect_format(filename, content)

# Get converter
converter = ConverterFactory.create(omics_type)

# Convert to unified
unified_data = await converter.to_unified(
    file_path=tmp_path,
    sample_id=sample_id,
    source_format=format,
    organism=organism,
    reference_genome=reference_genome
)

# Serialize and store
json_bytes = FileProcessor.serialize_unified_data(unified_data)
```

---

## Adding New Converters

### 1. Create Converter Class

```python
# backend/app/converters/my_omics.py

from app.converters.base import BaseConverter, ConverterFactory
from app.schemas.unified_format import OmicsType

@ConverterFactory.register(OmicsType.MY_OMICS)
class MyOmicsConverter(BaseConverter):

    async def to_unified(self, file_path, sample_id, source_format, **kwargs):
        # Read file
        # Parse data
        # Create UnifiedData object
        return unified_data

    async def from_unified(self, unified_data, target_format, output_path, **kwargs):
        # Convert from unified to target format
        return output_path
```

### 2. Register Format Detection

```python
# Update FileProcessor.EXTENSION_TO_OMICS
EXTENSION_TO_OMICS = {
    ".myformat": OmicsType.MY_OMICS,
}
```

### 3. Import in **init**.py

```python
# backend/app/converters/__init__.py
from .my_omics import MyOmicsConverter
```

---

## Error Handling

### If Conversion Fails

- Raw file is still stored ✅
- Processing info includes error message
- File marked as "raw" (not processed)
- Can be reprocessed later

### Example Error Response

```json
{
  "id": 123,
  "filename": "data.unknown",
  "processing": {
    "processed": false,
    "reason": "Unsupported format for conversion"
  }
}
```

---

## Benefits

1. **Data Standardization**: All data in consistent format
2. **Easy Analysis**: Processed data ready for pipelines
3. **Format Flexibility**: Support for multiple input formats
4. **Metadata Tracking**: Comprehensive metadata stored
5. **Dual Access**: Access both raw and processed versions
6. **Error Tolerance**: Upload succeeds even if processing fails
7. **Extensibility**: Easy to add new format converters

---

## Frontend Integration

### Upload Form

```typescript
const formData = new FormData();
formData.append("project_id", projectId);
formData.append("file", file);
formData.append("process_file", "true");
formData.append("organism", "Homo sapiens");
formData.append("reference_genome", "hg38");

const response = await api.post("/data/upload", formData);

if (response.data.processing.processed) {
  console.log("File processed successfully!");
  console.log(
    `Processed file ID: ${response.data.processing.processed_file_id}`
  );
}
```

### View Processed Data

```typescript
const response = await api.get(`/data/${fileId}/processed`);
const unifiedData = response.data.unified_data;

// Access records
unifiedData.records.forEach((record) => {
  console.log(record.id, record.values);
});

// View statistics
console.log(unifiedData.statistics);
```

---

## Testing

### Test Script

```bash
# Upload VCF file
curl -X POST "http://localhost:8000/data/upload" \
  -H "Authorization: Bearer $TOKEN" \
  -F "project_id=1" \
  -F "file=@test.vcf" \
  -F "organism=Homo sapiens" \
  -F "reference_genome=hg38"

# Get processed data
curl -H "Authorization: Bearer $TOKEN" \
  "http://localhost:8000/data/{file_id}/processed"
```

---

## Configuration

### Enable/Disable Processing

```python
# Per-upload
process_file=False  # Skip processing, store raw only

# Global (future enhancement)
# app/settings.py
ENABLE_FILE_PROCESSING = True
```

---

## Performance Considerations

- **Large Files**: Processing happens asynchronously (future enhancement)
- **Memory**: Files loaded in chunks where possible
- **Storage**: Processed files typically smaller (JSON compression)
- **Caching**: Processed data can be cached for faster access

---

## Future Enhancements

1. **Async Processing**: Background job queue for large files
2. **Progress Tracking**: Real-time processing progress
3. **Format Validation**: Pre-upload format validation
4. **Batch Processing**: Process multiple files at once
5. **Format Conversion**: Convert between formats via API
6. **Quality Control**: Automated QC checks during processing
7. **Visualization**: Preview processed data in UI

---

## Summary

The file processing module provides:

- ✅ Automatic format detection
- ✅ Multi-format support
- ✅ Unified data structure
- ✅ Dual storage (raw + processed)
- ✅ Easy extensibility
- ✅ Error tolerance
- ✅ Rich metadata

All uploaded files are automatically processed to enable seamless analysis workflows!

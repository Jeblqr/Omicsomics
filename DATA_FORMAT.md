# Unified Data Format & Converters

## Overview

Omicsomics uses a **unified JSON-based format** for all omics data types. This enables:

- **Consistent data handling** across different omics
- **Easy format conversion** between tools
- **Reproducible data processing** with complete metadata
- **Cross-platform compatibility**

## Unified Format Structure

All omics data follows this standard structure:

```json
{
  "metadata": {
    "omics_type": "genomics|transcriptomics|proteomics|metabolomics|epigenomics",
    "format_version": "1.0.0",
    "source_format": "vcf|fastq|csv|...",
    "sample_id": "sample_123",
    "organism": "Homo sapiens",
    "reference_genome": "GRCh38",
    "created_at": "2025-01-09T14:00:00Z",
    "updated_at": "2025-01-09T14:00:00Z",
    "processing_steps": [
      {
        "step_name": "quality_control",
        "tool": "FastQC",
        "version": "0.11.9",
        "parameters": {"threads": 4},
        "timestamp": "2025-01-09T14:05:00Z"
      }
    ],
    "custom_fields": {}
  },
  "headers": ["CHROM", "POS", "REF", "ALT", "QUAL"],
  "records": [
    {
      "id": "chr1:1000:A:G",
      "values": {
        "CHROM": "chr1",
        "POS": 1000,
        "REF": "A",
        "ALT": "G",
        "QUAL": 30.5
      }
    }
  ],
  "statistics": {
    "total_records": 1,
    "file_size_bytes": 12345
  }
}
```

## Supported Formats

### Genomics (WGS/WES)

**Input Formats:**
- VCF/VCF.GZ - Variant Call Format
- BAM/SAM - Binary/Sequence Alignment Map
- FASTA/FASTQ - Sequence files
- BED - Browser Extensible Data
- GFF/GTF/GFF3 - Gene annotation

**Export Formats:**
- VCF - For GATK, bcftools
- BED - For bedtools, UCSC
- PLINK binary (.bed/.bim/.fam) - For PLINK analysis

### Transcriptomics (RNA-seq)

**Input Formats:**
- CSV/TSV/TXT - Count matrices, expression tables
- FASTQ - Raw sequencing reads
- BAM - Aligned reads
- H5/H5AD - HDF5 formats
- Matrix Market (.mtx) - Sparse matrices

**Export Formats:**
- CSV/TSV - For DESeq2, edgeR
- H5AD - For Scanpy
- RDS - For Seurat

### Proteomics (LC-MS/MS)

**Input Formats:**
- mzML/mzXML - Mass spectrometry data
- MGF - Mascot Generic Format
- CSV/TSV - Protein/peptide tables
- MaxQuant outputs (peptides.txt, proteinGroups.txt)

**Export Formats:**
- mzML - For OpenMS, MSFragger
- CSV - For custom analysis
- MaxQuant format - For reanalysis

### Metabolomics (LC-MS/GC-MS)

**Input Formats:**
- mzML/mzXML - Mass spectrometry data
- NetCDF/CDF - Chromatography data
- CSV/TSV - Feature tables
- mzData - Mass spec data

**Export Formats:**
- mzML - For XCMS, MZmine
- CSV - For MetaboAnalyst
- NetCDF - For AMDIS

### Epigenomics (ChIP-seq, ATAC-seq)

**Input Formats:**
- BED/narrowPeak/broadPeak - Peak files
- BigWig/bedGraph - Signal tracks
- BAM - Aligned reads
- CSV/TSV - Methylation tables

**Export Formats:**
- BED - For MACS2, Homer
- BigWig - For IGV, UCSC
- BedGraph - For visualization

### Single-Cell (scRNA-seq)

**Input Formats:**
- H5/H5AD - AnnData format
- Matrix Market (.mtx) - 10x format
- CSV/TSV - Count matrices
- Loom - Single-cell format
- RDS - Seurat objects

**Export Formats:**
- H5AD - For Scanpy
- RDS - For Seurat
- CSV - For custom analysis

## Using the Converters

### API Endpoints

#### Convert to Unified Format

```bash
# Upload and convert a VCF file
curl -X POST "http://localhost:8001/api/v1/converters/convert/to-unified" \
  -H "Authorization: Bearer $TOKEN" \
  -F "file=@variants.vcf" \
  -F "omics_type=genomics" \
  -F "sample_id=sample_001" \
  -F "source_format=vcf" \
  -F "organism=Homo sapiens" \
  -F "reference_genome=GRCh38"
```

Response:
```json
{
  "metadata": {...},
  "headers": [...],
  "records": [...],
  "statistics": {...}
}
```

#### Convert from Unified Format

```bash
# Convert unified format to BED
curl -X POST "http://localhost:8001/api/v1/converters/convert/from-unified" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "unified_data": {...},
    "target_format": "bed"
  }'
```

#### Get Supported Formats

```bash
# Get supported formats for genomics
curl "http://localhost:8001/api/v1/converters/supported-formats/genomics" \
  -H "Authorization: Bearer $TOKEN"
```

Response:
```json
{
  "omics_type": "genomics",
  "supported_formats": ["vcf", "vcf.gz", "bam", "sam", "fasta", "fastq", "bed", "gff", "gtf"]
}
```

#### Get Export Formats for Tools

```bash
# Get export formats for PLINK
curl "http://localhost:8001/api/v1/converters/export-formats/plink" \
  -H "Authorization: Bearer $TOKEN"
```

Response:
```json
{
  "tool": "plink",
  "export_formats": ["bed", "bim", "fam"]
}
```

## Python API

### Converting Data

```python
from app.converters import ConverterFactory
from app.schemas.unified_format import OmicsType
from pathlib import Path

# Create converter
converter = ConverterFactory.create(OmicsType.GENOMICS)

# Convert to unified format
unified_data = await converter.to_unified(
    file_path=Path("data.vcf"),
    sample_id="sample_001",
    source_format="vcf",
    organism="Homo sapiens",
    reference_genome="GRCh38"
)

# Save unified format
await converter.save_unified(
    unified_data,
    Path("output/unified_data.json")
)

# Convert from unified to target format
output_path = await converter.from_unified(
    unified_data=unified_data,
    target_format="bed",
    output_path=Path("output/data.bed")
)
```

### Registering Custom Converters

```python
from app.converters.base import BaseConverter, ConverterFactory
from app.schemas.unified_format import OmicsType, UnifiedData

@ConverterFactory.register(OmicsType.METABOLOMICS)
class MetabolomicsConverter(BaseConverter):
    
    async def to_unified(self, file_path, sample_id, source_format, **kwargs):
        # Your conversion logic here
        metadata = self.create_metadata(
            sample_id=sample_id,
            source_format=source_format,
            **kwargs
        )
        
        # Parse file and create records
        records = []
        # ... parsing logic ...
        
        return UnifiedData(
            metadata=metadata,
            headers=headers,
            records=records,
            statistics=statistics
        )
    
    async def from_unified(self, unified_data, target_format, output_path, **kwargs):
        # Your export logic here
        # ... write to output_path ...
        return output_path
```

## Processing Steps Tracking

Every conversion and processing step is tracked in metadata:

```python
# Add processing step
unified_data = converter.add_processing_step(
    unified_data,
    step_name="variant_calling",
    tool="GATK",
    version="4.2.0.0",
    min_quality=30,
    ploidy=2
)
```

This creates:
```json
{
  "processing_steps": [
    {
      "step_name": "variant_calling",
      "tool": "GATK",
      "version": "4.2.0.0",
      "parameters": {
        "min_quality": 30,
        "ploidy": 2
      },
      "timestamp": "2025-01-09T14:10:00Z"
    }
  ]
}
```

## Benefits

### 1. Reproducibility
- Complete processing history
- Tool versions tracked
- Parameters preserved

### 2. Interoperability
- Easy data exchange between tools
- Consistent format across platforms
- Standard API for all omics types

### 3. Data Integration
- Unified format enables multi-omics analysis
- Common structure for all data types
- Standardized field names

### 4. Quality Control
- Embedded statistics
- Format validation
- Metadata completeness checks

## Extension Points

### Adding New Omics Types

1. Define in `OmicsType` enum
2. Add to `SUPPORTED_FORMATS` mapping
3. Create converter class
4. Register with `@ConverterFactory.register()`

### Adding New File Formats

1. Update `SUPPORTED_FORMATS` for omics type
2. Implement parser in converter class
3. Add export method if needed
4. Update documentation

### Custom Metadata Fields

Use `custom_fields` dict in metadata:

```python
metadata.custom_fields['sequencing_platform'] = 'Illumina NovaSeq'
metadata.custom_fields['library_prep'] = 'TruSeq'
metadata.custom_fields['read_length'] = 150
```

## Future Enhancements

- [ ] Streaming conversion for large files
- [ ] Parallel processing for multiple files
- [ ] Cloud storage integration (S3, GCS, Azure)
- [ ] Format auto-detection
- [ ] Compression support
- [ ] Checksum validation
- [ ] Schema versioning and migration
- [ ] REST API for bulk conversion
- [ ] CLI tool for local conversion
- [ ] Web UI for format conversion

## See Also

- [Pipeline Format](PIPELINE_FORMAT.md) - Unified pipeline definition
- [API Documentation](http://localhost:8001/docs) - Full API reference
- [User Guide](docs/PIPELINE_BUILDER.md) - Using the platform

"""Unified data format schemas for different omics types."""

from enum import Enum
from typing import Any, Dict, List, Optional
from datetime import datetime
from pydantic import BaseModel, Field


class OmicsType(str, Enum):
    """Supported omics types."""
    GENOMICS = "genomics"
    TRANSCRIPTOMICS = "transcriptomics"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    EPIGENOMICS = "epigenomics"
    SINGLE_CELL = "single_cell"
    MULTI_OMICS = "multi_omics"


class ProcessingStep(BaseModel):
    """Record of a processing step applied to the data."""
    step_name: str
    tool: str
    version: Optional[str] = None
    parameters: Dict[str, Any] = Field(default_factory=dict)
    timestamp: datetime = Field(default_factory=datetime.utcnow)


class UnifiedMetadata(BaseModel):
    """Common metadata for all omics data."""
    omics_type: OmicsType
    format_version: str = "1.0.0"
    source_format: str
    sample_id: str
    organism: Optional[str] = None
    reference_genome: Optional[str] = None
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    processing_steps: List[ProcessingStep] = Field(default_factory=list)
    custom_fields: Dict[str, Any] = Field(default_factory=dict)


class UnifiedDataRecord(BaseModel):
    """A single data record in the unified format."""
    id: str
    values: Dict[str, Any]
    

class UnifiedData(BaseModel):
    """Unified data format for all omics types."""
    metadata: UnifiedMetadata
    headers: List[str]
    records: List[UnifiedDataRecord]
    statistics: Dict[str, Any] = Field(default_factory=dict)
    

# Format-specific extensions

class GenomicsData(UnifiedData):
    """Unified format for genomics data (VCF, BAM, etc.)."""
    pass


class TranscriptomicsData(UnifiedData):
    """Unified format for transcriptomics data (count matrix, expression)."""
    pass


class ProteomicsData(UnifiedData):
    """Unified format for proteomics data (peptides, proteins)."""
    pass


class MetabolomicsData(UnifiedData):
    """Unified format for metabolomics data (features, compounds)."""
    pass


class EpigenomicsData(UnifiedData):
    """Unified format for epigenomics data (peaks, methylation)."""
    pass


class SingleCellData(UnifiedData):
    """Unified format for single-cell data (cell x gene matrix)."""
    pass


# Supported input formats mapping
SUPPORTED_FORMATS = {
    OmicsType.GENOMICS: [
        "vcf", "vcf.gz", "bam", "sam", "fasta", "fastq", 
        "fastq.gz", "bed", "gff", "gtf", "gff3"
    ],
    OmicsType.TRANSCRIPTOMICS: [
        "csv", "tsv", "txt", "count", "matrix", "h5", 
        "h5ad", "fastq", "fastq.gz", "bam"
    ],
    OmicsType.PROTEOMICS: [
        "mzml", "mzxml", "mgf", "csv", "tsv", "txt",
        "maxquant", "peptides.txt", "proteinGroups.txt"
    ],
    OmicsType.METABOLOMICS: [
        "mzml", "mzxml", "csv", "tsv", "txt", "mzdata",
        "netcdf", "cdf"
    ],
    OmicsType.EPIGENOMICS: [
        "bed", "narrowPeak", "broadPeak", "bigWig", "bw",
        "wig", "bedGraph", "bam", "csv", "tsv"
    ],
    OmicsType.SINGLE_CELL: [
        "h5", "h5ad", "mtx", "matrix.mtx", "csv", "tsv",
        "loom", "rds"
    ]
}


# Target formats for export (tool-specific)
EXPORT_FORMATS = {
    "plink": ["bed", "bim", "fam"],
    "gatk": ["vcf", "vcf.gz"],
    "star": ["fastq", "fastq.gz"],
    "salmon": ["fastq", "fastq.gz"],
    "deseq2": ["csv", "txt"],
    "seurat": ["csv", "h5"],
    "scanpy": ["h5ad"],
    "maxquant": ["mzml", "raw"],
    "xcms": ["mzml", "mzxml", "netcdf"],
    "macs2": ["bed", "bam"]
}

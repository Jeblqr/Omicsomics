"""
File Processing Service
Processes uploaded raw files and converts them to unified format
"""

import io
import tempfile
from pathlib import Path
from typing import BinaryIO, Dict, Optional, Tuple, Any
import json
import logging

from app.converters.base import ConverterFactory
from app.schemas.unified_format import OmicsType, UnifiedData

logger = logging.getLogger(__name__)


class FileProcessor:
    """
    Handles processing of uploaded files:
    1. Detect file type and omics type
    2. Convert to unified format
    3. Store both raw and processed versions
    """

    # File extension to omics type mapping
    EXTENSION_TO_OMICS = {
        # Genomics
        ".vcf": OmicsType.GENOMICS,
        ".vcf.gz": OmicsType.GENOMICS,
        ".bed": OmicsType.GENOMICS,
        ".gtf": OmicsType.GENOMICS,
        ".gff": OmicsType.GENOMICS,
        ".gff3": OmicsType.GENOMICS,
        ".fasta": OmicsType.GENOMICS,
        ".fa": OmicsType.GENOMICS,
        ".fastq": OmicsType.GENOMICS,
        ".fq": OmicsType.GENOMICS,
        # Transcriptomics
        ".counts": OmicsType.TRANSCRIPTOMICS,
        ".tsv": OmicsType.TRANSCRIPTOMICS,  # Could be multiple types
        ".txt": OmicsType.TRANSCRIPTOMICS,  # Could be multiple types
        # Proteomics
        ".mzml": OmicsType.PROTEOMICS,
        ".mzxml": OmicsType.PROTEOMICS,
        ".mgf": OmicsType.PROTEOMICS,
        ".raw": OmicsType.PROTEOMICS,  # Thermo RAW
        # Metabolomics
        ".mzdata": OmicsType.METABOLOMICS,
        ".cdf": OmicsType.METABOLOMICS,
        ".netcdf": OmicsType.METABOLOMICS,
        ".cdf": OmicsType.METABOLOMICS,
        # Generic data
        ".csv": OmicsType.TRANSCRIPTOMICS,  # Default to transcriptomics for CSV
        ".xlsx": OmicsType.TRANSCRIPTOMICS,
        ".xls": OmicsType.TRANSCRIPTOMICS,
    }

    # Format detection patterns
    FORMAT_PATTERNS = {
        "vcf": [b"##fileformat=VCF", b"#CHROM"],
        "bed": [b"chr", b"\t"],  # Basic BED format detection
        "gtf": [b"\t", b"gene_id"],
        "counts": [b"gene", b"sample", b"\t"],
        "fasta": [b">"],
        "fastq": [b"@"],
    }

    @classmethod
    def detect_omics_type(
        cls, filename: str, file_content: Optional[bytes] = None
    ) -> OmicsType:
        """
        Detect omics type from filename and optionally file content.

        Args:
            filename: Original filename
            file_content: First few KB of file content for content-based detection

        Returns:
            OmicsType enum value
        """
        filename_lower = filename.lower()

        # For CSV files, use content-based detection first
        if filename_lower.endswith(".csv") and file_content:
            content_str = file_content[:2000].decode("utf-8", errors="ignore").lower()

            # Check for proteomics indicators
            if any(
                keyword in content_str
                for keyword in ["mz", "protein_id", "peptide", "charge", "intensity"]
            ):
                return OmicsType.PROTEOMICS

            # Check for metabolomics indicators
            elif any(
                keyword in content_str
                for keyword in ["metabolite", "compound", "hmdb", "kegg"]
            ):
                return OmicsType.METABOLOMICS

            # Check for genomics indicators
            elif any(
                keyword in content_str
                for keyword in ["chr", "chrom", "position", "variant", "ref", "alt"]
            ):
                return OmicsType.GENOMICS

            # Default to transcriptomics for CSV with gene indicators or generic sample data
            # (gene_id, gene_name, or sample columns)
            else:
                return OmicsType.TRANSCRIPTOMICS

        # Try extension-based detection for non-CSV files
        for ext, omics_type in cls.EXTENSION_TO_OMICS.items():
            if filename_lower.endswith(ext):
                return omics_type

        # If content provided, try content-based detection
        if file_content:
            if b"##fileformat=VCF" in file_content[:1000]:
                return OmicsType.GENOMICS
            elif file_content.startswith(b">"):
                return OmicsType.GENOMICS  # FASTA
            elif file_content.startswith(b"@"):
                return OmicsType.GENOMICS  # FASTQ
            elif b"mzML" in file_content[:1000] or b"mzXML" in file_content[:1000]:
                return OmicsType.PROTEOMICS

        # Default fallback
        logger.warning(
            f"Could not detect omics type for {filename}, defaulting to TRANSCRIPTOMICS"
        )
        return OmicsType.TRANSCRIPTOMICS

    @classmethod
    def detect_format(cls, filename: str, file_content: Optional[bytes] = None) -> str:
        """
        Detect file format from filename and content.

        Args:
            filename: Original filename
            file_content: First few KB of file content

        Returns:
            Format string (e.g., "vcf", "bed", "csv")
        """
        filename_lower = filename.lower()

        # Extension-based detection
        if filename_lower.endswith(".vcf.gz") or filename_lower.endswith(".vcf"):
            return "vcf"
        elif filename_lower.endswith(".bed"):
            return "bed"
        elif filename_lower.endswith(".gtf"):
            return "gtf"
        elif filename_lower.endswith(".gff") or filename_lower.endswith(".gff3"):
            return "gff"
        elif filename_lower.endswith(".csv"):
            return "csv"
        elif filename_lower.endswith(".tsv") or filename_lower.endswith(".txt"):
            return "tsv"
        elif filename_lower.endswith(".fasta") or filename_lower.endswith(".fa"):
            return "fasta"
        elif filename_lower.endswith(".fastq") or filename_lower.endswith(".fq"):
            return "fastq"
        elif filename_lower.endswith(".mzml"):
            return "mzml"
        elif filename_lower.endswith(".mzxml"):
            return "mzxml"

        # Content-based detection
        if file_content:
            if b"##fileformat=VCF" in file_content[:1000]:
                return "vcf"
            elif file_content.startswith(b">"):
                return "fasta"
            elif file_content.startswith(b"@"):
                return "fastq"

        # Default
        return "unknown"

    @classmethod
    async def process_file(
        cls,
        file_data: BinaryIO,
        filename: str,
        project_id: int,
        sample_id: Optional[str] = None,
        omics_type: Optional[OmicsType] = None,
        metadata: Optional[Dict] = None,
    ) -> Tuple[UnifiedData, Dict[str, Any]]:
        """
        Process uploaded file and convert to unified format.

        Args:
            file_data: File binary data
            filename: Original filename
            project_id: Project ID
            sample_id: Sample identifier (optional, will generate if not provided)
            omics_type: Force specific omics type (optional, will auto-detect)
            metadata: Additional metadata (organism, reference_genome, etc.)

        Returns:
            Tuple of (UnifiedData object, processing_info dict)
        """
        # Read content for detection
        file_data.seek(0)
        content_sample = file_data.read(10000)  # Read first 10KB for detection
        file_data.seek(0)

        # Detect omics type if not provided
        if omics_type is None:
            omics_type = cls.detect_omics_type(filename, content_sample)

        # Detect format
        source_format = cls.detect_format(filename, content_sample)

        # Generate sample_id if not provided
        if sample_id is None:
            sample_id = f"sample_{project_id}_{Path(filename).stem}"

        # Prepare metadata
        conversion_metadata = metadata or {}

        processing_info: Dict[str, Any] = {
            "omics_type": omics_type.value,
            "source_format": source_format,
            "sample_id": sample_id,
            "original_filename": filename,
        }

        try:
            # Create temporary file for conversion
            with tempfile.NamedTemporaryFile(
                mode="wb", delete=False, suffix=Path(filename).suffix
            ) as tmp_file:
                file_data.seek(0)
                tmp_file.write(file_data.read())
                tmp_path = Path(tmp_file.name)

            try:
                # Get appropriate converter
                converter = ConverterFactory.create(omics_type)

                # Check if format is supported
                if not converter.is_supported_format(source_format):
                    logger.warning(
                        f"Format {source_format} not supported for {omics_type}, "
                        f"storing as raw without conversion"
                    )
                    processing_info["converted"] = False
                    processing_info["reason"] = "Unsupported format for conversion"

                    # Return minimal unified data structure
                    unified_data = cls._create_raw_data_structure(
                        sample_id=sample_id,
                        omics_type=omics_type,
                        source_format=source_format,
                        filename=filename,
                        metadata=conversion_metadata,
                    )
                    return unified_data, processing_info

                # Convert to unified format
                unified_data = await converter.to_unified(
                    file_path=tmp_path,
                    sample_id=sample_id,
                    source_format=source_format,
                    **conversion_metadata,
                )

                processing_info["converted"] = True
                processing_info["record_count"] = len(unified_data.records)
                processing_info["statistics"] = unified_data.statistics

                logger.info(
                    f"Successfully converted {filename} to unified format: "
                    f"{len(unified_data.records)} records"
                )

                return unified_data, processing_info

            finally:
                # Clean up temporary file
                tmp_path.unlink(missing_ok=True)

        except Exception as e:
            logger.error(f"Error processing file {filename}: {str(e)}")
            processing_info["converted"] = False
            processing_info["error"] = str(e)

            # Return raw data structure on error
            unified_data = cls._create_raw_data_structure(
                sample_id=sample_id,
                omics_type=omics_type,
                source_format=source_format,
                filename=filename,
                metadata=conversion_metadata,
            )
            return unified_data, processing_info

    @classmethod
    def _create_raw_data_structure(
        cls,
        sample_id: str,
        omics_type: OmicsType,
        source_format: str,
        filename: str,
        metadata: Dict,
    ) -> UnifiedData:
        """
        Create a minimal unified data structure for raw (unconverted) files.
        """
        from app.schemas.unified_format import UnifiedMetadata

        unified_metadata = UnifiedMetadata(
            omics_type=omics_type,
            source_format=source_format,
            sample_id=sample_id,
            organism=metadata.get("organism"),
            reference_genome=metadata.get("reference_genome"),
            custom_fields={
                "original_filename": filename,
                "processing_status": "raw",
                **metadata,
            },
        )

        return UnifiedData(
            metadata=unified_metadata,
            headers=[],
            records=[],
            statistics={"processing_status": "raw", "original_filename": filename},
        )

    @classmethod
    def serialize_unified_data(cls, unified_data: UnifiedData) -> bytes:
        """
        Serialize UnifiedData object to JSON bytes for storage.

        Args:
            unified_data: UnifiedData object

        Returns:
            JSON bytes
        """
        json_str = json.dumps(
            unified_data.model_dump(mode="json"),
            indent=2,
            default=str,
        )
        return json_str.encode("utf-8")

    @classmethod
    def deserialize_unified_data(cls, data_bytes: bytes) -> UnifiedData:
        """
        Deserialize JSON bytes to UnifiedData object.

        Args:
            data_bytes: JSON bytes

        Returns:
            UnifiedData object
        """
        json_str = data_bytes.decode("utf-8")
        data_dict = json.loads(json_str)
        return UnifiedData(**data_dict)

    @classmethod
    async def process_file_async(
        cls,
        file_content: bytes,
        filename: str,
        sample_id: Optional[str] = None,
        omics_type: Optional[OmicsType] = None,
        metadata: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        """
        Process file from byte content (for Celery tasks).

        Args:
            file_content: File bytes
            filename: Original filename
            sample_id: Optional sample ID
            omics_type: Optional forced omics type
            metadata: Additional metadata

        Returns:
            Dict with success, unified_data, omics_type, and processing_info
        """
        try:
            # Create BytesIO object from content
            file_io = io.BytesIO(file_content)

            # Use existing process_file method
            unified_data, processing_info = await cls.process_file(
                file_data=file_io,
                filename=filename,
                project_id=0,  # Not needed for async processing
                sample_id=sample_id,
                omics_type=omics_type,
                metadata=metadata,
            )

            return {
                "success": True,
                "unified_data": unified_data,
                "omics_type": processing_info.get("omics_type", "unknown"),
                "processing_info": processing_info,
            }

        except Exception as e:
            logger.error(f"Error in process_file_async: {e}", exc_info=True)
            return {"success": False, "error": str(e), "omics_type": "unknown"}

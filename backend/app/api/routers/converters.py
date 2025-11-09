"""API endpoints for data format conversion."""

from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, Form
from sqlalchemy.ext.asyncio import AsyncSession
from pathlib import Path
import tempfile
import shutil
from typing import Optional

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.converters import ConverterFactory
from app.schemas.unified_format import OmicsType, UnifiedData

router = APIRouter()


@router.post("/convert/to-unified")
async def convert_to_unified(
    file: UploadFile = File(...),
    omics_type: str = Form(...),
    sample_id: str = Form(...),
    source_format: str = Form(...),
    organism: Optional[str] = Form(None),
    reference_genome: Optional[str] = Form(None),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Convert uploaded file to unified format.

    Args:
        file: Uploaded data file
        omics_type: Type of omics data (genomics, transcriptomics, etc.)
        sample_id: Sample identifier
        source_format: Original file format (vcf, bed, etc.)
        organism: Optional organism name
        reference_genome: Optional reference genome

    Returns:
        UnifiedData object in JSON format
    """

    try:
        # Validate omics type
        try:
            omics_enum = OmicsType(omics_type.lower())
        except ValueError:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid omics type: {omics_type}. Must be one of: {[e.value for e in OmicsType]}",
            )

        # Create converter
        converter = ConverterFactory.create(omics_enum)

        # Check if format is supported
        if not converter.is_supported_format(source_format):
            raise HTTPException(
                status_code=400,
                detail=f"Unsupported format {source_format} for {omics_type}. Supported: {converter.supported_inputs}",
            )

        # Save uploaded file temporarily
        with tempfile.NamedTemporaryFile(
            delete=False, suffix=f".{source_format}"
        ) as tmp_file:
            shutil.copyfileobj(file.file, tmp_file)
            tmp_path = Path(tmp_file.name)

        try:
            # Convert to unified format
            unified_data = await converter.to_unified(
                file_path=tmp_path,
                sample_id=sample_id,
                source_format=source_format,
                organism=organism,
                reference_genome=reference_genome,
            )

            return unified_data

        finally:
            # Clean up temp file
            tmp_path.unlink(missing_ok=True)

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Conversion failed: {str(e)}")


@router.post("/convert/from-unified")
async def convert_from_unified(
    unified_data: UnifiedData,
    target_format: str,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Convert unified format data to target format.

    Args:
        unified_data: Data in unified format
        target_format: Target format (vcf, bed, plink-bed, etc.)

    Returns:
        Download URL or file path for converted data
    """

    try:
        # Create converter for the omics type
        converter = ConverterFactory.create(unified_data.metadata.omics_type)

        # Create temporary output path
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / f"converted.{target_format}"

            # Convert from unified format
            result_path = await converter.from_unified(
                unified_data=unified_data,
                target_format=target_format,
                output_path=output_path,
            )

            # Read converted file
            with open(result_path, "rb") as f:
                content = f.read()

            return {
                "format": target_format,
                "size_bytes": len(content),
                "sample_id": unified_data.metadata.sample_id,
                "message": "Conversion successful. Use download endpoint to retrieve file.",
            }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Conversion failed: {str(e)}")


@router.get("/supported-formats/{omics_type}")
async def get_supported_formats(
    omics_type: str, current_user: User = Depends(get_current_user)
):
    """
    Get supported input formats for an omics type.

    Args:
        omics_type: Type of omics data

    Returns:
        List of supported formats
    """

    try:
        omics_enum = OmicsType(omics_type.lower())
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid omics type: {omics_type}")

    formats = ConverterFactory.get_supported_formats(omics_enum)
    return {"omics_type": omics_type, "supported_formats": formats}


@router.get("/export-formats/{tool}")
async def get_export_formats(tool: str, current_user: User = Depends(get_current_user)):
    """
    Get supported export formats for a specific tool.

    Args:
        tool: Tool name (plink, gatk, star, etc.)

    Returns:
        List of supported export formats
    """

    formats = ConverterFactory.get_export_formats(tool.lower())
    if not formats:
        raise HTTPException(
            status_code=404, detail=f"No export formats registered for tool: {tool}"
        )

    return {"tool": tool, "export_formats": formats}

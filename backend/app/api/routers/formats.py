"""API routes for data format conversion."""

from fastapi import APIRouter, Depends, HTTPException, BackgroundTasks
from sqlalchemy.orm import Session
from typing import List, Optional
from pydantic import BaseModel
from datetime import datetime
import os

from app.database import get_db
from app.models.format_conversion import FormatConversion, ConversionRule
from app.converters.format_converter import get_format_converter
from app.core.auth import get_current_user
from app.models.user import User

router = APIRouter(prefix="/api/formats", tags=["formats"])


# Schemas
class ConversionRequest(BaseModel):
    """Request for manual format conversion."""
    source_path: str
    target_path: str
    from_format: Optional[str] = None
    to_format: Optional[str] = None
    parameters: Optional[dict] = None


class ConversionPathRequest(BaseModel):
    """Request to get conversion path."""
    from_format: str
    to_format: str


class FormatDetectionRequest(BaseModel):
    """Request to detect file format."""
    file_path: str


class ConversionEstimateRequest(BaseModel):
    """Request for conversion time estimate."""
    file_path: str
    from_format: str
    to_format: str


class ConversionResponse(BaseModel):
    """Response for conversion operations."""
    id: int
    source_path: str
    target_path: str
    from_format: str
    to_format: str
    conversion_path: List[str]
    status: str
    duration_seconds: Optional[float] = None
    conversion_mode: str
    created_at: datetime
    
    class Config:
        from_attributes = True


class FormatInfo(BaseModel):
    """Information about a supported format."""
    format_id: str
    extensions: List[str]
    mime_type: str
    can_convert_to: List[str]
    can_convert_from: List[str]


@router.get("/supported", response_model=List[FormatInfo])
async def get_supported_formats():
    """
    Get list of all supported formats and their conversion capabilities.
    """
    converter = get_format_converter()
    formats = []
    
    for format_id, info in converter.FORMATS.items():
        # Determine conversion capabilities
        can_convert_to = []
        can_convert_from = []
        
        for other_format in converter.FORMATS.keys():
            if format_id == other_format:
                continue
            
            try:
                # Check if conversion is possible
                path = converter.get_conversion_path(format_id, other_format)
                if len(path) > 0:
                    can_convert_to.append(other_format)
                
                path = converter.get_conversion_path(other_format, format_id)
                if len(path) > 0:
                    can_convert_from.append(other_format)
            except ValueError:
                pass
        
        formats.append(FormatInfo(
            format_id=format_id,
            extensions=info['extensions'],
            mime_type=info['mime_type'],
            can_convert_to=can_convert_to,
            can_convert_from=can_convert_from
        ))
    
    return formats


@router.post("/detect")
async def detect_format(request: FormatDetectionRequest):
    """
    Detect format of a file.
    """
    converter = get_format_converter()
    
    if not os.path.exists(request.file_path):
        raise HTTPException(status_code=404, detail="File not found")
    
    format_id = converter.detect_format(request.file_path)
    
    if format_id is None:
        raise HTTPException(status_code=400, detail="Cannot detect format")
    
    return {
        "file_path": request.file_path,
        "format": format_id,
        "extensions": converter.FORMATS[format_id]['extensions'],
        "mime_type": converter.FORMATS[format_id]['mime_type']
    }


@router.post("/conversion-path")
async def get_conversion_path(request: ConversionPathRequest):
    """
    Get conversion path between two formats.
    """
    converter = get_format_converter()
    
    try:
        path = converter.get_conversion_path(request.from_format, request.to_format)
        
        return {
            "from_format": request.from_format,
            "to_format": request.to_format,
            "conversion_path": path,
            "num_steps": len(path) - 1
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/estimate")
async def estimate_conversion_time(request: ConversionEstimateRequest):
    """
    Estimate conversion time for a file.
    """
    converter = get_format_converter()
    
    if not os.path.exists(request.file_path):
        raise HTTPException(status_code=404, detail="File not found")
    
    try:
        estimated_time = converter.estimate_conversion_time(
            request.file_path,
            request.from_format,
            request.to_format
        )
        
        file_size = os.path.getsize(request.file_path)
        
        return {
            "file_path": request.file_path,
            "file_size_bytes": file_size,
            "file_size_mb": file_size / (1024 * 1024),
            "from_format": request.from_format,
            "to_format": request.to_format,
            "estimated_seconds": estimated_time,
            "estimated_minutes": estimated_time / 60
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/convert", response_model=ConversionResponse)
async def convert_format(
    request: ConversionRequest,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Convert a file from one format to another (manual conversion).
    """
    converter = get_format_converter()
    
    # Validate source file exists
    if not os.path.exists(request.source_path):
        raise HTTPException(status_code=404, detail="Source file not found")
    
    # Create database record
    conversion = FormatConversion(
        source_path=request.source_path,
        target_path=request.target_path,
        source_format=request.from_format or "",
        target_format=request.to_format or "",
        conversion_mode="manual",
        status="pending",
        created_by=current_user.id,
        parameters=request.parameters
    )
    
    db.add(conversion)
    db.commit()
    db.refresh(conversion)
    
    # Start conversion in background
    background_tasks.add_task(
        perform_conversion,
        conversion.id,
        request.source_path,
        request.target_path,
        request.from_format,
        request.to_format,
        request.parameters or {}
    )
    
    return conversion


@router.get("/conversions", response_model=List[ConversionResponse])
async def list_conversions(
    skip: int = 0,
    limit: int = 100,
    status: Optional[str] = None,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    List format conversions for current user.
    """
    query = db.query(FormatConversion).filter(
        FormatConversion.created_by == current_user.id
    )
    
    if status:
        query = query.filter(FormatConversion.status == status)
    
    conversions = query.order_by(
        FormatConversion.created_at.desc()
    ).offset(skip).limit(limit).all()
    
    return conversions


@router.get("/conversions/{conversion_id}", response_model=ConversionResponse)
async def get_conversion(
    conversion_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get details of a specific conversion.
    """
    conversion = db.query(FormatConversion).filter(
        FormatConversion.id == conversion_id,
        FormatConversion.created_by == current_user.id
    ).first()
    
    if not conversion:
        raise HTTPException(status_code=404, detail="Conversion not found")
    
    return conversion


@router.delete("/conversions/{conversion_id}")
async def delete_conversion(
    conversion_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Delete a conversion record.
    """
    conversion = db.query(FormatConversion).filter(
        FormatConversion.id == conversion_id,
        FormatConversion.created_by == current_user.id
    ).first()
    
    if not conversion:
        raise HTTPException(status_code=404, detail="Conversion not found")
    
    db.delete(conversion)
    db.commit()
    
    return {"message": "Conversion deleted"}


# Background task function
def perform_conversion(
    conversion_id: int,
    source_path: str,
    target_path: str,
    from_format: Optional[str],
    to_format: Optional[str],
    parameters: dict
):
    """
    Perform actual conversion in background.
    """
    from app.database import SessionLocal
    
    db = SessionLocal()
    converter = get_format_converter()
    
    try:
        # Update status to running
        conversion = db.query(FormatConversion).filter(
            FormatConversion.id == conversion_id
        ).first()
        
        if not conversion:
            return
        
        conversion.status = "running"
        db.commit()
        
        # Perform conversion
        result = converter.convert(
            source_path=source_path,
            target_path=target_path,
            from_format=from_format,
            to_format=to_format,
            **parameters
        )
        
        # Update record with results
        conversion.source_format = result['from_format']
        conversion.target_format = result['to_format']
        conversion.conversion_path = result['conversion_path']
        conversion.duration_seconds = result['duration_seconds']
        conversion.target_size_bytes = result.get('target_size_bytes')
        conversion.status = "completed"
        conversion.completed_at = datetime.utcnow()
        
        db.commit()
        
    except Exception as e:
        # Update status to failed
        conversion = db.query(FormatConversion).filter(
            FormatConversion.id == conversion_id
        ).first()
        
        if conversion:
            conversion.status = "failed"
            conversion.error_message = str(e)
            conversion.completed_at = datetime.utcnow()
            db.commit()
    
    finally:
        db.close()

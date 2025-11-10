"""
API endpoints for data export functionality.

Provides batch export with format conversion, metadata inclusion,
and progress tracking.
"""

from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from ..database import get_db
from ..services.data_export import DataExportService
from ..models.data_export import ExportStatus, ExportFormat
from ..models.user import User
from ..core.auth import get_current_user


router = APIRouter(prefix="/api/data-export", tags=["data-export"])


# Pydantic models
class ExportJobCreate(BaseModel):
    """Request model for creating export job"""

    name: str = Field(..., min_length=1, max_length=255)
    description: Optional[str] = None
    project_id: Optional[int] = None
    file_ids: List[int] = Field(..., min_items=1)
    dataset_ids: Optional[List[int]] = None
    export_format: ExportFormat = ExportFormat.CSV
    include_metadata: bool = True
    include_lineage: bool = False
    compress: bool = True
    export_options: Optional[dict] = None
    ttl_hours: int = Field(default=48, ge=1, le=168)  # 1 hour to 1 week


class ExportJobResponse(BaseModel):
    """Response model for export job"""

    id: int
    job_key: str
    name: str
    description: Optional[str]
    user_id: int
    project_id: Optional[int]
    export_format: str
    include_metadata: bool
    include_lineage: bool
    compress: bool
    file_ids: List[int]
    dataset_ids: Optional[List[int]]
    export_options: Optional[dict]
    output_path: Optional[str]
    output_size: Optional[int]
    download_url: Optional[str]
    status: str
    progress: int
    error_message: Optional[str]
    total_files: int
    processed_files: int
    failed_files: int
    created_at: str
    started_at: Optional[str]
    completed_at: Optional[str]
    expires_at: Optional[str]

    class Config:
        from_attributes = True


class ExportJobListItem(BaseModel):
    """Simplified export job for list view"""

    id: int
    job_key: str
    name: str
    user_id: int
    project_id: Optional[int]
    export_format: str
    status: str
    progress: int
    total_files: int
    processed_files: int
    output_size: Optional[int]
    created_at: str
    completed_at: Optional[str]


class ExportStatistics(BaseModel):
    """Export statistics"""

    total_jobs: int
    by_status: dict
    by_format: dict
    total_size_bytes: int
    total_files_exported: int


# Endpoints
@router.post("/jobs", response_model=ExportJobResponse)
def create_export_job(
    data: ExportJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Create a new export job.

    The job will be processed asynchronously. Use the job_key to check status.

    - **name**: Job name
    - **file_ids**: List of file IDs to export (required)
    - **export_format**: Output format (csv, tsv, json, excel, etc.)
    - **compress**: Create ZIP archive
    - **include_metadata**: Include metadata JSON file
    - **include_lineage**: Include data lineage information
    - **ttl_hours**: Time to live before auto-cleanup (default 48 hours)
    """
    service = DataExportService(db)

    job = service.create_export_job(
        user_id=current_user.id,
        name=data.name,
        file_ids=data.file_ids,
        export_format=data.export_format,
        project_id=data.project_id,
        dataset_ids=data.dataset_ids,
        description=data.description,
        include_metadata=data.include_metadata,
        include_lineage=data.include_lineage,
        compress=data.compress,
        export_options=data.export_options,
        ttl_hours=data.ttl_hours,
    )

    return ExportJobResponse(**job.to_dict())


@router.get("/jobs/{job_id}", response_model=ExportJobResponse)
def get_export_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Get export job by ID.

    Returns detailed information including status, progress, and download URL.
    """
    service = DataExportService(db)
    job = service.get_export_job(job_id, user_id=current_user.id)

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Export job not found"
        )

    return ExportJobResponse(**job.to_dict())


@router.get("/jobs/key/{job_key}", response_model=ExportJobResponse)
def get_export_job_by_key(
    job_key: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get export job by job key"""
    service = DataExportService(db)
    job = service.get_export_job_by_key(job_key, user_id=current_user.id)

    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Export job not found"
        )

    return ExportJobResponse(**job.to_dict())


@router.get("/jobs", response_model=List[ExportJobListItem])
def list_export_jobs(
    project_id: Optional[int] = None,
    status: Optional[ExportStatus] = None,
    limit: int = 50,
    offset: int = 0,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    List export jobs for current user.

    - **project_id**: Filter by project
    - **status**: Filter by status (pending, processing, completed, failed, cancelled)
    - **limit**: Max results (default 50)
    - **offset**: Pagination offset
    """
    service = DataExportService(db)
    jobs = service.list_export_jobs(
        user_id=current_user.id,
        project_id=project_id,
        status=status,
        limit=limit,
        offset=offset,
    )

    return [
        ExportJobListItem(
            id=job.id,
            job_key=job.job_key,
            name=job.name,
            user_id=job.user_id,
            project_id=job.project_id,
            export_format=job.export_format.value if job.export_format else "unknown",
            status=job.status.value if job.status else "unknown",
            progress=job.progress,
            total_files=job.total_files,
            processed_files=job.processed_files,
            output_size=job.output_size,
            created_at=job.created_at.isoformat() if job.created_at else None,
            completed_at=job.completed_at.isoformat() if job.completed_at else None,
        )
        for job in jobs
    ]


@router.post("/jobs/{job_id}/cancel")
def cancel_export_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Cancel a pending or processing export job.

    Already completed, failed, or cancelled jobs cannot be cancelled.
    """
    service = DataExportService(db)

    if not service.cancel_export_job(job_id, user_id=current_user.id):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot cancel this job (not found or already finished)",
        )

    return {"message": "Export job cancelled", "job_id": job_id}


@router.delete("/jobs/{job_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_export_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Delete an export job and its output files.

    This will permanently remove the job record and delete the exported files
    from storage.
    """
    service = DataExportService(db)

    if not service.delete_export_job(job_id, user_id=current_user.id):
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Export job not found"
        )


@router.get("/statistics", response_model=ExportStatistics)
def get_export_statistics(
    db: Session = Depends(get_db), current_user: User = Depends(get_current_user)
):
    """
    Get export statistics for current user.

    Returns counts by status and format, total size, and total files exported.
    """
    service = DataExportService(db)
    stats = service.get_export_statistics(user_id=current_user.id)

    return ExportStatistics(**stats)


@router.post("/cleanup")
def cleanup_expired_jobs(
    db: Session = Depends(get_db), current_user: User = Depends(get_current_user)
):
    """
    Clean up expired export jobs.

    This endpoint is typically called by scheduled tasks to remove
    old export files and job records.

    Note: Only admin users should have access to this endpoint in production.
    """
    service = DataExportService(db)
    count = service.cleanup_expired_jobs()

    return {"message": f"Cleaned up {count} expired export jobs", "count": count}


@router.get("/")
def get_export_info():
    """Get export service information"""
    return {
        "service": "Data Export",
        "version": "1.0",
        "features": [
            "Batch file export",
            "Multiple format support",
            "Metadata inclusion",
            "Data lineage tracking",
            "ZIP compression",
            "Progress tracking",
            "Async processing",
            "Auto-cleanup",
        ],
        "supported_formats": [f.value for f in ExportFormat],
        "endpoints": {
            "create_job": "POST /api/data-export/jobs",
            "get_job": "GET /api/data-export/jobs/{job_id}",
            "list_jobs": "GET /api/data-export/jobs",
            "cancel_job": "POST /api/data-export/jobs/{job_id}/cancel",
            "delete_job": "DELETE /api/data-export/jobs/{job_id}",
            "statistics": "GET /api/data-export/statistics",
            "cleanup": "POST /api/data-export/cleanup",
        },
    }

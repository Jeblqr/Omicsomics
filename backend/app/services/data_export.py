"""
Service layer for data export functionality.

Handles batch export operations with format conversion, metadata inclusion,
and async processing.
"""

import os
import json
import zipfile
import tempfile
import shutil
from datetime import datetime, timedelta
from typing import List, Optional, Dict, Any
from pathlib import Path
from sqlalchemy.orm import Session
from sqlalchemy import and_

from ..models.data_export import ExportJob, ExportStatus, ExportFormat
from ..models.file import File
from ..models.dataset import Dataset
from ..core.storage import storage_manager
from ..celery_app import celery_app


class DataExportService:
    """Service for managing data export operations"""

    def __init__(self, db: Session):
        self.db = db

    def create_export_job(
        self,
        user_id: int,
        name: str,
        file_ids: List[int],
        export_format: ExportFormat,
        project_id: Optional[int] = None,
        dataset_ids: Optional[List[int]] = None,
        description: Optional[str] = None,
        include_metadata: bool = True,
        include_lineage: bool = False,
        compress: bool = True,
        export_options: Optional[Dict[str, Any]] = None,
        ttl_hours: int = 48,
    ) -> ExportJob:
        """
        Create a new export job.

        Args:
            user_id: User creating the export
            name: Export job name
            file_ids: List of file IDs to export
            export_format: Output format
            project_id: Optional project ID
            dataset_ids: Optional dataset IDs to export
            description: Job description
            include_metadata: Include metadata files
            include_lineage: Include data lineage information
            compress: Compress output (ZIP)
            export_options: Format-specific options
            ttl_hours: Time to live in hours before auto-cleanup

        Returns:
            Created ExportJob
        """
        # Generate unique job key
        timestamp = datetime.utcnow().strftime("%Y%m%d%H%M%S")
        job_key = f"export_{user_id}_{timestamp}"

        # Calculate expiry
        expires_at = datetime.utcnow() + timedelta(hours=ttl_hours)

        # Create job
        job = ExportJob(
            job_key=job_key,
            name=name,
            description=description,
            user_id=user_id,
            project_id=project_id,
            export_format=export_format,
            include_metadata=include_metadata,
            include_lineage=include_lineage,
            compress=compress,
            file_ids=file_ids,
            dataset_ids=dataset_ids,
            export_options=export_options or {},
            status=ExportStatus.PENDING,
            total_files=len(file_ids),
            expires_at=expires_at,
        )

        self.db.add(job)
        self.db.commit()
        self.db.refresh(job)

        # Trigger async processing
        process_export_job.delay(job.id)

        return job

    def get_export_job(
        self, job_id: int, user_id: Optional[int] = None
    ) -> Optional[ExportJob]:
        """Get export job by ID"""
        query = self.db.query(ExportJob).filter(ExportJob.id == job_id)

        if user_id is not None:
            query = query.filter(ExportJob.user_id == user_id)

        return query.first()

    def get_export_job_by_key(
        self, job_key: str, user_id: Optional[int] = None
    ) -> Optional[ExportJob]:
        """Get export job by job key"""
        query = self.db.query(ExportJob).filter(ExportJob.job_key == job_key)

        if user_id is not None:
            query = query.filter(ExportJob.user_id == user_id)

        return query.first()

    def list_export_jobs(
        self,
        user_id: Optional[int] = None,
        project_id: Optional[int] = None,
        status: Optional[ExportStatus] = None,
        limit: int = 50,
        offset: int = 0,
    ) -> List[ExportJob]:
        """List export jobs with filters"""
        query = self.db.query(ExportJob)

        if user_id is not None:
            query = query.filter(ExportJob.user_id == user_id)

        if project_id is not None:
            query = query.filter(ExportJob.project_id == project_id)

        if status is not None:
            query = query.filter(ExportJob.status == status)

        query = query.order_by(ExportJob.created_at.desc())
        query = query.limit(limit).offset(offset)

        return query.all()

    def update_job_status(
        self,
        job_id: int,
        status: ExportStatus,
        progress: Optional[int] = None,
        error_message: Optional[str] = None,
        output_path: Optional[str] = None,
        output_size: Optional[int] = None,
        download_url: Optional[str] = None,
    ) -> Optional[ExportJob]:
        """Update export job status and progress"""
        job = self.get_export_job(job_id)

        if not job:
            return None

        job.status = status

        if progress is not None:
            job.progress = progress

        if error_message is not None:
            job.error_message = error_message

        if output_path is not None:
            job.output_path = output_path

        if output_size is not None:
            job.output_size = output_size

        if download_url is not None:
            job.download_url = download_url

        # Update timestamps
        if status == ExportStatus.PROCESSING and not job.started_at:
            job.started_at = datetime.utcnow()

        if status in [
            ExportStatus.COMPLETED,
            ExportStatus.FAILED,
            ExportStatus.CANCELLED,
        ]:
            job.completed_at = datetime.utcnow()

        self.db.commit()
        self.db.refresh(job)

        return job

    def update_job_progress(
        self, job_id: int, processed_files: int, failed_files: int = 0
    ) -> Optional[ExportJob]:
        """Update job progress counters"""
        job = self.get_export_job(job_id)

        if not job:
            return None

        job.processed_files = processed_files
        job.failed_files = failed_files

        # Calculate progress percentage
        if job.total_files > 0:
            job.progress = int((processed_files / job.total_files) * 100)

        self.db.commit()
        self.db.refresh(job)

        return job

    def cancel_export_job(self, job_id: int, user_id: Optional[int] = None) -> bool:
        """Cancel a pending or processing export job"""
        job = self.get_export_job(job_id, user_id)

        if not job:
            return False

        if job.status in [
            ExportStatus.COMPLETED,
            ExportStatus.FAILED,
            ExportStatus.CANCELLED,
        ]:
            return False  # Already finished

        job.status = ExportStatus.CANCELLED
        job.completed_at = datetime.utcnow()

        self.db.commit()

        return True

    def delete_export_job(self, job_id: int, user_id: Optional[int] = None) -> bool:
        """Delete an export job and its output files"""
        job = self.get_export_job(job_id, user_id)

        if not job:
            return False

        # Delete output file if exists
        if job.output_path:
            try:
                storage_manager.delete(job.output_path)
            except Exception as e:
                print(f"Failed to delete export output: {e}")

        # Delete job record
        self.db.delete(job)
        self.db.commit()

        return True

    def cleanup_expired_jobs(self) -> int:
        """Clean up expired export jobs (auto-maintenance)"""
        expired_jobs = (
            self.db.query(ExportJob)
            .filter(
                and_(
                    ExportJob.expires_at < datetime.utcnow(),
                    ExportJob.status.in_(
                        [
                            ExportStatus.COMPLETED,
                            ExportStatus.FAILED,
                            ExportStatus.CANCELLED,
                        ]
                    ),
                )
            )
            .all()
        )

        count = 0
        for job in expired_jobs:
            if self.delete_export_job(job.id):
                count += 1

        return count

    def get_export_statistics(self, user_id: Optional[int] = None) -> Dict[str, Any]:
        """Get export statistics"""
        query = self.db.query(ExportJob)

        if user_id is not None:
            query = query.filter(ExportJob.user_id == user_id)

        jobs = query.all()

        stats = {
            "total_jobs": len(jobs),
            "by_status": {
                "pending": len([j for j in jobs if j.status == ExportStatus.PENDING]),
                "processing": len(
                    [j for j in jobs if j.status == ExportStatus.PROCESSING]
                ),
                "completed": len(
                    [j for j in jobs if j.status == ExportStatus.COMPLETED]
                ),
                "failed": len([j for j in jobs if j.status == ExportStatus.FAILED]),
                "cancelled": len(
                    [j for j in jobs if j.status == ExportStatus.CANCELLED]
                ),
            },
            "by_format": {},
            "total_size_bytes": sum(j.output_size or 0 for j in jobs),
            "total_files_exported": sum(j.processed_files for j in jobs),
        }

        # Count by format
        for job in jobs:
            fmt = job.export_format.value if job.export_format else "unknown"
            stats["by_format"][fmt] = stats["by_format"].get(fmt, 0) + 1

        return stats


@celery_app.task(name="process_export_job")
def process_export_job(job_id: int):
    """
    Async task to process export job.

    This runs in Celery worker and performs:
    1. Fetch files from storage
    2. Apply format conversions
    3. Package with metadata
    4. Compress if requested
    5. Upload to storage
    6. Update job status
    """
    from ..database import SessionLocal

    db = SessionLocal()
    try:
        service = DataExportService(db)
        job = service.get_export_job(job_id)

        if not job:
            return {"error": "Job not found"}

        # Update to processing
        service.update_job_status(job_id, ExportStatus.PROCESSING)

        # Create temp directory for processing
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Fetch files
            files = db.query(File).filter(File.id.in_(job.file_ids)).all()

            if not files:
                service.update_job_status(
                    job_id, ExportStatus.FAILED, error_message="No files found"
                )
                return {"error": "No files found"}

            processed = 0
            failed = 0

            # Process each file
            for file in files:
                try:
                    # Download file from storage
                    file_data = storage_manager.download(file.storage_path)

                    # Save to temp directory
                    output_file = temp_path / file.filename
                    with open(output_file, "wb") as f:
                        f.write(file_data)

                    processed += 1
                    service.update_job_progress(job_id, processed, failed)

                except Exception as e:
                    print(f"Failed to export file {file.id}: {e}")
                    failed += 1
                    service.update_job_progress(job_id, processed, failed)

            # Add metadata if requested
            if job.include_metadata:
                metadata = {
                    "export_info": {
                        "job_key": job.job_key,
                        "name": job.name,
                        "created_at": job.created_at.isoformat(),
                        "format": job.export_format.value,
                    },
                    "files": [
                        {
                            "id": f.id,
                            "filename": f.filename,
                            "size": f.size,
                            "uploaded_at": (
                                f.uploaded_at.isoformat() if f.uploaded_at else None
                            ),
                        }
                        for f in files
                    ],
                }

                metadata_file = temp_path / "export_metadata.json"
                with open(metadata_file, "w") as f:
                    json.dump(metadata, f, indent=2)

            # Compress if requested
            if job.compress:
                zip_path = temp_path.parent / f"{job.job_key}.zip"

                with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
                    for file in temp_path.iterdir():
                        zipf.write(file, arcname=file.name)

                final_path = zip_path
            else:
                # For single file export
                final_path = temp_path / files[0].filename

            # Upload to storage
            with open(final_path, "rb") as f:
                file_data = f.read()

            storage_path = f"exports/{job.user_id}/{job.job_key}/{final_path.name}"
            storage_manager.upload(file_data, storage_path)

            # Generate download URL
            download_url = storage_manager.generate_presigned_url(
                storage_path, expiry=job.expires_at
            )

            # Update job as completed
            service.update_job_status(
                job_id,
                ExportStatus.COMPLETED,
                progress=100,
                output_path=storage_path,
                output_size=len(file_data),
                download_url=download_url,
            )

        return {"status": "completed", "job_id": job_id}

    except Exception as e:
        # Update job as failed
        service = DataExportService(db)
        service.update_job_status(job_id, ExportStatus.FAILED, error_message=str(e))
        return {"error": str(e)}

    finally:
        db.close()

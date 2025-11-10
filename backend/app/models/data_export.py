"""
Database models for data export functionality.

Provides batch export capabilities with format conversion, metadata inclusion,
and progress tracking.
"""

from datetime import datetime
from typing import Optional
from sqlalchemy import (
    Column,
    Integer,
    String,
    DateTime,
    Boolean,
    Text,
    JSON,
    ForeignKey,
    Enum as SQLEnum,
)
from sqlalchemy.orm import relationship
import enum

from ..database import Base


class ExportStatus(str, enum.Enum):
    """Export job status"""

    PENDING = "pending"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ExportFormat(str, enum.Enum):
    """Supported export formats"""

    CSV = "csv"
    TSV = "tsv"
    JSON = "json"
    EXCEL = "excel"
    PARQUET = "parquet"
    HDF5 = "hdf5"
    ZIP = "zip"


class ExportJob(Base):
    """
    Export job tracking table.

    Manages batch export operations with format conversion, metadata inclusion,
    and progress monitoring.
    """

    __tablename__ = "export_jobs"

    id = Column(Integer, primary_key=True, index=True)
    job_key = Column(String(100), unique=True, index=True, nullable=False)

    # Job metadata
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)

    # User and project
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False, index=True)
    project_id = Column(Integer, ForeignKey("projects.id"), nullable=True, index=True)

    # Export configuration
    export_format = Column(
        SQLEnum(ExportFormat), nullable=False, default=ExportFormat.CSV
    )
    include_metadata = Column(Boolean, default=True)
    include_lineage = Column(Boolean, default=False)
    compress = Column(Boolean, default=True)

    # File selection
    file_ids = Column(JSON, nullable=False)  # List of file IDs to export
    dataset_ids = Column(JSON, nullable=True)  # Optional: export entire datasets

    # Export options
    export_options = Column(JSON, nullable=True)  # Format-specific options

    # Output
    output_path = Column(String(500), nullable=True)  # Path in storage
    output_size = Column(Integer, nullable=True)  # Size in bytes
    download_url = Column(String(500), nullable=True)

    # Status tracking
    status = Column(SQLEnum(ExportStatus), nullable=False, default=ExportStatus.PENDING)
    progress = Column(Integer, default=0)  # 0-100
    error_message = Column(Text, nullable=True)

    # Processing info
    total_files = Column(Integer, default=0)
    processed_files = Column(Integer, default=0)
    failed_files = Column(Integer, default=0)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    expires_at = Column(DateTime, nullable=True)  # Auto-cleanup after expiry

    # Relationships
    user = relationship("User", back_populates="export_jobs")
    project = relationship("Project", back_populates="export_jobs")

    def __repr__(self):
        return (
            f"<ExportJob(id={self.id}, job_key={self.job_key}, status={self.status})>"
        )

    def to_dict(self):
        """Convert to dictionary representation"""
        return {
            "id": self.id,
            "job_key": self.job_key,
            "name": self.name,
            "description": self.description,
            "user_id": self.user_id,
            "project_id": self.project_id,
            "export_format": self.export_format.value if self.export_format else None,
            "include_metadata": self.include_metadata,
            "include_lineage": self.include_lineage,
            "compress": self.compress,
            "file_ids": self.file_ids,
            "dataset_ids": self.dataset_ids,
            "export_options": self.export_options,
            "output_path": self.output_path,
            "output_size": self.output_size,
            "download_url": self.download_url,
            "status": self.status.value if self.status else None,
            "progress": self.progress,
            "error_message": self.error_message,
            "total_files": self.total_files,
            "processed_files": self.processed_files,
            "failed_files": self.failed_files,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": (
                self.completed_at.isoformat() if self.completed_at else None
            ),
            "expires_at": self.expires_at.isoformat() if self.expires_at else None,
        }

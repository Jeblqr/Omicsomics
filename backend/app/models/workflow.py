"""Workflow execution model."""

from datetime import datetime
from enum import Enum

from sqlalchemy import DateTime, Enum as SQLEnum, ForeignKey, JSON, String, Text, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base


class WorkflowStatus(str, Enum):
    """Workflow execution status."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class Workflow(Base):
    """Workflow execution tracking."""

    __tablename__ = "workflows"

    id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), index=True)
    workflow_type: Mapped[str] = mapped_column(
        String(50)
    )  # e.g., 'nextflow', 'cwl', 'snakemake'
    status: Mapped[WorkflowStatus] = mapped_column(
        SQLEnum(WorkflowStatus), default=WorkflowStatus.PENDING, index=True
    )
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"), nullable=False)
    input_files: Mapped[dict] = mapped_column(JSON, default=dict)  # File IDs and paths
    output_files: Mapped[dict] = mapped_column(
        JSON, default=dict
    )  # Generated file paths
    parameters: Mapped[dict] = mapped_column(JSON, default=dict)  # Workflow parameters
    logs: Mapped[str] = mapped_column(Text, default="")
    error_message: Mapped[str | None] = mapped_column(Text, nullable=True)
    started_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    completed_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
    )

    sample = relationship("Sample")

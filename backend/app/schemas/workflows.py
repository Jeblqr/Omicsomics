"""Workflow schemas."""

from datetime import datetime
from typing import Any

from pydantic import BaseModel, ConfigDict

from app.models.workflow import WorkflowStatus


class WorkflowBase(BaseModel):
    """Base workflow schema."""

    name: str
    workflow_type: str
    sample_id: int
    parameters: dict[str, Any] | None = None


class WorkflowCreate(WorkflowBase):
    """Create workflow schema."""

    input_files: dict[str, Any] | None = None


class WorkflowUpdate(BaseModel):
    """Update workflow schema."""

    status: WorkflowStatus | None = None
    logs: str | None = None
    error_message: str | None = None
    output_files: dict[str, Any] | None = None


class Workflow(WorkflowBase):
    """Workflow response schema."""

    id: int
    status: WorkflowStatus
    input_files: dict[str, Any]
    output_files: dict[str, Any]
    logs: str
    error_message: str | None
    started_at: datetime | None
    completed_at: datetime | None
    created_at: datetime
    updated_at: datetime

    model_config = ConfigDict(from_attributes=True)

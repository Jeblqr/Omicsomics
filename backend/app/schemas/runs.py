"""Enhanced Run schemas with pipeline and data selection."""

from datetime import datetime
from pydantic import BaseModel, Field


class RunCreate(BaseModel):
    """Schema for creating a new run."""

    name: str = Field(..., description="Run name")
    description: str | None = Field("", description="Run description")
    project_id: int = Field(..., description="Project ID")

    # Pipeline selection
    pipeline_type: str = Field(
        ..., description="Type: 'template', 'custom', or 'merged'"
    )
    pipeline_template_id: str | None = Field(
        None, description="Template pipeline ID (e.g., 'rna-seq-basic')"
    )
    custom_pipeline_id: int | None = Field(
        None, description="Custom pipeline ID from database"
    )
    pipeline_definition: dict | None = Field(
        None, description="Inline pipeline definition for merged workflows"
    )

    # Pipeline parameters and configuration
    parameters: dict = Field(default_factory=dict, description="Pipeline parameters")
    tool_versions: dict = Field(
        default_factory=dict, description="Specific tool versions to use"
    )

    # Data file associations
    input_files: list[int] = Field(
        default_factory=list, description="List of DataFile IDs to process"
    )
    input_mapping: dict = Field(
        default_factory=dict, description="Map pipeline inputs to files"
    )

    # Resource configuration
    resources: dict | None = Field(None, description="Resource requirements override")

    # Execution options
    auto_start: bool = Field(
        False, description="Automatically start execution after creation"
    )
    priority: int = Field(0, description="Execution priority (0-10)")


class RunUpdate(BaseModel):
    """Schema for updating a run."""

    name: str | None = None
    description: str | None = None
    status: str | None = None
    parameters: dict | None = None
    started_at: datetime | None = None
    finished_at: datetime | None = None


class RunResponse(BaseModel):
    """Schema for run response."""

    id: int
    name: str
    description: str
    status: str

    # Pipeline info
    pipeline_type: str | None = None
    pipeline_template_id: str | None = None
    custom_pipeline_id: int | None = None
    pipeline_config: dict | None = None

    # Execution info
    parameters: dict = {}
    input_files: list[int] = []
    output_files: list[int] = []
    input_mapping: dict = {}
    progress: float = 0.0

    # Metadata
    project_id: int
    owner_id: int
    started_at: datetime | None = None
    finished_at: datetime | None = None
    created_at: datetime
    updated_at: datetime

    # Execution results
    logs: str | None = None
    error_message: str | None = None

    model_config = {"from_attributes": True}


class RunExecutionRequest(BaseModel):
    """Schema for requesting run execution."""

    run_id: int
    force_restart: bool = Field(False, description="Force restart if already running")
    dry_run: bool = Field(False, description="Validate without executing")


class RunExecutionStatus(BaseModel):
    """Schema for run execution status."""

    run_id: int
    status: str
    progress: float = Field(0.0, description="Progress percentage (0-100)")
    current_step: str | None = None
    completed_steps: int = 0
    failed_steps: int = 0
    logs: str | None = None
    started_at: datetime | None = None
    estimated_completion: datetime | None = None

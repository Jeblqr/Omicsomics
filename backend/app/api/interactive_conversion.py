"""
API endpoints for interactive format conversion.

Provides REST API for user-guided format conversions
that require parameter configuration.
"""

from fastapi import APIRouter, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
from pathlib import Path
import shutil
import uuid
from datetime import datetime

from app.converters.interactive_converter import (
    get_interactive_converter,
    ConversionProgress,
)

router = APIRouter(
    prefix="/api/interactive-conversion", tags=["Interactive Conversion"]
)

# Storage for conversion jobs
conversion_jobs: Dict[str, Dict[str, Any]] = {}


class ScenarioInfo(BaseModel):
    """Scenario information response."""

    id: str
    name: str
    description: str
    parameters: List[Dict[str, Any]]


class DetectionRequest(BaseModel):
    """File format detection request."""

    file_path: str


class DetectionResponse(BaseModel):
    """File format detection response."""

    detected_scenario: Optional[str] = None
    confidence: float = Field(ge=0, le=1)
    message: str


class ValidationRequest(BaseModel):
    """Validation request."""

    file_path: str
    scenario_id: str


class ValidationResponse(BaseModel):
    """Validation response."""

    valid: bool
    messages: List[Dict[str, Any]]


class PreviewRequest(BaseModel):
    """Conversion preview request."""

    file_path: str
    scenario_id: str
    parameters: Dict[str, Any]


class PreviewResponse(BaseModel):
    """Conversion preview response."""

    input_info: Dict[str, Any]
    output_info: Dict[str, Any]
    sample_data: Optional[List[Dict[str, Any]]] = None
    validation_messages: List[Dict[str, Any]]
    estimated_time: Optional[float] = None


class ConversionRequest(BaseModel):
    """Conversion execution request."""

    file_path: str
    scenario_id: str
    parameters: Dict[str, Any]
    output_filename: Optional[str] = None


class ConversionResponse(BaseModel):
    """Conversion execution response."""

    job_id: str
    status: str
    message: str


class JobStatusResponse(BaseModel):
    """Conversion job status response."""

    job_id: str
    status: str
    progress: Optional[Dict[str, Any]] = None
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


@router.get("/scenarios", response_model=List[ScenarioInfo])
async def list_scenarios():
    """
    List all available interactive conversion scenarios.

    Returns:
        List of scenario information
    """
    converter = get_interactive_converter()
    scenarios = converter.list_scenarios()

    return [
        ScenarioInfo(
            id=scenario["id"],
            name=scenario["name"],
            description=scenario.get("description", ""),
            parameters=scenario.get("parameters", []),
        )
        for scenario in scenarios
    ]


@router.post("/detect", response_model=DetectionResponse)
async def detect_scenario(request: DetectionRequest):
    """
    Detect which conversion scenario matches the input file.

    Args:
        request: Detection request with file path

    Returns:
        Detection response with scenario ID
    """
    converter = get_interactive_converter()

    # Check if file exists
    if not Path(request.file_path).exists():
        raise HTTPException(status_code=404, detail="File not found")

    # Detect scenario
    scenario_id = converter.detect_scenario(request.file_path)

    if scenario_id:
        return DetectionResponse(
            detected_scenario=scenario_id,
            confidence=0.9,  # High confidence for detected scenarios
            message=f"Detected scenario: {scenario_id}",
        )
    else:
        return DetectionResponse(
            detected_scenario=None,
            confidence=0.0,
            message="No matching scenario detected",
        )


@router.post("/validate", response_model=ValidationResponse)
async def validate_file(request: ValidationRequest):
    """
    Validate an input file for a specific scenario.

    Args:
        request: Validation request

    Returns:
        Validation response with messages
    """
    converter = get_interactive_converter()

    # Check if file exists
    if not Path(request.file_path).exists():
        raise HTTPException(status_code=404, detail="File not found")

    # Validate file
    messages = converter.validate_file(request.file_path, request.scenario_id)

    # Check if any errors
    has_errors = any(msg.level.value == "error" for msg in messages)

    return ValidationResponse(
        valid=not has_errors, messages=[msg.to_dict() for msg in messages]
    )


@router.post("/preview", response_model=PreviewResponse)
async def preview_conversion(request: PreviewRequest):
    """
    Generate a preview of conversion results.

    Args:
        request: Preview request

    Returns:
        Preview response with sample data
    """
    converter = get_interactive_converter()

    # Check if file exists
    if not Path(request.file_path).exists():
        raise HTTPException(status_code=404, detail="File not found")

    try:
        # Generate preview
        preview = converter.preview_conversion(
            request.file_path, request.scenario_id, request.parameters
        )

        return PreviewResponse(
            input_info=preview.input_info,
            output_info=preview.output_info,
            sample_data=preview.sample_data,
            validation_messages=[msg.to_dict() for msg in preview.validation_messages],
            estimated_time=preview.estimated_time,
        )

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Preview generation failed: {str(e)}"
        )


def execute_conversion(
    job_id: str,
    source_path: str,
    target_path: str,
    scenario_id: str,
    parameters: Dict[str, Any],
):
    """
    Background task to execute conversion.

    Args:
        job_id: Job identifier
        source_path: Input file path
        target_path: Output file path
        scenario_id: Scenario identifier
        parameters: Conversion parameters
    """
    converter = get_interactive_converter()

    try:
        # Update job status
        conversion_jobs[job_id]["status"] = "running"
        conversion_jobs[job_id]["started_at"] = datetime.now().isoformat()

        # Progress callback
        def progress_callback(progress: ConversionProgress):
            conversion_jobs[job_id]["progress"] = progress.to_dict()

        # Execute conversion
        result = converter.convert(
            source_path, target_path, scenario_id, parameters, progress_callback
        )

        # Update job with result
        conversion_jobs[job_id]["status"] = "completed"
        conversion_jobs[job_id]["result"] = result
        conversion_jobs[job_id]["completed_at"] = datetime.now().isoformat()

    except Exception as e:
        # Update job with error
        conversion_jobs[job_id]["status"] = "failed"
        conversion_jobs[job_id]["error"] = str(e)
        conversion_jobs[job_id]["failed_at"] = datetime.now().isoformat()


@router.post("/convert", response_model=ConversionResponse)
async def start_conversion(
    request: ConversionRequest, background_tasks: BackgroundTasks
):
    """
    Start an interactive conversion job.

    Args:
        request: Conversion request
        background_tasks: FastAPI background tasks

    Returns:
        Conversion response with job ID
    """
    # Check if file exists
    if not Path(request.file_path).exists():
        raise HTTPException(status_code=404, detail="File not found")

    # Generate job ID
    job_id = str(uuid.uuid4())

    # Determine output path
    if request.output_filename:
        output_filename = request.output_filename
    else:
        input_path = Path(request.file_path)
        output_filename = f"{input_path.stem}_standardized{input_path.suffix}"

    output_path = Path("/tmp/format_conversions") / output_filename
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create job record
    conversion_jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "source_path": request.file_path,
        "target_path": str(output_path),
        "scenario_id": request.scenario_id,
        "parameters": request.parameters,
        "created_at": datetime.now().isoformat(),
    }

    # Start background task
    background_tasks.add_task(
        execute_conversion,
        job_id,
        request.file_path,
        str(output_path),
        request.scenario_id,
        request.parameters,
    )

    return ConversionResponse(
        job_id=job_id, status="pending", message="Conversion job started"
    )


@router.get("/jobs/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """
    Get the status of a conversion job.

    Args:
        job_id: Job identifier

    Returns:
        Job status response
    """
    if job_id not in conversion_jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = conversion_jobs[job_id]

    return JobStatusResponse(
        job_id=job_id,
        status=job["status"],
        progress=job.get("progress"),
        result=job.get("result"),
        error=job.get("error"),
    )


@router.get("/jobs/{job_id}/download")
async def download_result(job_id: str):
    """
    Download the result file of a completed conversion job.

    Args:
        job_id: Job identifier

    Returns:
        File response with converted file
    """
    if job_id not in conversion_jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = conversion_jobs[job_id]

    if job["status"] != "completed":
        raise HTTPException(
            status_code=400, detail=f"Job not completed (status: {job['status']})"
        )

    output_path = job["target_path"]

    if not Path(output_path).exists():
        raise HTTPException(status_code=404, detail="Output file not found")

    return FileResponse(
        output_path,
        filename=Path(output_path).name,
        media_type="application/octet-stream",
    )


@router.delete("/jobs/{job_id}")
async def delete_job(job_id: str):
    """
    Delete a conversion job and its output file.

    Args:
        job_id: Job identifier

    Returns:
        Success message
    """
    if job_id not in conversion_jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = conversion_jobs[job_id]

    # Delete output file if exists
    if "target_path" in job:
        output_path = Path(job["target_path"])
        if output_path.exists():
            output_path.unlink()

    # Delete job record
    del conversion_jobs[job_id]

    return {"message": "Job deleted successfully"}

"""
API endpoints for data merging.

Provides REST API for merging multiple data files with
different strategies (vertical, horizontal, smart).
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
from pathlib import Path
from enum import Enum
import uuid
from datetime import datetime

from app.services.data_merger import (
    get_data_merger,
    MergeMode,
    JoinType,
    DuplicateHandling,
    MergeConfig,
)

router = APIRouter(prefix="/api/data-merger", tags=["Data Merger"])

# Storage for merge jobs
merge_jobs: Dict[str, Dict[str, Any]] = {}


class MergeModeEnum(str, Enum):
    """Merge modes for API"""

    vertical = "vertical"
    horizontal = "horizontal"
    smart = "smart"


class JoinTypeEnum(str, Enum):
    """Join types for API"""

    inner = "inner"
    left = "left"
    right = "right"
    outer = "outer"


class DuplicateHandlingEnum(str, Enum):
    """Duplicate handling for API"""

    keep_first = "keep_first"
    keep_last = "keep_last"
    keep_all = "keep_all"
    error = "error"


class MergeRequest(BaseModel):
    """Request to merge files"""

    file_paths: List[str] = Field(..., description="List of file paths to merge")
    output_filename: str = Field(
        default="merged_data.csv", description="Output filename"
    )
    mode: MergeModeEnum = Field(default=MergeModeEnum.smart, description="Merge mode")
    join_type: Optional[JoinTypeEnum] = Field(
        default=JoinTypeEnum.outer, description="Join type for horizontal merge"
    )
    key_columns: Optional[List[str]] = Field(
        default=None, description="Key columns for horizontal merge"
    )
    duplicate_handling: DuplicateHandlingEnum = Field(
        default=DuplicateHandlingEnum.keep_all,
        description="How to handle duplicates",
    )
    add_source_column: bool = Field(
        default=True, description="Add column tracking source file"
    )
    ignore_index: bool = Field(default=False, description="Reset index after merge")
    sort_keys: bool = Field(default=False, description="Sort by key columns")
    fill_missing: Optional[Any] = Field(
        default=None, description="Value to fill missing data"
    )


class PreviewRequest(BaseModel):
    """Request to preview merge"""

    file_paths: List[str]
    mode: MergeModeEnum = MergeModeEnum.smart
    join_type: Optional[JoinTypeEnum] = JoinTypeEnum.outer
    key_columns: Optional[List[str]] = None
    add_source_column: bool = True
    n_rows: int = Field(default=10, ge=1, le=100)


class MergeResponse(BaseModel):
    """Response from merge operation"""

    job_id: str
    status: str
    message: str


class JobStatusResponse(BaseModel):
    """Job status response"""

    job_id: str
    status: str
    progress: int = Field(ge=0, le=100)
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    created_at: str
    completed_at: Optional[str] = None


class PreviewResponse(BaseModel):
    """Preview response"""

    success: bool
    preview_data: Optional[Dict[str, Any]] = None
    detected_mode: Optional[str] = None
    key_columns: Optional[List[str]] = None
    estimated_rows: int = 0
    estimated_columns: int = 0
    error: Optional[str] = None


@router.post("/merge", response_model=MergeResponse)
async def merge_files(request: MergeRequest, background_tasks: BackgroundTasks):
    """
    Merge multiple files.

    Starts a background job to merge files according to configuration.
    Returns job_id for status polling.
    """
    try:
        # Generate job ID
        job_id = str(uuid.uuid4())

        # Convert paths
        file_paths = [Path(p) for p in request.file_paths]

        # Validate files exist
        for path in file_paths:
            if not path.exists():
                raise HTTPException(status_code=404, detail=f"File not found: {path}")

        # Create output path
        output_dir = Path("/tmp/merged_data")
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / f"{job_id}_{request.output_filename}"

        # Create merge config
        config = MergeConfig(
            mode=MergeMode(request.mode.value),
            join_type=JoinType(request.join_type.value) if request.join_type else None,
            key_columns=request.key_columns,
            duplicate_handling=DuplicateHandling(request.duplicate_handling.value),
            add_source_column=request.add_source_column,
            ignore_index=request.ignore_index,
            sort_keys=request.sort_keys,
            fill_missing=request.fill_missing,
        )

        # Create job record
        merge_jobs[job_id] = {
            "status": "pending",
            "progress": 0,
            "created_at": datetime.now().isoformat(),
            "file_paths": request.file_paths,
            "output_path": str(output_path),
            "config": request.dict(),
        }

        # Start background task
        background_tasks.add_task(
            execute_merge, job_id, file_paths, output_path, config
        )

        return MergeResponse(
            job_id=job_id,
            status="pending",
            message=f"Merge job started for {len(file_paths)} files",
        )

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


async def execute_merge(
    job_id: str, file_paths: List[Path], output_path: Path, config: MergeConfig
):
    """Execute merge in background"""
    try:
        merge_jobs[job_id]["status"] = "running"
        merge_jobs[job_id]["progress"] = 10

        # Get merger service
        merger = get_data_merger()

        # Perform merge
        result = merger.merge_files(file_paths, output_path, config)

        merge_jobs[job_id]["progress"] = 100

        if result.success:
            merge_jobs[job_id]["status"] = "completed"
            merge_jobs[job_id]["result"] = {
                "output_path": str(result.output_path),
                "n_rows": result.n_rows,
                "n_columns": result.n_columns,
                "n_files_merged": result.n_files_merged,
                "statistics": result.statistics,
            }
        else:
            merge_jobs[job_id]["status"] = "failed"
            merge_jobs[job_id]["error"] = ", ".join(result.warnings)

        merge_jobs[job_id]["completed_at"] = datetime.now().isoformat()

    except Exception as e:
        merge_jobs[job_id]["status"] = "failed"
        merge_jobs[job_id]["error"] = str(e)
        merge_jobs[job_id]["completed_at"] = datetime.now().isoformat()


@router.get("/jobs/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """Get status of merge job"""
    if job_id not in merge_jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = merge_jobs[job_id]

    return JobStatusResponse(
        job_id=job_id,
        status=job["status"],
        progress=job["progress"],
        result=job.get("result"),
        error=job.get("error"),
        created_at=job["created_at"],
        completed_at=job.get("completed_at"),
    )


@router.get("/jobs/{job_id}/download")
async def download_result(job_id: str):
    """Download merged file"""
    if job_id not in merge_jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = merge_jobs[job_id]

    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Job not completed")

    output_path = Path(job["output_path"])

    if not output_path.exists():
        raise HTTPException(status_code=404, detail="Output file not found")

    return FileResponse(
        path=output_path,
        filename=output_path.name,
        media_type="application/octet-stream",
    )


@router.delete("/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete merge job and output file"""
    if job_id not in merge_jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = merge_jobs[job_id]

    # Delete output file if exists
    if "output_path" in job:
        output_path = Path(job["output_path"])
        if output_path.exists():
            output_path.unlink()

    # Remove job from storage
    del merge_jobs[job_id]

    return {"message": "Job deleted successfully"}


@router.post("/preview", response_model=PreviewResponse)
async def preview_merge(request: PreviewRequest):
    """
    Preview merge operation without saving.

    Shows sample of merged data and estimated dimensions.
    """
    try:
        # Convert paths
        file_paths = [Path(p) for p in request.file_paths]

        # Validate files exist
        for path in file_paths:
            if not path.exists():
                raise HTTPException(status_code=404, detail=f"File not found: {path}")

        # Create merge config
        config = MergeConfig(
            mode=MergeMode(request.mode.value),
            join_type=JoinType(request.join_type.value) if request.join_type else None,
            key_columns=request.key_columns,
            add_source_column=request.add_source_column,
        )

        # Get merger service
        merger = get_data_merger()

        # Generate preview
        result = merger.preview_merge(file_paths, config, request.n_rows)

        if result["success"]:
            return PreviewResponse(
                success=True,
                preview_data=result["preview_data"],
                detected_mode=result.get("detected_mode"),
                key_columns=result.get("key_columns"),
                estimated_rows=result.get("estimated_rows", 0),
                estimated_columns=result.get("estimated_columns", 0),
            )
        else:
            return PreviewResponse(success=False, error=result.get("error"))

    except Exception as e:
        return PreviewResponse(success=False, error=str(e))


@router.get("/supported-formats")
async def get_supported_formats():
    """Get list of supported file formats"""
    merger = get_data_merger()
    return {
        "formats": merger.supported_formats,
        "description": "Supported file formats for data merging",
    }


@router.get("/merge-modes")
async def get_merge_modes():
    """Get available merge modes with descriptions"""
    return {
        "modes": [
            {
                "value": "vertical",
                "label": "Vertical Merge (Stack Rows)",
                "description": "Concatenate rows from multiple files. Good for appending data with same structure.",
            },
            {
                "value": "horizontal",
                "label": "Horizontal Merge (Join Columns)",
                "description": "Join files by common key columns. Good for adding features to existing records.",
            },
            {
                "value": "smart",
                "label": "Smart Merge (Auto-detect)",
                "description": "Automatically detect the best merge strategy based on data structure.",
            },
        ]
    }


@router.get("/join-types")
async def get_join_types():
    """Get available join types for horizontal merge"""
    return {
        "join_types": [
            {
                "value": "inner",
                "label": "Inner Join",
                "description": "Only keep records with matching keys in all files",
            },
            {
                "value": "left",
                "label": "Left Join",
                "description": "Keep all records from first file, add matching from others",
            },
            {
                "value": "right",
                "label": "Right Join",
                "description": "Keep all records from last file, add matching from others",
            },
            {
                "value": "outer",
                "label": "Outer Join (Full)",
                "description": "Keep all records from all files, fill missing with NaN",
            },
        ]
    }

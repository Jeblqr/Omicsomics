"""
Pipeline Dataset Integration API

Provides endpoints for integrating datasets with pipeline execution:
- Use datasets as pipeline inputs
- Auto-create datasets from pipeline outputs
- View dataset usage in pipelines
"""

from typing import Dict, List, Optional
from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from app.database import get_db
from app.services.pipeline_dataset_service import PipelineDatasetService
from app.models.run import Run

router = APIRouter(prefix="/api/pipeline-datasets")


# Pydantic Models
class DatasetFileInfo(BaseModel):
    """Information about a dataset file."""

    id: int
    path: str
    name: str
    type: str
    role: str
    size: int


class DatasetInputRequest(BaseModel):
    """Request to get dataset files for pipeline input."""

    dataset_id: int
    file_role: Optional[str] = None


class DatasetInputResponse(BaseModel):
    """Response with dataset files."""

    dataset_id: int
    files: List[DatasetFileInfo]


class CreateDatasetFromRunRequest(BaseModel):
    """Request to create dataset from run outputs."""

    run_id: int
    dataset_name: Optional[str] = None
    dataset_description: Optional[str] = None
    data_type: Optional[str] = None
    file_roles: Optional[Dict[str, str]] = Field(
        default=None, description="Map file IDs to roles"
    )
    tags: Optional[List[str]] = None


class DatasetInfo(BaseModel):
    """Basic dataset information."""

    id: int
    name: str
    description: str
    data_type: str
    status: str
    file_count: int
    created_at: str


class CreateDatasetFromRunResponse(BaseModel):
    """Response with created dataset."""

    dataset: DatasetInfo
    message: str


class LinkDatasetsToRunRequest(BaseModel):
    """Request to link input datasets to run."""

    run_id: int
    dataset_ids: List[int]
    input_mapping: Optional[Dict] = Field(
        default=None, description="Map node inputs to dataset files"
    )


class RunInfo(BaseModel):
    """Pipeline run information."""

    run_id: int
    run_name: str
    status: str
    started_at: Optional[str]
    finished_at: Optional[str]


class DatasetUsageResponse(BaseModel):
    """Response with dataset usage information."""

    dataset_id: int
    runs: List[RunInfo]


class ValidationResult(BaseModel):
    """Dataset validation result."""

    is_valid: bool
    errors: List[str]


class ValidateDatasetRequest(BaseModel):
    """Request to validate dataset for pipeline input."""

    dataset_id: int
    required_file_types: Optional[List[str]] = None


class AutoCreateDatasetRequest(BaseModel):
    """Request to auto-create dataset from run."""

    run_id: int
    auto_tags: bool = True


# API Endpoints


@router.post("/input-files", response_model=DatasetInputResponse)
def get_dataset_input_files(
    request: DatasetInputRequest,
    db: Session = Depends(get_db),
):
    """
    Get files from a dataset for pipeline input.

    Returns file information including paths, names, types, and roles.
    """
    service = PipelineDatasetService(db)

    try:
        files = service.get_dataset_files_for_input(
            dataset_id=request.dataset_id, file_role=request.file_role
        )

        return DatasetInputResponse(
            dataset_id=request.dataset_id,
            files=[DatasetFileInfo(**f) for f in files],
        )
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.post("/create-from-run", response_model=CreateDatasetFromRunResponse)
def create_dataset_from_run(
    request: CreateDatasetFromRunRequest,
    db: Session = Depends(get_db),
):
    """
    Create a dataset from pipeline run outputs.

    Automatically adds all output files from the run to the new dataset,
    records lineage, and applies tags.
    """
    service = PipelineDatasetService(db)

    # Get run
    run = db.query(Run).filter(Run.id == request.run_id).first()
    if not run:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Run {request.run_id} not found",
        )

    try:
        dataset = service.create_dataset_from_run(
            run=run,
            dataset_name=request.dataset_name,
            dataset_description=request.dataset_description,
            data_type=request.data_type,
            file_roles=request.file_roles,
            tags=request.tags,
        )

        # Record lineage
        service.record_pipeline_lineage(dataset_id=dataset.id, run=run)

        return CreateDatasetFromRunResponse(
            dataset=DatasetInfo(
                id=dataset.id,
                name=dataset.name,
                description=dataset.description,
                data_type=dataset.data_type,
                status=dataset.status,
                file_count=len(dataset.files),
                created_at=dataset.created_at.isoformat(),
            ),
            message=f"Dataset created successfully with {len(dataset.files)} files",
        )
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.post("/link-to-run")
def link_datasets_to_run(
    request: LinkDatasetsToRunRequest,
    db: Session = Depends(get_db),
):
    """
    Link input datasets to a pipeline run.

    Stores dataset references in run metadata for tracking and lineage.
    """
    service = PipelineDatasetService(db)

    try:
        run = service.link_input_datasets_to_run(
            run_id=request.run_id,
            dataset_ids=request.dataset_ids,
            input_mapping=request.input_mapping,
        )

        return {
            "message": f"Linked {len(request.dataset_ids)} datasets to run {run.id}",
            "run_id": run.id,
            "dataset_ids": request.dataset_ids,
        }
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.get("/run/{run_id}/output-datasets", response_model=List[DatasetInfo])
def get_run_output_datasets(
    run_id: int,
    db: Session = Depends(get_db),
):
    """
    Get all datasets created from a pipeline run's outputs.
    """
    service = PipelineDatasetService(db)

    try:
        datasets = service.get_run_output_datasets(run_id)

        return [
            DatasetInfo(
                id=ds.id,
                name=ds.name,
                description=ds.description,
                data_type=ds.data_type,
                status=ds.status,
                file_count=len(ds.files),
                created_at=ds.created_at.isoformat(),
            )
            for ds in datasets
        ]
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.get("/dataset/{dataset_id}/usage", response_model=DatasetUsageResponse)
def get_dataset_usage_in_runs(
    dataset_id: int,
    db: Session = Depends(get_db),
):
    """
    Get all pipeline runs that used this dataset as input.
    """
    service = PipelineDatasetService(db)

    try:
        runs = service.get_dataset_usage_in_runs(dataset_id)

        return DatasetUsageResponse(
            dataset_id=dataset_id, runs=[RunInfo(**r) for r in runs]
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.post("/validate", response_model=ValidationResult)
def validate_dataset_for_pipeline(
    request: ValidateDatasetRequest,
    db: Session = Depends(get_db),
):
    """
    Validate that a dataset is suitable for pipeline input.

    Checks dataset status, file existence, and required file types.
    """
    service = PipelineDatasetService(db)

    try:
        is_valid, errors = service.validate_dataset_for_pipeline_input(
            dataset_id=request.dataset_id,
            required_file_types=request.required_file_types,
        )

        return ValidationResult(is_valid=is_valid, errors=errors)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.post("/auto-create", response_model=CreateDatasetFromRunResponse)
def auto_create_output_dataset(
    request: AutoCreateDatasetRequest,
    db: Session = Depends(get_db),
):
    """
    Automatically create output dataset for a completed pipeline run.

    Infers data type, adds automatic tags, and records lineage.
    """
    service = PipelineDatasetService(db)

    # Get run
    run = db.query(Run).filter(Run.id == request.run_id).first()
    if not run:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Run {request.run_id} not found",
        )

    try:
        dataset = service.auto_create_output_dataset(
            run=run, auto_tags=request.auto_tags
        )

        if not dataset:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Failed to create dataset (run may not be completed or has no outputs)",
            )

        return CreateDatasetFromRunResponse(
            dataset=DatasetInfo(
                id=dataset.id,
                name=dataset.name,
                description=dataset.description,
                data_type=dataset.data_type,
                status=dataset.status,
                file_count=len(dataset.files),
                created_at=dataset.created_at.isoformat(),
            ),
            message=f"Dataset auto-created successfully with {len(dataset.files)} files",
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e)
        )


@router.get("/")
def pipeline_datasets_info():
    """Get information about pipeline dataset integration."""
    return {
        "service": "Pipeline Dataset Integration",
        "version": "1.0",
        "description": "Integrate datasets with pipeline execution for seamless data flow",
        "features": [
            "Use datasets as pipeline inputs",
            "Auto-create datasets from pipeline outputs",
            "Record lineage for pipeline operations",
            "Track dataset usage across runs",
            "Validate datasets for pipeline compatibility",
        ],
        "endpoints": {
            "POST /input-files": "Get dataset files for pipeline input",
            "POST /create-from-run": "Create dataset from run outputs",
            "POST /link-to-run": "Link input datasets to run",
            "GET /run/{run_id}/output-datasets": "Get datasets created by run",
            "GET /dataset/{dataset_id}/usage": "Get runs that used dataset",
            "POST /validate": "Validate dataset for pipeline input",
            "POST /auto-create": "Auto-create dataset from run",
        },
    }

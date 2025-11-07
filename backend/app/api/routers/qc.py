"""Quality Control (QC) API endpoints."""

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.schemas import workflows as workflow_schema
from app.services import files as file_service
from app.services import samples as sample_service
from app.services import workflows as workflow_service
from app.workflows.executor import workflow_executor

router = APIRouter()


class QCRequest(BaseModel):
    """Request to run QC on files."""

    sample_id: int
    file_ids: list[int]


class QCResponse(BaseModel):
    """QC workflow response."""

    workflow_id: int
    status: str
    message: str


@router.post("/fastqc", response_model=QCResponse)
async def run_fastqc(
    qc_request: QCRequest,
    background_tasks: BackgroundTasks,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> QCResponse:
    """
    Run FastQC on the specified files.

    Args:
        qc_request: QC request with sample_id and file_ids
        background_tasks: FastAPI background tasks
        current_user: Authenticated user
        db: Database session

    Returns:
        QC workflow information
    """
    # Verify sample exists
    sample = await sample_service.get_sample(db, qc_request.sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    # Verify files exist and belong to the sample
    input_files = []
    for file_id in qc_request.file_ids:
        file = await file_service.get_file(db, file_id)
        if file is None:
            raise HTTPException(status_code=404, detail=f"File {file_id} not found")
        if file.sample_id != qc_request.sample_id:
            raise HTTPException(
                status_code=400,
                detail=f"File {file_id} does not belong to sample {qc_request.sample_id}",
            )
        # For now, use placeholder paths - in production, download from S3
        input_files.append(f"/tmp/file_{file_id}.fastq")

    # Create workflow record
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"FastQC for Sample {qc_request.sample_id}",
            workflow_type="fastqc",
            sample_id=qc_request.sample_id,
            input_files={"file_ids": qc_request.file_ids},
            parameters={
                "input_files": input_files,
                "output_dir": f"/tmp/qc_{qc_request.sample_id}",
            },
        ),
    )

    # Schedule FastQC execution
    background_tasks.add_task(
        workflow_executor.execute_fastqc,
        workflow.id,
        input_files,
        f"/tmp/qc_{qc_request.sample_id}",
        db,
    )

    return QCResponse(
        workflow_id=workflow.id,
        status="pending",
        message="FastQC analysis has been queued",
    )


@router.get("/results/{workflow_id}")
async def get_qc_results(
    workflow_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Get QC results for a workflow.

    Args:
        workflow_id: Workflow ID
        current_user: Authenticated user
        db: Database session

    Returns:
        QC results including status, logs, and output files
    """
    workflow = await workflow_service.get_workflow(db, workflow_id)
    if workflow is None:
        raise HTTPException(status_code=404, detail="Workflow not found")

    if workflow.workflow_type != "fastqc":
        raise HTTPException(status_code=400, detail="Workflow is not a FastQC workflow")

    return {
        "workflow_id": workflow.id,
        "status": workflow.status.value,
        "started_at": workflow.started_at,
        "completed_at": workflow.completed_at,
        "logs": workflow.logs,
        "output_files": workflow.output_files,
        "error_message": workflow.error_message,
    }

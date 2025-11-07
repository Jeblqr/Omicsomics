"""Workflow API endpoints."""

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import samples as sample_service
from app.services import workflows as workflow_service
from app.workflows.executor import workflow_executor

router = APIRouter()


@router.get("/", response_model=list[workflow_schema.Workflow])
async def list_workflows(
    sample_id: int | None = None,
    status: WorkflowStatus | None = None,
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> list[workflow_schema.Workflow]:
    """List workflows, optionally filtered by sample and status."""
    workflows = await workflow_service.get_workflows(
        db, skip=skip, limit=limit, sample_id=sample_id, status=status
    )
    return workflows


@router.post(
    "/", response_model=workflow_schema.Workflow, status_code=status.HTTP_201_CREATED
)
async def create_workflow(
    workflow: workflow_schema.WorkflowCreate,
    background_tasks: BackgroundTasks,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> workflow_schema.Workflow:
    """Create and optionally execute a workflow."""
    # Verify sample exists and user has access
    sample = await sample_service.get_sample(db, workflow.sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    # Create workflow record
    db_workflow = await workflow_service.create_workflow(db, workflow)

    # Schedule execution in background if parameters provided
    if workflow.parameters:
        if workflow.workflow_type == "nextflow":
            pipeline = workflow.parameters.get("pipeline")
            if pipeline:
                background_tasks.add_task(
                    workflow_executor.execute_nextflow,
                    db_workflow.id,
                    pipeline,
                    workflow.parameters,
                    db,
                )
        elif workflow.workflow_type == "fastqc":
            input_files = workflow.parameters.get("input_files", [])
            output_dir = workflow.parameters.get(
                "output_dir", f"/tmp/qc_{db_workflow.id}"
            )
            if input_files:
                background_tasks.add_task(
                    workflow_executor.execute_fastqc,
                    db_workflow.id,
                    input_files,
                    output_dir,
                    db,
                )

    return db_workflow


@router.get("/{workflow_id}", response_model=workflow_schema.Workflow)
async def get_workflow(
    workflow_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> workflow_schema.Workflow:
    """Get a specific workflow by ID."""
    workflow = await workflow_service.get_workflow(db, workflow_id)
    if workflow is None:
        raise HTTPException(status_code=404, detail="Workflow not found")
    return workflow


@router.put("/{workflow_id}", response_model=workflow_schema.Workflow)
async def update_workflow(
    workflow_id: int,
    workflow_update: workflow_schema.WorkflowUpdate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> workflow_schema.Workflow:
    """Update a workflow (e.g., cancel it)."""
    workflow = await workflow_service.get_workflow(db, workflow_id)
    if workflow is None:
        raise HTTPException(status_code=404, detail="Workflow not found")

    updated_workflow = await workflow_service.update_workflow(
        db, workflow_id, workflow_update
    )
    return updated_workflow


@router.delete("/{workflow_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_workflow(
    workflow_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Delete a workflow."""
    workflow = await workflow_service.get_workflow(db, workflow_id)
    if workflow is None:
        raise HTTPException(status_code=404, detail="Workflow not found")

    await workflow_service.delete_workflow(db, workflow_id)
    return None

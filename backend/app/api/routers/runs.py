from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.services import runs as runs_service
from app.services import projects as project_service
from app.schemas.runs import (
    RunCreate,
    RunUpdate,
    RunResponse,
    RunExecutionRequest,
    RunExecutionStatus,
)

router = APIRouter()


@router.post("/", status_code=status.HTTP_201_CREATED, response_model=RunResponse)
async def create_run(
    run_in: RunCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    project = await project_service.get_project(db, run_in.project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    run = await runs_service.create_run(
        db,
        run_data=run_in,
        owner_id=current_user.id,
    )

    return RunResponse(
        id=run.id,
        name=run.name,
        description=run.description,
        status=run.status,
        pipeline_type=run.pipeline_type,
        pipeline_template_id=run.pipeline_template_id,
        custom_pipeline_id=run.custom_pipeline_id,
        pipeline_config=run.pipeline_config or {},
        parameters=run.parameters or {},
        input_files=run.input_files or [],
        output_files=run.output_files or [],
        input_mapping=run.input_mapping or {},
        logs=run.logs,
        error_message=run.error_message,
        progress=run.progress or 0.0,
        project_id=run.project_id,
        owner_id=run.owner_id,
        started_at=run.started_at,
        finished_at=run.finished_at,
        created_at=run.created_at,
        updated_at=run.updated_at,
    )


@router.get("/", response_model=list[RunResponse])
async def list_runs(
    project_id: int | None = None,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    if project_id is not None:
        project = await project_service.get_project(db, project_id)
        if project is None:
            raise HTTPException(status_code=404, detail="Project not found")
        if project.owner_id != current_user.id:
            raise HTTPException(status_code=403, detail="Not authorized")

    rows = await runs_service.list_runs(db, project_id)
    return [
        RunResponse(
            id=r.id,
            name=r.name,
            description=r.description,
            status=r.status,
            pipeline_type=r.pipeline_type,
            pipeline_template_id=r.pipeline_template_id,
            custom_pipeline_id=r.custom_pipeline_id,
            pipeline_config=r.pipeline_config or {},
            parameters=r.parameters or {},
            input_files=r.input_files or [],
            output_files=r.output_files or [],
            input_mapping=r.input_mapping or {},
            logs=r.logs,
            error_message=r.error_message,
            progress=r.progress or 0.0,
            project_id=r.project_id,
            owner_id=r.owner_id,
            started_at=r.started_at,
            finished_at=r.finished_at,
            created_at=r.created_at,
            updated_at=r.updated_at,
        )
        for r in rows
    ]


@router.get("/{run_id}", response_model=RunResponse)
async def get_run(
    run_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    return RunResponse(
        id=run.id,
        name=run.name,
        description=run.description,
        status=run.status,
        pipeline_type=run.pipeline_type,
        pipeline_template_id=run.pipeline_template_id,
        custom_pipeline_id=run.custom_pipeline_id,
        pipeline_config=run.pipeline_config or {},
        parameters=run.parameters or {},
        input_files=run.input_files or [],
        output_files=run.output_files or [],
        input_mapping=run.input_mapping or {},
        logs=run.logs,
        error_message=run.error_message,
        progress=run.progress or 0.0,
        project_id=run.project_id,
        owner_id=run.owner_id,
        started_at=run.started_at,
        finished_at=run.finished_at,
        created_at=run.created_at,
        updated_at=run.updated_at,
    )


@router.patch("/{run_id}", response_model=RunResponse)
async def update_run(
    run_id: int,
    run_update: RunUpdate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Update a run's information or status."""
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Update fields if provided
    update_data = run_update.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(run, field, value)

    await db.commit()
    await db.refresh(run)

    return RunResponse(
        id=run.id,
        name=run.name,
        description=run.description,
        status=run.status,
        pipeline_type=run.pipeline_type,
        pipeline_template_id=run.pipeline_template_id,
        custom_pipeline_id=run.custom_pipeline_id,
        pipeline_config=run.pipeline_config or {},
        parameters=run.parameters or {},
        input_files=run.input_files or [],
        output_files=run.output_files or [],
        input_mapping=run.input_mapping or {},
        logs=run.logs,
        error_message=run.error_message,
        progress=run.progress or 0.0,
        project_id=run.project_id,
        owner_id=run.owner_id,
        started_at=run.started_at,
        finished_at=run.finished_at,
        created_at=run.created_at,
        updated_at=run.updated_at,
    )


@router.post("/{run_id}/start")
async def start_run(
    run_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Start execution of a run."""
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    if run.status == "running":
        raise HTTPException(status_code=400, detail="Run is already running")

    # Start pipeline execution in background
    import asyncio
    from app.services.pipeline_executor import execute_run_async

    # Create a new task to execute the run in background (with its own db session)
    asyncio.create_task(execute_run_async(run_id))

    return {"message": "Run started", "run_id": run_id, "status": "running"}


@router.post("/{run_id}/stop")
async def stop_run(
    run_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Stop execution of a running run."""
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    if run.status != "running":
        raise HTTPException(status_code=400, detail="Run is not running")

    run.status = "cancelled"
    from datetime import datetime, timezone

    run.finished_at = datetime.now(timezone.utc)
    await db.commit()

    # TODO: Integrate with execution engine to actually stop

    return {"message": "Run stopped", "run_id": run_id, "status": "cancelled"}


@router.get("/{run_id}/status", response_model=RunExecutionStatus)
async def get_run_status(
    run_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Get the execution status and progress of a run."""
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # TODO: Get detailed step information from execution engine
    return RunExecutionStatus(
        run_id=run.id,
        status=run.status,
        progress=run.progress or 0.0,
        current_step=None,
        completed_steps=0,
        failed_steps=0,
        logs=run.logs,
        started_at=run.started_at,
        estimated_completion=None,
    )


@router.get("/{run_id}/logs")
async def get_run_logs(
    run_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Get the execution logs of a run."""
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    return {
        "run_id": run.id,
        "logs": run.logs or "",
        "error_message": run.error_message,
    }


@router.delete("/{run_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_run(
    run_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Delete a run."""
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    await runs_service.delete_run(db, run_id)
    return None

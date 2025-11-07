"""Workflow service for managing workflow executions."""

from datetime import datetime

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.workflow import Workflow, WorkflowStatus
from app.schemas import workflows as workflow_schema


async def get_workflow(db: AsyncSession, workflow_id: int) -> Workflow | None:
    """Get a workflow by ID."""
    result = await db.execute(select(Workflow).where(Workflow.id == workflow_id))
    return result.scalars().first()


async def get_workflows(
    db: AsyncSession,
    skip: int = 0,
    limit: int = 100,
    sample_id: int | None = None,
    status: WorkflowStatus | None = None,
) -> list[Workflow]:
    """Get workflows, optionally filtered by sample and status."""
    query = select(Workflow).order_by(Workflow.created_at.desc())
    if sample_id is not None:
        query = query.where(Workflow.sample_id == sample_id)
    if status is not None:
        query = query.where(Workflow.status == status)
    query = query.offset(skip).limit(limit)
    result = await db.execute(query)
    return list(result.scalars().all())


async def create_workflow(
    db: AsyncSession, workflow: workflow_schema.WorkflowCreate
) -> Workflow:
    """Create a new workflow."""
    db_workflow = Workflow(
        name=workflow.name,
        workflow_type=workflow.workflow_type,
        sample_id=workflow.sample_id,
        input_files=workflow.input_files or {},
        parameters=workflow.parameters or {},
        status=WorkflowStatus.PENDING,
    )
    db.add(db_workflow)
    await db.commit()
    await db.refresh(db_workflow)
    return db_workflow


async def update_workflow(
    db: AsyncSession, workflow_id: int, workflow: workflow_schema.WorkflowUpdate
) -> Workflow | None:
    """Update a workflow."""
    db_workflow = await get_workflow(db, workflow_id)
    if db_workflow is None:
        return None

    update_data = workflow.model_dump(exclude_unset=True)
    
    # Handle status transitions
    if "status" in update_data:
        new_status = update_data["status"]
        if new_status == WorkflowStatus.RUNNING and db_workflow.started_at is None:
            db_workflow.started_at = datetime.utcnow()
        elif new_status in (WorkflowStatus.COMPLETED, WorkflowStatus.FAILED, WorkflowStatus.CANCELLED):
            if db_workflow.completed_at is None:
                db_workflow.completed_at = datetime.utcnow()
    
    for field, value in update_data.items():
        if field != "status" or value is not None:
            setattr(db_workflow, field, value)

    await db.commit()
    await db.refresh(db_workflow)
    return db_workflow


async def delete_workflow(db: AsyncSession, workflow_id: int) -> bool:
    """Delete a workflow."""
    db_workflow = await get_workflow(db, workflow_id)
    if db_workflow is None:
        return False

    await db.delete(db_workflow)
    await db.commit()
    return True

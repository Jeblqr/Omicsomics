from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from app.models.run import Run
from app.schemas.runs import RunCreate


async def create_run(
    db: AsyncSession,
    run_data: RunCreate,
    owner_id: int,
) -> Run:
    """Create a new run with pipeline configuration."""

    # Build pipeline_config from template or custom pipeline
    pipeline_config = {}
    if run_data.pipeline_type == "template" and run_data.pipeline_template_id:
        pipeline_config["type"] = "template"
        pipeline_config["template_id"] = run_data.pipeline_template_id
    elif run_data.pipeline_type == "custom" and run_data.custom_pipeline_id:
        pipeline_config["type"] = "custom"
        pipeline_config["custom_id"] = run_data.custom_pipeline_id
    elif run_data.pipeline_type == "merged":
        pipeline_config["type"] = "merged"
        if run_data.pipeline_template_id:
            pipeline_config["template_id"] = run_data.pipeline_template_id
        if run_data.custom_pipeline_id:
            pipeline_config["custom_id"] = run_data.custom_pipeline_id

    # Store pipeline definition if provided
    if run_data.pipeline_definition:
        pipeline_config["definition"] = run_data.pipeline_definition

    # Store tool versions if provided
    if run_data.tool_versions:
        pipeline_config["tool_versions"] = run_data.tool_versions

    run = Run(
        name=run_data.name,
        description=run_data.description or "",
        project_id=run_data.project_id,
        owner_id=owner_id,
        pipeline_type=run_data.pipeline_type,
        pipeline_template_id=run_data.pipeline_template_id,
        custom_pipeline_id=run_data.custom_pipeline_id,
        pipeline_config=pipeline_config,
        parameters=run_data.parameters or {},
        input_files=run_data.input_files or [],
        input_mapping=run_data.input_mapping or {},
        progress=0.0,
    )
    db.add(run)
    await db.commit()
    await db.refresh(run)
    return run


async def get_run(db: AsyncSession, run_id: int) -> Run | None:
    r = await db.execute(select(Run).where(Run.id == run_id))
    return r.scalars().first()


async def list_runs(
    db: AsyncSession, project_id: int | None = None, skip: int = 0, limit: int = 100
):
    q = select(Run).order_by(Run.created_at.desc())
    if project_id is not None:
        q = q.where(Run.project_id == project_id)
    q = q.offset(skip).limit(limit)
    r = await db.execute(q)
    return list(r.scalars().all())


async def delete_run(db: AsyncSession, run_id: int) -> None:
    """Delete a run."""
    run = await get_run(db, run_id)
    if run:
        await db.delete(run)
        await db.commit()

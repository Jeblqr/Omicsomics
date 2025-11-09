from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from app.models.run import Run


async def create_run(
    db: AsyncSession,
    name: str,
    description: str,
    project_id: int,
    owner_id: int,
    pipeline_config: dict | None = None,
) -> Run:
    run = Run(
        name=name,
        description=description,
        project_id=project_id,
        owner_id=owner_id,
        pipeline_config=pipeline_config,
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

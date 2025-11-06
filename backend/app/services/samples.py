from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.sample import Sample
from app.schemas import samples as sample_schema


async def get_sample(db: AsyncSession, sample_id: int) -> Sample | None:
    """Get a sample by ID."""
    result = await db.execute(select(Sample).where(Sample.id == sample_id))
    return result.scalars().first()


async def get_samples(
    db: AsyncSession, skip: int = 0, limit: int = 100, project_id: int | None = None
) -> list[Sample]:
    """Get samples, optionally filtered by project."""
    query = select(Sample).order_by(Sample.created_at.desc())
    if project_id is not None:
        query = query.where(Sample.project_id == project_id)
    query = query.offset(skip).limit(limit)
    result = await db.execute(query)
    return list(result.scalars().all())


async def create_sample(
    db: AsyncSession, sample: sample_schema.SampleCreate
) -> Sample:
    """Create a new sample."""
    db_sample = Sample(
        name=sample.name,
        description=sample.description or "",
        project_id=sample.project_id,
        metadata_=sample.metadata_ or {},
    )
    db.add(db_sample)
    await db.commit()
    await db.refresh(db_sample)
    return db_sample


async def update_sample(
    db: AsyncSession, sample_id: int, sample: sample_schema.SampleUpdate
) -> Sample | None:
    """Update a sample."""
    db_sample = await get_sample(db, sample_id)
    if db_sample is None:
        return None

    update_data = sample.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(db_sample, field, value)

    await db.commit()
    await db.refresh(db_sample)
    return db_sample


async def delete_sample(db: AsyncSession, sample_id: int) -> bool:
    """Delete a sample."""
    db_sample = await get_sample(db, sample_id)
    if db_sample is None:
        return False

    await db.delete(db_sample)
    await db.commit()
    return True

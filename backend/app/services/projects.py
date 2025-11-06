from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.project import Project
from app.schemas import projects as project_schema


async def get_project(db: AsyncSession, project_id: int) -> Project | None:
    """Get a project by ID."""
    result = await db.execute(select(Project).where(Project.id == project_id))
    return result.scalars().first()


async def get_projects(
    db: AsyncSession, skip: int = 0, limit: int = 100, owner_id: int | None = None
) -> list[Project]:
    """Get projects, optionally filtered by owner."""
    query = select(Project).order_by(Project.created_at.desc())
    if owner_id is not None:
        query = query.where(Project.owner_id == owner_id)
    query = query.offset(skip).limit(limit)
    result = await db.execute(query)
    return list(result.scalars().all())


async def create_project(
    db: AsyncSession, project: project_schema.ProjectCreate, owner_id: int
) -> Project:
    """Create a new project."""
    db_project = Project(
        name=project.name,
        description=project.description or "",
        owner_id=owner_id,
    )
    db.add(db_project)
    await db.commit()
    await db.refresh(db_project)
    return db_project


async def update_project(
    db: AsyncSession, project_id: int, project: project_schema.ProjectUpdate
) -> Project | None:
    """Update a project."""
    db_project = await get_project(db, project_id)
    if db_project is None:
        return None

    update_data = project.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(db_project, field, value)

    await db.commit()
    await db.refresh(db_project)
    return db_project


async def delete_project(db: AsyncSession, project_id: int) -> bool:
    """Delete a project."""
    db_project = await get_project(db, project_id)
    if db_project is None:
        return False

    await db.delete(db_project)
    await db.commit()
    return True

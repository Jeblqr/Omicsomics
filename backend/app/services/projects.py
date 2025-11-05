from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.project import Project
from app.schemas.projects import ProjectRead


async def list_projects(session: AsyncSession) -> list[ProjectRead]:
  result = await session.execute(select(Project))
  projects = result.scalars().all()
  return [ProjectRead.model_validate(project) for project in projects]

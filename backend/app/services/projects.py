from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.project import Project
from app.schemas.projects import ProjectCreate, ProjectRead


async def create_project(session: AsyncSession, payload: ProjectCreate) -> ProjectRead:
  project = Project(name=payload.name, description=payload.description or "")
  session.add(project)
  await session.commit()
  await session.refresh(project)
  return ProjectRead.model_validate(project)


async def list_projects(session: AsyncSession) -> list[ProjectRead]:
  result = await session.execute(select(Project).order_by(Project.id))
  projects = result.scalars().all()
  return [ProjectRead.model_validate(project) for project in projects]

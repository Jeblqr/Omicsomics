from sqlalchemy.ext.asyncio import AsyncSession

from fastapi import APIRouter, Depends, status

from app.api.dependencies.database import get_db
from app.schemas.projects import ProjectCreate, ProjectRead
from app.services.projects import create_project as create_project_service
from app.services.projects import list_projects as list_projects_service

router = APIRouter()


@router.get("/", response_model=list[ProjectRead], summary="List projects")
async def list_projects(session: AsyncSession = Depends(get_db)) -> list[ProjectRead]:
  """List registered projects."""
  return await list_projects_service(session)


@router.post("/", response_model=ProjectRead, status_code=status.HTTP_201_CREATED)
async def create_project(
  payload: ProjectCreate,
  session: AsyncSession = Depends(get_db),
) -> ProjectRead:
  """Create a new project for tracking samples and workflows."""
  return await create_project_service(session, payload)

from fastapi import APIRouter, Depends

from app.api.dependencies.database import get_db
from app.schemas.projects import ProjectRead
from app.services.projects import list_projects as list_projects_service

router = APIRouter()


@router.get("/", response_model=list[ProjectRead], summary="List projects")
async def list_projects(db=Depends(get_db)) -> list[ProjectRead]:
  """List registered projects."""
  return await list_projects_service(db)

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.schemas import projects as project_schema
from app.services import projects as project_service

router = APIRouter()


@router.get("/", response_model=list[project_schema.Project], summary="List projects")
async def list_projects(
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> list[project_schema.Project]:
    """List projects owned by the current user."""
    projects = await project_service.get_projects(
        db, skip=skip, limit=limit, owner_id=current_user.id
    )
    return projects


@router.post(
    "/", response_model=project_schema.Project, status_code=status.HTTP_200_OK
)
async def create_project(
    project: project_schema.ProjectCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> project_schema.Project:
    """Create a new project for tracking samples and workflows."""
    return await project_service.create_project(db, project, current_user.id)


@router.get("/{project_id}", response_model=project_schema.Project)
async def get_project(
    project_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> project_schema.Project:
    """Get a specific project by ID."""
    project = await project_service.get_project(db, project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=404, detail="Project not found")
    return project


@router.put("/{project_id}", response_model=project_schema.Project)
async def update_project(
    project_id: int,
    project_update: project_schema.ProjectUpdate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> project_schema.Project:
    """Update a project."""
    project = await project_service.get_project(db, project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=404, detail="Project not found")
    
    updated_project = await project_service.update_project(db, project_id, project_update)
    return updated_project


@router.delete("/{project_id}", status_code=status.HTTP_200_OK)
async def delete_project(
    project_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Delete a project."""
    project = await project_service.get_project(db, project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=404, detail="Project not found")
    
    await project_service.delete_project(db, project_id)
    return {"detail": "Project deleted"}

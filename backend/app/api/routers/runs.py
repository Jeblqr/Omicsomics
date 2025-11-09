from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.services import runs as runs_service
from app.services import projects as project_service

router = APIRouter()


class RunCreate(BaseModel):
    name: str
    description: str | None = ""
    project_id: int


@router.post("/", status_code=status.HTTP_201_CREATED)
async def create_run(
    run_in: RunCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    project = await project_service.get_project(db, run_in.project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    run = await runs_service.create_run(
        db, run_in.name, run_in.description or "", run_in.project_id, current_user.id
    )
    return {"id": run.id, "name": run.name, "status": run.status}


@router.get("/")
async def list_runs(
    project_id: int | None = None,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    if project_id is not None:
        project = await project_service.get_project(db, project_id)
        if project is None:
            raise HTTPException(status_code=404, detail="Project not found")
        if project.owner_id != current_user.id:
            raise HTTPException(status_code=403, detail="Not authorized")

    rows = await runs_service.list_runs(db, project_id)
    return [
        {"id": r.id, "name": r.name, "status": r.status, "created_at": r.created_at}
        for r in rows
    ]


@router.get("/{run_id}")
async def get_run(
    run_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    run = await runs_service.get_run(db, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    project = await project_service.get_project(db, run.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")
    return {
        "id": run.id,
        "name": run.name,
        "status": run.status,
        "description": run.description,
    }

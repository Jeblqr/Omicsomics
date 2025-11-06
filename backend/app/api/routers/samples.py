from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.schemas import samples as sample_schema
from app.services import projects as project_service
from app.services import samples as sample_service

router = APIRouter()


@router.get("/", response_model=list[sample_schema.Sample], summary="List samples")
async def list_samples(
    project_id: int | None = None,
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> list[sample_schema.Sample]:
    """List samples, optionally filtered by project."""
    # If project_id is provided, verify user has access to that project
    if project_id is not None:
        project = await project_service.get_project(db, project_id)
        if project is None:
            raise HTTPException(status_code=404, detail="Project not found")
        if project.owner_id != current_user.id:
            raise HTTPException(
                status_code=403, detail="Not authorized to access this project"
            )

    samples = await sample_service.get_samples(
        db, skip=skip, limit=limit, project_id=project_id
    )
    return samples


@router.post(
    "/", response_model=sample_schema.Sample, status_code=status.HTTP_201_CREATED
)
async def create_sample(
    sample: sample_schema.SampleCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> sample_schema.Sample:
    """Create a new sample within a project."""
    # Verify user owns the project
    project = await project_service.get_project(db, sample.project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to add samples to this project"
        )

    return await sample_service.create_sample(db, sample)


@router.get("/{sample_id}", response_model=sample_schema.Sample)
async def get_sample(
    sample_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> sample_schema.Sample:
    """Get a specific sample by ID."""
    sample = await sample_service.get_sample(db, sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    # Verify user owns the project
    project = await project_service.get_project(db, sample.project_id)
    if project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to access this sample"
        )

    return sample


@router.put("/{sample_id}", response_model=sample_schema.Sample)
async def update_sample(
    sample_id: int,
    sample_update: sample_schema.SampleUpdate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> sample_schema.Sample:
    """Update a sample."""
    sample = await sample_service.get_sample(db, sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    # Verify user owns the project
    project = await project_service.get_project(db, sample.project_id)
    if project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to modify this sample"
        )

    updated_sample = await sample_service.update_sample(db, sample_id, sample_update)
    return updated_sample


@router.delete("/{sample_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_sample(
    sample_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Delete a sample."""
    sample = await sample_service.get_sample(db, sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    # Verify user owns the project
    project = await project_service.get_project(db, sample.project_id)
    if project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to delete this sample"
        )

    await sample_service.delete_sample(db, sample_id)
    return None

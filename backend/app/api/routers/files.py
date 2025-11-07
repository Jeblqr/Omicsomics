import io
import mimetypes

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.schemas import files as file_schema
from app.services import files as file_service
from app.services import samples as sample_service

router = APIRouter()


@router.post(
    "/upload", response_model=file_schema.File, status_code=status.HTTP_201_CREATED
)
async def upload_file(
    file: UploadFile = File(...),
    sample_id: int = Form(...),
    file_type: str = Form(...),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> file_schema.File:
    """
    Upload a file to S3 and create a database record.

    Args:
        file: The uploaded file
        sample_id: ID of the sample this file belongs to
        file_type: Type of file (e.g., 'fastq', 'bam', 'vcf', 'mzml')
        current_user: Authenticated user
        db: Database session
    """
    # Verify sample exists and user has access
    sample = await sample_service.get_sample(db, sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    # Verify user owns the project (through sample)
    from app.services import projects as project_service

    project = await project_service.get_project(db, sample.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to upload files to this sample"
        )

    # Determine MIME type
    mime_type = (
        file.content_type
        or mimetypes.guess_type(file.filename)[0]
        or "application/octet-stream"
    )

    # Read file content
    content = await file.read()
    file_data = io.BytesIO(content)

    # Create file record and upload to S3
    db_file = await file_service.create_file(
        db=db,
        file_data=file_data,
        filename=file.filename,
        file_type=file_type,
        mime_type=mime_type,
        sample_id=sample_id,
        uploaded_by_id=current_user.id,
    )

    return db_file


@router.get("/", response_model=list[file_schema.File])
async def list_files(
    sample_id: int | None = None,
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> list[file_schema.File]:
    """List files, optionally filtered by sample."""
    # If sample_id provided, verify access
    if sample_id is not None:
        sample = await sample_service.get_sample(db, sample_id)
        if sample is None:
            raise HTTPException(status_code=404, detail="Sample not found")

        from app.services import projects as project_service

        project = await project_service.get_project(db, sample.project_id)
        if project is None or project.owner_id != current_user.id:
            raise HTTPException(
                status_code=403, detail="Not authorized to access this sample"
            )

    files = await file_service.get_files(
        db, skip=skip, limit=limit, sample_id=sample_id
    )
    return files


@router.get("/{file_id}", response_model=file_schema.File)
async def get_file(
    file_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> file_schema.File:
    """Get file metadata by ID."""
    file = await file_service.get_file(db, file_id)
    if file is None:
        raise HTTPException(status_code=404, detail="File not found")

    # Verify access through sample -> project
    sample = await sample_service.get_sample(db, file.sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    from app.services import projects as project_service

    project = await project_service.get_project(db, sample.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to access this file"
        )

    return file


@router.get("/{file_id}/download")
async def get_download_url(
    file_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
) -> dict[str, str]:
    """Get a presigned URL for downloading the file."""
    file = await file_service.get_file(db, file_id)
    if file is None:
        raise HTTPException(status_code=404, detail="File not found")

    # Verify access
    sample = await sample_service.get_sample(db, file.sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    from app.services import projects as project_service

    project = await project_service.get_project(db, sample.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to download this file"
        )

    url = await file_service.get_download_url(db, file_id)
    if url is None:
        raise HTTPException(status_code=500, detail="Failed to generate download URL")

    return {"download_url": url}


@router.delete("/{file_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_file(
    file_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Delete a file."""
    file = await file_service.get_file(db, file_id)
    if file is None:
        raise HTTPException(status_code=404, detail="File not found")

    # Verify access
    sample = await sample_service.get_sample(db, file.sample_id)
    if sample is None:
        raise HTTPException(status_code=404, detail="Sample not found")

    from app.services import projects as project_service

    project = await project_service.get_project(db, sample.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(
            status_code=403, detail="Not authorized to delete this file"
        )

    await file_service.delete_file(db, file_id)
    return None

import io
import mimetypes

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.services import datafiles as datafile_service
from app.services import projects as project_service

router = APIRouter()


@router.post("/upload", status_code=status.HTTP_201_CREATED)
async def upload_datafile(
    project_id: int = Form(...),
    file: UploadFile = File(...),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    # Verify project exists and user is owner
    project = await project_service.get_project(db, project_id)
    if project is None:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    filename = file.filename or "unnamed"
    mime_type = (
        file.content_type
        or mimetypes.guess_type(filename)[0]
        or "application/octet-stream"
    )
    content = await file.read()
    file_obj = io.BytesIO(content)

    df = await datafile_service.create_datafile(
        db=db,
        file_data=file_obj,
        filename=filename,
        project_id=project_id,
        user_id=current_user.id,
        content_type=mime_type,
    )
    return {
        "id": df.id,
        "filename": df.filename,
        "object_key": df.object_key,
        "metadata_": df.metadata_,
        "size": df.size,
        "checksum": df.checksum,
        "project_id": df.project_id,
        "run_id": df.run_id,
        "uploaded_by_id": df.uploaded_by_id,
        "created_at": df.created_at,
        "updated_at": df.updated_at,
    }


@router.get("/")
async def list_datafiles(
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

    rows = await datafile_service.list_datafiles(db, project_id)
    return [
        {
            "id": r.id,
            "filename": r.filename,
            "object_key": r.object_key,
            "metadata_": r.metadata_,
            "size": r.size,
            "checksum": r.checksum,
            "project_id": r.project_id,
            "run_id": r.run_id,
            "uploaded_by_id": r.uploaded_by_id,
            "created_at": r.created_at,
            "updated_at": r.updated_at,
        }
        for r in rows
    ]


@router.get("/{datafile_id}/download")
async def download_datafile(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Return presigned URL to encrypted object (client must decrypt or use download-decrypted endpoint)."""
    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    from app.storage.s3 import s3_client

    url = await s3_client.get_presigned_url(df.object_key)
    return {"download_url": url}


@router.get("/{datafile_id}/download-decrypted")
async def download_datafile_decrypted(
    datafile_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Server-side decrypt and stream the file content."""
    from fastapi.responses import StreamingResponse
    from app.services import storage_service

    df = await datafile_service.get_datafile(db, datafile_id)
    if df is None:
        raise HTTPException(status_code=404, detail="DataFile not found")

    project = await project_service.get_project(db, df.project_id)
    if project is None or project.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    # Download and decrypt
    plaintext = await storage_service.download_and_decrypt(
        current_user.id, df.object_key
    )

    # Stream back to client
    return StreamingResponse(
        io.BytesIO(plaintext),
        media_type="application/octet-stream",
        headers={"Content-Disposition": f'attachment; filename="{df.filename}"'},
    )

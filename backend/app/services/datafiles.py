import hashlib
import io
from typing import BinaryIO

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.datafile import DataFile
from app.services import storage_service


async def calculate_checksum(file_data: BinaryIO) -> str:
    file_data.seek(0)
    md = hashlib.sha256()
    for chunk in iter(lambda: file_data.read(8192), b""):
        md.update(chunk)
    file_data.seek(0)
    return md.hexdigest()


async def create_datafile(
    db: AsyncSession,
    file_data: BinaryIO,
    filename: str,
    project_id: int,
    user_id: int,
    content_type: str = "application/octet-stream",
) -> DataFile:
    """Encrypt, upload and create DataFile record."""
    # Determine size
    file_data.seek(0, 2)
    size = file_data.tell()
    file_data.seek(0)

    # Checksum
    checksum = await calculate_checksum(file_data)

    # Read bytes for upload (we will re-use file_data)
    content = file_data.read()

    object_key, metadata = await storage_service.upload_encrypted_bytes(
        user_id=user_id,
        project_id=project_id,
        filename=filename,
        content=content,
        content_type=content_type,
    )

    db_obj = DataFile(
        filename=filename,
        object_key=object_key,
        metadata_=metadata,
        size=size,
        checksum=checksum,
        project_id=project_id,
        uploaded_by_id=user_id,
    )
    db.add(db_obj)
    await db.commit()
    await db.refresh(db_obj)
    return db_obj


async def get_datafile(db: AsyncSession, datafile_id: int) -> DataFile | None:
    r = await db.execute(select(DataFile).where(DataFile.id == datafile_id))
    return r.scalars().first()


async def list_datafiles(db: AsyncSession, project_id: int | None = None):
    q = select(DataFile).order_by(DataFile.created_at.desc())
    if project_id is not None:
        q = q.where(DataFile.project_id == project_id)
    r = await db.execute(q)
    return list(r.scalars().all())

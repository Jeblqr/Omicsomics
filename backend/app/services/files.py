import hashlib
from typing import BinaryIO

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.file import File
from app.schemas import files as file_schema
from app.storage.s3 import s3_client


async def get_file(db: AsyncSession, file_id: int) -> File | None:
    """Get a file by ID."""
    result = await db.execute(select(File).where(File.id == file_id))
    return result.scalars().first()


async def get_files(
    db: AsyncSession, skip: int = 0, limit: int = 100, sample_id: int | None = None
) -> list[File]:
    """Get files, optionally filtered by sample."""
    query = select(File).order_by(File.created_at.desc())
    if sample_id is not None:
        query = query.where(File.sample_id == sample_id)
    query = query.offset(skip).limit(limit)
    result = await db.execute(query)
    return list(result.scalars().all())


def calculate_checksum(file_data: BinaryIO) -> str:
    """Calculate MD5 checksum of file."""
    md5_hash = hashlib.md5()
    # Read in chunks to handle large files
    for chunk in iter(lambda: file_data.read(8192), b""):
        md5_hash.update(chunk)
    file_data.seek(0)  # Reset file pointer
    return md5_hash.hexdigest()


async def create_file(
    db: AsyncSession,
    file_data: BinaryIO,
    filename: str,
    file_type: str,
    mime_type: str,
    sample_id: int,
    uploaded_by_id: int,
) -> File:
    """Create a new file record and upload to S3."""
    # Calculate file size
    file_data.seek(0, 2)  # Seek to end
    size = file_data.tell()
    file_data.seek(0)  # Reset to beginning

    # Calculate checksum
    checksum = calculate_checksum(file_data)

    # Generate S3 key (filepath)
    filepath = f"samples/{sample_id}/{checksum}_{filename}"

    # Upload to S3
    await s3_client.upload_file(file_data, filepath, mime_type)

    # Create database record
    db_file = File(
        filename=filename,
        filepath=filepath,
        file_type=file_type,
        mime_type=mime_type,
        size=size,
        checksum=checksum,
        sample_id=sample_id,
        uploaded_by_id=uploaded_by_id,
    )
    db.add(db_file)
    await db.commit()
    await db.refresh(db_file)
    return db_file


async def delete_file(db: AsyncSession, file_id: int) -> bool:
    """Delete a file from database and S3."""
    db_file = await get_file(db, file_id)
    if db_file is None:
        return False

    # Delete from S3
    try:
        await s3_client.delete_file(db_file.filepath)
    except Exception as e:
        print(f"Warning: Failed to delete file from S3: {e}")

    # Delete from database
    await db.delete(db_file)
    await db.commit()
    return True


async def get_download_url(
    db: AsyncSession, file_id: int, expiration: int = 3600
) -> str | None:
    """Generate a presigned URL for file download."""
    db_file = await get_file(db, file_id)
    if db_file is None:
        return None

    return await s3_client.get_presigned_url(db_file.filepath, expiration)

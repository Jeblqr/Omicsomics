from datetime import datetime

from pydantic import BaseModel, ConfigDict


class FileBase(BaseModel):
    filename: str
    file_type: str
    mime_type: str | None = "application/octet-stream"


class FileCreate(FileBase):
    filepath: str
    size: int
    checksum: str
    sample_id: int
    uploaded_by_id: int


class FileUpdate(BaseModel):
    filename: str | None = None


class File(FileBase):
    id: int
    filepath: str
    size: int
    checksum: str
    sample_id: int
    uploaded_by_id: int
    created_at: datetime
    updated_at: datetime

    model_config = ConfigDict(from_attributes=True)

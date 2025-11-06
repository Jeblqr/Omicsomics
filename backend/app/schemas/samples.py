from datetime import datetime
from typing import Any

from pydantic import BaseModel, ConfigDict


class SampleBase(BaseModel):
    name: str
    description: str | None = None
    metadata_: dict[str, Any] | None = None


class SampleCreate(SampleBase):
    project_id: int


class SampleUpdate(BaseModel):
    name: str | None = None
    description: str | None = None
    metadata_: dict[str, Any] | None = None


class Sample(SampleBase):
    id: int
    project_id: int
    created_at: datetime
    updated_at: datetime

    model_config = ConfigDict(from_attributes=True)

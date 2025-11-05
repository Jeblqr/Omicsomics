from pydantic import BaseModel


class ProjectBase(BaseModel):
  name: str
  description: str | None = None


class ProjectRead(ProjectBase):
  id: int

  class Config:
    orm_mode = True

from datetime import datetime

from sqlalchemy import DateTime, String, Text, func
from sqlalchemy.orm import Mapped, mapped_column

from .base import Base


class Project(Base):
  __tablename__ = "projects"

  id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
  name: Mapped[str] = mapped_column(String(255), unique=True, index=True)
  description: Mapped[str] = mapped_column(Text, default="")
  created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), server_default=func.now())
  updated_at: Mapped[datetime] = mapped_column(
    DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
  )

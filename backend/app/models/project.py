from sqlalchemy.orm import Mapped, mapped_column

from .base import Base


class Project(Base):
  __tablename__ = "projects"

  id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
  name: Mapped[str] = mapped_column(unique=True, index=True)
  description: Mapped[str] = mapped_column(default="")

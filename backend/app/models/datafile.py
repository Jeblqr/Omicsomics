from datetime import datetime

from sqlalchemy import BigInteger, DateTime, ForeignKey, JSON, String, Text, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base


class DataFile(Base):
    __tablename__ = "data_files"

    id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
    filename: Mapped[str] = mapped_column(String(512), index=True, nullable=False)
    object_key: Mapped[str] = mapped_column(String(1024), unique=True, nullable=False)
    # 'metadata' is a reserved attribute name on Declarative base, use attribute name
    # 'metadata_' while keeping the DB column name as 'metadata'.
    metadata_: Mapped[dict] = mapped_column(
        "metadata", JSON, nullable=False
    )  # JSON package (nonce, alg, version, size, etc)
    size: Mapped[int] = mapped_column(BigInteger, nullable=False)
    checksum: Mapped[str] = mapped_column(String(128), nullable=False)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id"), nullable=False)
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.id"), nullable=True)
    uploaded_by_id: Mapped[int] = mapped_column(ForeignKey("users.id"), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
    )

    project = relationship("Project")
    run = relationship("Run")
    uploaded_by = relationship("User")

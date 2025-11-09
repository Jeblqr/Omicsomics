"""Custom user-defined pipeline model."""

from datetime import datetime

from sqlalchemy import DateTime, ForeignKey, String, Text, JSON, Boolean, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base


class CustomPipeline(Base):
    """User-defined custom pipelines."""

    __tablename__ = "custom_pipelines"

    id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    description: Mapped[str] = mapped_column(Text, default="")

    # Pipeline definition as a visual graph
    # Contains: nodes (steps), edges (connections), parameters, validation_rules
    definition: Mapped[dict] = mapped_column(JSON, nullable=False)

    # Category for organization
    category: Mapped[str] = mapped_column(String(100), default="custom")

    # Whether this pipeline is public/shared with other users
    is_public: Mapped[bool] = mapped_column(Boolean, default=False)

    # Owner of the pipeline
    owner_id: Mapped[int] = mapped_column(ForeignKey("users.id"), nullable=False)

    # Optional: based on template pipeline
    template_id: Mapped[str] = mapped_column(String(100), nullable=True)

    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
    )

    owner = relationship("User")

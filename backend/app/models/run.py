from datetime import datetime

from sqlalchemy import DateTime, ForeignKey, String, Text, JSON, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base


class Run(Base):
    __tablename__ = "runs"

    id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    description: Mapped[str] = mapped_column(Text, default="")
    status: Mapped[str] = mapped_column(String(50), default="pending", index=True)

    # Pipeline selection
    pipeline_type: Mapped[str] = mapped_column(
        String(50), nullable=True
    )  # 'template', 'custom', 'merged'
    pipeline_template_id: Mapped[str] = mapped_column(String(100), nullable=True)
    custom_pipeline_id: Mapped[int] = mapped_column(nullable=True)

    # Pipeline configuration stored as JSON
    # Contains: pipeline_definition, execution_order, tool_versions
    pipeline_config: Mapped[dict] = mapped_column(JSON, nullable=True)

    # Parameters for pipeline execution
    parameters: Mapped[dict] = mapped_column(JSON, default={})

    # Input/output file associations
    input_files: Mapped[list] = mapped_column(JSON, default=[])  # List of DataFile IDs
    output_files: Mapped[list] = mapped_column(JSON, default=[])  # List of DataFile IDs
    input_mapping: Mapped[dict] = mapped_column(JSON, default={})  # Map inputs to files

    # Execution results
    logs: Mapped[str] = mapped_column(Text, nullable=True)
    error_message: Mapped[str] = mapped_column(Text, nullable=True)
    progress: Mapped[float] = mapped_column(default=0.0)  # 0-100

    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id"), nullable=False)
    owner_id: Mapped[int] = mapped_column(ForeignKey("users.id"), nullable=False)
    started_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=True)
    finished_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
    )

    project = relationship("Project")
    owner = relationship("User")

"""
Visualization Workspace Models

Database models for storing dashboard configurations:
- Dashboard: Top-level container for visualization workspace
- Panel: Individual visualization panel within a dashboard
"""

from datetime import datetime
from typing import Optional

from sqlalchemy import DateTime, ForeignKey, String, Text, JSON, Boolean, Integer, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base


class Dashboard(Base):
    """Dashboard workspace for multi-panel visualizations."""

    __tablename__ = "dashboards"

    id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    description: Mapped[str] = mapped_column(Text, default="")

    # Layout configuration as JSON
    # Contains: grid layout, panel positions, sizes, responsive breakpoints
    layout: Mapped[dict] = mapped_column(JSON, nullable=False, default=dict)

    # Dashboard metadata
    # Contains: tags, categories, custom properties
    metadata: Mapped[dict] = mapped_column(JSON, default=dict)

    # Sharing and visibility
    is_public: Mapped[bool] = mapped_column(Boolean, default=False)
    is_template: Mapped[bool] = mapped_column(Boolean, default=False)

    # Ownership
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id"), nullable=False)
    owner_id: Mapped[int] = mapped_column(ForeignKey("users.id"), nullable=False)

    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
    )

    # Relationships
    project = relationship("Project")
    owner = relationship("User")
    panels = relationship(
        "Panel", back_populates="dashboard", cascade="all, delete-orphan"
    )


class Panel(Base):
    """Individual visualization panel within a dashboard."""

    __tablename__ = "panels"

    id: Mapped[int] = mapped_column(primary_key=True, index=True, autoincrement=True)
    dashboard_id: Mapped[int] = mapped_column(
        ForeignKey("dashboards.id"), nullable=False
    )

    # Panel identification
    panel_key: Mapped[str] = mapped_column(
        String(100), nullable=False
    )  # Unique within dashboard
    title: Mapped[str] = mapped_column(String(255), nullable=False)
    description: Mapped[str] = mapped_column(Text, default="")

    # Visualization type
    viz_type: Mapped[str] = mapped_column(String(50), nullable=False)
    # e.g., 'line', 'bar', 'scatter', 'heatmap', 'table', 'metric', 'text'

    # Data source configuration
    # Contains: file_id, dataset_id, data_url, static_data
    data_source: Mapped[dict] = mapped_column(JSON, nullable=False)

    # Visualization configuration
    # Contains: chart options, axes, colors, filters, etc.
    viz_config: Mapped[dict] = mapped_column(JSON, nullable=False, default=dict)

    # Layout position (managed by dashboard layout)
    position: Mapped[dict] = mapped_column(JSON, nullable=False)
    # Contains: x, y, w, h (grid units)

    # Refresh and update settings
    auto_refresh: Mapped[bool] = mapped_column(Boolean, default=False)
    refresh_interval: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    # Interval in seconds

    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now()
    )

    # Relationships
    dashboard = relationship("Dashboard", back_populates="panels")

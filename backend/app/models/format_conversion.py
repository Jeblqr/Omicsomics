"""Database models for format conversions."""

from sqlalchemy import (
    Column,
    Integer,
    String,
    Float,
    DateTime,
    JSON,
    Boolean,
    ForeignKey,
)
from sqlalchemy.orm import relationship
from datetime import datetime
from app.database import Base


class FormatConversion(Base):
    """Track format conversion operations."""

    __tablename__ = "format_conversions"

    id = Column(Integer, primary_key=True, index=True)

    # Source information
    source_file_id = Column(Integer, ForeignKey("files.id"), nullable=True)
    source_path = Column(String, nullable=False)
    source_format = Column(String, nullable=False)
    source_size_bytes = Column(Integer)

    # Target information
    target_file_id = Column(Integer, ForeignKey("files.id"), nullable=True)
    target_path = Column(String, nullable=False)
    target_format = Column(String, nullable=False)
    target_size_bytes = Column(Integer)

    # Conversion details
    conversion_path = Column(JSON)  # List of intermediate formats
    conversion_mode = Column(String, nullable=False)  # 'auto' or 'manual'

    # Execution info
    status = Column(String, default="pending")  # pending, running, completed, failed
    duration_seconds = Column(Float)
    error_message = Column(String, nullable=True)

    # Metadata
    created_by = Column(Integer, ForeignKey("users.id"))
    created_at = Column(DateTime, default=datetime.utcnow)
    completed_at = Column(DateTime, nullable=True)

    # Additional parameters used for conversion
    parameters = Column(JSON, nullable=True)

    # Relationships
    creator = relationship("User", foreign_keys=[created_by])
    source_file = relationship("File", foreign_keys=[source_file_id])
    target_file = relationship("File", foreign_keys=[target_file_id])


class ConversionRule(Base):
    """Store conversion rules and scripts."""

    __tablename__ = "conversion_rules"

    id = Column(Integer, primary_key=True, index=True)

    # Format pair
    from_format = Column(String, nullable=False)
    to_format = Column(String, nullable=False)

    # Conversion method
    method = Column(
        String, nullable=False
    )  # 'pandas', 'r_script', 'python_script', 'binary'
    script_path = Column(String, nullable=True)  # Path to conversion script
    command_template = Column(
        String, nullable=True
    )  # Command template for binary tools

    # Performance metrics
    avg_time_per_gb = Column(Float)  # Average conversion time per GB
    success_rate = Column(Float, default=1.0)  # Success rate (0-1)
    usage_count = Column(Integer, default=0)  # Number of times used

    # Requirements
    requires_runtime = Column(String, nullable=True)  # 'r', 'python', 'binary', None
    requires_packages = Column(JSON, nullable=True)  # List of required packages

    # Status
    is_active = Column(Boolean, default=True)

    # Metadata
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    # Notes
    description = Column(String, nullable=True)
    notes = Column(String, nullable=True)

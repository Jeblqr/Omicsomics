"""
Database models for custom script tools.

Allows users to upload and integrate custom Python/R scripts with
parameter configuration and execution tracking.
"""

from datetime import datetime
from typing import Optional
from sqlalchemy import (
    Column,
    Integer,
    String,
    DateTime,
    Boolean,
    Text,
    JSON,
    ForeignKey,
    Enum as SQLEnum,
)
from sqlalchemy.orm import relationship
import enum

from ..database import Base


class ScriptLanguage(str, enum.Enum):
    """Supported script languages"""

    PYTHON = "python"
    R = "r"
    BASH = "bash"


class ScriptStatus(str, enum.Enum):
    """Script execution status"""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ScriptVisibility(str, enum.Enum):
    """Script visibility"""

    PRIVATE = "private"
    PROJECT = "project"
    PUBLIC = "public"


class CustomScript(Base):
    """
    Custom script definition table.

    Stores user-uploaded scripts with metadata and configuration.
    """

    __tablename__ = "custom_scripts"

    id = Column(Integer, primary_key=True, index=True)
    script_key = Column(String(100), unique=True, index=True, nullable=False)

    # Script metadata
    name = Column(String(255), nullable=False, index=True)
    description = Column(Text, nullable=True)
    version = Column(String(50), default="1.0.0")

    # User and project
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False, index=True)
    project_id = Column(Integer, ForeignKey("projects.id"), nullable=True, index=True)

    # Script details
    language = Column(SQLEnum(ScriptLanguage), nullable=False)
    script_content = Column(Text, nullable=False)
    entry_point = Column(String(255), nullable=True)  # Function/class to call

    # Parameters definition
    parameters_schema = Column(JSON, nullable=True)  # JSON Schema for parameters
    default_parameters = Column(JSON, nullable=True)

    # Dependencies
    requirements = Column(JSON, nullable=True)  # List of package requirements

    # Execution configuration
    timeout_seconds = Column(Integer, default=300)  # 5 minutes default
    max_memory_mb = Column(Integer, default=1024)  # 1GB default

    # Visibility and sharing
    visibility = Column(SQLEnum(ScriptVisibility), default=ScriptVisibility.PRIVATE)
    is_verified = Column(Boolean, default=False)  # Admin-verified scripts
    is_template = Column(Boolean, default=False)

    # Storage
    storage_path = Column(String(500), nullable=True)  # Path in storage

    # Usage statistics
    execution_count = Column(Integer, default=0)
    success_count = Column(Integer, default=0)
    failure_count = Column(Integer, default=0)
    average_duration_seconds = Column(Integer, nullable=True)

    # Tags and categories
    tags = Column(JSON, nullable=True)
    category = Column(String(100), nullable=True, index=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    last_executed_at = Column(DateTime, nullable=True)

    # Relationships
    user = relationship("User", back_populates="custom_scripts")
    project = relationship("Project", back_populates="custom_scripts")
    executions = relationship(
        "ScriptExecution", back_populates="script", cascade="all, delete-orphan"
    )

    def __repr__(self):
        return (
            f"<CustomScript(id={self.id}, name={self.name}, language={self.language})>"
        )

    def to_dict(self):
        """Convert to dictionary representation"""
        return {
            "id": self.id,
            "script_key": self.script_key,
            "name": self.name,
            "description": self.description,
            "version": self.version,
            "user_id": self.user_id,
            "project_id": self.project_id,
            "language": self.language.value if self.language else None,
            "entry_point": self.entry_point,
            "parameters_schema": self.parameters_schema,
            "default_parameters": self.default_parameters,
            "requirements": self.requirements,
            "timeout_seconds": self.timeout_seconds,
            "max_memory_mb": self.max_memory_mb,
            "visibility": self.visibility.value if self.visibility else None,
            "is_verified": self.is_verified,
            "is_template": self.is_template,
            "execution_count": self.execution_count,
            "success_count": self.success_count,
            "failure_count": self.failure_count,
            "average_duration_seconds": self.average_duration_seconds,
            "tags": self.tags,
            "category": self.category,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
            "last_executed_at": (
                self.last_executed_at.isoformat() if self.last_executed_at else None
            ),
        }


class ScriptExecution(Base):
    """
    Script execution tracking table.

    Records each execution with parameters, results, and logs.
    """

    __tablename__ = "script_executions"

    id = Column(Integer, primary_key=True, index=True)
    execution_key = Column(String(100), unique=True, index=True, nullable=False)

    # Script reference
    script_id = Column(
        Integer, ForeignKey("custom_scripts.id"), nullable=False, index=True
    )

    # User and project
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False, index=True)
    project_id = Column(Integer, ForeignKey("projects.id"), nullable=True, index=True)

    # Execution parameters
    parameters = Column(JSON, nullable=True)

    # Input/output
    input_files = Column(JSON, nullable=True)  # List of input file IDs
    output_files = Column(JSON, nullable=True)  # List of output file IDs

    # Execution status
    status = Column(SQLEnum(ScriptStatus), nullable=False, default=ScriptStatus.PENDING)

    # Results
    result_data = Column(JSON, nullable=True)  # Structured results
    output_text = Column(Text, nullable=True)  # stdout
    error_text = Column(Text, nullable=True)  # stderr
    exit_code = Column(Integer, nullable=True)

    # Performance metrics
    duration_seconds = Column(Integer, nullable=True)
    memory_usage_mb = Column(Integer, nullable=True)
    cpu_usage_percent = Column(Integer, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

    # Relationships
    script = relationship("CustomScript", back_populates="executions")
    user = relationship("User")
    project = relationship("Project")

    def __repr__(self):
        return f"<ScriptExecution(id={self.id}, script_id={self.script_id}, status={self.status})>"

    def to_dict(self):
        """Convert to dictionary representation"""
        return {
            "id": self.id,
            "execution_key": self.execution_key,
            "script_id": self.script_id,
            "user_id": self.user_id,
            "project_id": self.project_id,
            "parameters": self.parameters,
            "input_files": self.input_files,
            "output_files": self.output_files,
            "status": self.status.value if self.status else None,
            "result_data": self.result_data,
            "output_text": self.output_text,
            "error_text": self.error_text,
            "exit_code": self.exit_code,
            "duration_seconds": self.duration_seconds,
            "memory_usage_mb": self.memory_usage_mb,
            "cpu_usage_percent": self.cpu_usage_percent,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": (
                self.completed_at.isoformat() if self.completed_at else None
            ),
        }

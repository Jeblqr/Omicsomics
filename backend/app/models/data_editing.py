"""
Database models for data editing functionality.

Provides in-place data operations with preview, validation, undo/redo,
and operation history tracking.
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


class EditOperationType(str, enum.Enum):
    """Types of edit operations"""

    # Row operations
    DELETE_ROWS = "delete_rows"
    FILTER_ROWS = "filter_rows"
    SORT_ROWS = "sort_rows"
    DEDUPLICATE = "deduplicate"

    # Column operations
    ADD_COLUMN = "add_column"
    DELETE_COLUMN = "delete_column"
    RENAME_COLUMN = "rename_column"
    REORDER_COLUMNS = "reorder_columns"

    # Value operations
    REPLACE_VALUES = "replace_values"
    FILL_MISSING = "fill_missing"
    TRANSFORM_VALUES = "transform_values"
    NORMALIZE = "normalize"

    # Type operations
    CONVERT_TYPE = "convert_type"
    CAST_COLUMN = "cast_column"

    # Merge operations
    JOIN_DATA = "join_data"
    CONCAT_DATA = "concat_data"

    # Custom
    CUSTOM_SCRIPT = "custom_script"


class EditStatus(str, enum.Enum):
    """Edit session status"""

    DRAFT = "draft"
    PREVIEWING = "previewing"
    APPLYING = "applying"
    APPLIED = "applied"
    FAILED = "failed"
    REVERTED = "reverted"


class EditSession(Base):
    """
    Edit session tracking table.

    Manages in-place data editing with preview, validation, and undo capabilities.
    """

    __tablename__ = "edit_sessions"

    id = Column(Integer, primary_key=True, index=True)
    session_key = Column(String(100), unique=True, index=True, nullable=False)

    # Session metadata
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)

    # User and project
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False, index=True)
    project_id = Column(Integer, ForeignKey("projects.id"), nullable=True, index=True)

    # Target data
    file_id = Column(Integer, ForeignKey("files.id"), nullable=False, index=True)
    dataset_id = Column(Integer, ForeignKey("datasets.id"), nullable=True, index=True)

    # Session state
    status = Column(SQLEnum(EditStatus), nullable=False, default=EditStatus.DRAFT)

    # Operations
    operations = Column(JSON, nullable=False, default=list)  # List of operations
    current_operation_index = Column(Integer, default=-1)  # For undo/redo

    # Preview data
    preview_data = Column(JSON, nullable=True)  # Sample of result
    preview_summary = Column(JSON, nullable=True)  # Statistics

    # Validation
    validation_errors = Column(JSON, nullable=True)  # List of errors
    validation_warnings = Column(JSON, nullable=True)  # List of warnings
    is_valid = Column(Boolean, default=False)

    # Backup
    backup_path = Column(String(500), nullable=True)  # Original file backup

    # Result
    output_file_id = Column(Integer, ForeignKey("files.id"), nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    applied_at = Column(DateTime, nullable=True)

    # Relationships
    user = relationship("User", foreign_keys=[user_id])
    project = relationship("Project", back_populates="edit_sessions")
    file = relationship("File", foreign_keys=[file_id])
    output_file = relationship("File", foreign_keys=[output_file_id])
    dataset = relationship("Dataset", back_populates="edit_sessions")

    def __repr__(self):
        return f"<EditSession(id={self.id}, session_key={self.session_key}, status={self.status})>"

    def to_dict(self):
        """Convert to dictionary representation"""
        return {
            "id": self.id,
            "session_key": self.session_key,
            "name": self.name,
            "description": self.description,
            "user_id": self.user_id,
            "project_id": self.project_id,
            "file_id": self.file_id,
            "dataset_id": self.dataset_id,
            "status": self.status.value if self.status else None,
            "operations": self.operations,
            "current_operation_index": self.current_operation_index,
            "preview_data": self.preview_data,
            "preview_summary": self.preview_summary,
            "validation_errors": self.validation_errors,
            "validation_warnings": self.validation_warnings,
            "is_valid": self.is_valid,
            "backup_path": self.backup_path,
            "output_file_id": self.output_file_id,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
            "applied_at": self.applied_at.isoformat() if self.applied_at else None,
        }


class EditOperation(Base):
    """
    Individual edit operation record.

    Stores detailed information about each operation for audit trail.
    """

    __tablename__ = "edit_operations"

    id = Column(Integer, primary_key=True, index=True)

    # Session reference
    session_id = Column(
        Integer, ForeignKey("edit_sessions.id"), nullable=False, index=True
    )

    # Operation details
    operation_type = Column(SQLEnum(EditOperationType), nullable=False)
    operation_index = Column(Integer, nullable=False)  # Order in session

    # Operation parameters
    parameters = Column(JSON, nullable=False)

    # Operation metadata
    description = Column(Text, nullable=True)

    # Impact
    rows_affected = Column(Integer, nullable=True)
    columns_affected = Column(Integer, nullable=True)

    # Undo information
    undo_data = Column(JSON, nullable=True)  # Data needed to revert
    can_undo = Column(Boolean, default=True)

    # Execution
    executed = Column(Boolean, default=False)
    execution_time_ms = Column(Integer, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    executed_at = Column(DateTime, nullable=True)

    # Relationships
    session = relationship("EditSession", backref="operation_records")

    def __repr__(self):
        return f"<EditOperation(id={self.id}, type={self.operation_type}, executed={self.executed})>"

    def to_dict(self):
        """Convert to dictionary representation"""
        return {
            "id": self.id,
            "session_id": self.session_id,
            "operation_type": (
                self.operation_type.value if self.operation_type else None
            ),
            "operation_index": self.operation_index,
            "parameters": self.parameters,
            "description": self.description,
            "rows_affected": self.rows_affected,
            "columns_affected": self.columns_affected,
            "undo_data": self.undo_data,
            "can_undo": self.can_undo,
            "executed": self.executed,
            "execution_time_ms": self.execution_time_ms,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "executed_at": self.executed_at.isoformat() if self.executed_at else None,
        }

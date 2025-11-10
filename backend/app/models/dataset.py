"""
Dataset Models

Database models for dataset management with version control,
file lineage tracking, and metadata management.
"""

from sqlalchemy import (
    Column,
    String,
    Integer,
    DateTime,
    Text,
    JSON,
    ForeignKey,
    Boolean,
    Table,
)
from sqlalchemy.orm import relationship
from datetime import datetime
import uuid

from app.database import Base


# Association table for dataset files
dataset_files = Table(
    "dataset_files",
    Base.metadata,
    Column("dataset_id", String, ForeignKey("datasets.id", ondelete="CASCADE")),
    Column(
        "file_id", String, ForeignKey("dataset_file_entries.id", ondelete="CASCADE")
    ),
)

# Association table for dataset tags
dataset_tags = Table(
    "dataset_tags",
    Base.metadata,
    Column("dataset_id", String, ForeignKey("datasets.id", ondelete="CASCADE")),
    Column("tag_id", String, ForeignKey("tags.id", ondelete="CASCADE")),
)


class Dataset(Base):
    """Dataset model with version control"""

    __tablename__ = "datasets"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    project_id = Column(
        String, ForeignKey("projects.id", ondelete="CASCADE"), nullable=False
    )

    # Basic info
    name = Column(String, nullable=False, index=True)
    description = Column(Text)
    version = Column(Integer, default=1, nullable=False)

    # Metadata
    data_type = Column(String)  # e.g., "genomics", "proteomics", "transcriptomics"
    file_format = Column(String)  # Primary file format
    metadata = Column(JSON, default=dict)  # Flexible metadata storage

    # Status
    status = Column(String, default="active")  # active, archived, deprecated
    is_public = Column(Boolean, default=False)

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    # User info
    created_by = Column(String, ForeignKey("users.id"))

    # Version control
    parent_id = Column(String, ForeignKey("datasets.id"))  # For versioning

    # Relationships
    project = relationship("Project", back_populates="datasets")
    files = relationship(
        "DatasetFileEntry",
        secondary=dataset_files,
        back_populates="datasets",
        lazy="dynamic",
    )
    tags = relationship(
        "Tag", secondary=dataset_tags, back_populates="datasets", lazy="dynamic"
    )
    versions = relationship(
        "Dataset", backref="parent", remote_side=[id], lazy="dynamic"
    )
    lineage_records = relationship(
        "DatasetLineage",
        foreign_keys="DatasetLineage.dataset_id",
        back_populates="dataset",
        lazy="dynamic",
    )

    def __repr__(self):
        return f"<Dataset {self.name} v{self.version}>"


class DatasetFileEntry(Base):
    """File entry in a dataset"""

    __tablename__ = "dataset_file_entries"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))

    # File info
    file_path = Column(String, nullable=False)  # Absolute or relative path
    file_name = Column(String, nullable=False)
    file_type = Column(String)  # File extension
    file_size = Column(Integer)  # Size in bytes

    # Role in dataset
    role = Column(String)  # e.g., "primary", "metadata", "supplementary"
    description = Column(Text)

    # Checksums for integrity
    md5_hash = Column(String)
    sha256_hash = Column(String)

    # Timestamps
    added_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    datasets = relationship("Dataset", secondary=dataset_files, back_populates="files")
    lineage_as_source = relationship(
        "DatasetLineage",
        foreign_keys="DatasetLineage.source_file_id",
        back_populates="source_file",
        lazy="dynamic",
    )
    lineage_as_output = relationship(
        "DatasetLineage",
        foreign_keys="DatasetLineage.output_file_id",
        back_populates="output_file",
        lazy="dynamic",
    )

    def __repr__(self):
        return f"<DatasetFileEntry {self.file_name}>"


class DatasetLineage(Base):
    """Track data lineage and provenance"""

    __tablename__ = "dataset_lineage"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))

    # Dataset
    dataset_id = Column(String, ForeignKey("datasets.id", ondelete="CASCADE"))

    # Lineage tracking
    operation_type = Column(String)  # e.g., "merge", "conversion", "pipeline", "manual"
    operation_id = Column(String)  # ID of the operation (run_id, job_id, etc.)

    # Source and output files
    source_file_id = Column(String, ForeignKey("dataset_file_entries.id"))
    output_file_id = Column(String, ForeignKey("dataset_file_entries.id"))

    # Operation details
    operation_params = Column(JSON)  # Parameters used in the operation
    tool_name = Column(String)  # Tool or pipeline used
    tool_version = Column(String)

    # Timestamps
    executed_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    dataset = relationship("Dataset", back_populates="lineage_records")
    source_file = relationship(
        "DatasetFileEntry",
        foreign_keys=[source_file_id],
        back_populates="lineage_as_source",
    )
    output_file = relationship(
        "DatasetFileEntry",
        foreign_keys=[output_file_id],
        back_populates="lineage_as_output",
    )

    def __repr__(self):
        return f"<DatasetLineage {self.operation_type}>"


class Tag(Base):
    """Tags for organizing datasets"""

    __tablename__ = "tags"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String, unique=True, nullable=False, index=True)
    color = Column(String, default="#3b82f6")  # Hex color code
    description = Column(Text)

    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    datasets = relationship("Dataset", secondary=dataset_tags, back_populates="tags")

    def __repr__(self):
        return f"<Tag {self.name}>"

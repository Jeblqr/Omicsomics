"""
Dataset Manager Service

Manages datasets with version control, file lineage tracking,
and metadata management.
"""

from typing import List, Optional, Dict, Any
from pathlib import Path
from datetime import datetime
import hashlib
import logging

from sqlalchemy.orm import Session
from sqlalchemy import or_, and_

from app.models.dataset import (
    Dataset,
    DatasetFileEntry,
    DatasetLineage,
    Tag,
    dataset_files,
    dataset_tags,
)

logger = logging.getLogger(__name__)


class DatasetManager:
    """Service for managing datasets"""

    def __init__(self, db: Session):
        self.db = db

    # Dataset operations
    def create_dataset(
        self,
        project_id: str,
        name: str,
        description: Optional[str] = None,
        data_type: Optional[str] = None,
        file_format: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        created_by: Optional[str] = None,
        tags: Optional[List[str]] = None,
    ) -> Dataset:
        """Create a new dataset"""

        dataset = Dataset(
            project_id=project_id,
            name=name,
            description=description,
            data_type=data_type,
            file_format=file_format,
            metadata=metadata or {},
            created_by=created_by,
        )

        # Add tags
        if tags:
            for tag_name in tags:
                tag = self.get_or_create_tag(tag_name)
                dataset.tags.append(tag)

        self.db.add(dataset)
        self.db.commit()
        self.db.refresh(dataset)

        logger.info(f"Created dataset: {dataset.name} (ID: {dataset.id})")
        return dataset

    def get_dataset(self, dataset_id: str) -> Optional[Dataset]:
        """Get dataset by ID"""
        return self.db.query(Dataset).filter(Dataset.id == dataset_id).first()

    def list_datasets(
        self,
        project_id: Optional[str] = None,
        data_type: Optional[str] = None,
        status: Optional[str] = None,
        tags: Optional[List[str]] = None,
        search_query: Optional[str] = None,
        skip: int = 0,
        limit: int = 100,
    ) -> List[Dataset]:
        """List datasets with filters"""

        query = self.db.query(Dataset)

        if project_id:
            query = query.filter(Dataset.project_id == project_id)

        if data_type:
            query = query.filter(Dataset.data_type == data_type)

        if status:
            query = query.filter(Dataset.status == status)

        if tags:
            # Filter by tags
            for tag_name in tags:
                query = query.filter(Dataset.tags.any(Tag.name == tag_name))

        if search_query:
            # Search in name and description
            search_pattern = f"%{search_query}%"
            query = query.filter(
                or_(
                    Dataset.name.ilike(search_pattern),
                    Dataset.description.ilike(search_pattern),
                )
            )

        query = query.order_by(Dataset.created_at.desc())
        query = query.offset(skip).limit(limit)

        return query.all()

    def update_dataset(
        self,
        dataset_id: str,
        name: Optional[str] = None,
        description: Optional[str] = None,
        data_type: Optional[str] = None,
        file_format: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        status: Optional[str] = None,
    ) -> Optional[Dataset]:
        """Update dataset metadata"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return None

        if name is not None:
            dataset.name = name
        if description is not None:
            dataset.description = description
        if data_type is not None:
            dataset.data_type = data_type
        if file_format is not None:
            dataset.file_format = file_format
        if metadata is not None:
            dataset.metadata = metadata
        if status is not None:
            dataset.status = status

        dataset.updated_at = datetime.utcnow()

        self.db.commit()
        self.db.refresh(dataset)

        logger.info(f"Updated dataset: {dataset.name}")
        return dataset

    def delete_dataset(self, dataset_id: str) -> bool:
        """Delete dataset (soft delete by setting status to 'deleted')"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return False

        dataset.status = "deleted"
        dataset.updated_at = datetime.utcnow()

        self.db.commit()
        logger.info(f"Deleted dataset: {dataset.name}")
        return True

    def archive_dataset(self, dataset_id: str) -> bool:
        """Archive dataset"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return False

        dataset.status = "archived"
        dataset.updated_at = datetime.utcnow()

        self.db.commit()
        logger.info(f"Archived dataset: {dataset.name}")
        return True

    # Version control
    def create_version(
        self,
        dataset_id: str,
        description: Optional[str] = None,
        created_by: Optional[str] = None,
    ) -> Optional[Dataset]:
        """Create a new version of a dataset"""

        parent = self.get_dataset(dataset_id)
        if not parent:
            return None

        # Create new version
        new_version = Dataset(
            project_id=parent.project_id,
            name=parent.name,
            description=description or f"Version {parent.version + 1} of {parent.name}",
            version=parent.version + 1,
            data_type=parent.data_type,
            file_format=parent.file_format,
            metadata=parent.metadata.copy(),
            parent_id=parent.id,
            created_by=created_by,
        )

        # Copy files from parent
        for file in parent.files:
            new_version.files.append(file)

        # Copy tags from parent
        for tag in parent.tags:
            new_version.tags.append(tag)

        self.db.add(new_version)
        self.db.commit()
        self.db.refresh(new_version)

        logger.info(f"Created version {new_version.version} of dataset: {parent.name}")
        return new_version

    def get_versions(self, dataset_id: str) -> List[Dataset]:
        """Get all versions of a dataset"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return []

        # If this is a version, get the parent
        root_id = dataset.parent_id if dataset.parent_id else dataset.id

        # Get all versions
        versions = (
            self.db.query(Dataset)
            .filter(or_(Dataset.id == root_id, Dataset.parent_id == root_id))
            .order_by(Dataset.version)
            .all()
        )

        return versions

    # File operations
    def add_file(
        self,
        dataset_id: str,
        file_path: str,
        role: str = "primary",
        description: Optional[str] = None,
        compute_hash: bool = True,
    ) -> Optional[DatasetFileEntry]:
        """Add a file to dataset"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return None

        path = Path(file_path)

        # Check if file exists
        if not path.exists():
            logger.error(f"File not found: {file_path}")
            return None

        # Compute hashes if requested
        md5_hash = None
        sha256_hash = None
        if compute_hash:
            try:
                md5_hash = self._compute_md5(path)
                sha256_hash = self._compute_sha256(path)
            except Exception as e:
                logger.warning(f"Failed to compute hash for {file_path}: {e}")

        # Create file entry
        file_entry = DatasetFileEntry(
            file_path=str(path.absolute()),
            file_name=path.name,
            file_type=path.suffix,
            file_size=path.stat().st_size,
            role=role,
            description=description,
            md5_hash=md5_hash,
            sha256_hash=sha256_hash,
        )

        self.db.add(file_entry)
        dataset.files.append(file_entry)
        dataset.updated_at = datetime.utcnow()

        self.db.commit()
        self.db.refresh(file_entry)

        logger.info(f"Added file {path.name} to dataset {dataset.name}")
        return file_entry

    def remove_file(self, dataset_id: str, file_id: str) -> bool:
        """Remove a file from dataset"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return False

        file_entry = (
            self.db.query(DatasetFileEntry)
            .filter(DatasetFileEntry.id == file_id)
            .first()
        )

        if not file_entry:
            return False

        dataset.files.remove(file_entry)
        dataset.updated_at = datetime.utcnow()

        self.db.commit()
        logger.info(f"Removed file {file_entry.file_name} from dataset {dataset.name}")
        return True

    def list_files(self, dataset_id: str) -> List[DatasetFileEntry]:
        """List all files in dataset"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return []

        return list(dataset.files.all())

    # Lineage tracking
    def add_lineage(
        self,
        dataset_id: str,
        operation_type: str,
        source_file_id: Optional[str] = None,
        output_file_id: Optional[str] = None,
        operation_id: Optional[str] = None,
        operation_params: Optional[Dict[str, Any]] = None,
        tool_name: Optional[str] = None,
        tool_version: Optional[str] = None,
    ) -> Optional[DatasetLineage]:
        """Add lineage record"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return None

        lineage = DatasetLineage(
            dataset_id=dataset_id,
            operation_type=operation_type,
            operation_id=operation_id,
            source_file_id=source_file_id,
            output_file_id=output_file_id,
            operation_params=operation_params or {},
            tool_name=tool_name,
            tool_version=tool_version,
        )

        self.db.add(lineage)
        self.db.commit()
        self.db.refresh(lineage)

        logger.info(
            f"Added lineage record: {operation_type} for dataset {dataset.name}"
        )
        return lineage

    def get_lineage(self, dataset_id: str) -> List[DatasetLineage]:
        """Get lineage history for dataset"""

        return (
            self.db.query(DatasetLineage)
            .filter(DatasetLineage.dataset_id == dataset_id)
            .order_by(DatasetLineage.executed_at.desc())
            .all()
        )

    def trace_file_lineage(self, file_id: str) -> Dict[str, Any]:
        """Trace the complete lineage of a file"""

        file_entry = (
            self.db.query(DatasetFileEntry)
            .filter(DatasetFileEntry.id == file_id)
            .first()
        )

        if not file_entry:
            return {}

        # Get operations where this file is output
        as_output = list(file_entry.lineage_as_output.all())

        # Get operations where this file is source
        as_source = list(file_entry.lineage_as_source.all())

        # Build lineage tree
        lineage_tree = {
            "file": {
                "id": file_entry.id,
                "name": file_entry.file_name,
                "path": file_entry.file_path,
            },
            "created_by": [
                {
                    "operation": op.operation_type,
                    "tool": op.tool_name,
                    "executed_at": op.executed_at.isoformat(),
                    "source_files": (
                        [{"id": op.source_file.id, "name": op.source_file.file_name}]
                        if op.source_file
                        else []
                    ),
                }
                for op in as_output
            ],
            "used_in": [
                {
                    "operation": op.operation_type,
                    "tool": op.tool_name,
                    "executed_at": op.executed_at.isoformat(),
                    "output_files": (
                        [{"id": op.output_file.id, "name": op.output_file.file_name}]
                        if op.output_file
                        else []
                    ),
                }
                for op in as_source
            ],
        }

        return lineage_tree

    # Tag operations
    def get_or_create_tag(self, tag_name: str, color: Optional[str] = None) -> Tag:
        """Get existing tag or create new one"""

        tag = self.db.query(Tag).filter(Tag.name == tag_name).first()

        if not tag:
            tag = Tag(name=tag_name, color=color or "#3b82f6")
            self.db.add(tag)
            self.db.commit()
            self.db.refresh(tag)
            logger.info(f"Created tag: {tag_name}")

        return tag

    def add_tag_to_dataset(self, dataset_id: str, tag_name: str) -> bool:
        """Add tag to dataset"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return False

        tag = self.get_or_create_tag(tag_name)

        if tag not in dataset.tags:
            dataset.tags.append(tag)
            self.db.commit()
            logger.info(f"Added tag {tag_name} to dataset {dataset.name}")

        return True

    def remove_tag_from_dataset(self, dataset_id: str, tag_id: str) -> bool:
        """Remove tag from dataset"""

        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return False

        tag = self.db.query(Tag).filter(Tag.id == tag_id).first()
        if not tag:
            return False

        if tag in dataset.tags:
            dataset.tags.remove(tag)
            self.db.commit()
            logger.info(f"Removed tag {tag.name} from dataset {dataset.name}")

        return True

    def list_tags(self) -> List[Tag]:
        """List all tags"""
        return self.db.query(Tag).order_by(Tag.name).all()

    # Utility methods
    @staticmethod
    def _compute_md5(file_path: Path) -> str:
        """Compute MD5 hash of file"""
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    @staticmethod
    def _compute_sha256(file_path: Path) -> str:
        """Compute SHA256 hash of file"""
        hash_sha256 = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_sha256.update(chunk)
        return hash_sha256.hexdigest()

    # Statistics
    def get_statistics(self, project_id: Optional[str] = None) -> Dict[str, Any]:
        """Get dataset statistics"""

        query = self.db.query(Dataset)
        if project_id:
            query = query.filter(Dataset.project_id == project_id)

        total = query.count()
        active = query.filter(Dataset.status == "active").count()
        archived = query.filter(Dataset.status == "archived").count()

        # Count by data type
        data_types = {}
        for dataset in query.all():
            dt = dataset.data_type or "unknown"
            data_types[dt] = data_types.get(dt, 0) + 1

        return {
            "total_datasets": total,
            "active": active,
            "archived": archived,
            "by_data_type": data_types,
        }

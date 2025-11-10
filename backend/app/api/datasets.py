"""
Dataset Manager API

REST API for dataset management with version control,
file lineage tracking, and metadata management.
"""

from fastapi import APIRouter, HTTPException, Depends, Query
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime

from sqlalchemy.orm import Session
from app.database import get_db
from app.services.dataset_manager import DatasetManager

router = APIRouter(prefix="/api/datasets", tags=["datasets"])


# Pydantic models
class DatasetCreate(BaseModel):
    """Request to create dataset"""

    project_id: str
    name: str
    description: Optional[str] = None
    data_type: Optional[str] = None
    file_format: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = Field(default_factory=dict)
    tags: Optional[List[str]] = Field(default_factory=list)


class DatasetUpdate(BaseModel):
    """Request to update dataset"""

    name: Optional[str] = None
    description: Optional[str] = None
    data_type: Optional[str] = None
    file_format: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
    status: Optional[str] = None


class DatasetResponse(BaseModel):
    """Dataset response"""

    id: str
    project_id: str
    name: str
    description: Optional[str]
    version: int
    data_type: Optional[str]
    file_format: Optional[str]
    metadata: Dict[str, Any]
    status: str
    is_public: bool
    created_at: datetime
    updated_at: datetime
    created_by: Optional[str]
    parent_id: Optional[str]
    file_count: int
    tag_names: List[str]

    class Config:
        from_attributes = True


class FileAddRequest(BaseModel):
    """Request to add file to dataset"""

    file_path: str
    role: str = "primary"
    description: Optional[str] = None
    compute_hash: bool = True


class FileResponse(BaseModel):
    """File response"""

    id: str
    file_path: str
    file_name: str
    file_type: Optional[str]
    file_size: Optional[int]
    role: Optional[str]
    description: Optional[str]
    md5_hash: Optional[str]
    sha256_hash: Optional[str]
    added_at: datetime

    class Config:
        from_attributes = True


class LineageAddRequest(BaseModel):
    """Request to add lineage record"""

    operation_type: str
    source_file_id: Optional[str] = None
    output_file_id: Optional[str] = None
    operation_id: Optional[str] = None
    operation_params: Optional[Dict[str, Any]] = Field(default_factory=dict)
    tool_name: Optional[str] = None
    tool_version: Optional[str] = None


class LineageResponse(BaseModel):
    """Lineage response"""

    id: str
    dataset_id: str
    operation_type: str
    operation_id: Optional[str]
    source_file_id: Optional[str]
    output_file_id: Optional[str]
    operation_params: Dict[str, Any]
    tool_name: Optional[str]
    tool_version: Optional[str]
    executed_at: datetime

    class Config:
        from_attributes = True


class TagCreate(BaseModel):
    """Request to create tag"""

    name: str
    color: Optional[str] = "#3b82f6"
    description: Optional[str] = None


class TagResponse(BaseModel):
    """Tag response"""

    id: str
    name: str
    color: str
    description: Optional[str]
    created_at: datetime
    dataset_count: int = 0

    class Config:
        from_attributes = True


# Helper functions
def format_dataset_response(dataset, db: Session) -> Dict[str, Any]:
    """Format dataset for response"""
    return {
        "id": dataset.id,
        "project_id": dataset.project_id,
        "name": dataset.name,
        "description": dataset.description,
        "version": dataset.version,
        "data_type": dataset.data_type,
        "file_format": dataset.file_format,
        "metadata": dataset.metadata,
        "status": dataset.status,
        "is_public": dataset.is_public,
        "created_at": dataset.created_at,
        "updated_at": dataset.updated_at,
        "created_by": dataset.created_by,
        "parent_id": dataset.parent_id,
        "file_count": dataset.files.count(),
        "tag_names": [tag.name for tag in dataset.tags.all()],
    }


def format_file_response(file_entry) -> Dict[str, Any]:
    """Format file entry for response"""
    return {
        "id": file_entry.id,
        "file_path": file_entry.file_path,
        "file_name": file_entry.file_name,
        "file_type": file_entry.file_type,
        "file_size": file_entry.file_size,
        "role": file_entry.role,
        "description": file_entry.description,
        "md5_hash": file_entry.md5_hash,
        "sha256_hash": file_entry.sha256_hash,
        "added_at": file_entry.added_at,
    }


def format_tag_response(tag, db: Session) -> Dict[str, Any]:
    """Format tag for response"""
    return {
        "id": tag.id,
        "name": tag.name,
        "color": tag.color,
        "description": tag.description,
        "created_at": tag.created_at,
        "dataset_count": tag.datasets.count(),
    }


# Endpoints
@router.post("/", response_model=DatasetResponse)
async def create_dataset(request: DatasetCreate, db: Session = Depends(get_db)):
    """
    Create a new dataset

    Creates a dataset with metadata, tags, and initial configuration.
    """
    try:
        manager = DatasetManager(db)
        dataset = manager.create_dataset(
            project_id=request.project_id,
            name=request.name,
            description=request.description,
            data_type=request.data_type,
            file_format=request.file_format,
            metadata=request.metadata,
            tags=request.tags,
        )

        return format_dataset_response(dataset, db)

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{dataset_id}", response_model=DatasetResponse)
async def get_dataset(dataset_id: str, db: Session = Depends(get_db)):
    """
    Get dataset by ID

    Returns complete dataset information including file count and tags.
    """
    manager = DatasetManager(db)
    dataset = manager.get_dataset(dataset_id)

    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    return format_dataset_response(dataset, db)


@router.get("/", response_model=List[DatasetResponse])
async def list_datasets(
    project_id: Optional[str] = Query(None),
    data_type: Optional[str] = Query(None),
    status: Optional[str] = Query(None),
    tags: Optional[str] = Query(None, description="Comma-separated tag names"),
    search: Optional[str] = Query(None),
    skip: int = Query(0, ge=0),
    limit: int = Query(100, ge=1, le=1000),
    db: Session = Depends(get_db),
):
    """
    List datasets with filters

    Supports filtering by project, data type, status, tags, and search query.
    """
    manager = DatasetManager(db)

    tag_list = tags.split(",") if tags else None

    datasets = manager.list_datasets(
        project_id=project_id,
        data_type=data_type,
        status=status,
        tags=tag_list,
        search_query=search,
        skip=skip,
        limit=limit,
    )

    return [format_dataset_response(ds, db) for ds in datasets]


@router.put("/{dataset_id}", response_model=DatasetResponse)
async def update_dataset(
    dataset_id: str, request: DatasetUpdate, db: Session = Depends(get_db)
):
    """
    Update dataset metadata

    Updates basic information and metadata fields.
    """
    manager = DatasetManager(db)
    dataset = manager.update_dataset(
        dataset_id=dataset_id,
        name=request.name,
        description=request.description,
        data_type=request.data_type,
        file_format=request.file_format,
        metadata=request.metadata,
        status=request.status,
    )

    if not dataset:
        raise HTTPException(status_code=404, detail="Dataset not found")

    return format_dataset_response(dataset, db)


@router.delete("/{dataset_id}")
async def delete_dataset(dataset_id: str, db: Session = Depends(get_db)):
    """
    Delete dataset (soft delete)

    Sets dataset status to 'deleted' without removing from database.
    """
    manager = DatasetManager(db)
    success = manager.delete_dataset(dataset_id)

    if not success:
        raise HTTPException(status_code=404, detail="Dataset not found")

    return {"message": "Dataset deleted successfully"}


@router.post("/{dataset_id}/archive")
async def archive_dataset(dataset_id: str, db: Session = Depends(get_db)):
    """
    Archive dataset

    Sets dataset status to 'archived'.
    """
    manager = DatasetManager(db)
    success = manager.archive_dataset(dataset_id)

    if not success:
        raise HTTPException(status_code=404, detail="Dataset not found")

    return {"message": "Dataset archived successfully"}


# Version control endpoints
@router.post("/{dataset_id}/versions", response_model=DatasetResponse)
async def create_version(
    dataset_id: str, description: Optional[str] = None, db: Session = Depends(get_db)
):
    """
    Create a new version of dataset

    Creates a new version with incremented version number,
    copying files and tags from parent.
    """
    manager = DatasetManager(db)
    new_version = manager.create_version(
        dataset_id=dataset_id,
        description=description,
    )

    if not new_version:
        raise HTTPException(status_code=404, detail="Dataset not found")

    return format_dataset_response(new_version, db)


@router.get("/{dataset_id}/versions", response_model=List[DatasetResponse])
async def get_versions(dataset_id: str, db: Session = Depends(get_db)):
    """
    Get all versions of dataset

    Returns all versions ordered by version number.
    """
    manager = DatasetManager(db)
    versions = manager.get_versions(dataset_id)

    return [format_dataset_response(v, db) for v in versions]


# File management endpoints
@router.post("/{dataset_id}/files", response_model=FileResponse)
async def add_file(
    dataset_id: str, request: FileAddRequest, db: Session = Depends(get_db)
):
    """
    Add file to dataset

    Adds a file with optional hash computation for integrity checking.
    """
    manager = DatasetManager(db)
    file_entry = manager.add_file(
        dataset_id=dataset_id,
        file_path=request.file_path,
        role=request.role,
        description=request.description,
        compute_hash=request.compute_hash,
    )

    if not file_entry:
        raise HTTPException(
            status_code=400,
            detail="Failed to add file (dataset not found or file doesn't exist)",
        )

    return format_file_response(file_entry)


@router.delete("/{dataset_id}/files/{file_id}")
async def remove_file(dataset_id: str, file_id: str, db: Session = Depends(get_db)):
    """
    Remove file from dataset

    Removes file association without deleting the physical file.
    """
    manager = DatasetManager(db)
    success = manager.remove_file(dataset_id, file_id)

    if not success:
        raise HTTPException(status_code=404, detail="Dataset or file not found")

    return {"message": "File removed from dataset"}


@router.get("/{dataset_id}/files", response_model=List[FileResponse])
async def list_files(dataset_id: str, db: Session = Depends(get_db)):
    """
    List all files in dataset

    Returns all files with their metadata and hashes.
    """
    manager = DatasetManager(db)
    files = manager.list_files(dataset_id)

    return [format_file_response(f) for f in files]


# Lineage tracking endpoints
@router.post("/{dataset_id}/lineage", response_model=LineageResponse)
async def add_lineage(
    dataset_id: str, request: LineageAddRequest, db: Session = Depends(get_db)
):
    """
    Add lineage record

    Records a data processing operation for provenance tracking.
    """
    manager = DatasetManager(db)
    lineage = manager.add_lineage(
        dataset_id=dataset_id,
        operation_type=request.operation_type,
        source_file_id=request.source_file_id,
        output_file_id=request.output_file_id,
        operation_id=request.operation_id,
        operation_params=request.operation_params,
        tool_name=request.tool_name,
        tool_version=request.tool_version,
    )

    if not lineage:
        raise HTTPException(status_code=404, detail="Dataset not found")

    return LineageResponse.from_orm(lineage)


@router.get("/{dataset_id}/lineage", response_model=List[LineageResponse])
async def get_lineage(dataset_id: str, db: Session = Depends(get_db)):
    """
    Get lineage history

    Returns all lineage records for the dataset ordered by time.
    """
    manager = DatasetManager(db)
    lineage_records = manager.get_lineage(dataset_id)

    return [LineageResponse.from_orm(lr) for lr in lineage_records]


@router.get("/files/{file_id}/lineage")
async def trace_file_lineage(file_id: str, db: Session = Depends(get_db)):
    """
    Trace file lineage

    Returns complete lineage tree showing how file was created
    and where it was used.
    """
    manager = DatasetManager(db)
    lineage_tree = manager.trace_file_lineage(file_id)

    if not lineage_tree:
        raise HTTPException(status_code=404, detail="File not found")

    return lineage_tree


# Tag management endpoints
@router.post("/{dataset_id}/tags")
async def add_tag(
    dataset_id: str,
    tag_name: str = Query(..., description="Tag name to add"),
    db: Session = Depends(get_db),
):
    """
    Add tag to dataset

    Creates tag if it doesn't exist and adds it to dataset.
    """
    manager = DatasetManager(db)
    success = manager.add_tag_to_dataset(dataset_id, tag_name)

    if not success:
        raise HTTPException(status_code=404, detail="Dataset not found")

    return {"message": f"Tag '{tag_name}' added to dataset"}


@router.delete("/{dataset_id}/tags/{tag_id}")
async def remove_tag(dataset_id: str, tag_id: str, db: Session = Depends(get_db)):
    """
    Remove tag from dataset

    Removes tag association without deleting the tag.
    """
    manager = DatasetManager(db)
    success = manager.remove_tag_from_dataset(dataset_id, tag_id)

    if not success:
        raise HTTPException(status_code=404, detail="Dataset or tag not found")

    return {"message": "Tag removed from dataset"}


@router.get("/tags/list", response_model=List[TagResponse])
async def list_tags(db: Session = Depends(get_db)):
    """
    List all tags

    Returns all tags with dataset counts.
    """
    manager = DatasetManager(db)
    tags = manager.list_tags()

    return [format_tag_response(tag, db) for tag in tags]


@router.post("/tags", response_model=TagResponse)
async def create_tag(request: TagCreate, db: Session = Depends(get_db)):
    """
    Create a new tag

    Creates a tag with specified color and description.
    """
    manager = DatasetManager(db)
    tag = manager.get_or_create_tag(request.name, request.color)

    if request.description:
        tag.description = request.description
        db.commit()
        db.refresh(tag)

    return format_tag_response(tag, db)


# Statistics endpoint
@router.get("/stats/summary")
async def get_statistics(
    project_id: Optional[str] = Query(None), db: Session = Depends(get_db)
):
    """
    Get dataset statistics

    Returns counts and breakdowns by data type.
    """
    manager = DatasetManager(db)
    stats = manager.get_statistics(project_id)

    return stats

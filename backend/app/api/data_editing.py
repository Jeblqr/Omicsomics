"""
API endpoints for data editing functionality.

Provides in-place data operations with preview, validation, and undo/redo.
"""

from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from ..database import get_db
from ..services.data_editing import DataEditingService
from ..models.data_editing import EditOperationType, EditStatus
from ..models.user import User
from ..core.auth import get_current_user


router = APIRouter(prefix="/api/data-editing", tags=["data-editing"])


# Pydantic models
class EditSessionCreate(BaseModel):
    """Request model for creating edit session"""

    name: str = Field(..., min_length=1, max_length=255)
    description: Optional[str] = None
    file_id: int
    project_id: Optional[int] = None
    dataset_id: Optional[int] = None


class OperationAdd(BaseModel):
    """Request model for adding operation"""

    operation_type: EditOperationType
    parameters: dict
    description: Optional[str] = None


class PreviewRequest(BaseModel):
    """Request model for preview"""

    sample_size: int = Field(default=100, ge=10, le=1000)


class EditSessionResponse(BaseModel):
    """Response model for edit session"""

    id: int
    session_key: str
    name: str
    description: Optional[str]
    user_id: int
    project_id: Optional[int]
    file_id: int
    dataset_id: Optional[int]
    status: str
    operations: List[dict]
    current_operation_index: int
    preview_data: Optional[dict]
    preview_summary: Optional[dict]
    validation_errors: Optional[List[dict]]
    validation_warnings: Optional[List[dict]]
    is_valid: bool
    backup_path: Optional[str]
    output_file_id: Optional[int]
    created_at: str
    updated_at: Optional[str]
    applied_at: Optional[str]

    class Config:
        from_attributes = True


class EditSessionListItem(BaseModel):
    """Simplified edit session for list view"""

    id: int
    session_key: str
    name: str
    user_id: int
    file_id: int
    status: str
    operation_count: int
    is_valid: bool
    created_at: str
    applied_at: Optional[str]


class PreviewResponse(BaseModel):
    """Preview response"""

    preview: dict
    summary: dict
    validation: dict


# Endpoints
@router.post("/sessions", response_model=EditSessionResponse)
def create_edit_session(
    data: EditSessionCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Create a new edit session.

    This creates a backup of the original file and initializes an edit session
    for in-place data operations.

    - **name**: Session name
    - **file_id**: ID of file to edit (required)
    - **project_id**: Optional project ID
    - **dataset_id**: Optional dataset ID
    """
    service = DataEditingService(db)

    try:
        session = service.create_edit_session(
            user_id=current_user.id,
            file_id=data.file_id,
            name=data.name,
            project_id=data.project_id,
            dataset_id=data.dataset_id,
            description=data.description,
        )

        return EditSessionResponse(**session.to_dict())

    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get("/sessions/{session_id}", response_model=EditSessionResponse)
def get_edit_session(
    session_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Get edit session by ID.

    Returns detailed information including all operations, preview data,
    and validation results.
    """
    service = DataEditingService(db)
    session = service.get_edit_session(session_id, user_id=current_user.id)

    if not session:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Edit session not found"
        )

    return EditSessionResponse(**session.to_dict())


@router.get("/sessions/key/{session_key}", response_model=EditSessionResponse)
def get_edit_session_by_key(
    session_key: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Get edit session by session key"""
    service = DataEditingService(db)
    session = service.get_edit_session_by_key(session_key, user_id=current_user.id)

    if not session:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Edit session not found"
        )

    return EditSessionResponse(**session.to_dict())


@router.get("/sessions", response_model=List[EditSessionListItem])
def list_edit_sessions(
    project_id: Optional[int] = None,
    status: Optional[EditStatus] = None,
    file_id: Optional[int] = None,
    limit: int = 50,
    offset: int = 0,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    List edit sessions for current user.

    - **project_id**: Filter by project
    - **status**: Filter by status (draft, previewing, applying, applied, failed, reverted)
    - **file_id**: Filter by file
    - **limit**: Max results (default 50)
    - **offset**: Pagination offset
    """
    service = DataEditingService(db)
    sessions = service.list_edit_sessions(
        user_id=current_user.id,
        project_id=project_id,
        status=status,
        file_id=file_id,
        limit=limit,
        offset=offset,
    )

    return [
        EditSessionListItem(
            id=session.id,
            session_key=session.session_key,
            name=session.name,
            user_id=session.user_id,
            file_id=session.file_id,
            status=session.status.value if session.status else "unknown",
            operation_count=len(session.operations or []),
            is_valid=session.is_valid,
            created_at=session.created_at.isoformat() if session.created_at else None,
            applied_at=session.applied_at.isoformat() if session.applied_at else None,
        )
        for session in sessions
    ]


@router.post("/sessions/{session_id}/operations", response_model=EditSessionResponse)
def add_operation(
    session_id: int,
    data: OperationAdd,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Add an operation to the edit session.

    Operations are queued and can be previewed before applying.

    **Operation Types**:
    - Row operations: delete_rows, filter_rows, sort_rows, deduplicate
    - Column operations: add_column, delete_column, rename_column, reorder_columns
    - Value operations: replace_values, fill_missing, transform_values, normalize
    - Type operations: convert_type, cast_column
    - Merge operations: join_data, concat_data
    - Custom: custom_script

    **Example - Filter Rows**:
    ```json
    {
      "operation_type": "filter_rows",
      "parameters": {
        "column": "expression",
        "operator": "greater_than",
        "value": 2.0
      }
    }
    ```

    **Example - Fill Missing**:
    ```json
    {
      "operation_type": "fill_missing",
      "parameters": {
        "column": "sample_1",
        "method": "mean"
      }
    }
    ```
    """
    service = DataEditingService(db)

    try:
        session = service.add_operation(
            session_id=session_id,
            operation_type=data.operation_type,
            parameters=data.parameters,
            description=data.description,
        )

        return EditSessionResponse(**session.to_dict())

    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.delete("/sessions/{session_id}/operations/{operation_index}")
def remove_operation(
    session_id: int,
    operation_index: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Remove an operation from the session.

    This allows you to undo operations before applying them.
    """
    service = DataEditingService(db)

    try:
        session = service.remove_operation(session_id, operation_index)
        return EditSessionResponse(**session.to_dict())

    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.post("/sessions/{session_id}/preview", response_model=PreviewResponse)
def preview_operations(
    session_id: int,
    data: PreviewRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Preview the result of applying all operations.

    This loads the original file, applies all operations, and returns:
    - Sample of result data (configurable size)
    - Summary statistics (shape, types, missing values)
    - Validation errors and warnings

    **No data is modified** - this is a dry-run preview only.

    - **sample_size**: Number of rows to preview (default 100, max 1000)
    """
    service = DataEditingService(db)

    try:
        result = service.preview_operations(
            session_id=session_id, sample_size=data.sample_size
        )

        return PreviewResponse(**result)

    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.post("/sessions/{session_id}/apply", response_model=EditSessionResponse)
def apply_operations(
    session_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Apply all operations and save the result.

    This:
    1. Validates the operations (must run preview first)
    2. Applies all operations to the original file
    3. Creates a new file with the edited data
    4. Updates the session status to 'applied'

    **Prerequisites**:
    - Session must have been previewed (is_valid = true)
    - No validation errors

    **Result**:
    - Original file is preserved (backup created)
    - New file is created with edited data (output_file_id)
    """
    service = DataEditingService(db)

    try:
        session = service.apply_operations(session_id)
        return EditSessionResponse(**session.to_dict())

    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.post("/sessions/{session_id}/revert", response_model=EditSessionResponse)
def revert_operations(
    session_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Revert session to original file.

    This restores the file from backup and marks the session as reverted.
    """
    service = DataEditingService(db)

    try:
        session = service.revert_operations(session_id)
        return EditSessionResponse(**session.to_dict())

    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.delete("/sessions/{session_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_edit_session(
    session_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Delete an edit session and cleanup resources.

    This will:
    - Delete the backup file
    - Delete all operation records
    - Delete the session record

    **Note**: The output file (if applied) is NOT deleted.
    """
    service = DataEditingService(db)

    if not service.delete_edit_session(session_id, user_id=current_user.id):
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Edit session not found"
        )


@router.get("/operation-types")
def get_operation_types():
    """
    Get available operation types and their parameters.

    Returns a list of all supported operation types with:
    - Operation type name
    - Description
    - Required parameters
    - Optional parameters
    - Examples
    """
    return {
        "operation_types": [
            {
                "type": "delete_rows",
                "category": "row",
                "description": "Delete specific rows by index",
                "parameters": {"required": ["indices"], "optional": []},
                "example": {
                    "operation_type": "delete_rows",
                    "parameters": {"indices": [0, 5, 10]},
                },
            },
            {
                "type": "filter_rows",
                "category": "row",
                "description": "Filter rows based on condition",
                "parameters": {
                    "required": ["column", "operator", "value"],
                    "optional": [],
                },
                "operators": [
                    "equals",
                    "not_equals",
                    "greater_than",
                    "less_than",
                    "contains",
                ],
                "example": {
                    "operation_type": "filter_rows",
                    "parameters": {
                        "column": "expression",
                        "operator": "greater_than",
                        "value": 2.0,
                    },
                },
            },
            {
                "type": "sort_rows",
                "category": "row",
                "description": "Sort rows by columns",
                "parameters": {"required": ["columns"], "optional": ["ascending"]},
                "example": {
                    "operation_type": "sort_rows",
                    "parameters": {"columns": ["gene_name"], "ascending": True},
                },
            },
            {
                "type": "deduplicate",
                "category": "row",
                "description": "Remove duplicate rows",
                "parameters": {"required": [], "optional": ["columns"]},
                "example": {
                    "operation_type": "deduplicate",
                    "parameters": {"columns": ["gene_id"]},
                },
            },
            {
                "type": "add_column",
                "category": "column",
                "description": "Add a new column",
                "parameters": {
                    "required": ["column_name"],
                    "optional": ["default_value"],
                },
                "example": {
                    "operation_type": "add_column",
                    "parameters": {"column_name": "log2_fc", "default_value": 0.0},
                },
            },
            {
                "type": "delete_column",
                "category": "column",
                "description": "Delete columns",
                "parameters": {"required": ["columns"], "optional": []},
                "example": {
                    "operation_type": "delete_column",
                    "parameters": {"columns": ["temp_col"]},
                },
            },
            {
                "type": "rename_column",
                "category": "column",
                "description": "Rename a column",
                "parameters": {"required": ["old_name", "new_name"], "optional": []},
                "example": {
                    "operation_type": "rename_column",
                    "parameters": {"old_name": "gene", "new_name": "gene_symbol"},
                },
            },
            {
                "type": "replace_values",
                "category": "value",
                "description": "Replace specific values in a column",
                "parameters": {
                    "required": ["column", "old_value", "new_value"],
                    "optional": [],
                },
                "example": {
                    "operation_type": "replace_values",
                    "parameters": {
                        "column": "status",
                        "old_value": "NA",
                        "new_value": "Unknown",
                    },
                },
            },
            {
                "type": "fill_missing",
                "category": "value",
                "description": "Fill missing values",
                "parameters": {
                    "required": ["column"],
                    "optional": ["fill_value", "method"],
                },
                "methods": ["value", "mean", "median", "forward_fill", "backward_fill"],
                "example": {
                    "operation_type": "fill_missing",
                    "parameters": {"column": "expression", "method": "mean"},
                },
            },
            {
                "type": "convert_type",
                "category": "type",
                "description": "Convert column data type",
                "parameters": {"required": ["column", "target_type"], "optional": []},
                "types": ["int", "float", "string", "datetime", "bool"],
                "example": {
                    "operation_type": "convert_type",
                    "parameters": {"column": "count", "target_type": "int"},
                },
            },
        ]
    }


@router.get("/")
def get_editing_info():
    """Get data editing service information"""
    return {
        "service": "Data Editing",
        "version": "1.0",
        "features": [
            "In-place data operations",
            "Preview before apply",
            "Operation validation",
            "Undo/redo support",
            "Automatic backup",
            "Multiple operation types",
            "Type conversion",
            "Missing value handling",
        ],
        "operation_categories": [
            "row_operations",
            "column_operations",
            "value_operations",
            "type_operations",
        ],
        "endpoints": {
            "create_session": "POST /api/data-editing/sessions",
            "get_session": "GET /api/data-editing/sessions/{session_id}",
            "list_sessions": "GET /api/data-editing/sessions",
            "add_operation": "POST /api/data-editing/sessions/{session_id}/operations",
            "remove_operation": "DELETE /api/data-editing/sessions/{session_id}/operations/{index}",
            "preview": "POST /api/data-editing/sessions/{session_id}/preview",
            "apply": "POST /api/data-editing/sessions/{session_id}/apply",
            "revert": "POST /api/data-editing/sessions/{session_id}/revert",
            "delete_session": "DELETE /api/data-editing/sessions/{session_id}",
            "operation_types": "GET /api/data-editing/operation-types",
        },
    }

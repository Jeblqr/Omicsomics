"""
Custom Script Tools API

Provides REST API endpoints for managing and executing custom scripts.
"""

from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from ..core.deps import get_current_user, get_db
from ..models.custom_scripts import ScriptLanguage, ScriptStatus, ScriptVisibility
from ..models.user import User
from ..services.custom_scripts import CustomScriptService

router = APIRouter(prefix="/api/custom-scripts")


# Pydantic Models


class ScriptCreate(BaseModel):
    """Model for creating a script."""

    name: str = Field(..., min_length=1, max_length=200)
    description: str = Field(..., min_length=1)
    language: ScriptLanguage
    script_content: str = Field(..., min_length=1)
    entry_point: Optional[str] = None
    parameters_schema: Optional[Dict[str, Any]] = None
    requirements: Optional[List[str]] = None
    timeout: int = Field(default=300, ge=1, le=3600)
    max_memory: Optional[int] = Field(default=None, ge=100, le=8192)
    visibility: ScriptVisibility = ScriptVisibility.PRIVATE
    category: Optional[str] = None
    tags: Optional[List[str]] = None


class ScriptUpdate(BaseModel):
    """Model for updating a script."""

    name: Optional[str] = Field(None, min_length=1, max_length=200)
    description: Optional[str] = None
    script_content: Optional[str] = None
    entry_point: Optional[str] = None
    parameters_schema: Optional[Dict[str, Any]] = None
    requirements: Optional[List[str]] = None
    timeout: Optional[int] = Field(None, ge=1, le=3600)
    max_memory: Optional[int] = Field(None, ge=100, le=8192)
    visibility: Optional[ScriptVisibility] = None
    category: Optional[str] = None
    tags: Optional[List[str]] = None


class ScriptExecuteRequest(BaseModel):
    """Model for executing a script."""

    parameters: Optional[Dict[str, Any]] = None
    input_file_ids: Optional[List[int]] = None
    description: Optional[str] = None


class ScriptResponse(BaseModel):
    """Model for script response."""

    id: int
    script_key: str
    name: str
    description: str
    user_id: int
    language: ScriptLanguage
    script_content: str
    entry_point: Optional[str]
    parameters_schema: Optional[Dict[str, Any]]
    requirements: List[str]
    timeout: int
    max_memory: Optional[int]
    visibility: ScriptVisibility
    category: Optional[str]
    tags: List[str]
    is_verified: bool
    total_executions: int
    successful_executions: int
    failed_executions: int
    avg_duration: Optional[float]
    last_executed_at: Optional[str]
    created_at: str
    updated_at: str

    class Config:
        from_attributes = True


class ScriptListItem(BaseModel):
    """Model for script list item (without script content)."""

    id: int
    script_key: str
    name: str
    description: str
    user_id: int
    language: ScriptLanguage
    visibility: ScriptVisibility
    category: Optional[str]
    tags: List[str]
    is_verified: bool
    total_executions: int
    successful_executions: int
    last_executed_at: Optional[str]
    created_at: str

    class Config:
        from_attributes = True


class ExecutionResponse(BaseModel):
    """Model for execution response."""

    id: int
    execution_key: str
    script_id: int
    user_id: int
    parameters: Dict[str, Any]
    input_file_ids: List[int]
    output_file_ids: List[int]
    status: ScriptStatus
    description: Optional[str]
    result_data: Optional[Dict[str, Any]]
    output_text: Optional[str]
    error_text: Optional[str]
    exit_code: Optional[int]
    duration: Optional[float]
    memory_usage: Optional[int]
    started_at: Optional[str]
    completed_at: Optional[str]
    created_at: str

    class Config:
        from_attributes = True


class ScriptListResponse(BaseModel):
    """Model for paginated script list response."""

    items: List[ScriptListItem]
    total: int
    skip: int
    limit: int


class ExecutionListResponse(BaseModel):
    """Model for paginated execution list response."""

    items: List[ExecutionResponse]
    total: int
    skip: int
    limit: int


class ParameterValidationRequest(BaseModel):
    """Model for parameter validation request."""

    parameters: Dict[str, Any]


class ParameterValidationResponse(BaseModel):
    """Model for parameter validation response."""

    valid: bool
    errors: Optional[List[str]] = None


# API Endpoints


@router.post(
    "/scripts", response_model=ScriptResponse, status_code=status.HTTP_201_CREATED
)
def create_script(
    script_data: ScriptCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Create a new custom script.

    - **name**: Script name (1-200 characters)
    - **description**: Script description
    - **language**: Script language (python, r, bash)
    - **script_content**: The actual script code
    - **entry_point**: Entry point function/method (optional)
    - **parameters_schema**: JSON Schema for parameters validation (optional)
    - **requirements**: List of dependencies (optional)
    - **timeout**: Execution timeout in seconds (1-3600, default: 300)
    - **max_memory**: Maximum memory in MB (100-8192, optional)
    - **visibility**: Script visibility (private, project, public, default: private)
    - **category**: Script category (optional)
    - **tags**: Script tags (optional)
    """
    service = CustomScriptService(db)

    try:
        script = service.create_script(
            user_id=current_user.id,
            name=script_data.name,
            description=script_data.description,
            language=script_data.language,
            script_content=script_data.script_content,
            entry_point=script_data.entry_point,
            parameters_schema=script_data.parameters_schema,
            requirements=script_data.requirements,
            timeout=script_data.timeout,
            max_memory=script_data.max_memory,
            visibility=script_data.visibility,
            category=script_data.category,
            tags=script_data.tags,
        )
        return script
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get("/scripts/{script_id}", response_model=ScriptResponse)
def get_script(
    script_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Get a script by ID.

    Returns full script details including script content.
    Only accessible if user is owner or script is public.
    """
    service = CustomScriptService(db)
    script = service.get_script(script_id, current_user.id)

    if not script:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Script not found: {script_id}",
        )

    return script


@router.get("/scripts", response_model=ScriptListResponse)
def list_scripts(
    language: Optional[ScriptLanguage] = None,
    category: Optional[str] = None,
    visibility: Optional[ScriptVisibility] = None,
    tags: Optional[str] = None,
    search: Optional[str] = None,
    verified_only: bool = False,
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    List scripts with filters.

    - **language**: Filter by language (python, r, bash)
    - **category**: Filter by category
    - **visibility**: Filter by visibility (private, project, public)
    - **tags**: Comma-separated tags to filter by
    - **search**: Search in name and description
    - **verified_only**: Show only verified scripts
    - **skip**: Number of records to skip (pagination)
    - **limit**: Maximum number of records to return (max: 100)
    """
    if limit > 100:
        limit = 100

    service = CustomScriptService(db)
    tag_list = [t.strip() for t in tags.split(",")] if tags else None

    scripts, total = service.list_scripts(
        user_id=current_user.id,
        language=language,
        category=category,
        visibility=visibility,
        tags=tag_list,
        search=search,
        verified_only=verified_only,
        skip=skip,
        limit=limit,
    )

    return ScriptListResponse(
        items=scripts,
        total=total,
        skip=skip,
        limit=limit,
    )


@router.put("/scripts/{script_id}", response_model=ScriptResponse)
def update_script(
    script_id: int,
    script_data: ScriptUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Update a script (owner only).

    Only the script owner can update the script.
    All fields are optional - only provided fields will be updated.
    """
    service = CustomScriptService(db)

    # Convert to dict and remove None values
    update_fields = {k: v for k, v in script_data.dict().items() if v is not None}

    try:
        script = service.update_script(script_id, current_user.id, **update_fields)

        if not script:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Script not found or access denied: {script_id}",
            )

        return script
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.delete("/scripts/{script_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_script(
    script_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Delete a script (owner only).

    Only the script owner can delete the script.
    This will also delete all associated executions.
    """
    service = CustomScriptService(db)
    success = service.delete_script(script_id, current_user.id)

    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Script not found or access denied: {script_id}",
        )


@router.post(
    "/scripts/{script_id}/execute",
    response_model=ExecutionResponse,
    status_code=status.HTTP_202_ACCEPTED,
)
def execute_script(
    script_id: int,
    exec_request: ScriptExecuteRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Execute a script with given parameters.

    - **parameters**: Script parameters (must match parameters_schema if defined)
    - **input_file_ids**: List of input file IDs
    - **description**: Execution description (optional)

    Returns execution record with status "pending" or "running".
    Use GET /executions/{execution_id} to check execution status and results.
    """
    service = CustomScriptService(db)

    try:
        execution = service.execute_script(
            script_id=script_id,
            user_id=current_user.id,
            parameters=exec_request.parameters,
            input_file_ids=exec_request.input_file_ids,
            description=exec_request.description,
        )
        return execution
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.post(
    "/scripts/{script_id}/validate-parameters",
    response_model=ParameterValidationResponse,
)
def validate_parameters(
    script_id: int,
    validation_request: ParameterValidationRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Validate parameters against script's parameters_schema.

    Useful for checking parameters before execution.
    Returns validation result with any errors.
    """
    service = CustomScriptService(db)
    script = service.get_script(script_id, current_user.id)

    if not script:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Script not found: {script_id}",
        )

    if not script.parameters_schema:
        return ParameterValidationResponse(valid=True)

    try:
        service._validate_parameters(
            validation_request.parameters, script.parameters_schema
        )
        return ParameterValidationResponse(valid=True)
    except ValueError as e:
        return ParameterValidationResponse(valid=False, errors=[str(e)])


@router.get("/executions/{execution_id}", response_model=ExecutionResponse)
def get_execution(
    execution_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Get execution by ID.

    Returns execution details including status, results, output, and errors.
    """
    service = CustomScriptService(db)
    execution = service.get_execution(execution_id, current_user.id)

    if not execution:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Execution not found: {execution_id}",
        )

    return execution


@router.get("/executions", response_model=ExecutionListResponse)
def list_executions(
    script_id: Optional[int] = None,
    status_filter: Optional[ScriptStatus] = None,
    skip: int = 0,
    limit: int = 50,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    List executions with filters.

    - **script_id**: Filter by script ID
    - **status**: Filter by status (pending, running, completed, failed, cancelled)
    - **skip**: Number of records to skip (pagination)
    - **limit**: Maximum number of records to return (max: 100)
    """
    if limit > 100:
        limit = 100

    service = CustomScriptService(db)
    executions, total = service.list_executions(
        script_id=script_id,
        user_id=current_user.id,
        status=status_filter,
        skip=skip,
        limit=limit,
    )

    return ExecutionListResponse(
        items=executions,
        total=total,
        skip=skip,
        limit=limit,
    )


@router.post(
    "/executions/{execution_id}/cancel", status_code=status.HTTP_204_NO_CONTENT
)
def cancel_execution(
    execution_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Cancel a running execution.

    Only pending or running executions can be cancelled.
    """
    service = CustomScriptService(db)
    success = service.cancel_execution(execution_id, current_user.id)

    if not success:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Execution cannot be cancelled (not found, access denied, or already completed)",
        )


@router.get("/languages", response_model=List[str])
def get_supported_languages():
    """
    Get list of supported script languages.
    """
    return [lang.value for lang in ScriptLanguage]


@router.get("/categories", response_model=List[str])
def get_categories(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """
    Get list of all script categories in use.
    """
    from sqlalchemy import func
    from ..models.custom_scripts import CustomScript

    categories = (
        db.query(CustomScript.category)
        .filter(CustomScript.category.isnot(None))
        .distinct()
        .all()
    )

    return [cat[0] for cat in categories if cat[0]]


@router.get("/templates", response_model=List[Dict[str, Any]])
def get_script_templates():
    """
    Get script templates for each language.

    Returns example scripts with parameter schemas to help users get started.
    """
    templates = [
        {
            "language": "python",
            "name": "Data Analysis Template",
            "description": "Basic data analysis script with pandas",
            "script_content": """import pandas as pd
import json
import sys
import os

def main(params):
    # Load parameters
    input_file = params.get('input_file')
    output_file = params.get('output_file', 'output.csv')
    
    # Read data
    df = pd.read_csv(input_file)
    
    # Perform analysis
    summary = df.describe()
    
    # Save results
    os.makedirs('output', exist_ok=True)
    summary.to_csv(f'output/{output_file}')
    
    # Return results
    return {
        'rows': len(df),
        'columns': len(df.columns),
        'output_file': output_file
    }

if __name__ == '__main__':
    # Load parameters from JSON file
    with open(sys.argv[2], 'r') as f:
        params = json.load(f)
    
    result = main(params)
    print(json.dumps(result))
""",
            "parameters_schema": {
                "type": "object",
                "properties": {
                    "input_file": {
                        "type": "string",
                        "description": "Input CSV file path",
                    },
                    "output_file": {"type": "string", "default": "output.csv"},
                },
                "required": ["input_file"],
            },
        },
        {
            "language": "r",
            "name": "Statistical Analysis Template",
            "description": "Basic statistical analysis script",
            "script_content": """library(jsonlite)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(args[2])

# Load data
data <- read.csv(params$input_file)

# Perform analysis
summary_stats <- summary(data)

# Save results
dir.create('output', showWarnings = FALSE)
write.csv(summary_stats, file.path('output', params$output_file))

# Return results
result <- list(
  rows = nrow(data),
  columns = ncol(data),
  output_file = params$output_file
)

cat(toJSON(result))
""",
            "parameters_schema": {
                "type": "object",
                "properties": {
                    "input_file": {"type": "string"},
                    "output_file": {"type": "string", "default": "output.csv"},
                },
                "required": ["input_file"],
            },
        },
        {
            "language": "bash",
            "name": "File Processing Template",
            "description": "Basic file processing script",
            "script_content": """#!/bin/bash

# Load parameters
PARAMS_FILE=$2
INPUT_FILE=$(jq -r '.input_file' $PARAMS_FILE)
OUTPUT_FILE=$(jq -r '.output_file // "output.txt"' $PARAMS_FILE)

# Create output directory
mkdir -p output

# Process file
wc -l $INPUT_FILE > output/$OUTPUT_FILE

# Return results
echo "{\"output_file\": \"$OUTPUT_FILE\"}"
""",
            "parameters_schema": {
                "type": "object",
                "properties": {
                    "input_file": {"type": "string"},
                    "output_file": {"type": "string", "default": "output.txt"},
                },
                "required": ["input_file"],
            },
        },
    ]

    return templates

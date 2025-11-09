"""Pipeline template API endpoints."""

from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel

from app.api.dependencies.auth import get_current_user
from app.models.user import User
from app.services import pipeline_templates
from app.services.pipeline_templates import TemplateCategory

router = APIRouter()


class PipelineTemplateResponse(BaseModel):
    """Response model for pipeline template."""

    id: str
    name: str
    description: str
    category: str
    steps: List[dict]
    parameters: dict
    inputs: List[str]
    outputs: List[str]


class TemplateCategoryResponse(BaseModel):
    """Response model for template categories."""

    categories: List[str]


@router.get(
    "/",
    summary="List all pipeline templates",
    response_model=List[PipelineTemplateResponse],
)
async def list_pipeline_templates(
    category: Optional[str] = Query(None, description="Filter by category"),
    search: Optional[str] = Query(None, description="Search in name/description"),
    current_user: User = Depends(get_current_user),
):
    """
    List all available pipeline templates.
    These are common, pre-configured workflows available to all users.

    Supports filtering by category and text search.
    """
    if category:
        templates = pipeline_templates.get_pipelines_by_category(category)
    else:
        templates = pipeline_templates.get_all_pipeline_templates()

    # Apply search filter if provided
    if search:
        search_lower = search.lower()
        templates = [
            t
            for t in templates
            if search_lower in t["name"].lower()
            or search_lower in t["description"].lower()
        ]

    return templates


@router.get(
    "/categories",
    summary="Get available template categories",
    response_model=TemplateCategoryResponse,
)
async def get_template_categories(
    current_user: User = Depends(get_current_user),
):
    """Get list of all template categories."""
    categories = [cat.value for cat in TemplateCategory]
    return {"categories": categories}


@router.get(
    "/{pipeline_id}",
    summary="Get pipeline template details",
    response_model=PipelineTemplateResponse,
)
async def get_pipeline_template(
    pipeline_id: str,
    current_user: User = Depends(get_current_user),
):
    """Get detailed information about a specific pipeline template."""
    template = pipeline_templates.get_pipeline_template(pipeline_id)
    if template is None:
        raise HTTPException(status_code=404, detail="Pipeline template not found")

    return template


@router.get("/{pipeline_id}/validate", summary="Validate template parameters")
async def validate_template_parameters(
    pipeline_id: str,
    current_user: User = Depends(get_current_user),
):
    """
    Validate that all required parameters for a template are understood.
    Returns parameter schema and validation rules.
    """
    template = pipeline_templates.get_pipeline_template(pipeline_id)
    if template is None:
        raise HTTPException(status_code=404, detail="Pipeline template not found")

    parameters = template.get("parameters", {})

    # Check required parameters
    required = {k: v for k, v in parameters.items() if v.get("required", False)}
    optional = {k: v for k, v in parameters.items() if not v.get("required", False)}

    return {
        "template_id": pipeline_id,
        "required_parameters": required,
        "optional_parameters": optional,
        "total_parameters": len(parameters),
    }

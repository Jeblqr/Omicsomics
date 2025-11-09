"""Pipeline template API endpoints."""

from fastapi import APIRouter, Depends, HTTPException

from app.api.dependencies.auth import get_current_user
from app.models.user import User
from app.services import pipeline_templates

router = APIRouter()


@router.get("/", summary="List all pipeline templates")
async def list_pipeline_templates(
    category: str | None = None,
    current_user: User = Depends(get_current_user),
):
    """
    List all available pipeline templates.
    These are common, pre-configured workflows available to all users.
    """
    if category:
        templates = pipeline_templates.get_pipelines_by_category(category)
    else:
        templates = pipeline_templates.get_all_pipeline_templates()
    
    return templates


@router.get("/{pipeline_id}", summary="Get pipeline template details")
async def get_pipeline_template(
    pipeline_id: str,
    current_user: User = Depends(get_current_user),
):
    """Get detailed information about a specific pipeline template."""
    template = pipeline_templates.get_pipeline_template(pipeline_id)
    if template is None:
        raise HTTPException(status_code=404, detail="Pipeline template not found")
    
    return template

"""Custom pipeline API endpoints."""

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.auth import get_current_user
from app.api.dependencies.database import get_async_db
from app.models.user import User
from app.services import custom_pipelines as pipeline_service

router = APIRouter()


class PipelineNode(BaseModel):
    """Pipeline node/step definition."""

    id: str
    type: str
    label: str
    data: dict
    position: dict | None = None


class PipelineEdge(BaseModel):
    """Connection between pipeline nodes."""

    id: str
    source: str
    target: str
    sourceHandle: str | None = None
    targetHandle: str | None = None


class PipelineDefinition(BaseModel):
    """Complete pipeline definition."""

    nodes: list[PipelineNode]
    edges: list[PipelineEdge]
    parameters: dict = {}


class CustomPipelineCreate(BaseModel):
    name: str
    description: str
    definition: PipelineDefinition
    category: str = "custom"
    is_public: bool = False
    template_id: str | None = None


class CustomPipelineUpdate(BaseModel):
    name: str | None = None
    description: str | None = None
    definition: PipelineDefinition | None = None
    category: str | None = None
    is_public: bool | None = None


class MergePipelinesRequest(BaseModel):
    """Request to merge multiple pipelines."""

    pipeline_ids: list[int] = []
    template_ids: list[str] = []
    custom_definitions: list[PipelineDefinition] = []


@router.post("/", status_code=status.HTTP_201_CREATED)
async def create_custom_pipeline(
    pipeline_in: CustomPipelineCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Create a new custom pipeline."""
    pipeline = await pipeline_service.create_custom_pipeline(
        db=db,
        name=pipeline_in.name,
        description=pipeline_in.description,
        definition=pipeline_in.definition.model_dump(),
        owner_id=current_user.id,
        category=pipeline_in.category,
        is_public=pipeline_in.is_public,
        template_id=pipeline_in.template_id,
    )

    return {
        "id": pipeline.id,
        "name": pipeline.name,
        "description": pipeline.description,
        "definition": pipeline.definition,
        "category": pipeline.category,
        "is_public": pipeline.is_public,
        "template_id": pipeline.template_id,
        "owner_id": pipeline.owner_id,
        "created_at": pipeline.created_at,
        "updated_at": pipeline.updated_at,
    }


@router.get("/")
async def list_custom_pipelines(
    category: str | None = None,
    include_public: bool = True,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """List custom pipelines for the current user."""
    pipelines = await pipeline_service.list_custom_pipelines(
        db=db,
        owner_id=current_user.id,
        include_public=include_public,
        category=category,
    )

    return [
        {
            "id": p.id,
            "name": p.name,
            "description": p.description,
            "definition": p.definition,
            "category": p.category,
            "is_public": p.is_public,
            "template_id": p.template_id,
            "owner_id": p.owner_id,
            "created_at": p.created_at,
            "updated_at": p.updated_at,
        }
        for p in pipelines
    ]


@router.get("/{pipeline_id}")
async def get_custom_pipeline(
    pipeline_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Get a specific custom pipeline."""
    pipeline = await pipeline_service.get_custom_pipeline(db, pipeline_id)

    if pipeline is None:
        raise HTTPException(status_code=404, detail="Pipeline not found")

    # Check access: owner or public
    if pipeline.owner_id != current_user.id and not pipeline.is_public:
        raise HTTPException(status_code=403, detail="Not authorized")

    return {
        "id": pipeline.id,
        "name": pipeline.name,
        "description": pipeline.description,
        "definition": pipeline.definition,
        "category": pipeline.category,
        "is_public": pipeline.is_public,
        "template_id": pipeline.template_id,
        "owner_id": pipeline.owner_id,
        "created_at": pipeline.created_at,
        "updated_at": pipeline.updated_at,
    }


@router.put("/{pipeline_id}")
async def update_custom_pipeline(
    pipeline_id: int,
    pipeline_update: CustomPipelineUpdate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Update a custom pipeline."""
    pipeline = await pipeline_service.get_custom_pipeline(db, pipeline_id)

    if pipeline is None:
        raise HTTPException(status_code=404, detail="Pipeline not found")

    if pipeline.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    definition_dict = (
        pipeline_update.definition.model_dump() if pipeline_update.definition else None
    )

    updated = await pipeline_service.update_custom_pipeline(
        db=db,
        pipeline_id=pipeline_id,
        name=pipeline_update.name,
        description=pipeline_update.description,
        definition=definition_dict,
        category=pipeline_update.category,
        is_public=pipeline_update.is_public,
    )

    return {
        "id": updated.id,
        "name": updated.name,
        "description": updated.description,
        "definition": updated.definition,
        "category": updated.category,
        "is_public": updated.is_public,
        "template_id": updated.template_id,
        "owner_id": updated.owner_id,
        "created_at": updated.created_at,
        "updated_at": updated.updated_at,
    }


@router.delete("/{pipeline_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_custom_pipeline(
    pipeline_id: int,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """Delete a custom pipeline."""
    pipeline = await pipeline_service.get_custom_pipeline(db, pipeline_id)

    if pipeline is None:
        raise HTTPException(status_code=404, detail="Pipeline not found")

    if pipeline.owner_id != current_user.id:
        raise HTTPException(status_code=403, detail="Not authorized")

    await pipeline_service.delete_custom_pipeline(db, pipeline_id)
    return None


@router.post("/merge")
async def merge_pipelines(
    request: MergePipelinesRequest,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_db),
):
    """
    Merge multiple pipelines into a single workflow.
    Can merge custom pipelines, template pipelines, or inline definitions.
    """
    from app.services import pipeline_templates

    pipeline_defs = []

    # Collect custom pipelines
    for pid in request.pipeline_ids:
        pipeline = await pipeline_service.get_custom_pipeline(db, pid)
        if pipeline is None:
            raise HTTPException(status_code=404, detail=f"Pipeline {pid} not found")
        if pipeline.owner_id != current_user.id and not pipeline.is_public:
            raise HTTPException(
                status_code=403, detail=f"Not authorized to access pipeline {pid}"
            )
        pipeline_defs.append(pipeline.definition)

    # Collect template pipelines
    for template_id in request.template_ids:
        template = pipeline_templates.get_pipeline_template(template_id)
        if template is None:
            raise HTTPException(
                status_code=404, detail=f"Template {template_id} not found"
            )

        # Convert template to pipeline definition format
        template_def = {
            "id": template["id"],
            "name": template["name"],
            "nodes": [
                {
                    "id": f"step_{i}",
                    "type": "process",
                    "label": step,
                    "data": {"step": step},
                }
                for i, step in enumerate(template["steps"])
            ],
            "edges": [
                {
                    "id": f"edge_{i}",
                    "source": f"step_{i}",
                    "target": f"step_{i+1}",
                }
                for i in range(len(template["steps"]) - 1)
            ],
            "parameters": {
                "inputs": template["inputs"],
                "outputs": template["outputs"],
                "category": template["category"],
            },
        }
        pipeline_defs.append(template_def)

    # Add custom inline definitions
    for custom_def in request.custom_definitions:
        pipeline_defs.append(custom_def.model_dump())

    if not pipeline_defs:
        raise HTTPException(status_code=400, detail="No pipelines provided to merge")

    # Perform merge
    merged_definition = pipeline_service.merge_pipelines(pipeline_defs)

    return {
        "merged_definition": merged_definition,
        "source_count": len(pipeline_defs),
        "total_nodes": len(merged_definition["nodes"]),
        "total_edges": len(merged_definition["edges"]),
    }

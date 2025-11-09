"""Custom pipeline service for CRUD operations."""

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.custom_pipeline import CustomPipeline


async def create_custom_pipeline(
    db: AsyncSession,
    name: str,
    description: str,
    definition: dict,
    owner_id: int,
    category: str = "custom",
    is_public: bool = False,
    template_id: str | None = None,
) -> CustomPipeline:
    """Create a new custom pipeline."""
    pipeline = CustomPipeline(
        name=name,
        description=description,
        definition=definition,
        owner_id=owner_id,
        category=category,
        is_public=is_public,
        template_id=template_id,
    )
    db.add(pipeline)
    await db.commit()
    await db.refresh(pipeline)
    return pipeline


async def get_custom_pipeline(
    db: AsyncSession, pipeline_id: int
) -> CustomPipeline | None:
    """Get a custom pipeline by ID."""
    result = await db.execute(
        select(CustomPipeline).where(CustomPipeline.id == pipeline_id)
    )
    return result.scalar_one_or_none()


async def list_custom_pipelines(
    db: AsyncSession,
    owner_id: int | None = None,
    include_public: bool = True,
    category: str | None = None,
) -> list[CustomPipeline]:
    """List custom pipelines with optional filtering."""
    query = select(CustomPipeline)

    if owner_id is not None:
        if include_public:
            # Get user's own pipelines and public ones
            query = query.where(
                (CustomPipeline.owner_id == owner_id)
                | (CustomPipeline.is_public == True)
            )
        else:
            # Only user's own pipelines
            query = query.where(CustomPipeline.owner_id == owner_id)
    elif include_public:
        # Only public pipelines
        query = query.where(CustomPipeline.is_public == True)

    if category is not None:
        query = query.where(CustomPipeline.category == category)

    query = query.order_by(CustomPipeline.created_at.desc())

    result = await db.execute(query)
    return list(result.scalars().all())


async def update_custom_pipeline(
    db: AsyncSession,
    pipeline_id: int,
    name: str | None = None,
    description: str | None = None,
    definition: dict | None = None,
    category: str | None = None,
    is_public: bool | None = None,
) -> CustomPipeline | None:
    """Update a custom pipeline."""
    pipeline = await get_custom_pipeline(db, pipeline_id)
    if pipeline is None:
        return None

    if name is not None:
        pipeline.name = name
    if description is not None:
        pipeline.description = description
    if definition is not None:
        pipeline.definition = definition
    if category is not None:
        pipeline.category = category
    if is_public is not None:
        pipeline.is_public = is_public

    await db.commit()
    await db.refresh(pipeline)
    return pipeline


async def delete_custom_pipeline(db: AsyncSession, pipeline_id: int) -> bool:
    """Delete a custom pipeline."""
    pipeline = await get_custom_pipeline(db, pipeline_id)
    if pipeline is None:
        return False

    await db.delete(pipeline)
    await db.commit()
    return True


def merge_pipelines(pipeline_definitions: list[dict]) -> dict:
    """
    Merge multiple pipeline definitions into a single workflow.

    Returns a merged pipeline definition with:
    - Combined nodes from all pipelines
    - Connected edges between pipelines
    - Merged parameters
    """
    if not pipeline_definitions:
        return {"nodes": [], "edges": [], "parameters": {}}

    merged = {
        "nodes": [],
        "edges": [],
        "parameters": {},
        "metadata": {
            "merged_from": [p.get("id") or p.get("name") for p in pipeline_definitions],
            "merge_strategy": "sequential",
        },
    }

    node_id_offset = 0
    previous_output_nodes = []

    for i, pipeline_def in enumerate(pipeline_definitions):
        nodes = pipeline_def.get("nodes", [])
        edges = pipeline_def.get("edges", [])
        params = pipeline_def.get("parameters", {})

        # Re-number node IDs to avoid conflicts
        id_mapping = {}
        for node in nodes:
            old_id = node.get("id")
            new_id = f"p{i}_n{node_id_offset}"
            id_mapping[old_id] = new_id
            node_id_offset += 1

            new_node = node.copy()
            new_node["id"] = new_id
            new_node["pipeline_index"] = i
            merged["nodes"].append(new_node)

        # Re-map edges with new node IDs
        for edge in edges:
            new_edge = {
                "id": f"p{i}_e{len(merged['edges'])}",
                "source": id_mapping.get(edge.get("source"), edge.get("source")),
                "target": id_mapping.get(edge.get("target"), edge.get("target")),
                "sourceHandle": edge.get("sourceHandle"),
                "targetHandle": edge.get("targetHandle"),
            }
            merged["edges"].append(new_edge)

        # Connect to previous pipeline's output nodes
        if i > 0 and previous_output_nodes:
            # Find input nodes of current pipeline
            input_node_ids = [
                n["id"]
                for n in merged["nodes"]
                if n.get("pipeline_index") == i and n.get("type") == "input"
            ]

            if not input_node_ids:
                # If no explicit input nodes, connect to first node
                input_node_ids = [
                    merged["nodes"][j]["id"]
                    for j, n in enumerate(merged["nodes"])
                    if n.get("pipeline_index") == i
                ][:1]

            # Create bridge edges
            for out_node in previous_output_nodes:
                for in_node in input_node_ids:
                    merged["edges"].append(
                        {
                            "id": f"bridge_p{i-1}_to_p{i}_{len(merged['edges'])}",
                            "source": out_node,
                            "target": in_node,
                            "type": "bridge",
                        }
                    )

        # Find output nodes for next iteration
        previous_output_nodes = [
            n["id"]
            for n in merged["nodes"]
            if n.get("pipeline_index") == i and n.get("type") == "output"
        ]
        if not previous_output_nodes:
            # If no explicit output nodes, use last node
            pipeline_nodes = [
                n["id"] for n in merged["nodes"] if n.get("pipeline_index") == i
            ]
            if pipeline_nodes:
                previous_output_nodes = [pipeline_nodes[-1]]

        # Merge parameters with namespace
        for key, value in params.items():
            namespaced_key = f"pipeline_{i}_{key}"
            merged["parameters"][namespaced_key] = value

    return merged

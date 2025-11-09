"""Custom pipeline service for CRUD operations."""

from typing import Optional
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.custom_pipeline import CustomPipeline
from app.core.tools_registry import get_tools_registry
from app.schemas.tools import NodeToolConfig, ValidationResult


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


def validate_pipeline_definition(definition: dict) -> ValidationResult:
    """
    Validate a pipeline definition including tool configurations.

    Checks:
    - Pipeline structure (nodes, edges)
    - Tool configurations for each node
    - Edge connectivity
    - Parameter dependencies
    """
    errors = []
    warnings = []

    nodes = definition.get("nodes", [])
    edges = definition.get("edges", [])

    if not nodes:
        errors.append("Pipeline must contain at least one node")
        return ValidationResult(valid=False, errors=errors, warnings=warnings)

    # Validate nodes
    node_ids = set()
    tool_nodes = []

    for i, node in enumerate(nodes):
        node_id = node.get("id")
        if not node_id:
            errors.append(f"Node at index {i} is missing 'id' field")
            continue

        if node_id in node_ids:
            errors.append(f"Duplicate node ID: {node_id}")
        node_ids.add(node_id)

        node_type = node.get("type")
        if not node_type:
            warnings.append(f"Node {node_id} is missing 'type' field")

        # Check for tool configuration
        node_data = node.get("data", {})
        tool_config_data = node_data.get("tool_config")

        if tool_config_data and node_type in ["process", "tool", "analysis"]:
            tool_nodes.append({"node_id": node_id, "config": tool_config_data})

    # Validate tool configurations
    registry = get_tools_registry()

    for tool_node in tool_nodes:
        node_id = tool_node["node_id"]
        config_data = tool_node["config"]

        try:
            # Convert dict to NodeToolConfig
            config = NodeToolConfig(**config_data)

            # Validate tool configuration
            tool = registry.get_tool(config.tool_id)
            if not tool:
                errors.append(
                    f"Node {node_id}: Tool '{config.tool_id}' not found in registry"
                )
                continue

            # Check required parameters
            required_params = [p.name for p in tool.parameters if p.required]
            for param_name in required_params:
                if param_name not in config.parameters:
                    errors.append(
                        f"Node {node_id}: Required parameter '{param_name}' is missing"
                    )

            # Check required inputs are mapped
            required_inputs = [i.name for i in tool.inputs if i.required]
            for input_name in required_inputs:
                if input_name not in config.input_mappings:
                    warnings.append(
                        f"Node {node_id}: Required input '{input_name}' is not mapped. "
                        "This may be connected via edges."
                    )

        except Exception as e:
            errors.append(f"Node {node_id}: Invalid tool configuration - {str(e)}")

    # Validate edges
    for i, edge in enumerate(edges):
        edge_id = edge.get("id")
        source = edge.get("source")
        target = edge.get("target")

        if not source:
            errors.append(f"Edge at index {i} is missing 'source' field")
        elif source not in node_ids:
            errors.append(f"Edge {edge_id or i}: Source node '{source}' not found")

        if not target:
            errors.append(f"Edge at index {i} is missing 'target' field")
        elif target not in node_ids:
            errors.append(f"Edge {edge_id or i}: Target node '{target}' not found")

    # Check for cycles (simple check)
    if has_cycle(nodes, edges):
        warnings.append("Pipeline contains cycles. This may cause execution issues.")

    # Check for disconnected nodes
    connected_nodes = set()
    for edge in edges:
        connected_nodes.add(edge.get("source"))
        connected_nodes.add(edge.get("target"))

    disconnected = node_ids - connected_nodes
    if len(disconnected) > 1:  # Allow one disconnected node (could be start/end)
        warnings.append(f"Pipeline has disconnected nodes: {', '.join(disconnected)}")

    return ValidationResult(valid=len(errors) == 0, errors=errors, warnings=warnings)


def has_cycle(nodes: list[dict], edges: list[dict]) -> bool:
    """
    Check if the pipeline graph has cycles using DFS.
    """
    # Build adjacency list
    graph = {}
    for node in nodes:
        graph[node.get("id")] = []

    for edge in edges:
        source = edge.get("source")
        target = edge.get("target")
        if source and target:
            if source not in graph:
                graph[source] = []
            graph[source].append(target)

    # DFS cycle detection
    visited = set()
    rec_stack = set()

    def dfs(node_id):
        visited.add(node_id)
        rec_stack.add(node_id)

        for neighbor in graph.get(node_id, []):
            if neighbor not in visited:
                if dfs(neighbor):
                    return True
            elif neighbor in rec_stack:
                return True

        rec_stack.remove(node_id)
        return False

    for node_id in graph:
        if node_id not in visited:
            if dfs(node_id):
                return True

    return False


def get_pipeline_dependencies(definition: dict) -> dict:
    """
    Analyze pipeline dependencies including tools and data formats.

    Returns:
    - tools: List of required tools with versions
    - formats: List of required data formats
    - resources: Aggregated resource requirements
    """
    registry = get_tools_registry()

    tools_required = {}
    formats_required = set()
    total_resources = {
        "min_cpu": 0,
        "max_cpu": 0,
        "min_memory_gb": 0,
        "max_memory_gb": 0,
        "gpu_required": False,
    }

    nodes = definition.get("nodes", [])

    for node in nodes:
        node_data = node.get("data", {})
        tool_config_data = node_data.get("tool_config")

        if tool_config_data:
            tool_id = tool_config_data.get("tool_id")
            tool_version = tool_config_data.get("tool_version")

            if tool_id:
                tool = registry.get_tool(tool_id)
                if tool:
                    tools_required[tool_id] = {
                        "name": tool.name,
                        "version": tool_version or tool.current_version,
                        "docker_image": tool.docker_image,
                        "category": tool.category.value,
                    }

                    # Collect input/output formats
                    for input_def in tool.inputs:
                        if input_def.format:
                            formats_required.add(input_def.format)
                    for output_def in tool.outputs:
                        if output_def.format:
                            formats_required.add(output_def.format)

                    # Aggregate resources
                    if tool.resources.min_cpu:
                        total_resources["min_cpu"] += tool.resources.min_cpu
                    if tool.resources.max_cpu:
                        total_resources["max_cpu"] += tool.resources.max_cpu
                    if tool.resources.min_memory_gb:
                        total_resources["min_memory_gb"] += tool.resources.min_memory_gb
                    if tool.resources.max_memory_gb:
                        total_resources["max_memory_gb"] += tool.resources.max_memory_gb
                    if tool.resources.gpu_required:
                        total_resources["gpu_required"] = True

    return {
        "tools": list(tools_required.values()),
        "formats": sorted(list(formats_required)),
        "resources": total_resources,
        "node_count": len(nodes),
    }

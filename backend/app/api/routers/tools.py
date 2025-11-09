"""Tools management API endpoints."""

from typing import Optional
from fastapi import APIRouter, HTTPException, Query
from app.core.tools_registry import get_tools_registry
from app.schemas.tools import (
    ToolSpecification,
    ToolCategory,
    OmicsType,
    ToolSearchQuery,
    NodeToolConfig,
    ValidationResult,
)

router = APIRouter()


@router.get("/", response_model=list[ToolSpecification])
async def list_tools(
    category: Optional[ToolCategory] = Query(
        None, description="Filter by tool category"
    ),
    omics_type: Optional[OmicsType] = Query(None, description="Filter by omics type"),
    tags: Optional[str] = Query(None, description="Filter by tags (comma-separated)"),
):
    """
    List all available tools with optional filters.

    - **category**: Filter tools by category (alignment, variant_calling, etc.)
    - **omics_type**: Filter tools by supported omics type
    - **tags**: Filter tools by tags (comma-separated list)
    """
    registry = get_tools_registry()

    tag_list = None
    if tags:
        tag_list = [t.strip() for t in tags.split(",")]

    tools = registry.list_tools(
        category=category,
        omics_type=omics_type,
        tags=tag_list,
    )

    return tools


@router.get("/search", response_model=list[ToolSpecification])
async def search_tools(
    q: str = Query(..., description="Search query string"),
    category: Optional[ToolCategory] = Query(None, description="Filter by category"),
    omics_type: Optional[OmicsType] = Query(None, description="Filter by omics type"),
    limit: int = Query(50, ge=1, le=1000, description="Maximum number of results"),
):
    """
    Search tools by name, description, or tags.

    - **q**: Search query (searches in tool name, description, and tags)
    - **category**: Optional category filter
    - **omics_type**: Optional omics type filter
    - **limit**: Maximum number of results to return
    """
    registry = get_tools_registry()

    # First, search by query
    results = registry.search_tools(q)

    # Apply additional filters
    if category:
        results = [t for t in results if t.category == category]

    if omics_type:
        results = [t for t in results if omics_type in t.omics_types]

    # Limit results
    return results[:limit]


@router.get("/categories", response_model=list[str])
async def list_categories():
    """
    List all available tool categories.
    """
    return [category.value for category in ToolCategory]


@router.get("/omics-types", response_model=list[str])
async def list_omics_types():
    """
    List all supported omics types.
    """
    return [omics_type.value for omics_type in OmicsType]


@router.get("/{tool_id}", response_model=ToolSpecification)
async def get_tool(tool_id: str):
    """
    Get detailed information about a specific tool.

    Returns the complete tool specification including:
    - Input/output definitions
    - Parameter specifications with validation rules
    - Resource requirements
    - Documentation and citation information
    """
    registry = get_tools_registry()
    tool = registry.get_tool(tool_id)

    if tool is None:
        raise HTTPException(status_code=404, detail=f"Tool '{tool_id}' not found")

    return tool


@router.get("/{tool_id}/parameters")
async def get_tool_parameters(tool_id: str):
    """
    Get parameter definitions for a specific tool.

    Returns structured parameter information for building UI forms.
    """
    registry = get_tools_registry()
    tool = registry.get_tool(tool_id)

    if tool is None:
        raise HTTPException(status_code=404, detail=f"Tool '{tool_id}' not found")

    return {
        "tool_id": tool.id,
        "tool_name": tool.name,
        "parameters": [
            {
                "name": p.name,
                "type": p.type,
                "label": p.label,
                "description": p.description,
                "required": p.required,
                "default": p.default,
                "enum_values": p.enum_values,
                "min_value": p.min_value,
                "max_value": p.max_value,
                "ui_widget": p.ui_widget,
                "group": p.group,
                "advanced": p.advanced,
            }
            for p in tool.parameters
        ],
    }


@router.get("/{tool_id}/inputs")
async def get_tool_inputs(tool_id: str):
    """
    Get input definitions for a specific tool.
    """
    registry = get_tools_registry()
    tool = registry.get_tool(tool_id)

    if tool is None:
        raise HTTPException(status_code=404, detail=f"Tool '{tool_id}' not found")

    return {
        "tool_id": tool.id,
        "tool_name": tool.name,
        "inputs": [
            {
                "name": i.name,
                "label": i.label,
                "description": i.description,
                "type": i.type,
                "required": i.required,
                "multiple": i.multiple,
                "file_extensions": i.file_extensions,
                "format": i.format,
            }
            for i in tool.inputs
        ],
    }


@router.get("/{tool_id}/outputs")
async def get_tool_outputs(tool_id: str):
    """
    Get output definitions for a specific tool.
    """
    registry = get_tools_registry()
    tool = registry.get_tool(tool_id)

    if tool is None:
        raise HTTPException(status_code=404, detail=f"Tool '{tool_id}' not found")

    return {
        "tool_id": tool.id,
        "tool_name": tool.name,
        "outputs": [
            {
                "name": o.name,
                "label": o.label,
                "description": o.description,
                "type": o.type,
                "file_extension": o.file_extension,
                "format": o.format,
                "optional": o.optional,
            }
            for o in tool.outputs
        ],
    }


@router.post("/validate", response_model=ValidationResult)
async def validate_tool_config(config: NodeToolConfig):
    """
    Validate a tool configuration for a pipeline node.

    Checks:
    - Tool exists in registry
    - Required parameters are provided
    - Parameter values match type constraints
    - Required inputs are mapped
    - Resource requirements are reasonable
    """
    registry = get_tools_registry()
    tool = registry.get_tool(config.tool_id)

    if tool is None:
        return ValidationResult(
            valid=False, errors=[f"Tool '{config.tool_id}' not found in registry"]
        )

    errors = []
    warnings = []

    # Check tool version
    if config.tool_version:
        valid_versions = [v.version for v in tool.versions]
        if config.tool_version not in valid_versions:
            warnings.append(
                f"Tool version '{config.tool_version}' not found. "
                f"Available versions: {', '.join(valid_versions)}"
            )

    # Validate required parameters
    required_params = [p.name for p in tool.parameters if p.required]
    for param_name in required_params:
        if param_name not in config.parameters:
            errors.append(f"Required parameter '{param_name}' is missing")

    # Validate parameter types and constraints
    param_map = {p.name: p for p in tool.parameters}
    for param_name, param_value in config.parameters.items():
        if param_name not in param_map:
            warnings.append(f"Unknown parameter '{param_name}'")
            continue

        param_def = param_map[param_name]

        # Type checking
        if param_def.type == "integer" and not isinstance(param_value, int):
            errors.append(f"Parameter '{param_name}' must be an integer")
        elif param_def.type == "float" and not isinstance(param_value, (int, float)):
            errors.append(f"Parameter '{param_name}' must be a number")
        elif param_def.type == "boolean" and not isinstance(param_value, bool):
            errors.append(f"Parameter '{param_name}' must be a boolean")
        elif param_def.type == "string" and not isinstance(param_value, str):
            errors.append(f"Parameter '{param_name}' must be a string")

        # Range checking for numeric types
        if param_def.type in ["integer", "float"] and isinstance(
            param_value, (int, float)
        ):
            if param_def.min_value is not None and param_value < param_def.min_value:
                errors.append(
                    f"Parameter '{param_name}' value {param_value} is below minimum {param_def.min_value}"
                )
            if param_def.max_value is not None and param_value > param_def.max_value:
                errors.append(
                    f"Parameter '{param_name}' value {param_value} exceeds maximum {param_def.max_value}"
                )

        # Enum validation
        if param_def.type == "enum" and param_def.enum_values:
            if param_value not in param_def.enum_values:
                errors.append(
                    f"Parameter '{param_name}' value '{param_value}' not in allowed values: "
                    f"{', '.join(param_def.enum_values)}"
                )

    # Validate required inputs are mapped
    required_inputs = [i.name for i in tool.inputs if i.required]
    for input_name in required_inputs:
        if input_name not in config.input_mappings:
            errors.append(f"Required input '{input_name}' is not mapped")

    # Check for unknown input mappings
    valid_inputs = {i.name for i in tool.inputs}
    for input_name in config.input_mappings:
        if input_name not in valid_inputs:
            warnings.append(f"Unknown input '{input_name}' in mappings")

    # Validate resource requirements
    if config.resources:
        if config.resources.min_cpu and config.resources.max_cpu:
            if config.resources.min_cpu > config.resources.max_cpu:
                errors.append("Minimum CPU cannot exceed maximum CPU")

        if config.resources.min_memory_gb and config.resources.max_memory_gb:
            if config.resources.min_memory_gb > config.resources.max_memory_gb:
                errors.append("Minimum memory cannot exceed maximum memory")

    return ValidationResult(valid=len(errors) == 0, errors=errors, warnings=warnings)


@router.get("/{tool_id}/versions")
async def get_tool_versions(tool_id: str):
    """
    Get available versions for a specific tool.
    """
    registry = get_tools_registry()
    tool = registry.get_tool(tool_id)

    if tool is None:
        raise HTTPException(status_code=404, detail=f"Tool '{tool_id}' not found")

    return {
        "tool_id": tool.id,
        "current_version": tool.current_version,
        "versions": [
            {
                "version": v.version,
                "release_date": v.release_date,
                "deprecated": v.deprecated,
                "docker_image": v.docker_image,
                "conda_package": v.conda_package,
            }
            for v in tool.versions
        ],
    }


@router.get("/{tool_id}/resources")
async def get_tool_resources(tool_id: str):
    """
    Get resource requirements for a specific tool.
    """
    registry = get_tools_registry()
    tool = registry.get_tool(tool_id)

    if tool is None:
        raise HTTPException(status_code=404, detail=f"Tool '{tool_id}' not found")

    return {
        "tool_id": tool.id,
        "tool_name": tool.name,
        "resources": {
            "min_cpu": tool.resources.min_cpu,
            "max_cpu": tool.resources.max_cpu,
            "min_memory_gb": tool.resources.min_memory_gb,
            "max_memory_gb": tool.resources.max_memory_gb,
            "min_disk_gb": tool.resources.min_disk_gb,
            "gpu_required": tool.resources.gpu_required,
            "gpu_count": tool.resources.gpu_count,
        },
    }

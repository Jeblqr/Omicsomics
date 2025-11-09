"""Tool specification schemas for pipeline builder."""

from enum import Enum
from typing import Any, Literal
from pydantic import BaseModel, Field


class ParameterType(str, Enum):
    """Parameter data types."""

    STRING = "string"
    INTEGER = "integer"
    FLOAT = "float"
    BOOLEAN = "boolean"
    FILE = "file"
    DIRECTORY = "directory"
    ENUM = "enum"
    ARRAY = "array"
    OBJECT = "object"


class ResourceRequirement(BaseModel):
    """Computational resource requirements."""

    min_cpu: int | None = Field(None, description="Minimum CPU cores")
    max_cpu: int | None = Field(None, description="Maximum CPU cores")
    min_memory_gb: float | None = Field(None, description="Minimum memory in GB")
    max_memory_gb: float | None = Field(None, description="Maximum memory in GB")
    min_disk_gb: float | None = Field(None, description="Minimum disk space in GB")
    gpu_required: bool = Field(False, description="Whether GPU is required")
    gpu_count: int | None = Field(None, description="Number of GPUs needed")


class ParameterDefinition(BaseModel):
    """Definition of a tool parameter."""

    name: str = Field(..., description="Parameter name")
    type: ParameterType = Field(..., description="Parameter data type")
    label: str = Field(..., description="Human-readable label")
    description: str = Field("", description="Parameter description")
    required: bool = Field(False, description="Whether parameter is required")
    default: Any | None = Field(None, description="Default value")

    # For enum type
    enum_values: list[str] | None = Field(
        None, description="Allowed values for enum type"
    )

    # For numeric types
    min_value: float | None = Field(None, description="Minimum value")
    max_value: float | None = Field(None, description="Maximum value")

    # For string types
    pattern: str | None = Field(None, description="Regex pattern for validation")
    min_length: int | None = Field(None, description="Minimum string length")
    max_length: int | None = Field(None, description="Maximum string length")

    # For array types
    item_type: ParameterType | None = Field(None, description="Type of array items")
    min_items: int | None = Field(None, description="Minimum number of items")
    max_items: int | None = Field(None, description="Maximum number of items")

    # For file types
    file_extensions: list[str] | None = Field(
        None, description="Allowed file extensions"
    )
    mime_types: list[str] | None = Field(None, description="Allowed MIME types")

    # UI hints
    ui_widget: str | None = Field(
        None, description="Suggested UI widget (textbox, slider, dropdown, etc.)"
    )
    group: str | None = Field(None, description="Parameter group for organization")
    order: int = Field(0, description="Display order within group")
    advanced: bool = Field(False, description="Whether this is an advanced parameter")


class InputDefinition(BaseModel):
    """Definition of a tool input."""

    name: str = Field(..., description="Input name")
    label: str = Field(..., description="Human-readable label")
    description: str = Field("", description="Input description")
    type: ParameterType = Field(..., description="Input data type")
    required: bool = Field(True, description="Whether input is required")
    multiple: bool = Field(False, description="Whether multiple inputs are allowed")
    file_extensions: list[str] | None = Field(
        None, description="Allowed file extensions"
    )
    mime_types: list[str] | None = Field(None, description="Allowed MIME types")
    omics_types: list[str] | None = Field(None, description="Accepted omics data types")
    format: str | None = Field(
        None, description="Expected data format (e.g., 'vcf', 'bam', 'fastq')"
    )


class OutputDefinition(BaseModel):
    """Definition of a tool output."""

    name: str = Field(..., description="Output name")
    label: str = Field(..., description="Human-readable label")
    description: str = Field("", description="Output description")
    type: ParameterType = Field(..., description="Output data type")
    file_extension: str | None = Field(None, description="Output file extension")
    mime_type: str | None = Field(None, description="Output MIME type")
    format: str | None = Field(None, description="Output data format")
    optional: bool = Field(False, description="Whether output is optional")


class ToolVersion(BaseModel):
    """Tool version information."""

    version: str = Field(..., description="Version string")
    release_date: str | None = Field(None, description="Release date")
    changelog: str | None = Field(None, description="Version changelog")
    deprecated: bool = Field(False, description="Whether this version is deprecated")
    docker_image: str | None = Field(None, description="Docker image for this version")
    conda_package: str | None = Field(
        None, description="Conda package for this version"
    )


class ToolCategory(str, Enum):
    """Tool categories."""

    ALIGNMENT = "alignment"
    VARIANT_CALLING = "variant_calling"
    QUALITY_CONTROL = "quality_control"
    QUANTIFICATION = "quantification"
    DIFFERENTIAL_EXPRESSION = "differential_expression"
    FUNCTIONAL_ANALYSIS = "functional_analysis"
    VISUALIZATION = "visualization"
    PREPROCESSING = "preprocessing"
    NORMALIZATION = "normalization"
    STATISTICAL_ANALYSIS = "statistical_analysis"
    CLUSTERING = "clustering"
    DIMENSION_REDUCTION = "dimension_reduction"
    CELL_TYPE_ANNOTATION = "cell_type_annotation"
    PEAK_CALLING = "peak_calling"
    MOTIF_ANALYSIS = "motif_analysis"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    CUSTOM = "custom"


class OmicsType(str, Enum):
    """Omics data types."""

    GENOMICS = "genomics"
    TRANSCRIPTOMICS = "transcriptomics"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    EPIGENOMICS = "epigenomics"
    SINGLE_CELL = "single_cell"
    MULTIOMICS = "multiomics"


class ToolSpecification(BaseModel):
    """Complete tool specification."""

    id: str = Field(..., description="Unique tool identifier")
    name: str = Field(..., description="Tool name")
    description: str = Field(..., description="Tool description")
    category: ToolCategory = Field(..., description="Tool category")
    omics_types: list[OmicsType] = Field(..., description="Supported omics types")

    # Version information
    current_version: str = Field(..., description="Current/default version")
    versions: list[ToolVersion] = Field(
        default_factory=list, description="Available versions"
    )

    # Tool definition
    inputs: list[InputDefinition] = Field(
        default_factory=list, description="Input definitions"
    )
    outputs: list[OutputDefinition] = Field(
        default_factory=list, description="Output definitions"
    )
    parameters: list[ParameterDefinition] = Field(
        default_factory=list, description="Parameter definitions"
    )

    # Requirements
    resources: ResourceRequirement = Field(
        default_factory=ResourceRequirement, description="Resource requirements"
    )
    dependencies: list[str] = Field(
        default_factory=list, description="Tool dependencies"
    )

    # Documentation
    documentation_url: str | None = Field(
        None, description="Official documentation URL"
    )
    citation: str | None = Field(None, description="Citation information")
    license: str | None = Field(None, description="Software license")

    # Execution
    command_template: str | None = Field(None, description="Command line template")
    docker_image: str | None = Field(None, description="Docker image")
    conda_package: str | None = Field(None, description="Conda package")

    # Metadata
    author: str | None = Field(None, description="Tool author/organization")
    homepage: str | None = Field(None, description="Tool homepage URL")
    tags: list[str] = Field(default_factory=list, description="Search tags")


class NodeToolConfig(BaseModel):
    """Tool configuration for a pipeline node."""

    tool_id: str = Field(..., description="Tool identifier from registry")
    tool_version: str | None = Field(None, description="Specific tool version")
    parameters: dict[str, Any] = Field(
        default_factory=dict, description="Tool parameter values"
    )
    input_mappings: dict[str, str] = Field(
        default_factory=dict, description="Input port to data source mappings"
    )
    output_mappings: dict[str, str] = Field(
        default_factory=dict, description="Output port to destination mappings"
    )
    resources: ResourceRequirement | None = Field(
        None, description="Custom resource requirements"
    )
    environment: dict[str, str] = Field(
        default_factory=dict, description="Environment variables"
    )
    enabled: bool = Field(
        True, description="Whether this tool is enabled in the pipeline"
    )


class ValidationResult(BaseModel):
    """Result of parameter validation."""

    valid: bool = Field(..., description="Whether validation passed")
    errors: list[str] = Field(
        default_factory=list, description="Validation error messages"
    )
    warnings: list[str] = Field(default_factory=list, description="Validation warnings")


class ToolSearchQuery(BaseModel):
    """Tool search query parameters."""

    query: str | None = Field(None, description="Search query string")
    category: ToolCategory | None = Field(None, description="Filter by category")
    omics_type: OmicsType | None = Field(None, description="Filter by omics type")
    tags: list[str] = Field(default_factory=list, description="Filter by tags")
    limit: int = Field(50, ge=1, le=1000, description="Maximum results")
    offset: int = Field(0, ge=0, description="Results offset for pagination")

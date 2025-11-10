"""
Quick Visualizer API

REST API for quick data visualization from Data Browser.
Provides endpoints for:
- Analyzing files and suggesting visualizations
- Generating visualizations with custom options
- Listing available chart types
"""

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any, List
from pathlib import Path
import logging

from app.services.quick_visualizer import get_quick_visualizer, ChartType, DataType

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/quick-visualizer", tags=["quick-visualizer"])


# Pydantic models
class AnalyzeRequest(BaseModel):
    """Request to analyze a file"""

    file_path: str = Field(..., description="Absolute path to data file")


class AnalyzeResponse(BaseModel):
    """Analysis result"""

    success: bool
    error: Optional[str] = None
    data_type: Optional[str] = None
    n_rows: Optional[int] = None
    n_columns: Optional[int] = None
    limited: Optional[bool] = Field(
        None, description="Whether data was limited for performance"
    )
    column_info: Optional[List[Dict[str, Any]]] = None
    suggestions: Optional[List[Dict[str, Any]]] = Field(
        None, description="Suggested visualizations with options"
    )


class VisualizeRequest(BaseModel):
    """Request to generate visualization"""

    file_path: str = Field(..., description="Absolute path to data file")
    chart_type: str = Field(..., description="Type of chart to generate")
    options: Optional[Dict[str, Any]] = Field(
        default_factory=dict, description="Chart-specific options"
    )


class VisualizeResponse(BaseModel):
    """Visualization result"""

    success: bool
    error: Optional[str] = None
    image: Optional[str] = Field(None, description="Base64-encoded PNG image")
    chart_type: Optional[str] = None


class ChartTypeInfo(BaseModel):
    """Chart type information"""

    type: str
    name: str
    description: str
    required_options: List[str]
    optional_options: List[str]


# Endpoints
@router.post("/analyze", response_model=AnalyzeResponse)
async def analyze_file(request: AnalyzeRequest):
    """
    Analyze a data file and suggest visualizations

    Detects data type, analyzes columns, and suggests appropriate
    visualization methods.

    Args:
        request: File path to analyze

    Returns:
        Analysis results with suggestions
    """
    try:
        file_path = Path(request.file_path)

        if not file_path.exists():
            raise HTTPException(status_code=404, detail="File not found")

        visualizer = get_quick_visualizer()
        result = visualizer.analyze_file(file_path)

        return AnalyzeResponse(**result)

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to analyze file: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/visualize", response_model=VisualizeResponse)
async def generate_visualization(request: VisualizeRequest):
    """
    Generate a visualization for a data file

    Creates a chart based on the specified type and options.
    Returns the chart as a base64-encoded PNG image.

    Args:
        request: Visualization request with file path, chart type, and options

    Returns:
        Base64-encoded image data
    """
    try:
        file_path = Path(request.file_path)

        if not file_path.exists():
            raise HTTPException(status_code=404, detail="File not found")

        # Validate chart type
        try:
            ChartType(request.chart_type)
        except ValueError:
            raise HTTPException(
                status_code=400, detail=f"Invalid chart type: {request.chart_type}"
            )

        visualizer = get_quick_visualizer()
        result = visualizer.generate_visualization(
            file_path, request.chart_type, request.options
        )

        if not result.get("success"):
            raise HTTPException(
                status_code=400, detail=result.get("error", "Visualization failed")
            )

        return VisualizeResponse(**result)

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to generate visualization: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/chart-types", response_model=List[ChartTypeInfo])
async def list_chart_types():
    """
    List all available chart types

    Returns information about each chart type, including
    required and optional parameters.

    Returns:
        List of chart type definitions
    """
    chart_types = [
        {
            "type": ChartType.HISTOGRAM.value,
            "name": "Histogram",
            "description": "Distribution of a numeric variable",
            "required_options": [],
            "optional_options": ["column", "bins"],
        },
        {
            "type": ChartType.SCATTER.value,
            "name": "Scatter Plot",
            "description": "Relationship between two numeric variables",
            "required_options": [],
            "optional_options": ["x_column", "y_column"],
        },
        {
            "type": ChartType.LINE.value,
            "name": "Line Chart",
            "description": "Trends over time or sequence",
            "required_options": [],
            "optional_options": ["columns"],
        },
        {
            "type": ChartType.BAR.value,
            "name": "Bar Chart",
            "description": "Count or comparison of categories",
            "required_options": [],
            "optional_options": ["column"],
        },
        {
            "type": ChartType.BOXPLOT.value,
            "name": "Box Plot",
            "description": "Distribution and outliers of numeric variables",
            "required_options": [],
            "optional_options": ["columns"],
        },
        {
            "type": ChartType.HEATMAP.value,
            "name": "Heatmap",
            "description": "Visualize matrix data with colors",
            "required_options": [],
            "optional_options": [],
        },
        {
            "type": ChartType.VIOLIN.value,
            "name": "Violin Plot",
            "description": "Distribution shape with density estimation",
            "required_options": [],
            "optional_options": ["columns"],
        },
        {
            "type": ChartType.CORRELATION.value,
            "name": "Correlation Matrix",
            "description": "Pairwise correlations between numeric variables",
            "required_options": [],
            "optional_options": [],
        },
    ]

    return [ChartTypeInfo(**ct) for ct in chart_types]


@router.get("/data-types")
async def list_data_types():
    """
    List all detectable data types

    Returns:
        List of data type definitions
    """
    return {
        "data_types": [
            {
                "type": DataType.NUMERIC_MATRIX.value,
                "name": "Numeric Matrix",
                "description": "Matrix with mostly numeric values",
            },
            {
                "type": DataType.TIME_SERIES.value,
                "name": "Time Series",
                "description": "Data with temporal dimension",
            },
            {
                "type": DataType.CATEGORICAL.value,
                "name": "Categorical",
                "description": "Data with categorical variables",
            },
            {
                "type": DataType.MIXED.value,
                "name": "Mixed",
                "description": "Mix of numeric and categorical data",
            },
            {
                "type": DataType.GENOMIC_INTERVALS.value,
                "name": "Genomic Intervals",
                "description": "Genomic coordinates (chr, start, end)",
            },
            {
                "type": DataType.EXPRESSION_MATRIX.value,
                "name": "Expression Matrix",
                "description": "Gene expression or similar omics data",
            },
            {
                "type": DataType.UNKNOWN.value,
                "name": "Unknown",
                "description": "Data type could not be determined",
            },
        ]
    }

"""Visualization API endpoints."""

import logging
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.visualizations.generator import viz_generator
from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/visualizations", tags=["visualizations"])


# Request Models
class VolcanoPlotRequest(BaseModel):
    """Request for volcano plot generation."""

    results_file: str = Field(..., description="Path to DE results file")
    log2fc_col: str = Field(default="log2FoldChange", description="Log2FC column name")
    pvalue_col: str = Field(default="pvalue", description="P-value column name")
    padj_col: str = Field(default="padj", description="Adjusted p-value column name")
    gene_col: str = Field(default="gene", description="Gene identifier column")
    fc_threshold: float = Field(default=1.0, description="Log2FC threshold")
    pval_threshold: float = Field(default=0.05, description="Adjusted p-value threshold")


class UMAPPlotRequest(BaseModel):
    """Request for UMAP plot generation."""

    h5ad_file: str = Field(..., description="Path to AnnData H5AD file")
    color_by: str = Field(default="leiden", description="Column to color by")


class HeatmapRequest(BaseModel):
    """Request for heatmap generation."""

    h5ad_file: str = Field(..., description="Path to AnnData H5AD file")
    gene_list: list[str] = Field(..., description="List of genes to include")
    groupby: str = Field(default="leiden", description="Column to group by")
    standard_scale: str = Field(default="var", description="Scaling method")


# Endpoints
@router.post("/volcano")
async def get_volcano_plot(
    request: VolcanoPlotRequest,
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Generate volcano plot data from differential expression results.
    
    Returns JSON data formatted for Plotly visualization including:
    - X coordinates (log2 fold change)
    - Y coordinates (-log10 p-value)
    - Gene identifiers
    - Significance categories (up, down, not_sig)
    - Summary counts
    """
    result = viz_generator.generate_volcano_plot_data(
        results_file=request.results_file,
        log2fc_col=request.log2fc_col,
        pvalue_col=request.pvalue_col,
        padj_col=request.padj_col,
        gene_col=request.gene_col,
        fc_threshold=request.fc_threshold,
        pval_threshold=request.pval_threshold,
    )

    if result["status"] == "error":
        raise HTTPException(status_code=400, detail=result["error"])

    return result["data"]


@router.post("/pca")
async def get_pca_plot(
    h5ad_file: str = Query(..., description="Path to AnnData H5AD file"),
    color_by: str = Query(default="leiden", description="Column to color by"),
    n_components: int = Query(default=2, description="Number of components (2 or 3)"),
    current_user: User = Depends(get_current_user),
):
    """
    Generate PCA plot data from single-cell AnnData object.
    
    Returns PCA coordinates and metadata for interactive plotting.
    Supports both 2D and 3D visualization.
    """
    result = viz_generator.generate_pca_plot_data(
        h5ad_file=h5ad_file,
        color_by=color_by,
        n_components=n_components,
    )

    if result["status"] == "error":
        raise HTTPException(status_code=400, detail=result["error"])

    return result["data"]


@router.post("/umap")
async def get_umap_plot(
    request: UMAPPlotRequest,
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Generate UMAP plot data from single-cell AnnData object.
    
    Returns UMAP coordinates, cell metadata, and cluster counts
    for interactive visualization with hover tooltips and selection.
    """
    result = viz_generator.generate_umap_plot_data(
        h5ad_file=request.h5ad_file,
        color_by=request.color_by,
    )

    if result["status"] == "error":
        raise HTTPException(status_code=400, detail=result["error"])

    return result["data"]


@router.post("/heatmap")
async def get_heatmap(
    request: HeatmapRequest,
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Generate heatmap data for gene expression.
    
    Returns expression matrix (groups x genes) suitable for
    heatmap visualization with Plotly or similar libraries.
    """
    result = viz_generator.generate_heatmap_data(
        h5ad_file=request.h5ad_file,
        gene_list=request.gene_list,
        groupby=request.groupby,
        standard_scale=request.standard_scale,
    )

    if result["status"] == "error":
        raise HTTPException(status_code=400, detail=result["error"])

    return result["data"]


@router.get("/igv")
async def get_igv_config(
    bam_file: str = Query(..., description="Path to BAM file"),
    region: str = Query(..., description="Genomic region (chr:start-end)"),
    reference_genome: str = Query(default="hg38", description="Reference genome"),
    current_user: User = Depends(get_current_user),
):
    """
    Get IGV.js configuration for genome browser visualization.
    
    Returns configuration object that can be directly used to
    initialize IGV.js browser in the frontend.
    """
    result = viz_generator.generate_igv_data(
        bam_file=bam_file,
        region=region,
        reference_genome=reference_genome,
    )

    if result["status"] == "error":
        raise HTTPException(status_code=400, detail=result["error"])

    return result["config"]


@router.get("/quality-metrics")
async def get_quality_metrics(
    h5ad_file: str = Query(..., description="Path to AnnData H5AD file"),
    metrics: list[str] = Query(
        default=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        description="Metrics to retrieve",
    ),
    current_user: User = Depends(get_current_user),
):
    """
    Get quality control metrics for single-cell data.
    
    Returns distribution statistics for QC metrics including
    mean, median, standard deviation for visualization.
    """
    result = viz_generator.generate_quality_metrics_plot(
        h5ad_file=h5ad_file,
        metrics=metrics,
    )

    if result["status"] == "error":
        raise HTTPException(status_code=400, detail=result["error"])

    return result["data"]


@router.get("/formats")
async def get_supported_formats():
    """
    Get list of supported file formats for visualization.
    
    Returns information about file types that can be visualized
    and the corresponding visualization types available.
    """
    return {
        "formats": {
            "h5ad": {
                "description": "AnnData format for single-cell data",
                "visualizations": ["umap", "pca", "heatmap", "quality-metrics"],
            },
            "csv": {
                "description": "CSV/TSV tables for differential expression",
                "visualizations": ["volcano", "scatter"],
            },
            "bam": {
                "description": "Binary alignment format",
                "visualizations": ["igv", "coverage"],
            },
            "vcf": {
                "description": "Variant call format",
                "visualizations": ["igv", "variant-table"],
            },
        },
        "plot_types": [
            {"type": "volcano", "description": "Volcano plot for differential expression"},
            {"type": "umap", "description": "UMAP dimensionality reduction plot"},
            {"type": "pca", "description": "PCA plot"},
            {"type": "heatmap", "description": "Gene expression heatmap"},
            {"type": "igv", "description": "Genome browser (IGV.js)"},
            {"type": "quality-metrics", "description": "QC metrics distributions"},
        ],
    }
